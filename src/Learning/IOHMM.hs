{-# LANGUAGE RecordWildCards #-}

module Learning.IOHMM
  ( IOHMM (..)
  , LogLikelihood
  , init
  , withEmission
  , euclideanDistance
  , viterbi
  , baumWelch
  , baumWelch'
  , simulate
  ) where

import           Control.Applicative                         ( (<$>) )
import           Control.Arrow                               ( first )
import           Data.List                                   ( elemIndex )
import           Data.Maybe                                  ( fromJust )
import           Data.Random.Distribution                    ( rvar )
import qualified Data.Random.Distribution.Categorical as C   ( Categorical, fromList, normalizeCategoricalPs )
import           Data.Random.Distribution.Extra              ( pmf )
import           Data.Random.RVar                            ( RVar )
import qualified Data.Vector                          as V   ( elemIndex, fromList, map, toList, unsafeIndex )
import qualified Data.Vector.Generic                  as G   ( convert )
import qualified Data.Vector.Unboxed                  as U   ( fromList, zip )
import qualified Numeric.LinearAlgebra.Data           as H   ( (!), fromList, fromLists, toList )
import qualified Numeric.LinearAlgebra.HMatrix        as H   ( tr )
import           Learning.IOHMM.Internal                     ( LogLikelihood )
import qualified Learning.IOHMM.Internal              as I
import           Prelude                              hiding ( init )

-- | Parameter set of the input-output hidden Markov model with discrete emission.
--   This 'IOHMM' assumes that the inputs affect only the transition
--   probabilities. The model schema is as follows.
--
--   @
--       x_0    x_1           x_n
--               |             |
--               v             v
--       z_0 -> z_1 -> ... -> z_n
--        |      |             |
--        v      v             v
--       y_0    y_1           y_n
--   @
--
--   Here, @[x_0, x_1, ..., x_n]@ are given inputs, @[z_0, z_1, ..., z_n]@
--   are hidden states, and @[y_0, y_1, ..., y_n]@ are observed outputs.
--   @z_0@ is determined by the 'initialStateDist'.
--   For @i = 1, ..., n@, @z_i@ is determined by the 'transitionDist'
--   conditioned by @x_i@ and @z_{i-1}@.
--   For @i = 0, ..., n@, @y_i@ is determined by the 'emissionDist'
--   conditioned by @z_i@.
data IOHMM i s o = IOHMM { inputs  :: [i]
                         , states  :: [s]
                         , outputs :: [o]
                         , initialStateDist :: C.Categorical Double s
                           -- ^ Categorical distribution of initial state
                         , transitionDist :: i -> s -> C.Categorical Double s
                           -- ^ Categorical distribution of next state
                           --   conditioned by the input and previous state
                         , emissionDist :: s -> C.Categorical Double o
                           -- ^ Categorical distribution of output conditioned
                           --   by the hidden state
                         }

instance (Show i, Show s, Show o) => Show (IOHMM i s o) where
  show IOHMM {..} = "IOHMM {inputs = "           ++ show inputs
                      ++ ", states = "           ++ show states
                      ++ ", outputs = "          ++ show outputs
                      ++ ", initialStateDist = " ++ show initialStateDist
                      ++ ", transitionDist = "   ++ show [(transitionDist i s, (i, s)) | i <- inputs, s <- states]
                      ++ ", emissionDist = "     ++ show [(emissionDist s, s) | s <- states]
                      ++ "}"

-- | @init inputs states outputs@ returns a random variable of models with the
--   @inputs@, @states@, and @outputs@, wherein parameters are sampled from uniform
--   distributions.
init :: (Eq i, Eq s, Eq o) => [i] -> [s] -> [o] -> RVar (IOHMM i s o)
init is ss os = fromInternal is ss os <$> I.init (length is) (length ss) (length os)

-- | @withEmission model xs ys@ returns a model in which the
--   'emissionDist' is updated by re-estimations using the inputs @xs@ and
--   outputs @ys@. The 'emissionDist' is set to be normalized histograms
--   each of which is calculated from segumentations of @ys@ based on the
--   Viterbi state path.
--   If the lengths of @xs@ and @ys@ are different, the longer one is cut
--   by the length of the shorter one.
withEmission :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> IOHMM i s o
withEmission (model @ IOHMM {..}) xs ys = fromInternal inputs states outputs $ I.withEmission model' $ U.zip xs' ys'
  where
    inputs'  = V.fromList inputs
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` inputs') xs
    ys'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') ys

-- | Return the Euclidean distance between two models that have the same
--   inputs, states, and outputs.
euclideanDistance :: (Eq i, Eq s, Eq o) => IOHMM i s o -> IOHMM i s o -> Double
euclideanDistance model1 model2 =
  checkTwoModelsIn "euclideanDistance" model1 model2 `seq`
  I.euclideanDistance model1' model2'
  where
    model1' = toInternal model1
    model2' = toInternal model2

-- | @viterbi model xs ys@ performs the Viterbi algorithm using the inputs
--   @xs@ and outputs @ys@, and returns the most likely state path and its
--   log likelihood.
--   If the lengths of @xs@ and @ys@ are different, the longer one is cut
--   by the length of the shorter one.
viterbi :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> ([s], LogLikelihood)
viterbi (model @ IOHMM {..}) xs ys =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs ys `seq`
  first toStates $ I.viterbi model' $ U.zip xs' ys'
  where
    inputs'  = V.fromList inputs
    states'  = V.fromList states
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` inputs') xs
    ys'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') ys
    toStates = V.toList . V.map (V.unsafeIndex states') . G.convert

-- | @baumWelch model xs ys@ iteratively performs the Baum-Welch algorithm
--   using the inputs @xs@ and outputs @ys@, and returns a list of updated
--   models and their corresponding log likelihoods.
--   If the lengths of @xs@ and @ys@ are different, the longer one is cut
--   by the length of the shorter one.
baumWelch :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> [(IOHMM i s o, LogLikelihood)]
baumWelch (model @ IOHMM {..}) xs ys =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs ys `seq`
  map (first $ fromInternal inputs states outputs) $ I.baumWelch model' $ U.zip xs' ys'
  where
    inputs'  = V.fromList inputs
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` inputs') xs
    ys'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') ys

-- | @baumWelch' model xs@ performs the Baum-Welch algorithm using the
--   inputs @xs@ and outputs @ys@, and returns a model locally maximizing
--   its log likelihood.
baumWelch' :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> (IOHMM i s o, LogLikelihood)
baumWelch' (model @ IOHMM {..}) xs ys =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs ys `seq`
  first (fromInternal inputs states outputs) $ I.baumWelch' model' $ U.zip xs' ys'
  where
    inputs'  = V.fromList inputs
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` inputs') xs
    ys'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') ys

-- | @simulate model xs@ generates a Markov process coinciding with the
--   inputs @xs@ using the @model@, and returns its state path and observed
--   outputs.
simulate :: IOHMM i s o -> [i] -> RVar ([s], [o])
simulate IOHMM {..} xs
  | null xs   = return ([], [])
  | otherwise = do s0 <- rvar initialStateDist
                   y0 <- rvar $ emissionDist s0
                   unzip . ((s0, y0) :) <$> sim s0 (tail xs)
  where
    sim _ []      = return []
    sim s (x:xs') = do s' <- rvar $ transitionDist x s
                       y' <- rvar $ emissionDist s'
                       ((s', y') :) <$> sim s' xs'

-- | Check if the model is valid in the sense of whether the 'states' and
--   'outputs' are not empty.
checkModelIn :: String -> IOHMM i s o -> ()
checkModelIn fun IOHMM {..}
  | null inputs  = errorIn fun "empty inputs"
  | null states  = errorIn fun "empty states"
  | null outputs = errorIn fun "empty outputs"
  | otherwise    = ()

-- | Check if the two models have the same inputs, states, and outputs.
checkTwoModelsIn :: (Eq i, Eq s, Eq o) => String -> IOHMM i s o -> IOHMM i s o -> ()
checkTwoModelsIn fun model model'
  | is /= is' = errorIn fun "inputs disagree"
  | ss /= ss' = errorIn fun "states disagree"
  | os /= os' = errorIn fun "outputs disagree"
  | otherwise = ()
  where
    is  = inputs model
    is' = inputs model'
    ss  = states model
    ss' = states model'
    os  = outputs model
    os' = outputs model'

-- | Check if all the elements of the given inputs (outputs) are contained
--   in the 'inputs' ('outputs') of the model.
checkDataIn :: (Eq i, Eq o) => String -> IOHMM i s o -> [i] -> [o] -> ()
checkDataIn fun IOHMM {..} xs ys
  | any (`notElem` inputs) xs  = errorIn fun "illegal input data"
  | any (`notElem` outputs) ys = errorIn fun "illegal output data"
  | any (`notElem` xs) inputs  = errorIn fun "insufficient input data"
  | otherwise                  = ()

-- | Convert internal 'IOHMM' to 'IOHMM'.
fromInternal :: (Eq i, Eq s, Eq o) => [i] -> [s] -> [o] -> I.IOHMM -> IOHMM i s o
fromInternal is ss os I.IOHMM {..} = IOHMM { inputs           = is
                                           , states           = ss
                                           , outputs          = os
                                           , initialStateDist = C.fromList pi0'
                                           , transitionDist   = \i s -> case (elemIndex i is, elemIndex s ss) of
                                                                          (Nothing, _)     -> C.fromList []
                                                                          (_, Nothing)     -> C.fromList []
                                                                          (Just j, Just k) -> C.fromList $ w' j k
                                           , emissionDist     = \s -> case elemIndex s ss of
                                                                        Nothing -> C.fromList []
                                                                        Just i  -> C.fromList $ phi' i
                                           }
  where
    pi0'   = zip (H.toList initialStateDist) ss
    w' j k = zip (H.toList $ V.unsafeIndex transitionDist j H.! k) ss
    phi' i = zip (H.toList $ H.tr emissionDistT H.! i) os

-- | Convert 'IOHMM' to internal 'IOHMM'. The 'initialStateDist'',
--   'transitionDist'', and 'emissionDistT'' are normalized.
toInternal :: (Eq i, Eq s, Eq o) => IOHMM i s o -> I.IOHMM
toInternal IOHMM {..} = I.IOHMM { I.nInputs          = length inputs
                                , I.nStates          = length states
                                , I.nOutputs         = length outputs
                                , I.initialStateDist = pi0
                                , I.transitionDist   = w
                                , I.emissionDistT    = phi'
                                }
  where
    pi0_ = C.normalizeCategoricalPs initialStateDist
    w_ i = C.normalizeCategoricalPs . transitionDist i
    phi_ = C.normalizeCategoricalPs . emissionDist
    pi0  = H.fromList [pmf pi0_ s | s <- states]
    w    = V.fromList $ map (\i -> H.fromLists [[pmf (w_ i s) s' | s' <- states] | s <- states]) inputs
    phi' = H.fromLists [[pmf (phi_ s) o | s <- states] | o <- outputs]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.IOHMM." ++ fun ++ ": " ++ msg
