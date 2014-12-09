module Learning.IOHMM (
    IOHMM (..)
  , LogLikelihood
  , init
  , withEmission
  , viterbi
  , baumWelch
  , simulate
  ) where

import Prelude hiding (init)
import Control.Applicative ((<$>))
import Control.Arrow (first)
import Data.List (elemIndex)
import Data.Maybe (fromJust)
import Data.Random.Distribution (pdf, rvar)
import Data.Random.Distribution.Categorical (Categorical)
import qualified Data.Random.Distribution.Categorical as C (
    fromList, normalizeCategoricalPs
  )
import Data.Random.Distribution.Categorical.Util ()
import Data.Random.RVar (RVar)
import Data.Random.Sample (sample)
import qualified Data.Vector as V (
    elemIndex, fromList, map, toList, unsafeIndex
  )
import qualified Data.Vector.Generic as G (convert)
import qualified Data.Vector.Unboxed as U (fromList, zip)
import qualified Numeric.LinearAlgebra.Data as H (
    (!), fromList, fromLists, toList
  )
import qualified Numeric.LinearAlgebra.HMatrix as H (tr)
import Learning.IOHMM.Internal (LogLikelihood)
import qualified Learning.IOHMM.Internal as I

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
                         , initialStateDist :: Categorical Double s
                           -- ^ Categorical distribution of initial state
                         , transitionDist :: i -> s -> Categorical Double s
                           -- ^ Categorical distribution of next state
                           --   conditioned by the input and previous state
                         , emissionDist :: s -> Categorical Double o
                           -- ^ Categorical distribution of output conditioned
                           --   by the hidden state
                         }

instance (Show i, Show s, Show o) => Show (IOHMM i s o) where
  show = showIOHMM

showIOHMM :: (Show i, Show s, Show o) => IOHMM i s o -> String
showIOHMM hmm = "IOHMM {inputs = "           ++ show is
                  ++ ", states = "           ++ show ss
                  ++ ", outputs = "          ++ show os
                  ++ ", initialStateDist = " ++ show pi0
                  ++ ", transitionDist = "   ++ show [(w i s, (i, s)) | i <- is, s <- ss]
                  ++ ", emissionDist = "     ++ show [(phi s, s) | s <- ss]
                  ++ "}"
  where
    is  = inputs hmm
    ss  = states hmm
    os  = outputs hmm
    pi0 = initialStateDist hmm
    w   = transitionDist hmm
    phi = emissionDist hmm

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
withEmission model xs ys = fromInternal is ss os $ I.withEmission model' $ U.zip xs' ys'
  where
    is     = inputs model
    is'    = V.fromList is
    ss     = states model
    os     = outputs model
    os'    = V.fromList os
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` is') xs
    ys'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') ys

-- | @viterbi model xs ys@ performs the Viterbi algorithm using the inputs
--   @xs@ and outputs @ys@, and returns the most likely state path and its
--   log likelihood.
--   If the lengths of @xs@ and @ys@ are different, the longer one is cut
--   by the length of the shorter one.
viterbi :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> ([s], LogLikelihood)
viterbi model xs ys =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs ys `seq`
  first toStates $ I.viterbi model' $ U.zip xs' ys'
  where
    is'    = V.fromList $ inputs model
    ss'    = V.fromList $ states model
    os'    = V.fromList $ outputs model
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` is') xs
    ys'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') ys
    toStates = V.toList . V.map (V.unsafeIndex ss') . G.convert

-- | @baumWelch model xs ys@ iteratively performs the Baum-Welch algorithm
--   using the inputs @xs@ and outputs @ys@, and returns a list of updated
--   models and their corresponding log likelihoods.
--   If the lengths of @xs@ and @ys@ are different, the longer one is cut
--   by the length of the shorter one.
baumWelch :: (Eq i, Eq s, Eq o) => IOHMM i s o -> [i] -> [o] -> [(IOHMM i s o, LogLikelihood)]
baumWelch model xs ys =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs ys `seq`
  map (first $ fromInternal is ss os) $ I.baumWelch model' $ U.zip xs' ys'
  where
    is     = inputs model
    is'    = V.fromList is
    ss     = states model
    os     = outputs model
    os'    = V.fromList os
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` is') xs
    ys'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') ys

-- | @simulate model xs@ generates a Markov process coinciding with the
--   inputs @xs@ using the @model@, and returns its state path and observed
--   outputs.
simulate :: IOHMM i s o -> [i] -> RVar ([s], [o])
simulate model xs
  | null xs   = return ([], [])
  | otherwise = do s0 <- sample $ rvar pi0
                   y0 <- sample $ rvar $ phi s0
                   unzip . ((s0, y0) :) <$> sim s0 (tail xs)
  where
    sim _ []     = return []
    sim s (x:xs') = do s' <- sample $ rvar $ w x s
                       y' <- sample $ rvar $ phi s'
                       ((s', y') :) <$> sim s' xs'
    pi0 = initialStateDist model
    w   = transitionDist model
    phi = emissionDist model

-- | Check if the model is valid in the sense of whether the 'states' and
--   'outputs' are not empty.
checkModelIn :: String -> IOHMM i s o -> ()
checkModelIn fun hmm
  | null is   = err "empty inputs"
  | null ss   = err "empty states"
  | null os   = err "empty outputs"
  | otherwise = ()
  where
    is = inputs hmm
    ss = states hmm
    os = outputs hmm
    err = errorIn fun

-- | Check if all the elements of the given inputs (outputs) are contained
--   in the 'inputs' ('outputs') of the model.
checkDataIn :: (Eq i, Eq o) => String -> IOHMM i s o -> [i] -> [o] -> ()
checkDataIn fun hmm xs ys
  | all (`elem` is) xs && all (`elem` os) ys = ()
  | otherwise                                = err "illegal data"
  where
    is = inputs hmm
    os = outputs hmm
    err = errorIn fun

-- | Convert internal 'IOHMM' to 'IOHMM'.
fromInternal :: (Eq i, Eq s, Eq o) => [i] -> [s] -> [o] -> I.IOHMM -> IOHMM i s o
fromInternal is ss os hmm' = IOHMM { inputs           = is
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
    pi0 = I.initialStateDist hmm'
    w   = I.transitionDist hmm'
    phi = H.tr $ I.emissionDistT hmm'
    pi0'   = zip (H.toList pi0) ss
    w' j k = zip (H.toList $ V.unsafeIndex w j H.! k) ss
    phi' i = zip (H.toList $ phi H.! i) os

-- | Convert 'IOHMM' to internal 'IOHMM'. The 'initialStateDist'',
--   'transitionDist'', and 'emissionDistT'' are normalized.
toInternal :: (Eq i, Eq s, Eq o) => IOHMM i s o -> I.IOHMM
toInternal hmm = I.IOHMM { I.nInputs          = length is
                         , I.nStates          = length ss
                         , I.nOutputs         = length os
                         , I.initialStateDist = pi0
                         , I.transitionDist   = w
                         , I.emissionDistT    = phi'
                         }
  where
    is   = inputs hmm
    ss   = states hmm
    os   = outputs hmm
    pi0_ = C.normalizeCategoricalPs $ initialStateDist hmm
    w_ i = C.normalizeCategoricalPs . (transitionDist hmm) i
    phi_ = C.normalizeCategoricalPs . emissionDist hmm
    pi0  = H.fromList [pdf pi0_ s | s <- ss]
    w    = V.fromList $ map (\i -> H.fromLists [[pdf (w_ i s) s' | s' <- ss] | s <- ss]) is
    phi' = H.fromLists [[pdf (phi_ s) o | s <- ss] | o <- os]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.IOHMM." ++ fun ++ ": " ++ msg
