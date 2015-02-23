{-# LANGUAGE RecordWildCards #-}

module Learning.HMM
  ( HMM (..)
  , LogLikelihood
  , init
  , withEmission
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
import qualified Data.Vector.Unboxed                  as U   ( fromList )
import           Learning.HMM.Internal                       ( LogLikelihood )
import qualified Learning.HMM.Internal                as I
import qualified Numeric.LinearAlgebra.Data           as H   ( (!), fromList, fromLists, toList )
import qualified Numeric.LinearAlgebra.HMatrix        as H   ( tr )
import           Prelude                              hiding ( init )

-- | Parameter set of the hidden Markov model with discrete emission.
--   The model schema is as follows.
--
--   @
--       z_0 -> z_1 -> ... -> z_n
--        |      |             |
--        v      v             v
--       x_0    x_1           x_n
--   @
--
--   Here, @[z_0, z_1, ..., z_n]@ are hidden states and @[x_0, x_1, ..., x_n]@
--   are observed outputs. @z_0@ is determined by the 'initialStateDist'.
--   For @i = 1, ..., n@, @z_i@ is determined by the 'transitionDist'
--   conditioned by @z_{i-1}@.
--   For @i = 0, ..., n@, @x_i@ is determined by the 'emissionDist'
--   conditioned by @z_i@.
data HMM s o = HMM { states  :: [s]
                   , outputs :: [o]
                   , initialStateDist :: C.Categorical Double s
                     -- ^ Categorical distribution of initial state
                   , transitionDist :: s -> C.Categorical Double s
                     -- ^ Categorical distribution of next state
                     --   conditioned by the previous state
                   , emissionDist :: s -> C.Categorical Double o
                     -- ^ Categorical distribution of output conditioned
                     --   by the hidden state
                   }

instance (Show s, Show o) => Show (HMM s o) where
  show HMM {..} = "HMM {states = "           ++ show states
                  ++ ", outputs = "          ++ show outputs
                  ++ ", initialStateDist = " ++ show initialStateDist
                  ++ ", transitionDist = "   ++ show [(transitionDist s, s) | s <- states]
                  ++ ", emissionDist = "     ++ show [(emissionDist s, s) | s <- states]
                  ++ "}"

-- | @init states outputs@ returns a random variable of models with the
--   @states@ and @outputs@, wherein parameters are sampled from uniform
--   distributions.
init :: (Eq s, Eq o) => [s] -> [o] -> RVar (HMM s o)
init ss os = fromInternal ss os <$> I.init (length ss) (length os)

-- | @model \`withEmission\` xs@ returns a model in which the
--   'emissionDist' is updated by re-estimations using the observed outputs
--   @xs@. The 'emissionDist' is set to be normalized histograms each of
--   which is calculated from segumentations of @xs@ based on the Viterbi
--   state path.
withEmission :: (Eq s, Eq o) => HMM s o -> [o] -> HMM s o
withEmission (model @ HMM {..}) xs = fromInternal states outputs $ I.withEmission model' xs'
  where
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') xs

-- | @viterbi model xs@ performs the Viterbi algorithm using the observed
--   outputs @xs@, and returns the most likely state path and its log
--   likelihood.
viterbi :: (Eq s, Eq o) => HMM s o -> [o] -> ([s], LogLikelihood)
viterbi (model @ HMM {..}) xs =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs `seq`
  first toStates $ I.viterbi model' xs'
  where
    states'  = V.fromList states
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') xs
    toStates = V.toList . V.map (V.unsafeIndex states') . G.convert

-- | @baumWelch model xs@ iteratively performs the Baum-Welch algorithm
--   using the observed outputs @xs@, and returns a list of updated models
--   and their corresponding log likelihoods.
baumWelch :: (Eq s, Eq o) => HMM s o -> [o] -> [(HMM s o, LogLikelihood)]
baumWelch (model @ HMM {..}) xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  map (first $ fromInternal states outputs) $ I.baumWelch model' xs'
  where
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') xs

-- | @baumWelch' model xs@ performs the Baum-Welch algorithm using the
--   observed outputs @xs@, and returns a model locally maximizing its log
--   likelihood.
baumWelch' :: (Eq s, Eq o) => HMM s o -> [o] -> (HMM s o, LogLikelihood)
baumWelch' (model @ HMM {..}) xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  first (fromInternal states outputs) $ I.baumWelch' model' xs'
  where
    outputs' = V.fromList outputs
    model'   = toInternal model
    xs'      = U.fromList $ fromJust $ mapM (`V.elemIndex` outputs') xs

-- | @simulate model t@ generates a Markov process of length @t@ using the
--   @model@, and returns its state path and outputs.
simulate :: HMM s o -> Int -> RVar ([s], [o])
simulate HMM {..} step
  | step < 1  = return ([], [])
  | otherwise = do s0 <- rvar initialStateDist
                   x0 <- rvar $ emissionDist s0
                   unzip . ((s0, x0) :) <$> sim s0 (step - 1)
  where
    sim _ 0 = return []
    sim s t = do s' <- rvar $ transitionDist s
                 x' <- rvar $ emissionDist s'
                 ((s', x') :) <$> sim s' (t - 1)

-- | Check if the model is valid in the sense of whether the 'states' and
--   'outputs' are not empty.
checkModelIn :: String -> HMM s o -> ()
checkModelIn fun HMM {..}
  | null states  = errorIn fun "empty states"
  | null outputs = errorIn fun "empty outputs"
  | otherwise    = ()

-- | Check if all the elements of the observed outputs are contained in the
--   'outputs' of the model.
checkDataIn :: Eq o => String -> HMM s o -> [o] -> ()
checkDataIn fun HMM {..} xs
  | all (`elem` outputs) xs = ()
  | otherwise               = errorIn fun "illegal data"

-- | Convert internal 'HMM' to 'HMM'.
fromInternal :: (Eq s, Eq o) => [s] -> [o] -> I.HMM -> HMM s o
fromInternal ss os I.HMM {..} = HMM { states           = ss
                                    , outputs          = os
                                    , initialStateDist = C.fromList pi0'
                                    , transitionDist   = \s -> case elemIndex s ss of
                                                                 Nothing -> C.fromList []
                                                                 Just i  -> C.fromList $ w' i
                                    , emissionDist     = \s -> case elemIndex s ss of
                                                                 Nothing -> C.fromList []
                                                                 Just i  -> C.fromList $ phi' i
                                    }
  where
    pi0'   = zip (H.toList initialStateDist) ss
    w' i   = zip (H.toList $ transitionDist H.! i) ss
    phi' i = zip (H.toList $ H.tr emissionDistT H.! i) os

-- | Convert 'HMM' to internal 'HMM'. The 'initialStateDist'',
--   'transitionDist'', and 'emissionDistT'' are normalized.
toInternal :: (Eq s, Eq o) => HMM s o -> I.HMM
toInternal HMM {..} = I.HMM { I.nStates          = length states
                            , I.nOutputs         = length outputs
                            , I.initialStateDist = pi0
                            , I.transitionDist   = w
                            , I.emissionDistT    = phi'
                            }
  where
    pi0_ = C.normalizeCategoricalPs initialStateDist
    w_   = C.normalizeCategoricalPs . transitionDist
    phi_ = C.normalizeCategoricalPs . emissionDist
    pi0  = H.fromList [pmf pi0_ s | s <- states]
    w    = H.fromLists [[pmf (w_ s) s' | s' <- states] | s <- states]
    phi' = H.fromLists [[pmf (phi_ s) o | s <- states] | o <- outputs]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.HMM." ++ fun ++ ": " ++ msg
