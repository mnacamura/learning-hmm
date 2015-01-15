module Learning.HMM (
    HMM (..)
  , LogLikelihood
  , init
  , withEmission
  , viterbi
  , baumWelch
  , baumWelch'
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
import qualified Data.Vector.Unboxed as U (fromList)
import qualified Numeric.LinearAlgebra.Data as H (
    (!), fromList, fromLists, toList
  )
import qualified Numeric.LinearAlgebra.HMatrix as H (tr)
import Learning.HMM.Internal (LogLikelihood)
import qualified Learning.HMM.Internal as I

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
                   , initialStateDist :: Categorical Double s
                     -- ^ Categorical distribution of initial state
                   , transitionDist :: s -> Categorical Double s
                     -- ^ Categorical distribution of next state
                     --   conditioned by the previous state
                   , emissionDist :: s -> Categorical Double o
                     -- ^ Categorical distribution of output conditioned
                     --   by the hidden state
                   }

instance (Show s, Show o) => Show (HMM s o) where
  show = showHMM

showHMM :: (Show s, Show o) => HMM s o -> String
showHMM hmm = "HMM {states = "           ++ show ss
              ++ ", outputs = "          ++ show os
              ++ ", initialStateDist = " ++ show pi0
              ++ ", transitionDist = "   ++ show [(w s, s) | s <- ss]
              ++ ", emissionDist = "     ++ show [(phi s, s) | s <- ss]
              ++ "}"
  where
    ss  = states hmm
    os  = outputs hmm
    pi0 = initialStateDist hmm
    w   = transitionDist hmm
    phi = emissionDist hmm

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
withEmission model xs = fromInternal ss os $ I.withEmission model' xs'
  where
    ss     = states model
    os     = outputs model
    os'    = V.fromList os
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') xs

-- | @viterbi model xs@ performs the Viterbi algorithm using the observed
--   outputs @xs@, and returns the most likely state path and its log
--   likelihood.
viterbi :: (Eq s, Eq o) => HMM s o -> [o] -> ([s], LogLikelihood)
viterbi model xs =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs `seq`
  first toStates $ I.viterbi model' xs'
  where
    ss'    = V.fromList $ states model
    os'    = V.fromList $ outputs model
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') xs
    toStates = V.toList . V.map (V.unsafeIndex ss') . G.convert

-- | @baumWelch model xs@ iteratively performs the Baum-Welch algorithm
--   using the observed outputs @xs@, and returns a list of updated models
--   and their corresponding log likelihoods.
baumWelch :: (Eq s, Eq o) => HMM s o -> [o] -> [(HMM s o, LogLikelihood)]
baumWelch model xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  map (first $ fromInternal ss os) $ I.baumWelch model' xs'
  where
    ss     = states model
    os     = outputs model
    os'    = V.fromList os
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') xs

-- | @baumWelch' model xs@ performs the Baum-Welch algorithm using the
--   observed outputs @xs@, and returns a model locally maximizing its log
--   likelihood.
baumWelch' :: (Eq s, Eq o) => HMM s o -> [o] -> (HMM s o, LogLikelihood)
baumWelch' model xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  first (fromInternal ss os) $ I.baumWelch' model' xs'
  where
    ss     = states model
    os     = outputs model
    os'    = V.fromList os
    model' = toInternal model
    xs'    = U.fromList $ fromJust $ mapM (`V.elemIndex` os') xs

-- | @simulate model t@ generates a Markov process of length @t@ using the
--   @model@, and returns its state path and outputs.
simulate :: HMM s o -> Int -> RVar ([s], [o])
simulate model step
  | step < 1  = return ([], [])
  | otherwise = do s0 <- sample $ rvar pi0
                   x0 <- sample $ rvar $ phi s0
                   unzip . ((s0, x0) :) <$> sim s0 (step - 1)
  where
    sim _ 0 = return []
    sim s t = do s' <- sample $ rvar $ w s
                 x' <- sample $ rvar $ phi s'
                 ((s', x') :) <$> sim s' (t - 1)
    pi0 = initialStateDist model
    w   = transitionDist model
    phi = emissionDist model

-- | Check if the model is valid in the sense of whether the 'states' and
--   'outputs' are not empty.
checkModelIn :: String -> HMM s o -> ()
checkModelIn fun hmm
  | null ss   = err "empty states"
  | null os   = err "empty outputs"
  | otherwise = ()
  where
    ss = states hmm
    os = outputs hmm
    err = errorIn fun

-- | Check if all the elements of the observed outputs are contained in the
--   'outputs' of the model.
checkDataIn :: Eq o => String -> HMM s o -> [o] -> ()
checkDataIn fun hmm xs
  | all (`elem` os) xs = ()
  | otherwise          = err "illegal data"
  where
    os = outputs hmm
    err = errorIn fun

-- | Convert internal 'HMM' to 'HMM'.
fromInternal :: (Eq s, Eq o) => [s] -> [o] -> I.HMM -> HMM s o
fromInternal ss os hmm' = HMM { states           = ss
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
    pi0 = I.initialStateDist hmm'
    w   = I.transitionDist hmm'
    phi = H.tr $ I.emissionDistT hmm'
    pi0'   = zip (H.toList pi0) ss
    w' i   = zip (H.toList $ w H.! i) ss
    phi' i = zip (H.toList $ phi H.! i) os

-- | Convert 'HMM' to internal 'HMM'. The 'initialStateDist'',
--   'transitionDist'', and 'emissionDistT'' are normalized.
toInternal :: (Eq s, Eq o) => HMM s o -> I.HMM
toInternal hmm = I.HMM { I.nStates          = length ss
                       , I.nOutputs         = length os
                       , I.initialStateDist = pi0
                       , I.transitionDist   = w
                       , I.emissionDistT    = phi'
                       }
  where
    ss   = states hmm
    os   = outputs hmm
    pi0_ = C.normalizeCategoricalPs $ initialStateDist hmm
    w_   = C.normalizeCategoricalPs . transitionDist hmm
    phi_ = C.normalizeCategoricalPs . emissionDist hmm
    pi0  = H.fromList [pdf pi0_ s | s <- ss]
    w    = H.fromLists [[pdf (w_ s) s' | s' <- ss] | s <- ss]
    phi' = H.fromLists [[pdf (phi_ s) o | s <- ss] | o <- os]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.HMM." ++ fun ++ ": " ++ msg
