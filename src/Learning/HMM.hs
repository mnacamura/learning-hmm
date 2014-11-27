module Learning.HMM (
    HMM (..)
  , LogLikelihood
  , new
  , init
  , withEmission
  , viterbi
  , baumWelch
  , simulate
  ) where

import Prelude hiding (init)
import Control.Applicative ((<$>))
import Control.Arrow (first)
import Data.List (elemIndex, genericLength)
import Data.Maybe (fromJust)
import Data.Random.Distribution (pdf, rvar)
import Data.Random.Distribution.Categorical (Categorical)
import qualified Data.Random.Distribution.Categorical as C (
    fromList, fromWeightedList, normalizeCategoricalPs
  )
import Data.Random.Distribution.Categorical.Util ()
import Data.Random.RVar (RVar)
import Data.Random.Sample (sample)
import Data.Vector ((!))
import qualified Data.Vector as V (elemIndex, fromList, map, mapM, toList)
import qualified Data.Vector.Generic.Util.LinearAlgebra as V (transpose)
import Learning.HMM.Internal

-- | Parameter set of the hidden Markov model. Direct use of the
--   constructor is not recommended. Instead, call 'new' or 'init'.
data HMM s o = HMM { states  :: [s] -- ^ Hidden states
                   , outputs :: [o] -- ^ Outputs
                   , initialStateDist :: Categorical Double s
                     -- ^ Categorical distribution of initial states
                   , transitionDist :: s -> Categorical Double s
                     -- ^ Categorical distribution of next states
                     --   conditioned by the previous states
                   , emissionDist :: s -> Categorical Double o
                     -- ^ Categorical distribution of outputs conditioned
                     --   by the hidden states
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

-- | @new states outputs@ returns a model from the @states@ and @outputs@.
--   The 'initialStateDist' and 'emissionDist' are set to be uniform
--   distributions. The 'transitionDist' is specified as follows: with
--   probability 1/2, move to the same state, otherwise, move to a random
--   state (which might be the same state).
--
--   >>> new [1, 2 :: Int] ['C', 'D']
--   HMM {states = [1,2], outputs = "CD", initialStateDist = fromList [(0.5,1),(0.5,2)], transitionDist = [(fromList [(0.75,1),(0.25,2)],1),(fromList [(0.25,1),(0.75,2)],2)], emissionDist = [(fromList [(0.5,'C'),(0.5,'D')],1),(fromList [(0.5,'C'),(0.5,'D')],2)]}
new :: (Ord s, Ord o) => [s] -> [o] -> HMM s o
new ss os = HMM { states           = ss
                , outputs          = os
                , initialStateDist = pi0
                , transitionDist   = w
                , emissionDist     = phi
                }
  where
    pi0 = C.fromWeightedList [(1, s) | s <- ss]
    w s | s `elem` ss = C.fromList [(p s', s') | s' <- ss]
        | otherwise   = C.fromList []
      where
        k = genericLength ss
        p s' | s' == s   = 1/2 * (1 + 1/k)
             | otherwise = 1/2 / k
    phi s | s `elem` ss = C.fromWeightedList [(1, o) | o <- os]
          | otherwise   = C.fromList []

-- | @init states outputs@ returns a random variable of the model with
--   @states@ and @outputs@, wherein parameters are sampled from uniform
--   distributions.
init :: (Eq s, Eq o) => [s] -> [o] -> RVar (HMM s o)
init ss os = fromHMM' ss os <$> init' (length ss) (length os)

-- | @model \`withEmission\` xs@ returns a model in which the
--   'emissionDist' is updated by using the observed outputs @xs@. The
--   'emissionDist' is set to be normalized histograms each of which is
--   calculated from a partial set of @xs@ for each state. The partition is
--   based on the most likely state path obtained by the Viterbi algorithm.
withEmission :: (Eq s, Eq o) => HMM s o -> [o] -> HMM s o
withEmission model xs = fromHMM' ss os $ withEmission' model' xs'
  where
    ss     = states model
    os     = outputs model
    model' = toHMM' model
    xs'    = V.fromList $ fromJust $ mapM (`elemIndex` os) xs

-- | @viterbi model xs@ performs the Viterbi algorithm using the observed
--   outputs @xs@, and returns the most likely state path and its log
--   likelihood.
viterbi :: (Eq s, Eq o) => HMM s o -> [o] -> ([s], LogLikelihood)
viterbi model xs =
  checkModelIn "viterbi" model `seq`
  checkDataIn "viterbi" model xs `seq`
  first (V.toList . V.map (ss !)) $ viterbi' model' xs'
  where
    ss     = V.fromList $ states model
    os     = V.fromList $ outputs model
    model' = toHMM' model
    xs'    = fromJust $ V.mapM (`V.elemIndex` os) $ V.fromList xs

-- | @baumWelch model xs@ iteratively performs the Baum-Welch algorithm
--   using the observed outputs @xs@, and returns a list of updated models
--   and their corresponding log likelihoods.
baumWelch :: (Eq s, Eq o) => HMM s o -> [o] -> [(HMM s o, LogLikelihood)]
baumWelch model xs =
  checkModelIn "baumWelch" model `seq`
  checkDataIn "baumWelch" model xs `seq`
  map (first $ fromHMM' ss os) $ baumWelch' model' xs'
  where
    ss     = states model
    os     = outputs model
    model' = toHMM' model
    xs'    = V.fromList $ fromJust $ mapM (`elemIndex` os) xs

-- | @simulate model t@ generates a Markov process of length @t@ using the
--   @model@, and returns its state path and observed outputs.
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

-- | Convert 'HMM'' to 'HMM'.
fromHMM' :: (Eq s, Eq o) => [s] -> [o] -> HMM' -> HMM s o
fromHMM' ss os hmm' = HMM { states           = ss
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
    pi0 = initialStateDist' hmm'
    w   = transitionDist' hmm'
    phi = V.transpose $ emissionDistT' hmm'
    pi0'   = zip (V.toList pi0) ss
    w' i   = zip (V.toList $ w ! i) ss
    phi' i = zip (V.toList $ phi ! i) os

-- | Convert 'HMM' to 'HMM''. The 'initialStateDist'', 'transitionDist'',
--   and 'emissionDistT'' are normalized.
toHMM' :: (Eq s, Eq o) => HMM s o -> HMM'
toHMM' hmm = HMM' { nStates'          = length ss
                  , nOutputs'         = length os
                  , initialStateDist' = V.fromList pi0'
                  , transitionDist'   = V.fromList w'
                  , emissionDistT'    = V.fromList phi'
                  }
  where
    ss  = states hmm
    os  = outputs hmm
    pi0 = C.normalizeCategoricalPs $ initialStateDist hmm
    w   = C.normalizeCategoricalPs . transitionDist hmm
    phi = C.normalizeCategoricalPs . emissionDist hmm
    pi0' = [pdf pi0 s | s <- ss]
    w'   = [V.fromList [pdf (w s) s' | s' <- ss] | s <- ss]
    phi' = [V.fromList [pdf (phi s) o | s <- ss] | o <- os]

errorIn :: String -> String -> a
errorIn fun msg = error $ "Learning.HMM." ++ fun ++ ": " ++ msg
