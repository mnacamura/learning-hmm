module Learning.HMM.Internal (
    HMM' (..)
  , LogLikelihood
  , Probability
  , init'
  , withEmission'
  , viterbi'
  , baumWelch'
  -- , baumWelch1'
  -- , forward'
  -- , backward'
  ) where

import Control.Applicative ((<$>))
import Control.Monad (forM_, replicateM)
import Control.Monad.ST (runST)
import qualified Data.Map.Strict as M (findWithDefault)
import Data.Random.RVar (RVar)
import Data.Random.Distribution.Simplex (stdSimplex)
import Data.Vector (Vector, (!))
import qualified Data.Vector as V (
    filter, foldl1', freeze, fromList, generate, length, map, maximum
  , maxIndex, replicate, sum, tail, zip, zipWith, zipWith3
  )
import qualified Data.Vector.Mutable as MV (new, read, write)
import qualified Data.Vector.Util as V (frequencies)
import Data.Vector.Util.LinearAlgebra (
    (>+>), (>.>), (>/>), (#+#), (.>), (>/), (#.>), (<.#)
  )
import qualified Data.Vector.Util.LinearAlgebra as V (transpose)

type LogLikelihood = Double
type Probability   = Double

-- | More efficient data structure of the 'HMM' model. The 'states' and
--   'outputs' in 'HMM' are represented by their indices. The
--   'initialStateDist', 'transitionDist', and 'emissionDist' are
--   represented by matrices. The 'emissionDistT'' is a transposed matrix
--   in order to simplify the calculation.
data HMM' = HMM' { nStates'          :: Int -- ^ Number of states
                 , nOutputs'         :: Int -- ^ Number of outputs
                 , initialStateDist' :: Vector Probability
                 , transitionDist'   :: Vector (Vector Probability)
                 , emissionDistT'    :: Vector (Vector Probability)
                 }

init' :: Int -> Int -> RVar HMM'
init' n m = do
  pi0 <- V.fromList <$> stdSimplex (n-1)
  w   <- V.fromList <$> replicateM n (V.fromList <$> stdSimplex (n-1))
  phi <- V.fromList <$> replicateM n (V.fromList <$> stdSimplex (m-1))
  return HMM' { nStates'          = n
              , nOutputs'         = m
              , initialStateDist' = pi0
              , transitionDist'   = w
              , emissionDistT'    = V.transpose phi
              }

withEmission' :: HMM' -> Vector Int -> HMM'
withEmission' model xs = model { emissionDistT' = phi' }
  where
    ss   = V.generate (nStates' model) id
    os   = V.generate (nOutputs' model) id
    phi' = let (path, _) = viterbi' model xs
               freqs     = V.frequencies $ V.zip path xs
               hists     = V.map (\s -> V.map (\o ->
                                 M.findWithDefault 0 (s, o) freqs) os) ss
           in V.transpose $ V.map (\f -> f >/ V.sum f) hists

viterbi' :: HMM' -> Vector Int -> (Vector Int, LogLikelihood)
viterbi' model xs = (path, logL)
  where
    n = V.length xs

    -- First, we calculate the value function and the state maximizing it
    -- for each time.
    deltas :: Vector (Vector LogLikelihood)
    psis   :: Vector (Vector Int)
    (deltas, psis) = runST $ do
      ds <- MV.new n
      ps <- MV.new n
      ds `MV.write` 0 $ V.map log (phi' ! (xs ! 0)) >+> V.map log pi0
      ps `MV.write` 0 $ V.replicate k 0
      forM_ [1..(n-1)] $ \t -> do
        d <- ds `MV.read` (t-1)
        let dws = V.map (\wj -> d >+> V.map log wj) w'
        ds `MV.write` t $ V.map log (phi' ! (xs ! t)) >+> V.map V.maximum dws
        ps `MV.write` t $ V.map V.maxIndex dws
      ds' <- V.freeze ds
      ps' <- V.freeze ps
      return (ds', ps')
      where
        k    = nStates' model
        pi0  = initialStateDist' model
        w'   = V.transpose $ transitionDist' model
        phi' = emissionDistT' model

    -- The most likely path and corresponding log likelihood is as follows.
    path = runST $ do
      ix <- MV.new n
      ix `MV.write` (n-1) $ V.maxIndex $ deltas ! (n-1)
      forM_ (reverse [0..(n-2)]) $ \t -> do
        i <- ix `MV.read` (t+1)
        ix `MV.write` t $ psis ! (t+1) ! i
      V.freeze ix
    logL = V.maximum $ deltas ! (n-1)

baumWelch' :: HMM' -> Vector Int -> [(HMM', LogLikelihood)]
baumWelch' model xs = zip models (tail logLs)
  where
    n = V.length xs
    step (m, _)     = m `seq` baumWelch1' m n xs
    (models, logLs) = unzip $ iterate step (model, undefined)

-- | Perform one step of the Baum-Welch algorithm and return the updated
--   model and the likelihood of the old model.
baumWelch1' :: HMM' -> Int -> Vector Int -> (HMM', LogLikelihood)
baumWelch1' model n xs = (model', logL)
  where
    -- First, we calculate the alpha, beta, and scaling values using the
    -- forward-backward algorithm.
    (alphas, cs) = forward' model n xs
    betas        = backward' model n xs cs

    -- Based on the alpha, beta, and scaling values, we calculate the
    -- gamma and xi values.
    gammas = V.zipWith3 (\a b c -> a >.> b >/ c) alphas betas cs
    xis    = V.zipWith3 (\a b x -> let w1 = V.zipWith (.>) a w0
                                   in V.map ((phi0 ! x) >.> b >.>) w1)
               alphas (V.tail betas) (V.tail xs)
      where
        w0   = transitionDist' model
        phi0 = emissionDistT' model

    -- Using the gamma and xi values, we obtain the optimal initial state
    -- probability vector, transition probability matrix, and emission
    -- probability matrix.
    pi0  = gammas ! 0
    w    = let ds = V.foldl1' (#+#) xis -- denominators
               ns = V.map V.sum ds      -- numerators
           in V.zipWith (>/) ds ns
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip xs gammas
               ds    = V.foldl1' (>+>) . gs'  -- denominators
               ns    = V.foldl1' (>+>) gammas -- numerators
           in V.map (\o -> ds o >/> ns) os
      where
        os = V.generate (nOutputs' model) id

    -- We finally obtain the new model and the likelihood for the old model.
    model' = model { initialStateDist' = pi0
                   , transitionDist'   = w
                   , emissionDistT'    = phi'
                   }
    logL = - (V.sum $ V.map log cs)

-- | Return alphas and scaling variables.
forward' :: HMM' -> Int -> Vector Int -> (Vector (Vector Probability), Vector Probability)
{-# INLINE forward' #-}
forward' model n xs = runST $ do
  as <- MV.new n
  cs <- MV.new n
  let a0 = (phi' ! (xs ! 0)) >.> pi0
      c0 = 1 / V.sum a0
  MV.write as 0 (c0 .> a0)
  MV.write cs 0 c0
  forM_ [1..(n-1)] $ \t -> do
    a <- as `MV.read` (t-1)
    let a' = (phi' ! (xs ! t)) >.> (a <.# w)
        c' = 1 / V.sum a'
    MV.write as t (c' .> a')
    MV.write cs t c'
  as' <- V.freeze as
  cs' <- V.freeze cs
  return (as', cs')
  where
    pi0  = initialStateDist' model
    w    = transitionDist' model
    phi' = emissionDistT' model

-- | Return betas using scaling variables.
backward' :: HMM' -> Int -> Vector Int -> Vector Probability -> Vector (Vector Probability)
{-# INLINE backward' #-}
backward' model n xs cs = runST $ do
  bs <- MV.new n
  let bE = V.replicate k 1
      cE = cs ! (n-1)
  MV.write bs (n-1) $ cE .> bE
  forM_ (reverse [0..(n-2)]) $ \t -> do
    b <- MV.read bs (t+1)
    let b' = w #.> ((phi' ! (xs ! (t+1))) >.> b)
        c' = cs ! t
    MV.write bs t $ c' .> b'
  V.freeze bs
  where
    k    = nStates' model
    w    = transitionDist' model
    phi' = emissionDistT' model
