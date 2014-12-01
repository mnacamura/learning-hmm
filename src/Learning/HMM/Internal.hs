module Learning.HMM.Internal (
    HMM' (..)
  , LogLikelihood
  , init'
  , withEmission'
  , viterbi'
  , baumWelch'
  -- , baumWelch1'
  -- , forward'
  -- , backward'
  ) where

import Control.Applicative ((<$>))
import Control.DeepSeq (NFData, force, rnf)
import Control.Monad (forM_, replicateM)
import Control.Monad.ST (runST)
import qualified Data.Map.Strict as M (findWithDefault)
import Data.Random.RVar (RVar)
import Data.Random.Distribution.Simplex (stdSimplex)
import qualified Data.Vector as V (
    Vector, filter, foldl1', fromList, generate, map
  , unsafeFreeze, unsafeIndex, unsafeTail, zip, zipWith, zipWith3
  )
import qualified Data.Vector.Generic as G (convert)
import qualified Data.Vector.Generic.Util as G (frequencies)
import Data.Vector.Generic.Util.LinearAlgebra (
    (>+>), (>.>), (>/>), (#+#), (.>), (>/), (#.>), (<.#)
  )
import qualified Data.Vector.Generic.Util.LinearAlgebra as G (transpose)
import qualified Data.Vector.Mutable as MV (
    unsafeNew, unsafeRead, unsafeWrite
  )
import qualified Data.Vector.Unboxed as U (
    Vector, fromList, generate, length, map, maxIndex, maximum
  , replicate, sum, unsafeFreeze, unsafeIndex, unsafeTail, zip
  )
import qualified Data.Vector.Unboxed.Mutable as MU (
    unsafeNew, unsafeRead, unsafeWrite
  )

type LogLikelihood = Double

-- | More efficient data structure of the 'HMM' model. The 'states' and
--   'outputs' in 'HMM' are represented by their indices. The
--   'initialStateDist', 'transitionDist', and 'emissionDist' are
--   represented by matrices. The 'emissionDistT'' is a transposed matrix
--   in order to simplify the calculation.
data HMM' = HMM' { nStates'          :: Int -- ^ Number of states
                 , nOutputs'         :: Int -- ^ Number of outputs
                 , initialStateDist' :: U.Vector Double
                 , transitionDist'   :: V.Vector (U.Vector Double)
                 , emissionDistT'    :: V.Vector (U.Vector Double)
                 }

instance NFData HMM' where
  rnf hmm' = rnf n `seq` rnf m `seq` rnf pi0 `seq` rnf w `seq` rnf phi'
    where
      n    = nStates' hmm'
      m    = nOutputs' hmm'
      pi0  = initialStateDist' hmm'
      w    = transitionDist' hmm'
      phi' = emissionDistT' hmm'

init' :: Int -> Int -> RVar HMM'
init' n m = do
  pi0 <- U.fromList <$> stdSimplex (n-1)
  w   <- V.fromList <$> replicateM n (U.fromList <$> stdSimplex (n-1))
  phi <- V.fromList <$> replicateM n (U.fromList <$> stdSimplex (m-1))
  return HMM' { nStates'          = n
              , nOutputs'         = m
              , initialStateDist' = pi0
              , transitionDist'   = w
              , emissionDistT'    = G.transpose phi
              }

withEmission' :: HMM' -> U.Vector Int -> HMM'
withEmission' model xs = model { emissionDistT' = G.transpose phi }
  where
    ss  = V.generate (nStates' model) id
    os  = U.generate (nOutputs' model) id
    phi = let (path, _) = viterbi' model xs
              freqs     = G.frequencies $ U.zip path xs
              hists     = V.map (\s -> U.map (\o ->
                                M.findWithDefault 0 (s, o) freqs) os) ss
          in V.map (\f -> f >/ U.sum f) hists

viterbi' :: HMM' -> U.Vector Int -> (U.Vector Int, LogLikelihood)
viterbi' model xs = (path, logL)
  where
    n = U.length xs

    -- First, we calculate the value function and the state maximizing it
    -- for each time.
    deltas :: V.Vector (U.Vector Double)
    psis   :: V.Vector (U.Vector Int)
    (deltas, psis) = runST $ do
      ds <- MV.unsafeNew n
      ps <- MV.unsafeNew n
      MV.unsafeWrite ds 0 $ U.map log (V.unsafeIndex phi' (U.unsafeIndex xs 0)) >+> U.map log pi0
      MV.unsafeWrite ps 0 $ U.replicate k 0
      forM_ [1..(n-1)] $ \t -> do
        d <- MV.unsafeRead ds (t-1)
        let dws = V.map (\wj -> d >+> U.map log wj) w'
        MV.unsafeWrite ds t $ U.map log (V.unsafeIndex phi' (U.unsafeIndex xs t)) >+> G.convert (V.map U.maximum dws)
        MV.unsafeWrite ps t $ G.convert (V.map U.maxIndex dws)
      ds' <- V.unsafeFreeze ds
      ps' <- V.unsafeFreeze ps
      return (ds', ps')
      where
        k    = nStates' model
        pi0  = initialStateDist' model
        w'   = G.transpose $ transitionDist' model
        phi' = emissionDistT' model

    -- The most likely path and corresponding log likelihood are as follows.
    path = runST $ do
      ix <- MU.unsafeNew n
      MU.unsafeWrite ix (n-1) $ U.maxIndex (V.unsafeIndex deltas (n-1))
      forM_ (reverse [0..(n-2)]) $ \t -> do
        i <- MU.unsafeRead ix (t+1)
        MU.unsafeWrite ix t $ U.unsafeIndex (V.unsafeIndex psis (t+1)) i
      U.unsafeFreeze ix
    logL = U.maximum $ V.unsafeIndex deltas (n-1)

baumWelch' :: HMM' -> U.Vector Int -> [(HMM', LogLikelihood)]
baumWelch' model xs = zip models (tail logLs)
  where
    n = U.length xs
    step (m, _)     = baumWelch1' m n xs
    (models, logLs) = unzip $ iterate step (model, undefined)

-- | Perform one step of the Baum-Welch algorithm and return the updated
--   model and the likelihood of the old model.
baumWelch1' :: HMM' -> Int -> U.Vector Int -> (HMM', LogLikelihood)
baumWelch1' model n xs = force (model', logL)
  where
    -- First, we calculate the alpha, beta, and scaling values using the
    -- forward-backward algorithm.
    (alphas, cs) = forward' model n xs
    betas        = backward' model n xs cs

    -- Based on the alpha, beta, and scaling values, we calculate the
    -- posterior distribution, i.e., gamma and xi values.
    (gammas, xis) = posterior' model n xs alphas betas cs

    -- Using the gamma and xi values, we obtain the optimal initial state
    -- probability vector, transition probability matrix, and emission
    -- probability matrix.
    pi0  = V.unsafeIndex gammas 0
    w    = let ds = V.foldl1' (#+#) xis -- denominators
               ns = V.map U.sum ds      -- numerators
           in V.zipWith (>/) ds ns
    phi' = let gs' o = V.map snd $ V.filter ((== o) . fst) $ V.zip (G.convert xs) gammas
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
    logL = - (U.sum $ U.map log cs)

-- | Return alphas and scaling variables.
forward' :: HMM' -> Int -> U.Vector Int -> (V.Vector (U.Vector Double), U.Vector Double)
{-# INLINE forward' #-}
forward' model n xs = runST $ do
  as <- MV.unsafeNew n
  cs <- MU.unsafeNew n
  let a0 = V.unsafeIndex phi' (U.unsafeIndex xs 0) >.> pi0
      c0 = 1 / U.sum a0
  MV.unsafeWrite as 0 (c0 .> a0)
  MU.unsafeWrite cs 0 c0
  forM_ [1..(n-1)] $ \t -> do
    a <- MV.unsafeRead as (t-1)
    let a' = V.unsafeIndex phi' (U.unsafeIndex xs t) >.> (a <.# w)
        c' = 1 / U.sum a'
    MV.unsafeWrite as t (c' .> a')
    MU.unsafeWrite cs t c'
  as' <- V.unsafeFreeze as
  cs' <- U.unsafeFreeze cs
  return (as', cs')
  where
    pi0  = initialStateDist' model
    w    = transitionDist' model
    phi' = emissionDistT' model

-- | Return betas using scaling variables.
backward' :: HMM' -> Int -> U.Vector Int -> U.Vector Double -> V.Vector (U.Vector Double)
{-# INLINE backward' #-}
backward' model n xs cs = runST $ do
  bs <- MV.unsafeNew n
  let bE = U.replicate k 1
      cE = U.unsafeIndex cs (n-1)
  MV.unsafeWrite bs (n-1) $ cE .> bE
  forM_ (reverse [0..(n-2)]) $ \t -> do
    b <- MV.unsafeRead bs (t+1)
    let b' = w #.> (V.unsafeIndex phi' (U.unsafeIndex xs (t+1)) >.> b)
        c' = U.unsafeIndex cs t
    MV.unsafeWrite bs t $ c' .> b'
  V.unsafeFreeze bs
  where
    k    = nStates' model
    w    = transitionDist' model
    phi' = emissionDistT' model

-- | Return the posterior distribution.
posterior' :: HMM' -> Int -> U.Vector Int -> V.Vector (U.Vector Double) -> V.Vector (U.Vector Double) -> U.Vector Double -> (V.Vector (U.Vector Double), V.Vector (V.Vector (U.Vector Double)))
{-# INLINE posterior' #-}
posterior' model _ xs alphas betas cs = (gammas, xis)
  where
    gammas = V.zipWith3 (\a b c -> a >.> b >/ c) alphas betas (G.convert cs)
    xis    = V.zipWith3 (\a b x -> let w' = V.zipWith (.>) (G.convert a) w
                                   in V.map (V.unsafeIndex phi' x >.> b >.>) w')
               alphas (V.unsafeTail betas) (G.convert $ U.unsafeTail xs)
    w    = transitionDist' model
    phi' = emissionDistT' model
