{-# LANGUAGE
    MultiParamTypeClasses,
    FlexibleContexts, FlexibleInstances,
    UndecidableInstances
  #-}

module Data.Random.Distribution.Extra
  ( PMF
  , pmf
  ) where

import           Data.Maybe                                ( fromMaybe )
import           Data.Tuple                                ( swap )
import           Data.Random.Distribution                  ( Distribution )
import qualified Data.Random.Distribution.Categorical as C ( Categorical, toList )

class Distribution d t => PMF d t where
  -- | Probability mass function for discrete distributions
  pmf :: d t -> t -> Double

instance (Real p, Eq a, Distribution (C.Categorical p) a) => PMF (C.Categorical p) a where
  pmf d a = realToFrac $ fromMaybe 0 $ lookup a dict
    where dict = map swap $ C.toList d
