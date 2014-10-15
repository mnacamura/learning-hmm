{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

module Data.Random.Distribution.Categorical.Util () where

import Data.Random.Distribution (PDF, pdf)
import Data.Random.Distribution.Categorical (Categorical, toList, totalWeight)
import Data.Tuple (swap)

instance Eq a => PDF (Categorical Double) a where
  pdf cat a = case lookup a $ map swap $ toList cat of
                Nothing -> 0
                Just p  -> p / totalWeight cat
