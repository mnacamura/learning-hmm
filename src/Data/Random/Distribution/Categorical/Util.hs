{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

module Data.Random.Distribution.Categorical.Util () where

import Data.Maybe (fromMaybe)
import Data.Random.Distribution (PDF, pdf)
import Data.Random.Distribution.Categorical (Categorical, toList)
import Data.Tuple (swap)

instance Eq a => PDF (Categorical Double) a where
  pdf cat a = fromMaybe 0 (lookup a $ map swap $ toList cat)
