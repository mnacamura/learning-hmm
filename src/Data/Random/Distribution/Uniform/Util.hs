{-# LANGUAGE FlexibleContexts, FlexibleInstances, MultiParamTypeClasses #-}

module Data.Random.Distribution.Uniform.Util () where

import Data.Number.LogFloat (LogFloat, logFloat, fromLogFloat)
import Data.Random.Distribution (Distribution)
import Data.Random.Distribution.Uniform (
    StdUniform (..), Uniform (..), doubleUniform
  )
import Data.Random (rvarT)
import Data.Random.Source (getRandomDouble)

instance Distribution Uniform LogFloat where
  rvarT (Uniform a b) = fmap logFloat (doubleUniform a' b')
    where a' = fromLogFloat a
          b' = fromLogFloat b

instance Distribution StdUniform LogFloat where
  rvarT _ = fmap logFloat getRandomDouble
