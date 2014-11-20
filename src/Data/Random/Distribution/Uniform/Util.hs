{-# LANGUAGE FlexibleContexts, FlexibleInstances, MultiParamTypeClasses #-}

module Data.Random.Distribution.Uniform.Util () where

import Control.Applicative ((<$>))
import Data.Number.LogFloat (LogFloat, logFloat, fromLogFloat)
import Data.Random.Distribution (Distribution)
import Data.Random.Distribution.Uniform -- (StdUniform(..), Uniform(..), doubleUniform)
import Data.Random (rvarT)
import Data.Random.Source (getRandomDouble)

instance Distribution Uniform LogFloat where
  rvarT (Uniform a b) = do x <- doubleUniform (fromLogFloat a) (fromLogFloat b)
                           return $ logFloat x

instance Distribution StdUniform LogFloat where
  rvarT _ = logFloat <$> getRandomDouble
