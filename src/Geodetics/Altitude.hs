module Geodetics.Altitude (
  HasAltitude (..)
) where


-- | All geographical coordinate systems need the concept of
-- altitude above a reference point, usually associated with
-- local sea level.
-- 
-- Minimum definition: altitude, setAltitude.
class HasAltitude a where
   altitude :: a -> Double
   setAltitude :: Double -> a -> a
   -- | Set altitude to zero.
   groundPosition :: a -> a
   groundPosition = setAltitude 0
