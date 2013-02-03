{- | An Ellipsoid is a reasonable best fit for the surface of the 
Earth over some defined area. WGS84 is the standard used for the whole
of the Earth. Other Ellipsoids are considered a best fit for some
specific area.
-}

module Geodetics.Ellipsoids where

import Data.Monoid
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P


-- | Earth-centred, Earth-fixed coordinates as @(X,Y,Z)@. The origin and axes are
-- not defined: use with caution.
type ECEF = (Length Double, Length Double, Length Double)


-- | The 7 parameter Helmert transformation. The monoid instance allows composition.
data Helmert = Helmert {
   cX, cY, cZ :: Length Double,
   helmertScale :: Dimensionless Double,  -- ^ Parts per million
   rX, rY, rZ :: Dimensionless Double } deriving (Eq, Show)

instance Monoid Helmert where
   mempty = Helmert (0 *~ meter) (0 *~ meter) (0 *~ meter) _1 _0 _0 _0
   mappend h1 h2 = Helmert (cX h1 + cX h2) (cY h1 + cY h2) (cZ h1 + cZ h2)
                           (helmertScale h1 + helmertScale h2)
                           (rX h1 + rX h2) (rY h1 + rY h2) (rZ h1 + rZ h2)


-- | The inverse of a Helmert transformation.
inverseHelmert :: Helmert -> Helmert
inverseHelmert h = Helmert (negate $ cX h) (negate $ cY h) (negate $ cZ h) 
                           (negate $ helmertScale h) 
                           (negate $ rX h) (negate $ rY h) (negate $ rZ h)


-- | Apply a Helmert transformation to earth-centered coordinates.
applyHelmert:: Helmert -> ECEF -> ECEF
applyHelmert h (x,y,z) = (
      cX h + s * (                  x - rZ h * y + rY h * z),
      cY h + s * (          rZ h  * x +        y - rX h * z),
      cZ h + s * ((negate $ rY h) * x + rX h * y +        z))
   where
      s = _1 + helmertScale h * (1e-6 *~ one)


-- | An Ellipsoid is defined by the major radius and the inverse flattening (which define its shape), 
-- and its Helmert transform relative to WGS84 (which defines its position and orientation).
--
-- The inclusion of the Helmert parameters relative to WGS84 actually make this a Terrestrial 
-- Reference Frame (TRF), but the term "Ellipsoid" will be used in this library for readability.
--
-- Minimum definition: @majorRadius@, @flatR@ & @helmert@.
-- 
-- Laws:
-- 
-- > helmertToWGS84 = applyHelmert . helmert
-- > helmertFromWGS84 e . helmertToWGS84 e = id
class (Show a, Eq a) => Ellipsoid a where
   majorRadius :: a -> Length Double
   flatR :: a -> Dimensionless Double
      -- ^ Inverse of the flattening.
   helmert :: a -> Helmert
   helmertToWSG84 :: a -> ECEF -> ECEF
      -- ^ The Helmert transform that will convert a position wrt 
      -- this ellipsoid into a position wrt WGS84.
   helmertToWSG84 e ec = applyHelmert (helmert e) ec
   helmertFromWSG84 :: a -> ECEF -> ECEF
      -- ^ And its inverse.
   helmertFromWSG84 e ec = applyHelmert (inverseHelmert $ helmert e) ec


-- | The WGS84 geoid, major radius 6378137.0 meters, flattening = 1 / 298.257223563
-- as defined in \"Technical Manual DMA TM 8358.1 - Datums, Ellipsoids, Grids, and 
-- Grid Reference Systems\" at the National Geospatial-Intelligence Agency (NGA).
-- 
-- The WGS84 has a special place in this library as the standard Ellipsoid against
-- which all others are defined.
data WGS84 = WGS84

instance Eq WGS84 where _ == _ = True

instance Show WGS84 where
   show _ = "WGS84"
   
instance Ellipsoid WGS84 where
   majorRadius _ = 6378137.0 *~ meter
   flatR _ = 298.257223563 *~ one
   helmert _ = mempty
   helmertToWSG84 _ = id
   helmertFromWSG84 _ = id
   
   
-- | Ellipsoids other than WGS84, used within a defined geographical area where
-- they are a better fit to the local geoid. Can also be used for historical ellipsoids.
--
-- The @Show@ instance just shows the name.
-- Creating two different local ellipsoids with the same name is a Bad Thing.
data LocalEllipsoid = LocalE {
   nameLocal :: String,
   majorRadiusLocal :: Length Double,
   flatRLocal :: Dimensionless Double,
   helmertLocal :: Helmert } deriving (Eq)

instance Show LocalEllipsoid where
    show = nameLocal  

instance Ellipsoid LocalEllipsoid where
   majorRadius = majorRadiusLocal
   flatR = flatRLocal
   helmert = helmertLocal


-- | Flattening (f) of an ellipsoid.
flattening :: (Ellipsoid e) => e -> Dimensionless Double
flattening e = _1 / flatR e

-- | The minor radius of an ellipsoid.
minorRadius :: (Ellipsoid e) => e -> Length Double
minorRadius e = majorRadius e * (_1 - flattening e)


-- | The eccentricity squared of an ellipsoid.
eccentricity2 :: (Ellipsoid e) => e -> Dimensionless Double
eccentricity2 e = _2 * f - (f * f) where f = flattening e

-- | The second eccentricity squared of an ellipsoid.
eccentricity'2 :: (Ellipsoid e) => e -> Dimensionless Double
eccentricity'2 e = (f * (_2 - f)) / (_1 - f * f) where f = flattening e

