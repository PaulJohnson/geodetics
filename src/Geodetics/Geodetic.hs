module Geodetics.Geodetic (
   Geodetic (..),
   toLocal,
   toWGS84,
   antipode,
   geometricalDistance,
   geometricalDistanceSq
) where


import Data.Char (chr)
import Data.Function
import Data.Monoid
import Geodetics.Altitude
import Geodetics.Ellipsoids
import Numeric
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P

-- | Defines a three-D position on or around the Earth using latitude, longitude and altitude
-- with respect to a specified ellipsoid, with positive directions being North and East.
-- The default "show" instance gives position in degrees to 5 decimal
-- places, which is a resolution of about 1m on the Earth's surface. Internally latitude
-- and longitude are stored as double precision radians. Convert to degrees using e.g. 
-- @latitude g /~ degree@.
-- 
-- The altitude assumes that the local height datum is always co-incident with the ellipsoid.
-- In practice the "mean sea level" (the usual height datum) can be tens of meters above or below
-- the ellipsoid, and two ellipsoids can differ by similar amounts. However in practice the altitude
-- is usually known with reference to a local datum regardless of the ellipsoid in use.
-- 
-- There is no "Eq" instance because comparing two arbitrary co-ordinates on the Earth
-- is a non-trivial exercise. Clearly if all the parameters are equal on the same ellipsoid
-- then they are indeed in the same place. However if different ellipsoids are used then
-- two co-ordinates with different numbers can still refer to the same physical location. 
-- If you want to find out if two co-ordinates are the same to within a given tolerance then
-- use "geometricDistance" (or its squared variant to avoid an extra @sqrt@ operation).
data (Ellipsoid e) => Geodetic e = Geodetic {
   latitude, longitude :: Angle Double,
   geoAlt :: Length Double,
   ellipsoid :: e
}

instance (Ellipsoid e) => Show (Geodetic e) where
   show g = concat [
      letter "SN" (latitude g), " ", showAngle (abs $ latitude g), ", ",
      letter "WE" (longitude g), " ", showAngle (abs $ longitude g), ", ", 
      show (altitude g), " ", show (ellipsoid g)]
      where letter s n = [s !! (if n < _0 then 0 else 1)] 

         

-- | Show an angle as degrees, minutes and seconds.
showAngle :: Dimensionless Double -> String
showAngle a = 
   show d ++ [chr 0xB0, ' '] ++ show (P.abs m1) ++ "' " ++ showFFloat (Just 2) (P.abs s) "''"
   where
      d, m1 :: Int
      (d, f1) = P.properFraction $ P.abs $ a /~ degree
      (m1, f2) = P.properFraction $ f1 P.* 60
      s = f2 P.* 60
         

instance (Ellipsoid e) => HasAltitude (Geodetic e) where
   altitude = geoAlt
   setAltitude h g = g{geoAlt = h}

   
   
-- | The point on the Earth diametrically opposite the argument, with the same altitude.
antipode :: (Ellipsoid e) => Geodetic e -> Geodetic e
antipode g = Geodetic lat long (geoAlt g) (ellipsoid g)
   where
      lat = negate $ latitude g
      long' = longitude g - 180 *~ degree
      long | long' < _0  = long' + 360 *~ degree
           | otherwise  = long' 


-- | Distance from the surface to the Z axis straight down.
normal :: (Ellipsoid e) => e -> Angle Double -> Length Double
normal e lat = majorRadius e / sqrt (_1 - eccentricity2 e * sin lat ^ pos2)

   
   
-- | Convert a geodetic coordinate into earth centred, relative to the ellipsoid in use.
geoToEarth :: (Ellipsoid e) => Geodetic e -> ECEF
geoToEarth geo = (
      (n + h) * coslat * coslong,
      (n + h) * coslat * sinlong,
      (n * (_1 - eccentricity2 e) + h) * sinlat)
   where 
      n = normal e $ latitude geo
      e = ellipsoid geo
      coslat = cos $ latitude geo
      coslong = cos $ longitude geo
      sinlat = sin $ latitude geo
      sinlong = sin $ longitude geo
      h = altitude geo


-- | Convert an earth centred coordinate into a geodetic coordinate on 
-- the specified geoid.
--
-- Uses the closed form solution of H. Vermeille: Direct transformation from 
-- geocentric coordinates to geodetic coordinates.  Journal of Geodesy
-- Volume 76, Number 8 (2002), 451-454
earthToGeo :: (Ellipsoid e) => e -> ECEF -> (Angle Double, Angle Double, Length Double)
earthToGeo e (x,y,z) = (phi, atan2 y x, sqrt (l ^ pos2 + p2) - norm)
   where
      -- Naming: numeric suffix inicates power. Hence x2 = x * x, x3 = x2 * x, etc.
      p2 = x ^ pos2 + y ^ pos2
      a = majorRadius e
      a2 = a ^ pos2
      e2 = eccentricity2 e
      e4 = e2 ^ pos2
      zeta = (_1-e2) * (z ^ pos2 / a2)
      rho = (p2 / a2 + zeta - e4) / _6
      rho2 = rho ^ pos2
      rho3 = rho * rho2
      s = e4 * zeta * p2 / (_4 * a2)
      t = cbrt (s + rho3 + sqrt (s * (s + _2 * rho3)))
      u = rho + t + rho2 / t
      v = sqrt (u ^ pos2 + e4 * zeta)
      w = e2 * (u + v - zeta) / (_2 * v)
      kappa = _1 + e2 * (sqrt (u + v + w ^ pos2) + w) / (u + v)
      phi = atan (kappa * z / sqrt p2)
      norm = normal e phi
      l = z + e2 * norm * sin phi


-- | Convert a position from any geodetic to a another one, assuming local altitude stays constant.
toLocal :: (Ellipsoid e1, Ellipsoid e2) => e2 -> Geodetic e1 -> Geodetic e2
toLocal e2 g = Geodetic lat lon alt e2
   where
      alt = altitude g
      (lat, lon, _) = earthToGeo e2 $ applyHelmert h $ geoToEarth g
      h = helmert (ellipsoid g) `mappend` inverseHelmert (helmert e2)

-- | Convert a position from any geodetic to WGS84, assuming local altitude stays constant.
toWGS84 :: (Ellipsoid e) => Geodetic e -> Geodetic WGS84
toWGS84 g = Geodetic lat lon alt WGS84
   where
      alt = altitude g
      (lat, lon, _) = earthToGeo WGS84 $ applyHelmert h $ geoToEarth g
      h = helmert (ellipsoid g)


-- | The absolute distance in a straight line between two geodetic 
-- points. They must be on the same ellipsoid.
-- Note that this is not the Great Circle distance taken by following 
-- the curvature of the earth.
geometricalDistance :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Length Double
geometricalDistance g1 g2 = sqrt $ geometricalDistanceSq g1 g2

-- | The square of the absolute distance. Comes out as "Area" type of course.
geometricalDistanceSq :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Area Double
geometricalDistanceSq g1 g2 = (x1-x2) ^ pos2 + (y1-y2) ^ pos2 + (z1-z2) ^ pos2
   where
      (x1,y1,z1) = geoToEarth g1
      (x2,y2,z2) = geoToEarth g2
