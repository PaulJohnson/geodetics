module Geodetics.Geodetic (
   -- * Geodetic Coordinates
   Geodetic (..),
   readGroundPosition,
   toLocal,
   toWGS84,
   antipode,
   geometricalDistance,
   geometricalDistanceSq,
   groundDistance,
   properAngle,
   showAngle,
   -- * Earth Centred Earth Fixed Coordinates
   ECEF,
   geoToEarth,
   earthToGeo,
   -- * Re-exported for convenience
   WGS84 (..)
) where


import Data.Char (chr)
import Data.Maybe
import Geodetics.Altitude
import Geodetics.Ellipsoids
import Geodetics.LatLongParser
import Text.ParserCombinators.ReadP

-- | Defines a three-D position on or around the Earth using latitude,
-- longitude and altitude with respect to a specified ellipsoid, with
-- positive directions being North and East.  The default "show"
-- instance gives position in degrees, minutes and seconds to 5 decimal
-- places, which is a
-- resolution of about 1m on the Earth's surface. Internally latitude
-- and longitude are stored as double precision radians. Convert to
-- degrees using e.g.  @latitude g /~ degree@.
--
-- The functions here deal with altitude by assuming that the local
-- height datum is always co-incident with the ellipsoid in use,
-- even though the \"mean sea level\" (the usual height datum) can be tens
-- of meters above or below the ellipsoid, and two ellipsoids can
-- differ by similar amounts. This is because the altitude is
-- usually known with reference to a local datum regardless of the
-- ellipsoid in use, so it is simpler to preserve the altitude across
-- all operations. However if
-- you are working with ECEF coordinates from some other source then
-- this may give you the wrong results, depending on the altitude
-- correction your source has used.
--
-- There is no "Eq" instance because comparing two arbitrary
-- co-ordinates on the Earth is a non-trivial exercise. Clearly if all
-- the parameters are equal on the same ellipsoid then they are indeed
-- in the same place. However if different ellipsoids are used then
-- two co-ordinates with different numbers can still refer to the same
-- physical location.  If you want to find out if two co-ordinates are
-- the same to within a given tolerance then use "geometricDistance"
-- (or its squared variant to avoid an extra @sqrt@ operation).
data Geodetic e = Geodetic {
   latitude, longitude :: Double,  -- ^ In radians.
   geoAlt :: Double,  -- ^ In meters.
   ellipsoid :: e
}

instance (Ellipsoid e) => Show (Geodetic e) where
   show g = concat [
      showAngle (abs $ latitude g),  " ", letter "SN" (latitude g),  ", ",
      showAngle (abs $ longitude g), " ", letter "WE" (longitude g), ", ",
      show (altitude g), " ", show (ellipsoid g)]
      where letter s n = [s !! (if n < 0 then 0 else 1)]


-- | Read the latitude and longitude of a ground position and
-- return a Geodetic position on the specified ellipsoid.
--
-- The latitude and longitude may be in any of the following formats.
-- The comma between latitude and longitude is optional in all cases.
-- Latitude must always be first.
--
-- * Signed decimal degrees: 34.52327, -46.23234
--
-- * Decimal degrees NSEW: 34.52327N, 46.23234W
--
-- * Degrees and decimal minutes (units optional): 34° 31.43' N, 46° 13.92'
--
-- * Degrees, minutes and seconds (units optional): 34° 31' 23.52\" N, 46° 13' 56.43\" W
--
-- * DDDMMSS format with optional leading zeros: 343123.52N, 0461356.43W
readGroundPosition :: (Ellipsoid e) => e -> String -> Maybe (Geodetic e)
readGroundPosition e str =
   case map fst $ filter (null . snd) $ readP_to_S latLong str of
      [] -> Nothing
      (lat,long) : _ -> Just $ groundPosition $ Geodetic (lat * degree) (long * degree) undefined e


-- | Show an angle as degrees, minutes and seconds to two decimal places.
showAngle :: Double -> String
showAngle a
   | isNaN a        = "NaN"  -- Not a Nangle
   | isInfinite a   = sgn ++ "Infinity"
   | otherwise      = concat [sgn, show d, [chr 0xB0, ' '],
                              show m, "\8242 ",
                              show s, ".", dstr, "\8243" ]
   where
      sgn = if a < 0 then "-" else ""
      centisecs :: Integer
      centisecs = abs $ round $ (a / (arcsecond / 100))
      (d, m1) = centisecs `divMod` 360000
      (m, s1) = m1 `divMod` 6000   -- hundredths of arcsec per arcmin
      (s, ds) = s1 `divMod` 100
      dstr = reverse $ take 2 $ reverse (show ds) ++ "00" -- Decimal fraction with zero padding.


instance (Ellipsoid e) => HasAltitude (Geodetic e) where
   altitude = geoAlt
   setAltitude h g = g{geoAlt = h}



-- | The point on the Earth diametrically opposite the argument, with
-- the same altitude.
antipode :: Geodetic e -> Geodetic e
antipode g = Geodetic lat long (geoAlt g) (ellipsoid g)
   where
      lat = negate $ latitude g
      long' = longitude g - 180 * degree
      long | long' < 0  = long' + 360 * degree
           | otherwise  = long'



-- | Convert a geodetic coordinate into earth centered, relative to the
-- ellipsoid in use.
geoToEarth :: (Ellipsoid e) => Geodetic e -> ECEF
geoToEarth geo = (
      (n + h) * coslat * coslong,
      (n + h) * coslat * sinlong,
      (n * (1 - eccentricity2 e) + h) * sinlat)
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
-- Uses the closed form solution of H. Vermeille: Direct
-- transformation from geocentric coordinates to geodetic coordinates.
-- Journal of Geodesy Volume 76, Number 8 (2002), 451-454. Result is in the form
-- @(latitude, longitude, altitude)@.
earthToGeo :: (Ellipsoid e) => e -> ECEF -> (Double, Double, Double)
earthToGeo e (x,y,z) = (phi, atan2 y x, sqrt (l ^ _2 + p2) - norm)
   where
      -- Naming: numeric suffix inicates power. Hence x2 = x * x, x3 = x2 * x, etc.
      p2 = x * x + y * y
      a = majorRadius e
      a2 = a * a
      e2 = eccentricity2 e
      e4 = e2 * e2
      zeta = (1-e2) * (z * z / a2)
      rho = (p2 / a2 + zeta - e4) / 6
      rho2 = rho * rho
      rho3 = rho * rho2
      s = e4 * zeta * p2 / (4 * a2)
      t = (s + rho3 + sqrt (s * (s + 2 * rho3))) ** (1/3) -- Cube root
      u = rho + t + rho2 / t
      v = sqrt (u * u + e4 * zeta)
      w = e2 * (u + v - zeta) / (2 * v)
      kappa = 1 + e2 * (sqrt (u + v + w * w) + w) / (u + v)
      phi = atan (kappa * z / sqrt p2)
      norm = normal e phi
      l = z + e2 * norm * sin phi


-- | Convert a position from any geodetic to another one, assuming local altitude stays constant.
toLocal :: (Ellipsoid e1, Ellipsoid e2) => e2 -> Geodetic e1 -> Geodetic e2
toLocal e2 g = Geodetic lat lon alt e2
   where
      alt = altitude g
      (lat, lon, _) = earthToGeo e2 $ applyHelmert h $ geoToEarth g
      h = helmert (ellipsoid g) `mappend` inverseHelmert (helmert e2)

-- | Convert a position from any geodetic to WGS84, assuming local
-- altitude stays constant.
toWGS84 :: (Ellipsoid e) => Geodetic e -> Geodetic WGS84
toWGS84 g = Geodetic lat lon alt WGS84
   where
      alt = altitude g
      (lat, lon, _) = earthToGeo WGS84 $ applyHelmert h $ geoToEarth g
      h = helmert (ellipsoid g)


-- | The absolute distance in a straight line between two geodetic
-- points. They must be on the same ellipsoid.
-- Note that this is not the geodetic distance taken by following
-- the curvature of the earth.
geometricalDistance :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Double
geometricalDistance g1 g2 = sqrt $ geometricalDistanceSq g1 g2

-- | The square of the absolute distance.
geometricalDistanceSq :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Double
geometricalDistanceSq g1 g2 = (x1-x2) ^ _2 + (y1-y2) ^ _2 + (z1-z2) ^ _2
   where
      (x1,y1,z1) = geoToEarth g1
      (x2,y2,z2) = geoToEarth g2


-- | The shortest ellipsoidal distance between two points on the
-- ground with reference to the same ellipsoid. Altitude is ignored.
--
-- The results are the distance between the points, the bearing of
-- the second point from the first, and (180 degrees - the bearing
-- of the first point from the second).
--
-- The algorithm can fail to converge where the arguments are near to
-- antipodal. In this case it returns @Nothing@.
--
-- Uses Vincenty's formula. \"Direct and inverse solutions of
-- geodesics on the ellipsoid with application of nested
-- equations\". T. Vincenty. Survey Review XXII 176, April
-- 1975. <http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf>
groundDistance :: (Ellipsoid e) => Geodetic e -> Geodetic e ->
                  Maybe (Double, Double, Double)
groundDistance p1 p2 = do
     (_, (lambda, (cos2Alpha, delta, sinDelta, cosDelta, cos2DeltaM))) <-
       listToMaybe $ dropWhile converging $ take 100 $ zip lambdas $ drop 1 lambdas
     let
       uSq = cos2Alpha * (a^ _2 - b^ _2) / b^ _2
       bigA = 1 + uSq/16384 * (4096 + uSq * ((-768) + uSq * ((320 - 175*uSq))))
       bigB =     uSq/1024  * (256  + uSq * ((-128) + uSq * ((74 -  47* uSq))))
       deltaDelta =
         bigB * sinDelta * (cos2DeltaM +
                             bigB/4 * (cosDelta * (2 * cos2DeltaM^ _2 - 1)
                                       - bigB/6 * cos2DeltaM * (4 * sinDelta^ _2 - 3)
                                          * (4 * cos2DeltaM - 3)))
       s = b * bigA * (delta - deltaDelta)
       alpha1 = atan2(cosU2 * sin lambda) (cosU1 * sinU2 - sinU1 * cosU2 * cos lambda)
       alpha2 = atan2(cosU1 * sin lambda) (cosU1 * sinU2 * cos lambda - sinU1 * cosU2)
     return (s, alpha1, alpha2)
  where
    f = flattening $ ellipsoid p1
    a = majorRadius $ ellipsoid p1
    b = minorRadius $ ellipsoid p1
    l = abs $ longitude p1 - longitude p2
    u1 = atan ((1-f) * tan (latitude p1))
    u2 = atan ((1-f) * tan (latitude p2))
    sinU1 = sin u1
    cosU1 = cos u1
    sinU2 = sin u2
    cosU2 = cos u2

    nextLambda lambda = (lambda1, (cos2Alpha, delta, sinDelta, cosDelta, cos2DeltaM))
      where
        sinLambda = sin lambda
        cosLambda = cos lambda
        sinDelta = sqrt((cosU2 * sinLambda) ^ _2 +
                        (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ^ _2)
        cosDelta = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        delta = atan2 sinDelta cosDelta
        sinAlpha = if sinDelta == 0 then 0 else cosU1 * cosU2 * sinLambda / sinDelta
        cos2Alpha = 1 - sinAlpha ^ _2
        cos2DeltaM = if cos2Alpha == 0
                     then 0
                     else cosDelta - 2 * sinU1 * sinU2 / cos2Alpha
        c = (f/16) * cos2Alpha * (4 + f * (4 - 3 * cos2Alpha))
        lambda1 = l + (1-c) * f * sinAlpha
                  * (delta + c * sinDelta
                     * (cos2DeltaM + c * cosDelta *(2 * cos2DeltaM ^ _2 - 1)))
    lambdas = iterate (nextLambda . fst) (l, undefined)
    converging ((l1,_),(l2,_)) = abs (l1 - l2) > 1e-14


-- | Add or subtract multiples of 2*pi so that for all @t@, @-pi < properAngle t < pi@.
properAngle :: Double -> Double
properAngle t
   | r1 <= negate pi    = r1 + pi2
   | r1 > pi            = r1 - pi2
   | otherwise          = r1
   where
      pf :: Double -> (Int, Double)
      pf = properFraction  -- Shut up GHC warning about defaulting to Integer.
      (_,r) = pf (t/pi2)
      r1 = r * pi2
      pi2 = pi * 2
