{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE RoleAnnotations #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}

{- | An Ellipsoid is a reasonable best fit for the surface of the
Earth over some defined area. WGS84 is the standard used for the whole
of the Earth. Other Ellipsoids are considered a best fit for some
specific area.
-}

module Geodetics.Ellipsoids (
   -- * Conversion constants
   degree,
   arcminute,
   arcsecond,
   kilometer,
   -- * Helmert transform between geodetic reference systems
   Helmert (..),
   inverseHelmert,
   ECEF,
   applyHelmert,
   -- * Ellipsoid models of the Geoid
   Ellipsoid (..),
   WGS84 (..),
   LocalEllipsoid (..),
   flattening,
   minorRadius,
   eccentricity2,
   eccentricity'2,
   -- * Auxiliary latitudes and related Values
   normal,
   latitudeRadius,
   meridianRadius,
   primeVerticalRadius,
   isometricLatitude,
   -- * Tiny linear algebra library for 3D vectors
   Vec3,
   Matrix3,
   add3,
   scale3,
   negate3,
   transform3,
   invert3,
   trans3,
   dot3,
   cross3
) where


-- | All angles in this library are in radians. This is one degree in radians.
degree :: Double
degree = pi/180

-- | One arc-minute in radians.
arcminute :: Double
arcminute = degree / 60

-- | One arc-second in radians.
arcsecond :: Double
arcsecond = arcminute / 60


-- | All distances in this library are in meters. This is one kilometer in meters.
kilometer :: Double
kilometer = 1000


-- | Small integers, specialised to @Int@, used for raising to powers.
--
-- If you say @x^2@ then Haskell complains that the @2@ has ambiguous type, so you
-- need to say @x^(2::Int)@ to disambiguate it. This gets tedious in complex formulae.

-- | 3d vector as @(X,Y,Z)@.
type Vec3 a = (a,a,a)

-- | 3x3 transform matrix for Vec3.
type Matrix3 a = Vec3 (Vec3 a)


-- | Multiply a vector by a scalar.
scale3 :: (Num a) =>  Vec3 a -> a -> Vec3 a
scale3 (x,y,z) s = (x*s, y*s, z*s)


-- | Negation of a vector.
negate3 :: (Num a) => Vec3 a -> Vec3 a
negate3 (x,y,z) = (negate x, negate y, negate z)

-- | Add two vectors
add3 :: (Num a) => Vec3 a -> Vec3 a -> Vec3 a
add3 (x1,y1,z1) (x2,y2,z2) = (x1+x2, y1+y2, z1+z2)


-- | Multiply a matrix by a vector in the Dimensional type system.
transform3 :: (Num a) =>
   Matrix3 a -> Vec3 a -> Vec3 a
transform3 (tx,ty,tz) v = (t tx v, t ty v, t tz v)
   where
      t (x1,y1,z1) (x2,y2,z2) = x1*x2 + y1*y2 + z1*z2


-- | Inverse of a 3x3 matrix.
invert3 :: (Fractional a) => Matrix3 a -> Matrix3 a
invert3 ((x1,y1,z1),
         (x2,y2,z2),
         (x3,y3,z3)) =
      ((det2 y2 z2 y3 z3 / det, det2 z1 y1 z3 y3 / det, det2 y1 z1 y2 z2 / det),
       (det2 z2 x2 z3 x3 / det, det2 x1 z1 x3 z3 / det, det2 z1 x1 z2 x2 / det),
       (det2 x2 y2 x3 y3 / det, det2 y1 x1 y3 x3 / det, det2 x1 y1 x2 y2 / det))
   where
      det = (x1*y2*z3 + y1*z2*x3 + z1*x2*y3) - (z1*y2*x3 + y1*x2*z3 + x1*z2*y3)
      det2 a b c d = a*d - b*c

-- | Transpose of a 3x3 matrix.
trans3 :: Matrix3 a -> Matrix3 a
trans3 ((x1,y1,z1),(x2,y2,z2),(x3,y3,z3)) = ((x1,x2,x3),(y1,y2,y3),(z1,z2,z3))


-- | Dot product of two vectors
dot3 :: (Num a) => Vec3 a -> Vec3 a -> a
dot3 (x1,y1,z1) (x2,y2,z2) = x1*x2 + y1*y2 + z1*z2

-- | Cross product of two vectors
cross3 :: (Num a) => Vec3 a -> Vec3 a -> Vec3 a
cross3 (x1,y1,z1) (x2,y2,z2) = (y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2)


-- | The 7 parameter Helmert transformation. The monoid instance allows composition.
data Helmert = Helmert {
   cX, cY, cZ :: Double,  -- ^ Offset in meters
   helmertScale :: Double,  -- ^ Parts per million
   rX, rY, rZ :: Double  -- ^ Rotation around each axis in radians.
} deriving (Eq, Show)

instance Semigroup Helmert where
    h1 <> h2 = Helmert (cX h1 + cX h2) (cY h1 + cY h2) (cZ h1 + cZ h2)
                       (helmertScale h1 + helmertScale h2)
                       (rX h1 + rX h2) (rY h1 + rY h2) (rZ h1 + rZ h2)

instance Monoid Helmert where
   mempty = Helmert 0 0 0 0 0 0 0
   mappend = (<>)

-- | The inverse of a Helmert transformation.
inverseHelmert :: Helmert -> Helmert
inverseHelmert h = Helmert (negate $ cX h) (negate $ cY h) (negate $ cZ h)
                           (negate $ helmertScale h)
                           (negate $ rX h) (negate $ rY h) (negate $ rZ h)


-- | Earth-centred, Earth-fixed coordinates as a vector. The origin and axes are
-- not defined: use with caution.
type ECEF = Vec3 Double

-- | Apply a Helmert transformation to earth-centered coordinates.
applyHelmert:: Helmert -> ECEF -> ECEF
applyHelmert h (x,y,z) = (
      cX h + s * (                x - rZ h * y + rY h * z),
      cY h + s * (        rZ h  * x +        y - rX h * z),
      cZ h + s * (negate (rY h) * x + rX h * y +        z))
   where
      s = 1 + helmertScale h * 1e-6


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
   majorRadius :: a -> Double
   flatR :: a -> Double
      -- ^ Inverse of the flattening.
   helmert :: a -> Helmert  -- ^ The Helmert parameters relative to WGS84,
   helmertToWGS84 :: a -> ECEF -> ECEF
      -- ^ The Helmert transform that will convert a position wrt
      -- this ellipsoid into a position wrt WGS84.
   helmertToWGS84 e = applyHelmert (helmert e)
   helmertFromWGS84 :: a -> ECEF -> ECEF
      -- ^ And its inverse.
   helmertFromWGS84 e = applyHelmert (inverseHelmert $ helmert e)


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
   majorRadius _ = 6378137.0
   flatR _ = 298.257223563
   helmert _ = mempty
   helmertToWGS84 _ = id
   helmertFromWGS84 _ = id


-- | Ellipsoids other than WGS84, used within a defined geographical area where
-- they are a better fit to the local geoid. Can also be used for historical ellipsoids.
--
-- The @Show@ instance just returns the name.
-- Creating two different local ellipsoids with the same name is a Bad Thing.
data LocalEllipsoid = LocalEllipsoid {
   nameLocal :: String,
   majorRadiusLocal :: Double,
   flatRLocal :: Double,
   helmertLocal :: Helmert
} deriving (Eq)

instance Show LocalEllipsoid where
    show = nameLocal

instance Ellipsoid LocalEllipsoid where
   majorRadius = majorRadiusLocal
   flatR = flatRLocal
   helmert = helmertLocal


-- | Flattening (f) of an ellipsoid.
flattening :: (Ellipsoid e) => e -> Double
flattening e = 1 / flatR e

-- | The minor radius of an ellipsoid in meters.
minorRadius :: (Ellipsoid e) => e -> Double
minorRadius e = majorRadius e * (1 - flattening e)


-- | The eccentricity squared of an ellipsoid.
eccentricity2 :: (Ellipsoid e) => e -> Double
eccentricity2 e = 2 * f - (f * f) where f = flattening e

-- | The second eccentricity squared of an ellipsoid.
eccentricity'2 :: (Ellipsoid e) => e -> Double
eccentricity'2 e = (f * (2 - f)) / (1 - f * f) where f = flattening e


-- | Distance in meters from the surface at the specified latitude to the
-- axis of the Earth straight down. Also known as the radius of
-- curvature in the prime vertical, and often denoted @N@.
normal :: (Ellipsoid e) => e -> Double -> Double
normal e lat = majorRadius e / sqrt (1 - eccentricity2 e * sin lat ^ (2 :: Int))


-- | Radius of the circle of latitude: the distance from a point
-- at that latitude to the axis of the Earth, in meters.
latitudeRadius :: (Ellipsoid e) => e -> Double -> Double
latitudeRadius e lat = normal e lat * cos lat


-- | Radius of curvature in the meridian at the specified latitude, in meters
-- Often denoted @M@.
meridianRadius :: (Ellipsoid e) => e -> Double -> Double
meridianRadius e lat =
   majorRadius e * (1 - eccentricity2 e)
   / sqrt ((1 - eccentricity2 e * sin lat ** 2) ** 3)


-- | Radius of curvature of the ellipsoid perpendicular to the meridian at the specified latitude, in meters.
primeVerticalRadius :: (Ellipsoid e) => e -> Double -> Double
primeVerticalRadius e lat =
   majorRadius e / sqrt (1 - eccentricity2 e * sin lat ^ (2 :: Int))


-- | The isometric latitude. The isometric latitude is conventionally denoted by ψ
-- (not to be confused with the geocentric latitude): it is used in the development
-- of the ellipsoidal versions of the normal Mercator projection and the Transverse
-- Mercator projection. The name "isometric" arises from the fact that at any point
-- on the ellipsoid equal increments of ψ and longitude λ give rise to equal distance
-- displacements along the meridians and parallels respectively.
isometricLatitude :: (Ellipsoid e) => e -> Double -> Double
isometricLatitude ellipse lat = atanh sinLat - e * atanh (e * sinLat)
   where
      sinLat = sin lat
      e = sqrt $ eccentricity2 ellipse
