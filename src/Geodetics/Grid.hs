{-# LANGUAGE FunctionalDependencies #-}

module Geodetics.Grid (
   -- * Grid types
   GridClass (..),
   GridPoint (..),
   GridOffset (..),
   -- * Grid operations
   polarOffset,
   offsetScale,
   offsetNegate,
   applyOffset,
   offsetDistance,
   offsetDistanceSq,
   offsetBearing,
   gridOffset,
   -- * Unsafe conversion
   unsafeGridCoerce,
   -- * Utility functions for grid references
   fromGridDigits,
   toGridDigits,
   GridUnit (..)
) where

import Data.Char
import Geodetics.Altitude
import Geodetics.Ellipsoids
import Geodetics.Geodetic

-- | A Grid is a two-dimensional projection of the ellipsoid onto a plane. Any given type of grid can
-- usually be instantiated with parameters such as a tangential point or line, and these parameters
-- will include the terrestrial reference frame ("Ellipsoid" in this library) used as a foundation.
-- Hence conversion from a geodetic to a grid point requires the \"basis\" for the grid in question,
-- and grid points carry that basis with them because without it there is no defined relationship
-- between the grid points and terrestrial positions.
class GridClass r e | r->e where
   fromGrid :: GridPoint r -> Geodetic e
   toGrid :: r -> Geodetic e -> GridPoint r
   gridEllipsoid :: r -> e


-- | A point on the specified grid.
data GridPoint r = GridPoint {
   eastings, northings, altGP :: Double,
   gridBasis :: r
} deriving (Eq, Show)

instance HasAltitude (GridPoint g) where
   altitude = altGP
   setAltitude h gp = gp{altGP = h}


-- | A vector relative to a point on a grid. All distances are in meters.
-- Operations that use offsets will only give
-- meaningful results if all the points come from the same grid.
--
-- The monoid instance is the sum of offsets.
data GridOffset = GridOffset {
   deltaEast, deltaNorth, deltaAltitude :: Double
} deriving (Eq, Show)

instance Semigroup GridOffset where
  g1 <> g2 = GridOffset (deltaEast g1 + deltaEast g2)
                        (deltaNorth g1 + deltaNorth g2)
                        (deltaAltitude g1 + deltaAltitude g2)

instance Monoid GridOffset where
   mempty = GridOffset 0 0 0
   mappend = (<>)


-- | An offset defined by a distance (m) and a bearing (radians) to the right of North.
--
-- There is no elevation parameter because we are using a plane to approximate an ellipsoid,
-- so elevation would not provide a useful result.  If you want to work with elevations
-- then 'Geodetics.Path.rayPath' will give meaningful results.
polarOffset :: Double -> Double -> GridOffset
polarOffset r d = GridOffset (r * sin d) (r * cos d) 0


-- | Scale an offset by a scalar.
offsetScale :: Double -> GridOffset -> GridOffset
offsetScale s off = GridOffset (deltaEast off * s)
                               (deltaNorth off * s)
                               (deltaAltitude off * s)

-- | Invert an offset.
offsetNegate :: GridOffset -> GridOffset
offsetNegate off = GridOffset (negate $ deltaEast off)
                              (negate $ deltaNorth off)
                              (negate $ deltaAltitude off)


-- Add an offset on to a point to get another point.
applyOffset :: GridOffset -> GridPoint g -> GridPoint g
applyOffset off p = GridPoint (eastings p + deltaEast off)
                           (northings p + deltaNorth off)
                           (altitude p + deltaAltitude off)
                           (gridBasis p)


-- | The distance represented by an offset.
offsetDistance :: GridOffset -> Double
offsetDistance = sqrt . offsetDistanceSq


-- | The square of the distance represented by an offset.
offsetDistanceSq :: GridOffset -> Double
offsetDistanceSq off =
   deltaEast off ^ _2 + deltaNorth off ^ _2 + deltaAltitude off ^ _2


-- | The direction represented by an offset, as bearing to the right of North.
offsetBearing :: GridOffset -> Double
offsetBearing off = atan2 (deltaEast off) (deltaNorth off)


-- | The offset required to move from p1 to p2.
gridOffset :: GridPoint g -> GridPoint g -> GridOffset
gridOffset p1 p2 = GridOffset (eastings p2 - eastings p1)
                              (northings p2 - northings p1)
                              (altitude p2 - altitude p1)


-- | Coerce a grid point of one type into a grid point of a different type,
-- but with the same easting, northing and altitude. This is unsafe because it
-- will produce a different position unless the two grids are actually equal.
--
-- This function should be used only to convert between distinguished grids
-- (e.g. "UkNationalGrid") and their equivalent numerical definitions.
unsafeGridCoerce :: b -> GridPoint a -> GridPoint b
unsafeGridCoerce base p = GridPoint (eastings p) (northings p) (altitude p) base


-- | Convert a list of digits to a distance. The first argument is the size of the
-- grid square within which these digits specify a position. The first digit is
-- in units of one tenth of the grid square, the second one hundredth, and so on.
-- The first result is the lower limit of the result, and the second is the size
-- of the specified offset.
--
-- For instance @fromGridDigits (100 * kilometer) "237"@ will return
--
-- > Just (23700, 100)
--
-- If there are any non-digits in the string then the function returns @Nothing@.
fromGridDigits :: Double -> String -> Maybe (Double, Double)
fromGridDigits sq ds = if all isDigit ds then Just (d, p) else Nothing
   where
      n :: Integer
      n = fromIntegral $ length ds
      d = sum $ zipWith (*)
         (map (fromIntegral . digitToInt) ds)
         (drop 1 $ iterate (/ 10) sq)
      p = sq / fromIntegral ((10 :: Integer) ^ n)


-- | Convert a distance into a digit string suitable for printing as part
-- of a grid reference. The result is the south or west side of the enclosing grid square,
-- where the size of the square is defined by the number of digits.
-- The result is expressed as an integer count of squares and a string of digits.
-- If any arguments are invalid then @Nothing@ is returned.
toGridDigits ::
   Double           -- ^ Size of enclosing grid square. Must be at least 1000m.
   -> Int           -- ^ Number of digits to return. Must be positive.
   -> Double        -- ^ Offset to convert into grid (m).
   -> Maybe (Integer, String)
toGridDigits sq n d =
   if sq < 1000 || n < 0 || d < 0
   then Nothing
   else
      Just (sqs, pad)
   where
      p :: Integer
      p = 10 ^ n
      unit :: Double
      unit = sq / fromIntegral p
      u = floor (d / unit)
      (sqs, d1) = u `divMod` p
      s = show d1
      pad = if n == 0 then "" else replicate (n - length s) '0' ++ s


-- | Some grids (notably UTM and UPS) can have optional units in their grid references.
data GridUnit = GridMeters | GridKilometers deriving Eq

instance Show GridUnit where
   show GridMeters = "m"
   show GridKilometers = "km"
