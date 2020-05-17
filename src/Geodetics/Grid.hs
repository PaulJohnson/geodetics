{-# LANGUAGE FunctionalDependencies #-}

module Geodetics.Grid (
   -- ** Grid types
   GridClass (..),
   GridPoint (..),
   GridOffset (..),
   -- ** Grid operations
   polarOffset,
   offsetScale,
   offsetNegate,
   applyOffset,
   offsetDistance,
   offsetDistanceSq,
   offsetBearing,
   gridOffset,
   -- ** Unsafe conversion
   unsafeGridCoerce,
   -- ** Utility functions for grid references
   fromGridDigits,
   toGridDigits
) where

import Data.Char
import Data.Function
import Geodetics.Altitude
import Geodetics.Geodetic
import Numeric.Units.Dimensional.Prelude hiding ((.))
import qualified Prelude as P

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
   eastings, northings, altGP :: Length Double,
   gridBasis :: r
} deriving (Show)


instance Eq (GridPoint r) where
   p1 == p2  =
      eastings p1 == eastings p2 &&
      northings p1 == northings p2 &&
      altGP p1 == altGP p2

instance HasAltitude (GridPoint g) where
   altitude = altGP
   setAltitude h gp = gp{altGP = h}



-- | A vector relative to a point on a grid.
-- Operations that use offsets will only give
-- meaningful results if all the points come from the same grid.
--
-- The monoid instance is the sum of offsets.
data GridOffset = GridOffset {
   deltaEast, deltaNorth, deltaAltitude :: Length Double
} deriving (Eq, Show)

instance Semigroup GridOffset where
  g1 <> g2 = GridOffset (deltaEast g1 + deltaEast g2)
                        (deltaNorth g1 + deltaNorth g2)
                        (deltaAltitude g1 + deltaAltitude g2)

instance Monoid GridOffset where
   mempty = GridOffset _0 _0 _0
   mappend = (<>)

-- | An offset defined by a distance and a bearing to the right of North.
--
-- There is no elevation parameter because we are using a plane to approximate an ellipsoid,
-- so elevation would not provide a useful result.  If you want to work with elevations
-- then "rayPath" will give meaningful results.
polarOffset :: Length Double -> Angle Double -> GridOffset
polarOffset r d = GridOffset (r * sin d) (r * cos d) _0


-- | Scale an offset by a scalar.
offsetScale :: Dimensionless Double -> GridOffset -> GridOffset
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
offsetDistance :: GridOffset -> Length Double
offsetDistance = sqrt . offsetDistanceSq


-- | The square of the distance represented by an offset.
offsetDistanceSq :: GridOffset -> Area Double
offsetDistanceSq off =
   deltaEast off ^ pos2 + deltaNorth off ^ pos2 + deltaAltitude off ^ pos2


-- | The direction represented by an offset, as bearing to the right of North.
offsetBearing :: GridOffset -> Angle Double
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
-- It should be used only to convert between distinguished grids (e.g. "UkNationalGrid") and
-- their equivalent numerical definitions.
unsafeGridCoerce :: b -> GridPoint a -> GridPoint b
unsafeGridCoerce base p = GridPoint (eastings p) (northings p) (altitude p) base



-- | Convert a list of digits to a distance. The first argument is the size of the
-- grid square within which these digits specify a position. The first digit is
-- in units of one tenth of the grid square, the second one hundredth, and so on.
-- The first result is the lower limit of the result, and the second is the size
-- of the specified offset.
-- So for instance @fromGridDigits (100 *~ kilo meter) "237"@ will return
--
-- > Just (23700 meters, 100 meters)
--
-- If there are any non-digits in the string then the function returns @Nothing@.
fromGridDigits :: Length Double -> String -> Maybe (Length Double, Length Double)
fromGridDigits sq ds = if all isDigit ds then Just (d, p) else Nothing
   where
      n = length ds
      d = sum $ zipWith (*)
         (map ((*~ one) . fromIntegral . digitToInt) ds)
         (tail $ iterate (/ (10 *~ one)) sq)
      p = sq / ((10 *~ one) ** (fromIntegral n *~ one))

-- | Convert a distance into a digit string suitable for printing as part
-- of a grid reference. The result is the nearest position to the specified
-- number of digits, expressed as an integer count of squares and a string of digits.
-- If any arguments are invalid then @Nothing@ is returned.
toGridDigits ::
   Length Double    -- ^ Size of enclosing grid square. Must be at least 1 km.
   -> Int           -- ^ Number of digits to return. Must be positive.
   -> Length Double -- ^ Offset to convert into grid.
   -> Maybe (Integer, String)
toGridDigits sq n d =
   if sq < (1 *~ kilo meter) || n < 0 || d < _0
   then Nothing
   else
      Just (sqs, pad)
   where
      p :: Integer
      p = 10 P.^ n
      unit :: Length Double
      unit = sq / (fromIntegral p *~ one)
      u = round ((d / unit) /~ one)
      (sqs, d1) = u `divMod` p
      s = show d1
      pad = if n == 0 then "" else replicate (n P.- length s) '0' ++ s
