{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies #-}

module Geodetics.Grid (
   -- ** Grid Types
   GridClass (..),
   GridPoint (..),
   GridOffset (..),
   -- ** Grid Operations
   polarOffset,
   offsetScale,
   offsetNegate,
   applyOffset,
   offsetDistance,
   offsetDistanceSq,
   offsetBearing,
   gridOffset,
   -- ** Unsafe Conversion
   unsafeGridCoerce
) where

import Data.Function
import Data.Monoid
import Geodetics.Altitude
import Geodetics.Geodetic
import Numeric.Units.Dimensional.Prelude
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
   altitude gp = altGP gp
   setAltitude h gp = gp{altGP = h}



-- | A vector relative to a point on a grid.
-- Operations that use offsets will only give
-- meaningful results if all the points come from the same grid.
data GridOffset = GridOffset {
   deltaEast, deltaNorth, deltaAltitude :: Length Double
} deriving (Eq, Show)

instance Monoid GridOffset where
   mempty = GridOffset (0 *~ meter) (0 *~ meter) (0 *~ meter)
   mappend g1 g2 = GridOffset (deltaEast g1 + deltaEast g2) 
                              (deltaNorth g1 + deltaNorth g2) 
                              (deltaAltitude g1 + deltaAltitude g2) 

-- | An offset defined by a distance and a bearing to the right of North.
--
-- There is no elevation parameter because this does not provide
-- useful results because we are using a plane to approximate a sphere,
-- and the 
polarOffset :: Length Double -> Dimensionless Double -> GridOffset
polarOffset r d = GridOffset (r * sin d) (r * cos d) (0 *~ meter)


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

              
-- | The direction represented by an offset, as radians to the right of North.
offsetBearing :: GridOffset -> Dimensionless Double
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
-- It should be used only to convert between priviledged grids (e.g. "UkNationalGrid") and
-- their equivalent numerical definitions.
unsafeGridCoerce :: b -> GridPoint a -> GridPoint b
unsafeGridCoerce base p = GridPoint (eastings p) (northings p) (altitude p) base
