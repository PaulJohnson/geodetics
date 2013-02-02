{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleInstances #-}

{- |
This module defines coordinate systems for positions on planet Earth and conversions between them. 

The question \"where am I?\" on the Earth is not as simple as it sounds, and the greater the degree of precision
required, the more difficult the question becomes. An accessible introduction to the topic can be found at
http://www.ordnancesurvey.co.uk/oswebsite/gps/docs/A_Guide_to_Coordinate_Systems_in_Great_Britain.pdf.

This library is designed to give answers accurate to within a few meters, apart from conversion between 
local ellipsoids, which can be a few tens of meters out. 


This library does not take account of the difference between an ellipsoid and the actual geoid. Instead
vertical position is treated as \"altitude\"; that is, the height above some known datum (such as the 
Ordnance Data Newlyn, which is the official measurement of \"mean sea level\" in the UK). All calculations
assume that this is the ellipsoid height, so when converting between ellipsoids the altitude value
is held constant and only the latitude and longitude change.

For more information on the geoid, see the Guide referenced above or http://en.wikipedia.org/wiki/Geoid.

Despite the number of getter and setter functions this library does not use lenses; there is no accepted
standard lens library at present, so it is better if a client wishing to use one creates the relevant
instances.
-}


module Geodetics.Coordinates (
   -- ** Geodetic Coordinates
   HasAltitude (..),
   Geodetic (..),
   toLocal,
   toWGS84,
   antipodes,
   geometricalDistance,
   geometricalDistanceSq,
   -- ** General Grid Coordinates
   GridClass (..),
   GridPoint (..),
   GridOffset (..),
   polarOffset,
   offsetScale,
   offsetNegate,
   applyOffset,
   offsetDistance,
   offsetDistanceSq,
   offsetBearing,
   gridOffset,
   -- ** Transverse Mercator Grid
   GridTM (trueOrigin, falseOrigin, gridScale),
   mkGridTM
) where

import Data.Char (chr)
import Data.Function
import Data.Monoid
import Geodetics.Ellipsoids
import Numeric
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P

import Debug.Trace   -- Uncomment to enable tracing.


-- | All geographical coordinate systems need the concept of altitude above a point on the ground.
-- 
-- Minimum definition: altitude, setAltitude.
class HasAltitude a where
   altitude :: a -> Length Double
   setAltitude :: Length Double -> a -> a
   groundPosition :: a -> a  -- ^ Set altitude to zero.
   groundPosition = setAltitude (0 *~ meter)
   
   
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
-- use "geometricDistance" (or its squared variant).
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
antipodes :: (Ellipsoid e) => Geodetic e -> Geodetic e
antipodes g = Geodetic lat long (geoAlt g) (ellipsoid g)
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

-- ---------------------------


-- | A Grid is a two-dimensional projection of the ellipsoid onto a plane.
class GridClass r e | r->e where
   fromGrid :: GridPoint r -> Geodetic e
   toGrid :: r -> Geodetic e -> GridPoint r
   gridEllipsoid :: r -> e


-- | A point on the specified grid. 
data GridPoint r = GridPoint {
   eastings, northings, altGP :: Length Double,
   gridBase :: r
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
                           (gridBase p)


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





-- | A Transverse Mercator projection gives an approximate mapping of the ellipsoid on to a 2-D grid. It models
-- a sheet curved around the ellipsoid so that it touches it at one north-south line (hence making it part of
-- a slightly elliptical cylinder).
data GridTM e = GridTM {
   trueOrigin :: Geodetic e,
      -- ^ A point on the line where the projection touches the ellipsoid (altitude is ignored).
   falseOrigin :: GridOffset,
      -- ^ An offset from the true origin expressed as a vector from the true origin, 
      -- used to avoid negative coordinates over the area of interest.
      -- The altitude gives a vertical offset from the ellipsoid.
   gridScale :: Dimensionless Double,
      -- ^ A scaling factor that balances the distortion between the east & west edges and the middle 
      -- of the projection.
      
   -- Remaining elements are memoised parameters computed from the ellipsoid underlying the true origin.
   grid_n :: Dimensionless Double,
   gridN1, gridN2, gridN3, gridN4 :: Dimensionless Double
   -- grid_A :: Length Double,
   -- grid_alpha :: [Dimensionless Double],
   -- grid_beta :: [Dimensionless Double],
   -- grid_delta :: [Dimensionless Double]
} deriving (Show)



mkGridTM :: (Ellipsoid e) => Geodetic e -> GridOffset -> Dimensionless Double -> GridTM e
mkGridTM origin offset sf =
   GridTM {trueOrigin = origin,
           falseOrigin = offset,
           gridScale = sf,
           grid_n = n,
           gridN1 = _1 + n + (_5/_4) * n^pos2 + (_5/_4) * n^pos3,
           gridN2 = _3 * n + _3 * n^pos2 + ((21*~one)/_8) * n^pos3,
           gridN3 = ((15*~one)/_8) * (n^pos2 + n^pos3),
           gridN4 = ((35*~one)/(24*~one)) * n^pos3
        }
    where 
       f = flattening $ ellipsoid origin
       n = f / (_2-f)  -- Equivalent to (a-b)/(a+b) where b = (1-f)*a



-- [1]: "A Guide to Coordinate Systems in Great Britain", 
--      v2.1, Oct 2010, by the Ordnance Survey.

-- | Equation C3 from reference [1].
m :: (Ellipsoid e) => GridTM e -> Dimensionless Double -> Length Double
m grid lat = bF0 * (gridN1 grid * dLat 
                    - gridN2 grid * sin dLat * cos sLat
                    + gridN3 grid * sin (_2 * dLat) * cos (_2 * sLat) 
                    - gridN4 grid * sin (_3 * dLat) * cos (_3 * sLat))
   where
      dLat = lat - latitude (trueOrigin grid)
      sLat = lat + latitude (trueOrigin grid)
      bF0 = (minorRadius $ gridEllipsoid grid) * gridScale grid


instance (Ellipsoid e) => GridClass (GridTM e) e where
   fromGrid p = Geodetic
      (lat' - east' ^ pos2 * tanLat / (_2 * rho * v)  -- Term VII
            + east' ^ pos4 * (tanLat / ((24 *~ one) * rho * v ^ pos3)) 
                           * (_5 + _3 * tanLat ^ pos2 + eta2 - _9 * tanLat ^ pos2 * eta2)  -- Term VIII
            - east' * east' ^ pos5 * (tanLat / ((720 *~ one) * rho * v ^ pos5))
                           * (61 *~ one + (90 *~ one) * tanLat ^ pos2 + (45 *~ one) * tanLat ^ pos4)) -- Term IX
      (longitude (trueOrigin grid) 
            + east' / (cosLat * v)  -- Term X
            - (east' ^ pos3 / (_6 * cosLat * v ^ pos3)) * (v / rho + _2 * tanLat ^ pos2)  -- Term XI
            + (east' ^ pos5 / ((120 *~ one) * cosLat * v ^ pos5)) 
                 * (_5 + (28 *~ one) * tanLat ^ pos2  + (24 *~ one) * tanLat ^ pos4)  -- Term XII
            - (east' ^ pos5 * east' ^ pos2 / ((5040 *~ one) * cosLat * v * v * v ^ pos5))
                 * ((61 *~ one) + (662 *~ one) * tanLat ^ pos2 + (1320 *~ one) * tanLat ^ pos4 + (720 *~ one) * tanLat * tanLat ^ pos5)) -- Term XIIa
     (0 *~ meter) (gridEllipsoid grid)
            
            
      where
         GridPoint east' north' _ _ = (offsetNegate $ falseOrigin grid) `applyOffset` p
         lat' = fst $ head $ dropWhile ((> 0.01 *~ milli meter) . snd) 
               $ tail $ iterate next (latitude $ trueOrigin grid, 1 *~ meter) 
            where
               next (phi, _) = let delta = north' - m grid phi in (phi + delta / aF0, delta) 
               -- head and tail are safe because iterate returns an infinite list.
          
         sinLat = sin lat'
         cosLat = cos lat'
         tanLat = tan lat'
         sinLat2 = sinLat ^ pos2
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
               
               
         aF0 = (majorRadius $ gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         grid = gridBase p
         
   toGrid grid geo = trace traceMsg $ applyOffset (off  `mappend` falseOrigin grid) $ 
                     GridPoint (0 *~ metre) (0 *~ metre) (0 *~ metre) grid
      where
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
         off = GridOffset
                  (dLong * term_IV
                   + dLong ^ pos3 * term_V
                   + dLong ^ pos5 * term_VI)
                  (m grid lat + dLong ^ pos2 * term_II
                     + dLong ^ pos4 * term_III 
                     + dLong * dLong ^ pos5 * term_IIIa)
                  (0 *~ meter)
         -- Terms defined in [1].
         term_II   = (v/_2) * sinLat * cosLat
         term_III  = (v/(24*~one)) * sinLat * cosLat ^ pos3 
                     * (_5 - tanLat ^ pos2 + _9 * eta2)
         term_IIIa = (v/(720*~one)) * sinLat * cosLat ^ pos5 
                     * ((61 *~ one) - (58 *~ one) * tanLat ^ pos2 + tanLat ^ pos4)
         term_IV   = v * cosLat
         term_V    = (v/_6) * cosLat ^ pos3 * (v/rho - tanLat ^ pos2)
         term_VI   = (v/(120*~one)) * cosLat ^ pos5 
                     * (_5 - (18*~one) * tanLat ^ pos2 
                              + tanLat ^ pos4 + (14*~one) * eta2
                              - (58*~one) * tanLat ^ pos2 * eta2)
         -- Trace message for debugging. Uncomment this code if necessary
         traceMsg = concat [
            "v    = ", show v, "\n",
            "rho  = ", show rho, "\n",
            "eta2 = ", show eta2, "\n",
            "M    = ", show $ m grid lat, "\n",
            "I    = ", show $ m grid lat + deltaNorth (falseOrigin grid), "\n",
            "II   = ", show term_II, "\n",
            "III  = ", show term_III, "\n",
            "IIIa = ", show term_IIIa, "\n",
            "IV   = ", show term_IV, "\n",
            "V    = ", show term_V, "\n",
            "VI   = ", show term_VI, "\n"]
         -- Common subexpressions
         lat = latitude geo
         long = longitude geo
         dLong = long - longitude (trueOrigin grid)
         sinLat = sin lat
         cosLat = cos lat
         tanLat = tan lat
         sinLat2 = sinLat ^ pos2
         aF0 = (majorRadius $ gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         
   gridEllipsoid = ellipsoid . trueOrigin

{- This code does not work. Rather than debug it I'll try using the OSGB formulae.

mkGridTM :: (Ellipsoid e) => Geodetic e -> GridOffset -> Dimensionless Double -> GridTM e
mkGridTM origin offset sf =
   GridTM {trueOrigin = origin,
           falseOrigin = offset,
           gridScale = sf,
           grid_n = n,
           grid_A = (a / (_1+n)) * (_1 + n ^ pos2 / _4 + n ^ pos4 / (64 *~ one)),
           grid_alpha = [polyN [240, (-320),    150], 
                         polyN [0,      130, -(288)], 
                         polyN [0,        0,    122]],
           grid_beta =  [polyN [240, (-320),    185],
                         polyN [0,       10,     32],
                         polyN [0,        0,     17]],
           grid_delta = [polyN [960, (-320), (-960)],
                         polyN [0,     1120, (-768)],
                         polyN [0,        0,   1792]]
          }
   where
      -- The magic numbers and formulae used here and in the Grid instance are taken from 
      -- http://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#Simplified_formulas
      -- retrieved 22 Nov 2012.
      polyN = sum . zipWith (*) nSeries . map ((/ (480 *~ one)) . (*~ one))
         -- Compute a polynomial of 'n'.
      n = f / (_2 - f)
      nSeries = iterate (n*) n
      a = majorRadius $ ellipsoid origin
      f = flattening $ ellipsoid origin
      
      
instance (Ellipsoid e) => GridClass (GridTM e) e where
   fromGrid p = Geodetic lat long alt (gridEllipsoid $ gridBase p) 
      where
         trueP = applyOffset (offsetNegate $ falseOrigin $ gridBase p) p
         lat = chi + sum (zipWith (*) delta (map (sin . (chi *)) [_2,_4..]))
         long = longitude (trueOrigin $ gridBase p) + atan2 (sinh eta') (cos xi')
         alt = altitude trueP
         xi = northings trueP / scaleA
         eta = eastings trueP / scaleA
         xi' = xi - sum (zipWith3 prod3 beta (map (sin . (xi *)) [_2,_4..]) (map (cosh . (eta *)) [_2,_4..]))
         eta' = eta - sum (zipWith3 prod3 beta (map (cos . (xi *)) [_2,_4..]) (map (sinh . (eta *)) [_2,_4..]))
         prod3 x1 x2 x3 = x1 * x2 * x3 
         chi = asin $ sin xi' / cosh eta'
         scaleA = gridScale (gridBase p) * grid_A (gridBase p)
         beta = grid_beta $ gridBase p
         delta = grid_delta $ gridBase p
   toGrid grid geo = applyOffset (off  `mappend` falseOrigin grid) $ 
                     GridPoint (0 *~ metre) (0 *~ metre) (0 *~ metre) grid
      where
         off = GridOffset 
            (scaleA * (eta' + sum (zipWith3 prod3 alpha (map (cos . (xi' *)) [_2,_4..]) (map (sinh . (eta' *)) [_2,_4..]))))
            (scaleA * (xi'  + sum (zipWith3 prod3 alpha (map (sin . (xi' *)) [_2,_4..]) (map (cosh . (eta' *)) [_2,_4..]))))
            (altitude geo)
         t = sinh $ atanh sinLat - n' * atanh (n' * sinLat)
         xi' = atan $ t / cos long
         eta' = atanh $ sin long / sqrt (_1+t^pos2)
         prod3 x1 x2 x3 = x1 * x2 * x3
         n' = _2 * sqrt n / (_1+n)
         n = grid_n grid
         sinLat = sin lat
         scaleA = gridScale grid * grid_A grid
         alpha = grid_alpha grid
         lat = latitude geo
         long = longitude geo - longitude (trueOrigin grid)
   gridEllipsoid = ellipsoid . trueOrigin
   
   -}