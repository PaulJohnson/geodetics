{-# LANGUAGE TypeOperators, TypeFamilies, FlexibleContexts #-}
-- | The implementation assumes IEEE 754 arithmetic.

module Geodetics.Path where

import Control.Monad
import Geodetics.Ellipsoids
import Geodetics.Geodetic


-- | Lower and upper exclusive distance bounds within which a path is valid. 
type PathValidity = (Double, Double)

-- | A path is a parametric function of distance along the path. The result is the
-- position, and the direction of the path at that point as heading and elevation angles.
--
-- A well-behaved path must be continuous and monotonic (that is, if the distance increases
-- then the result is further along the path) for all distances within its validity range.
-- Ideally the physical distance along the path from the origin to the result should be equal
-- to the argument. If this is not practical then a reasonably close approximation may be used,
-- such as a spheroid instead of an ellipsoid, provided that this is documented. 
-- Outside its validity the path function may
-- return anything or bottom.
data Path e = Path {
      pathFunc :: Double -> (Geodetic e, Double, Double),
         -- ^ Takes a length and returns a position, and direction as heading and elevation angles.
      pathValidity :: PathValidity
   }
   
-- | Convenience value for paths that are valid for all distances.
alwaysValid :: PathValidity
alwaysValid = (negate inf, inf) where
   inf = 1.0 / 0  -- Assumes IEE arithmetic.


-- | True if the path is valid at that distance.
pathValidAt :: Path e -> Double -> Bool
pathValidAt path d = d > x1 && d < x2
   where (x1,x2) = pathValidity path


-- | Find where a path meets a given condition using bisection. This is not the
-- most efficient algorithm, but it will always produce a result if there is one
-- within the initial bounds. If there is more than one result then an arbitrary
-- one will be returned.
-- 
-- The initial bounds must return one GT or EQ value and one LT or EQ value. If they
-- do not then @Nothing@ is returned.
bisect :: 
   Path e 
   -> (Geodetic e -> Ordering)        -- ^ Evaluation function.
   -> Double                          -- ^ Required accuracy in terms of distance along the path.
   -> Double -> Double                -- ^ Initial bounds.
   -> Maybe Double
bisect path f t b1 b2 = do
      guard $ pathValidAt path b1
      guard $ pathValidAt path b2
      let r = pairResults b1 b2
      guard $ hasRoot r
      bisect1 $ sortPair r
   where
      f' d = let (p, _, _) = pathFunc path d in f p 
      pairResults d1 d2 = ((d1, f' d1), (d2, f' d2))
      hasRoot (v1, v2) = snd v1 <= EQ && EQ <= snd v2
      sortPair (v1, v2) = if snd v1 <= snd v2 then (v1, v2) else (v2, v1)
      bisect1 ((d1, r1), (d2, r2)) =
         let d3 = (d1 + d2) / 2
             r3 = f' d3
             c1 = ((d1, r1), (d3, r3))
             c2 = ((d3, r3), (d2, r2))
         in if abs (d1 - d2) <= t 
            then return d3
            else bisect1 $ if hasRoot c1 then c1 else c2


-- | Try to find the intersection point of two ground paths (i.e. ignoring 
-- altitude). Returns the distance of the intersection point along each path
-- using a modified Newton-Raphson method. If the two paths are 
-- fairly straight and not close to parallel then this will converge rapidly.
--
-- The algorithm projects great-circle paths forwards using the bearing at the 
-- estimate to find the estimated intersection, and then uses the distances to this 
-- intersection as the next estimates.
--
-- If either estimate departs from its path validity then @Nothing@ is returned.
intersect :: (Ellipsoid e) =>
   Double -> Double                   -- ^ Starting estimates.
   -> Double                          -- ^ Required accuracy.
   -> Int                             -- ^ Iteration limit. Returns @Nothing@ if this is reached.  
   -> Path e -> Path e                -- ^ Paths to intersect.
   -> Maybe (Double, Double)
intersect d1 d2 accuracy n path1 path2
   | not $ pathValidAt path1 d1     = Nothing
   | not $ pathValidAt path2 d2     = Nothing
   | n <= 0                         = Nothing
   | mag < 1e-15                    = Nothing
   | mag3 (nv1 `cross3` nv2) * r <= accuracy = Just (d1, d2)
       -- Assumes that sin (accuracy/r) == accuracy/r
   | otherwise = 
      if abs d1a + abs d2a < abs d1b + abs d2b
         then intersect (d1 + d1a) (d2 + d2a) accuracy (pred n) path1 path2
         else intersect (d1 + d1b) (d2 + d2b) accuracy (pred n) path1 path2
   where
      (pt1, h1, _) = pathFunc path1 d1
      (pt2, h2, _) = pathFunc path2 d2
      vectors :: Double -> Double -> Double 
                 -> (Vec3 Double, Vec3 Double)
      vectors lat lon b = (
          -- Unit vector of normal to surface at (lat,lon)
         (cosLat*cosLon, cosLat*sinLon, sinLat),
         -- Normal of great circle defined by bearing b at (lat,lon)
         (sinLon * cosB - sinLat * cosLon * sinB,
          negate cosLon * cosB - sinLat * sinLon * sinB,
           cosLat * sinB))
         where
            sinLon = sin lon
            sinLat = sin lat
            cosLon = cos lon
            cosLat = cos lat
            sinB = sin b
            cosB = cos b
      mag3 (x,y,z) = sqrt $ x*x + y*y + z*z
      (nv1, gc1) = vectors (latitude pt1) (longitude pt1) h1
      (nv2, gc2) = vectors (latitude pt2) (longitude pt2) h2
      nv3 = gc1 `cross3` gc2         -- Intersection of the great circles
      mag = mag3 nv3
      nv3a = scale3 nv3 (1 / mag)   -- Scale to unit. See outer function for case when mag3 == 0
      nv3b = negate3 nv3a            -- Antipodal result. Take the closest.
      -- Find "nearest" intersection, defined as smaller of sum of distances to current points.
      d1a = gcDist gc1 nv1 nv3a * r
      d2a = gcDist gc2 nv2 nv3a * r
      d1b = gcDist gc1 nv1 nv3b * r
      d2b = gcDist gc2 nv2 nv3b * r
      -- Signed angle between v1 and v2, 
      gcDist norm v1 v2 = 
         let c = v1 `cross3` v2 
         in (if c `dot3` norm < 0 then negate else id) $ atan2 (mag3 c) (v1 `dot3` v2) 
      r = majorRadius $ ellipsoid pt1
          
{- Note on derivation of "intersect"

The algorithm is a variant of the Newton-Raphson method, and shares its advantage
of rapid convergence in many useful cases. Each path has a current approximation to
the intersection, and the next approximation is computed by projecting both paths 
along great circles from the current approximation and finding the point where those
great circles intersect. A spherical Earth is assumed for simplicity. 

The Great Circle calculations use a vector method rather than spherical trigonometry.
This avoids a lot of transcendental functions and also the singularities inherent in 
polar coordinate systems. This implementation is based on formulae from 
<http://www.movable-type.co.uk/scripts/latlong-vectors.html>, which in turn is based on 
"A Non-singular Horizontal Position Representation" by Kenneth Gade, THE JOURNAL OF 
NAVIGATION (2010), 63, 395–417.
<http://www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf>


"pt1" is the current approximation for the result on "path1".  The vector "nv1" is the 
unit vector pointing "pt1". "gc1" is the normal to the great circle tangent to 
"path1" at "pt1". "nv2" and "gc2" are similarly derived from "path2".

The intersection of the planes of "gc1" and "gc2" is given by their cross product, "nv3".
This is scaled to a unit vector "nv3a", and the other intersection of the great circles is
at "nv3b", opposite. The distances from nv1 and nv2 to nv3a and nv3b are computed using
the corresponding great-circle normals to determine whether the distances are positive or 
negative along the paths. The nearest solution (defined in terms of great-circle distance) 
is taken as the basis for the next approximation.
-}

-- | A ray from a point heading in a straight line in 3 dimensions. 
rayPath :: (Ellipsoid e) => 
   Geodetic e          -- ^ Start point.
   -> Double           -- ^ Bearing.
   -> Double           -- ^ Elevation.
   -> Path e
rayPath pt1 bearing elevation = Path ray alwaysValid
   where
      ray distance = (Geodetic lat long alt (ellipsoid pt1), bearing2, elevation2)
         where
            pt2' = pt1' `add3` (delta `scale3` distance)      -- ECEF of result point.
            (lat,long,alt) = earthToGeo (ellipsoid pt1) pt2'  -- Geodetic of result point.
            (dE,dN,dU) = transform3 (trans3 $ ecefMatrix lat long) delta  -- Direction of ray at result point.
            elevation2 = asin dU
            bearing2 = if dE == 0 && dN == 0 then bearing else atan2 dE dN  -- Allow for vertical elevation.
            
      ecefMatrix lat long =   -- Transform matrix for vectors from (East, North, Up) to (X,Y,Z).
         ((negate sinLong, negate cosLong*sinLat, cosLong*cosLat),
              --    East X      North X               Up X
          (       cosLong, negate sinLong*sinLat, sinLong*cosLat),
              --    East Y      North Y               Up Y
          (             0,      cosLat         , sinLat))
              --    East Z      North Z               Up Z
         where
            sinLong = sin long
            cosLong = cos long
            sinLat = sin lat
            cosLat = cos lat
      
      direction = (sinB*cosE, cosB*cosE, sinE)  -- Direction of ray in ENU
      delta = transform3 (ecefMatrix (latitude pt1) (longitude pt1)) direction  -- Convert to ECEF
      pt1' = geoToEarth pt1    -- ECEF of origin point.
      sinB = sin bearing
      cosB = cos bearing
      sinE = sin elevation
      cosE = cos elevation
      
-- | Rhumb line: path following a constant course. Also known as a loxodrome.
--
-- The valid range stops a few arc-minutes short of the poles to ensure that the 
-- polar singularities are not included. Anyone using a rhumb line that close to a pole
-- must be going round the twist anyway.
--
-- Based on /Practical Sailing Formulas for Rhumb-Line Tracks on an Oblate Earth/
-- by G.H. Kaplan, U.S. Naval Observatory. Except for points close to the poles 
-- the approximation is accurate to within a few meters over 1000km.
rhumbPath :: (Ellipsoid e) =>
   Geodetic e            -- ^ Start point.
   -> Double             -- ^ Course.
   -> Path e
rhumbPath pt course = Path rhumb validity
   where
      rhumb distance = (Geodetic lat (properAngle lon) 0 (ellipsoid pt), course, 0)
         where
            lat' = lat0 + distance * cosC / m0   -- Kaplan Eq 13.
            lat = lat0 + (m0 / (a*(1-e2))) * ((1-3*e2/4)*(lat'-lat0)
                                            + (3*e2/8)*(sin (2*lat') - sin (2*lat0)))
            lon | abs cosC > 1e-7
                     = lon0 + tanC * (q lat - q0)     -- Kaplan Eq 16.
                | otherwise
                     = lon0 + distance * sinC / latitudeRadius (ellipsoid pt) ((lat0 + lat')/2)
      validity
         | cosC > 0  = ((negate pi/2 - latitude pt) * b / cosC, (pi/2 - latitude pt) * b / cosC)
         | otherwise  = ((pi/2 - latitude pt) * b / cosC, (negate pi/2 - latitude pt) * b / cosC)
      q0 = q lat0
      q phi = log (tan (pi/4+phi/2)) + e * log ((1-eSinPhi)/(1+eSinPhi)) / 2
         where                                -- Factor out expression from Eq 16 of Kaplan
            eSinPhi = e * sin phi
      sinC = sin course
      cosC = cos course
      tanC = tan course
      lat0 = latitude pt
      lon0 = longitude pt
      e2 = eccentricity2 $ ellipsoid pt
      e = sqrt e2
      m0 = meridianRadius (ellipsoid pt) lat0
      a = majorRadius $ ellipsoid pt
      b = minorRadius $ ellipsoid pt
   

-- | A path following the line of latitude around the Earth eastwards.
--
-- This is equivalent to @rhumbPath pt (pi/2)@
latitudePath :: (Ellipsoid e) =>
   Geodetic e           -- ^ Start point.
   -> Path e
latitudePath pt = Path line alwaysValid
   where
      line distance = (pt2, pi/2, 0)
         where
            pt2 = Geodetic 
               (latitude pt) (longitude pt + distance / r)
               0 (ellipsoid pt)
      r = latitudeRadius (ellipsoid pt) (latitude pt)


-- | A path from the specified point to the North Pole. Use negative distances
-- for the southward path.
--
-- This is equivalent to @rhumbPath pt _0@
longitudePath :: (Ellipsoid e) =>
   Geodetic e    -- ^ Start point.
   -> Path e
longitudePath pt = rhumbPath pt 0
