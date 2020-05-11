{-# OPTIONS_GHC -fno-warn-orphans #-}

module Main where

import Data.Maybe
import Data.Monoid
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P
import Test.Framework (Test, defaultMainWithOpts, testGroup)
import Test.Framework.Options (TestOptions, TestOptions'(..))
import Test.Framework.Runners.Options (RunnerOptions, RunnerOptions'(..))
import Test.Framework.Providers.HUnit
import Test.Framework.Providers.QuickCheck2 (testProperty)
import qualified Test.HUnit as HU
import Test.QuickCheck
import Test.QuickCheck.Checkers (EqProp, eq, (=-=), unbatch)
import Test.QuickCheck.Classes (monoid)

import ArbitraryInstances
import Geodetics.Altitude
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.Path
import Geodetics.Stereographic
import Geodetics.TransverseMercator
import Geodetics.UK


main :: IO ()
main = do
   let empty_test_opts = mempty :: TestOptions
   let my_test_opts = empty_test_opts {
     topt_maximum_generated_tests = Just 1000
   }

   let empty_runner_opts = mempty :: RunnerOptions
   let my_runner_opts = empty_runner_opts {
     ropt_test_options = Just my_test_opts
   }

   defaultMainWithOpts tests my_runner_opts

instance EqProp GridOffset where
  (GridOffset a b c) =-= (GridOffset a' b' c') =
    eq True $ a ≈ a' && b ≈ b' && c ≈ c'
    where x ≈ y = abs (x - y) < 0.00001 *~ meter

instance EqProp Helmert where
  (Helmert cX' cY' cZ' s rX' rY' rZ') =-= (Helmert cX'' cY'' cZ'' s' rX'' rY'' rZ'') =
    eq True $ and [cX' ≈ cX'', cY' ≈ cY'', cZ' ≈ cZ'',
                   s ≈- s',
                   rX' ≈- rX'', rY' ≈- rY'', rZ' ≈- rZ'']

    where x ≈ y = abs (x - y) < 0.00001 *~ meter
          x ≈- y = abs (x - y) < (_1 / (_5 * _2) ** (_5))

tests :: [Test]
tests = [
   testGroup "Geodetic" [
      testProperty "WGS84 and back" prop_WGS84_and_back,
      testProperty "Zero ground distance" prop_zero_ground,
      testGroup "UK Points" $ map pointTest ukPoints],
      testGroup "World lines" $ map worldLineTests worldLines,
   testGroup "Grid" [
      testProperty "Grid Offset 1" prop_offset1,
      testProperty "Grid Offset 2" prop_offset2,
      testProperty "Grid Offset 3" prop_offset3,
      testProperty "Grid 1" prop_grid1 ],
   testGroup "TransverseMercator" [
      testCase "fromGrid . toGrid == id" $ HU.assertBool "" prop_tmGridInverse
      ],
   testGroup "UK" [
      testProperty "UK Grid 1" prop_ukGrid1,
      testGroup "UK Grid 2" $ map ukGridTest2 ukSampleGrid,
      testGroup "UK Grid 3" $ map ukGridTest3 ukSampleGrid,
      testGroup "UK Grid 4" $ map ukGridTest4 ukSampleGrid,
      testGroup "UK Grid 5" $ map ukGridTest5 ukSampleGrid
      ],
   testGroup "Stereographic" [
      testCase "toGrid north" $ HU.assertBool "" stereographicToGridN,
      testCase "fromGrid north" $ HU.assertBool "" stereographicFromGridN,
      testCase "toGrid south" $ HU.assertBool "" stereographicToGridS,
      testCase "fromGrid south" $ HU.assertBool "" stereographicFromGridS,
      testProperty "Stereographic round trip" prop_stereographic
      ],
   testGroup "Paths" [
      testProperty "Ray Path 1" prop_rayPath1,
      testProperty "Ray Continuity" prop_rayContinuity,
      testProperty "Ray Bisection" prop_rayBisect,
      testProperty "Rhumb Continuity" prop_rhumbContinuity,
      testProperty "Rhumb Intersection" prop_rhumbIntersect
      ],
   testGroup "GridOffset" $ map (uncurry testProperty) $ unbatch $ monoid (mempty :: GridOffset),
   testGroup "Helmert" $ map (uncurry testProperty) $ unbatch $ monoid (mempty :: Helmert)
   ]


-- | The positions are within 30 cm.
samePlace :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
samePlace p1 p2 = geometricalDistance p1 p2 < 0.3 *~ meter


-- | The positions are within 10 m.
closeEnough :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
closeEnough p1 p2 = geometricalDistance p1 p2 < 10 *~ meter


-- | The angles are within 0.01 arcsec
sameAngle :: Angle Double -> Angle Double -> Bool
sameAngle v1 v2 = abs (properAngle (v1 - v2)) < 0.01 *~ arcsecond

-- | The grid positions are within 1mm
sameGrid :: (GridClass r e) => GridPoint r -> GridPoint r -> Bool
sameGrid p1 p2 = check eastings && check northings && check altitude
   where check f = f p1 - f p2 < 1 *~ milli meter


-- | Grid offsets are within 1mm.
sameOffset :: GridOffset -> GridOffset -> Bool
sameOffset go1 go2 = check deltaNorth && check deltaEast && check deltaAltitude
   where check f = f go1 - f go2 < 1 *~ milli meter


-- | The grid X and Y are both within 1 meter
closeGrid :: (GridClass r e) => GridPoint r -> GridPoint r -> Bool
closeGrid p1 p2 = check eastings && check northings && check altitude
   where check f = f p1 - f p2 < 1 *~ meter

-- | Degrees, minutes and seconds into radians.
dms :: Int -> Int -> Double -> Dimensionless Double
dms d m s = fromIntegral d *~ degree + fromIntegral m *~ arcminute + s *~ arcsecond

-- | Round-trip from local to WGS84 and back is identity (approximately)
prop_WGS84_and_back :: Geodetic LocalEllipsoid -> Bool
prop_WGS84_and_back p = samePlace p $ toLocal (ellipsoid p) $ toWGS84 p


-- | Test that for all points p, the ground distance from p to p is zero.
prop_zero_ground :: Geodetic WGS84 -> Bool
prop_zero_ground p =
   case groundDistance p p of
      Nothing -> False
      Just (d, _, _) -> abs d < 1 *~ milli meter


-- | Sample pairs of points with bearings and distances.
-- The Oracle for these values is the @FORWARD@ program from
--  <http://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html>
worldLines :: [(String, Geodetic WGS84, Geodetic WGS84, Length Double, Dimensionless Double, Dimensionless Double)]
worldLines = [
   ("Ordinary", Geodetic (40*~degree) (30*~degree) _0 WGS84, Geodetic (30*~degree) (50*~degree) _0 WGS84,
      2128852.999*~meter, 115.19596706*~degree, 126.79044315*~degree),
   ("Over Pole", Geodetic (60*~degree) (0*~degree) _0 WGS84, Geodetic (60*~degree) (180*~degree) _0 WGS84,
      6695785.820*~meter, 0*~degree, 180*~degree),
   ("Equator to Pole", Geodetic (0*~degree) (0*~degree) _0 WGS84, Geodetic (90*~degree) (180*~degree) _0 WGS84,
      10001965.729*~meter, 0*~degree, 180*~degree)]


worldLineTests :: (String, Geodetic WGS84, Geodetic WGS84, Length Double, Dimensionless Double, Dimensionless Double) -> Test
worldLineTests (str, g1, g2, d, a, b) = testCase str $ HU.assertBool "" $ ok $ groundDistance g1 g2
   where
      ok Nothing = False
      ok (Just (d1, a1, b1)) =
         abs (d - d1) < 0.01 *~ meter
         && abs (a - a1) < 0.01 *~ arcsecond
         && abs (b - b1) < 0.01 *~ arcsecond

-- | Sample points for UK tests. The oracle for these values is the script at
-- <http://www.movable-type.co.uk/scripts/latlong-convert-coords.html>, which uses
-- the same Helmert transform as this library. Hence the results should match to within 30 cm.
ukPoints :: [(String, Geodetic WGS84, Geodetic OSGB36)]
ukPoints = [
   ("Greenwich",        Geodetic (dms 51 28 40.86) (dms 0 0 (-5.83)) _0 WGS84,
                        Geodetic (dms 51 28 39.00) (dms 0 0 0) _0 OSGB36),
   ("Edinburgh Castle", Geodetic (dms 55 56 56.30) (dms (-3) (-12) (-2.73)) _0 WGS84,
                        Geodetic (dms 55 56 56.51) (dms (-3) (-11) (-57.61)) _0 OSGB36),
   ("Lands End",        Geodetic (dms 50 03 56.68) (dms (-5) (-42) (-51.20)) _0 WGS84,
                        Geodetic (dms 50 03 54.51) (dms (-5) (-42) (-47.87)) _0 OSGB36),
   ("Gt. Yarmouth Pier",Geodetic (dms 52 36 29.33) (dms 1 44 27.79) _0 WGS84,
                        Geodetic (dms 52 36 27.84) (dms 1 44 34.52) _0 OSGB36),
   ("Stanhope",         Geodetic (dms 54 44 49.08) (dms (-2) 0 (-19.89)) _0 WGS84,
                        Geodetic (dms 54 44 48.71) (dms (-2) 0 (-14.41)) _0 OSGB36) ]



-- Convert a named point into a test
pointTest :: (Ellipsoid e2) => (String, Geodetic WGS84, Geodetic e2) -> Test
pointTest (testName, wgs84, local) =  testCase testName $ HU.assertBool "" $ samePlace wgs84 (toWGS84 local)


-- The negation of the sum of a list of offsets is equal to the sum of the negated items.
prop_offset1 :: [GridOffset] -> Bool
prop_offset1 offsets = sameOffset (offsetNegate $ mconcat offsets) (mconcat $ map offsetNegate offsets)

-- A polar offset multiplied by a scalar is equal to an offset in the same direction with the length multiplied.
prop_offset2 :: Distance -> Bearing -> Scalar -> Bool
prop_offset2 (Distance d) (Bearing h) (Scalar s) = sameOffset go1 go2
   where
      go1 = offsetScale s $ polarOffset d h
      go2 = polarOffset (d * s) h

-- | A polar offset has the offset distance and bearing of its arguments.
prop_offset3 :: GridOffset -> Bool
prop_offset3 delta = sameOffset delta0
                                (polarOffset (offsetDistance delta0) (offsetBearing delta))
   where delta0 = delta {deltaAltitude = 0 *~ meter}

-- | Given a grid point and an offset, applying the offset to the point gives a new point which
-- is offset from the first point by the argument offset.
prop_grid1 :: GridPoint (GridTM LocalEllipsoid) -> GridOffset -> Bool
prop_grid1 p d = sameOffset d $ p `gridOffset` applyOffset d p

-- | Check that using toGrid/fromGrid for TransverseMercator projection are inverses
-- | for negative latitudes near the coordinates 0,0
prop_tmGridInverse :: Bool
prop_tmGridInverse = 
   let origin = Geodetic 
         { latitude = 0 *~ degree
         , longitude = 0 *~ degree
         , geoAlt = 0 *~ meter
         , ellipsoid = WGS84
         }
       g = mkGridTM origin mempty (1 *~ one)
       testPoint = origin { latitude = (-0.02) *~ degree }
   in fromGrid (toGrid g testPoint) `closeEnough` testPoint
   
-- | Converting a UK grid reference to a GridPoint and back is a null operation.
prop_ukGrid1 :: GridRef -> Bool
prop_ukGrid1 (GridRef str) =
   str ==
   (fromJust $ toUkGridReference ((length str P.- 2) `div` 2) $ fst $ fromJust $ fromUkGridReference str)

-- | UK Grid Reference points. The oracle for these points was the
-- UK Grid Reference Finder (gridreferencefinder.com), retrieved on 26 Jan 2013.
ukSampleGrid :: [(String, GridPoint UkNationalGrid, Geodetic WGS84, String)]
ukSampleGrid = map convert [
 -- Grid Reference, X,      Y,      Latitude,  Longitude,    Description
   ("SW3425625070", 134256, 025070, 50.066230, -5.7148278,   "Lands End"),
   ("TR3302139945", 633021, 139945, 51.111396,  1.3277159,   "Dover Harbour"),
   ("TQ3001980417", 530019, 180417, 51.507736, -0.12793230,  "Nelsons Column"),
   ("TA2542370644", 525423, 470644, 54.116376, -0.082668990, "Flamborough Lighthouse"),
   ("NK1354745166", 413547, 845166, 57.496512, -1.7756310,   "Peterhead harbour"),
   ("ND3804872787", 338048, 972787, 58.638518, -3.0688688,   "John O Groats"),
   ("SC3915875189", 239158, 475189, 54.147275, -4.4641148,   "Douglas Harbour"),
   ("ST1922474591", 319224, 174591, 51.464505, -3.1641741,   "Torchwood HQ"),
   ("SK3520736502", 435207, 336502, 52.924784, -1.4777486,   "Derby Cathedral")]
   where
      convert (grid, x, y, lat, long, desc) =
         (grid, GridPoint (x *~ meter) (y *~ meter) (0 *~ meter) UkNationalGrid,
          Geodetic (lat *~ degree) (long *~ degree) (0 *~ meter) WGS84, desc)

type GridPointTest = (String, GridPoint UkNationalGrid, Geodetic WGS84, String) -> Test

-- | Check that grid reference to grid point works for sample points.
ukGridTest2 :: GridPointTest
ukGridTest2 (gridRef, gp, _, testName) = testCase testName $ HU.assertBool ""
   $ (fst $ fromJust $ fromUkGridReference gridRef) == gp

-- | Check that grid point to grid reference works for sample points.
ukGridTest3 :: GridPointTest
ukGridTest3 (gridRef, gp, _, testName) = testCase testName $ HU.assertBool ""
   $ toUkGridReference 5 gp == Just gridRef

-- | Check that grid point to WGS84 works close enough for sample points.
ukGridTest4 :: GridPointTest
ukGridTest4 (_, gp, geo, testName) = testCase testName $ HU.assertBool ""
   $ closeEnough geo $ toWGS84 $ fromGrid gp

-- | Check that WGS84 to grid point works close enough for sample points.
ukGridTest5 :: GridPointTest
ukGridTest5 (_, gp, geo, testName) = testCase testName $ HU.assertBool ""
   $ offsetDistance (gridOffset gp $ toGrid UkNationalGrid $ toLocal OSGB36 geo) < 1 *~ meter


-- | Worked example for UK Geodetic to GridPoint, taken from "A Guide to Coordinate Systems in Great Britain" [1]
ukTest :: Geodetic OSGB36
ukTest = Geodetic (dms 52 39 27.2531) (dms 1 43 4.5177) (0 *~ meter) OSGB36

{-
   v = 6.3885023333E+06
   rho = 6.3727564399E+06
   eta2 = 2.4708136169E-03
   m = 4.0668829596E+05
   I = 3.0668829596E+05
   II = 1.5404079092E+06
   III = 1.5606875424E+05
   IIIa = -2.0671123011E+04
   IV = 3.8751205749E+06
   V = -1.7000078208E+05
   VI = -1.0134470432E+05
   E = 651409.903 m
   N = 313177.270 m
-}


-- | Standard stereographic grid for point tests in the Northern Hemisphere.
stereoGridN :: GridStereo LocalEllipsoid
stereoGridN = mkGridStereo tangent origin (0.9999079 *~ one)
   where
      ellipse = LocalEllipsoid "Bessel 1841" (6377397.155 *~ metre) (299.15281 *~ one) mempty
      tangent = Geodetic (dms 52 9 22.178) (dms 5 23 15.500) (0 *~ meter) ellipse
      origin = GridOffset (155000 *~ metre) (463000 *~ metre) (0 *~ meter)


-- | Standard steregraphic grid for point tests in the Southern Hemisphere.
--
-- This is the same as stereoGridN but with the tangent latitude and the false origin northings negated.
stereoGridS :: GridStereo LocalEllipsoid
stereoGridS = mkGridStereo tangent origin (0.9999079 *~ one)
   where
      ellipse = LocalEllipsoid "Bessel 1841" (6377397.155 *~ metre) (299.15281 *~ one) mempty
      tangent = Geodetic (negate $ dms 52 9 22.178) (dms 5 23 15.500) (0 *~ meter) ellipse
      origin = GridOffset ((-155000) *~ metre) (463000 *~ metre) (0 *~ meter)


-- | Data for the stereographic tests taken from
-- <http://ftp.stu.edu.tw/BSD/NetBSD/pkgsrc/distfiles/epsg-6.11/G7-2.pdf>
stereographicToGridN :: Bool
stereographicToGridN = sameGrid g1 g1'
   where
      p1 = Geodetic (dms 53 0 0) (dms 6 0 0) (0 *~ meter) $ gridEllipsoid stereoGridN
      g1 = GridPoint (196105.283 *~ meter) (557057.739 *~ meter) (0 *~ meter) stereoGridN
      g1' = toGrid stereoGridN p1

stereographicFromGridN :: Bool
stereographicFromGridN = samePlace p1 p1'
   where
      p1 = Geodetic (dms 53 0 0) (dms 6 0 0) (0 *~ meter) $ gridEllipsoid stereoGridN
      g1 = GridPoint (196105.283 *~ meter) (557057.739 *~ meter) (0 *~ meter) stereoGridN
      p1' = fromGrid g1


stereographicToGridS :: Bool
stereographicToGridS = sameGrid g1 g1'
   where
      p1 = Geodetic (negate $ dms 53 0 0) (dms 6 0 0) (0 *~ meter) $ gridEllipsoid stereoGridS
      g1 = GridPoint ((-196105.283) *~ meter) (557057.739 *~ meter) (0 *~ meter) stereoGridS
      g1' = toGrid stereoGridS p1


stereographicFromGridS :: Bool
stereographicFromGridS = samePlace p1 p1'
   where
      p1 = Geodetic (negate $ dms 53 0 0) (dms 6 0 0) (0 *~ meter) $ gridEllipsoid stereoGridS
      g1 = GridPoint ((-196105.283) *~ meter) (557057.739 *~ meter) (0 *~ meter) stereoGridS
      p1' = fromGrid g1


-- | Check the round trip for a stereographic projection.
prop_stereographic :: GridPoint (GridStereo LocalEllipsoid) -> Property
prop_stereographic p =
   let g = fromGrid p
       r = toGrid (gridBasis p) g
   in counterexample ("p = " ++ show p ++ "\ng = " ++ show g ++ "\nr = " ++ show r) $
     closeGrid p r



-- | A ray at distance zero returns its original arguments.
prop_rayPath1 :: Ray WGS84 -> Bool
prop_rayPath1 r@(Ray pt b e) =
      samePlace pt pt1 && sameAngle b b1 && sameAngle e e1
   where (pt1,b1,e1) = pathFunc (getRay r) _0


type ContinuityTest e = Geodetic e -> Bearing -> Azimuth -> Distance -> Distance -> Property

type ContinuityTest1 e = Geodetic e -> Bearing -> Distance2 -> Distance2 -> Property

-- | Many paths can be specified by a start point, bearing and azimuth,
-- and have the property that any (point,bearing,azimuth) triple on
-- the path will specify the same path with a distance offset.
prop_pathContinuity :: (Ellipsoid e) =>
   (Geodetic e -> Angle Double -> Angle Double -> Path e) -> ContinuityTest e
prop_pathContinuity pf pt0 (Bearing b0) (Azimuth a0) (Distance d1) (Distance d2) =
   counterexample (show ((pt2, Bearing b2, Azimuth a2), (pt3, Bearing b3, Azimuth a3))) $
      pathValidAt path0 d1 && pathValidAt path0 d2 && pathValidAt path0 (d1+d2) ==>
      closeEnough pt2 pt3 && sameAngle b2 b3 && sameAngle a2 a3
   where
      path0 = pf pt0 b0 a0
      (pt1, b1, a1) = pathFunc path0 d1
      path1 = pf pt1 b1 a1
      (pt2, b2, a2) = pathFunc path1 d2
      (pt3, b3, a3) = pathFunc path0 (d1 + d2)  -- Points 2 and 3 should be the same.


-- | For continuity testing of ground-based paths (azimuth & altitude always zero)
-- where lower accuracy is required.
prop_pathContinuity1 :: (Ellipsoid e) => (Geodetic e -> Angle Double -> Path e) -> ContinuityTest1 e
prop_pathContinuity1 pf pt0 (Bearing b0) (Distance2 d1) (Distance2 d2) =
   counterexample (show ((pt2, Bearing b2), (pt3, Bearing b3))) $
      pathValidAt path0 d1 && pathValidAt path0 d2 && pathValidAt path0 (d1+d2) ==>
      closeEnough pt2 pt3 && sameAngle b2 b3
   where
      path0 = pf pt0 b0
      (pt1, b1, _) = pathFunc path0 d1
      path1 = pf pt1 b1
      (pt2, b2, _) = pathFunc path1 d2
      (pt3, b3, _) = pathFunc path0 (d1 + d2)  -- Points 2 and 3 should be the same.


-- | A point on a ray will continue along the same ray, and hence give the same points.
prop_rayContinuity :: ContinuityTest WGS84
prop_rayContinuity = prop_pathContinuity rayPath


-- | A ray bisected to an altitude will give that altitude.
-- This is a test of bisection rather than rays.
prop_rayBisect :: Ray WGS84 -> Altitude -> Bool
prop_rayBisect r (Altitude height) =
   case bisect ray0 f (1 *~ centi meter) (0 *~ meter) (1000 *~ kilo meter) of
      Nothing -> False
      Just d -> let (g, _, _) = pathFunc ray0 d in abs (altitude g - height) < 1 *~ centi meter
   where
      f g = compare (altitude g) height
      ray0 = getRay r


-- | A point on a rhumb line will continue along the same rhumb.
prop_rhumbContinuity :: ContinuityTest1 WGS84
prop_rhumbContinuity = prop_pathContinuity1 rhumbPath


-- | Two rhumb paths intersect at the same place.
prop_rhumbIntersect :: RhumbPaths2 -> Property
prop_rhumbIntersect rp =
   case intersect _0 _0 (10.0 *~ centi meter) 100 path1 path2 of
      Just (d1, d2) ->
         let (pt1, _, _) = pathFunc path1 d1
             (pt2, _, _) = pathFunc path2 d2
         in counterexample (show (pt1, pt2)) $ label "Intersection" $ samePlace pt1 pt2
      Nothing -> label "No intersection" True
   where
      (path1, path2) = mk2RhumbPaths rp
