{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}
{-# HLINT ignore "Redundant bracket" #-}

module Main where

import Control.Monad
import Data.Char
import Data.Either
import Data.Maybe
import Test.Hspec
import Test.Hspec.QuickCheck
import Test.HUnit (assertFailure)
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
import Geodetics.UTM
import LatLongParser (parserTests)

main :: IO ()
main = hspec $ do
   describe "Geodetic" $ do
      prop "WGS84 and back" prop_WGS84_and_back
      prop "Zero ground distance" prop_zero_ground
      describe "UK Points" $ mapM_ pointTest ukPoints
      describe "World line" $ mapM_ worldLineTests worldLines
   parserTests
   describe "Grid" $ do
      prop "Grid Offset 1" prop_offset1
      prop "Grid Offset 2" prop_offset2
      prop "Grid Offset 3" prop_offset3
      prop "Grid 1" prop_grid1
   describe "TransverseMercator" $ prop "fromGrid . toGrid == id" prop_tmGridInverse
   describe "UK" $ do
      prop "UK Grid 1" prop_ukGrid1
      describe "UK Grid 2" $ mapM_ ukGridTest2 ukSampleGrid
      describe "UK Grid 3" $ mapM_ ukGridTest3 ukSampleGrid
      describe "UK Grid 4" $ mapM_ ukGridTest4 ukSampleGrid
      describe "UK Grid 5" $ mapM_ ukGridTest5 ukSampleGrid
   describe "UTM" $ do
      prop "UTM Grid 1" prop_utmGridTest1
      describe "UTM Grid 2" $ mapM_ utmGridTest2 utmSampleGrid
      describe "UTM Grid 3" $ mapM_ utmGridTest3 utmSampleGrid
      describe "UTM Grid 4" $ mapM_ utmGridTest4 utmSampleGrid
      describe "UTM Grid 5" $ mapM_ utmGridTest5 utmSampleGrid
   describe "MGRS" $ do
      prop "MGRS Grid 1" prop_mgrs_gridTest1
      prop "MGRS Grid 2" prop_mgrs_gridTest2
      describe "MGRS Grid 3" $ mapM_ mgrsGridTest3 utmSampleGrid
      describe "MGRS Grid 4" $ mapM_ mgrsGridTest4 utmSampleGrid
   describe "Stereographic" $ do
      it "toGrid north" stereographicToGridN
      it "fromGrid north" stereographicFromGridN
      it "toGrid south" stereographicToGridS
      it "fromGrid south" stereographicFromGridS
      prop "Stereographic round trip" prop_stereographic
   describe "Paths" $ do
      prop "Ray Path 1" prop_rayPath1
      prop "Ray Continuity" prop_rayContinuity
      prop "Ray Bisection" prop_rayBisect
      prop "Rhumb Continuity" prop_rhumbContinuity
      prop "Rhumb Intersection" prop_rhumbIntersect
   describe "GridOffset monoid" $ mapM_ (uncurry prop) $ unbatch $ monoid (mempty :: GridOffset)
   describe "Helmert monoid" $ mapM_ (uncurry prop) $ unbatch $ monoid (mempty :: Helmert)

instance EqProp GridOffset where
  (GridOffset a b c) =-= (GridOffset a' b' c') =
    eq True $ a ≈ a' && b ≈ b' && c ≈ c'
    where x ≈ y = abs (x - y) < 0.00001

instance EqProp Helmert where
  (Helmert cX' cY' cZ' s rX' rY' rZ') =-= (Helmert cX'' cY'' cZ'' s' rX'' rY'' rZ'') =
    eq True $ and [cX' ≈ cX'', cY' ≈ cY'', cZ' ≈ cZ'',
                   s ≈- s',
                   rX' ≈- rX'', rY' ≈- rY'', rZ' ≈- rZ'']

    where x ≈ y = abs (x - y) < 0.00001
          x ≈- y = abs (x - y) < 1 / (5 * 2) ^ _5


-- | The positions are within 30 cm.
samePlace :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Expectation
samePlace p1 p2 = expectTrue msg $ geometricalDistance p1 p2 < 0.3
   where
      msg = "location " <> show p2 <> " is > 30cm from expected " <> show p1

samePlace' :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
samePlace' p1 p2 = geometricalDistance p1 p2 < 0.3

-- | The positions are within 10 m.
closeEnough :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Expectation
closeEnough p1 p2 = expectTrue msg $ geometricalDistance p1 p2 < 10
   where
      msg = "location " <> show p2 <> " is > 10m from expected " <> show p1

closeEnough' :: (Ellipsoid e) => Geodetic e -> Geodetic e -> Bool
closeEnough' p1 p2 = geometricalDistance p1 p2 < 10

-- | The angles are within 0.01 arcsec
sameAngle :: Double -> Double -> Expectation
sameAngle v1 v2 = expectTrue msg $ abs (properAngle (v1 - v2)) < 0.01 * arcsecond
   where
      msg = "expected angle " <> show (v1 / degree) <> ", got " <> show (v2 / degree)

sameAngle' :: Double -> Double -> Bool
sameAngle' v1 v2 = abs (properAngle (v1 - v2)) < 0.01 * arcsecond

-- | The grid positions are within 1mm
sameGrid :: (Show r) => GridPoint r -> GridPoint r -> Expectation
sameGrid p1 p2 = expectTrue msg $ check eastings && check northings && check altitude
   where
      msg = "expected " <> show p1 <> ", got " <> show p2
      check f = f p1 - f p2 < 1e-3

-- | Grid offsets are within 1mm.
sameOffset :: GridOffset -> GridOffset -> Expectation
sameOffset go1 go2 = expectTrue msg $ check deltaNorth && check deltaEast && check deltaAltitude
   where
      msg = "expected " <> show go1 <> ", got " <> show go2
      check f = f go1 - f go2 < 1e-3


-- | The grid X and Y are both within 1 meter
closeGrid :: (Show r) => GridPoint r -> GridPoint r -> Expectation
closeGrid p1 p2 = expectTrue msg $ check eastings && check northings && check altitude
   where
      msg = "expected " <> show p1 <> ", got " <> show p2
      check f = f p1 - f p2 < 1


-- | Degrees, minutes and seconds into radians.
dms :: Int -> Int -> Double -> Double
dms d m s = fromIntegral d * degree + fromIntegral m * arcminute + s * arcsecond

-- | Round-trip from local to WGS84 and back is identity (approximately)
prop_WGS84_and_back :: Geodetic LocalEllipsoid -> Expectation
prop_WGS84_and_back p = samePlace p $ toLocal (ellipsoid p) $ toWGS84 p


-- | Test that for all points p, the ground distance from p to p is zero.
prop_zero_ground :: Geodetic WGS84 -> Bool
prop_zero_ground p =
   case groundDistance p p of
      Nothing -> False
      Just (d, _, _) -> abs d < 1e-3


-- | Sample pairs of points with bearings and distances.
-- The Oracle for these values is the @FORWARD@ program from
--  <http://www.ngs.noaa.gov/TOOLS/Inv_Fwd/Inv_Fwd.html>
worldLines :: [(String, Geodetic WGS84, Geodetic WGS84, {-Length-} Double, {-Angle-} Double, {-Angle-} Double)]
worldLines = [
   ("Ordinary", Geodetic (40 * degree) (30 * degree) 0 WGS84, Geodetic (30 * degree) (50 * degree) 0 WGS84,
      2128852.999, 115.19596706 * degree, 126.79044315 * degree),
   ("Over Pole", Geodetic (60 * degree) (0 * degree) 0 WGS84, Geodetic (60 * degree) (180 * degree) 0 WGS84,
      6695785.820, 0 * degree, 180 * degree),
   ("Equator to Pole", Geodetic (0 * degree) (0 * degree) 0 WGS84, Geodetic (90 * degree) (180 * degree) 0 WGS84,
      10001965.729, 0 * degree, 180 * degree)]


worldLineTests :: (String, Geodetic WGS84, Geodetic WGS84, Double, Double, Double) -> SpecWith (Arg Expectation)
worldLineTests (str, g1, g2, d, a, b) = it str $ ok $ groundDistance g1 g2
   where
      ok Nothing = False
      ok (Just (d1, a1, b1)) =
         abs (d - d1) < 0.01
         && abs (a - a1) < 0.01 * arcsecond
         && abs (b - b1) < 0.01 * arcsecond

-- | Sample points for UK tests. The oracle for these values is the script at
-- <http://www.movable-type.co.uk/scripts/latlong-convert-coords.html>, which uses
-- the same Helmert transform as this library. Hence the results should match to within 30 cm.
ukPoints :: [(String, Geodetic WGS84, Geodetic OSGB36)]
ukPoints = [
   ("Greenwich",        Geodetic (dms 51 28 40.86) (dms 0 0 (-5.83)) 0 WGS84,
                        Geodetic (dms 51 28 39.00) (dms 0 0 0) 0 OSGB36),
   ("Edinburgh Castle", Geodetic (dms 55 56 56.30) (dms (-3) (-12) (-2.73)) 0 WGS84,
                        Geodetic (dms 55 56 56.51) (dms (-3) (-11) (-57.61)) 0 OSGB36),
   ("Lands End",        Geodetic (dms 50 03 56.68) (dms (-5) (-42) (-51.20)) 0 WGS84,
                        Geodetic (dms 50 03 54.51) (dms (-5) (-42) (-47.87)) 0 OSGB36),
   ("Gt. Yarmouth Pier",Geodetic (dms 52 36 29.33) (dms 1 44 27.79) 0 WGS84,
                        Geodetic (dms 52 36 27.84) (dms 1 44 34.52) 0 OSGB36),
   ("Stanhope",         Geodetic (dms 54 44 49.08) (dms (-2) 0 (-19.89)) 0 WGS84,
                        Geodetic (dms 54 44 48.71) (dms (-2) 0 (-14.41)) 0 OSGB36) ]



-- Convert a named point into a test
pointTest :: (Ellipsoid e2) => (String, Geodetic WGS84, Geodetic e2) -> SpecWith (Arg Expectation)
pointTest (testName, wgs84, local) =  it testName $ wgs84 `samePlace` toWGS84 local


-- The negation of the sum of a list of offsets is equal to the sum of the negated items.
prop_offset1 :: [GridOffset] -> Expectation
prop_offset1 offsets = sameOffset (offsetNegate $ mconcat offsets) (mconcat $ map offsetNegate offsets)

-- A polar offset multiplied by a scalar is equal to an offset in the same direction with the length multiplied.
prop_offset2 :: Distance -> Bearing -> Scalar -> Expectation
prop_offset2 (Distance d) (Bearing h) (Scalar s) = sameOffset go1 go2
   where
      go1 = offsetScale s $ polarOffset d h
      go2 = polarOffset (d * s) h

-- | A polar offset has the offset distance and bearing of its arguments.
prop_offset3 :: GridOffset -> Expectation
prop_offset3 delta = sameOffset delta0
                                (polarOffset (offsetDistance delta0) (offsetBearing delta))
   where delta0 = delta {deltaAltitude = 0}

-- | Given a grid point and an offset, applying the offset to the point gives a new point which
-- is offset from the first point by the argument offset.
prop_grid1 :: GridPoint (GridTM LocalEllipsoid) -> GridOffset -> Expectation
prop_grid1 p d = sameOffset d $ p `gridOffset` applyOffset d p

-- | Check that using toGrid/fromGrid for TransverseMercator projection are inverses
-- | for negative latitudes near the coordinates 0,0
prop_tmGridInverse :: Expectation
prop_tmGridInverse =
   let origin = Geodetic
         { latitude = 0 * degree
         , longitude = 0 * degree
         , geoAlt = 0
         , ellipsoid = WGS84
         }
       g = mkGridTM origin mempty 1
       testPoint = origin { latitude = (-1) * arcminute }
       tp1 = toGrid g testPoint
       tp2 = fromGrid tp1
   in tp2 `closeEnough` testPoint

-- | Converting a UK grid reference to a GridPoint and back is a null operation.
prop_ukGrid1 :: UkGridRef -> Expectation
prop_ukGrid1 (UkGridRef str) =
   str `shouldBe`
   fromJust (toUkGridReference ((length str - 2) `div` 2) $ fst $ fromJust $ fromUkGridReference str)

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
   ("SK3520736502", 435207, 336502, 52.924784, -1.4777486,   "Derby Cathedral"),
   ("TG5141013177", 651410, 313177, 52.657979 , 1.7160519,   "Caister Water Tower"),
   ("TG2623802646", 626238, 302646, 52.574548 , 1.3373749,   "Framingham")]
   -- Caister and Framingham are taken from Ordnance Survey worked examples.
   where
      convert (grid, x, y, lat, long, desc) =
         (grid, GridPoint x y 0 UkNationalGrid,
          Geodetic (lat * degree) (long * degree) 0 WGS84, desc)

type UkGridPointTest = (String, GridPoint UkNationalGrid, Geodetic WGS84, String) -> SpecWith (Arg Expectation)

-- | Check that grid reference to grid point works for sample points.
ukGridTest2 :: UkGridPointTest
ukGridTest2 (gridRef, gp, _, testName) =
   it testName $ fst (fromJust $ fromUkGridReference gridRef) `shouldBe` gp

-- | Check that grid point to grid reference works for sample points.
ukGridTest3 :: UkGridPointTest
ukGridTest3 (gridRef, gp, _, testName) =
   it testName $ toUkGridReference 5 gp `shouldBe` Just gridRef

-- | Check that grid point to WGS84 works close enough for sample points.
ukGridTest4 :: UkGridPointTest
ukGridTest4 (_, gp, geo, testName) =
   it testName $ geo `closeEnough` toWGS84 (fromGrid gp)

-- | Check that WGS84 to grid point works close enough for sample points.
ukGridTest5 :: UkGridPointTest
ukGridTest5 (_, gp, geo, testName) =
   it testName $ gp `closeGrid` toGrid UkNationalGrid (toLocal OSGB36 geo)


-- | Worked example for UK Geodetic to GridPoint, taken from "A Guide to Coordinate Systems in Great Britain" [1]
ukTest :: Geodetic OSGB36
ukTest = Geodetic (dms 52 39 27.2531) (dms 1 43 4.5177) 0 OSGB36

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


-- | Check that a UTM grid reference round-trips to a gridpoint and back.
prop_utmGridTest1 :: UtmGridRef -> Expectation
prop_utmGridTest1 (UtmGridRef str) =
   str `shouldBe`
   toUtmGridReference Nothing True 0 (fromRight (error str) $ fromUtmGridReference str)

-- | Sample points for UTM, in both UTM and MGRS formats. The oracle for these points was the Montana
-- State University converter. http://rcn.montana.edu/resources/Converter.aspx.
-- Retrieved on 1st Feb 2025.
utmSampleGrid :: [(String, String, GridPoint UtmZone, Geodetic WGS84, String)]
utmSampleGrid = map convert [
 -- UTM Reference,        MGRS Reference,       X,      Y,        Latitude,    Longitude,   Description  
   ("30N 699304 5710208", "30U XC 99304 10208", 699304, 5710208,   51.5078064,  -0.1279388, "Nelson's Column"),
   ("35S 379101 8017747", "35K LA 79101 17747", 379101, 8017747,  -17.9249501,  25.8585071, "Victoria Falls"),
   ("27N 454366 7111715", "27W VM 54366 11715", 454366, 7111715,   64.1289075, -21.9373288, "Reykjavik Airport"),
   ("32N 297697 6700532", "32V KN 97697 00532", 297697, 6700532,   60.3904298,   5.3284215, "Bergen"),
   ("23S 683473 7460697", "23K PQ 83473 60697", 683473, 7460697,  -22.9518122, -43.2105383, "Christ the Redeemer"),
   ("21S 439699 4272868", "21F VC 39699 72868", 439699, 4272868,  -51.6919073, -57.8724006, "Port Stanley"),
   ("59S 461474 4919822", "59G MK 61474 19822", 461474, 4919822,  -45.8740880, 170.5035807, "Dunedin"),
   ("31N 166022 0"      , "31N AA 66022 00000", 166022,       0,    0.0      ,   0.0      , "The Origin"),
   ("31S 166022 9999999", "31M AV 66022 99999", 166022, 9999999,  -0.00000903,   0.0      , "1 meter south")]
   where
      convert (utm, mgrs, x, y, lat, long, desc) =
         (utm, mgrs, GridPoint x y 0 zone, geo, desc)
         where
            geo = Geodetic (lat * degree) (long * degree) 0 WGS84
            znum = fromJust $ utmZoneNumber geo
            hemi = if lat < 0 then UtmSouth else UtmNorth
            zone = fromJust $ mkUtmZone hemi znum

type UtmGridPointTest = (String, String, GridPoint UtmZone, Geodetic WGS84, String) -> SpecWith (Arg Expectation)


-- | Check that UTM reference to grid point works for sample points.
utmGridTest2 :: UtmGridPointTest
utmGridTest2 (gridRef, _, gp, _, testName) =
   it testName $ fromUtmGridReference gridRef `shouldBe` Right gp

-- | Check that Grid point to UTM reference works for sample points.
utmGridTest3 :: UtmGridPointTest
utmGridTest3 (gridRef, _, gp, _, testName) =
   it testName $ toUtmGridReference Nothing False 0 gp `shouldBe` gridRef

-- | Check that grid point to WGS84 works close enough for sample points.
utmGridTest4 :: UtmGridPointTest
utmGridTest4 (_, _, gp, geo, testName) =
   it testName $ geo `closeEnough` fromGrid gp

-- | Check that WGS84 to grid point works close enough for sample points.
utmGridTest5 :: UtmGridPointTest
utmGridTest5 (_, _, gp, geo, testName) =
   it testName $ gp `closeGrid` toGrid (fromJust $ utmZone geo) geo


-- | Check that a UTM grid point round-trips to MGRS and back, with spaces.
prop_mgrs_gridTest1 :: UtmGridRef -> Expectation
prop_mgrs_gridTest1 (UtmGridRef str) =
   case fromUtmGridReference str of
      Left msg -> assertFailure $ "gridTest1 bogus UTM ref: " <> str <> ". Messages = " <> show msg
      Right gp ->
         fromMgrsGridReference <$> toMgrsGridReference True 5 gp `shouldBe` Just (Right (gp, GridOffset 0.5 0.5 0))

-- | Check that a UTM grid point round-trips to MGRS and back, without spaces.
prop_mgrs_gridTest2 :: UtmGridRef -> Expectation
prop_mgrs_gridTest2 (UtmGridRef str) =
   case fromUtmGridReference str of
      Left msg -> assertFailure $ "gridTest2 bogus UTM ref: " <> str <> ". Messages = " <> show msg
      Right gp ->
         fromMgrsGridReference <$> toMgrsGridReference False 5 gp `shouldBe` Just (Right (gp, GridOffset 0.5 0.5 0))

-- | Check that MGRS reference to grid point works for sample points.
mgrsGridTest3 :: UtmGridPointTest
mgrsGridTest3 (_, mgrs, gp, _, testName) =
   it testName $ fromMgrsGridReference mgrs `shouldBe` Right (gp, GridOffset 0.5 0.5 0)

-- | Check that grid point to MGRS reference works for sample points.
mgrsGridTest4 :: UtmGridPointTest
mgrsGridTest4 (_, mgrs, gp, _, testName) = do
   it (testName <> " with spaces") $ toMgrsGridReference True 5 gp `shouldBe` Just mgrs
   it (testName <> " without spaces") $ toMgrsGridReference False 5 gp `shouldBe` Just (filter (not . isSpace) mgrs)


-- | Standard stereographic grid for point tests in the Northern Hemisphere.
stereoGridN :: GridStereo LocalEllipsoid
stereoGridN = mkGridStereo tangent origin 0.9999079
   where
      ellipse = LocalEllipsoid "Bessel 1841" 6377397.155 299.15281 mempty
      tangent = Geodetic (dms 52 9 22.178) (dms 5 23 15.500) 0 ellipse
      origin = GridOffset 155000 463000 0


-- | Standard steregraphic grid for point tests in the Southern Hemisphere.
--
-- This is the same as stereoGridN but with the tangent latitude and the false origin northings negated.
stereoGridS :: GridStereo LocalEllipsoid
stereoGridS = mkGridStereo tangent origin 0.9999079
   where
      ellipse = LocalEllipsoid "Bessel 1841" 6377397.155 299.15281 mempty
      tangent = Geodetic (negate $ dms 52 9 22.178) (dms 5 23 15.500) 0 ellipse
      origin = GridOffset (-155000) 463000 0


-- | Data for the stereographic tests taken from
-- <http://ftp.stu.edu.tw/BSD/NetBSD/pkgsrc/distfiles/epsg-6.11/G7-2.pdf>
stereographicToGridN :: Expectation
stereographicToGridN = sameGrid g1 g1'
   where
      p1 = Geodetic (dms 53 0 0) (dms 6 0 0) 0 $ gridEllipsoid stereoGridN
      g1 = GridPoint 196105.283 557057.739 0 stereoGridN
      g1' = toGrid stereoGridN p1

stereographicFromGridN :: Expectation
stereographicFromGridN = samePlace p1 p1'
   where
      p1 = Geodetic (dms 53 0 0) (dms 6 0 0) 0 $ gridEllipsoid stereoGridN
      g1 = GridPoint 196105.283 557057.739 0 stereoGridN
      p1' = fromGrid g1


stereographicToGridS :: Expectation
stereographicToGridS = sameGrid g1 g1'
   where
      p1 = Geodetic (negate $ dms 53 0 0) (dms 6 0 0) 0 $ gridEllipsoid stereoGridS
      g1 = GridPoint (-196105.283) 557057.739 0 stereoGridS
      g1' = toGrid stereoGridS p1


stereographicFromGridS :: Expectation
stereographicFromGridS = samePlace p1 p1'
   where
      p1 = Geodetic (negate $ dms 53 0 0) (dms 6 0 0) 0 $ gridEllipsoid stereoGridS
      g1 = GridPoint (-196105.283) 557057.739 0 stereoGridS
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
      samePlace' pt pt1 && sameAngle' b b1 && sameAngle' e e1
   where (pt1,b1,e1) = pathFunc (getRay r) 0


type ContinuityTest e = Geodetic e -> Bearing -> Azimuth -> Distance -> Distance -> Property

type ContinuityTest1 e = Geodetic e -> Bearing -> Distance2 -> Distance2 -> Property

-- | Many paths can be specified by a start point, bearing and azimuth,
-- and have the property that any (point,bearing,azimuth) triple on
-- the path will specify the same path with a distance offset.
prop_pathContinuity :: (Ellipsoid e) =>
   (Geodetic e -> Double -> Double -> Path e) -> ContinuityTest e
prop_pathContinuity pf pt0 (Bearing b0) (Azimuth a0) (Distance d1) (Distance d2) =
   counterexample (show ((pt2, Bearing b2, Azimuth a2), (pt3, Bearing b3, Azimuth a3))) $
      pathValidAt path0 d1 && pathValidAt path0 d2 && pathValidAt path0 (d1+d2) ==>
      closeEnough' pt2 pt3 && sameAngle' b2 b3 && sameAngle' a2 a3
   where
      path0 = pf pt0 b0 a0
      (pt1, b1, a1) = pathFunc path0 d1
      path1 = pf pt1 b1 a1
      (pt2, b2, a2) = pathFunc path1 d2
      (pt3, b3, a3) = pathFunc path0 (d1 + d2)  -- Points 2 and 3 should be the same.


-- | For continuity testing of ground-based paths (azimuth & altitude always zero)
-- where lower accuracy is required.
prop_pathContinuity1 :: (Ellipsoid e) => (Geodetic e -> Double -> Path e) -> ContinuityTest1 e
prop_pathContinuity1 pf pt0 (Bearing b0) (Distance2 d1) (Distance2 d2) =
   counterexample (show ((pt2, Bearing b2), (pt3, Bearing b3))) $
      pathValidAt path0 d1 && pathValidAt path0 d2 && pathValidAt path0 (d1+d2) ==>
      closeEnough' pt2 pt3 && sameAngle' b2 b3
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
   case bisect ray0 f 1e-2 0 (1000 * kilometer) of
      Nothing -> False
      Just d -> let (g, _, _) = pathFunc ray0 d in abs (altitude g - height) < 1e-2
   where
      f g = compare (altitude g) height
      ray0 = getRay r


-- | A point on a rhumb line will continue along the same rhumb.
prop_rhumbContinuity :: ContinuityTest1 WGS84
prop_rhumbContinuity = prop_pathContinuity1 rhumbPath


-- | Two rhumb paths intersect at the same place.
prop_rhumbIntersect :: RhumbPaths2 -> Property
prop_rhumbIntersect rp =
   case intersect 0 0 0.1 100 path1 path2 of
      Just (d1, d2) ->
         let (pt1, _, _) = pathFunc path1 d1
             (pt2, _, _) = pathFunc path2 d2
         in counterexample (show (pt1, pt2)) $ label "Intersection" $ samePlace pt1 pt2
      Nothing -> label "No intersection" True
   where
      (path1, path2) = mk2RhumbPaths rp


-- | Copied from `Test.Hspec.Expectations` source. It ought to be exported from there.
expectTrue :: HasCallStack => String -> Bool -> Expectation
expectTrue msg b = unless b (expectationFailure msg)