{-# LANGUAGE FlexibleInstances, RankNTypes, KindSignatures, DataKinds, NumericUnderscores #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}
{-# HLINT ignore "Redundant <$>" #-}
{-# HLINT ignore "Use underscore" #-}
{-# HLINT ignore "Redundant bracket" #-}

-- | Orphan "Arbitrary" and related instances for testing purposes. 

module ArbitraryInstances where

import Control.Monad
import Geodetics.Altitude
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.Ellipsoids
import Geodetics.Path
import Geodetics.Stereographic as SG
import Geodetics.TransverseMercator as TM
import Test.QuickCheck
import Text.Printf


-- | Shrink an angle so that shrunk values are round numbers of degrees.
shrinkAngle :: Double -> [Double]
shrinkAngle v = (degree *) <$> shrink (v/degree)

-- | Shrink a distance so that shrunk values are round numbers of kilometers.
shrinkDistance :: Double -> [Double]
shrinkDistance v = (kilometer *) <$> shrink (v/kilometer)

-- | Wrapper for arbitrary angles.
newtype Bearing = Bearing Double

instance Show Bearing where
   show (Bearing b) = "Bearing " ++ showAngle b

instance Arbitrary Bearing where
   arbitrary = Bearing . (* degree) <$> choose (-180,180)
   shrink (Bearing b) = Bearing <$> shrinkAngle b


newtype Azimuth = Azimuth Double

instance Show Azimuth where
   show (Azimuth a) = "Azimuth " ++ showAngle a

instance Arbitrary Azimuth where
   arbitrary = Azimuth . (* degree) <$> choose (0,90)
   shrink (Azimuth a) = Azimuth <$> shrinkAngle a


-- | Wrapper for arbitrary distances up to 10,000 km
newtype Distance = Distance Double deriving (Show)

instance Arbitrary Distance where
   arbitrary = Distance . (* kilometer) <$> choose (0,10_000)
   shrink (Distance d) = Distance <$> shrinkDistance d


-- | Wrapper for arbitrary distances up to 1,000 km
newtype Distance2 = Distance2 Double deriving (Show)

instance Arbitrary Distance2 where
   arbitrary = Distance2 . (* kilometer) <$> choose (0,1_000)
   shrink (Distance2 d) = Distance2 <$> shrinkDistance d

-- | Wrapper for arbitrary altitudes up to 10 km
newtype Altitude = Altitude Double deriving (Show)

instance Arbitrary Altitude where
   arbitrary = Altitude . (* kilometer) <$> choose (0,10)
   shrink (Altitude h) = Altitude <$> shrinkDistance h


-- | Wrapper for arbitrary dimensionless numbers (-10 .. 10)
newtype Scalar = Scalar Double deriving (Show)

instance Arbitrary Scalar where
   arbitrary = Scalar <$> choose (-10,10)
   shrink (Scalar s) = Scalar <$> shrink s


-- | Wrapper for arbitrary UK grid references.
newtype UkGridRef = UkGridRef String deriving Show

instance Arbitrary UkGridRef where
   arbitrary = do
      n <- choose (0,4)
      c1 <- elements "HJNOST" -- General vicinity of UK
      c2 <- elements $ ['A'..'H'] ++ ['J'..'Z']
      dx <- vectorOf n $ choose ('0','9')
      dy <- vectorOf n $ choose ('0','9')
      return $ UkGridRef $ c1 : c2 : (dx ++ dy)
   shrink = shrinkNothing


-- | Wrapper for arbitrary UTM grid references.
newtype UtmGridRef = UtmGridRef String deriving Show

instance Arbitrary UtmGridRef where
   arbitrary = do
      zone <- choose (1,60 :: Int)
      hemi <- elements "NS"
      e <- choose (100_000, 899_999 :: Integer)
      n <- case hemi of
         'N' -> choose (0, 9_328_195 :: Integer) -- Equator to 84 degrees north.
         _   -> choose (1_118_248, 10_000_000)  -- 80 degrees south to Equator.
      return $ UtmGridRef $ printf "%02d" zone <> [hemi] <> " " <> show e <> "E " <> show n <> "N"


-- | Generate in range +/- <arg> m.
genOffset :: Double -> Gen Double
genOffset d = choose (-d, d)

genAlt :: Gen Double
genAlt = choose (0,10000)


genLatitude :: Gen Double
genLatitude = (* degree) <$> choose (-90,90)

genLongitude :: Gen Double
genLongitude = (* degree) <$> choose (-180,180)

genSeconds :: Gen Double
genSeconds = (* arcsecond) <$> choose (-10,10)


-- | Shrinking with the original value preserved. Used for shrinking records.  See 
-- http://stackoverflow.com/questions/14006005/idiomatic-way-to-shrink-a-record-in-quickcheck for details.
shrink' :: (Arbitrary a) => a -> [a]
shrink' x = x : shrink x

shrinkAngle' :: Double -> [Double]
shrinkAngle' a = a : shrinkAngle a


instance Arbitrary Helmert where
   arbitrary =
      Helmert <$> genOffset 300 <*> genOffset 300 <*> genOffset 300 <*>
         (choose (-5,10)) <*> genSeconds <*> genSeconds <*> genSeconds
   shrink h =
      drop 1 $ Helmert <$> shrink' (cX h) <*> shrink' (cY h) <*> shrink' (cZ h) <*>
         shrink' (helmertScale h) <*>
         shrink' (rX h) <*> shrink' (rY h) <*> shrink' (rZ h)


instance Arbitrary WGS84 where
   arbitrary = return WGS84
   shrink = shrinkNothing


instance Arbitrary LocalEllipsoid where
   arbitrary =
      (LocalEllipsoid . ("Local_" ++) <$> replicateM 3 (choose ('A','Z'))) <*>  -- name
         (choose (6378100, 6378400)) <*>                    -- majorRadius
         (choose (297,300)) <*>                             -- flatR
         arbitrary                                          -- helmert
   shrink e = drop 1 $ LocalEllipsoid (nameLocal e) (majorRadius e) (flatR e) <$> shrink' (helmert e)


instance (Ellipsoid e, Arbitrary e) => Arbitrary (Geodetic e) where
   arbitrary =
      Geodetic <$> genLatitude <*> genLongitude <*> genOffset 1 <*> arbitrary
   shrink g =
      drop 1 $ Geodetic <$> shrinkAngle' (latitude g) <*> shrinkAngle' (longitude g) <*>
         shrink' (altitude g) <*> shrink' (ellipsoid g)

instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridPoint (GridTM e)) where
   arbitrary = GridPoint <$> genOffset 100000 <*> genOffset 100000 <*> genOffset 1 <*> arbitrary
   shrink p = drop 1 $ GridPoint <$>
      shrink' (eastings p) <*>
      shrink' (northings p) <*>
      shrink' (altitude p) <*>
      shrink' (gridBasis p)


instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridPoint (GridStereo e)) where
   arbitrary = GridPoint <$> genOffset 100000 <*> genOffset 100000 <*> genOffset 1 <*> arbitrary
   shrink p = drop 1 $ GridPoint <$>
      shrink' (eastings p) <*>
      shrink' (northings p) <*>
      shrink' (altitude p) <*>
      shrink' (gridBasis p)


instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridTM e) where
   arbitrary = mkGridTM <$> arbitrary <*> arbitrary <*> choose (0.95,1.0)
   shrink tm = drop 1 $ mkGridTM <$> shrink' (trueOrigin tm) <*> shrink' (falseOrigin tm) <*> [TM.gridScale tm]


instance Arbitrary GridOffset where
   arbitrary = GridOffset <$> genOffset 100000 <*> genOffset 100000 <*> genAlt
   shrink d = drop 1 $ GridOffset <$>
      shrink' (deltaEast d) <*> shrink' (deltaNorth d) <*> shrink' (deltaAltitude d)


instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridStereo e) where
   arbitrary = mkGridStereo <$> arbitrary <*> arbitrary <*> choose (0.95,1.0)
   shrink sg = drop 1 $ mkGridStereo <$> shrink' (gridTangent sg) <*> shrink' (gridOrigin sg) <*> [SG.gridScale sg]


-- | Wrapper for arbitrary rays, along with creation parameters for printing and shrinking.
data Ray e = Ray (Geodetic e) Double Double

instance (Ellipsoid e) => Show (Ray e) where
   show (Ray p0 b e ) = "(Ray " ++ show p0 ++ ", " ++ showAngle b ++ ", " ++ showAngle e ++ ")"

getRay :: (Ellipsoid e) => Ray e -> Path e
getRay (Ray p0 b e) = rayPath p0 b e

instance (Ellipsoid e, Arbitrary e) => Arbitrary (Ray e) where
   arbitrary = do
      p0 <- arbitrary
      b <- (* degree) <$> choose (-180,180)
      e <- (* degree) <$> choose (0,90)
      return $ Ray p0 b e
   shrink (Ray p0 b e) = drop 1 $ do
      p0' <- shrink' p0
      b' <- shrinkAngle' b
      e' <- shrinkAngle' e
      return $ Ray p0' b' e'


-- | Two rhumb paths starting not more than 1000 km apart.
data RhumbPaths2 = RP2 {
      rp2Point0 :: Geodetic WGS84,
      rp2Bearing0 :: Bearing,
      rp2Distance :: Distance2,
      rp2Bearing1 :: Bearing,
      rp2Bearing2 :: Bearing
   }

instance Show RhumbPaths2 where
   show rp2 = show (pt1, Bearing b1) ++ show (pt2, Bearing b2)
      where
         (p1, p2) = mk2RhumbPaths rp2
         (pt1, b1, _) = pathFunc p1 0
         (pt2, b2, _) = pathFunc p2 0

instance Arbitrary RhumbPaths2 where
   arbitrary = RP2
      <$> arbitrary `suchThat` ((< (70 * degree)) . abs . latitude)
      <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
   shrink rp =
      drop 1 $ RP2 <$>
         shrink' (rp2Point0 rp) <*>
         shrink' (rp2Bearing0 rp) <*>
         shrink' (rp2Distance rp) <*>
         shrink' (rp2Bearing1 rp) <*>
         shrink' (rp2Bearing2 rp)

mk2RhumbPaths :: RhumbPaths2 -> (Path WGS84, Path WGS84)
mk2RhumbPaths (RP2 pt0 (Bearing b0) (Distance2 d) (Bearing b1) (Bearing b2)) =
   (path1, path2)
   where
      path0 = rhumbPath pt0 b0
      path1 = rhumbPath pt0 b1
      (pt2, _, _) = pathFunc path0 d
      path2 = rhumbPath pt2 b2

