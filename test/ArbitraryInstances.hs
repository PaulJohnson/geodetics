{-# LANGUAGE FlexibleInstances #-}

-- | Orphan "Arbitrary" and related instances for testing purposes. 

module ArbitraryInstances where

import Control.Applicative
import Control.Monad
import Test.QuickCheck
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P

import Geodetics.Coordinates
import Geodetics.Ellipsoids


-- | Wrapper for arbitrary angles.
newtype Bearing = Bearing (Dimensionless Double) deriving (Show)

instance Arbitrary Bearing where
   arbitrary = Bearing <$> (*~ degree) <$> choose (0,360)
   shrink = shrinkNothing
   
-- | Wrapper for arbitrary distances up to 1000km
newtype Distance = Distance (Length Double) deriving (Show)

instance Arbitrary Distance where
   arbitrary = Distance <$> (*~ kilo meter) <$> choose (0,1000)
   shrink (Distance d) = Distance <$> (*~ kilo meter) <$> shrink (d /~ kilo meter)

-- | Wrapper for arbitrary dimensionless numbers (-10 .. 10)
newtype Scalar = Scalar (Dimensionless Double) deriving (Show)

instance Arbitrary Scalar where
   arbitrary = Scalar <$> (*~ one) <$> choose (-10,10)
   shrink (Scalar s) = Scalar <$> (*~ one) <$> shrink (s /~ one)

-- | Wrapper for arbitrary grid references.
newtype GridRef = GridRef String deriving Show

instance Arbitrary GridRef where
   arbitrary = do
      n <- choose (0,4)
      c1 <- elements "HJNOST" -- General vicinity of UK
      c2 <- elements $ ['A'..'H'] ++ ['J'..'Z']
      dx <- vectorOf n $ choose ('0','9')
      dy <- vectorOf n $ choose ('0','9')
      return $ GridRef $ c1 : c2 : (dx ++ dy)
   shrink = shrinkNothing

-- | Generate in range +/- <arg> m.
genOffset :: Double -> Gen (Length Double)
genOffset d = (*~ meter) <$> choose (-d, d)

genAlt :: Gen (Length Double)
genAlt = (*~ meter) <$> choose (0,10000)


genLatitude :: Gen (Dimensionless Double)
genLatitude = (*~ degree) <$> choose (-90,90)

genLongitude :: Gen (Dimensionless Double)
genLongitude = (*~ degree) <$> choose (-180,180)

genSeconds :: Gen (Dimensionless Double)
genSeconds = (*~ arcsecond) <$> choose (-10,10)
    

-- | Shrinking with the original value preserved. Used for shrinking records.  See 
-- http://stackoverflow.com/questions/14006005/idiomatic-way-to-shrink-a-record-in-quickcheck for details.
shrink' :: (Arbitrary a) => a -> [a]
shrink' x = x : shrink x

-- | Shrink a quantity in the given units.
shrinkQuantity :: (Arbitrary a, Fractional a) => Unit d a -> Quantity d a -> [Quantity d a]
shrinkQuantity u q = map (*~ u) $ shrink' $ q /~ u

shrinkLength :: (Arbitrary a, Fractional a) => Length a -> [Length a]
shrinkLength = shrinkQuantity meter

shrinkUnit :: (Arbitrary a, Fractional a) => Dimensionless a -> [Dimensionless a]
shrinkUnit = shrinkQuantity one



instance Arbitrary Helmert where
   arbitrary = 
      Helmert <$> genOffset 300 <*> genOffset 300 <*> genOffset 300 <*> 
         ((*~ one) <$> choose (-5,10)) <*>
         genSeconds <*> genSeconds <*> genSeconds
   shrink h = 
      tail $ Helmert <$> shrinkLength (cX h) <*> shrinkLength (cY h) <*> shrinkLength (cZ h) <*>
         shrinkUnit (helmertScale h) <*>
         shrinkUnit (rX h) <*> shrinkUnit (rY h) <*> shrinkUnit (rZ h)      


instance Arbitrary WGS84 where
   arbitrary = return WGS84
   shrink = shrinkNothing
   

instance Arbitrary LocalEllipsoid where
   arbitrary =
      LocalE <$> (("Local_" ++) <$> replicateM 3 (choose ('A','Z'))) <*>  -- name
         ((*~ meter) <$> choose (6378100, 6378400)) <*>                  -- majorRadius
         ((*~ one) <$> choose (297,300)) <*>                             -- flatR
         arbitrary                                                       -- helmert
   shrink e = tail $ LocalE (nameLocal e) (majorRadius e) (flatR e) <$> shrink' (helmert e)

        
instance (Ellipsoid e, Arbitrary e) => Arbitrary (Geodetic e) where
   arbitrary = 
      Geodetic <$> genLatitude <*> genLongitude <*> genOffset 1 <*> arbitrary
   shrink g = 
      tail $ Geodetic <$> shrinkUnit (latitude g) <*> shrinkUnit (longitude g) <*> 
         shrinkLength (altitude g) <*> shrink' (ellipsoid g)

instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridPoint (GridTM e)) where
   arbitrary = GridPoint <$> genOffset 100000 <*> genOffset 100000 <*> genOffset 1 <*> arbitrary
   shrink p = tail $ GridPoint <$> 
      shrinkLength (eastings p) <*> 
      shrinkLength (northings p) <*> 
      shrinkLength (altitude p) <*> 
      shrink' (gridBase p)


instance (Ellipsoid e, Arbitrary e) => Arbitrary (GridTM e) where
   arbitrary = mkGridTM <$> arbitrary <*> arbitrary <*> ((*~ one) <$> choose (0.95,1.0))
   shrink tm = tail $ mkGridTM <$> shrink' (trueOrigin tm) <*> shrink' (falseOrigin tm) <*> [gridScale tm]
   
instance Arbitrary GridOffset where
   arbitrary = GridOffset <$> genOffset 100000 <*> genOffset 100000 <*> genAlt
   shrink d = tail $ GridOffset <$> 
      shrinkLength (deltaEast d) <*> shrinkLength (deltaNorth d) <*> shrinkLength (deltaAltitude d)

         