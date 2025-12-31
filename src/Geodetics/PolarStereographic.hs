{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleContexts #-}

module Geodetics.PolarStereographic (
  Pole (..),
  PolarStereographic (trueOrigin, falseOrigin, polarEllipsoid, gridScale),
  mkGridPolarStereographic,
  UpsGrid,
  upsNorth,
  upsSouth,
  fromUpsGridReference,
  parseUpsGridReference,
  toUpsGridReference
) where

import Control.Monad
import Data.Char
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Text.Parsec
import Text.Parsec.Error
import Text.Read (readMaybe)
import Text.Printf

-- | Polar stereographic grids are defined for true origins at the north and south poles.
data Pole = NorthPole | SouthPole deriving (Show, Ord, Eq, Enum, Bounded)


{- | Polar Stereographic Grids

Formulae are taken from
/The Universal Grids: Univerersal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS)/
DMA Technical Manual 8358.2, Defense Mapping Agency, Fairfax, VA. https://apps.dtic.mil/sti/tr/pdf/ADA266497.pdf

When working with polar grids all directions are relative to the grid rather than the actual pole.
So in the Arctic \"North\" on the Universal Polar Stereographic grid means towards the Bering Sea
rather than towards the North Pole.
-}
data PolarStereographic e = PolarStereographic {
  trueOrigin :: Pole,
  falseOrigin :: GridOffset,
    -- ^ The negation of the grid position of the true origin. Used to avoid negative coordinates over the area
    -- of interest. The altitude gives a vertical offset from the ellipsoid.
  polarEllipsoid :: e,
    -- ^ The ellipsoid for the projection. Arguments passed to `toGrid` *must* use this ellipsoid.
    -- The type system cannot verify this for `LocalEllipsoid`.
  gridScale :: Double,
    -- ^ The scaling factor applied at the pole. This balances the distortion between the center
    -- and edges of the projection.

  -- Remaining elements are memoised parameters computed from the ellipsoid.
  gridA, gridB, gridC, gridD, gridC0 :: !Double
} deriving (Show)

instance (Eq e) => Eq (PolarStereographic e) where
  g1 == g2 =
    trueOrigin g1 == trueOrigin g2 &&
    falseOrigin g1 == falseOrigin g2 &&
    polarEllipsoid g1 == polarEllipsoid g2 &&
    gridScale g1 == gridScale g2

instance (Ellipsoid e) => GridClass (PolarStereographic e) e where
  fromGrid p = Geodetic lat long (altGP p) (polarEllipsoid gb)
    where
      gridZero = GridPoint 0 0 0 gb
      gb = gridBasis p
      p' = gridZero `gridOffset` (falseOrigin gb `applyOffset` p)
      radius = offsetDistance p'
      isoColat = 2 * atan (radius / (gridScale gb * gridC0 gb))
      isoLat = pi/2 - isoColat
      lat1 = isoLat +
        gridA gb * sin (2*isoLat) +
        gridB gb * sin (4*isoLat) +
        gridC gb * sin (6*isoLat) +
        gridD gb * sin (8*isoLat)
      lat = case trueOrigin gb of
        NorthPole -> lat1
        SouthPole -> negate lat1
      long = case trueOrigin gb of
        NorthPole -> offsetBearing p' { deltaNorth = negate $ deltaNorth p'}
        SouthPole -> offsetBearing p'

  toGrid r geo = offsetNegate (falseOrigin r) `applyOffset` GridPoint east north 0 r
    where
      absLat = abs $ latitude geo
      e = sqrt (eccentricity2 $ polarEllipsoid r)
      eSinLat = e * sin absLat
      tz2 = ((1 + eSinLat)/(1-eSinLat))**(e/2) * tan (pi/4 - absLat / 2)
      radius = gridScale r * gridC0 r * tz2
      north = case trueOrigin r of
        NorthPole -> negate $ radius * cos (longitude geo)
        SouthPole -> radius * cos (longitude geo)
      east = radius * sin (longitude geo)
  gridEllipsoid = polarEllipsoid


mkGridPolarStereographic :: (Ellipsoid e) =>
  Pole          -- ^ True origin at north or south pole.
  -> e          -- ^ The ellipsoid used for the projection.
  -> GridOffset -- ^ Vector from true origin to the false origin.
  -> Double     -- ^ Scale factor.
  -> PolarStereographic e
mkGridPolarStereographic pole ellip offset scale =
  PolarStereographic {
    trueOrigin = pole,
    falseOrigin = offset,
    polarEllipsoid = ellip,
    gridScale = scale,
    gridA = e2/2 + (5/24)*e4 + e6/12 + (13/360)*e8,
    gridB = (7/48)*e4 + (29/240)*e6 + (811/11520)*e8,
    gridC = (7/120)*e6 + (81/1120)*e8,
    gridD = (4279/161280)*e8,
    gridC0 = (2 * majorRadius ellip / sqrt (1 - e2)) * ((1-e1)/(1+e1))**(e1/2)
  }
  where
    e1 = sqrt $ eccentricity2 ellip
    e2 = eccentricity2 ellip
    e4 = e2^_2
    e6 = e2^_3
    e8 = e2^_4


-- | The Universal Polar Stereographic (UPS) grids for north and south poles.
type UpsGrid = PolarStereographic WGS84


-- | UPS grid for the North Pole.
upsNorth :: UpsGrid
upsNorth = mkGridPolarStereographic
  NorthPole
  WGS84
  (GridOffset { deltaEast = -(2000 * kilometer), deltaNorth = -(2000 * kilometer), deltaAltitude = 0 })
  0.994  -- Scale factor


-- | UPS grid for the South Pole.
upsSouth :: UpsGrid
upsSouth = mkGridPolarStereographic
  SouthPole
  WGS84
  (GridOffset { deltaEast = -(2000 * kilometer), deltaNorth = -(2000 * kilometer), deltaAltitude = 0 })
  0.994   -- Scale factor


-- | Convert a grid reference into a UPS grid location.
--
-- There doesn't appear to be any conventional representation for polar grid references,
-- so this is an attempt to cover as many bases as possible. It takes an Easting followed
-- by a Northing with spaces in between. Both can have optional units of m or km,
-- and be optionally followed by an \"N\" or \"E\" as appropriate.
--
-- The choice of pole is provided in an extra argument rather than within the string because humans
-- will normally assume this from the context and so not provide it.
--
-- If the string cannot be parsed then one or more error messages are returned.
fromUpsGridReference :: Pole -> String -> Either [String] (GridPoint UpsGrid)
fromUpsGridReference pole str = case parse (parseUpsGridReference pole) "" str of
    Left err -> Left $ lines $ showErrorMessages
      "or" "unknown parse error" "expecting" "unexpected" "end of input"
      (errorMessages err)
    Right r -> Right r


parseUpsGridReference :: Stream s m Char => Pole -> ParsecT s u m (GridPoint UpsGrid)
parseUpsGridReference pole = do
      spaces1
      (eastings1, eastUnit) <- readDistance
      spaces
      optional (oneOf "Ee" <?> "E")
      spaces1
      (northings1, northUnit) <- readDistance
      unless (eastUnit == northUnit) $ fail "Northings and Eastings units don't match."
      spaces1
      optional (oneOf "Nn" <?> "N")
      spaces1
      eof
      return $ GridPoint eastings1 northings1 0 $
        case pole of
          NorthPole -> upsNorth
          SouthPole -> upsSouth
  where
    readDistance = do  -- Returns (distance in meters, unit from input)
      digits <- many1 (digit <|> char '.' <?> "number")
      spaces1
      when (length digits > 20) $ fail "Too many digits."
      (multiplier, unit) <- do
        unit <- optionMaybe (string1' "m" <|> string1' "km" <?> "units (m or km)")
        case unit of
          Just "km" -> pure (1000, Just GridKilometers)
          Just _    -> pure (1, Just GridMeters)
          Nothing   -> pure (1, Nothing)
      case readMaybe digits of
        Just d -> return (d * multiplier, unit)
        Nothing -> fail $ "Cannot read number: " <> digits
    string1' target = try $ do  -- Case-insensitive version of string'
      cs <- count (length target) anyToken
      if map toLower target == map toLower cs then return cs else unexpected cs
    spaces1 = void $ many (char ' ' <?> "space")  -- Other white space not permitted.


toUpsGridReference :: 
  Maybe GridUnit  -- ^ Include explicit units in the output. @Nothing@ means meters without units.
  -> Bool  -- ^ Include \"E\" and \"N\" in the output.
  -> Int  -- ^ Digits of resolution. 0 = 1m resolution, 1 = 10m, 2 = 100m etc. (-2) = 1cm.
  -> GridPoint UpsGrid
  -> String
toUpsGridReference unit letters res gp =
  dist (eastings gp)  <> (if letters then "E " else " ") <>
    dist (northings gp) <> (if letters then "N" else "")
  where
    res1 :: Double
    res1 = 10 ** fromIntegral res   -- Resolution in meters.
    floorRes :: Double -> Double
    floorRes d = res1 * fromIntegral (floor (d/res1) :: Integer)
    dist d = case unit of
      Nothing            -> printf "%.*f" (-res) $ floorRes d
      Just GridMeters     -> printf "%.*fm" (-res) $ floorRes d
      Just GridKilometers -> printf "%.*fkm" (3-res) $ floorRes d / 1000
