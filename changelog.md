Release Notes
-------------

Version 0.0.2: Tided up cabal file and removed spurious dependency on Parsec.

Version 0.0.3: Updated for Haskell Platform 2014.2.0.0 and GHC 7.8.3. Fixed
   some minor documentation issues.

Version 0.0.4: Updated for Dimensional 1.0.

Version 0.0.5: Fixed bug in Monoid instance for Helmert. Created Semigroup
   instance for Helmert.

Version 0.0.6: Prevent attempted building on GHC 7.8 (it doesn't work)
   and fix the build on 7.10 with a conditional semigroups dependency

Version 0.1.0: Updated for Dimensional 1.3 and GHC 8.6.

Version 0.1.1: Fixed bug #15: for a point p, groundDistance p p returned NaN

Version 0.1.2: Fixed bugs #16 and #17: Unicode PRIME and DOUBLE PRIME now allowed in
   position strings, and the degree symbol is allowed for decimal degrees.

Version 1.0.0: Removed dependency on Dimensional library. This is a breaking change:
   hence the major version bump. Also fixed bug #18 (and #19).

## Version 1.1.0

* Dropped a redundant `Ellipsoid` constraint on `antipode`.
* Added functions `showGeodeticLatLong`, `showGeodeticSignedDecimal`,
  `showGeodeticNSEWDecimal`, and `showGeodeticDDDMMSS`.
