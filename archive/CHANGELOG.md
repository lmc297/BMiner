# BMiner CHANGELOG

All noteable changes to BMiner will be documented in this file

## [2.1.0] 2018-09-17
### Added
- Added ANIb tab and parsing to construct a bar graph using ANIb results from BTyper version 2.3.0; all BTyper results from BTyper versions 2.0.0 and up can still be used with this version of BMiner

### Changed
-Updated README to reflect addition of ANIb results parsing

## [2.0.2] 2018-01-17
### Changed
- Added simplify=FALSE to ensure sapply returns a list of virulence or AMR genes when parsing virulence and AMR sections of files (lines 279 and 281; in a case where final results files were provided that all possessed identical virulence and AMR profiles, sapply returned a matrix, and BMiner could not parse the files)

## [2.0.1] 2017-11-14
### Changed
- Moved older versions of BMiner (1.0.0 and 2.0.0) to their own directories within the archive directory so they can be run easily using shiny's runGitHub command with the subdir option in R 
- Updated README to clarify how to launch the current version of BMiner, as well as the older versions
- Updated README to clarify that results files created using BTyper version 1.0.0 need to be run with BMiner version 1.0.0
- Updated README to include the library(shiny) command as a step in launching BMiner

## [2.0.1] 2017-08-17
### Added
- Added ggtree as a dependency

### Changed
- Fixed bugs in UI to make app compatible with newest version of R (3.4.1 Single Candle), R Studio, and required packages
- Changed README to reflect ggtree's addition as a dependency

## [2.0.0] 2017-06-29

### Added
- Antimicrobial resistance gene analysis tab (barplot, NMDS, PCA, presence/absence matrix) 

### Changed
- Renamed "Samples" to "Genomes" for all plots
- Updated README to reflect addition of antimicrobial resistance gene detection function; updated disclaimer

## [1.0.0] 2017-05-14

### Added
- BMiner initial commit (virulence, panC clade, rpoB allelic type, MLST, and 16S aggregation)
