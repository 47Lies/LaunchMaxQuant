# Change Log
All notable changes to MaxQuant will be documented in this file.

## [1.6.12.0] - 2020-03-13
 - Isobaric matching between runs
 - Isobaric PSM-level weighted ration normalization
 - Bugfixes

## [1.6.11.0] - 2020-03-01
 - TIMS performance improvements
 - Bugfixes
 
## [1.6.10.43] - 2019-10-17
- Fixed [MaxQuant-348](https://maxquant.myjetbrains.com/youtrack/issue/MaxQuant-348)

## [1.6.10.0] - 2019-10-17
 - Fixed bug with retention time calibration [MaxQuant-318](https://maxquant.myjetbrains.com/youtrack/issue/MaxQuant-318)

## [1.6.9.0] - 2019-10-14
 - Command line improvements
   - Folders of raw and fasta files can be changed in mqpar file. 
   - Job indices for partial processing in MaxQuantCmd are now 1-based as they always were in the GUI.
 - Bug fixes, performance improvements

## [1.6.8.0] - 2019-09-27
### Changed
 - Bug fixes
 - TIMS: Performance improvements

## [1.6.7.0] - 2019-07-10
### Changed
 - PTMs in sequence strings are specified by their full name in round brackets, 
 not bei their abbreviations any more.
### Added 
 - Sciex Linux support
### Removed
 - Matching from and to options removed
### Fixed 
 - TIMS bug fixing and performance optimisations

## [1.6.6.0] - 2019-05-13
### Fixed
 - TIMS data: fixed NaN values in evidence table
 - TIMS data: Performance optimization in MSMS preparation step

## [1.6.5.0] - 2019-01-27
### Added
 - Bruker raw files: Visual C++ Redistributable package. Installation check and added vcredist folder
## Fixed
 - Bruker raw files: Linux fixes for RawFile plugins

## [1.6.4.1] - 2018-12-13
### Fixed
 - Support for the mzTab file format at TIMS data sets

## [1.6.4.0] - 2018-12-10
### Added 
 - TIMS: Added Match between runs
 - TIMS: Added ion mobility columns to the evidence.txt table
 - TIMS: Performance fixes
 - new dependent peptides output table: Enhances the analysis of dependent (modified) peptides with comprehensive qualitative and quantitative information. Important: Always use the new dependent peptide search in conjunction with match between runs -> match unidentifed features and experiment specifications to enable optimal information gain.

## [1.6.3.4] - 2018-11-08
### Changed
 - TIMS/TMT stability issues fixed

## [1.6.3.3] - 2018-10-24
### Added
 - Support for the mzTab file format

### Changed
 - Parameter for isobaric labels works different now. Correction factors are not configured with the modification any more, but they are specified in the parameter. This allows for analyzing multiple batches with different correction factors together.
 - The channel index for isobaric labels is now one-based. It used to be zero-based.
 - Mass deficit is only calculated for peptides lighter than 2000Da
 - Stability improvements for large number of threads (>=120)

## [1.6.2.10] - 2018-07-26
### Fixed
 - TIMS-TOF workflow stable

## [1.6.2.6] - 2018-07-04
### Added
 - "Checking fasta files" step at the beginning
 - "Variation search" is fully enabled


## [1.6.2.3] - 2018-06-11
### Fixed
- MaxQuant-247 MQ is crashing if more than one parameter group is specified

## [1.6.2.2] - 2018-06-10
### Added
- TMT support for timsTOF

### Fixed
- MaxQuant-246 some searches crash after main search (error reading search results)
- MaxQuant-243 MaxQuant 1.6.2.1 creates invalid XML in mqpar file (spaces in encoding specification) resolved

 ## [1.6.2.1] - 2018-06-01
 ### Added
 - A new commandline interface for AndromedaCmd.exe

 ### Fixed
  - "mqpar.xml conversion error" is resolved (MaxQuant-242)
  - Bruker TIMS files are loadable

## [1.6.2.0] - 2018-05-25
### Added
- Fasta file input parameter: parse rules and taxonomy can now be specified
  in the parameter tabs
- Bruker TIMS Tof PASEF is supported  
- Additional parameters for command line call

### Fixed
 - "Could not find file tmpPeptideFiles" is resolved (MaxQuant-192)

 ## [1.6.1.0] - 2017-12-08
 ### Added

 ### Fixed
  - MaxQuant running smoothly on Ubuntu Linux
  - MaxQuant-222 Fragmentation type always identified as CID

 ## [1.6.0.16] - 2017-08-29
 ### Fixed
  - Added missing button icons/tool tips to Andromeda configuration parse rule testing GUI [[MaxQuant-197](https://maxquant.myjetbrains.com/youtrack/issue/MaxQuant-197)].
  - MaxQuant-191 Export of MS/MS spectrum to pdf in MaxQuant 1.6.0.1 viewer
  - Equal results between mono and .NET framework

## [1.6.0.13] - 2017-08-11
### Added
  - Raw data access API is in open source code on GitHub
  - MzXml raw data reader is in open source code on GitHub

### Fixed
  - 3D spectrum viewer is working again.
  - MS3 based TMT quantification is working again.
  - Special characters as german umlauts can be used in file paths for Bruker.

## [1.5.8.0] - 2017-03-13
### Added
 - Winforms based GUI.

## [1.5.7.4] - 2017-01-27
### Fixed
 - Group-specific parameters are loaded properly from the mqpar.xml file.
 - Command line version works again.

## [1.5.7.0] - 2017-01-12
### Fixed
 - Support for NeuCode labeling added.

## [1.5.6.5] - 2016-11-22
### Fixed
 - Removed possibility to load zipped fasta files since it is not working.
 - Specifying fractions does not lead to a crash any more upon clicking starrt.

## [1.5.6.0] - 2016-10-31
### Added
 - This version corresponds to the feature set described in the the publication Tyanova et al. (2016) Nature Protocols 11, 2301-19.

