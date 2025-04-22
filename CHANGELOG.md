# Changelog

## [1.2.0] - 2025-04-22
### Added
- Cross-platform support improved for Linux and Mac
- Consistent version handling across the package
- Expanded documentation in README

### Changed
- Updated author contact information
- Improved dependency specifications in setup.py and conda files
- Switched repository from GitHub to Bitbucket

### Fixed
- Version synchronization between version.py and setup.py
- Package metadata in setup.py to include template and asset files
- Added missing dependencies (scipy, numpy) to requirements

## [1.1.0] - 2025-02-28
### Added
- Data type-specific quality assessment
  - Specialized thresholds for bacterial, viral, and metagenomic data
  - Added `get_thresholds()` function to select appropriate thresholds by data type
  - Improved reports with dataset type indicators

### Enhanced
- QC Report Generation
  - Fixed PDF layout issues with first page
  - Added dataset type indicator in reports
  - Tailored metrics display based on data type
  - Customized assessment text for different data types

### Changed
- Modularized threshold configuration
  - Separated bacterial, viral, and metagenomic thresholds

## [1.0.1] - 2025-02-20
- Fixed package structure to resolve entry point issues
- Added proper Python package structure with bin/__init__.py
- Updated setup.py to properly handle script installation
- Added KMA dependency information in documentation
- Fixed conda configuration

## [1.0.0] - 2025-02-15
- Initial project structure for cgeqc