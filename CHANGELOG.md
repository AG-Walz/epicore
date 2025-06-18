# AG-Walz/epicore: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.6 - [dev]

### Changed
- Change output format from csv to tsv.
- Increase the landscape value for all non-repetitive peptides, not just at the first occurrence.

### Fixed
- Remove modification of peptides for peptide length histogram.
- Adjust the epicore intensity plot range to make sure that all epitopes are included.
- Remove redundant intensity information. 
- Fix aggregation method for peptides with the same sequence but different properties.

## v0.1.5 - 2025/06/16

### Changed
- Replace path in html file with svg content.

## v0.1.4 - 2025/05/28

### Added
- Add optional --html flag for html images.

### Fixed 
- Fix bin range for number of peptides mapped to one core in plot.
- Fix landscape computation.
- Fix core epitope end position computation.

### Changed
- Move location of log file from current directory to result directory. 
- Change output format.
- Change input parameters file to command line parameters. 
- Add proteome occurrence and intensity columns to output.


## v0.1.3 - 2025/04/24

### Added
- Add CHANGELOG.md

### Fixed
- Fix handling of provided start and end position of peptides in the proteome sequence. 
- Fix mapping error of core epitopes and input peptides
- Decrease runtime

### Changed
- Update epicore commands in README.
- Update conda installation command in README.   

### Removed
- Workflow for automatic version update for each release

## v0.1.2 - 2025/04/07

## v0.1.1 - 2025/04/07

### Changed
- Update github workflow to automatically update the version to the current release tag. 

### Removed
- setup.py
- requirements.txt

## v0.1.0 - 2025/03/27

- Initial release of epicore