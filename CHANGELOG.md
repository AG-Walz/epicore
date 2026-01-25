# AG-Walz/epicore: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.8 - [dev]

### Added
- Enable epicore for cohorts.
- The flag --strict, enabling a strict peptide grouping mode.
- Add quality control plots to output.

### Changed 
- Sort peptides with descending end positions and ascending start positions.
- Split peptide groups at landscape minima in default mode.
- Add all peptides completely included in a group to it.
- Change pandas to polars for input parsing.
- Use multiprocessing.

### Fixed
- Include modification information in output.
- Fix mapping of consensus sequences to input peptides.
- Fix peptide sorting for repetitive regions.

### Removed
- Support for intensity column.


## v0.1.7 - 2025/10/14

### Changed
- Remove index column from pep_cores_mapping.tsv.
- Set required python version from >=3.12 to >=3.10.

## v0.1.6 - 2025/07/30

### Changed
- Change output format from csv to tsv.
- Increase the landscape value for all non-repetitive peptides, not just at the first occurrence.
- Plot core epitope lengths instead of the length of the entire consensus epitope. 
- Change the calculation of the total intensity of the evidence file by counting each peptide intensity only once.
- Adjust column output order to column input order.

### Fixed
- Remove modification of peptides for peptide length histogram.
- Adjust the epicore intensity plot range to make sure that all epitopes are included.
- Remove redundant intensity information. 
- Fix aggregation method for peptides with the same sequence but different properties.
- Fix mapping of cores and input peptides for peptide sequences that occur multiple times in the evidence file.
- Fix proteome_occurrence column for peptides that are part of a repetitive region.

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