# Changelog

Main changes to this code/ project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Center output nifti files of AIBL
- Remove FSL library depency for OASIS-to-bids conversion

### Changed
- Fix minor typos in documentation (#205)
- Replace exception by warning when CAPs folder not recognized
- In AIBL-to-bids, center output nifti files of AIBL

### Fixed

- Add SPM12 dependency check for pet_surface pipelines

## Clinica 0.4.0

### Added

- pet-linear pipeline: spatial normalization to the MNI space and intensity normalization of PET images
- pet-surface-longitudinal pipeline: Surface-based longitudinal processing of PET images 
- check-missing-processing tool allows creating a TSV file containing information about the pipelines executed into a specific CAPS folder
- Conversion information is added once the converter is run to facilitate traceability.
- Add new keywords available in ADNI3 to the adni-2-bids converter

### Changed

- Harmonize output message display when running a pipeline.
- Automatically ignore an image to process if it is found in the CAPS folder.
- 


### Fixed

- Fix a minor bug when the pet-volume pipeline performs PVC.
- Code source was completely reformatted using the Black code style.
- Documentation for the project is now versioned (versions from 0.3.8 are publicly available).
- Functions used for multiple pipelines are now mutualized (e.g container_from_filename function).
- f-strings are used massively.

## Clinica 0.3.8

### Added

- Add option to run deeplearning-prepare-data in output of t1-extension pipeline
  and custom pipelines (PR #150).
- Add Build and publish documentation with CI (PR #146).
- Add CHANGELOG.md file

### Changed

- Harmonize PET tracers handling (ML/DL/Stats) (PR #137).
- Behaviour of ADNI converter: some minor bugs and updated wrt ADNI3. E.g., the
  field age_bl was added. (PR #139, #138, #140, #142)

### Fixed

- Add DataDictionary_NIFD_.xlsx file when using NIFD2BIDS.


## Clinica 0.3.7 - FreeSurfer-Longitudinal

### Changes

#### Clinica Core:

- [New] Remove CAT12 from dependencies
- [New] Add checksum for volume atlases
- [Fix] Remove duplicated lines in AICHA ROI file
- [New] Integrate clinica.wiki repository into clinica repository: New versions
  of Clinica will now have their own version of documentation

#### Pipelines:

- [New] `t1-freesurfer-longitudinal` pipeline: FreeSurfer-based longitudinal
  processing of T1-weighted MR images [[Reuter et al.,
  2012](http://dx.doi.org/10.1016/j.neuroimage.2012.02.084)]. More info in the
  wiki: https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer_Longitudinal/
- [Change] The `fmri-preprocessing` pipeline is removed from the Clinica
  software as we will not actively maintain it. It is now in a separate
  repository: https://github.com/aramis-lab/clinica_pipeline_fmri_preprocessing

#### Converters:

- [Enh] Improve how dependencies are checked for converters
- [Fix] Add diagnosis conversion for ADNI3
- [Fix] Avoid creation of sessions `ses-V01` in `*_sessions.tsv` files


##Clinica 0.3.6 - Pip

### Changes

#### Clinica Core:

- [Change] Pip is the main way to install Clinica. Conda packages are not
  available for the new versions.
- [Update] Set the minimal version of Python to 3.7.
- [Fix] Remove non-breaking spaces.

#### Pipelines:

- [New] Display failed image(s) when running `t1-linear` pipeline.

#### Converters:

- [Fix] The `aibl-2-bids` converter now handles new version of clinical data.
- [Update] The `oasis-2-bids` converter now uses NiBabel instead of FreeSurfer
  to convert OASIS dataset.


## Clinica 0.3.5 - DL-Prepare-Data

### Changes

#### Clinica Core:

- [CI] Improve Jenkins configuration (Automatic generation of testing reports
  in order to be displayed in the CI interface; Recreate Python environment if
  `requirements.txt` changes)

#### Pipelines:

- [New] `deeplearning-prepare-data` pipeline: Prepare input data for deep
  learning with PyTorch. More info on the Wiki:
  https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DeepLearning_PrepareData/
- [Change] `t1-linear` pipeline now crops image on default. If
  `--uncropped_image` is added to the command line, the image is not cropped.
- [Change] Refactor machine learning modules. Main changes involve use of
  CamelCase convention for classes and parameters used dictionaries.


##Clinica 0.3.4 - Bugfixes

### Changes

#### Clinica Core:

- [Improvement] Remove Clinica dependencies while updating and unfreezing some
  of them

#### Pipelines:

- [Enh] Improve how template files are downloaded for `t1-linear` and
  `statistics-volume` pipelines

