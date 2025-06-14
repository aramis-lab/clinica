<!-- markdownlint-disable MD024 MD026 -->
# Changelog

Main changes to this code/ project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Clinica 0.10.0

### Fixed

- [IOTools] Fix broken import of `bids.dataset` module (PR [#1484](https://github.com/aramis-lab/clinica/issues/1484))
- [IOTools] `Center-Nifti` now runs by default on T1 as it was expected (PR [#1418](https://github.com/aramis-lab/clinica/pull/1418))
- [DWI] Fix broken import of `bids.dataset` module (PR[#1523](https://github.com/aramis-lab/clinica/pull/1523))
- [Converters] GENFI : Fix validation of optional path to clinical data (PR [#1513](https://github.com/aramis-lab/clinica/pull/1513))
- [Pipelines] Fix SPM standalone version check (PR [#1488](https://github.com/aramis-lab/clinica/pull/1488))

### Enhanced

- [Pipelines] Search for SPM standalone version is not hardcoded to SPM12 anymore (PR [#1490](https://github.com/aramis-lab/clinica/pull/1490))
- Improve architecture for BIDS and CAPS datasets logic (PR [#1453](https://github.com/aramis-lab/clinica/pull/1453))

### Added

- [IOTools] When using `Center-Nifti` it is now possible to choose the centering threshold (PR [#1433](https://github.com/aramis-lab/clinica/pull/1433))
- [DOC] Warning about SPM standalone version (PR [#1496](https://github.com/aramis-lab/clinica/pull/1496))
- [T1Volume] Skip processed visits (PR [#1403](https://github.com/aramis-lab/clinica/pull/1403))
- [AIBL-to-BIDS] Now has a json constructor based on dicoms to retrieve information required to fit BIDS specifications (PR [1500](https://github.com/aramis-lab/clinica/pull/1500))

### Removed

- [Pipelines] Remove interactivity when centering images is required by the pipeline (PET-Volume, T1-Volume, PET-Surface, T1-Freesurfer) (PR [#1419](https://github.com/aramis-lab/clinica/pull/1419))

## Clinica 0.9.4

### Fixed

- [PET] Fix query pattern when a tracer is specified but not the reconstruction method (see PR [#1470](https://github.com/aramis-lab/clinica/pull/1470))
- [Converters] Fix type mismatch in `clinica.iotools.bids_utils` (see PR [#1376](https://github.com/aramis-lab/clinica/pull/1376))

## Clinica 0.9.3

### Fixed

- [PETLinear] Fix the SUVR normalization step (see PR [#1348](https://github.com/aramis-lab/clinica/pull/1348))

## Clinica 0.9.2

### Fixed

- [Converters] Fix KeyError "APGEN" in adni-to-bids (see PR [#1342](https://github.com/aramis-lab/clinica/pull/1342))
- [Converters] Fix function `load_clinical_csv` to avoid hidden files (see PR [#1343](https://github.com/aramis-lab/clinica/pull/1343))
- [Converters] Fix logging in oasis-to-bids and aibl-to-bids (see PR [#1340](https://github.com/aramis-lab/clinica/pull/1340))
- [Converters] Fix method `BIDSSubjectID.to_study_id()` (see PR [#1335](https://github.com/aramis-lab/clinica/pull/1335))
- [T1Linear] Fix missing import in utility function (see PR [#1329](https://github.com/aramis-lab/clinica/pull/1329))
- [Converters] Fix oasis-to-bids to work with new clinical data (see PR [#1321](https://github.com/aramis-lab/clinica/pull/1321))
- [Converters] Fix function `define_participants` in ixi-to-bids, which is now deterministic (see PR [#1318](https://github.com/aramis-lab/clinica/pull/1318))
- [Converters] Fix Oasis1 `BIDSSubjectID` method (see PR [#1312](https://github.com/aramis-lab/clinica/pull/1312))

## Clinica 0.9.1

### Fixed

- Fix broken install of release `0.9.0` (see PR [#1304](https://github.com/aramis-lab/clinica/pull/1304))

## Clinica 0.9.0

### Added

- [ixi-to-bids] There is a new converter to convert IXI to BIDS (see PR [#1239](https://github.com/aramis-lab/clinica/pull/1239))
- [iotools] There is a new tool `clinica iotools describe` to nicely display metadata in `dataset_description.json` files in the console (see PR [#1287](https://github.com/aramis-lab/clinica/pull/1287))
- [t1-linear] It is now possible to use `t1-linear` with ANTsPy (see PR [#1244](https://github.com/aramis-lab/clinica/pull/1244))
- [ADNI2BIDS] The converter now support FMAP (see PR [#1119](https://github.com/aramis-lab/clinica/pull/1119))
- [Doc] Some pages have been refactored to improve readability and make information easier to find (see PR [#1284](https://github.com/aramis-lab/clinica/pull/1284), [#1288](https://github.com/aramis-lab/clinica/pull/1288), and [#1295](https://github.com/aramis-lab/clinica/pull/1295))
- [Doc] The online documentation now has a glossary page (see PR [#1110](https://github.com/aramis-lab/clinica/pull/1110))

### Breaking changes

- BIDS and CAPS datasets must have a `dataset_description.json` file at the root (see PR [#1127](https://github.com/aramis-lab/clinica/pull/1127) and PR [#1158](https://github.com/aramis-lab/clinica/pull/1158))

### Enhanced

- [Converters] Converters now implement the same API which make them easier to use from Python (see PR [#1140](https://github.com/aramis-lab/clinica/pull/1140))

### Fixed

- [Converters] A few issues with logging have been fixed (see PR [#1181](https://github.com/aramis-lab/clinica/pull/1181))
- [t1-linear] Fix the `crop_nifti` function which was doing unnecessary resampling (see PR [#1215](https://github.com/aramis-lab/clinica/pull/1215))

## Clinica 0.8.4

### Fixed

- [Pipelines] Fix 'SPM not in matlab path' error when using SPM12 with Matlab (see PR [#1261](https://github.com/aramis-lab/clinica/pull/1261))

## Clinica 0.8.3

### Fixed

- [IOTools] Fix `merge-tsv` when providing a `T1Freesurfer` CAPS directory (see PR [#1240](https://github.com/aramis-lab/clinica/pull/1240))

## Clinica 0.8.2

### Fixed

- [ADNI2BIDS] Fix wrong preprocessing sequence in FDG PET Uniform (see PR [#1159](https://github.com/aramis-lab/clinica/pull/1159))
- [NIFD2BIDS] Fix check on type in pandas dataframe (see PR [#1230](https://github.com/aramis-lab/clinica/pull/1227))
- [DWI] `DWIPreprocessingUsingT1` pipeline skips problematic images instead of crashing (see PR [#1169](https://github.com/aramis-lab/clinica/pull/1169))
- Fix compatibility issues with matplotlib 3.9 (see PR [#1185](https://github.com/aramis-lab/clinica/pull/1185))

## Clinica 0.8.1

### Fixed

- [ADNI2BIDS] Fix `DXSUM_PDXCONV_ADNIALL` unknown clinical file for adni-to-bids converter (see PR [#1144](https://github.com/aramis-lab/clinica/pull/1144))
- [DWI] Fix DWI preprocessing using T1 rename to caps node (see PR [#1146](https://github.com/aramis-lab/clinica/pull/1146))

## Clinica 0.8.0

### Added

- Python 3.12 is officially supported (see PR [#1065](https://github.com/aramis-lab/clinica/pull/1065))

### Breaking changes

- [DWI] Improve BIDS compliance for DWI preprocessing pipelines (see PR [#1050](https://github.com/aramis-lab/clinica/pull/1050))
- Python 3.8 is not supported anymore (see PR [#1065](https://github.com/aramis-lab/clinica/pull/1065))

### Enhanced

- [ADNI2BIDS] The converter now takes fmri multiband into account (see PR [#1041](https://github.com/aramis-lab/clinica/pull/1041))
- [DWI] Set the random seed in pipelines using ANTs (requires having ANTs >= 2.5) (see PR [#1071](https://github.com/aramis-lab/clinica/pull/1071))
- [Pipelines] Using SPM standalone only requires to set `$SPMSTANDALONE_HOME` and `$MCR_HOME` (see PR [#1105](https://github.com/aramis-lab/clinica/pull/1105))
- [DOC] The Third-party documentation page has been improved (see PR [#1099](https://github.com/aramis-lab/clinica/pull/1099))
- [DOC] The documentation page for `StatisticsSurface` has been updated (see PR [#1100](https://github.com/aramis-lab/clinica/pull/1100))

### Fixed

- [ADNI2BIDS] Ignore real and imaginary scans outputted by dcm2niix (see PR [#1029](https://github.com/aramis-lab/clinica/pull/1029))
- [ADNI2BIDS] Ignore ADC DWI images outputted by dcm2niix (see PR [#1061](https://github.com/aramis-lab/clinica/pull/1061))

## Clinica 0.7.7

### Enhanced

- [ADNI2BIDS] Enable converter to read CSV files with new naming convention from ADNI (see PR #1016)
- [Converters] Expose `n_procs` option to converters using multiprocessing (see PR #1009)
- [GENFI2BIDS] Enable converter to extract more clinical data (see PR #1005)

### Fixed

- [IOTools] Fix regression in `merge-tsv` command (see PR #1013 and #1019)
- [GENFI2BIDS] Fix issue with link and readme data not found by the converter (see PR #1015)
- [PETLinear] Fix the CLI of PETLinear due to previously added option (see PR #1004)
- [T1Linear] Fix bad session labels in output file names (see PR #1000)
- [StatisticsSurface] Some fixes to the pipeline (see PR #840)

## Clinica 0.7.6

### Enhanced

- [Pipelines] The clinica file reader now supports the `run` entity and picks the latest run available by default (see PR #943)
- [ADNI2BIDS] The ADNI2BIDS converter now handles PET uniform data (see PR #703)
- [PETLinear] Improve the BIDS compliance of the CAPS output (see PR #935)
- [T1Linear & FLAIRLinear] Improve the BIDS compliance of the CAPS output (see PR #940)

### Fixed

- [DWI] Fix regression in pipeline DWIPreprocessingUsingT1 which was not generating bvec file since release `0.7.4` (see PR #967)
- [ADNI2BIDS] Fix broken assertion in the metadata reader for derived images (see PR #957)
- [ADNI2BIDS] Fix BIDS violations in converter's output (see PR #959)

## Clinica 0.7.5

### Enhanced

- [Converters] Improvements to GENFI to BIDS converter (see PR #909)
- [Converters] Improvements to Oasis3 to BIDS converter (see PR #911)
- [Pipelines] Add the possibility to specify a random seed to be used with ANTs (see PR #916)

### Fixed

- [Converters] GENFI to BIDS various fixes (see PRs #912, #923, #927, #930, and #934)
- [Converters] Fix ADNI to BIDS issue with unsupported visit code "uns1" (see PRs #914 and #920)
- [MachineLearning] Fix bug in CLI of pipeline `machinelearning-prepare-spatial-svm` (see PR #910)

## Clinica 0.7.4

### Enhanced

- [DWI] Read total readout time and phase encoding direction from JSON files
- [IOTools] Session numbers are now represented using three digits instead of two
- [PET] Handle 4D images with dummy fourth dimension
- [Dependencies] Upgrade Poetry lock files to version 2

### Added

- [Converters] New GENFI to BIDS converter
- Support for Python 3.11

### Fixed

- [BUG] Fixed issue with NetworkX 3
- [BUG] Fix run statistics-volume-correction assertion error
- [DOC] Fix command in statistics volume documentation
- [IOTools] Fix `create_scan_dict` row extraction condition

## Clinica 0.7.3

### Enhanced

- [CI] Add caching support for unit tests
- [CI] Refactor testing tools
- [Dependencies] Bump lxml from 4.9.0 to 4.9.1
- [Dependencies] Upgrade joblib to 1.2.0
- [Dependencies] build: Install nipype up to version 1.8.2
- [SurfStat] Pure python implementation
- [IOTools] Fix warnings in merge-tsv
- [Adni2BIDS] Deal with new data from ADNI3
- [DWIPreprocessingUsingT1] Optimized disk usage of Pipeline DWIPreprocessingUsingT1
- [IOTools] Allow setting a custom logging directory via environment variable
- [IOTools] Center all modalities if no modality is specified
- [Pipelines] Report uncompliant BIDS subjects

### Added

- [Converters] Add support for BIDS Readme 
- [IOTools] Extend the create-subjects-sessions iotool to CAPS directories
- [IOTools] Add pet-linear to checks for missing processing

### Fixed

- [UKB2BIDS] Add error if data is not found or filtered 
- [DWIPreprocessingUsingT1] Add missing `out_file` parameter to DWIBiasCorrect
- [Converters] UKB2BIDS drop directories labeled as unusable
- [Adni2BIDS] Handle empty lines in `create_subs_sess_list`
- [IOTools] Fix `vox_to_world_space_method_1`

## Clinica 0.7.2

### Fixed

- [Pipelines] Fix bug introduced in previous version with the use of the gunzip interface
- [DWIConnectome] Use ConstrainedSphericalDeconvolution instead of buggy EstimateFOD 

### Enhanced 

- [Adni2Bids] Add compatibility for edge cases introduced in Adni3

## Clinica 0.7.1

### Added
- [Doc] add ukbiobank documentation
- [DWIConnectome] Fetch meta data directly from MRtrix github repository

### Changed
- [Core] Enable parallelization when grabbing files

### Fixed
- [Converters] Fix several warnings


## Clinica 0.7.0

### Added
- [flair-linear] new pipeline to affinely align FLAIR images to the MNI space
- [Ukbiobank] new converter to modify T1W/T2/DWI/SWI/tfmri/rsfMRI UKBiobank data into BIDS standard

## Clinica 0.6.0

### Changed
- [PET*]   Use `trc`instead of `acq` for BIDS compliance
- [Converters] Remove superfluous use of `acq` entity in filenames for BIDS compliance

### Added
- [adni-to-bids] allow extraction of metadata from xml
- [CI] Initiate use of unit tests

### Fixed
- [adni-to-bids] fix edge case for supporting `nan` session-ids

## Clinica 0.5.6

### Fixed
- [DWIPreprocessUsingT1] Updated call to antsApplyTransform
- [Utils] Replace deprecated call to pandas `append` by `concat`

### Changed

- Upgrade minimum Python version to 3.8 and upgrade dependencies
- Set BIDS version to 1.7.0 by default (overwritten for some converters)

## Clinica 0.5.5

### Fixed
- [`pet-linear`] fix bug in `pet-linear`which had the pipeline not terminate 


## Clinica 0.5.4

### Added

- [`merge-tsv`] Add `t1-freesurfer-longitudinal` and `dwi-dti` results

### Changed

- [`t1-freesurfer`] Enable t1-freesurfer to run with missing files
- [all converters] Normalize subprocess calls to `dcm2niix`

### Fixed

- [`OASIS3`/`NIFD`/`HABS`]] add data_description file to BIDS
- [`DWI-DTI`] Remove thresholding for DECFA
- [`adni-to-bids`]  Tighten check on `session-id` values
- [`adni-to-bids`] Fix bug related to multiple conversions


## Clinica 0.5.3

### Added

- [`t1-freesurfer`] Add option to t1-freesurfer to project the results of `recon-all` onto another atlas
- [`CI`] Use `poetry` for dependency management

### Changed

- [`CI`] Refactor non-regression tests for easier parallelization
- [`Atlas`] Update checksum to make pipelines compatible with `fsl 6.0.5`

### Fixed

- [`t1-volume*/pet*`] Add command line argument `yes` for turning interactivity off
- [`t1-volume-existing-template] Fix chained invocation
- [`t1-volume*/pet-volume*] Fix default value of `--smooth` parameter for click compatibility
- [`dwi-connectome`] Set `--n_tracks`'s type for click compatibility 
- [`dwi-preprocessing*`] Change type of `initrand` and `use_cuda` to bool 
- [`t1-freesurfer-longitudinal`] Fix broken pipeline due to typo in code
- [`Documentation`] Update OASIS3_to_bids instructions for conversion
- [`StatisticsSurface`] Fix type in `covariate` argument
- [`StatisticsVolume`] Fix bug in `feature` argument 

## Clinica 0.5.2

### Changed

- [`dwi-preprocessing*`] Rewrite of `dwi-preprocessing*` pipelines using FSL's `eddy` tool

### Removed

- [`deeplearning-prepare-data`] Migration of pipeline to ClinicaDL

### Fixed

- [`oasis3-to-bids/nifd-to-bids`] Change code for backward compatibility with pandas 1.1.x
- [T1-FreeSurfer/DWI] Remove Typing for compatibility with Nipype

## Clinica 0.5.1

### Added

- [`oasis3-to-bids`] Add converter
- [GitHub] Add citation file

### Changed

- [`adni-to-bids`] Improve fetching of participants
- [`adni-to-bids`] Image path finder more robust
- [Doc] Update the OASIS3 documentation
- [CI] Code refactoring/cleanup

### Fixed

- [Atlas] Fix ROI index for left amygdala in AAL2 atlas
- [`adni-to-bids`] Prevent crash when files exists
- [`adni-to-bids`] Revert behavior to encode Dementia as AD
- [`adni-to-bids`] Remove entries with incoherent session names
- [`nifd-to-bids`] Several bugfixes and enhancements
- [I/O tools] Fix bug on empty dataframe
- [CI] Fix bash instruction to init conda
- [Doc] Correct DWI-Connectome description paragraph

## Clinica 0.5.0

### Added

- [Doc] Add missing documentation on check-missing-processing iotool
- [`deeplearning-prepare-data`]: Add option to run pipeline with ROI for tensor_format option to extract region of interest according to a mask.

### Changed

- [Core] Improve Logging for Clinica
- [Core] Improve CLI through using Click
- [Core] Nibabel replace get_data() by get_fdata() method for dataobj_images (nibabel)
- [`adni-to-bids`] Optimization of `adni-to-bids` clinical data extraction
- [`adni-to-bids`] Replace xlsx by tsv files for clinical data specification

### Fixed

- [Core] Fix bug in `write_scan_tsv` function
- [Doc] Add documentation for `check-missing-processing` command
- [Doc] Fix several small typos
- [Doc] Instructions for installing SPM dependency on MacOs Big Sur
- [CI] Fix several small issues with non-regression tests
- [CI] Fix typo in Jenkins script
- [CI] Automatically delete conda environments after PR is merged
- [ML] Fix unresolved reference in SVM pipeline
- [ML] Fix typo in parameters for SVC pipeline

## Clinica 0.4.1

### Added

- [`deeplearning-prepare-data`]: Add option to run pipeline with pet-linear outputs

### Changed

- [`oasis-to-bids`]: Remove FSL library dependency for OASIS-to-bids conversion.
- [Clinica]: Replace exception by warning when CAPs folder not recognized.
- [`aibl-to-bids`]: Center output nifti files of AIBL.
- [`aibl-to-bids`]: Extracts DICOM metadata in JSON files.
- [`merge-tsv`]: Fetch subcortical volumes generated by `t1-freesurfer` pipelines and JSON files (if exists).
- Fix minor typos in documentation.

### Fixed

- [`pet-surface`]: Verify SPM12 installation when running pipeline

## Clinica 0.4.0

### Added

- `pet-linear` pipeline: spatial normalization to the MNI space and intensity normalization of PET images
- `pet-surface-longitudinal` pipeline: Surface-based longitudinal processing of PET images
- `check-missing-processing` tool allows creating a TSV file containing information about the pipelines executed into a specific CAPS folder
- Conversion information is added once the converter is run to facilitate traceability.
- Add new keywords available in ADNI3 to the `adni-to-bids` converter

### Changed

- Harmonize output message display when running a pipeline.
- Automatically ignore an image to process if it is found in the CAPS folder.

### Fixed

- Fix a minor bug when the pet-volume pipeline performs PVC.
- Code source was completely reformatted using the Black code style.
- Documentation for the project is now versioned (versions from 0.3.8 are publicly available).
- Functions used for multiple pipelines are now mutualized (e.g container_from_filename function).
- f-strings are used massively.

## Clinica 0.3.8

### Added

- Add option to run `deeplearning-prepare-data` in output of `t1-extension` pipeline
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
  wiki: <https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer_Longitudinal/>
- [Change] The `fmri-preprocessing` pipeline is removed from the Clinica
  software as we will not actively maintain it. It is now in a separate
  repository: <https://github.com/aramis-lab/clinica_pipeline_fmri_preprocessing/>

#### Converters:

- [Enh] Improve how dependencies are checked for converters
- [Fix] Add diagnosis conversion for ADNI3
- [Fix] Avoid creation of sessions `ses-V01` in `*_sessions.tsv` files

## Clinica 0.3.6 - Pip

### Changes

#### Clinica Core:

- [Change] Pip is the main way to install Clinica. Conda packages are not
  available for the new versions.
- [Update] Set the minimal version of Python to 3.7.
- [Fix] Remove non-breaking spaces.

#### Pipelines:

- [New] Display failed image(s) when running `t1-linear` pipeline.

#### Converters:

- [Fix] The `aibl-to-bids` converter now handles new version of clinical data.
- [Update] The `oasis-to-bids` converter now uses NiBabel instead of FreeSurfer
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
  <https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DeepLearning_PrepareData/>
- [Change] `t1-linear` pipeline now crops image on default. If
  `--uncropped_image` is added to the command line, the image is not cropped.
- [Change] Refactor machine learning modules. Main changes involve use of
  CamelCase convention for classes and parameters used dictionaries.

## Clinica 0.3.4 - Bugfixes

### Changes

#### Clinica Core:

- [Improvement] Remove Clinica dependencies while updating and unfreezing some
  of them

#### Pipelines:

- [Enh] Improve how template files are downloaded for `t1-linear` and
  `statistics-volume` pipelines
