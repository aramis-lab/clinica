# `statistics_volume_correction`` - Performs corrections after analysis in statistics-volume pipeline



## Dependencies
No dependency is required.

## Running the pipeline
The pipeline can be run with the following command line:

```
clinica run statistics_volume_correction caps_directory t_map height_threshold FWEp FDRp FWEc FDRc 
```
where:

  - `caps_directory` is the output folder containing the results in a [CAPS](../CAPS) hierarchy.
  - `t_map`: name of the T statistic map used for the correction
  - `height_threshold`: height threshold indicated in the SPM report in output of `statistics-volume` pipeline
  - `FWEp`: threshold indicated in the SPM report in output of `statistics-volume` pipeline
  - `FDRp`: indicated in the SPM report in output of `statistics-volume` pipeline
  - `FWEc`: indicated in the SPM report in output of `statistics-volume` pipeline
  - `FDRc`: indicated in the SPM report in output of `statistics-volume` pipeline

Optional parameters:

  - `n_cuts`: number of cuts in final visualization


## Outputs

Results are stored in the following folder of the [CAPS hierarchy](docs/CAPS): `subjects/sub-<participant_label>/ses-<session_label>/<some_folder>`.

The main output files are:

  - `<source_file>_labelname-<label>_mainouput1`: description main output 1.

  - `<source_file>_labelname-<label>_mainouput2`: description main output 2.

The full list of output files can be found in the [ClinicA Processed Structure (CAPS) Specification](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8/edit#heading=h.f4ddnk971gkn).


<!--## Visualization of the results-->

<!--!!! note-->
<!--    The visualization command is not available for the moment. Please come back later, this section will be updated ASAP.-->


## Describing this pipeline in your paper

<!--You can have a single version for your pipeline-->

!!! cite "Example of paragraph (short version):"
    These results have been obtained using the `statistics_volume_correction` pipeline of Clinica. This pipeline is a ...

!!! cite "Example of paragraph (long version):"
    These results have been obtained using the `statistics_volume_correction` pipeline of Clinica. More precisely,...

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/1517933/aramis_clinica/items/collectionKey/2DHP3WXH).
