# `statistics-volume` - Statistics for volume

Perform 2-sample t-test on your images

## Dependencies
If you installed the docker image of Clinica, nothing is required. If you only installed the core of Clinica, this pipeline needs the installation of **SPM** on its standalone form or with **Matlab** on your computer. You can find how to install this software on the [installation](../#installing-clinica-from-source) page.


## Running the pipeline
The pipeline can be run with the following command line:

```
clinica run statistics-volume caps_directory tsv_file contrast
```
where:

  - `caps_directory` is the output folder containing the results in a [CAPS](../CAPS) hierarchy.
  - `tsv_file` is the path to the TSV file with covariables (mandatory, unlike other Clinica pipelines)
  - `contrast` is the name of the contrast used to differntiate the 2 groups


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
    These results have been obtained using the `statistics_volume` pipeline of Clinica. This pipeline is a ...

!!! cite "Example of paragraph (long version):"
    These results have been obtained using the `statistics_volume` pipeline of Clinica. More precisely,...

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/1517933/aramis_clinica/items/collectionKey/2DHP3WXH).
