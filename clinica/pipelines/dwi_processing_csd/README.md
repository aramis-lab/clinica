# `tractography` - Tractography

<SHORT_DESCRIPTION>

## Prerequisites

This pipeline requires that the following pipelines have run on your BIDS dataset:

  - [`t1-freesurfer-cross-sectional`](T1_FreeSurfer)  in order to get the brain mask and atlas-based parcellations.

## Dependencies
<!-- If you installed the docker image of Clinica, nothing is required. -->
If you only installed the core of Clinica, this pipeline needs the installation of **FSL** and **MRTrix3** on your computer. You can find how to install this software on the [installation](../#installing-clinica-from-source) page.

## Running the pipeline
The pipeline can be run with the following command line:

```shell
clinica run tractography bids_directory caps_directory
```
where:

  - `bids_directory` is the input folder containing the dataset in a [BIDS](../BIDS) hierarchy.
  - `caps_directory` is the output folder containing the results in a [CAPS](../CAPS) hierarchy.

If you want to run the pipeline on a subset of your BIDS dataset, you can use the `-tsv` flag to specify in a TSV file the participants belonging to your subset.

### Options

<!-- I put all the options as header titles in order to be able to get the hyperlink ancher. -->

#### `--n_tracks`

MRtrix3 `tckgen`'s `-select` option: 

> set the desired number of streamlines to be selected by `tckgen`, after all selection criteria have been applied (i.e. inclusion/exclusion ROIs, min/max length, etc). tckgen will keep seeding streamlines until this number of streamlines have been selected, or the maximum allowed number of seeds has been exceeded (see -seeds option). By default, 5000 streamlines are to be selected.

### Outputs

Results are stored in the following folder of the [CAPS hierarchy](docs/CAPS): `subjects/sub-<participant_label>/ses-<session_label>/dwi/tractography`.

The main output files are:

  - `<source_file>_labelname-<label>_mainouput1`: description main output 1.
  - `<source_file>_labelname-<label>_mainouput2`: description main output 2.

The full list of output files can be found in the [ClinicA Processed Structure (CAPS) Specification](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8/edit#heading=h.f4ddnk971gkn).


### Visualization of the results

We advise you to use the following commands to visualize the tractography results of a given subject using `mrview` utility:

```shell
caps_directory= # Example: "MY_DATASET_CAPS"
subject_id= # Example: "01"
session_id= # Example: "M00"
atlas_id= # Example: "desikan"

mrview -mode 2 \
	-load 			${caps_directory}/subjects/sub-${subject_id}/ses-${session_id}/t1/freesurfer_cross_sectional/sub-${subject_id}_ses-${session_id}/mri/orig.mgz \
	-tractography.load 	${caps_directory}/subjects/sub-${subject_id}/ses-${session_id}/dwi/tractography/*.tck \
	-odf.load_sh		${caps_directory}/subjects/sub-${subject_id}/ses-${session_id}/dwi/tractography/*_fod.mif \
	-connectome.init 	${caps_directory}/subjects/sub-${subject_id}/ses-${session_id}/dwi/tractography/*_parcellation-${atlas_id}_nodes.mif \
	-connectome.load 	${caps_directory}/subjects/sub-${subject_id}/ses-${session_id}/dwi/tractography/*_parcellation-${atlas_id}_connectome.csv
```

Do not forget to fill in the missing information (after the `=` signs) and do not hesitate to remove lines of the `mrview` command that you may not be interested in or that may take to much time to load.

## Describing this pipeline in your paper

<!--You can have a single version for your pipeline-->

!!! cite "Example of paragraph (short version):"
    These results have been obtained using the `tractography` pipeline of Clinica. This pipeline is a ...

!!! cite "Example of paragraph (long version):"
    These results have been obtained using the `tractography` pipeline of Clinica. More precisely,...

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/1517933/aramis_clinica/items/collectionKey/2DHP3WXH).
