# `t1-spm-dartel-existing-template`` - <VERY_SHORT_DESCRIPTION>

<SHORT_DESCRIPTION>


#### Contents

- [Dependencies](#dependencies)
- [Running the pipeline (command line)](#running-the-pipeline-command-line)
- [Outputs](#outputs)
- [Visualization of the results](#visualization-of-the-results)
- [Describing this pipeline in your paper](#describing-this-pipeline-in-your-paper)
- [Appendix](#appendix)


## Dependencies

If you installed the docker image of Clinica, nothing is required.

If you only installed the core of Clinica, this pipeline needs the installation of **<software_package>** on your computer. You can find how to install this software on the [installation](docs/BeforeYouInstall) page.

## Running the pipeline (command line)
The pipeline can be run with the following command line:
```
clinica run t1-spm-dartel-existing-template bids_directory caps_directory
```
where:
- `bids_directory` is the input folder containing the dataset in a [BIDS](docs/BIDS) hierarchy
- `caps_directory` is the output folder containing the results in a [CAPS](docs/CAPS) hierarchy
- `<ARG_1>` <ARG_1_DESCRIPTION>
- `<ARG_2>` <ARG_2_DESCRIPTION>

## Outputs

Results are stored in the following folder of the [CAPS hierarchy](docs/CAPS): `subjects/sub-<participant_label>/ses-<session_label>/<some_folder>`.

The main output files are:
- `<source_file>_main_ouput_1`: description main output 1.
- `<source_file>_main_ouput_2`: description main output 2.

The full list of output files can be found in the [ClinicA Processed Structure (CAPS) Specification](https://docs.google.com/document/d/14mjXbqRceHK0fD0BIONniLK713zY7DbQHJEV7kxqsd8/edit#heading=h.f4ddnk971gkn).


## Visualization of the results
After the execution of the pipeline, you can check the outputs of a subject by running the command:
> **Notes:**
>
> _The visualization command is not available for the moment. Please come back later, this section will be updated ASAP._


## Describing this pipeline in your paper

> **Example of paragraph:**
>
>_These results have been obtained using the my-pipeline pipeline of Clinica. More precisely ..._

## Appendix
Further information can be found on [this supplementary page](docs/Pipelines/<My_Pipeline_Appendix>).