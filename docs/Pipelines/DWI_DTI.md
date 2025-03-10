<!-- markdownlint-disable MD007 MD046 -->
# `dwi-dti` – DTI-based processing of corrected DWI datasets

The `dwi-dti` pipeline computes diffusion tensor imaging (DTI) with extraction of DTI-based measures, namely the fractional anisotropy (FA), mean diffusivity (MD), axial diffusivity (AD) and radial diffusivity (RD).
Then, the DTI-derived scalar maps (FA, MD, AD, RD) are normalized onto an FA-atlas with labelled tracts.
Finally, TSV files containing a summary of the regional statistics (mean DTI-based measures) are generated to ease subsequent statistical analyses.

To that aim, it mainly relies on the **MRtrix3** [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)] software for DTI computations and on **ANTs** for the normalization aspects [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)].

## Prerequisites

You need [preprocessed DWI data](../DWI_Preprocessing) prior to running any of these pipelines.

## Dependencies
If you only installed the core of Clinica, this pipeline needs the installation of [ANTs v2.5.0](../Software/Third-party.md#ants), [FSL 6.0](../Software/Third-party.md#fsl), and [MRtrix3](../Software/Third-party.md#mrtrix3) on your computer.

## Running the pipeline

The `dwi-dti` pipeline can be run with the following command line:

```Text
clinica run dwi-dti CAPS_DIRECTORY
```

where:

- `CAPS_DIRECTORY` is the input/output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.

??? info "Optional parameters common to all pipelines"
    --8<-- "snippets/pipelines_options.md:all"

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications/#dwi-dti-dti-based-processing-of-corrected-dwi-datasets):
`subjects/<participant_id>/<session_id>/dwi/dti_based_processing/`.

The main output files are:

- native_space/:
    - `<source_file>_space-{b0|T1w}_model-DTI_diffmodel.nii.gz`:
    The diffusion tensor imaging (DTI) data of the subject.
    - `<source_file>_space-[b0|T1w]_{FA|MD|AD|RD}.nii.gz`:
    The DTI-based measures, namely the fractional anisotropy (`FA`), mean diffusivity (`MD`), axial diffusivity (`AD`) and radial diffusivity (`RD`).
- normalized_space/
    - `<source_file>_space-<space>_{FA|MD|AD|RD}.nii.gz`:
    The DTI-based measures registered to the space of an FA-atlas.
- atlas_statistics/
    - `<source_file>_space-<space>_map-{FA|MD|AD|RD}_statistics.tsv`:
    TSV files summarizing the regional statistics on the labelled atlas `<space>`.

!!! note "Atlases available for the DTI-based processing pipeline:"
    - [JHUDTI81](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases)
    [[Hua et al., 2008](https://doi.org/10.1016/j.neuroimage.2007.07.053);
    [Wakana et al., 2007](https://doi.org/10.1016/j.neuroimage.2007.02.049)]:
    This atlas contains 48 white matter tract labels that were created by manually segmenting a standard-space average of diffusion MRI tensor maps from 81 subjects.
    - [JHUTracts 0|25|50](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases)
    [[Mori et al., 2005](https://www.elsevier.com/books/mri-atlas-of-human-white-matter/mori/978-0-444-51741-8)].
    This atlas contains 20 white matter tract labels that were identified probabilistically by averaging the results of deterministic tractography run on 28 subjects.
    Several thresholds of these probabilistic tracts are proposed (0%, 25%, 50%).

!!! note
    The full list of output files can be found in [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#dwi-dti-dti-based-processing-of-corrected-dwi-datasets).

## Describing this pipeline in your paper

??? cite "Example of paragraph for the `dwi-dti` pipeline:"
    These results have been obtained using the `dwi-dti` pipeline of Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].
    A diffusion tensor imaging (DTI) model was fitted to each voxel to calculate the fractional anisotropy (FA), mean diffusivity (MD), radial diffusivity (RD) and axial diffusivity (AD) maps using **MRtrix** [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)].
    The FA map of each subject was then registered to the FA map of the JHU atlas template with the ANTs SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)], and the estimated non-linear deformation was applied to the MD, AD and RD maps to have, for each individual, all the DTI-based maps in the space of the JHU atlas.

    We then assessed the integrity of a set of anatomical white matter tracts defined in the:

    - (Description for JHUDTI81 atlas) DTI-81 white-matter atlas [[Hua et al., 2008](https://doi.org/10.1016/j.neuroimage.2007.07.053); [Wakana et al., 2007](https://doi.org/10.1016/j.neuroimage.2007.02.049)].
    This atlas contains 48 white matter tract labels that were created by manually segmenting a standard-space average of diffusion MRI tensor maps from 81 subjects.

    - (Description for JHUTracts[0|25|50] atlas) JHU white-matter tractography atlas [Mori et al., 2005].
    This atlas contains 20 white matter tract labels that were identified probabilistically by averaging the results of deterministic tractography run on 28 subjects.
    Several thresholds of these probabilistic tracts are proposed (0%, 25%, 50%).

    The warping of this atlas to each individual subject provides a parcellation of the subject’s white matter into anatomical tracts.
    The integrity of the tracts was assessed by analyzing the average FA, MD, AD and RD in each tract.
    To that purpose, the scalar maps of each subject were put into correspondence with the FA-map in the atlas space using the ANTs SyN algorithm [[Avants et al., 2008](https://doi.org/10.1016/j.media.2007.06.004)].
    Finally, for each subject, the mean scalar value in each tract was computed for each DTI-based measure.

!!! tip
    Easily access the papers cited on this page on Zotero: [DTI](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/9URIGJNJ).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
