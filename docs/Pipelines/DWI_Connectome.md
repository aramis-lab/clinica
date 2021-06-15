# `dwi-connectome` - Computation of structural connectome from corrected DWI datasets

The `dwi-connectome` pipeline computes a weighted graph encoding anatomical connections between a set of brain regions from corrected DWI datasets.

To aim that, it relies on the **MRtrix3** [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)] software to compute the constrained spherical deconvolution diffusion model, perform probabilistic tractography and computes a connectome using the Desikan & Destrieux atlases from **FreeSurfer**.

## Prerequisites

You need to [preprocess your DWI data](../DWI_Preprocessing) and run the [`t1-freesurfer`](../T1_FreeSurfer) pipeline on your T1-weighted MRI images prior to running this pipeline.

## Dependencies
<!-- If you installed the docker image of Clinica, nothing is required.-->

If you only installed the core of Clinica, this pipeline needs the installation of **FSL 6.0** and **MRtrix3** on your computer.

You can find how to install these software packages on the [third-party](../../Third-party) page.

## Running the pipeline

The pipeline can be run with the following command line:

```Text
clinica run dwi-connectome <caps_directory>
```

where:

- `caps_directory` is the input/output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.

If you want to run the pipeline on a subset of your CAPS dataset, you can use the `-tsv` flag to specify in a TSV file the participants belonging to your subset.

!!! note "Number of streamlines (`--n_tracks` parameters)"
    The quality of the tractography and, as a result, the connectome mainly depends on the number of streamlines you can generate (the more the better).
    However, increasing the number of streamlines increases the need for computational resources and space to store the results.
    On default, 1 million streamlines are generated which represents 1 hour of computation time.

## Outputs

Results are stored in the following folder of the
[CAPS hierarchy](../../CAPS/Specifications):
`subjects/<participant_id>/<session_id>/dwi/connectome_based_processing/`.

The main output files are:

- `<source_file>_space-{b0|T1w}_model-CSD_diffmodel.nii.gz`:
Constrained spherical deconvolution (CSD) diffusion model.
- `<source_file>_space-{b0|T1w}_model-CSD_tractography.tck`:
The whole-brain tractography.
- `<source_file>_space-{b0|T1w}_model-CSD_parcellation-{desikan|destrieux}_connectivity.tsv`:
The connectivity matrix based on the Desikan or Destrieux parcellation.

!!! note "Atlases available for the Connectome-based processing pipeline:"
    - [Desikan](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
    [[Desikan et al., 2006](https://doi.org/10.1016/j.neuroimage.2006.01.021)]:
    This atlas is a subdivision of the cerebral cortex into gyri and contains 34 regions per hemisphere.
    It was built using a dataset of 40 MRI scans from which 34 cortical ROIs were manually identified in each of the individual hemispheres.
    - [Destrieux](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
    [[Destrieux et al., 2010](https://dx.doi.org/10.1016%2Fj.neuroimage.2010.06.010)]:
    This atlas is a subdivision of the cerebral cortex into gyri and sulci, and contains 74 regions per hemisphere.
    It was built on anatomical MRI of 24 healthy subjects from which 74 cortical ROIs were manually identified in each of the individual hemispheres.

!!! note
    The full list of output files can be found in the [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications).

<!--## Visualization of the results-->

<!--
We advise you to use the following commands to visualize the tractography results of a given subject using `mrview` utility:

```shell
caps_directory= # Example: "MY_DATASET_CAPS"
participant_id= # Example: "sub-CLNC01"
session_id= # Example: "ses-M00"
atlas_label= # Example: "desikan"

mrview -mode 2 \
        -load                   ${caps_directory}/subjects/${participant_id}/${session_id}/dwi/preprocessing/${participant_id}_${session_id}_preproc.nii.gz \
        -tractography.load      ${caps_directory}/subjects/${participant_id}/${session_id}/dwi/connectome_based_processing/*_tractography.tck \
        -odf.load_sh            ${caps_directory}/subjects/${participant_id}/${session_id}/dwi/connectome_based_processing/*_FOD.mif \
        -connectome.init        ${caps_directory}/subjects/${participant_id}/${session_id}/dwi/connectome_based_processing/*_parcellation-${atlas_label}_node.nii.gz \
        -connectome.load        ${caps_directory}/subjects/${participant_id}/${session_id}/dwi/connectome_based_processing/*_parcellation-${atlas_label}_connectivity.tsv
```

Do not forget to fill in the missing information (after the `=` signs) and do not hesitate to remove lines of the `mrview` command that you may not be interested in or that may take to much time to load.
-->

## Describing this pipeline in your paper

!!! cite "Example of paragraph for the `dwi-connectome` pipeline"
    These results have been obtained using the `dwi-connectome` pipeline of Clinica
    [[Routier et al](https://hal.inria.fr/hal-02308126/)] relying on the **MRtrix3**
    [[Tournier et al., 2019](https://doi.org/10.1016/j.neuroimage.2019.116137)] software package.
    Fiber orientation distributions (FOD) at highly anisotropic voxels (FA >0.7) was computed to determine the response function, which was used for constrained spherical deconvolution to accurately estimate the FOD [[Tournier et al., 2007](https://doi.org/10.1016/j.neuroimage.2007.02.016)].
    Then, `<n_tracks>` fibers with a probabilistic tracking algorithm [[Tournier et al., 2010](https://cds.ismrm.org/protected/10MProceedings/files/1670_4298.pdf)] were generated.
    Default parameters included minimum length 20 mm, a step size of 0.2 mm, minimum radius of curvature of 1 mm and FOD cutoff of 0.1.
    All voxels in the 1-mm dilated white-matter mask were used as seeds and the tracking procedure was stopped if a fiber reached a voxel outside the mask or if a stopping criterion was met (high fiber curvature or low FOD).
    Finally, the connectome is estimated by counting the number of tracks connecting each pair of nodes according to the [Desikan|Destrieux] parcellation.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/UJRXE4AP).

## Support

- You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
- Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).
