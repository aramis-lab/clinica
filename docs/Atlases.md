<!-- markdownlint-disable MD007 -->
# Atlases used in Clinica

This page describes the different atlases used in the pipelines of Clinica.

## Volume atlases

### Processing of T1-weighted MRI & PET images

These atlases, all defined in MNI space, are mainly used when performing volumetric processing of T1 and PET images, as done in the  [`t1-volume-*`](../Pipelines/T1_Volume) and  [`pet-volume`](../Pipelines/PET_Volume) pipelines.

- [AAL2](http://www.gin.cnrs.fr/en/tools/aal-aal2/)
[Tzourio-Mazoyer et al., [2002](http://dx.doi.org/10.1006/nimg.2001.0978),
[2015](http://dx.doi.org/10.1016/j.neuroimage.2015.07.075)]
is anatomical atlas based on a single subject.
It is the updated version of AAL, which is probably the most widely used cortical parcellation map in the neuroimaging literature.
It was built using manual tracing on the spatially normalized single-subject high-resolution T1 volume in MNI space.
It is composed of 120 regions covering the whole cortex as well as the main subcortical structures.
- [AICHA](http://www.gin.cnrs.fr/en/tools/aicha/)
[[Joliot et al., 2015](http://dx.doi.org/10.1016/j.jneumeth.2015.07.013)]
is a functional atlas based on multiple subjects.
It was built using parcellation of group-level functional connectivity profiles computed from resting-state fMRI data of 281 healthy subjects.
It is composed of 384 regions covering the whole cortex as well as the main subcortical structures.
- [Hammers](http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD)
[[Hammers et al., 2003](http://dx.doi.org/10.1002/hbm.10123);
[Gousias et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.11.034)]
is an anatomical atlas based on multiple subjects.
It was built using manual tracing on anatomical MRI from 30 healthy subjects.
The individual subjects parcellations were then registered to MNI space to generate a probabilistic atlas as well as a maximum probability map.
The latter was used in the present work.
It is composed of 69 regions covering the whole cortex as well as he main subcortical structures.
- [LPBA40](http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD)
[[Shattuck et al., 2008](http://dx.doi.org/10.1016/j.neuroimage.2007.09.031)]
is an anatomical atlas based on multiple subjects.
It was built using manual tracing on anatomical MRI from 40 healthy subjects.
The individual subjects parcellations were then registered to MNI space to generate a maximum probability map.
It is composed of 56 regions covering the whole cortex as well as he main subcortical structures.
- [Neuromorphometrics](http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD)
is an anatomical atlas based on multiple subjects.
It was built using manual tracing on anatomical MRI from 30 healthy subjects.
The individual subjects parcellations were then registered to MNI space to generate a maximum probability map.
It is composed of 140 regions covering the whole cortex as well as he main subcortical structures.
Data were made available for the “[MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling](http://masiweb.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details)”.

The main difference between LBPA40, Hammers and Neuromorphometrics atlases is the degree of detail (i.e. the number of regions) of the anatomical parcellation.

### Processing of DWI images

These atlases, all defined in MNI space, are mainly used when processing DWI images, as done in the [`dwi-dti`](../Pipelines/DWI_DTI) pipeline.

- [JHUDTI81](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases)
[[Hua et al., 2008](https://doi.org/10.1016/j.neuroimage.2007.07.053);
[Wakana et al., 2007](https://doi.org/10.1016/j.neuroimage.2007.02.049)]:
This atlas contains 48 white matter tract labels that were created by manually segmenting a standard-space average of diffusion MRI tensor maps from 81 subjects.
- [JHUTracts[0|25|50]](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases)
[[Mori et al., 2005](https://www.elsevier.com/books/mri-atlas-of-human-white-matter/mori/978-0-444-51741-8)]:
This atlas contains 20 white matter tract labels that were identified probabilistically by averaging the results of deterministic tractography run on 28 subjects.
Several thresholds of these probabilistic tracts are proposed (0%, 25%, 50%).

## Surface atlases

These atlases are mainly used when processing T1-weighted images with the [`t1-freesurfer`](../Pipelines/T1_FreeSurfer) pipeline and PET images with the [`pet-surface`](../Pipelines/PET_Surface) pipeline.

- [Desikan](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
[[Desikan et al., 2006]](https://doi.org/10.1016/j.neuroimage.2006.01.021):
This atlas is a subdivision of the cerebral cortex into gyri and contains 34 regions per hemisphere.
It was built using a dataset of 40 MRI scans from which 34 cortical ROIs were manually identified in each of the individual hemispheres.
- [Destrieux](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation)
[[Destrieux et al., 2010]](https://dx.doi.org/10.1016%2Fj.neuroimage.2010.06.010):
This atlas is a subdivision of the cerebral cortex into gyri and sulci, and contains 74 regions per hemisphere.
It was built on anatomical MRI of 24 healthy subjects from which 74 cortical ROIs were manually identified in each of the individual hemispheres.

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/JPGDLCMZ).

## Tutorial: How to add a new volume atlas to Clinica?

It is possible to run the [`t1-volume`](../T1_Volume) and [`pet-volume`](../PET_Volume) <!--and [`dwi-dti`](../DWI_DTI)--> pipelines using a custom parcellation.
To do so:

- Install Clinica following the [developer instructions](../../Installation/#install-clinica);

- In the `<clinica>/clinica/utils/atlas.py` file, modify the following two elements:
    - The label of the volume atlas that will be stored in CAPS filename(s):

    ```python
    T1_VOLUME_ATLASES = [
        "AAL2",
        "AICHA",
        "Hammers",
        "LPBA40",
        "Neuromorphometrics",
    ]
    ```

    Simply define a new label that will be your new volume.
    `T1_VOLUME_ATLASES` is used by all the command-line interfaces using atlases from the [`t1-volume`](../T1_Volume) pipeline so you do not need to modify the pipelines' CLI to make this new region appear.
    The same rationale applies to `PET_VOLUME_ATLASES`<!-- and `DWI_DTI_VOLUME_ATLASES`-->.

    - Create a new class inherited from `AtlasAbstract` and fill the three compulsory methods.
    If we take for instance the AAL2 parcellation:

    ```python
    class AAL2(AtlasAbstract):
    def __init__(self):
        AtlasAbstract.__init__(self)

    @staticmethod
    def get_name_atlas():
        return "AAL2"

    @staticmethod
    def get_atlas_labels():
        from os.path import join, split, realpath

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AAL2_dseg.nii.gz",
        )

    @staticmethod
    def get_tsv_roi():
        from os.path import join, split, realpath

        return join(
            split(realpath(__file__))[0],
            "..",
            "resources",
            "atlases",
            "atlas-AAL2_dseg.tsv",
        )
    ```

    The string returned by the `get_name_atlas()` method must match the label given in the `{T1|PET}_VOLUME_ATLASES` list.
    The `get_atlas_labels()` method must return the path to the parcellation in NIfTI format while the `get_tsv_roi()` method must return the path to a TSV file.
    In this example, the labels and TSV files associated with the `AAL2` atlas are located at `<clinica>/resources/atlases/atlas-AAL2_dseg.{nii.gz|tsv}`.
    Finally, the TSV file must contain the `roi_value` and `roi_name` columns and looks like:

    ```text
    roi_value   roi_name
    0           Background
    2001        Precentral_L
    [...]       [...]
    9170        Vermis_10
    ```
