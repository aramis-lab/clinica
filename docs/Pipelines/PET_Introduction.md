<!-- markdownlint-disable MD007 -->
# Introduction

!!! note "Clinica & BIDS specifications for PET modality"
    Since Clinica `v0.6`, PET data following the official specifications in BIDS version 1.6.0 are now compatible with Clinica.
    See [BIDS](../../BIDS) page for more information.

## Partial volume correction (PVC)

To correct for partial volume effects, several [PVC](../glossary.md#pvc) algorithms exist and are implemented in the [PETPVC toolbox](https://github.com/UCL/PETPVC).

To perform PVC (compulsory for [`pet-surface`](../PET_Surface), optional for [`pet-volume`](../PET_Volume)), you will need to specify in a [TSV](../glossary.md#tsv) file the full width at half maximum (FWHM), in millimeters, of the [PSF](../glossary.md#psf) associated with your data, in the x, y and z directions.

For instance, if the FWHM of the PSF associated with your first image is 5 mm along the x and y axes, and 6 mm along the z axis, the first row of your [TSV](../glossary.md#tsv) file will look like this:

```Text
participant_id    session_id     acq_label     psf_x    psf_y    psf_z
sub-CLNC0001      ses-M000        18FFDG        5        5        6
sub-CLNC0001      ses-M000        18FAV45       4.5      4.5      5
sub-CLNC0002      ses-M000        18FFDG        5        5        6
sub-CLNC0003      ses-M000        18FFDG        7        7        7
```

Since the PSF depends on the [PET](../glossary.md#pet) tracer and scanner, the `participant_id`, `session_id`, `acq_label`, `psf_x`, `psf_y` and `psf_z` columns are compulsory.

The values in the column `acq_label` should match the value associated to the `trc` key in the [BIDS](../glossary.md#bids) dataset.

For example in the following [BIDS](../glossary.md#bids) layout the values associated would be `18FFDG` and `18FAV45`:

```text
bids
└─ sub-CLNC0001
   └─ ses-M000
      ├─ sub-CLNC001_ses-M000_trc-18FAV45_pet.nii.gz
      └─ sub-CLNC001_ses-M000_trc-18FFDG_pet.nii.gz
```

## Reference regions used for intensity normalization

In neurology, an approach widely used to allow inter- and intra-subject comparison of [PET](../glossary.md#pet) images is to compute standardized uptake value ratio ([SUVR](../glossary.md#suvr)) maps.
The images are intensity normalized by dividing each [voxel](../glossary.md#voxel) of the image by the average uptake in a reference region.
This region is chosen according to the tracer and disease studied as it must be unaffected by the disease.

Clinica `v0.3.8` introduces the possibility for the user to select the reference region for the [SUVR](../glossary.md#suvr) map computation.

Reference regions provided by Clinica come from the [Pick atlas](https://www.nitrc.org/projects/wfu_pickatlas) in MNI space and currently are:

- `pons`: 6 mm eroded version of the pons region
- `cerebellumPons`: 6 mm eroded version of the cerebellum + pons regions
- `pons2`: new in Clinica `v0.4`
- `cerebellumPons2`: new in Clinica `v0.4`

In Clinica `v0.4` two new versions of these masks have been introduced: `pons2` and `cerebellumPons2`.
Indeed, we wanted to improve the reference regions mask to better fit the [MNI152NLin2009cSym template](https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#template-based-coordinate-systems) used in linear processing pipelines.
The new masks still come from the [Pick atlas](https://www.nitrc.org/projects/wfu_pickatlas) but with a different processing: we decided to first truncate the mask using SPM12 tissue probability maps to remove voxels overlapping with regions outside the brain (bone, CSF, background...).
Then, we eroded the mask using scipy [`binary_erosion`](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.morphology.binary_erosion.html) with 3 iterations.

## Tutorial: How to add new SUVR reference regions to Clinica?

It is possible to run the [`pet-surface`](../PET_Surface) and [`pet-volume`](../PET_Volume) pipelines using a custom reference region.

- Install Clinica following the [developer instructions](../../Installation/#install-clinica) ;

- In the `<clinica>/clinica/utils/pet.py` file:
    - **Step 1:** Define the label of the [SUVR](../glossary.md#suvr) reference region that will be stored in CAPS filename(s).
      To do so, you need to add a variant to the `SUVRReferenceRegion` enumeration, which should look like this:

    ```python
    class SUVRReferenceRegion(str, Enum):
        PONS = "pons"
        CEREBELLUM_PONS = "cerebellumPons"
        PONS2 = "pons2"
        CEREBELLUM_PONS2 = "cerebellumPons2"
    ```

    Simply define a new label that will be your new [SUVR](../glossary.md#suvr) reference region.
    The `SUVRReferenceRegion` enumeration is used by all command-line interfaces, so you do not need to modify the pipelines' CLI to make this new region appear.

    - **Step 2:** Define the path of the [SUVR](../glossary.md#suvr) reference region that you will use.
      The function responsible to get the [SUVR](../glossary.md#suvr) mask is called `get_suvr_mask`, and it looks by default in the folder `<clinica>/resources/masks/`.
      You can put your mask in this folder and edit the following function (add an if statement to handle the enumeration variant you added in step 1): 

    ```python
    def _get_suvr_reference_region_labels_filename(region: SUVRReferenceRegion) -> str:
        if region == SUVRReferenceRegion.PONS:
            return "region-pons_eroded-6mm_mask.nii.gz"
        if region == SUVRReferenceRegion.CEREBELLUM_PONS:
            return "region-cerebellumPons_eroded-6mm_mask.nii.gz"
        if region == SUVRReferenceRegion.PONS2:
            return "region-pons_remove-extrabrain_eroded-2it_mask.nii.gz"
        if region == SUVRReferenceRegion.CEREBELLUM_PONS2:
            return "region-cerebellumPons_remove-extrabrain_eroded-3it_mask.nii.gz"
    ```
