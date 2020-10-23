# Introduction


## Partial volume correction (PVC)

To correct for [partial volume effects](http://www.turkupetcentre.net/petanalysis/image_pve.html), several PVC algorithms exist and are implemented in the [PETPVC toolbox](https://github.com/UCL/PETPVC).

To perform PVC (compulsory for [`pet-surface`](../PET_Surface), optional for [`pet-volume`](../PET_Volume)), you will need to specify in a TSV file the full width at half maximum (FWHM), in millimeters, of the [point spread function (PSF)](https://en.wikipedia.org/wiki/Point_spread_function) associated with your data, in the x, y and z directions.

For instance, if the FWHM of the PSF associated with your first image is 5 mm along the x and y axes, and 6 mm along the z axis, the first row of your TSV file will look like this:

```
participant_id    session_id     acq_label     psf_x    psf_y    psf_z
sub-CLNC0001      ses-M00        FDG           5        5        6
sub-CLNC0001      ses-M00        AV45          4.5      4.5      5
sub-CLNC0002      ses-M00        FDG           5        5        6
sub-CLNC0003      ses-M00        FDG           7        7        7
```

Since the PSF depends on the PET tracer and scanner, the `participant_id`, `session_id`, ` acq_label`, ` psf_x`, `psf_y` and `psf_z` columns are compulsory.


## Reference regions used for intensity normalization

In neurology, an approach widely used to allow inter- and intra-subject comparison of PET images is to compute standardized uptake value ratio (SUVR) maps. The images are intensity normalized by dividing each voxel of the image by the average uptake in a reference region. This region is chosen according to the tracer and disease studied as it must be unaffected by the disease.

Clinica `v0.3.8` introduces the possibility for the user to select the reference region for the SUVR map computation.

Reference regions provided by Clinica come from the [Pick atlas](https://www.nitrc.org/projects/wfu_pickatlas) in MNI space and currently are:

- `pons`: 6 mm eroded version of the pons region

- `cerebellumPons`:  6 mm eroded version of the cerebellum + pons regions

## Tutorial: How to add new SUVR reference regions to Clinica?

It is possible to run the [`pet-surface`](../PET_Surface) and [`pet-volume`](../PET_Volume) pipelines using a custom reference region.

- Install Clinica following the [developer instructions](../../Installation/#install-clinica) ;

- In the `<clinica>/clinica/utils/pet.py` file, modify the following two elements:
    - The label of the SUVR reference region that will be stored in CAPS filename(s):
    ```python
    LIST_SUVR_REFERENCE_REGIONS = [
        "pons",
        "cerebellumPons",
    ]
    ```
    Simply define a new label that will be your new SUVR reference region. `LIST_SUVR_REFERENCE_REGIONS` is used by all command-line interfaces so you do not need to modify the pipelines' CLI to make this new region appear.

    - The path of the SUVR reference region that you will use:
    ```python
    def get_suvr_mask(suvr_reference_region):
        """Get path of the SUVR mask from SUVR reference region label.

        Args:
            suvr_reference_region: Label of the SUVR reference region

        Returns:
            Path of the SUVR mask
        """
        import os

        suvr_reference_region_to_suvr = {
            "pons": os.path.join(
                os.path.split(os.path.realpath(__file__))[0],
                "..",
                "resources",
                "masks",
                "region-pons_eroded-6mm_mask.nii.gz",
            ),
            "cerebellumPons": os.path.join(
                os.path.split(os.path.realpath(__file__))[0],
                "..",
                "resources",
                "masks",
                "region-cerebellumPons_eroded-6mm_mask.nii.gz",
            ),
        }
        return suvr_reference_region_to_suvr[suvr_reference_region]
    ```
    In this example, the SUVR reference region associated with the `cerebellumPons` label is located at `<clinica>/resources/masks/region-cerebellumPons_eroded-6mm_mask.nii.gz`.
