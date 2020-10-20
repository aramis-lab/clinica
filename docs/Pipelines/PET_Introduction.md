# Introduction


## Partial volume correction (PVC)

To correct for [partial volume effects](http://www.turkupetcentre.net/petanalysis/image_pve.html), several PVC algorithms exist and are implemented in the [PETPVC toolbox](https://github.com/UCL/PETPVC).

To perform PVC (compulsory for [`pet-surface` pipeline](../PET_Surface), optional for [`pet-volume` pipeline](../PET_Volume)), you will need to specify in a TSV file the full width at half maximum (FWHM), in millimeters, of the [point spread function (PSF)](https://en.wikipedia.org/wiki/Point_spread_function) associated with your data, in the x, y and z directions.

For instance, if the FWHM of the PSF associated with your first image is 8 mm along the x axis, 9 mm along the y axis, and 10 mm along z axis, the first row of your TSV file will look like this:

```text
participant_id    session_id     acq_label     psf_x    psf_y    psf_z
sub-CLNC01        ses-M00        FDG           8        9        10
sub-CLNC01        ses-M18        FDG           8        9        10
sub-CLNC01        ses-M00        AV45          7        6        5
sub-CLNC02        ses-M00        FDG           8        9        10
sub-CLNC03        ses-M00        FDG           8        9        10
```

Since PSF information may differ according to the PET tracer, `participant_id`, `session_id`, ` acq_label`, ` psf_x`, `psf_y` and `psf_z` columns are compulsory columns.



## Reference regions for standardized uptake value ratio (SUVR) map

Clinica `v0.3.8` introduces the possibility for the user to select the reference region for the SUVR map computation.

Reference regions provided by Clinica come from the Pick atlas in MNI space and currently are:

- `pons`: 6 mm eroded version of the pons region

- `cerebellumPons`:  6 mm eroded version of the cerebellum + pons regions



## Tutorial: How to add new SUVR reference region to Clinica?

If you need to use a reference region not provided by Clinica but still want to use [`pet-surface`](../PET_Surface) or [`pet-volume`](../PET_Volume) pipelines, it is possible to easily extend the list of SUVR regions.

- You first need to install Clinica following [developer instructions](../../Installation/#install-clinica) ;

- Once done you will need to modify your `<clinica>/clinica/utils/pet.py` file in particular the following two elements:
    - The label of the SUVR reference region that will be stored in CAPS filename(s):
    ```python
    LIST_SUVR_REFERENCE_REGIONS = [
        "pons",
        "cerebellumPons",
    ]
    ```
    Simply define a new label that will be your new SUVR reference region. `LIST_SUVR_REFERENCE_REGIONS` is used by all command-line interfaces so you do need to modify the pipelines' CLI to make appear this new region.

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
    In this example, the SUVR reference region associated to `cerebellumPons` label is located at `<clinica>/resources/masks/region-cerebellumPons_eroded-6mm_mask.nii.gz`.
