# `center-nifti` - Center NIfTI files of a BIDS directory

Your [BIDS](http://bids.neuroimaging.io) dataset may contain NIfTI files where the origin does not correspond to the center of the image (i.e. the anterior commissure).
[SPM](../Software/Third-party.md#spm12) is especially sensitive to this case and segmentation procedures may result in blank images or even fail.
To mitigate this issue we offer a simple tool that generates from your BIDS a new dataset with centered NIfTI files for the selected modalities.

!!! warning "By default :"

    - This tool will only center **T1w** images.
    - Only NIfTI volumes whose center is at **more than 50 mm** from the origin of the world coordinate system are centered. This threshold has been chosen empirically after a set of experiments to determine at which distance from the origin SPM segmentation and coregistration procedures stop working properly.

```shell
clinica iotools center-nifti [OPTIONS] BIDS_DIRECTORY OUTPUT_BIDS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](http://bids.neuroimaging.io) hierarchy.
- `OUTPUT_BIDS_DIRECTORY` is the output path to the new version of your BIDS dataset, with faulty NIfTI centered.
This folder can be empty or nonexistent.

Optional arguments:

- `--modality` is a case-insensitive parameter that defines which modalities are converted.

    !!! tip "How to use :"
        - If you want to convert T1w images only, do not use the option.
        - If you want to convert all types of pet, use `--modality pet` or `-m pet`
        - If you want to convert 18FFDG_PET, use `-m 18ffdg_pet`
        - If you want to convert both pet and T1, use `-m T1 -m pet`

        Basically, the software searches for the modality key inside the filename. Understanding this, you can now center any modality you want!

- `--center_all_files` is an option that forces Clinica to center all the files of the modalities selected with the `--modality` flag.

!!! note
    The images contained in the input `bids_directory` folder that do not need to be centered will also be copied to the output folder `new_bids_directory`.

     The list of the converted files will appear in a text file in `new_bids_directory/centered_nifti_list_TIMESTAMP.txt`.
