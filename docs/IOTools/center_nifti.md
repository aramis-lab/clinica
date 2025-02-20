# `center-nifti` - Center NIfTI files of a BIDS directory

Your [BIDS](http://bids.neuroimaging.io) dataset may contain NIfTI files where the origin does not correspond to the center of the image (i.e. the anterior commissure).
[SPM](../Software/Third-party.md#spm12) is especially sensitive to this case and segmentation procedures may result in blank images or even fail.
To mitigate this issue we offer a simple tool that generates from your BIDS a new dataset with centered NIfTI files for the selected modalities.


```shell
clinica iotools center-nifti [OPTIONS] BIDS_DIRECTORY OUTPUT_BIDS_DIRECTORY
```

where:

- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](http://bids.neuroimaging.io) hierarchy.
- `OUTPUT_BIDS_DIRECTORY` is the output path to the new version of your BIDS dataset, with faulty NIfTI centered.
This folder can be empty or nonexistent.

Optional arguments:

- `-m/--modality` is a case-insensitive parameter that defines which modalities are converted. By default, the tool centers only **T1w** images.

    !!! tip "How to use :"
        - If you want to convert T1w images only, do not use the option.
        - If you want to convert all types of pet, use `--modality pet` or `-m pet`
        - If you want to convert 18FFDG_PET, use `-m 18ffdg_pet`
        - If you want to convert both pet and T1, use `-m T1 -m pet`

        Basically, the software searches for the modality key inside the filename. Understanding this, you can now center any modality you want!

- `-t/--threshold` allows choosing the critical distance (mm) from the origin of the world coordinate system to an image center above which images will be centered. By default, this threshold is set to **50**. To center images regardless of their original distance to the center, set it to 0.

    !!! info "Default centering threshold"
        This threshold of 50 mm was tailored for SPM. It was chosen empirically after a set of experiments to determine at which distance from the origin SPM segmentation and coregistration procedures stop working properly.
    

!!! note
    - The images contained in the input `bids_directory` folder that do not need to be centered will also be copied to the output folder `output_bids_directory`. 
    - In this output directory, a text file named `centered_nifti_list_TIMESTAMP.txt` with the list of the converted files will also be created.
