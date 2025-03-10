--8<-- [start:matlab]
## :warning: Known issue with Matlab and SPM12

[Matlab](https://www.mathworks.com/products/matlab.html) and [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (whose implementation is based on Matlab) can sometimes randomly crash, causing a rather unreadable error in the console.
Those events are unpredictable. In case it occurs to you, please do the following:

- Check that you have a valid Matlab license.
- Before relaunching the command line, be sure to remove the content of the working directory (if you specified one).
--8<-- [end:matlab]


--8<-- [start:center-nifti]
!! warning "Centering BIDS nifti"
    If the images from the `BIDS_DIRECTORY` are not centered, Clinica will give a **warning** because this can be an issue **if** later processing steps involve SPM (for instance if you are planning to run [pet-surface](./PET_Surface.md) afterwards).
    The warning message will contain a suggestion of a command to be run on your `BIDS_DIRECTORY` in order to generate a new BIDS dataset with images centered. This relies on the IOTool [center-nifti](../IOTools/center_nifti.md).
    It is highly recommended to follow this recommendation but Clinica won't force you to do so.
--8<-- [end:center-nifti]
