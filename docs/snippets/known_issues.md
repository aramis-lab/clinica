--8<-- [start:matlab]
## :warning: Known issue with Matlab and SPM12

[Matlab](https://www.mathworks.com/products/matlab.html) and [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (whose implementation is based on Matlab) can sometimes randomly crash, causing a rather unreadable error in the console.
Those events are unpredictable. In case it occurs to you, please do the following:

- Check that you have a valid Matlab license.
- Before relaunching the command line, be sure to remove the content of the working directory (if you specified one).
--8<-- [end:matlab]


--8<-- [start:center-nifti]
!!! warning "Centering BIDS nifti"
    If the images from the `BIDS_DIRECTORY` are not centered, Clinica will give a **warning** because this can be an issue **if** later processing steps involve SPM (for instance if you are planning to run [pet-surface](./PET_Surface.md) afterwards).
    The warning message will contain a suggestion of a command to be run on your `BIDS_DIRECTORY` in order to generate a new BIDS dataset with images centered. This relies on the IOTool [center-nifti](../IOTools/center_nifti.md).
    It is highly recommended to follow this recommendation but Clinica won't force you to do so.
--8<-- [end:center-nifti]


--8<-- [start:gtmseg]
!!! failure "Known error on macOS with gtmseg"
     If you are running `pet-surface` or `pet-surface-longitudinal` on macOS, we noticed that if the path to the CAPS is too long, the pipeline fails when the `gtmseg` command from FreeSurfer is executed.
    This generates crash files with `gtmseg` in the filename, for instance:

    ```console
    $ nipypecli crash crash-20210404-115414-sheldon.cooper-gtmseg-278e3a57-294f-4121-8a46-9975801f24aa.pklz
    [...]
    Abort
    ERROR: mri_gtmseg --s sub-ADNI011S4105_ses-M000 --usf 2 --o gtmseg.mgz --apas apas+head.mgz --no-subseg-wm --no-keep-cc --no-keep-hypo
    gtmseg exited with errors
    Standard error:
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M000/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53505/tmp.mgh' (type = MGH )                       [ ok ]
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M00/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53727/tmp.mgh' (type = MGH )                       [ ok ]
    Saving result to '<caps>/subjects/sub-ADNI011S4105/ses-M000/t1/freesurfer_cross_sectional/sub-ADNI011S4105_ses-M000/tmp/tmpdir.xcerebralseg.50819/tmpdir.fscalc.53946/tmp.mgh' (type = MGH )                       [ ok ]
    Return code: 1
    ```

    This is under investigation (see [Issue #119](https://github.com/aramis-lab/clinica/issues/119) for details) and will be solved as soon as possible.
--8<-- [end:gtmseg]
