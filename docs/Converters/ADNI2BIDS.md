# `adni-2-bids` – Conversion of the Alzheimer’s Disease Neuroimaging Initiative (ADNI) to BIDS

!!! quote "Description reproduced from the [ADNI's Study Design webpage](http://adni.loni.usc.edu/study-design/)"
    ADNI is a global research effort that actively supports the investigation and development of treatments that slow or stop the progression of AD. This multisite, longitudinal study assesses clinical, imaging, genetic and biospecimen biomarkers through the process of normal aging to early mild cognitive impairment (EMCI), to late mild cognitive impairment (LMCI), to dementia or AD. With established, standardized methods for imaging and biomarker collection and analysis, ADNI facilitates a way for scientists to conduct cohesive research and share compatible data with other researchers around the world.

    The ADNI study has three phases: ADNI1, ADNI GO and ADNI2. New participants were recruited across North America during each phase of the study and agreed to complete a variety of imaging and clinical assessments. Participants are followed and reassessed over time to track the pathology of the disease as it progresses. Results are then shared by ADNI through USC’s Laboratory of Neuro Imaging’s Image Data Archive (IDA).

    The table below summarizes the North American ADNI study target participant numbers, as well as the types of data taken at the different phases. In reality, the total number of study participants vary.

    <center>![](http://adni.loni.usc.edu/wp-content/uploads/2012/10/about-study-design_2_UPDATED_WIP2.png)</center>

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of the **dcm2nii** and **dcm2niix** DICOM to NIfTI converters, and of **FreeSurfer**. You can find how to install these software packages on the [installation](../#installing-clinica-from-source) page.

## Downloading ADNI

To download the ADNI dataset you first need to register to the [LONI Image & Data Archive (IDA)](https://ida.loni.usc.edu/login.jsp), a secure research data repository, and then request access to the ADNI dataset through the submission of an [online application form](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=ADNI).

In order to use the converter, you will need to download both the images and the clinical data. To do so, from the [main page](https://ida.loni.usc.edu/login.jsp?returnPage=UserManagement.jsp&project=) click on `PROJECTS` and `ADNI`. To download the imaging data, click on `Download` and choose `Image collections`. In the 'Advanced search' tab, untick ADNI 3 as we currently do not support it, pick the images you wish to download, for example tick `MRI` to download all the MR images, and then click on `SEARCH`. In the `Advanced search results` tab, click `Select All` and `Add To Collection`. Finally, in the `Data Collection` tab, select the collection you just created, tick `All` and click on `Advanced download`. We advise you to group files as 10 zip files. To download the clinical data, click on `Download` and choose `Study Data`. Select all the CSV files which are present in ALL by ticking `Select ALL tabular data` and click `Download`.

!!! note
    You do not have to modify the original folder name or rename the clinical data files before using the converter.


!!! warning
    We do not currently support the conversion of ADNI 3. This will be done in the future.

## Modalities supported
Currently, the modalities supported by our converter are:

  - T1-weighted MRI
  - Diffusion weighted imaging (DWI)
  - Fluorodeoxyglucose (FDG) PET
  - Florbetapir (AV45) PET
  - Clinical data

To convert the imaging data to BIDS, a list of subjects with their sessions is first obtained from the `ADNIMERGE` spreadsheet. This list is compared for each modality of interest to the list of scans available, as provided by modality-specific csv files (e.g. `MRILIST.csv`). If the modality was acquired for a specific pair of subject-session, and several scans and/or preprocessed images are available, only one is converted. Regarding the T1 scans, when several are available for a single session, the preferred scan (as identified in `MAYOADIRL_MRI_IMAGEQC_12_08_15.csv`) is chosen. If a preferred scan is not specified then the higher quality scan (as defined in `MRIQUALITY.csv`) is selected. If no quality control is found, then we choose the first scan. Gradwarp and B1-inhomogeneity corrected images are selected when available as these corrections can be performed in a clinical setting, otherwise the original image is selected. 1.5 T images are preferred for ADNI 1 since they are available for a larger number of patients. Regarding the DWI scans, we selected imaging following the 'Axial DTI' sequence. Regarding the FDG and AV45 PET scans, the images co-registered and averaged across time frames are selected. The scans failing quality control (if specified in `PETQC.csv`) are discarded. Data that do not change over time, such as the subject's sex, education level or diagnosis at baseline, are obtained from the ADNIMERGE spreadsheet and gathered in the participants.tsv file, located at the top of the BIDS folder hierarchy. The session-dependent data, such as the clinical scores, are obtained from specific csv files (e.g. `MMSE.csv`) and gathered in `<subjectID> _session.tsv` files in each participant subfolder. The clinical data being converted are defined in a spreadsheet (`clinical_specifications_adni.xlsx`) that is available with the code of the converter. The user can easily modify this file if he/she wants to convert additional clinical data.

## Using the converter

The converter can be run with the following command line:

```
clinica convert adni-to-bids dataset_directory clinical_data_directory bids_directory
```

where:

  - `dataset_directory` is the path to the original ADNI images's directory;
  - `clinical_data_directory` is the path to the directory where the csv file with the clinical data is located;
  - `bids_directory` is the path to the output directory, where the BIDS-converted version of ADNI will be stored.

### Optional parameters
The converter offers the possibility of converting only the clinical data (once that the images are in the BIDS format) using the optional parameter `-c`.

Due to the high computational time requested for converting all the modalities of the whole ADNI dataset, it is possible to convert a single modality at the time using the parameter `-m` with one of the following values:

*  `T1` for the T1-weighted MRI
*  `DWI` for Diffusion Weighted Images (DWI)
*  `PET_FDG` for the Fluorodeoxyglucose (FDG) PET
*  `PET_AV45` for the Florbetapir (AV45) PET

Is also possible to provide the path to a .txt file with the list of subjects to convert using the optional parameter `-s`.

For more information about the optional parameters, you can type:

```
clinica convert adni-to-bids -h
```

??? failure "Known errors"
    - T1
       ```
       Subject sub-ADNI031S0830 for session ses-M48 folder contains sub-ADNI031S0830_ses-M48_T1w_Eq_1.nii
       Subject sub-ADNI100S0995 for session ses-M18 folder contains sub-ADNI100S0995_ses-M18_T1w_Eq_1.nii
       Subject sub-ADNI031S0867 for session ses-M48 folder contains sub-ADNI031S0867_ses-M48_T1w_Eq_1.nii
       Subject sub-ADNI100S0892 for session ses-M18 folder contains sub-ADNI100S0892_ses-M18_T1w_Eq_1.nii

       Subject sub-ADNI029S0845 for session ses-M24 folder is empty
       Subject sub-ADNI094S1267 for session ses-M24 folder is empty
       Subject sub-ADNI029S0843 for session ses-M24 folder is empty
       Subject sub-ADNI027S0307 for session ses-M48 folder is empty
       Subject sub-ADNI057S1269 for session ses-M24 folder is empty
       Subject sub-ADNI036S4899 for session ses-M03 folder is empty
       ```

    - DWI
       ```
       Subject sub-ADNI016S4638 for session ses-M00 folder contains sub-ADNI016S4638_ses-M00_dwi_Eq_1.nii
       Subject sub-ADNI007S4611 for session ses-M03 folder contains sub-ADNI007S4611_ses-M03_dwi_Eq_1.nii
       Subject sub-ADNI027S5118 for session ses-M00 folder contains sub-ADNI027S5118_ses-M00_dwi_Eq_1.nii
       Subject sub-ADNI094S2238 for session ses-M48 folder contains sub-ADNI094S2238_ses-M48_dwi_Eq_1.nii
       Subject sub-ADNI129S4287 for session ses-M00 folder contains sub-ADNI129S4287_ses-M00_dwi_Eq_1.nii

       Subject sub-ADNI098S4018 for session ses-M00 folder contains sub-ADNI098S4018_ses-M00_dwia.nii
       Subject sub-ADNI098S4003 for session ses-M12 folder contains sub-ADNI098S4003_ses-M12_dwia.nii

       Subject sub-ADNI029S4585 for session ses-M48 has wrong b-val/b-vec file
       Subject sub-ADNI029S2395 for session ses-M60 has wrong b-val/b-vec file
       Subject sub-ADNI029S0824 for session ses-M108 has wrong b-val/b-vec file
       Subject sub-ADNI029S0914 for session ses-M108 has wrong b-val/b-vec file
       Subject sub-ADNI029S4384 for session ses-M48 has wrong b-val/b-vec file
       Subject sub-ADNI029S4385 for session ses-M48 has wrong b-val/b-vec file
       Subject sub-ADNI094S4630 for session ses-M06 has wrong b-val/b-vec file
       Subject sub-ADNI094S4649 for session ses-M06 has wrong b-val/b-vec file
       Subject sub-ADNI029S5219 for session ses-M24 has wrong b-val/b-vec file

       Subject sub-ADNI027S2219 for session ses-M36 has wrong dimensions (256 x 256 x 2013)
       Subject sub-ADNI129S2332 for session ses-M12 has wrong dimensions (256 x 256 x 1549)
       ```

    - FDG PET
       ```
       Subject sub-ADNI941S1195 for session ses-M48 folder is empty
       Subject sub-ADNI005S0223 for session ses-M12 folder is empty
       Subject sub-ADNI037S1421 for session ses-M36 folder contains NONAME.nii
       Subject sub-ADNI037S1078 for session ses-M36 folder contains NONAME.nii
       ```

    - AV45 PET
      ```
      Subject sub-ADNI128S2220 for session ses-M48 folder contains sub-ADNI128S2220_ses-M48_task-rest_acq-AV45_pet_Eq_1.nii
      ```

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    - **T1 and/or PET**:
    The ADNI data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https:// doi.org/10.1038/sdata.2016.44)] using Clinica [[Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].
    - **DWI**:
    The ADNI data have been curated and converted to the Brain Imaging Data Structure (BIDS) format [[Gorgolewski et al., 2016](https:// doi.org/10.1038/sdata.2016.44)] using Clinica [[Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042); [Wen et al., 2018](https://arxiv.org/abs/1812.11183)].    

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
