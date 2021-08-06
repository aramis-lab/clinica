<!-- markdownlint-disable MD046 -->
# `aibl-2-bids` – Conversion of the Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) to BIDS

!!! quote "Description reproduced from the [AIBL's Study Design webpage](http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/)"
    The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD.
    Although AIBL and ADNI have many of the same goals, there are differences between the two projects.
    Read more about the AIBL study from their [website](http://www.aibl.csiro.au/).

    **Study Participants** AIBL has enrolled 1100 participants and collected over 4.5 years worth of longitudinal data:

      - 211 AD patients
      - 133 MCI patients
      - 768 comparable healthy controls

    **Data** AIBL follows ADNI 1 protocols. Available data includes:

      - Clinical and cognitive data
      - Image data: MRI, PET
      - Biomarkers data: blood, genotype, ApoE
      - Dietary/lifestyle assessment

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of the **dcm2niix** DICOM to NIfTI converter.
You can find how to install these software packages on the [installation](../../#installing-clinica-from-source) page.

### Downloading AIBL

To download the AIBL dataset you first need to register to the [LONI Image & Data Archive (IDA)](https://ida.loni.usc.edu/login.jsp), a secure research data repository, and then request access to the AIBL dataset through the submission of an [online application form](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=AIBL).

In order to use the converter, you will need to download both the images and the clinical data.
To do so, from the [main page](https://ida.loni.usc.edu/login.jsp?returnPage=UserManagement.jsp&project=) click on `PROJECTS` and `AIBL`.
To download the imaging data, click on `Download` and choose `Image collections`.
In the `Advanced search` tab, pick the images you wish to download, for example tick `MRI` to download all the MR images, and then click on `SEARCH`.
In the `Advanced search results` tab, click `Select All` and `Add To Collection`.
Finally, in the `Data Collection` tab, select the collection you just created, tick `All` and click on `Advanced download`.
We advise you to group files as 10 zip files.
To download the clinical data, from the [main page of the AIBL project](https://ida.loni.usc.edu/home/projectPage.jsp?project=AIBL), click on `Download Clinical Data` and on the next page click on `DOWNLOAD`.

!!! note
    You do not have to modify the original folder name or rename the clinical data files before using the converter.

### Modalities supported

Currently, the modalities supported by our converter are:

- T1-weighted MRI
- Pittsburgh compound B (PiB) PET
- Florbetapir (AV45) PET
- Flutemetamol (FLUTE) PET
- Clinical data

The conversion of the imaging data to BIDS relies on modality-specific csv files that provide the list of scans available.
For each AIBL participant, the only T1w MR, PiB PET, Florbetapir PET, and Flutemetamol PET image available per session is converted.
For each image, the coordinates of the origin are set to the center of the box containing the image data.
This allows other image processing pipelines in Clinica (mainly SPM based) to run without needing further image preprocessing.

The conversion of the clinical data relies on the list of subjects and sessions obtained after the conversion of the imaging data and on the csv files containing the non-imaging data.
Data that do not change over time are gathered in the `participants.tsv` file, located at the top of the BIDS folder hierarchy, while the session-dependent data are gathered in `<subjectID>_session.tsv` files in each participant subfolder.
The clinical data being converted are defined in a spreadsheet (`clinical_specifications.xlsx`) available with the code of the converter, which the user can modify.

### Using the converter

The converter can be run with the following command line:

```shell
clinica convert aibl-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original AIBL images' directory;
- `CLINICAL_DATA_DIRECTORY` is the path to the directory where the csv file with the clinical data is located;
- `BIDS_DIRECTORY` is the path to the output directory, where the BIDS-converted version of AIBL will be stored.

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The AIBL data have been curated and converted to the Brain Imaging Data Structure (BIDS) format
    [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675);
    [Samper-González et al., 2018](https://doi.org/10.1016/j.neuroimage.2018.08.042)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/NASGJPVL).
