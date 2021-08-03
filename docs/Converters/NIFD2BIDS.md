# `nifd-2-bids` – Conversion of Neuroimaging in Frontotemporal Dementia (NIFD) to BIDS

!!! quote "Description reproduced from the [NIFD's LONI Image & Data Archive (IDA) webpage](https://ida.loni.usc.edu/home/projectPage.jsp?project=NIFD&page=HOME&subPage=OVERVIEW_PR#)"
    NIFD is the nickname for the frontotemporal lobar degeneration neuroimaging initiative (FTLDNI, AG032306), which was funded by the NIA and NINDS to characterize longitudinal clinical and imaging changes in FTLD.
    The imaging and clinical methods are the same for NIFD and for the 4-Repeat Tauopathy Neuroimaging Initiative (4RTNI), which is also available for download from LONI.
    Controls for NIFD are the same controls as those collected for 4RTNI.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of the **dcm2niix** DICOM to NIfTI converter and **FreeSurfer**.
You can find how to install these software packages on the [installation](../../#installing-clinica-from-source) page.

## Downloading NIFD

To download the NIFD dataset you first need to register to the [LONI Image & Data Archive (IDA)](https://ida.loni.usc.edu/login.jsp), a secure research data repository, and then request access to the NIFD dataset through the submission of an [online application form](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=NIFD).

In order to use the converter, you will need to download both the images and the clinical data.
To do so, from the [main page](https://ida.loni.usc.edu/login.jsp?returnPage=UserManagement.jsp&project=) click on `PROJECTS` and `NIFD`.
To download the imaging data, click on `Download` and choose `Image collections`.
In the `Advanced search` tab, pick the images you wish to download, for example tick `MRI` to download all the MR images, and then click on `SEARCH`.
In the `Advanced search results` tab, click `Select All` and `Add To Collection`.
Finally, in the `Data Collection` tab, select the collection you just created, tick `All` and click on `Advanced download`.
We advise you to group files as 10 zip files.

To download the clinical data, click on `Download` and choose `Study Data`.
Select all the csv files which are present by ticking `ALL` and click `Download`.
You should get two files: `NIFD_Clinical_Data_2017_final_updated.xlsx` and `DataDictionary_NIFD_2017.10.18.xlsx`.
One last file is needed, click on `Download` and choose `Image collections`.
In the `Advanced search` tab, tick all modalities and all boxes from the `Display in result` column.
Click `SEARCH` and `CSV Download`.
You should get a file named `idaSearch_<date>.csv`.
Please rename this file to `idaSearch_all.csv`.

## Modalities supported

Currently, the modalities supported by our converter are:

- T1-weighted MRI
- T2-weighted FLAIR MRI
- Fluorodeoxyglucose (FDG) PET
- Pittsburgh compound B (PiB) PET
- Clinical data

## Using the converter

The converter can be run with the following command line:

```Text
clinica convert nifd-to-bids [OPTIONS] DATASET_DIRECTORY CLINICAL_DATA_DIRECTORY BIDS_DIRECTORY 
```

where:

- `DATASET_DIRECTORY` is the path to the original NIFD images' directory;
- `CLINICAL_DATA_DIRECTORY` is the path to the directory where the following clinical data files are located: `NIFD_Clinical_Data_2017_final_updated.xlsx`, `DataDictionary_NIFD_2017.10.18.xlsx` and `idaSearch_all.csv`;
- `BIDS_DIRECTORY` is the path to the output directory, where the BIDS-converted version of NIFD will be stored.

## Citing this converter in your paper

!!! cite "Example of paragraph:"
    The NIFD data have been curated and converted to the Brain Imaging Data Structure (BIDS) format
    [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)] using Clinica
    [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/collections/NASGJPVL).
