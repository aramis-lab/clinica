# `nifd-to-bids` – Conversion of Neuroimaging in Frontotemporal Dementia (NIFD) to BIDS

!!! quote "Description reproduced from the [NIFD's LONI Image & Data Archive (IDA) webpage](https://ida.loni.usc.edu/home/projectPage.jsp?project=NIFD&page=HOME&subPage=OVERVIEW_PR#)"
    NIFD is the nickname for the frontotemporal lobar degeneration neuroimaging initiative (FTLDNI, AG032306), which was funded by the NIA and NINDS to characterize longitudinal clinical and imaging changes in FTLD.
    The imaging and clinical methods are the same for NIFD and for the 4-Repeat Tauopathy Neuroimaging Initiative (4RTNI), which is also available for download from LONI.
    Controls for NIFD are the same controls as those collected for 4RTNI.

## Dependencies

If you only installed the core of Clinica, this pipeline needs the installation of the **dcm2niix** DICOM to NIfTI converter.
You can find how to install these software packages on the [installation](../../#installing-clinica-from-source) page.

## Downloading NIFD

In order to use the converter, you will need to download both the images and the clinical data for NIFD.

First, you will have to register to the [LONI Image & Data Archive (IDA)](https://ida.loni.usc.edu/login.jsp), a secure
research data repository, through the submission of
an [online application form](https://ida.loni.usc.edu/collaboration/access/appApply.jsp?project=NIFD). Once your access
is granted, head to the [NIFD projects's page](https://ida.loni.usc.edu/home/projectPage.jsp?project=NIFD) and login
using your credentials.

### Downloading the imaging data

To download the imaging data:

1. Press `Download` and select `Image collections` in the menubar below.

2. Select the `Advanced Search` tab and configure the appropriate search criteria for your collection. Please, ensure
   only `NIFD` is selected under the `PROJECT/PHASE` section of the form. `MRI` should be selected in the `Modality`
   section by default. Select `PET` and leave the selector to `OR` to download PET imaging data in addition to MRI.

3. Once all search criteria have been selected, press `SEARCH` at the bottom of the form. Your search results should be
   presented in a new tab named `Advanced Search Results`.

4. Select the desired subjects and scans using the tick boxes displayed in the form or press `Select All` in the top
   right corner. Press `Add To Collection` to create a new collection from your selection and give it a name. This name
   will be used as a stem for future downloads.

5. Select the `Data Collections` tab and find your collection within the `My Collections` tree on the left-hand side.

6. Click the `CSV` button to download the collection metadata in tabular form.

7. Select the whole collection by ticking `All` and press `Advanced Download`. A download summary should be displayed
   with the name of the collection, the number of items selected and a dropdown menu to select different groups of file.
   Depending on the size of the collection, it is advised to download the collection in 5 or 10 files instead of 1 as
   default.

8. Upon completion of the download process, create a new folder and move the collection metadata file as well as the
   content of the archive(s) into it.

### Downloading the clinical data

To download the clinical data, press `Download`, select `Study Data`, tick `NIFD Clinical Data` and press download. Once
downloaded, you may move this file to the same location where the imaging data are stored or somewhere else.

## Supported modalities

Currently, the modalities supported by our converter are:

- T1-weighted MRI
- T2-weighted FLAIR MRI
- Fluorodeoxyglucose (FDG) PET
- Pittsburgh compound B (PiB) PET
- Clinical data and survey (MMSE, CDR, ...)

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
