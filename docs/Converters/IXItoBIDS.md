# `ixi-to-bids` – Conversion of Information eXtraction from Images (IXI) to BIDS format
??? quote "Dataset Description"
    The [Information eXtraction from Images](https://brain-development.org/ixi-dataset/)
    is a project which issued a dataset of nearly 600 nifti images from healthy subjects.
    The MR acquisition protocol includes T1,T2, PD weighted, MRA and diffusion-weighted images.
    Three hospitals in London were involved in data collection.

## Downloading the data
The IXI dataset can be downloaded freely from the [IXI webpage](https://brain-development.org/ixi-dataset/).
!!! danger "Organising data with the aim of using the converter"
    The folders and files downloaded from the website should be left **as they are** (name, inner organisation...) but can be
    placed where the user wants them.

## Using the converter
### Available Modalities

!!! info inline end ""
    DTI files are merged together to produce **one** DWI image.
The converter can convert to [BIDS](../glossary.md#bids) all the modalities offered by IXI :
DTI ; T1 ; T2 ; PD ; angiography.


### Dependencies
If you [installed clinica](../Software/Installation.md#install-clinica), this converter needs no further dependencies.

### Understanding the command line
```{ .bash .copy }
clinica convert ixi-to-bids DATASET_DIRECTORY BIDS_DIRECTORY CLINICAL_DATA_DIRECTORY [OPTIONS]
```
where :

- `DATASET_DIRECTORY` is the path to the raw IXI dataset directory, which should contain all the IXI folders previously downloaded :
  ```text title="DATASET_DIRECTORY Organisation"
  DATASET_DIRECTORY
  ├── IXI-T1
  │   ├── IXI002-Guys-0828-T1.nii.gz
  │   ├── ...
  ├── IXI-DTI
  │   └── IXI002-Guys-0828-DTI-00.nii.gz
  │   ├── ...
  ...
  ```

- `BIDS_DIRECTORY` is the path to the BIDS folder to be created
- `CLINICAL_DATA_DIRECTORY` is the path to the folder where the clinical data of the IXI dataset is stored.
```text title="CLINICAL_DATA_DIRECTORY Organisation"
    CLINICAL_DATA_DIRECTORY
    ├── IXI.xls
    ...
```

--8<-- "snippets/converters_options.md"
