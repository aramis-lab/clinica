# `ixi-to-bids` – Conversion of Information eXtraction from Images (IXI) to BIDS format

??? quote "Dataset Description"
    [IXI webpage](https://brain-development.org/ixi-dataset/)"
    The Information eXtraction from Images is a project which issued a dataset of nearly
    600 images from healthy subjects. The MR acquisition protocol includes T1,T2, PD weighted,
    MRA and diffusion-weighted images. Three hospitals in London were involved in data collection.


## Downloading IXI
The IXI dataset can be downloaded freely from the [IXI webpage](https://brain-development.org/ixi-dataset/){ data-preview }.
The converter can convert to [BIDS](../glossary.md#bids){ data-preview } all the modalities offered by IXI :
DTI, T1, T2, PD and angiography.

## Using the converter
### Available Modalities
for IXI basically renaming files except DTI where images
are merged into dwi

### Dependencies
If you [installed clinica](../Installation.md#install-clinica), this converter needs no further dependencies.

### Understanding the command line
```bash
clinica convert ixi-to-bids DATASET_DIRECTORY BIDS_DIRECTORY CLINICAL_DATA_DIRECTORY
```
where :

- `DATASET_DIRECTORY` is the path to the raw IXI dataset directory, which should contain all the IXI folders previously downloaded :

```
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
```
CLINICAL_DATA_DIRECTORY
├── IXI.xls
...
```

{!Converters/converters_options.md!}
