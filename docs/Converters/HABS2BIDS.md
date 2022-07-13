# `habs-2-bids` – Conversion of Harvard Aging Brain Study (HABS) to BIDS

## Dependencies

This converter does not require additional dependencies beyond Clinica.

## Downloading HABS

This converter requires both the clinical and imaging data available for download on the Harvard Aging Brain
Study [website].

Access to the dataset must be requested through the submission of an [online form]. Please review the _Requested Data_
section carefully and check that both the clinical data and desired imaging modalities are selected. Upon acceptance, a
list of personalized download links should be emailed to you. Keep in mind those links will be set with an expiry date
and the size of dataset is reasonably large (more than 100 GB), so it is advised to complete the download process as
early as possible.

## Supported modalities

At this stage, the HABS to BIDS converter supports the T1, FLAIR and PET modalities. Support for additional modalities
may be implemented in the future.

## Running the converter

Move the downloaded files to a common folder (say `sourcedata`). Do not alter or unzip the content of the original
files. Run the following command to convert the dataset to BIDS format in the target `rawdata` folder:

```sh
clinica convert habs-to-bids /path/to/sourcedata /path/to/rawdata
```

## Citation

!!! cite "Example paragraph:"

    The HABS dataset has been converted to BIDS [[Gorgolewski et al., 2016](https://doi.org/10.1038/sdata.2016.44)]
    using Clinica [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)].


[website]: https://habs.mgh.harvard.edu
[online form]: https://habs.mgh.harvard.edu/researchers/request-data