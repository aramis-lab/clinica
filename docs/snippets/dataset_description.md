## The dataset_description.json file

### Specifications

This file MUST be present at the root of a BIDS or CAPS dataset and MUST contain the following minimal information:

```json
{
	"Name": "name identifier for the dataset",
	"BIDSVersion": "1.7.0",
	"CAPSVersion": "1.0.0",
	"DatasetType": "derivative"
}
```

- `Name`: String identifier of the dataset. It can be the name of your study for example. By default Clinica generates a random UUID for this field. When running a pipeline which will create a new CAPS dataset, you can use the `--caps-name` option to provide a name (see [common options](/Software/InteractingWithClinica.md#common-options-for-pipelines)). If the CAPS dataset already exist, the existing name will be kept.
- `BIDSVersion`: The version number of the BIDS specifications that the BIDS input dataset is using when this CAPS dataset was generated.
- `CAPSVersion`: The version number of the CAPS specifications used for this dataset.
- `DatasetType`: Either "raw" or "derivative". For a CAPS dataset this should always be "derivative" as it contains processed data.

In addition, the `dataset_description.json` file MAY contain a `Processing` key which is a list of objects describing the different processing pipelines that were run on this CAPS.
Here is an example for a CAPS dataset containing the outputs of two pipelines: `t1-linear` and `pet-linear`:

```json
{
    "Name": "e6719ef6-2411-4ad2-8abd-da1fd8fbdf32",
    "BIDSVersion": "1.7.0",
    "CAPSVersion": "1.0.0",
    "DatasetType": "derivative",
    "Processing": [
        {
            "Name": "t1-linear",
            "Date": "2024-08-06T10:28:21.848950",
            "Author": "ci",
            "Machine": "ubuntu",
            "InputPath": "/mnt/data_ci/T1Linear/in/bids"
        },
        {
            "Name": "pet-linear",
            "Date": "2024-08-06T10:36:27.403373",
            "Author": "ci",
            "Machine": "ubuntu",
            "InputPath": "/mnt/data_ci/PETLinear/in/bids"
        }
    ]
}
```

A `Processing` is described with the following fields:

- `Name`: The name of the processing. For Clinica pipelines, this is the name of the pipeline.
- `Date`: This date is in iso-format and indicates when the processing was run.
- `Author`: This indicates the user name which triggered the processing.
- `Machine`: This indicates the name of the machine on which the processing was run.
- `InputPath`: This is the full path (on the machine on which the processing was run) to the input dataset of the processing.

### Potential problems

The `dataset_description.json` file for BIDS and CAPS datasets was introduced in Clinica `0.9.0`.

This means that results obtained with prior versions of Clinica do not have this file automatically generated.
Clinica will interpret this as a `<1.0.0` dataset and should error with a suggestion of a minimal `dataset_description.json` that you should add to your dataset.
In this situation, create this new file with the suggested content and re-start the pipeline.

You might also see the following error message:

```
Impossible to write the 'dataset_description.json' file in <FOLDER> because it already exists and it contains incompatible metadata.
```

This means that you have version mismatch for the BIDS and/or CAPS specifications.
That is, the versions indicated in the input (or output) dataset(s) does not match the versions currently used by Clinica.
If this happens, it is recommended to re-run the conversion or the pipeline which initially generated the dataset with the current version of Clinica.
