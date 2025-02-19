# `describe` - Describe a BIDS or CAPS dataset

This tool describes a BIDS or CAPS dataset, basically parsing the `dataset_description.json` file and displaying the information in the console in a nice way.

!!! note
    This tool has been added in Clinica `0.9.0` and is not available in older versions.

```shell
clinica iotools describe DATASET_PATH
```

where:

- `DATASET_PATH`: path to the BIDS or CAPS dataset to be described.

**Examples:**

With a BIDS dataset:

```shell
clinica iotools describe ./bids
                         Dataset information
╭───────────────────────────────┬──────┬─────────────────────────────╮
│                          Name │ Type │ BIDS Specifications Version │
├───────────────────────────────┼──────┼─────────────────────────────┤
│ The mother of all experiments │  raw │ 1.6.0                       │
╰───────────────────────────────┴──────┴─────────────────────────────╯
```

With a CAPS dataset:

```shell
clinica iotools describe ./caps
                                               Dataset information
╭──────────────────────────────────────┬────────────┬─────────────────────────────┬─────────────────────────────╮
│                                 Name │       Type │ BIDS Specifications Version │ CAPS Specifications Version │
├──────────────────────────────────────┼────────────┼─────────────────────────────┼─────────────────────────────┤
│ 05c5794b-2d20-4217-b727-c215d079ab35 │ derivative │ 1.9.0                       │                       1.0.0 │
╰──────────────────────────────────────┴────────────┴─────────────────────────────┴─────────────────────────────╯
                                                                                           Processing information
╭────────────────────────────┬────────────────────────────┬───────────────────┬─────────────────┬────────────────────────────────────────────┬────────────────┬────────────────────────────────────────────╮
│                       Name │                       Date │            Author │         Machine │                                  InputPath │ ClinicaVersion │                               Dependencies │
├────────────────────────────┼────────────────────────────┼───────────────────┼─────────────────┼────────────────────────────────────────────┼────────────────┼────────────────────────────────────────────┤
│ dwi-preprocessing-using-t1 │ 2024-09-01 13:41:05.529799 │ nicolas.gensollen │ UMR-COLLI-MP050 │ /Users/nicolas.gensollen/dwi_datasets/bids │          0.9.0 │                Dependencies                │
│                            │                            │                   │                 │                                            │                │ ╭───────────┬──────────────┬─────────────╮ │
│                            │                            │                   │                 │                                            │                │ │      Name │ VersionCons… │ InstalledV… │ │
│                            │                            │                   │                 │                                            │                │ ├───────────┼──────────────┼─────────────┤ │
│                            │                            │                   │                 │                                            │                │ │      ants │      >=2.3.0 │       2.5.0 │ │
│                            │                            │                   │                 │                                            │                │ │       fsl │      >=6.0.1 │       6.0.2 │ │
│                            │                            │                   │                 │                                            │                │ │    mrtrix │      >=0.0.0 │       3.0.3 │ │
│                            │                            │                   │                 │                                            │                │ │ convert3d │      >=0.0.0 │       1.0.0 │ │
│                            │                            │                   │                 │                                            │                │ ╰───────────┴──────────────┴─────────────╯ │
╰────────────────────────────┴────────────────────────────┴───────────────────┴─────────────────┴────────────────────────────────────────────┴────────────────┴────────────────────────────────────────────╯
```
