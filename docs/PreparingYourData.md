# Preparing your data

## How to prepare your data?

In order to use Clinica pipelines by command line, you will need to give as input:

1. A dataset that follows the [BIDS standard](http://bids.neuroimaging.io/).
2. A file (named by default ``subject_sessions_list.tsv``) with the list of all the available sessions for each subject.

This page explains how to obtain a [BIDS](glossary.md#bids) compliant dataset and the file `subject_sessions_list.tsv`.

## Dataset BIDS compliant

[BIDS](http://bids.neuroimaging.io/) is a standard for organizing neuroimaging data, and it has been also adopted by the Aramis team.

The pipelines used by command line work only if you give as input a BIDS compliant dataset.
So, if your dataset doesn't follow this standard, you will need to convert to it.

On BIDS website is possible to find some tools that can help you for the conversion like, for example, the [OpenfMRI to BIDS converter](https://github.com/INCF/openfmri2bids).
A quick overview of the BIDS structure can be found on the [BIDS](BIDS) page.
For more details you can also check the documentation.

> Note:
>
> This is a future development and is still work in progress.

Clinica offers the possibility to convert the ADNI dataset to BIDS standard using the following command:

```Text
clinica convert adni_to_bids dataset_directory output directory -c
```

## Generate subject_sessions_list.tsv file

Together with a BIDS compliant dataset, you will need to give as input a TSV file with two columns (participant_id and session_id) containing the list of visits for each subject:

```text
participant_id   session_id
sub-01           ses-M000
sub-01           ses-M024
sub-02           ses-M024
```

You can generate this file using Clinica with command `create-subjects-visit`:

```bash
clinica io create-subjects-visits bids_directory output_directory
```

More detail about the command can be found [here](./IOTools/create_subjects_visits.md).
