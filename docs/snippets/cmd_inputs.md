--8<-- [start:bids]
- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](../BIDS.md) hierarchy
--8<-- [end:bids]

--8<-- [start:caps]
- `CAPS_DIRECTORY` is the output folder containing the results in a [CAPS](../CAPS/Introduction.md) hierarchy
--8<-- [end:caps]

--8<-- [start:bids_caps]
- `BIDS_DIRECTORY` is the input folder containing the dataset in a [BIDS](../BIDS.md) hierarchy
- `CAPS_DIRECTORY` is the output folder containing the results in a [CAPS](../CAPS/Introduction.md) hierarchy
--8<-- [end:bids_caps]

--8<-- [start:acq]
- `ACQ_LABEL` is the label given to the PET acquisition, specifying the tracer used (`trc-<acq_label>`).
--8<-- [end:acq]

--8<-- [start:group]
- `GROUP_LABEL` is the label of the group that is associated to the DARTEL template that you had created when running the [`t1-volume`](./T1_Volume.md) pipeline.
--8<-- [end:group]

--8<-- [start:region]
- The reference region is used to perform intensity normalization (i.e. dividing each voxel of the image by the average uptake in this region) resulting in a standardized uptake value ratio ([SUVR](../glossary.md#suvr)) map.
  It can be `cerebellumPons` or `cerebellumPons2` (used for amyloid tracers) and `pons` or `pons2` (used for FDG).
  See [PET introduction](./PET_Introduction.md) for more details about masks versions.
--8<-- [end:region]
