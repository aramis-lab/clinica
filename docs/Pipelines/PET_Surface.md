# `pet-surface` - Surface-based processing of PET images

This pipeline performs several processing steps for the analysis of PET data on the cortical surface [[Marcoux et al., 2018](https://doi.org/10.3389/fninf.2018.00094)]:

- co-registration of PET and T1-w MRI (T1) images;
- intensity normalization;
- partial volume correction;
- robust projection of the PET signal onto the subject’s cortical surface;
- spatial normalization to a template;
- atlas statistics.

This pipeline relies mainly on tools from **[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)** and **[PETPVC](https://github.com/UCL/PETPVC)** [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)].

## Prerequisite
You need to have performed the [`t1-freesurfer`](../T1_FreeSurfer) pipeline on your T1 images.


## Dependencies
<!-- If you installed the docker image of Clinica, nothing is required.-->

If you only installed the core of Clinica, this pipeline needs the installation of **FreeSurfer 6.0**, **SPM12**, **FSL 6.0** and **PETPVC 1.2.4** (which depends on **ITK 4**) on your computer. You can find how to install these software packages on the [third-party](../../Third-party) page.

## Running the pipeline
The pipeline can be run with the following command line:

```Text
clinica run pet-surface <bids_directory> <caps_directory>
```
where:

  - `bids_directory` is the input folder containing the dataset in a [BIDS](../../BIDS) hierarchy.
  - `caps_directory` is the output folder containing the results in a [CAPS](../../CAPS/Introduction) hierarchy.

If you want to run the pipeline on a subset of your BIDS dataset, you can use the `-tsv` flag to specify in a TSV file the participants belonging to your subset.

Please note that next to each PET file in your BIDS folder, a `json` file must be added to specify the `EffectiveResolutionInPlane` and `EffectiveResolutionAxial` in mm relative to the point spread function (PSF).

Your bids hierarchy (for a given subject `sub-001`) must look like this :
```
bids
└── sub-001
    └── ses-M00
        ├── anat
        │   └── sub-001_ses-M00_T1w.nii.gz
        └── pet
            ├── sub-001_ses-M00_task-rest_acq-fdg_pet.json
            └── sub-001_ses-M00_task-rest_acq-fdg_pet.nii.gz
```

The `sub-001_ses-M00_task-rest_acq-fdg_pet.json` must look like this :

```
{
	"Psf":[
		{
			"EffectiveResolutionInPlane": 5.5,
			"EffectiveResolutionAxial": 5.5
		}
	]
}
```


Pipeline options

- `--pet_tracer`: type of PET image to process. Possible values are fdg and av45. Default value is fdg. This parameter affects the reference region used for the intensity normalization (FDG: pons, AV45: pons and cerebellum).

- `-np`: This parameter specifies the number of threads to run in parallel. We recommand using `your_number_of_cpu - 1`. Please note that PETPVC is extremely demanding in terms of resources and may cause the pipeline to crash if many subjects happen to be partial volume corrected at the same time (Error : `Failed to allocate memory for image`). To mitigate this issue, you can do the following:

**1)** Use a working directory when you launch clinica

**2)** If the pipeline crash, just launch again the command (while giving the same working directory)

**3)** The whole processing will continue where it left ! (you can lower the number of thread to run in parallel the second time)

!!! note
    The arguments common to all Clinica pipelines are described in [Interacting with clinica](../../InteractingWithClinica).

!!! tip
    Do not hesitate to type `clinica run pet-surface --help` to see the full list of parameters.

## Outputs

Results are stored in the following folder of the [CAPS hierarchy](../../CAPS/Specifications/#pet-surface-surface-based-processing-of-pet-images): `subjects/sub-<participant_label>/ses-<session_label>/pet/surface`

The files are (where `*` stands for `sub-<participant_label>_ses-<session_label>`):

- `atlas_statistics/*_task-<label>_acq-<label>_pet_space-<label>_pvc-iy_suvr-<label>_statistics.tsv`: TSV files summarizing the regional statistics on the labelled atlases (Desikan and Destrieux).
- `*_hemi-{left|right}_midcorticalsurface`: surface at equal distance between the white matter/gray matter interface and the pial surface (one per hemisphere).
- `*_task-rest_acq-<label>_pet_space-<label>_suvr-<label>_pvc-iy_hemi-<label>_fwhm-<value>_projection.mgh`: PET data that can be mapped onto meshes. If the `space` is `fsaverage`, it can be mapped either onto the white or pial surface of FsAverage. If the `space` is `native`, it can be mapped onto the white or pial surface of the subject’s surface (i.e. `{l|r}h.white`, `{l|r}h.pial` files from the `t1-freesurfer` pipeline).

!!! note
    The full list of output files from the pet-volume pipeline can be found in the [The ClinicA Processed Structure (CAPS) specifications](../../CAPS/Specifications/#pet-surface-surface-based-processing-of-pet-images).


<center>![PET surface results](../../img/PET_Surface.jpg)</center>
*<center><small>FDG PET SUVR projected onto the cortical surface (left hemisphere) for (from left to right) a cognitively normal subject (CN), a patient with Alzheimer’s disease (AD), a patient with semantic variant primary progressive aphasia (svPPA) and a patient with logopenic variant primary progressive aphasia (lvPPA). The first row is the projection in the subject’s space. The second row is the same signal for each subject, but warped to FsAverage after smoothing with a 20 mm Gaussian kernel.</small></center>*



## Going further

- You can use projected PET data to perform group comparison or correlation analysis with the [`statistics-surface` pipeline](../Stats_Surface).
- You can use projected PET data to perform classification based on [machine learning](../MachineLearning_Classification), as showcased in the [AD-ML framework](https://github.com/aramis-lab/AD-ML).


## Describing this pipeline in your paper

!!! cite "Example of paragraph:"
    These results have been obtained using the `pet-surface` pipeline of Clinica [[Routier et al](https://hal.inria.fr/hal-02308126/); [Marcoux et al., 2018](https://doi.org/10.3389/fninf.2018.00094)]. The subject’s PET image was registered to the T1 using spmregister ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)) and intensity normalized using the [pons | pons and cerebellum] from the Pick atlas in MNI space as reference region (registration to MNI space was performed using [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)). Partial volume correction was then performed using the iterative Yang algorithm implemented in [PETPVC](https://github.com/UCL/PETPVC) [[Thomas et al., 2016](https://doi.org/10.1088/0031-9155/61/22/7975)] with regions obtained from gtmseg ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)). Based on the subject’s white surface and cortical thickness, seven surfaces for each hemisphere were computed, ranging from 35% to 65% of the grey matter thickness. The partial volume corrected data were projected onto these meshes and the seven values were averaged, giving more weight to the vertices near the center of the cortex. Finally, the projected PET signal in the subject’s native space was spatially normalized to the standard space of FsAverage ([FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)).

!!! tip
    Easily access the papers cited on this page on [Zotero](https://www.zotero.org/groups/2240070/clinica_aramislab/items/collectionKey/RGVVHC5W).

## Support

-   You can use the [Clinica Google Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!
-   Report an issue on [GitHub](https://github.com/aramis-lab/clinica/issues).


## Appendix I: Diagram of the pipeline execution

<center>![Diagram of the pipeline execution](https://www.frontiersin.org/files/Articles/400789/fninf-12-00094-HTML-r1/image_m/fninf-12-00094-g001.jpg)</center>
*<center><small>The subject's T1-w MRI is coregistered with the PET image and the PET image is intensity normalized using the average uptake in a reference region. In parallel, cortical surfaces and a parcellation are generated from the subject's T1-w MRI. The PET image, after partial volume correction performed using the parcellation, is robustly projected onto the cortical surface. Finally, regional mean uptake values are extracted from the projected PET data, and the projected PET signal in the subject's native space is spatially normalized to the standard space of FsAverage.</small></center>*



## Appendix II: How to manipulate outputs

Outputs of pipeline are composed of two different type of file: surface files and MGH data that are to be overlaid onto a surface.  

### Surface file
They can be read using various software. You can open it using `freeview` (FreeSurfer viewer), with `freeview -f /path/to/your/surface/file`.

You can also open it in MATLAB, using SurfStat: `mysurface = SurfStatReadSurf('/path/to/your/surface/file')`. This will give you a structure with fields `coord` (for coordinates), a list of coordinates for each point of the mesh, and also field `tri` (for triangle), a list of triplet for each triangle, each number representing the Nth vertex of the `coord` list. Here is below an example to make things clearer (read with Matlab)

<center>![PET surface file](../../img/PET_Surface_File.png)</center>


### Data files

The data files wear the `.mgh` extension. It is composed of a single vector. This file contains a vector, where value at Nth position must be mapped into the Nth vertex of the `coord` list to be correctly represented. You can access them either in Matlab with the command:
```
mydata = SurfStatReadData('/path/to/your/file.mgh');
```
(you will get a single row vector)

Or in Python with the `nibabel` library:
```
import nibabel
mydata = nibabel.load('/path/to/your/mgh/file')
```
`mydata` will then be a `MGHImage`, more information [here](http://nipy.org/nibabel/reference/nibabel.freesurfer.html#mghimage). Keep in mind that if you want to manipulate the data vector within this object, you will need to transform it a bit. Indeed, if you do the following:

```
raw_data = mydata.get_data()
print(raw_data.shape)
```

The shape of your "raw" vector will probably look like this: `(163842, 1, 1)`. Use the `squeeze` function from `numpy` to get a `(163842,)` shape. (documentation [here](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.squeeze.html)). The reverse operation (`(163842,)` to a `(163842, 1, 1)` shape) can be achieved with the `atleast_3d` function from `numpy` (documentation [here](https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.atleast_3d.html)). This may come handy when you need to create a `MGHImage` from scratch.

### Visualization of the results
After the execution of the pipeline, you can check the outputs of one subject by running this command (subject moved into FsAverage):

```bash
freeview -f $SUBJECTS_DIR/fsaverage/surf/lh.pial:overlay=path/to/your/projected/pet/in/fsaverage/left/hemi \
 -f $SUBJECTS_DIR/fsaverage/surf/rh.pial:overlay=path/to/your/projected/pet/in/fsaverage/right/hemi
```
<center>![PET-Surface Freeview](../../img/PET_Surface_Freeview.png)</center>


But you can also visualize your subjects cortical projection directly into his native space:
```bash
freeview -f path/to/midcortical/surface/left:overlay=path/to/your/projected/pet/in/nativespace/left/hemi \
 -f path/to/midcortical/surface/right:overlay=path/to/your/projected/pet/in/nativespace/right/hemi
```

You will need to adjust the colormap using the `Configure` button in the left panel, just below the `Overlay` section.

You can also visualize your surface using the SurfStat tool. Once SurfStat installation folder is added to your MATLAB path, you can display your surfaces with the following commands:
```
mydata = SurfStatReadData({'/path/to/left/data', '/path/to/right/data'});
mysurfaces = SurfStatReadSurf({'/path/to/left/surface', '/path/to/right/surface'});
figure, SurfStatViewData(mydata, mysurfaces, 'Title of figure');
```
It will get you the following figure:
<center>![PET-Surface SurfStat](../../img/PET_Surface_SurfStat.png)</center>
