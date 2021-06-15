# What is Clinica?

Clinica is a software platform for clinical neuroscience research studies using multimodal data and most often longitudinal follow-up.

## What do you mean by "clinical neuroscience studies"?

While the meaning of these terms may vary, we use it to refer to studies involving human participants (i.e. patients with neurological and psychiatric diseases, and control subjects) explored with multimodal data (neuroimaging, clinical and cognitive evaluations, genetic data, etc.) and most often involving longitudinal follow-up.

Currently, a large part of Clinica is devoted to neuroimaging data analysis.
Future versions will include more extensive support of genotyping data and clinical/cognitive evaluations.

## What are the main features of Clinica?

- Neuroimaging data analysis
  - Anatomical MRI
  - Diffusion MRI
  - PET
- Statistical analysis
- Machine learning
- Data management made easy:
  - Standardized data structures for inputs
  - Standardized data structures for outputs
  - Conversion of datasets

## Which technologies underlie Clinica?

Clinica is written in Python and uses the [Nipype](http://nipype.readthedocs.io/en/latest/) system for pipelining.
It combines widely-used software packages for neuroimaging data analysis ([SPM12](http://www.fil.ion.ucl.ac.uk/spm/), [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/), [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/), [MRtrix3](http://www.mrtrix.org/), etc.), machine learning ([Scikit-learn](http://scikit-learn.org/)) and the [BIDS](http://bids.neuroimaging.io/) standard for data organization.

## What are the benefits of using Clinica and not simply those widely-used software packages?

In short: **to make your life easier**.

Specifically, Clinica provides:

- complex processing pipelines involving the combination of different analysis software tools;
- integration between feature extraction and statistics/machine learning;
- standardized file organization.

This should help you to:

- easily share data and results within your institution and with external collaborators;
- make your research more reproducible;
- spend less time on data management and processing.

## Who develops Clinica?

Clinica was initiated by the [Aramis Lab](http://www.aramislab.fr/) at the [Brain and Spine Institute (ICM)](https://icm-institute.org/) in Paris.
