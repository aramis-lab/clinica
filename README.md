<!--(http://www.clinica.run/img/clinica_brainweb.png)-->
<!-- markdownlint-disable MD033 -->

<h1 align="center">
  <a href="http://www.clinica.run">
    <img src="http://www.clinica.run/assets/images/clinica-icon-257x257.png" alt="Logo" width="120" height="120">
  </a>
  <br/>
  Clinica
</h1>

<p align="center"><strong>Software platform for clinical neuroimaging studies</strong></p>

<p align="center">
  <a href="https://ci.inria.fr/clinica-aramis/job/clinica/job/dev/">
    <img src="https://ci.inria.fr/clinica-aramis/buildStatus/icon?job=clinica%2Fdev" alt="Build Status">
  </a>
  <a href="https://badge.fury.io/py/clinica">
    <img src="https://badge.fury.io/py/clinica.svg" alt="PyPI version">
  </a>
  <a href="https://pypi.org/project/clinica">
    <img src="https://img.shields.io/pypi/pyversions/clinica" alt="Supported Python versions">
  </a>
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/Installation/">
  </a>
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/Installation/">
    <img src="https://anaconda.org/aramislab/clinica/badges/platforms.svg" alt="platform">
  </a>
  <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black">
  </a>
</p>

<p align="center">
  <a href="http://www.clinica.run">Homepage</a> |
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/">Documentation</a> |
  <a href="https://doi.org/10.3389/fninf.2021.689675">Paper</a> |
  <a href="https://github.com/aramis-lab/clinica/discussions">Forum</a> |
  See also:
  <a href="#related-repositories">AD-ML</a>,
  <a href="#related-repositories">AD-DL</a>,
  <a href="#related-repositories">ClinicaDL</a>
</p>

## About The Project

Clinica is a software platform for clinical research studies involving patients
with neurological and psychiatric diseases and the acquisition of multimodal
data (neuroimaging, clinical and cognitive evaluations, genetics...),
most often with longitudinal follow-up.

Clinica is command-line driven and written in Python.
It uses the [Nipype](https://nipype.readthedocs.io/) system for pipelining and combines
widely-used software packages for neuroimaging data analysis
([ANTs](http://stnava.github.io/ANTs/),
[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/),
[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki),
[MRtrix](https://www.mrtrix.org/),
[PETPVC](https://github.com/UCL/PETPVC),
[SPM](https://www.fil.ion.ucl.ac.uk/spm/)), machine learning
([Scikit-learn](https://scikit-learn.org/stable/)) and the [BIDS
standard](http://bids-specification.readthedocs.io/) for data organization.

Clinica provides tools to convert publicly available neuroimaging datasets into
BIDS, namely:

- [ADNI: Alzheimer’s Disease Neuroimaging Initiative](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/ADNI2BIDS/)
- [AIBL: Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/AIBL2BIDS/)
- [HABS: Harvard Aging Brain Study](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/HABS2BIDS/)
- [NIFD: Neuroimaging in Frontotemporal Dementia](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/NIFD2BIDS/)
- [OASIS: Open Access Series of Imaging Studies](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/OASIS2BIDS/)
- [OASIS-3: Longitudinal Neuroimaging, Clinical, and Cognitive Dataset for Normal Aging and Alzheimer’s Disease](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/OASIS3TOBIDS/)

Clinica can process any BIDS-compliant dataset with a set of complex processing
pipelines involving different software packages for the analysis of
neuroimaging data (T1-weighted MRI, diffusion MRI and PET data).
It also provides integration between feature extraction and statistics, machine
learning or deep learning.

![ClinicaPipelines](http://www.clinica.run/img/Clinica_Pipelines_A4_2021-04-02_75dpi.jpg)

Clinica is also showcased as a framework for the reproducible classification of
Alzheimer's disease using
[machine learning](https://github.com/aramis-lab/AD-ML) and
[deep learning](https://github.com/aramis-lab/clinicadl).

## Getting Started

> Full instructions for installation and additional information can be found in
the [user documentation](https://aramislab.paris.inria.fr/clinica/docs/public/latest/).

Clinica currently supports macOS and Linux.
It can be installed by typing the following command:

```sh
pip install clinica
```

To avoid conflicts with other versions of the dependency packages installed by pip, it is strongly recommended to create a virtual environment before the installation.
For example, use [Conda](https://docs.conda.io/en/latest/miniconda.html), to create a virtual
environment and activate it before installing clinica (you can also use
`virtualenv`):

```sh
conda create --name clinicaEnv python=3.8
conda activate clinicaEnv
```

Depending on the pipeline that you want to use, you need to install pipeline-specific interfaces.
Not all the dependencies are necessary to run Clinica.
Please refer to this [page](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Third-party/)
to determine which third-party libraries you need to install.

## Example

Diagram illustrating the Clinica pipelines involved when performing a group
comparison of FDG PET data projected on the cortical surface between patients
with Alzheimer's disease and healthy controls from the ADNI database:

![ClinicaExample](http://www.clinica.run/img/Clinica_Example_2021-04-02_75dpi.jpg)

1. Clinical and neuroimaging data are downloaded from the ADNI website and data
   are converted into BIDS with the [`adni-to-bids`
   converter](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/ADNI2BIDS/).
2. Estimation of the cortical and white surface is then produced by the
   [`t1-freesurfer`
   pipeline](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_FreeSurfer/).
3. FDG PET data can be projected on the subject’s cortical surface and
   normalized to the FsAverage template from FreeSurfer using the
   [`pet-surface` pipeline](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Surface/).
4. TSV file with demographic information of the population studied is given to
   the [`statistics-surface`
   pipeline](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/Stats_Surface/) to generate
   the results of the group comparison.

> For more examples and details, please refer to the
> [Documentation](https://aramislab.paris.inria.fr/clinica/docs/public/latest/).

## Support

- Check for [past answers](https://groups.google.com/forum/#!forum/clinica-user) in the old Clinica Google Group
- Start a [discussion](https://github.com/aramis-lab/clinica/discussions) on Github
- Report an [issue](https://github.com/aramis-lab/clinica/issues) on GitHub

## Contributing

We encourage you to contribute to Clinica!
Please check out the [Contributing to Clinica guide](CONTRIBUTING.md) for
guidelines about how to proceed.  Do not hesitate to ask questions if something
is not clear for you, report an issue, etc.

## License

This software is distributed under the MIT License.
See [license file](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt)
for more information.

## Citing us

- Routier, A., Burgos, N., Díaz, M., Bacci, M., Bottani, S., El-Rifai O., Fontanella, S., Gori, P., Guillon, J., Guyot, A., Hassanaly, R., Jacquemont, T.,  Lu, P., Marcoux, A.,  Moreau, T., Samper-González, J., Teichmann, M., Thibeau-Sutre, E., Vaillant G., Wen, J., Wild, A., Habert, M.-O., Durrleman, S., and Colliot, O.:
*Clinica: An Open Source Software Platform for Reproducible Clinical Neuroscience Studies* Frontiers in Neuroinformatics, 2021
[doi:10.3389/fninf.2021.689675](https://doi.org/10.3389/fninf.2021.689675)

## Related Repositories

- [AD-DL: Classification of Alzheimer's disease status with convolutional neural networks](https://github.com/aramis-lab/AD-DL).
- [AD-ML: Framework for the reproducible classification of Alzheimer's disease using
machine learning](https://github.com/aramis-lab/AD-ML).
- [ClinicaDL: Framework for the reproducible processing of neuroimaging data with deep learning methods](https://github.com/aramis-lab/clinicadl).
