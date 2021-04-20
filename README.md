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
  <a href="https://ci.inria.fr/clinica-aramis/job/clinica/job/master/">
    <img src="https://ci.inria.fr/clinica-aramis/buildStatus/icon?job=clinica%2Fmaster" alt="Build Status">
  </a>
  <a href="https://badge.fury.io/py/clinica">
    <img src="https://badge.fury.io/py/clinica.svg" alt="PyPI version">
  </a>
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/Installation/">
    <img src="https://anaconda.org/aramislab/clinica/badges/installer/conda.svg" alt="conda install">
  </a>
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/Installation/">
    <img src="https://anaconda.org/aramislab/clinica/badges/platforms.svg" alt="platform">
  </a>
</p>

<p align="center">
  <a href="http://www.clinica.run">Homepage</a> |
  <a href="https://aramislab.paris.inria.fr/clinica/docs/public/latest/">Documentation</a> |
  <a href="https://hal.inria.fr/hal-02308126">Preprint</a> |
  <a href="https://groups.google.com/forum/#!forum/clinica-user">Forum</a> |
  See also:
  <a href="#related-repositories">AD-ML</a>,
  <a href="#related-repositories">AD-DL</a>
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
- [NIFD: Neuroimaging in Frontotemporal Dementia](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/NIFD2BIDS/)
- [OASIS: Open Access Series of Imaging Studies](https://aramislab.paris.inria.fr/clinica/docs/public/latest/Converters/OASIS2BIDS/)

Clinica can process any BIDS-compliant dataset with a set of complex processing
pipelines involving different software packages for the analysis of
neuroimaging data (T1-weighted MRI, diffusion MRI and PET data).
It also provides integration between feature extraction and statistics, machine
learning or deep learning.

![ClinicaPipelines](http://www.clinica.run/img/Clinica_Pipelines_A4_2021-04-02_75dpi.jpg)

Clinica is also showcased as a framework for the reproducible classification of
Alzheimer's disease using
[machine learning](https://github.com/aramis-lab/AD-ML) and
[deep learning](https://github.com/aramis-lab/AD-DL).

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
conda create --name clinicaEnv python=3.7
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

- [Report an issue on GitHub](https://github.com/aramis-lab/clinica/issues)
- Use the [Clinica Google
  Group](https://groups.google.com/forum/#!forum/clinica-user) to ask for help!

<!--
## Contributing
We encourage you to contribute to Clinica! Please check out the [Contributing
to Clinica guide](Contributing.md) for guidelines about how to proceed. Do not
hesitate to ask questions if something is not clear for you, report an issue,
etc.
-->

## License

This software is distributed under the MIT License.
See [license file](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt)
for more information.

## Related Repositories

- [AD-DL: Framework for the reproducible classification of Alzheimer's disease using
deep learning](https://github.com/aramis-lab/AD-DL)
- [AD-ML: Framework for the reproducible classification of Alzheimer's disease using
machine learning](https://github.com/aramis-lab/AD-ML)
