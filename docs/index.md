<!-- markdownlint-disable MD007 -->
# Clinica Documentation

TEST

<div style="background-color: #fcf8e5; color: #8C6D46; padding: 16px 24px;">
 <div style="font-size: 16px; font-weight: bold;">⚠ Warning</div>
 <div style="font-size: 12px">
 Issue started on (DATE) and is currently active
 </div>
</div>

## What is Clinica ?

Clinica is a software platform for clinical neuroscience research studies using multimodal data and most often longitudinal follow-up.
You can learn more on [this page](./WhatIsClinica.md).

## Installation

Clinica can be installed on **MacOS** and **Linux** (CentOS or Debian/Ubuntu) machines, and possibly on Windows computers with a Linux Virtual Machine.
We assume that users installing and using Clinica are comfortable using the command line.

- [Installation](./Software/Installation.md)
- [Third-party software](./Software/Third-party.md)  
- [Interacting with Clinica](./Software/InteractingWithClinica.md)
- [Contributing to Clinica](./Software/contributing.md)

## User documentation

### Clinica environment

- [BIDS: the input data structure](./BIDS.md)
- [CAPS: the processed data structure](./CAPS/Introduction.md)

### Pipelines (`clinica run`)

--8<-- "snippets/inventory_pipelines.md"

### Dataset converters (`clinica convert`)

Clinica provides tools to curate several publicly available neuroimaging datasets and convert them to [BIDS](./BIDS.md) namely:

--8<-- "snippets/inventory_converters.md"

!!! note
    We provide converters for the datasets used in the [Aramis Lab](http://www.aramislab.fr/).
    Feel free to contact us if you are interested in another dataset or to contribute!

### I/O tools (`clinica iotools`)

Clinica provides data handling tools for [BIDS](http://bids.neuroimaging.io) and [CAPS](../CAPS/Introduction) compliant datasets:

--8<-- "snippets/inventory_iotools.md"


## Clinica at conferences

--8<-- "snippets/inventory_conferences.md"

## Support

- Check for [past answers](https://groups.google.com/forum/#!forum/clinica-user) in the old Clinica Google Group
- Start a [discussion](https://github.com/aramis-lab/clinica/discussions) on GitHub
- Report an [issue](https://github.com/aramis-lab/clinica/issues) on GitHub

## License

Clinica is distributed under the terms of the MIT license given [here](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).

## Citing Clinica

For publications or communications using Clinica, please cite [[Routier et al., 2021](https://doi.org/10.3389/fninf.2021.689675)]
as well as the references mentioned on the wiki page of the pipelines you used.
Each page includes text to cite the software packages that are used by Clinica (for example, citing SPM when using the `t1-volume` pipeline).

!!! warning "Disclaimer"
    Clinica is a software for research studies.
    It is not intended for use in medical routine.

---

![Clinica_Partners_Banner](./img/Clinica_Partners_Banner.png)
