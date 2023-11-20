<!-- markdownlint-disable MD046 -->
# Installation

You will find below the steps for installing Clinica on Linux or Mac.
Please do not hesitate to contact us on the
[forum](https://groups.google.com/forum/#!forum/clinica-user) or
[GitHub](https://github.com/aramis-lab/clinica/issues)
if you encounter any issues.

## Prepare your Python environment

You will need a Python environment to run Clinica.

We advise you to use [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
Miniconda allows you to install, run, and update Python packages and their dependencies.
It can also create environments to isolate your libraries.
To install Miniconda, open a new terminal and type the following commands:

- If you are on Linux:

```{.sourceCode .bash}
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda-installer.sh
bash /tmp/miniconda-installer.sh
```

- If you are on Mac:

```{.sourceCode .bash}
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o /tmp/miniconda-installer.sh
bash /tmp/miniconda-installer.sh
```

[Miniconda](https://docs.conda.io/en/latest/miniconda.html) will ask you where to install it.
Do not forget to copy the `export PATH` given at the end of the installation.
If everything went fine, open a new terminal and type `conda info`, it will verify if
Conda is installed, check the version and show your Miniconda path.

## Install Clinica

The latest release of Clinica can be installed by using the conventional
[PyPI package manager](https://pypi.org/project/clinica/) as follows:

```shell
conda create --name clinicaEnv python=3.10
conda activate clinicaEnv
pip install clinica
```

!!! info
    Since Clinica `v0.3.5`, Conda installation is no longer available (i.e.
    `conda create --name clinicaEnv python=3.6 clinica -c Aramislab -c conda-forge`
    will only install Clinica `v0.3.4`).
    Pip is now the only way to install the latest version of Clinica.

## Installation of the third-party software packages

Depending on the pipeline that you want to use, you need to install **pipeline-specific interfaces**.
Not all the dependencies are necessary to run Clinica.
Please refer to [this section](../Third-party) to determine which third-party libraries you need to install.

## Shell completion (optional)

Shell completion for Clinica is available for [Bash](https://www.gnu.org/software/bash/),
[Fish](https://fishshell.com/docs/current/), and [Zsh](https://zsh.sourceforge.io/Doc/).

For Bash, add this to `~/.bashrc`:

```shell
eval "$(_CLINICA_COMPLETE=source_bash clinica)"
```

For Fish, add this to `~/.config/fish/completions/clinica.fish`:

```shell
eval (env _CLINICA_COMPLETE=source_fish clinica)
```

For Zsh, add this to `~/.zshrc`:

```shell
eval "$(_CLINICA_COMPLETE=source_zsh clinica)"
```

Finally, open a new shell to enable completion.

## Run the Clinica environment

### Activation of the Clinica environment

Now that you have created the Clinica environment, you can activate it:

```{.sourceCode .bash}
conda activate clinicaEnv
```

!!! success
    Congratulations, you have installed Clinica! At this point, you can try the
    basic `clinica` command and get the following help screen:

    ```console
    (clinicaEnv)$ clinica
    Usage: clinica [OPTIONS] COMMAND [ARGS]...

    Options:
      --version      Show the version and exit.
      -v, --verbose  Increase logging verbosity.
      -h, --help     Show this message and exit.

    Commands:
      convert   Convert popular neuroimaging datasets to the BIDS format.
      generate  Instantiate a new pipeline from available templates.
      iotools   Tools to handle BIDS/CAPS datasets.
      run       Run pipelines on BIDS and CAPS datasets.
    ```

If you have successfully installed the third-party software packages, you are
ready to run any of the pipelines proposed by Clinica.

You can now learn how to [interact with Clinica](../InteractingWithClinica).

### Deactivation of the Clinica environment

At the end of your session, remember to deactivate your Conda environment:

```{.sourceCode .bash}
conda deactivate
```

## Developer instructions

This section is intended for users who plan to contribute to Clinica or test the current development version.

Clinica uses [Poetry](https://python-poetry.org) to manage its development environment.
Please follow these [installation instructions](https://python-poetry.org/docs/#installation) and verify the `poetry` command is correctly setup.

Clone the development branch of Clinica:

```shell
git clone --branch dev https://github.com/aramis-lab/clinica.git
cd clinica
```

Create an environment for development:

```shell
conda env create -f environment.yml
```

Install Clinica with the necessary development dependencies:

```shell
conda activate clinica_env
poetry install
```

In case you need to test the documentation locally, install the additional dependencies with:

```shell
poetry install --extras docs
```
