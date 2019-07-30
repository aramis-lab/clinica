# Installation

You will find below the steps for installing Clinica on Linux or Mac. Please do
not hesitate to contact us on the
[forum](https://groups.google.com/forum/#!forum/clinica-user)
if you encounter any issues.


## Quick start

### Installation of Clinica

#### Python environment
You will need a Python environment to run Clinica. We advise you to
use [Miniconda](http://conda.pydata.org/miniconda.html).
Miniconda allows you to install, run, and update Python packages and their
dependencies. It can also create environments to isolate your libraries.
To install Miniconda, open a new terminal and type the following commands:

- If you are on Linux:
```bash
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda-installer.sh
bash /tmp/miniconda-installer.sh
```

- If you are on Mac:
```bash
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o /tmp/miniconda-installer.sh
bash /tmp/miniconda-installer.sh
```

Miniconda will ask you where to install it. Do not forget to copy the `export
PATH` given at the end of the installation. If everything went
fine, open a new terminal and type `conda info`, it will verify if Conda is
installed, check the version and show your Miniconda path.

Clinica can be installed either by using the conventional
[PyPI package manager](https://pypi.org/project/clinica/) or by using the
[Conda package manager](https://conda.io/docs/). In both cases, it is a good idea
to create an isolated environment if you want to preserve other already installed
Python libraries.  

#### Using Conda

The latest release of Clinica can be installed using Conda as follows:

```bash
conda create --name clinicaEnv python=3.6 clinica -c Aramislab -c conda-forge
```

It will download all the needed packages and create the Clinica environment.

!!! info
    This command will create an environment called `clinicaEnv`. Clinica is
    already installed inside it. The option `-c Aramislab` tells Conda to search
    the package on the [Aramis conda channel](https://anaconda.org/Aramislab).

#### Using pip

The latest release of Clinica can be installed using pip as follows:

```bash
conda create --name clinicaEnv python=3.6
conda activate clinicaEnv
pip install clinica
```

### Installation of the third-party software packages
Depending on the pipeline that you want to use, you need to install
**Pipeline-specific interfaces**. Not all the dependencies are necessary to run
Clinica.
Please refer to [this section](../Third-party) to determine which third-party
libraries you need to install.


### Running the Clinica environment
#### Activation of the Clinica environment

Now that you have created the Clinica environment, you can activate it:

```bash
conda activate clinicaEnv
activate-global-python-argcomplete --user #Only the first time you activate the environment
eval "$(register-python-argcomplete clinica)"
```

!!! success
    Congratulations, you have installed Clinica! At this point, you can try the
    basic `clinica` command and get the help screen:
    ```bash
    (clinicaEnv)$ clinica
    usage: clinica [-v] [-l file.log]  ...

    clinica expects one of the following keywords:

        run                 To run pipelines on BIDS/CAPS datasets.
        convert             To convert unorganized datasets into a BIDS hierarchy.
        iotools             Tools to handle BIDS/CAPS datasets.
        generate            To generate pre-filled files when creating new
                            pipelines (for developers).

    Optional arguments:
      -v, --verbose         Verbose: print all messages to the console
      -l file.log, --logname file.log
                            Define the log file name (default: clinica.log)
    ```

    If you have successfully installed the third-party software, you are ready
    to run any of the pipelines proposed by Clinica.

    It should display a help describing the different categories of command line
    (see [Interacting with Clinica](../InteractingWithClinica) for further
    explanations).

#### Deactivation of the Clinica environment
At the end of your session, remember to deactivate your Conda environment:
```bash
conda deactivate
```

## Developer installation

If you plan to contribute to Clinica or if you want to have the current development
version, you can either:

* Download the tarball for a specific version from our
[repository](https://github.com/aramis-lab/clinica/releases).
Then decompress it.
* Clone Clinica's repository from GitHub:
```bash
git clone https://github.com/aramis-lab/clinica.git
```

We suggest creating a custom Conda environment and installing Clinica using the
provided YML file:

```bash
cd clinica
conda env create -f environment.yml
```

By default, the environment is named `clinica_env`. You can choose a different
name by adding the option `--name my_clinica_environment`.

Clinica is installed within the environment created. Remember to
activate the environment before proceeding:

```bash
conda activate clinica_env
pip install -e . # Only the first time you activate the environment
activate-global-python-argcomplete --user # Only the first time you activate the environment
eval "$(register-python-argcomplete clinica)"
```

If everything goes well, type `clinica` and you should see the help message which
is displayed above.

At the end of your session, you can deactivate your Conda environment:
```bash
conda deactivate
```

Remember that Clinica will be only available inside your Conda environment.
Further information for Clinica's contributors can be found
[here](./CodingForClinica).
