# archground-dist

Python package to perform the analysis of planetary spectra from ground-based facilities using the archnemesis algorithm

This package is currently maintained by [Juan Alday](https://gapt.iaa.es/content/juan-alday).

If interested users are missing key points in the documentation, would appreciate seeing jupyter notebooks for certain purposes, or want to report issues, please do so by contacting us.

## Installation

To install this package, please follow the next steps:

### 1. Installing the latest version of the code from GitHub

The latest version of code has to be downloaded from [Github](https://github.com/juanaldayparejo/archground-dist.git) under a [GNU General Public License v3](LICENSE). To do so, type in the command window:

```bash
git clone https://github.com/juanaldayparejo/archground-dist.git
```

Before installing the library, we recommend users to create and load a new Python [virtual environment](https://docs.python.org/3/library/venv.html) for a clean install:

```bash
python -m venv name_of_virtual_environment/
source name_of_virtual_environment/bin/activate
```

Then move into the package directory:

```bash
cd archground-dist
```

Finally, we need to install the library. Given that this is a highly dynamic package were new additions are frequently introduced, we recommend installing the package but keeping it editable by typing:

```bash
pip install --editable .
```

This will install **archground**, but with the ability to update any changes made to the code. In addition, it will install all the required libraries it depends on.

### 2. Installing the Mars Climate Database and Venus Climate Database

Most of the atmospheric and surface models created for Mars and Venus rely on the extraction of atmospheric vertical profiles from the [Mars Climate Database](https://www-mars.lmd.jussieu.fr/mars/access.html) and the [Venus Climate Database](), respectively.

To enable this functionality, the users must download these databases and compile the Fortran functions that enable the extraction of the profiles from the files. Here, we are going to show the steps of how to set up the MCD, but the steps for the VCD should be analogous.

The full version of the MCD can be downloaded from this [link](https://www-mars.lmd.jussieu.fr/MCD_pro/mcd_pro.html).

Once downloaded, we need to compile the fortran functions of the MCD using the *archground-dist/archground/gcm/fmcd_gfortran.sh* compilation file stored within the acs-dist distribution. Before compiling it, we must modify with following fields:

- NETCDF: Path to the directory where the NETCDF libraries are compiled.
- wheremcd: Path to the Mars Climate Database.
- version: Version of the MCD that will be used.

Once these fields are modified according to your local settings, we can compile it by writing in the terminal:

```bash
source fmcd_gfortran.sh
```

This will create a file starting with *fmcd* and with a *.so* extension. This will be automatically loaded when importing the **archground** library into our Python programs.


### 3. Downloading the archNEMESIS spectroscopic databases

Some of the radiative transfer calculations performed with **archground** rely on spectroscopic parameters from different databases, such as HITRAN, GEISA, HITEMP, etc. These databases need to be stored in a format readable by archNEMESIS.

Reference spectroscopic databases compatible with archNEMESIS are available from our dedicated collection on DIGITAL.CSIC:

<a href="https://digital.csic.es/handle/10261/435473" target="_blank">
    <img src="docs/_static/digital_csic_logo.png"
         alt="DIGITAL.CSIC"
         width="250">
</a>

Users wishing to work with spectroscopic databases that are not included in this collection can generate their own archNEMESIS-formatted databases. If you require assistance with this process, please contact us.



### 4. Updating relevant paths to match local settings

The **archground** library relies on using data from other tools that cannot be unfortunately stored on the same GitHub repository, like the atmospheric and spectroscopic databases.

Once we have downloaded them. We need to update the paths in the library according to our local settings. The file we need to change is *archground-dist/archground/paths.py*.
