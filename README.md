<p align="center">
  <a href="https://archnemesis.readthedocs.io/en/latest/">
    <img src="https://raw.githubusercontent.com/juanaldayparejo/archnemesis-dist/main/docs/images/archnemesis_logo_black_background.png"
         alt="archNEMESIS logo"
         width="400">
  </a>
</p>

<h1 align="center">archground</h1>

<p align="center">
  <strong>Analysis of planetary spectra from ground-based observatories</strong>
</p>

**archground** extends [archNEMESIS](https://archnemesis.readthedocs.io/en/latest/)
with tools tailored to ground-based observations of planetary atmospheres.
It brings together the archNEMESIS radiative transfer and retrieval algorithms,
instrument-specific functionality, and reference atmospheric models.

This package is maintained by [Juan Alday](https://gapt.iaa.es/content/juan-alday).

Questions, suggestions, and contributions are welcome. To report a problem,
request documentation, or suggest a Jupyter notebook example, please
[open an issue](https://github.com/juanaldayparejo/archground-dist/issues)
or contact the maintainer.


## Supported instruments

**archground** supports the analysis of planetary spectra obtained with the following instruments:

<table>
  <tr>
    <th align="center">NASA Infrared Telescope Facility</th>
    <th align="center">Calar Alto Observatory</th>
    <th align="center">SOFIA Stratospheric Observatory</th>
  </tr>
  <tr>
    <td align="center">
      <img src="https://irsa.ipac.caltech.edu/images/icons/missions/irtf_logo.png"
           alt="NASA Infrared Telescope Facility logo"
           height="90">
    </td>
    <td align="center">
      <img src="https://www.caha.es/APPS/CfP/public/images/CAHA-Logo.gif"
           alt="Calar Alto Observatory logo"
           height="90">
    </td>
    <td align="center">
      <img src="https://upload.wikimedia.org/wikipedia/commons/8/83/SOFIA_mission_patch.png?utm_source=en.wikipedia.org&utm_campaign=index&utm_content=original"
           alt="SOFIA logo"
           height="90">
    </td>
  </tr>
  <tr>
    <td align="center">
      <a href="https://irtfweb.ifa.hawaii.edu/~texes/Main%20Website/NewTEXES.html"><strong>TEXES</strong></a>
    </td>
    <td align="center">
      <a href="https://carmenes.caha.es/"><strong>CARMENES</strong></a>
    </td>
    <td align="center">
      <a href="https://irsa.ipac.caltech.edu/data/SOFIA/docs/instruments/exes/"><strong>EXES</strong></a>
    </td>
  </tr>
</table>




## Supported planets

**archground** currently includes atmospheric models for **Mars** and **Venus**.
Reference atmospheric profiles are extracted from the corresponding climate databases:

| Planet | Reference atmospheric profiles |
| :--- | :--- |
| 🔴 **Mars** | [Mars Climate Database (MCD)](https://www-mars.lmd.jussieu.fr/mars/access.html) |
| 🟡 **Venus** | [Venus Climate Database (VCD)](http://www-venus.lmd.jussieu.fr/) |

See [Installation](#installation) for instructions on downloading the climate
databases, compiling their Fortran interfaces, and configuring local paths.





## Installation

To install this package, please follow the next steps:

### 0. Install archNEMESIS

This python package is essentially an adaptation of the archNEMESIS python package prepared for the analysis of ground-based telescope observations of planetary atmospheres.

As such, it is expected that users will already have a working version of archNEMESIS in their virtual environments. If not, please visit the archNEMESIS [website](https://archnemesis.readthedocs.io/en/latest/).


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


## Citing archNEMESIS

If archground or archNEMESIS have been significant in your research, we suggest citing the following articles:

- archNEMESIS reference publication:
   - Alday, J., Penn, J., Irwin, P., Mason, J., Yang, J. and Dobinson, J. (2025) archNEMESIS: An Open-Source Python Package for Analysis of Planetary Atmospheric Spectra, *Journal of Open Research Software*, 13(1), p. 10. doi: [10.5334/jors.554](https://doi.org/10.5334/jors.554).

- NEMESIS reference publication:
   - Irwin, P. G. J., Teanby, N. A., De Kok, R., Fletcher, L. N., Howett, C. J. A., Tsang, C. C. C., ... & Parrish, P. D. (2008). The NEMESIS planetary atmosphere radiative transfer and retrieval tool. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 109(6), 1136-1150. doi: [10.1016/j.jqsrt.2007.11.006](https://doi.org/10.1016/j.jqsrt.2007.11.006).
