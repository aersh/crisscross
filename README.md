# Robust nucleation control via crisscross polymerization of highly coordinated DNA slats

Dionis Minev, Christopher M. Wintersinger, Anastasia Ershova, William M. Shih

License: MIT License

## Stochastic model simulations

### System requirements

The stochastic model simulation code is provided as Jupyter notebooks running Python 3.

The following versions of the dependencies were used:

* jupyter notebook 6.0.1
* matplotlib 3.1.1
* numpy 1.17.2
* pandas 0.25.1
* seaborn 0.9.0
* scipy 1.3.1

This code was tested on MacOS Mojave and Catalina, but should readily work on other operating systems that support Jupyter notebooks.

### Installation guide

We recommend using Anaconda to install all the required dependencies: https://docs.anaconda.com/anaconda/install/

Alternatively all the packages above are pip-installable `pip install matplotlib numpy pandas seaborn scipy`

### Demo and instructions for use

The stochastic model Jupyter notebooks are provided as stand-alone scripts and can be simply re-run to reproduce the simulation results.

## Sequence design

### System requirements

Sequence design code is provided as individual Python scripts for each design version.

The following versions of the dependencies were used:

* python 2.7.16
* unafold 3.8

This code was tested on MacOS Mojave and Catalina.

### Installation guide

The UNAFold nucleic acid folding binaries are no longer publicly available so running these scripts relies on UNAFold already being installed on the user's computer.

### Demo and instructions for use

The scripts that were used for designs in the manuscript are provided as is, and can be simply re-run using Python 2 to generate sequences for the same architectures with the same design parameters if UNAFold is installed.