# DNA_stability
- [Overview](#overview)
- [Documentation](#documentation)
- [System Requierments](#system-Requierments)
- [Installation Guide](#installation-guide)
- [Issues](https://github.com/Vishvaraj2601/DNA-Project/issues)

# Overview
Design of error-correcting codes for DNA-based data storage is based on an analysis of DNA stability.. 

# Documentation

As documentation for the softwarwe becomes available, it will be placed under the docs folder. We also generated one .exe file. It will generate the plots as a result direct in one click.

# System Requierments

## Hardware Specifications
Only a basic computer with appropriate RAM and processing capability is needed to run this code.

## Software Requierments
### OS Requirements
This package is supported for windows, macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3
+ Windows : Windows 8/10/11

Note that most OSes will support our software by using Docker.

### Library Dependences

For high precision floating point calculations, Python packages are dependent on libraries created in other programming languages. They are installable using popular package managers:
```
apt-get install libgmp-dev libmpc-dev libmpfr-dev
```
If you use Docker, the provided Dockerfile will install the needed libraries.

### Python Dependences

Python versions 3.6 to 3.8 were used for the development and testing of our programmes. It is dependent on the following things:

```
nose
sphinx
biopython
editdistance
statistics
matplotlib
numpy
scipy
gmpy2
```

# Installation Guide

## Use your local python environment
The easiest thing to do is download or checkout the code from GitHub if you already have Python 3 installed on your computer. Run the subsequent commands in the DNA stability directory after that:

    git clone https://github.com/Vishvaraj2601/DNA-Project
    cd DNA-Project
    pip install -r requirements.txt
    
## Execute Our Analysis

Run the main script this way:

    python3 stability.py
    
This command should run to completion within a few minutes.
