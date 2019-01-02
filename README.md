# PhageMiner
PhageMiner is a user-supervised semi-automated bioinformatics tool for prophage identification in bacterial genomes

Author: Reza Rezaei Javan

E-mail: reza.rezaeijavan@ndm.ox.ac.uk

Copyright (C) 2019 University of Oxford

## General remarks
This script can run on both Linux and macOS. It has been tested on the following operating systems:
Linux Ubuntu – 14.04.5 LTS
macOS Mojave – 10.14
macOS High Sierra – 10.13

## Usage
To use PhageMiner, put the main script (PhageMiner.py) and the annotated GenBank file of the bacterial genome (*.gb) in the same folder. 
Run the software using the following command:
```
python PhageMiner.py <input file>
```

## Installation
PhageMiner has the following dependencies:

### Required dependencies
* biopython
* pandas
* sklearn
* reportlab

There are a number of ways to install these dependencies and details are provided below for those not familiar with command line interfaces.

### biopython
To install biopython, insert the following command in the terminal:
```
pip install biopython
```
### pandas
To install pandas, insert the following command in the terminal:
```
pip install pandas
```
### sklearn
To install sklearn, insert the following command in the terminal:
```
pip install -U scikit-learn
```
### reportlab
To install reportlab, insert the following command in the terminal:
```
pip install reportlab
```
