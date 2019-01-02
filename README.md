# PhageMiner
## Introduction
**PhageMiner** is a user-supervised semi-automated bioinformatics tool for prophage identification in bacterial genomes.

The main objective of the PhageMiner software is to streamline the otherwise tedious manual curation process for prophage sequence discovery. 

PhageMiner uses a mean shift algorithm combined with annotation-based genome mining in
order to rapidly identify prophage sequences within complete or fragmented bacterial genomes that result from whole genome sequencing (genomes with multiple contigs). 

PhageMiner requires only minimal manual input from the user and produces various genome diagrams to aid a thorough discovery of previously unidentified prophages. 

PhageMiner can run on both Linux and macOS. 

**Author**: Reza Rezaei Javan

**E-mail**: reza.rezaeijavan@ndm.ox.ac.uk

Copyright (C) 2019 University of Oxford

## Usage
To use PhageMiner, put the main script (PhageMiner.py) and the input file (.gb) in the same folder. Next, run the software using the following command:
```
python PhageMiner.py <input file>
```

## Installation

### Required dependencies
PhageMiner has the following python dependencies:
* biopython
* pandas
* sklearn
* reportlab

There are a number of ways to install these dependencies and guidelines are provided below for those not familiar with command line interfaces.

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
### Operating Systems
PhageMiner can run on both **Linux** and **macOS**. It has been tested on the following operating systems:
* Linux Ubuntu – 14.04.5 LTS
* macOS Mojave – 10.14
* macOS High Sierra – 10.13

## Input file
PhageMiner only requires an annotated Genbank file as input. It can work with annotated GenBank (GBK) files automatically downloaded from the NCBI server. However, for the best result, it is highly recommended to freshly annotate your bacterial genome using the RasT server (http://rast.nmpdr.org/). PhageMiner can also work with files annotated using Prokka. An annotated *Streptococcus agalactiae* draft genome is provided as an example input file.

## Output files
PhageMiner produces GenBank files for each identified prophage region. It also produces genome diagrams of the prophage sequences within the bacterial genome. Based on user input during the discovery process, PhageMiner categorises identifed prophage regions into three groups: full-length prophages, satellite prophages and unknown phage-related regions. PhageMiner also produces a table containing the nucleotide sequences and annotations of all phage-related genes in the bacterial genome. Furthermore, a PDF diagram of the host genome is produced, in which the prophage regions and assembly gaps are highlighted.  Each of these features can be inactivated based on the user's preference.  

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/RezaRezaeiJavan/PhageMiner/issues) or email reza.rezaeijavan@ndm.ox.ac.uk.
