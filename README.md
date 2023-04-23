# Epistasis-Recommendation-Systems-

The Epistasis Recommendation Systems Program consists of a generalizable k-means clustering and recommendation system approach for detection of higher order coalitions of single nucleotide polymorphisms displaying epistasis. 

The goal is to detect statistically significant coalitions of SNPs that display epistasis. Due to the exponential number of possible groupings of mutations, investigations coincide with challenges navigating the expansive computational search space. Thus, we present a user-item recommendation system complemented by a k-means clustering algorithm to effectively reduce the dimensionality of the search space for higher order interactions. More specifically, after implementing k-means clustering, we consider two similarity metrics, their intersection, and a more general framework to identify groups of mutations responsible for an increased quantitative phenotype with an adjusted count metric.

Future focus will be on the investigation of coupling our methodology with neural networks and deep learning approaches. We also plan on addressing the exponential number of possible SNP groupings using an adjusted k-means algorithm.

For a full description of the program, visit the
[project page](https://github.com/eandrade4/Epistasis-Recommendation-Systems-).

Submit bug reports and feature suggestions, or track changes in the
[issue queue](https://github.com/eandrade4/Epistasis-Recommendation-Systems-/issues).

## Table of contents

- Requirements
- Installation
- Framework
- Example/Quick Start Guide
- Acknowledgements

## Requirements

This program requires the following:

- Coding environment compatible with Python. Suggested: Jupyter Notebooks, Google Colaboratory
- Google Drive
- A chosen data set consisting of genotype/phenotype data for a set of individuals
- Several imported packages

<img width="591" alt="Screenshot 2023-04-22 at 6 28 53 PM" src="https://user-images.githubusercontent.com/130394178/233814571-1e895803-90a9-462d-81aa-511ab1f77306.png">

## Installation

Install the files titled s2023_cluster_and_nerds_functions_on_logbwd70.py, and s2023_onlyfunctions.py to your desired coding environment. These can be found on the [project page](https://github.com/eandrade4/Epistasis-Recommendation-Systems-).

## Framework
We begin with a clustering algorithm on genotype/phenotype data to reduce the number of individuals considered within our study. We represent this data as an m by n matrix were m is the number of individuals and n is the number of mutations present. The entries of this matrix are either 0, 1, or 2, indicating how many copies of the given mutation is present within the individual. 

**Not finished writing this section yet**


## Example/Quick Start Guide
**Not finished writing this section yet**


## Acknowledgements 
This material is based upon research done with Dr. Mario Banuelos and Nicholas Tom as part of the MBG Research group. This work was also supported by CSU Fresno Undergraduate Research Grant Program, and CSU Fresno Department of Mathematics.

More information about this project is discussed in our paper titled, [User-Item Recommendation Approaches to Detect Genomic Variant Interactions](https://ieeexplore.ieee.org/abstract/document/9980252?casa_token=XbdrxWBJBmoAAAAA:APklFDdywKot5zCRE8KnQq2iagfVyz7QnnF5GPvYVfdH47ZsccDkWGtinTdLtYvvCX3IGnr8). 

- E. Andrade, N. Tom and M. Banuelos, "User-Item Recommendation Approaches to Detect Genomic Variant Interactions," 2022 Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC), Chiang Mai, Thailand, 2022, pp. 1457-1461, doi: 10.23919/APSIPAASC55919.2022.9980252.


