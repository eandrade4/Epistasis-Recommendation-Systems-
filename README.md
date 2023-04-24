# Epistasis-Recommendation-Systems-

The Epistasis Recommendation Systems Program consists of a generalizable k-means clustering and recommendation system approach for detection of higher order coalitions of single nucleotide polymorphisms (SNPs) displaying epistasis. 

Single nucleotide polymorphisms (SNPs) are genetic mutations affecting a single base pair. 

The goal is to detect statistically significant coalitions of SNPs that display epistasis. Due to the exponential number of possible groupings of mutations, investigations coincide with challenges navigating the expansive computational search space. Thus, we present a user-item recommendation system complemented by a k-means clustering algorithm to effectively reduce the dimensionality of the search space for higher order interactions. More specifically, after implementing k-means clustering, we consider two similarity metrics, their intersection, and a more general framework to identify groups of mutations responsible for an increased quantitative phenotype with an adjusted count metric.

Future focus will be on the investigation of coupling our methodology with neural networks and deep learning approaches. We also plan on addressing the exponential number of possible SNP groupings using an adjusted k-means algorithm.

For a full description of the program, visit the
[project page](https://github.com/eandrade4/Epistasis-Recommendation-Systems-).

Submit bug reports and feature suggestions, or track changes in the
[issue queue](https://github.com/eandrade4/Epistasis-Recommendation-Systems-/issues).

## Table of Contents

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

```ruby
pip install kneed

pip install kaleido

import tarfile
import urllib
import time
from datetime import date
import os
from datetime import datetime


import kaleido
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
from itertools import combinations
import math
import collections
from matplotlib_venn import venn2

from kneed import KneeLocator, DataGenerator as dg
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler, StandardScaler
```

## Installation

Install the files titled s2023_cluster_and_nerds_functions_on_logbwd70.py, and s2023_onlyfunctions.py to your desired coding environment. These can be found on the [project page](https://github.com/eandrade4/Epistasis-Recommendation-Systems-).

## Framework
We begin with a clustering algorithm on genotype/phenotype data to reduce the number of individuals considered within our study. We represent this data as an m by n matrix were m is the number of individuals and n is the number of mutations present. The entries of this matrix are either 0, 1, or 2, indicating how many copies of the given mutation is present within the individual. 

This data is filtered using an optimized k-means++ algorithm. We move forward using the cluster with the highest or lowest average phenotype, depending on if the data set is being analyzed for positive or negative epistasis. Singular value decomposition and low-rank approximation are used to reduce the data even further (keeping 95 % of the variance). Next, the individual with the most extreme phenotype within the cluster is compared to the remaining individuals within the cluster by cosine and Jaccard similarity. The intersection between the most similar individuals from each method is taken, and the mean occurrences of each SNP coalition is calculated. This adjusted counting method results in a ranking of recommended SNP groupings that are predicted to display epistasis. 


## Example/Quick Start Guide
```ruby
# Mount your google drive in google colab
from google.colab import drive
drive.mount('/content/drive')

import sys
sys.path.insert(0,'/content/drive/My Drive/Epistasis Project') #This is the path to your folder containing your data set. 

#Import S2023_OnlyFunctions so we can access all of the functions. 
from S2023_OnlyFunctions import *

#This creates a folder called "Code Outputs" in your Google Drive for figures to be saved to. 
images_dir = drive_setup("/content/drive/MyDrive/Epistasis Project/Code Outputs")


#Make sure that in the csv or excel file you are going to import that the phenotype column comes first, and the genotypes are next. 
#Your column labels should be in increasing order before listing them here. 
phenotype, dfphen = import_clean("Name_of_phenotype_column_here_as_string", genotype_column_labels_go_here,'/content/drive/MyDrive/Epistasis Project/Name_of_dataset_here.csv' ) #.csv or .xlsx can be used

#Normalize the data.
df, features, scaled_features = standardize_features(dfphen,phenotype)

#This runs an optimization scheme using the elbow and silhouette method to choose the number of clusters best for your specific data set.
k = choose_num_clusters(df, scaled_features, images_dir)

#Create the clusters.
kmeans = cluster_features(k,scaled_features)

#Visualize clusters using 2 components.
cluster_viz(k,"Clustering results from Karst Data with 2 Components", features, images_dir)

#Creates a dataframe for each of the clusters, and names them cluster0,...,clusteri for i many clusters.
clusters = df_clusters(k,dfphen,kmeans);

#Choose the cluster you will continue with. If analyzing for positive epistasis, the first parameter is "max" else, use "min" for negative epistasis.
finaldf, n = choose_finaldf("max_or_min_here_as_str", clusters, phenotype)

#The dataframe called finaldf is the chosen cluster.

#Count the occurrences of SNPs in finaldf.
max_data_fun(finaldf, finaldf.index, phenotype)

#Create a barchart of the occurrences of each SNP.
occurrences_viz(finaldf, images_dir, phenotype, n)

#Create a scatterplot with the individuals in the cluster compared to the rest of the individuals in the study with respect to their phenotype.
clustervsavg_viz(dfphen, finaldf, phenotype, images_dir, n)

#We see what recommendations we get for first through fourth order SNPs using the counting method on the clustered individuals *before* doing our dimension reduction and similarity metrics. This allows us to compare results later. 
NoNerdsRecs = recommend1to4order(dfphen, finaldf, number_of_first_order_recs_wanted_here, number_of_second_order_recs_wanted_here, number_of_third_order_recs_wanted_here, number_of_fourth_order_recs_wanted_here, dfphen, phenotype)

#NoNerdsRecs is the dataframe containing the top recommended SNP coalitions.

#Filter data to prepare for Neighborhood Epistasis Recommendation Detection with varying Similarity (NERDS) algorithm.
reduceddf, reducedarray = filter_for_NERDS(finaldf, dfphen, phenotype)

#Singular value decomposition is used here.
U,S,V = normalize_data(reduceddf, reducedarray)

#Determine number of principal components to use in order to keep 95 percent variance. 
c = solve_for_c(.95, S, reduceddf)

#Create the sliced data with c principal components.
person_id, sliced, normalsliced = similarity_setup("max_or_min_here_as_str", reduceddf, phenotype, V, c)

#Rank individuals based on cosine similarity.
cosindex = top_cosine_similarity(sliced, person_id, c)

#This function can print the highest ranked individuals.
print_similar_people(reducedarray, person_id, cosindex)

#Rank individuals based on Jaccard similarity.
jacindex = Jaccard(normalsliced,person_id,c)

#Find the intersection of the top cosine and Jaccard similar individuals.
duped = metric_overlap(cosindex,jacindex, reduceddf, images_dir)

#Use counting method to recommend first through fourth order groupings of SNPs with respect to reduceddf. Note that in this function, the parameter reduceddf can be changed to dfphen in order to calculate the mean occurrences of SNP groupings with respect to all the individuals in the study, rather than only the ones in the chosen cluster.
NerdsRecs = recommend1to4order(reduceddf, duped, number_of_first_order_recs_wanted_here, number_of_second_order_recs_wanted_here, number_of_third_order_recs_wanted_here, number_of_fourth_order_recs_wanted_here, dfphen, phenotype)

#NerdsRecs consists of the top recommended first through fourth order SNP groupings predicted to display epistasis.

```


## Acknowledgements 
This material is based upon research done with Dr. Mario Banuelos and Nicholas Tom as part of the MBG Research group. This work was also supported by CSU Fresno Undergraduate Research Grant Program, and CSU Fresno Department of Mathematics.

More information about this project is discussed in our paper titled, [User-Item Recommendation Approaches to Detect Genomic Variant Interactions](https://ieeexplore.ieee.org/abstract/document/9980252?casa_token=XbdrxWBJBmoAAAAA:APklFDdywKot5zCRE8KnQq2iagfVyz7QnnF5GPvYVfdH47ZsccDkWGtinTdLtYvvCX3IGnr8). 

- E. Andrade, N. Tom and M. Banuelos, "User-Item Recommendation Approaches to Detect Genomic Variant Interactions," 2022 Asia-Pacific Signal and Information Processing Association Annual Summit and Conference (APSIPA ASC), Chiang Mai, Thailand, 2022, pp. 1457-1461, doi: 10.23919/APSIPAASC55919.2022.9980252.


