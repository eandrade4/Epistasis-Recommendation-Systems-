# -*- coding: utf-8 -*-
"""S2023 - CLUSTER AND NERDS FUNCTIONS ON LOGBWD70.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1BdRVsZUD63Kk1z41aqKyqBb5e6peIvLl

The purpose of this program is to use the generalized code to analyze the Karst Mice data kmeans clustering and NERDS. 

The recommendations from our adjusted counting method are found before and after the implementation of NERDS so the outputs can be compared. 

To begin, packages are imported, functions from S2023_OnlyFunctions are imported, and google drive is mounted. Then, the following steps are completed: 


*   A directory is set up in Google drive for figures and outputs to be saved to. 

*   The genotype/phenotype data is imported, filtered, and standardized. 

*   An optimization algorithm is used to determine how many clusters should be implemented within the k-means clustering algorithm. 

*   Clusters are created, and the cluster with the highest or lowest average phenotype is extracted. This can be controlled by the user. 

*   First through fourth order coalitions of SNPs are recommended based on an adjusted counting method. 

*   The data is reduced using SVD and low rank approximation keeping 95% of the variance. 

*   NERDS (Neighborhood Epistasis Recommendation Detection with varying Similarity) is implemented on the final cluster, resulting in a list of recommended first through fourth order coalitions.
"""

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

# Mount your google drive in google colab
from google.colab import drive
drive.mount('/content/drive')

import sys
sys.path.insert(0,'/content/drive/My Drive/Epistasis Project')

#import S2023_OnlyFunctions as ep
from S2023_OnlyFunctions import *

images_dir = drive_setup("/content/drive/MyDrive/Epistasis Project/Emmas Code Outputs")

phenotype, dfphen = import_clean("LogBWD70", np.append([4],np.arange(6,20)),'/content/drive/MyDrive/Epistasis Project/K_final_females.csv' )
df, features, scaled_features = standardize_features(dfphen,phenotype)

k = choose_num_clusters(df, scaled_features, images_dir)

kmeans = cluster_features(k,scaled_features)

cluster_viz(k,"Clustering results from Karst Data with 2 Components", features, images_dir)

clusters = df_clusters(k,dfphen,kmeans);

finaldf, n = choose_finaldf("max", clusters, phenotype)
finaldf

max_data_fun(finaldf, finaldf.index, phenotype)

occurrences_viz(finaldf, images_dir, phenotype, n)

clustervsavg_viz(dfphen, finaldf, phenotype, images_dir, n)

NoNerdsRecs = recommend1to4order(dfphen, finaldf, 1, 2, 4, 6, dfphen, phenotype)
NoNerdsRecs

reduceddf, reducedarray = filter_for_NERDS(finaldf, dfphen, phenotype)

U,S,V = normalize_data(reduceddf, reducedarray)

c = solve_for_c(.95, S, reduceddf)

person_id, sliced, normalsliced = similarity_setup("max", reduceddf, phenotype, V, c)

cosindex = top_cosine_similarity(sliced, person_id, c)

print_similar_people(reducedarray, person_id, cosindex)

jacindex = Jaccard(normalsliced,person_id,c)

duped = metric_overlap(cosindex,jacindex, reduceddf, images_dir)

NerdsRecs = recommend1to4order(reduceddf, duped, 1, 2, 4, 6, dfphen, phenotype)
NerdsRecs

NerdsRecs = recommend1to4order(dfphen, duped, 1, 2, 4, 6, dfphen, phenotype)
NerdsRecs

dfphen

"""#Creating a Figure of [2,10,11]"""

dfnew = dfphen.loc[(dfphen['2'] != 0) & (dfphen['10'] != 0) & (dfphen['11'] != 0)]
dfnew

fig = go.Figure()

# add a scatterplot (Individuals vs Phenotype)
fig.add_trace(go.Scatter(x=dfphen.index, y=dfphen.loc[:,phenotype], mode='markers',name ='Individuals',  marker=dict(
            color='steelblue')))
fig.add_trace(go.Scatter(x=dfnew.index, y=dfnew.loc[:,phenotype], mode='markers', name=f'Individuals with [2,10,11]', marker=dict(
            color='coral')))
fig.add_trace(go.Scatter(x = [dfphen.index.min(), dfphen.index.max()],
                           y = [dfphen[phenotype].mean(), dfphen[phenotype].mean()],
                           mode = "lines", name='Average of all Individuals', line=dict(color='steelblue')))
fig.add_trace(go.Scatter(x = [dfphen.index.min(), dfphen.index.max()],
                           y = [dfnew.loc[:,phenotype].mean(), dfnew.loc[:,phenotype].mean()],
                           mode = "lines", name=f'Average of Individuals with [2,10,11]', line=dict(color='coral')))

  # add labels and title
fig.update_layout(title=f'Individuals with [2,10,11] vs Average', xaxis_title="Individuals", yaxis_title=f"Phenotype ({phenotype})")

fig.update_layout(legend=dict(
    yanchor="top",
    y=0.99,
    xanchor="left",
    x=0.01
))

fig.show()