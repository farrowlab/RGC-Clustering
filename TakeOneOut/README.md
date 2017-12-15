# RGC-Cluster-Evaluation
# Overview
Included here are scripts used to evaluate and determine the number of clusters in a dataset. 

The first script, callClusterSumbulData.m, runs on Sumbul's data, and returns the Cluster Look up Table, Rand Index and Similarity Index for each values in the tables. 

The second script, callClusterOurData.m, runs on our data, and also returns the Cluster Look up Table, Rand Index and Similarity Index for each values in the tables.

# Details
Implemented here is the take-one-out method, where each cell is taken out, and the remaining cells are reclustered. The taken out cell is then approach with the new clustering scheme, and we calculated the similarity index between the new cluster the cell fits into and the old cluster
