# RGC-Cluster-Evaluation
# Overview
- Included here are scripts used to evaluate and determine the number of clusters in a dataset. 

- The first script, callClusterSumbulData.m, runs on Sumbul's data, and returns the Cluster Look up Table, Rand Index and Similarity Index for each values in the tables. 

- The second script, callClusterOurData.m, runs on our data, and also returns the Cluster Look up Table, Rand Index and Similarity Index for each values in the tables.

# Basic Usage
- Run callClusterOurData.m to run the Take One Out method that will validate the number of clusters based on the chosen Clustering Scheme.
- The default clustering scheme is e-linkage.
- To change to a different clustering scheme, open clusterRGCarborDensitiesQD.m and make the necessary adjustment.
 

# Details
- Implemented here is the take-one-out method, where each cell is taken out, and the remaining cells are reclustered. 

- The taken out cell is then approach with the new clustering scheme, and we calculated the similarity index between the new cluster the cell fits into and the old cluster

- The resulting cluster and its similarity index are then combined and plotted for inspection in a histogram.

- User look at the bin with the highest vote and pick the cluster number accordingly.




