# RGC-Classification
Including Scripts to Classify Cells

To run the entire process, open Matlab and run the script RunMe.m

The script will create the train dataset from Sumbul Data, and a test dataset from our data.

The function CreateTrainDataset creates the Train Data

The function callClassifyOurData creates the Test Data

The parameters are chosen to be the 5 peaks and the widths of the 5 peaks, along with the area under curve of stratification profile.

Only 298 cells from our data are used. 

Extra Decision Trees Classifier with SMOTE Oversampling is used on the train data, and provide classification result from our data.

The result is then plotted for inspection 
