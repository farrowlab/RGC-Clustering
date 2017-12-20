# RGC-Classification
# Overview
Included here are Scripts to Classify Cells based on the 15 cell types from Sumbul

# Step by Step Instruction

To run the entire process, open Matlab and run the script RunMe.m

# Details 
The script will create the train dataset from Sumbul Data, and a test dataset from our data.

The function CreateTrainDataset creates the Train Data.

The function callClassifyOurData creates the Test Data.

The parameters are chosen to be Location of the 5 peaks and the widths of the 5 peaks obtained from matlab findpeak function, along with the area under curve of the  stratification profile of individual cell.

Our data are automatically pulled from the Database file, so any changes there will change the number of cells being classified. 

Extra Decision Trees Classifier with SMOTE Oversampling is used on the train data, and provide classification result from our data.

This classifier is done with a python script called EDTClassifier.py

Finally, the result is then plotted for inspection, and this is inside the callCheckClassification function

