# RGC-Classification
# Overview
Included here are Scripts to Classify Cells based on the 15 cell types from Sumbul

# Step by Step Instruction

- To run the entire process, open Matlab and run the script RunMe.m

- The result will be stored in Results folder, with the labels of each cell contained in a txt file called classesLABEL.txt

- A trained model will also be there, and stored as ExtraTree.pkl

- There is also another txt file called classesPROB.txt, saying how confident the algo with its classification. 



# Details 

- The script will create the train dataset from Sumbul Data, and a test dataset from our data.

- The function CreateTrainDataset creates the Train Data.

- The function callClassifyOurData creates the Test Data.

- The parameters are chosen to be Location of the 5 peaks and the widths of the 5 peaks obtained from matlab findpeak function, along with the area under curve of the  stratification profile of individual cell.

- Our data are automatically pulled from the Database file, so any changes there will change the number of cells being classified. 

- This classifier is called with a python script called EDTClassifier.py

- Extra Trees Classifier with SMOTE Oversampling is used on the train data, and provide classification result from our data. SMOTE is used to balance the undersampled and oversampled classed in Sumbul Dataset. As for the classifier choice, Extra Trees give best training and testing result, but Random Forest Classifier also give good result, and can also be used instead. In fact, it's simple to switch Classifier model using sklearn. Inside the python code are several commented out option that can be used instead.

- Finally, the result is then plotted for inspection, and this is inside the callCheckClassification function

