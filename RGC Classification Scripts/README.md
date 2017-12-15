# RGC-Classification
Including Scripts to Classify Cells

To run the entire process, open Matlab and run the script RunMe.m

The script will create the train dataset from Sumbul Data, and a test dataset from our data.

The function CreateTrainDataset creates the Train Data

The function CreateTestDataset creates the Test Data

The parameters are chosen to be the 5 peaks and the widths of the 5 peaks, along with the area under curve of stratification profile.

After the data is trained, two labeling schemes from Random Forest and Extra Trees Classifier are created and cross-validated, to give the final result.

The result is then plotted for inspection 
