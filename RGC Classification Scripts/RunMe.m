%%%This is a script to create the train and test dataset
%%%Create the classfier from the train dataset, and classify the test dataset%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Create Train Dataset from Sumbul's data
createTrainDataset
%2) Create Test Dataset from Our Data
callClassifyOurDataQD
%%%%%%%%%%traindata.txt and testdata.txt will be created%%%%%%%%
%3) Run Python Script to read train data and create Extra Trees Classifier
system('python EDTClassifier.py'); %%detect 

%4) Check Results 
callCheckClassification

%system('python ExtraDecisionTreeSMOTE.py');
%combineResult
