%%%This is a script to create the train and test dataset
%%%Create the classfier from the train dataset, and classify the test dataset%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Create Train Dataset from Sumbul's data
createTrainDataset
%2) Create Test Dataset from Our Data
createTestDataset
%%%%%%%%%%traindata.txt and testdata.txt will be created%%%%%%%%
%3) Run Python Script to read train data and create Random Forest Classifier
system('python RandomForestSMOTE.py'); %%detect ON
%4) Run Python Script to read train data and create Extra Trees Classifier
system('python ExtraDecisionTreeSMOTE.py');
%5) Combine two Classifying Schemes to get final result
combineResult
