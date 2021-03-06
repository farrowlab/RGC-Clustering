import numpy as np
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.model_selection import cross_val_score
from sklearn import linear_model
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
#from sklearn.ensemble import GradientBoostingClassifier
#from sklearn.neural_network import MLPClassifier
from sklearn.metrics.pairwise import chi2_kernel
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from glmnet import LogitNet
from sklearn.model_selection import StratifiedKFold
import pickle
f = open ('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Dataset/traindata.txt','r').readlines()

N = len(f) 
l = [];
for i in range(0,N):
	w = f[i].split(',')
	try:
		list1 = [float(x) for x in w]
		l.append(list1)
	except ValueError,e:
		print "error",e,"on line", i
data = np.array(l)
print data.shape

f2 = open('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Dataset/OurDataTest.txt','r').readlines()
N = len(f2) 
l = [];
for i in range(0,N):
	w = f2[i].split(',')
	try:
        	list1 = [float(x) for x in w]
        	l.append(list1)
	except ValueError,e:
        	print "error",e,"on line", i
testdata = np.array(l)
print testdata.shape
a,b = testdata.shape

X = data[:,0:b]
Y = data[:,b].astype(int)
from sklearn.model_selection import ShuffleSplit

#cv = ShuffleSplit(n_splits=10, test_size=0.25, random_state=0)
cv = StratifiedKFold(n_splits=10, random_state=None, shuffle=True)
print "NOw with SMOTE oversampling"
##########################NOw do resampling
from imblearn.over_sampling import SMOTE,ADASYN
X_resampled, Y_resampled = SMOTE(random_state=42).fit_sample(X,Y)
X = X_resampled
Y = Y_resampled

print "Extra Trees Classifier"
model = ExtraTreesClassifier(n_estimators=200, max_depth = None, min_samples_split=2, random_state=0)
#model3 = model3.fit(X,Y)
#model = make_pipeline(preprocessing.StandardScaler(), model)
score3 = cross_val_score(model,X,Y, cv = cv)
print score3.mean()
model = model.fit(X,Y)

classes = list()
problist = list()
#pfile = open('prob.txt','w')
for i in range(0,a):
	testsubject = testdata[i,:]
	testsubject = testsubject.reshape(1,-1)
	c = model.predict(testsubject)
#	print c 
	prob = model.predict_proba(testsubject)
	#maxprob = max(prob)	
	maxi = 0
	for item in prob:
	#	if item >= maxi:
		problist.append(item)		
	#		maxi = item
#		for stuff in item:
#			if stuff >= maxi:
#				maxi = stuff
      #	print maxi
      # for item in prob:
      #         pfile.write("%s\t" %str(item))
      # pfile.write("\n")       	       
#	classes.append((c,maxi))
	classes.append(c)
#	problist.append(maxi)	

from sklearn.externals import joblib
joblib.dump(model, '/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Results/ExtraTree.pkl')

tfile = open('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Results/classesLABEL.txt','w')
pfile = open('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Results/classesPROB.txt','w')
for item in classes:
	for stuff in item:	
		tfile.write("%s\n" % str(stuff)) 
#	tfile.write("\n")

for item in problist:
	for stuff in item:
		pfile.write("%s " % str(stuff))
	pfile.write("\n")


