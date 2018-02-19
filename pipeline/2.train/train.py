#cdarby@jhu.edu
#updated 01-24-18

#Python 2.7.13 (tested 01-24-18)
#Python 3.6.2 (tested 01-24-18)

from __future__ import print_function


from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import random,argparse
from sklearn.metrics import confusion_matrix

parser = argparse.ArgumentParser(description='python train.py')
parser.add_argument('--out', help='Output classifier name',required=False,default="clf.pkl")
parser.add_argument('--mosaic', help='Feature file for mosaic sites',required=True)
parser.add_argument('--germline', help='Feature file for germline sites' ,required=True)
parser.add_argument('--mindepth', help='Minimum depth to use as training example' ,required=False,default=16)
parser.add_argument('--nestimators', help='Number of estimators in RFC' ,required=False,default=100)
parser.add_argument('--maxleafnodes', help='Maximum leaf nodes per tree in RFC' ,required=False,default=50)

args = parser.parse_args()

DPTH=args.mindepth
NEST=args.nestimators
MXNODES=args.maxleafnodes

#0th item in feature vecot is depth.
mosaicsites = [line.strip().split() for line in open(args.mosaic, 'r') if float(line.strip().split()[0]) > DPTH]
germlinesites = [line.strip().split() for line in open(args.germline, 'r') if float(line.strip().split()[0]) > DPTH]

#Make classes equal size, sample randomly
if len(mosaicsites) > len(germlinesites):
	mosaicsites = random.sample(mosaicsites,len(germlinesites)) 
elif len(germlinesites) > len(mosaicsites):
	germlinesites = random.sample(germlinesites,len(mosaicsites))

clf = RandomForestClassifier(n_estimators=NEST,max_leaf_nodes=MXNODES)
X = mosaicsites+germlinesites
Y = [1 for i in range(len(mosaicsites))] + [0 for i in range(len(germlinesites))]
clf = clf.fit(X,Y)
print("Confusion Matrix (class 0=Germline, class 1=Mosaic)")
#cols are predicted, rows are actual
print(confusion_matrix(Y,clf.predict(X)))

#print zip(clf.feature_importances_,descriptions0)		
joblib.dump(clf, args.out) 

