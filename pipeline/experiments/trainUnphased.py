#cdarby@jhu.edu
#updated 01-24-18

#Python 2.7.13 (tested 01-24-18)
#Python 3.6.2 (tested 01-24-18)

from __future__ import print_function


from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import random,argparse
from sklearn.metrics import confusion_matrix

FEATURENAMES = ["depth","MAF","Mavgbq","Mavgpos","Mavgclip",
        "Mavgind","weightedMbq","Javgbq","Javgpos",
        "Javgclip","Javgind"]

parser = argparse.ArgumentParser(description='python trainUnphased.py --mosaic mosaic.features.tsv --het het.features.tsv --hom hom.features.tsv')
parser.add_argument('--out', help='Output classifier name',required=False,default="clf.pkl")
parser.add_argument('--mosaic', help='Feature file for mosaic sites',required=True)
parser.add_argument('--het', help='Feature file for het sites' ,required=True)
parser.add_argument('--hom', help='Feature file for hom sites' ,required=True)
parser.add_argument('--mindepth', help='Minimum depth to use as training example' ,required=False,default=16)
parser.add_argument('--nestimators', help='Number of estimators in RFC' ,required=False,default=100)
parser.add_argument('--maxleafnodes', help='Maximum leaf nodes per tree in RFC' ,required=False,default=50)
parser.add_argument('--modelsize', help='Number of het + hom; mosaic to use for model' ,required=False,default=None)

args = parser.parse_args()

DPTH=args.mindepth
NEST=args.nestimators
MXNODES=args.maxleafnodes
if args.modelsize is not None: modelsize = int(args.modelsize)

#0th item in feature vect is depth.
mosaicsites = [line.strip().split() for line in open(args.mosaic, 'r') if float(line.strip().split()[0]) > DPTH]

mosaicsites = [ [line[0]] + [line[3]] + line[22:26] + line[27:32] for line in mosaicsites] #Remove all features regarding phasing
germlinesites = [line.strip().split() for line in open(args.het, 'r') if float(line.strip().split()[0]) > DPTH] + [line.strip().split() for line in open(args.hom, 'r') if float(line.strip().split()[0]) > DPTH]
germlinesites = [ [line[0]] + [line[3]] + line[22:26] + line[27:32] for line in germlinesites] #Remove all features regarding phasing

if args.modelsize is not None and len(germlinesites) >= modelsize and len(mosaicsites) >= modelsize:
	mosaicsites = random.sample(mosaicsites,modelsize)
	germlinesites = random.sample(germlinesites,modelsize)
else: #Make classes equal size, sample randomly
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
print("Feature importances")
print(sorted(list(zip(clf.feature_importances_,FEATURENAMES))))	
joblib.dump(clf, args.out) 

