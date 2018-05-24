#cdarby@jhu.edu
#updated 05-23-18

#Python 2.7.10 / scikit-learn 0.16.1 (tested 05-24-18)
#Python 3.4.2 / scikit-learn 0.19.1 (tested 05-24-18)

from __future__ import print_function

import argparse,sys

if sys.version_info > (3,0):
	import _pickle as cPickle
else:
	import cPickle

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

NFEAT = 33
MIN_CLF_SCORE = 0.9

parser = argparse.ArgumentParser(description='python classify.py --clf outs/clf.pkl --vectors outs/vectors.txt')
parser.add_argument('--clf',help="random forest scikit-learn classifier",required=True)
parser.add_argument('--vectors',help="text file of lists output by filter",required=True)

args = parser.parse_args()
clf = joblib.load(args.clf)

with open(args.vectors) as F:
	lines = F.readlines()

predictionArray = np.empty([len(lines),NFEAT])
for (i,L) in enumerate(lines):
	predictionArray[i] = np.array(eval(L)[3:])#field 0 is chrom; field 1 is pos; field 2 is minor allele base; rest are features
	
p = clf.predict_proba(predictionArray)
for i in range(len(p)):
	if p[i][1] >= MIN_CLF_SCORE:
		L = eval(lines[i])
		print(str(L[0]) + "\t" + str(L[1]) + "\t" + str(L[1]+1) + "\t" + str(p[i][1]) + "\t" + str(L[3]) + "\t" + str(L[8]) + "\t" + str(L[6]) + "\t" + L[2])	#chrom, pos, pos+1, score, depth, nC, MAF