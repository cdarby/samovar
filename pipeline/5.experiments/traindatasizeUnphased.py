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
parser.add_argument('--mosaic', help='Feature file for mosaic sites',required=True)
parser.add_argument('--het', help='Feature file for germline sites (het)' ,required=True)
parser.add_argument('--hom', help='Feature file for germline sites (hom)' ,required=True)
parser.add_argument('--mindepth', help='Minimum depth to use as training example' ,required=False,default=16)
parser.add_argument('--nestimators', help='Number of estimators in RFC' ,required=False,default=100)
parser.add_argument('--maxleafnodes', help='Maximum leaf nodes per tree in RFC' ,required=False,default=50)
parser.add_argument('--outpfx', help='outfile prefix' ,required=False,default="")


args = parser.parse_args()

DPTH=int(args.mindepth)
NEST=int(args.nestimators)
MXNODES=int(args.maxleafnodes)

#0th item in feature vect is depth.
hetsites = [line.strip().split() for line in open(args.het, "r") if float(line.strip().split()[0]) > DPTH]
hetsites = [ [line[0]] + [line[3]] + line[22:26] + line[27:32] for line in hetsites] #Remove all features regarding phasing

homsites = [line.strip().split() for line in open(args.hom, "r") if float(line.strip().split()[0]) > DPTH]
homsites = [ [line[0]] + [line[3]] + line[22:26] + line[27:32] for line in homsites] #Remove all features regarding phasing

test_het = hetsites[:5000]
test_hom = homsites[:5000]

hetsites = hetsites[5000:]
homsites = homsites[5000:]

#Make germline class equally het/hom, sample randomly
if len(hetsites) > len(homsites):
	hetsites = random.sample(hetsites,len(homsites)) 
elif len(homsites) > len(hetsites):
	homsites = random.sample(homsites,len(hetsites))

germlinesites = hetsites + homsites
test_germline = test_het + test_hom

hetsites = None
homsites = None

mosaicsites = [line.strip().split() for line in open(args.mosaic, "r") if float(line.strip().split()[0]) > DPTH]
mosaicsites = [ [line[0]] + [line[3]] + line[22:26] + line[27:32] for line in mosaicsites] #Remove all features regarding phasing

test_mosaic = mosaicsites[:10000]
mosaicsites = mosaicsites[10000:]

#Make classes equal size, sample randomly
if len(mosaicsites) > len(germlinesites):
	mosaicsites = random.sample(mosaicsites,len(germlinesites)) 
elif len(germlinesites) > len(mosaicsites):
	germlinesites = random.sample(germlinesites,len(mosaicsites))

print(str(len(mosaicsites)))

clf = RandomForestClassifier(n_estimators=NEST,max_leaf_nodes=MXNODES)

for N in range(1000,min(21000,len(mosaicsites)),1000):
	X = random.sample(mosaicsites,N)+random.sample(germlinesites,N)
	Y = [1 for i in range(N)] + [0 for i in range(N)]
	clf = clf.fit(X,Y)
	print(N)
	p_m = clf.predict_proba(test_mosaic)
	p_g = clf.predict_proba(test_germline)
	with open(args.outpfx + str(N)+".txt","w") as F:
		for i in range(10000):
			F.write("1\t" + str(round(p_m[i][1],5)) + "\n")
			F.write("0\t" + str(round(p_g[i][1],5)) + "\n")
