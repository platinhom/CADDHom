#! /usr/bin/env python
import sys,heapq
import numpy as np

neednum=20
#weighs=[20,10,10,15,10,5,6,7,8,9,1,2,3];
weighs=[ 1 for i in range(13)]

# Function defined here

def similarity(vals1,vals2,maxs,mins,weights,coeff=2):
	sv=0.0
	for i in range(len(vals1)):
		if ((maxs[i]-mins[i])==0): 
			#sv += weights[i];
			continue;
		sv+=weights[i]*pow(1-(abs(vals1[i]-vals2[i])/float(maxs[i]-mins[i])), coeff)
	return sv

# Read data
fr=open(sys.argv[1])
datas=[]
for line in fr:
	tmp=line.strip().split()
	datas.append(map(float,tmp))
fr.close()

# get the features and data number
numfeature=0
numdata=0
if (len(datas)>0):
	numdata=len(datas)
	numfeature=len(datas[0])

# Get the max/min of features
features=[ [] for i in range(numfeature)]
for i in range(numdata):
	for j in range(numfeature):
		features[j].append(datas[i][j])
maxfeatures=map(max,features)
minfeatures=map(min,features)

# Calculate the real weight of each item
# When max-min ==0 , weight=0, ignore it!
sumweigh=0.0
ignorefeatureNum=[]
weighreals=[]
for i in range(numfeature):
	if ((maxfeatures[i]-minfeatures[i])==0):
		ignorefeatureNum.append(i)
	else:
		sumweigh+=weighs[i]
for i in range(len(weighs)):
	if i in ignorefeatureNum:
		weighreals.append(0.0)
	else:
		weighreals.append(float(weighs[i])/sumweigh)

# Calculate similarity matrix
simmatrix=[ [0.0]*numdata for i in range(numdata)];
for i in range(numdata):
	simmatrix[i][i]=similarity(datas[i],datas[i],maxfeatures,minfeatures,weighreals,1)
	for j in range(i+1, numdata):
		sim=similarity(datas[i],datas[j],maxfeatures,minfeatures,weighreals,1)
		simmatrix[j][i]=simmatrix[i][j]=sim

# Write out similar matrix
formatlambda=(lambda val: format(val, '6.4f'));
fw=open("similarity.matrix",'w')
format
for i in range(numdata):
	fw.write(" ".join(map(formatlambda,(simmatrix[i])))+"\n")
fw.close()

# find out top similar id, start from 1
dlen='<'+str(len(str(numdata))); # < used for left alignment
formatlambda=(lambda val: format(val, dlen+'d'));
for i in range(numdata):
	largestvalue=heapq.nlargest(neednum+1, simmatrix[i])
	minlim=largestvalue[-1]
	outnums=[]
	for j in range(numdata):
		if (simmatrix[i][j] >= minlim and i != j): outnums.append(j+1)
	print " ".join(map(formatlambda,outnums))
