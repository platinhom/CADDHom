#! /usr/bin/env python
import fr
import numpy as np

neednum=20
weighs=[2,10,10,15,10];

sweighs=sum(weighs)
weighreals=[float(i)/sweighs for i in weighs]

fr=open(sys.argv[1])
datas=[]
for line in fr:
	tmp=line.strip().split()
	datas.append(tmp)
numfeature=0
numdata=0
if (len(datas)>0):
	numdata=len(datas)
	numfeature=len(datas[0])

features=[ [] for i in range(numfeature)]
for i in range(numdata):
	for j in range(numfeature):
		features[j].append(datas[i][j])
maxfeatures=map(max,features)


