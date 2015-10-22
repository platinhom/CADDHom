#! /usr/bin/env python
# -*- coding: utf8 -*-

'''
File 1: y+x data for kplane
File 2: x data for training SVM (x here may not be same as above)
File 3: y+x test data for test in SVM 
'''

import os,sys
import HHKplane
import SVM_Data

if (len(sys.argv) < 4):
	print "Please give three input files: ";
	print __doc__;
	exit(1);

fk=sys.argv[1];
fs=sys.argv[2];
ft=sys.argv[3];

ngroup=2

### K-Plane

# create Kplane object
kp=HHKplane.KPlaneObj();
# Initialize data by file and set up program
# Data in file should be "y x1 x2 x3...."
kp.initializeByFile(fk,num_group=ngroup,regroup=False,randomgroup=False,print_in=False);

# Initialize data by file and set up program with initial group information.
# Data in file should be "y x1 x2 x3...."
#kp.initializeByFileWithGroup(fk,num_group=ngroup,regroup=False,randomgroup=False,print_in=False);

# Setup parameters. PS: initial will reset all default parameters, 
# you have to set them after initialization if you don't want to use default. 
# Default: klamda=0.001, rlamda=500,rmslimit=0.0001
kp.SetAllParameters(klamda=0.0001,rlamda=500, rmslimit=0.0001);
#Calculate the partitioning by k-plane regression
rms, coef, partition_groups = kp.kplane();
kp.outputGroupData(fk,OnlyGroup = True);
#print rms
#print coef
#print partition_groups

fklist=os.path.splitext(fk);
fgdata=fklist[0]+'_gdata'+fklist[1];
fpar=fklist[0]+'_par'+fklist[1];

### SVM_Data Preparation
outsets=1
trainradio=0.8
methodID=6;
pmethod=500;

dsvm=SVM_Data.SVM_Data();
# extract from kplane file..2-13 column
dsvm.extractData(fk,[x for x in range(2,14)], fs);
fname=fklist[0]+"_data"+fklist[1];
dsvm.combineDataFile(fgdata,fs,outf=fname);

dsvm.setParameter(outsets,trainradio);
dsvm.readdata(fname,svmformat=False, reset=True)
dsvm.writeTV2SVMfile(fname,methodID=methodID, outputsets=outsets,parameter=pmethod)
