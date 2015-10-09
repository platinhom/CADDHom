#! /usr/bin/env python
# -*- coding: utf8 -*-

'''
Author: Hom, Date: 2015.8.22, Update: 2015.10.8
To pre-process SVM input for training set and validation set.
User can set the training set radio, defaul is 0.8
User can also assign the output sets number, defalut is 1
Usage: python SVM_trainvalid.py input.txt [outsets trainingradio] 
'''

print __doc__

import os,sys,random
from math import *

## Default value
trainradio=0.8;
outsets=1;
newline=False;

## Method for output:
# 1. Output many sets with shuffle list (Origin and All data)
# 2. Output many sets with random lines in a class(Possible repeated lines)
# 3. Output many sets with random lines in all classes (probability depends on origin data)
# 4. Output many sets with random lines in all classes (probability of each class is equal)
method=4;

## Data is all the lines in list.

class SVM_Data:
	# Contain the datas and process method. 
	# ! setParameter function to set the Output set number and Training set radio
	# ! readata function to read data into object
	# ! RandomTrainData function to randomize and initial train set/valid set according to self.tradio
	# ! writeTV2SVMfile function:
	#    Write out train/valid sets (multi-outset)
	# ! writeTV2SVMfileRandom function: 
	#    Write out train/valid sets (multi-outset) based on randomly choose data from a class
	#    Can expand total number of data to self.ndata.
	# ! writeTV2SVMfileRandom: 

	# class id sets
	classid=[]
	# data in each class (dict)
	classdata={}
	# data in training set (dict)
	train={}
	# data in validation set (dict)
	valid={}
	# number of set for output
	nset=1
	# training set radio
	tradio=0.8
	# totoal number of data
	ndata=0
	# output data number of each set
	nout=0
	
	def __init__(self):
		pass
	
	# clear all the data
	def clear(self):
		del self.classid[:];
		self.classdata.clear()
		self.train.clear();
		self.valid.clear();
	
	# reset train and valid data set
	def resetTV(self):
		self.train.clear()
		self.valid.clear()

	# Recommend to use to set parameter before obtain train/valid set
	def setParameter(self,nset=1,trainingradio=0.8):
		if (nset<1):
			print "Output set number should >1! Now: "+str(nset);
			print "Reset Output set number to 1"
			nset=1;
		self.nset=nset;
		if (trainingradio>1 or trainingradio<0):
			print "Training set radio should [0,1]! Now: "+str(trainingradio);
			print "Reset Training set radio to 0.8"
			trainingradio=0.8;
		self.trainingradio=trainingradio;
	
	# Perceive SVM data line
	def PerceiveSVMline(self,line):
		tmp=line.strip().split()
		if (len(tmp)>1):
			data=[tmp[0]]
			for item in tmp[1:]:
				p=item.split(':')
				idx=int(p[0]);
				ld=len(data)
				if (idx==ld):
					data.append(p[1])
				# svm can 1:1 5:3 10:2, other 0
				elif (idx>ld):
					data=data+["0"]*(idx-ld)+p[1]
				# should never happen in normal case
				else : data[idx]=p[1]
			return data
		else: return []

	# input filename 
	# Whether svm input data
	# seperator for the data in line
	def readdata(self,filename,svmformat=False, reset=True, sep=""):
		# Default to clear all data
		if reset: self.clear()
		ndata=0
		fr=open(filename)
		for line in fr:
			tmp=[]
			# seperate the data in line
			if (svmformat): tmp=self.PerceiveSVMline(line)
			else:	tmp=line.strip().split(sep)
			
			#class name + data
			if (len(tmp)>1):
				classNum=tmp[0]
				if (classNum.isdigit()):
					 ndata+=1
					if (classNum in self.classdata):
						self.classdata[classNum].append(tmp[1:]);
					else:
						self.classdata[classNum]=[tmp[1:]];
						self.classid.append(classNum);
		self.ndata=ndata
		self.nout=ndata

	#string a data without class id (save in list) to given format
	#No blank in start/end of string!
	def strdata(self,data, svmformat=True,sep=" "):
		if (len(data)==0):return "";
		if (svmformat):
			index=1;
			out=str(data[0]);
			for item in data[1:]:
				out+=" "+str(index)+":"+str(item)
				index+=1;
			return out
		else:
			return sep.join(data)
	
	# write out all the data to a file in svm format	
	# Default also write out the class id, can set off						
	def writeAll2SVMfile(self, filename,writeClassID=True):
		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i)+" ");
				fw.write(strdata(data,svmformat=True));
				fw.write('\n')
		fw.close()

	# Output all the data with classid. Can read by Excel :)
	def writeAll2CSVfile(self, filename,writeClassID=True):
		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+",");
				fw.write(strdata(data,svmformat=False,sep=','));
				fw.write('\n')
		fw.close()

	# Output all the data with classid. Sep=" "
	def writeAll2PRNfile(self, filename,writeClassID=True):
		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+" ");
				fw.write(strdata(data,svmformat=False,sep=' '));
				fw.write('\n')
		fw.close()

	# Reset and Randomize the train set and valid set element
	def RandomTrainData(self):
		self.resetTV();
		# Just check the parameter
		setParameter(self,nset=self.nset,trainingradio=self.tradio)
		for c in self.classid:
			classNow=self.classdata[c];
			count=len(classNow);
			# training set number (tradio: trainsing set radio)
			Ntrain=int(count*self.tradio);
			self.train[c]=[]
			self.valid[c]=[]
			# Get a random list
			shufflelist=range(count);
			random.shuffle(shufflelist);
			# get the number of data being used in set
			tlist=shufflelist[:Ntrain];
			tlist.sort()
			vlist=shufflelist[Ntrain:];
			vlist.sort()
			for i in tlist:
				self.train[c].append(classNow[i]);
			for i in vlist:
				self.valid[c].append(classNow[i]);	

	# write out data in training/valid set to corresponding file in svm format	
	# if outset=0 (not given), use the nset attr. in obj.
	# Note: Don't expand the set or randomly select data! 
	def writeTV2SVMfile(self, filename, outputsets=0):
		fnamelist=os.path.splitext(fname);
		if (outputsets==0): outputsets=self.nset;

		for outset in range(outputsets):
			trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
			validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
			ftrain=open(trainfname,'w');
			fvalid=open(validfname,'w');
			for i in self.classid:
				for data in self.train[i]:
					ftrain.write(str(int(i))+strdata(data,svmformat=True));
					ftrain.write('\n')
				for data in self.valid[i]:
					fvalid.write(str(int(i))+strdata(data,svmformat=True));
					fvalid.write('\n')					
			ftrain.close()
			fvalid.close()


## Only retain class number is digit
## data here file handle or a list contains all lines
def filterdata(data):
	fdata=[];
	for line in data:
		classNum=line.strip().split()[0];
		if (classNum.isdigit()):
			fdata.append(line);
	return fdata;
				
## Analyse and classify the data
## Return: [ Class numbers in list , classified lines in dict ] 
def countclass(data):
	classdict={};
	classlist=[];
	for line in data:
		classNum=line.strip().split()[0];
		if (classNum.isdigit()):
			if (classNum in classdict):
				classdict[classNum].append(line);
			else:
				classdict[classNum]=[line];
				classlist.append(classNum);
	return [classlist,classdict];

## Write out the training and validation sets (multi-model) in a given radio
def writeTrainValid(classlist, classdict, fname='data.txt',outputsets=1, trainsetradio=0.8, newline=False):
	fnamelist=os.path.splitext(fname);
	for outset in range(outputsets):
		trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
		validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
		ftrain=open(trainfname,'w');
		fvalid=open(validfname,'w');
		for c in classlist:
			classNow=classdict[c];
			count=len(classNow);
			Ntrain=int(count*trainsetradio);
			#Nvalid=count-Ntrain;
			shufflelist=range(count);
			random.shuffle(shufflelist);
			if (newline):
				for i in shufflelist[:Ntrain]:
					ftrain.write(classNow[i]+'\n');
				for i in shufflelist[Ntrain:]:
					fvalid.write(classNow[i]+'\n');	
			else:
				for i in shufflelist[:Ntrain]:
					ftrain.write(classNow[i]);
				for i in shufflelist[Ntrain:]:
					fvalid.write(classNow[i]);	
		ftrain.close();
		fvalid.close();

## Write out the training and validation sets (multi-model) in a given radio
## Randomize the output lines in a class, meaning that some lines may be used many times.
def writeRandomTrainValid(classlist, classdict, fname='data.txt',outputsets=1, trainsetradio=0.8, newline=False):
	fnamelist=os.path.splitext(fname);
	for outset in range(outputsets):
		trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
		validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
		ftrain=open(trainfname,'w');
		fvalid=open(validfname,'w');
		j=0;
		for c in classlist:
			classNow=classdict[c];
			count=len(classNow);
			Ntrain=int(count*trainsetradio);
			Nvalid=count-Ntrain;
			if (newline):
				for i in range(Ntrain):
					ftrain.write(classNow[random.randint(0,count-1)]+'\n');
				for i in range(Nvalid):
					fvalid.write(classNow[random.randint(0,count-1)]+'\n');	
			else:
				for i in range(Ntrain):
					ftrain.write(classNow[random.randint(0,count-1)]);
				for i in range(Nvalid):
					fvalid.write(classNow[random.randint(0,count-1)]);	
		ftrain.close();
		fvalid.close();

# Get all datas first and random get data to train set/valid set
def writeRandomAllTrainValid(fdata, fname='data.txt',outputsets=1, trainsetradio=0.8, newline=False):
	fnamelist=os.path.splitext(fname);
	count=len(fdata);
	#Ntrain=int(count*trainsetradio);
	#Nvalid=count-Ntrain;		
	Ntrain=1000;
	Nvalid=250;
	for outset in range(outputsets):
		trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
		validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
		ftrain=open(trainfname,'w');
		fvalid=open(validfname,'w');
		if (newline):
			for i in range(Ntrain):
				ranv=random.randint(0,count-1);
				ftrain.write(fdata[ranv]+'\n');
			for i in range(Nvalid):
				ranv=random.randint(0,count-1);
				fvalid.write(fdata[ranv]+'\n');	
		else:
			for i in range(Ntrain):
				ranv=random.randint(0,count-1);
				ftrain.write(fdata[ranv]);
			for i in range(Nvalid):
				ranv=random.randint(0,count-1);
				fvalid.write(fdata[ranv]);	
		ftrain.close();
		fvalid.close();

##  largest common divisor
def gongyue(m,n): 
	# 穷举法
    ##return max([x for x in range(min(m, n),0,-1) if m%x==0 and n%x==0])

    # 辗转相除法
	if (n==0): 
		return m;
	else:
		return gongyue(n,m%n);

## least common multiple
def gongbei(m,n):
	return m*n/gongyue(m,n); 

## multi-number largest common divisor
def multigongyue(nums):
	if (len(nums)<=0):
		print "Error: No input for multi-number largest common divisor";
		exit(1);
	elif(len(nums)==1):
		return nums[0];
	value=nums[0];
	for num in nums[1:]:
		value=gongyue(value,num);
	return value;

## multi-number least common multiple
def multigongbei(nums):
	if (len(nums)<=0):
		print "Error: No input for multi-number least common multiple";
		exit(1);
	elif(len(nums)==1):
		return nums[0];
	value=nums[0];
	for num in nums[1:]:
		value=gongbei(value,num);
	return value;
	#return i0*i1*i2.../multigongyue(nums)^n  

# random get a class
# then random get a data from a class
def writeRandomEqualAllTrainValid(classlist, classdict, fname='data.txt',outputsets=1, trainsetradio=0.8, newline=False):
	fnamelist=os.path.splitext(fname);
	numlist=[];
	for i in classlist:
	 	numlist.append(len(classdict[i]));
	# lcm=multigongbei(numlist);
	# clcm=lcm*len(classlist);
	# if (clcm<=10):
	# 	lcm=lcm*100;
	# 	clcm*=100;
	# elif (clcm<=100):
	# 	lcm=lcm*10;
	# 	clcm*=10;

	count=len(classlist);
	#Ntrain=int(count*trainsetradio);
	#Nvalid=count-Ntrain;		
	Ntrain=1000;
	Nvalid=250;
	for outset in range(outputsets):
		trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
		validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
		ftrain=open(trainfname,'w');
		fvalid=open(validfname,'w');
		if (newline):
			for i in range(Ntrain):
				ranc=random.randint(0,count-1);
				ranv=random.randint(0,numlist[ranc]-1);
				ftrain.write(classdict[classlist[ranc]][ranv]+'\n');
			for i in range(Nvalid):
				ranc=random.randint(0,count-1);
				ranv=random.randint(0,numlist[ranc]-1);
				fvalid.write(classdict[classlist[ranc]][ranv]+'\n');
		else:
			for i in range(Ntrain):
				ranc=random.randint(0,count-1);
				ranv=random.randint(0,numlist[ranc]-1);
				ftrain.write(classdict[classlist[ranc]][ranv]);
			for i in range(Nvalid):
				ranc=random.randint(0,count-1);
				ranv=random.randint(0,numlist[ranc]-1);
				fvalid.write(classdict[classlist[ranc]][ranv]);	
		ftrain.close();
		fvalid.close();

if (__name__ == '__main__'):
	if (len(sys.argv) < 2):
		print "Please give an input SVM data file!";
		exit(1);
	elif(len(sys.argv)>=3):
		if (sys.argv[2].isdigit()):
			outsets=int(sys.argv[2]);
		else: 
			print "Output sets should be a digit!";
			exit(1);
		if (len(sys.argv)==4):
			if (sys.argv[3].isdigit()):
				trainradio=float(sys.argv[3]);
				if (trainradio>1.0 or trainradio<0):
					print "Training set radio should be [0,1]";
					exit(1);
			else: 
				print "Training set radio should be a digit!";
				exit(1);

	fname=sys.argv[1];
	fr=open(fname);
	data=fr.readlines();

	if (method==1):
		# Analyze classification of data
		countresult=countclass(data);
		# Write out training set and validation set
		# Output many sets of training/validation set based on outsets parameter
		writeTrainValid(countresult[0],countresult[1],fname,outsets,trainradio,False);
	elif (method==2):
		# Analyze classification of data
		countresult=countclass(data);
		# Write out training set and validation set with random lines in a class
		# Output many sets of training/validation set based on outsets parameter
		writeRandomTrainValid(countresult[0],countresult[1],fname,outsets,trainradio,False);
	elif (method==3):
		# Analyze classification of data
		fdata=filterdata(data);
		writeRandomAllTrainValid(fdata,fname,outsets,trainradio,False);
	elif (method==4):
		# Analyze classification of data
		countresult=countclass(data);
		writeRandomEqualAllTrainValid(countresult[0],countresult[1],fname,outsets,trainradio,False);	
	fr.close();
#end main