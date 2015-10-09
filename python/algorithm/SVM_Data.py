#! /usr/bin/env python
# -*- coding: utf8 -*-

'''
Author: Hom, Date: 2015.8.22, Update: 2015.10.8
To pre-process SVM input for training set and validation set.
User can set the training set radio, defaul is 0.8
User can also assign the output sets number, defalut is 1

## MethodID for output method(LCM means Least Common Multiple method):
# 1. Output by shuffle list (Origin train set, origin class probability)
# 2. Output by shuffle list (LCM on train set, equal class probability)
# 3. Output by random data in the class(Origin train set, origin class probability)
# 4. Output by random data in the class(LCM on train set, equal class probability)
# 5. Output by random data in all data(Given size from train set, origin class probability)
# 6. Output by random data in random class(Given size from train set, equal class probability)

Usage: python SVM_Data.py input.txt [outsetNumber trainingRadio methodID methodParameter] 
methodParameter for method 2,4 for LCM multiple; method 5,6 for total output data number.
'''

print __doc__

import os,sys,random
from math import *

##  largest common divisor
def LargestCommon(m,n): 
	# 穷举法
	##return max([x for x in range(min(m, n),0,-1) if m%x==0 and n%x==0])

	# 辗转相除法
	if (n==0): 
		return m;
	else:
		return LargestCommon(n,m%n);

## least common multiple
def LeastCommon(m,n):
	return m*n/LargestCommon(m,n); 

## multi-number largest common divisor
def multiLargestCommon(nums):
	if (len(nums)<=0):
		print "Error: No input for multi-number largest common divisor";
		exit(1);
	elif(len(nums)==1):
		return nums[0];
	value=nums[0];
	for num in nums[1:]:
		value=LargestCommon(value,num);
	return value;

## multi-number least common multiple
def multiLeastCommon(nums):
	if (len(nums)<=0):
		print "Error: No input for multi-number least common multiple";
		exit(1);
	elif(len(nums)==1):
		return nums[0];
	value=nums[0];
	for num in nums[1:]:
		value=LeastCommon(value,num);
	return value;
	#return i0*i1*i2.../multiLargestCommon(nums)^n

class SVM_Data:
	# Contain the datas and process method. 
	# ! setParameter function to set the Output set number and Training set radio
	# ! readata function to read data into object
	# ! RandomTrainData function to randomize and initial train set/valid set according to self.tradio
	# ! writeTV2SVMfile is method to choose following method:
	# ! writeTV2SVMfileOrigin function:
	#    Output train/valid sets as origin data or least common method (even data for class).
	# ! writeTV2SVMfileRandom function: 
	#    Output train/valid sets based on randomly choose data from a class as origin or least common method
	# ! writeTV2SVMfileRandomTotal function: 
	#    Output train/valid sets based on randomly choose a data total data or from a random class.
	#    Output total number can be assigned. Don't need least common method
	# !!! Thoes three methods can output data based on origin class size 
	# !!!   or even class size(least common method or random class methond).
	# !!!   The difference is , writeTV2SVMfile output data as origin train set, 
	# !!!   writeTV2SVMfileRandom output data randomly in a class, maybe not equal amount for each data.
	# !!!   writeTV2SVMfileRandomTotal randomly select data from total train set and output given total data.

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
	# multi-times based on Least Common for each class
	classLC={}
	
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

	#string a data without class id (save in list) to given format
	#No blank in start/end of string!
	def strdata(self,data, svmformat=True,sep=" "):
		if (len(data)==0):return "";
		if (svmformat):
			index=1;
			out=str(index)+":"+str(data[0]);
			index+=1;
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
				if writeClassID: 
					fw.write(str(int(i)+" "));
				fw.write(self.strdata(data,svmformat=True));
				fw.write('\n')
		fw.close()

	# Output all the data with classid. Can read by Excel :)
	def writeAll2CSVfile(self, filename,writeClassID=True):
		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+",");
				fw.write(self.strdata(data,svmformat=False,sep=','));
				fw.write('\n')
		fw.close()

	# Output all the data with classid. Sep=" "
	def writeAll2PRNfile(self, filename,writeClassID=True):
		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+" ");
				fw.write(self.strdata(data,svmformat=False,sep=' '));
				fw.write('\n')
		fw.close()

	# Reset and Randomize the train set and valid set element
	def RandomTrainData(self):
		self.resetTV();
		# Just check the parameter
		self.setParameter(nset=self.nset,trainingradio=self.tradio)
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

	def writeValidData(self,filename):
		fvalid=open(filename,'w');
		for i in self.classid:
			for data in self.valid[i]:
				fvalid.write(str(int(i))+" "+self.strdata(data,svmformat=True));
				fvalid.write('\n')	
		fvalid.close()

	# Calculate multi-times for each class amount based on least common
	# LeastCommon: 
	#   0: don't use LeastCommon method.
	#   1: use LeastCommon method.
	#   >1(int): use LeastCommon method and each set amount multiply the value
	def CalcLeastCommon(self,LeastCommon=0):
		self.classLC.clear()
		mlc=1;
		if (LeastCommon>0):
			classLens=[]
			for c in self.classid:
				l=len(self.classdata[c])
				#exclude len=0 (blank class?)
				if (l!=0):classLens.append(l);
			mlc=multiLeastCommon(classLens);
			for c in self.classid:
				self.classLC[c]=int(mlc*LeastCommon/len(self.classdata[c]));
		elif (LeastCommon==0):
			for c in self.classid:
				self.classLC[c]=1;
		else:
			raise ValueError("Least Common is not right: "+str(LeastCommon));

	# write out data in training/valid set to corresponding file in svm format	
	# if outset=0 (not given), use the nset attr. in obj.
	# multiLeastCommon: 
	#   0: don't use LeastCommon method. Amount in class as origin.
	#   1: use LeastCommon method. Amount in class is equal.
	#   >1(int): use LeastCommon method and each set amount multiply the value
	# Note: No randomly selecting data when write out! 
	def writeTV2SVMfileOrigin(self, filename, outputsets=0,LeastCommon=0):
		fnamelist=os.path.splitext(fname);
		if (outputsets==0): outputsets=self.nset;
		self.CalcLeastCommon(LeastCommon);

		for outset in range(outputsets):
			self.RandomTrainData();
			trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
			validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
			self.writeValidData(validfname);
			ftrain=open(trainfname,'w');
			for i in self.classid:
				for data in self.train[i]:
					#write multi times bases on Least Common
					for mlc in range(self.classLC[i]):
						ftrain.write(str(int(i))+" "+self.strdata(data,svmformat=True));
						ftrain.write('\n')					
			ftrain.close()

	# write out data in training/valid set to corresponding file in svm format	
	# if outset=0 (not given), use the nset attr. in obj.
	# multiLeastCommon: 
	#   0: don't use LeastCommon method. Amount in class as origin.
	#   1: use LeastCommon method. Amount in class is equal.
	#   >1(int): use LeastCommon method and each set amount multiply the value
	# Note: Data in a class is randomly selected! Recommend to expand the data amount by LC method! 
	def writeTV2SVMfileRandom(self, filename, outputsets=0,LeastCommon=0):
		fnamelist=os.path.splitext(fname);
		if (outputsets==0): outputsets=self.nset;
		self.CalcLeastCommon(LeastCommon);

		for outset in range(outputsets):
			self.RandomTrainData();
			trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
			validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
			self.writeValidData(validfname);
			ftrain=open(trainfname,'w');
			for i in self.classid:
				classNow=self.train[i];
				count=len(classNow);
				#write multi times bases on Least Common(classlen*classLC[i])
				for mlc in range(self.classLC[i]*count):
					data=classNow[random.randint(0,count-1)];
					ftrain.write(str(int(i))+" "+self.strdata(data,svmformat=True));
					ftrain.write('\n')					
			ftrain.close()

	# write out data in training/valid set to corresponding file in svm format	
	# if outset=0 (not given), use the nset attr. in obj.
	# TotalDataNum: Output total data number
	#   0: Use default data number as output data number
	# randomClass: Whether use random class method to make class probability equal.
	# Note: Data in a class is randomly selected!  
	def writeTV2SVMfileRandomTotal(self, filename, outputsets=0,TotalDataNum=0,randomClass=True):
		fnamelist=os.path.splitext(fname);
		if (outputsets==0): outputsets=self.nset;
		if (TotalDataNum==0): 
			TotalDataNum=self.ndata;
		TotalDataNum=int(TotalDataNum);
		nclass=len(self.classid);
		totaldata=[]
		totaldataID=[]
		numlist=[]
		for outset in range(outputsets):
			self.RandomTrainData();
			del totaldata[:]
			del totaldataID[:]
			del numlist[:]
			for c in self.classid:
				# Get total data
				if (not randomClass): 
					totaldataID=totaldataID+[c]*len(self.train[c]);
					totaldata=totaldata+self.train[c];
				# Get length of each class
				else: 
					numlist.append(len(self.train[c]));
			ndata=len(totaldata);
			trainfname=fnamelist[0]+"_train_"+str(outset+1)+fnamelist[1];
			validfname=fnamelist[0]+"_valid_"+str(outset+1)+fnamelist[1];
			self.writeValidData(validfname);
			ftrain=open(trainfname,'w');

			for i in range(TotalDataNum):
				data=[]
				if (randomClass):
					# random class number
					ranc=random.randint(0,nclass-1);
					# random data number
					ranv=random.randint(0,numlist[ranc]-1);
					data=self.train[self.classid[ranc]][ranv]
					ftrain.write(str(int(self.classid[ranc]))+" "+self.strdata(data,svmformat=True));
				else:
					idx=random.randint(0,ndata-1)
					data=totaldata[idx]
					ftrain.write(str(int(totaldataID[idx]))+" "+self.strdata(data,svmformat=True));
				ftrain.write('\n')					
			ftrain.close()

	def writeTV2SVMfile(self, filename, methodID=1, outputsets=0,paramete=0):
		if (methodID>6 or methodID<1):
			raise ValueError("Error Method ID!: "+str(methodID));
		if (methodID==1):
			self.writeTV2SVMfileOrigin(filename, outputsets=outputsets,LeastCommon=0);
		elif (methodID==2):
			self.writeTV2SVMfileOrigin(filename, outputsets=outputsets,LeastCommon=paramete);
		elif (methodID==3):
			self.writeTV2SVMfileRandom(filename, outputsets=outputsets,LeastCommon=0);
		elif (methodID==4):
			self.writeTV2SVMfileRandom(filename, outputsets=outputsets,LeastCommon=paramete);
		elif (methodID==5):
			self.writeTV2SVMfileRandomTotal(filename, outputsets=outputsets,TotalDataNum=paramete,randomClass=False);
		elif (methodID==6):
			self.writeTV2SVMfileRandomTotal(filename, outputsets=outputsets,TotalDataNum=paramete,randomClass=True);

if (__name__ == '__main__'):

	## Default value
	trainradio=0.8;
	outsets=1;
	methodID=1;
	pmethod=1;

	if (len(sys.argv) < 2):
		print "Please give an input SVM data file!";
		exit(1);
	elif(len(sys.argv)>=3):
		if (sys.argv[2].isdigit()):
			outsets=int(sys.argv[2]);
		else: 
			print "Output sets should be a digit!";
			exit(1);
		if (len(sys.argv)>=4):
			trainradio=float(sys.argv[3]);
			if (trainradio>1.0 or trainradio<0):
				print "Training set radio should be [0,1]";
				exit(1);
		if (len(sys.argv)>=5):
			if (sys.argv[4].isdigit()):
				methodID=int(sys.argv[4]);
				if (methodID>6 or methodID<1):
					print "Method should be 1 to 6";
					exit(1);
			else: 
				print "Method ID should be a digit!";
				exit(1);
		if (len(sys.argv)>=6):
			if (sys.argv[5].isdigit()):
				pmethod=float(sys.argv[5]);
				if (pmethod<=0):
					print "Parameter for method should be > 0";
					exit(1);
			else: 
				print "Parameter for method should be a digit!";
				exit(1);

	fname=sys.argv[1];
	dsvm=SVM_Data();
	dsvm.setParameter(outsets,trainradio);
	dsvm.readdata(fname,svmformat=True, reset=True)
	dsvm.writeTV2SVMfile(fname,methodID=methodID, outputsets=outsets,paramete=pmethod)
	
#end main