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
# 7. Output for N-fold cross-validation (Generate N model, one set as validation)

Usage: python SVM_Data.py input.txt [outsetNumber trainingRadio methodID methodParameter] 
methodParameter in method 2,4 for LCM multiple factor; method 5,6 for total output data number.
'''

__author__	=	"Zhixiong Zhao"
__date__	=	"2015.10.8"
__version__	=	"1.0"
__email__	=	"zhaozxcpu@hotmail.com"

import os,sys,random
import math
import numpy as np


def LargestCommon(m,n): 
	'largest common divisor'
	# 穷举法
	##return max([x for x in range(min(m, n),0,-1) if m%x==0 and n%x==0])

	# 辗转相除法
	if (n==0): 
		return m;
	else:
		return LargestCommon(n,m%n);

def LeastCommon(m,n):
	'least common multiple'
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

def multiLeastCommon(nums):
	'multi-number least common multiple'
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

	'''
	# Contain the datas and process method. 
	# ! *setParameter function to set the Output set number and Training set radio
	# ! *readdata function to read data into object
	# ! *RandomTrainData function to randomize and initial train set/valid set according to self.tradio
	# ! *writeTV2SVMfile is method to choose the following method:
	# !  writeTV2SVMfileOrigin function:
	#     Output train/valid sets as origin data or least common method (even data for class).
	# !  writeTV2SVMfileRandom function: 
	#     Output train/valid sets based on randomly choose data from a class as origin or least common method
	# !  writeTV2SVMfileRandomTotal function: 
	#     Output train/valid sets based on randomly choose a data total data or from a random class.
	#     Output total number can be assigned. Don't need least common method
	# !! Thoes three methods can output data based on origin class size 
	# !!   or even class size(least common method or random class methond).
	# !!   The difference is , writeTV2SVMfile output data as origin train set, 
	# !!   writeTV2SVMfileRandom output data randomly in a class, maybe not equal amount for each data.
	# !!   writeTV2SVMfileRandomTotal randomly select data from total train set and output given total data.
	'''
	# class id sets
	classid=[]

	# data in each class (dict)
	classdata={}
	
	# data in training set (dict)
	train={}
	
	# data in validation set (dict)
	valid={}
	
	# data in test set (in list, not dict! No class id for them)
	test=[]
	
	# number of set for output
	nset=1
	
	# training set radio
	tradio=0.8
	
	# totoal number of data
	ndata=0
	
	# multi-times based on Least Common for each class
	classLC={}
	
	def __init__(self):
		'''SVM_Data object'''
		pass
	
####### Basic setting function ##########
	
	def clear(self):
		# clear all the data

		del self.classid[:];
		del self.test[:];
		self.classdata.clear()
		self.train.clear();
		self.valid.clear();
		self.classLC.clear();
	
	def resetTV(self):
		# reset train and valid data set

		self.train.clear()
		self.valid.clear()

	def setParameter(self,nset=1,trainingradio=0.8):
		# Recommend to use to set parameter before obtain train/valid set

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

	def CalcLeastCommon(self,LeastCommon=0):
		# Calculate multi-times for each class amount based on least common
		# LeastCommon: 
		#   0: don't use LeastCommon method.
		#   1: use LeastCommon method.
		#   >1(int): use LeastCommon method and each set amount multiply the value

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

	def RandomTrainData(self):
		# Reset and Randomize the train set and valid set element

		self.resetTV();
		# Just check the parameter
		self.setParameter(nset=self.nset,trainingradio=self.tradio)
		for c in self.classid:
			classNow=self.classdata[c];
			count=len(classNow);
			self.train[c]=[]
			self.valid[c]=[]

			# Exception: only one data in a class. 
			if (count==1):
				# if only one data, put it to training set.
				self.train[c].append(classNow[0]);
				continue;

			# training set number (tradio: trainsing set radio)
			Ntrain=int(count*self.tradio);
			# Get a random list
			shufflelist=range(count);
			random.shuffle(shufflelist);

			# get the number of data being used in set
			tlist=shufflelist[:Ntrain];
			tlist.sort()
			vlist=shufflelist[Ntrain:];
			vlist.sort()
			# Set the train / valid set
			for i in tlist:
				self.train[c].append(classNow[i]);
			for i in vlist:
				self.valid[c].append(classNow[i]);	

####### Read Data related function ##########

	def PerceiveSVMline(self,line):
		# Perceive SVM data line to a list

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

	def readdata(self,filename,svmformat=False, reset=True, sep=""):
		'''Read data from a file into the object.
		 1: input filename, file should be classID and features. 
		 2: Whether svm input data
		 3: seperator for the data in line'''

		# Default to clear all data
		if reset: self.clear()
		ndata=0
		fr=open(filename)
		for line in fr:
			tmp=[]
			# seperate the data in line
			if (svmformat): tmp=self.PerceiveSVMline(line)
			else:	
				if(sep==""):
					tmp=line.strip().split()
				else:
					tmp=line.strip().split(sep)
			
			#class name + data
			if (len(tmp)>1):
				try:
					# int('6.0e+00') error
					classNum=str(int(float(tmp[0])))
					if (int(classNum)>0):
						ndata+=1
						if (classNum in self.classdata):
							self.classdata[classNum].append(tmp[1:]);
						else:
							self.classdata[classNum]=[tmp[1:]];
							self.classid.append(classNum);
				except ValueError:
					pass
		self.ndata=ndata
		fr.close()

####### Write Data related function ##########

	def strdata(self,data, svmformat=True,sep=" "):
		#string a data without class id (save in list) to given format
		#No blank in start/end of string!

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
					
	def writeAll2SVMfile(self, filename,writeClassID=True):
		# write out all the data to a file in svm format	
		# Default also write out the class id, can set off	

		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: 
					fw.write(str(int(i)+" "));
				fw.write(self.strdata(data,svmformat=True));
				fw.write('\n')
		fw.close()

	def writeAll2CSVfile(self, filename,writeClassID=True):
		# Output all the data with classid. Can read by Excel :)

		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+",");
				fw.write(self.strdata(data,svmformat=False,sep=','));
				fw.write('\n')
		fw.close()
	
	def writeAll2PRNfile(self, filename,writeClassID=True):
		# Output all the data with classid. Sep=" "

		fw=open(filename,'w');
		for i in self.classid:
			for data in self.classdata[i]:
				if writeClassID: fw.write(str(int(i))+" ");
				fw.write(self.strdata(data,svmformat=False,sep=' '));
				fw.write('\n')
		fw.close()

	def writeValidData(self,filename):
		# write valid set data to a file

		fvalid=open(filename,'w');
		for i in self.classid:
			for data in self.valid[i]:
				fvalid.write(str(int(i))+" "+self.strdata(data,svmformat=True));
				fvalid.write('\n')	
		fvalid.close()

############# Main function for output Train/Valid Set ################

	def writeTV2SVMfileOrigin(self, filename, outputsets=0,LeastCommon=0):
		# write out data in training/valid set to corresponding file in svm format	
		# if outset=0 (not given), use the nset attr. in obj.
		# multiLeastCommon: 
		#   0: don't use LeastCommon method. Amount in class as origin.
		#   1: use LeastCommon method. Amount in class is equal.
		#   >1(int): use LeastCommon method and each set amount multiply the value
		# Note: No randomly selecting data when write out! 

		fnamelist=os.path.splitext(filename);
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

	def writeTV2SVMfileRandom(self, filename, outputsets=0,LeastCommon=0):
		# write out data in training/valid set to corresponding file in svm format	
		# if outset=0 (not given), use the nset attr. in obj.
		# multiLeastCommon: 
		#   0: don't use LeastCommon method. Amount in class as origin.
		#   1: use LeastCommon method. Amount in class is equal.
		#   >1(int): use LeastCommon method and each set amount multiply the value
		# Note: Data in a class is randomly selected! Recommend to expand the data amount by LC method!

		fnamelist=os.path.splitext(filename);
		if (outputsets==0): outputsets=self.nset;
		self.CalcLeastCommon(LeastCommon);

		for outset in range(outputsets):
			# Reset and Randomize the train/valid set
			self.RandomTrainData();
			# Output valid set
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

	def writeTV2SVMfileRandomTotal(self, filename, outputsets=0,TotalDataNum=0,randomClass=True):
		# write out data in training/valid set to corresponding file in svm format	
		# if outset=0 (not given), use the nset attr. in obj.
		# TotalDataNum: Output total data number
		#   0: Use default data number as output data number
		# randomClass: Whether use random class method to make class probability equal.
		# Note: Data in a class is randomly selected!  

		fnamelist=os.path.splitext(filename);
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

	def writeTV2SVMfileNfold(self, filename, outputsets=10):
		# write out data in training/valid set to corresponding file in svm format
		# Use N-fold method	
		# if outset=0 (not given), use the nset attr. in obj.
		# Note: Data were randomly group first and then one group as validation, write out N model.

		fnamelist=os.path.splitext(filename);
		if (outputsets==0): outputsets=self.nset;

		for outset in range(outputsets):
			# Reset and Randomize the train/valid set
			self.RandomTrainData();
			# Output valid set
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

	def writeTV2SVMfile(self, filename, methodID=1, outputsets=0,parameter=0):
		# To control output train/valid set method based methodID
		# parameter for LCM is the multiple factor; for randomTotal is the training set size.

		if (methodID>7 or methodID<1):
			raise ValueError("Error Method ID!: "+str(methodID));

		if (methodID==1):
			self.writeTV2SVMfileOrigin(filename, outputsets=outputsets,LeastCommon=0);
		elif (methodID==2):
			self.writeTV2SVMfileOrigin(filename, outputsets=outputsets,LeastCommon=parameter);
		elif (methodID==3):
			self.writeTV2SVMfileRandom(filename, outputsets=outputsets,LeastCommon=0);
		elif (methodID==4):
			self.writeTV2SVMfileRandom(filename, outputsets=outputsets,LeastCommon=parameter);
		elif (methodID==5):
			self.writeTV2SVMfileRandomTotal(filename, outputsets=outputsets,TotalDataNum=parameter,randomClass=False);
		elif (methodID==6):
			self.writeTV2SVMfileRandomTotal(filename, outputsets=outputsets,TotalDataNum=parameter,randomClass=True);

########## Post/Pre-process function ##########

	def combineDataFile(self,fname1,fname2,outf="",sep=" "):
		# Combine two data files to a outfile.
		# sep for separator
		# NOTE: Only done when data size is the same
		fnamelist1=os.path.splitext(fname1);
		fnamelist2=os.path.splitext(os.path.basename(fname2));

		if (outf==""):
			outf=fnamelist1[0]+"_"+fnamelist2[0]+fnamelist1[1];
		fr1=open(fname1);
		fr2=open(fname2);
		fw=open(outf,'w');
		lines1=[]
		lines2=[]
		for line in fr1:
			if (line.strip()): lines1.append(line)
		for line in fr2:
			if (line.strip()): lines2.append(line)
		if (len(lines1)==len(lines2)):
			for i in range(len(lines1)):
				fw.write(lines1[i].strip()+sep+lines2[i]);
		else:
			print "Data in in file is not equal!: "+str(len(lines1))+","+str(len(lines2))
		fr1.close()
		fr2.close()
		fw.close()

	def extractData(self,fname,columnList,outf="",sep=" "):
		# Extract given column from a data file.
		# Combine use with combineDataFile() to set a suitable input data file
		# column list start from 1, and in list [1, 3, 5]
		fnamelist=os.path.splitext(fname);
		if (outf==""):
			outf=fnamelist[0]+"_extract"+fnamelist[1];
		fr=open(fname);
		fw=open(outf,'w');
		for line in fr:
			tmp=line.strip().split();
			outlist=[];
			for i in range(len(tmp)):
				if ((i+1) in columnList):
					outlist.append(tmp[i]);
			fw.write(sep.join(outlist)+"\n");
			del outlist[:];
		fr.close()
		fw.close()

	def CalcTestSetRMS(self,svmout, wtsout, yfile, svmcolumn=1, ycolumn=1):
		# To calculate the RMS. Test 
		# svmout include the class number of each data in test set
		# wtsout include the weight of each class
		
		if (not self.test): 
			raise ValueError("No test set!");

		if (svmcolumn<1):
			raise ValueError("Column number for SVM class should be >=1!");

		if (ycolumn<1):
			raise ValueError("Column number for y value should be >=1!");

		f2=open(wtsout)
		f1=open(svmout)
		samefile=(yfile==svmout)
		if (not samefile): f3=open(yfile)
		svmc=[]
		wts=[]
		y=[]

		for line in f1:
			tmp=line.strip().split()
			if (len(tmp)>0):
				svmc.append(tmp[svmcolumn-1])
				if samefile:
					y.append(tmp[ycolumn-1])

		for line in f2:
			tmp=line.strip().split()
			if (len(tmp)>0):
				wts.append(tmp);

		if not samefile:
			for line in f3:
				tmp=line.strip().split()
				if (len(tmp)>0):
					y.append(tmp[ycolumn-1])				

		# check data:
		if (len(svmc) is 0 or len(wts) is 0 or len(y) is 0):
			raise ValueError("Size of svm output or weight or y is zero!")

		if (len(svmc)!=len(test) or len(svmc)!=len(y)):	
			raise ValueError("Data size is not equal! test:"
				+str(len(test))+"; svm:"+str(len(svmc))
				+"; y:"+str(len(y)) )

		if (len(test[0])!=len(wts[0])):
			raise ValueError("Features number is not equal! test: "
				+str(len(test[0]))+"; weight:"+str(len(wts[0])) )

		nclass=len(wts)
		svmcs=set(svmc)
		for i in svmcs:
			if (int(i)>nclass):
				raise ValueError("Class number in svm is off weight border! Nweight: "
					+str(nclass)+"; Class number in svm:"+str(i))

		testR=np.array(self.test,dtype=np.float64)
		wtsR=np.array(wts,dtype=np.float64)
		errors=np.zero((len(svmc),1))

		for i in range(len(svmc)):
			tmpx=testR[i,:]
			tmpwts=wtsR[svmc[i],:]
			errors[i]=(np.abs(float(y[i])-(tmpwts[0, :].dot(tmpx.T)).T)).T;
		rms=np.sqrt((errors.T.dot(errors)).ravel()[0]/errors.shape[0])
		print rms
		return rms

if (__name__ == '__main__'):
	
	print __doc__

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
	dsvm.writeTV2SVMfile(fname,methodID=methodID, outputsets=outsets,parameter=pmethod)
	
#end main