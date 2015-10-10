#! /usr/bin/env python
# -*- coding: utf8 -*-
'''
#######################################################################################
This program implements the K-plane regression algorithm based on EM algorithm 
to give the optimal partitioning of the data set in the sense ofd minimize the RMS
error of the regression. Together with the SVM classifer, we can do the classification.

Return:
rms is the total rms error
coef is the regression coefficients for each group
partition_groups is the group information for each data

Input:
x is the training set; x(1:num_samples, 1:number_fearures)
y is the target value: y(1:num_samples)
dgroup is the initial group information for each data if given

Parameters:
num_group is the number of groups
regroup is whether use regroup method to make sure the result with starting group number.
randomgroup is whether use randomized (or well-ordering method) group information in initialization.
print_in is whether print out the intermediate information.
klamda is the regularization constant for k-plane
rlamda is the regularization constant for least sqaure regression
rmslimit is the convergence limit for RMS 

Author: Bao Wang, wangbaonj@gmail.com, 09/30/2015 based on Matlab: basical algorithm
Reviser: Zhixiong Zhao, 10/06/2015 based on python: objective, regroup algorithm
########################################################################################
'''

import os,sys,traceback
import numpy as np

#The number of groups
num_group=8;

class KPlaneObj():
	# x: x data
	# y: y data
	# dlen: data number
	# plen: parameter numbers
	# num_group: group number
	# dgroup: group number for each data(start)
	# partition_groups: final group number for each data
	# errors: error for each data
	# centers: center for each group
	# wts: weight for each group
	def __init__(self):
		pass

	def CheckData(self,data=True,group=True):
		# Check input data and convert them to valid data type
		if (data):
			if (not isinstance(self.x,np.ndarray)):
				self.x=np.array(self.x);
			if (not isinstance(self.y,np.ndarray)):
				self.y=np.array(self.y);
			if (self.x.ndim==1):
				self.dlen=self.x.size;
				self.plen=1
				self.x=self.x.T;
			elif (self.x.ndim==2):
				self.dlen, self.plen=self.x.shape;
			else:
				raise ValueError, "x dimension > 2!";
			if (self.y.ndim==1):
				my=self.y.size;
				if (my!=self.dlen):
					raise ValueError, "data size not equal in x and y!";
				self.y=self.y.T;
			elif (self.y.ndim==2):
				my, ny=self.y.shape;
				if (my!=self.dlen):
					raise ValueError, "data size not equal in x and y!";
			else:
				raise ValueError, "y dimension > 2!";
			# At least two data
			assert (self.dlen>1 and self.plen>0)
		# check group number
		if (group):
			assert (self.num_group>0)
			if (self.dlen<2*self.num_group):
				print "Group number is larger than data number/2"
				print "Automatically set group number to "+str(int(self.dlen/2))
				self.SetGroupNumber(int(self.dlen/2));

	def SetAllParameters(self,klamda=0.001,rlamda=500, rmslimit=0.0001):
		self.klamda=klamda
		self.rlamda=rlamda
		self.rmslimit=rmslimit

	def SetGroupNumber(self,number=1):
		self.num_group=number;
		# May need further dgroup assignment

	def SetGroupData(self,dgroup):
		if (isinstance(dgroup,np.ndarray)):
			self.dgroup=dgroup;
		else:
			self.dgroup=np.array(dgroup);
		if (self.dgroup.size !=self.dlen):
			raise ValueError,"amount of group number not equal to amount of data!"
		if (self.dgroup.ndim==1):
			self.dgroup.resize((self.dlen,1));

	### Grouping function can use only after initialize
	#   Grouping Based on sequence
	def AutoGroup(self):
		num_group=self.num_group
		dlen=self.dlen
		for k in range(0,num_group-1):
			self.dgroup[k*np.floor(dlen/num_group):(k+1)*np.floor(dlen/num_group)]=k+1;
		self.dgroup[(num_group-1)*np.floor(dlen/num_group):dlen]=num_group;

	#   Grouping Based on randomize
	def RandomGroup(self):
		for i in range(self.dlen):
			self.dgroup[i]=np.random.randint(1,num_group+1);

	### Initialization on different input data
	#   Initialize based on x and y data (array or list or tuple)
	def initializeByData(self,x,y,num_group=1,regroup=True,randomgroup=False,print_in=False):
		try:
			if (isinstance(x,np.ndarray)):
				self.x=x
			else:
				self.x=np.array(x)
			if (isinstance(y,np.ndarray)):
				self.y=y
			else:
				self.y=np.array(y)
			# Set default parameters in calculation
			self.SetAllParameters();
			# Set number of group, whether show intermediate information
			self.num_group=num_group
			self.print_in=print_in
			# Check and adjust the data.
			self.CheckData();

			self.dgroup=np.zeros((self.dlen,1),dtype='int32')
			self.centers=np.zeros((self.num_group,self.plen),dtype='float64')
			if (regroup):
				self.regroup=True;
			else: 
				self.regroup=False;
			if (randomgroup):
				self.RandomGroup();
			else:
				self.AutoGroup();

		except:
			traceback.print_exc()
			exit()

	#   Initialize based on data file (with y x data)
	def initializeByFile(self,filename='data.txt',num_group=1,regroup=True,randomgroup=False,print_in=False):
		#load data from file
		try:
			data=np.loadtxt(filename) #dlen(m) data * plen(n) para + 1
			dlen=data.shape[0]
			plen=data.shape[1]
			assert (dlen>1 and plen>1)
			ymat=data[:,0];
			xmat=data[:,1:]
			self.initializeByData(xmat,ymat,num_group,regroup,randomgroup,print_in)
		except:
			traceback.print_exc()
			exit()

	#   Initialize based on data file with group information (with group y x data)
	def initializeByFileWithGroup(self,filename='data.txt',num_group=1,regroup=True,print_in=False):
		#load data from file
		try:
			data=np.loadtxt(filename); #dlen(m) data * plen(n) para + 1 y + 1 group
			dlen=data.shape[0];
			plen=data.shape[1];
			assert (dlen>1 and plen>2)
			gmat=data[:,0];
			ymat=data[:,1];
			xmat=data[:,2:];
			# Use random group method (True here) to initial may be faster
			self.initializeByData(xmat,ymat,num_group,regroup,True,print_in)
			self.SetGroupData(gmat);
		except:
			traceback.print_exc();
			exit();

	def LeastSqaures(self,x,y,rlamda=500):
		try:
			if (x.ndim==1):
				mx=x.size;
				nx=1
			elif (x.ndim==2):
				mx, nx=x.shape;
			else:
				raise ValueError, "x dimension > 2!";
			if (y.ndim==1):
				my=y.size;
				ny=1
			else:
				raise ValueError, "y dimension > 2!";

			if (mx==my and ny==1):
				wts=np.zeros((1,nx))
				wts[0, :]=(np.linalg.inv(x.T.dot(x)+rlamda*np.eye(nx)).dot(x.T.dot(y))).T
				return wts
			elif( mx==ny and my==1):
				wts=np.zeros((1,nx))
				wts[0, :]=(np.linalg.inv(x.T.dot(x)+rlamda*np.eye(nx)).dot(x.T.dot(y.T))).T
			else:
				raise ValueError, "Not match array size for x and y!";
		except:
			traceback.print_exc()
			exit()	

	def CalcRMSArray(self,arr):
		return np.sqrt((arr.T.dot(arr)).ravel()[0]/arr.shape[0])

	#on test
	def FindMatchPoint(self,num):
		err1=self.errors[num];
		gn1=self.dgroup[num];
		y1=self.y[num];
		x1=self.x[num];
		tmpwts=np.zeros((1,self.plen));
		rmsmin=err1;
		nummin=num;
		err2min=0.0;
		absmin=np.array((0.0,0.0))
		gdataIndex=np.where(self.partition_groups.ravel()==gn1)[0]
		#for i in range(self.dlen):
		for j in range(len(gdataIndex)):
			i=gdataIndex[j]
			if (i!=num):
				err2=self.errors[i];
				gn2=self.dgroup[i];
				y2=self.y[i];
				x2=self.x[i];
				#rmsbefore=self.CalcRMSArray(np.array((err1,err2)));
				tmpx=np.array((x1,x2));
				tmpy=np.array((y1,y2));
				tmpwts[0, :]=self.LeastSqaures(tmpx,tmpy,self.rlamda);
				absnow=(np.abs(tmpy.T-(tmpwts[0, :].dot(tmpx.T)).T)).T;
				rmsnow=self.CalcRMSArray(absnow);
				if ( rmsnow<rmsmin): #and rmsnow < rmsbefore):
					print "find better: "+str(absnow)+str(rmsnow)
					err2min=err2;
					absmin=absnow;
					nummin=i;
					rmsmin=rmsnow;
		if (nummin==num):
			print "Not found a better point!"
			return num;
		else:
			print "Better point found: before: "+str(err1)+","+str(err2min)+"; now: "+str(absmin[0])+","+str(absmin[1])
			return nummin;

	def UpdateCenters(self):
		updateHere=False # Whether need to update rms
		for i in range(self.num_group):
			try:
				gdataIndex=np.where(self.partition_groups.ravel()==i+1)[0]
				#print gdataIndex
				if (len(gdataIndex)>=2):
					self.centers[i, :]=sum(self.x.take(gdataIndex,axis=0))/(len(gdataIndex)); 
				elif (len(gdataIndex)==1 and self.regroup):
					maxerrn=np.argmax(self.errors)
					self.errors[maxerrn]=0.0
					tmplist=gdataIndex.tolist()
					tmplist.append(maxerrn)
					self.centers[i, :]=sum(self.x.take(tmplist,axis=0))/2;
					# Whether need to update rms
					updateHere = True
				elif (len(gdataIndex)==0):
					if (self.regroup):
						print "Group disappear: "+ str(i+1)
						if (self.print_in): print "Group disappear: "+ str(i+1)
						maxerrn=np.argmax(self.errors)
						# print self.errors[maxerrn]
						self.errors[maxerrn]=0.0
						maxerrn2=np.argmax(self.errors)
						#print self.errors[maxerrn2]
						self.errors[maxerrn2]=0.0

						# maxerrn2=self.FindMatchPoint(maxerrn)
						# self.errors[maxerrn]=0.0
						# self.errors[maxerrn2]=0.0

						self.centers[i, :]=sum(self.x.take((maxerrn,maxerrn2),axis=0))/2;
						self.partition_groups[maxerrn]=i+1
						self.partition_groups[maxerrn2]=i+1
						# Whether need to update rms
						updateHere=True
					else:
						self.centers[i, :]=np.nan
			# no data in a group
			except ZeroDivisionError:
				print "Still error at centers calculation.."
		return updateHere

	def UpdateRegression(self):
		for i in range(self.num_group):
			#500*eye(n) is the regularization in the regression
			dgindex=np.where(self.partition_groups.ravel()==i+1)[0]
			tmpx=self.x.take(dgindex,axis=0)
			tmpy=self.y.take(dgindex,axis=0)
			# Regression
			#wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;
			self.wts[i, :]=self.LeastSqaures(tmpx,tmpy,self.rlamda); 

			#Regression error for each data
			#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
			# Instead method:
			absnow=(np.abs(tmpy.T-(self.wts[i, :].dot(tmpx.T)).T)).T
			for j in range(len(absnow)):
				#if (self.errors[dgindex[j]]==0.0):print absnow[j]
				self.errors[dgindex[j]]=absnow[j];
		rms=self.CalcRMSArray(self.errors);
		return rms;

	# The main kpane function here
	def kplane(self):
		rms1=10.0; rms2=0.0;
		coef=0.0;

		x=self.x;
		y=self.y;
		m=self.dlen;
		n=self.plen;

		# grouping information during calc.
		dgroup=self.dgroup.copy()
		# Saving partition group result
		self.partition_groups=dgroup;
		# Initialize the centers of the groups
		self.centers=np.zeros((self.num_group,self.plen),dtype='float64')
		# Error for each data
		self.errors=np.zeros((self.dlen,1),dtype='float64')
		# Weights for each group
		self.wts=np.zeros((self.num_group, self.plen));

		# Calculate the weights and errors for each group
		wts=self.wts
		err=self.errors;
		err_temp=np.zeros((self.num_group, 1));
		rms1=self.UpdateRegression();
		# Averge x center value for each group
		centers=self.centers;
		self.UpdateCenters();
		# Update rms if using regroup method
		if(self.regroup):rms1=self.UpdateRegression();

		if self.print_in: print "Starting RMS: " + str(rms1)

		minrms=9999.0;
		minrmscount=0;
		minstartcount=0;
		oldminrms=9999.0;
		# Update the partition and regressor via EM algorithm
		while( np.abs(rms1-rms2)>self.rmslimit or rms1<rms2 or \
			(np.abs(rms1-rms2)<=self.rmslimit and self.regroup )):

			rms1=rms2;
			#Update the partition
			for i in range(m):
				for j in range(self.num_group):
				#lambda*norm2(x(i, :)-centers(j, :)) is the regularization term
					if (self.regroup):
						err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)+self.klamda*np.linalg.norm(x[i, :]-centers[j, :]);
					else:
						err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)
						# discard contribution from regularization term because NaN errors.
				dgroup[i]=np.argmin(err_temp)+1;

			# Update regressor
			rms2=self.UpdateRegression();
			# Update the centers of each group
			rmsUpdate2=self.UpdateCenters();	
			# Update rms if regroup in update center
			if(rmsUpdate2):
				rms2=self.UpdateRegression();

			if self.print_in: print "Grouped state: "+ str(dgroup.T);
			if self.print_in: print "RMS now: "+str(rms2);

			if (rms2<minrms):
				if (minrms-rms2<self.rmslimit):
					minrmscount+=1;
				else: minrmscount=1;
				minrms=rms2;
				minstartcount=0;
				if self.print_in: print "New Min RMS: "+str(minrms);
			if (rms2-minrms<self.rmslimit and rms2>=minrms):
				minrmscount+=1;
				#if (minrmscount>10):
				#	print "Minrms: "+str(minrms)+" for "+str(minrmscount)
			if (minrmscount>10):
				break
			# repeat doing
			if (minstartcount>500):
				oldminrms=minrms;
				minrms=rms2;
				if self.print_in: print "Best Min RMS before changing MINRMS: "+str(oldminrms);
				minstartcount=0;
			minstartcount+=1;

		partition_groups=dgroup.copy();
		coef=wts.copy();
		finalgroup=list(set(dgroup.T[0].tolist()));
		finalgroup.sort();
		if self.print_in: print "Final group number: "+str(len(finalgroup));
		#print rms1,rms2
		#print coef
		return (rms2,coef,partition_groups);

	# Basic file name
	# OnlyGroup to output only the group id as data sequence 
	def outnow(self,filename,OnlyGroup=False):
		fnamelist=os.path.splitext(filename);
		fdata=fnamelist[0]+'_gdata'+fnamelist[1];
		fpar=fnamelist[0]+'_par'+fnamelist[1];
		data=np.zeros((self.dlen,self.plen+2));
		gmat=self.partition_groups.ravel()

		data[:,0]=gmat;
		if (OnlyGroup):
			data[:,1]=self.y;
			data[:,2:]=self.x;
		np.savetxt(fdata,data);
		np.savetxt(fpar,self.wts);

if (__name__=="__main__"):
	try:
		# create Kplane object
		kp=KPlaneObj();
		# Initialize data by file and set up program
		kp.initializeByFile("data.txt",num_group=8,regroup=False,randomgroup=False,print_in=False);
		# Setup parameters. PS: initial will reset all default parameters, 
		# you have to set them after initialization if you don't want to use default. 
		# Default: klamda=0.001, rlamda=500,rmslimit=0.0001
		kp.SetAllParameters(klamda=0.0001,rlamda=500, rmslimit=0.0001);
		#Calculate the partitioning by k-plane regression
		rms, coef, partition_groups = kp.kplane();
		kp.outnow('data.txt');
		print rms
		print coef
		print partition_groups

	except:
		traceback.print_exc()

#raw_input()