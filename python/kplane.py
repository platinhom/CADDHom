#! /usr/bin/env python
# -*- coding: utf8 -*-
'''
#######################################################################################
This program implements the K-plane regression algorithm based on EM algorithm 
to give the optimal partitioning of the data set in the sense ofd minimize the RMS

error of the regression. Together with the SVM classifer, we can do the classification.
rms is the total rms error
coef is the regression coefficients for each plane
x is the training set; x(1:num_samples, 1:number_fearures)
y is the target value: y(1:num_samples)
num_group is the number of groups
dgroup is the initial partitioning
lamda is the regularization constant

Author: Bao Wang, wangbaonj@gmail.com, 09/30/2015 based on Matlab
Reviser: Zhixiong Zhao, 10/06/2015 based on python
########################################################################################
'''

import os,sys,traceback
import numpy as np

#The number of groups
num_group=8;

class KPlaneObj():
	def __init__(self):
		pass

	def CheckData(self,data=True,group=True):
		# Check input data
		if (data):
			if (self.x.ndim==1):
				self.x=self.x.T
				self.dlen=self.x.size;
				self.plen=1
			elif (self.x.ndim==2):
				self.dlen, self.plen=self.x.shape;
			else:
				raise ValueError, "x dimension > 2!";
			if (self.y.ndim==1):
				my=self.y.size;
				if (my!=self.dlen):
					raise ValueError, "data size not equal in x and y!";
				self.y=self.y.T
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

	def SetGroupNumber(self,number=1):
		self.num_group=number;
		# May need further dgroup assignment

	# grouping function can use only after initialize
	# Grouping Based on sequence
	def AutoGroup(self):
		num_group=self.num_group
		dlen=self.dlen

		for k in range(0,num_group-1):
			self.dgroup[k*np.floor(dlen/num_group):(k+1)*np.floor(dlen/num_group)]=k+1;
		self.dgroup[(num_group-1)*np.floor(dlen/num_group):dlen]=num_group;
	# Grouping Based on randomize
	def RandomGroup(self):
		for i in range(self.dlen):
			self.dgroup[i]=np.random.randint(1,num_group+1);

	# Initialize based on x and y data
	def initializeByData(self,x,y,num_group=1,lamda=0.001,print_in=False,randomgroup=False):
		try:
			if (isinstance(x,np.ndarray)):
				self.x=x
			else:
				self.x=np.array(x)
			if (isinstance(y,np.ndarray)):
				self.y=y
			else:
				self.y=np.array(y)
			self.num_group=num_group
			self.lamda=lamda
			self.print_in=print_in
			self.CheckData();
			self.dgroup=np.zeros((self.dlen,1),dtype='int32')
			self.centers=np.zeros((self.num_group,self.plen),dtype='float64')
			if (randomgroup):
				self.RandomGroup();
			else:
				self.AutoGroup();
		except:
			traceback.print_exc()
			#exit()

	# Initialize based on data file
	def initializeByFile(self,filename='data.txt',num_group=1,lamda=0.001,print_in=False,randomgroup=False):
		#load data from file
		try:
			data=np.loadtxt(filename) #dlen(m) data * plen(n) para + 1
			dlen=data.shape[0]
			plen=data.shape[1]
			assert (dlen>1 and plen>1)
			ymat=data[:,0];
			xmat=data[:,1:]
			self.initializeByData(xmat,ymat,num_group,lamda,print_in,randomgroup)
		except:
			traceback.print_exc()
			exit()

	def LeastSqaures(self,x,y,lamda=500):
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
			elif (y.ndim==2):
				my, ny=y.shape;
			else:
				raise ValueError, "y dimension > 2!";

			if (mx==my and ny==1):
				wts=np.zeros((1,nx))
				wts[0, :]=(np.linalg.inv(x.T.dot(x)+lamda*np.eye(nx)).dot(x.T.dot(y))).T
				return wts
			elif( mx==ny and my==1):
				wts=np.zeros((1,nx))
				wts[0, :]=(np.linalg.inv(x.T.dot(x)+lamda*np.eye(nx)).dot(x.T.dot(y.T))).T
			else:
				raise ValueError, "Not match array size for x and y!";
		except:
			traceback.print_exc()
			exit()	

	def CalcRMSArray(self,arr):
		return np.sqrt((arr.T.dot(arr)).ravel()[0]/arr.shape[0])

	def UpdateCenters(self,ReGroup=False):
		updateHere=False # Whether need to update rms
		for i in range(self.num_group):
			try:
				gdataIndex=np.where(self.partition_groups.ravel()==i+1)[0]
				#print gdataIndex
				if (len(gdataIndex)>=2):
					self.centers[i, :]=sum(self.x.take(gdataIndex,axis=0))/(len(gdataIndex)); 
				elif (len(gdataIndex)==1 and ReGroup):
					maxerrn=np.argmax(self.errors)
					self.errors[maxerrn]=0.0
					tmplist=gdataIndex.tolist()
					tmplist.append(maxerrn)
					self.centers[i, :]=sum(self.x.take(tmplist,axis=0))/2;
					# Whether need to update rms
					updateHere = True
				elif (len(gdataIndex)==0):
					if (ReGroup):
						if (self.print_in): print "Group disappear: "+ str(i+1)
						maxerrn=np.argmax(self.errors)
						self.errors[maxerrn]=0.0
						maxerrn2=np.argmax(self.errors)
						self.errors[maxerrn2]=0.0
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

	def UpdateRegression(self,lamda=500):
		for i in range(self.num_group):
			#500*eye(n) is the regularization in the regression
			dgindex=np.where(self.partition_groups.ravel()==i+1)[0]
			tmpx=self.x.take(dgindex,axis=0)
			tmpy=self.y.take(dgindex,axis=0)
			# Regression
			#wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;
			self.wts[i, :]=self.LeastSqaures(tmpx,tmpy,lamda); 

			#Regression error for each data
			#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
			# Instead method:
			absnow=(np.abs(tmpy.T-(self.wts[i, :].dot(tmpx.T)).T)).T
			for j in range(len(absnow)):
				self.errors[dgindex[j]]=absnow[j];
		rms=self.CalcRMSArray(self.errors);
		return rms;

	def kplane(self,ReGroup=False):
		rms1=10.0; rms2=0.0;
		coef=0.0;
		LeastSqauresLamda=500

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
		rms1=self.UpdateRegression(LeastSqauresLamda);
		# Averge x center value for each group
		centers=self.centers;
		self.UpdateCenters(ReGroup);
		# Update rms if using regroup method
		if(ReGroup):rms1=self.UpdateRegression(LeastSqauresLamda);

		if self.print_in: print "Starting RMS: " + str(rms1)

		minrms=9999.0
		minrmscount=0
		minstartcount=0
		oldminrms=9999.0
		# Update the partition and regressor via EM algorithm
		while( np.abs(rms1-rms2)>0.0001 or rms1<rms2 or (np.abs(rms1-rms2)<=0.0001 and ReGroup )):
			rms1=rms2;
			#Update the partition
			for i in range(m):
				for j in range(self.num_group):
				#lambda*norm2(x(i, :)-centers(j, :)) is the regularization term
					if (ReGroup):
						err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)+self.lamda*np.linalg.norm(x[i, :]-centers[j, :]);
					else:
						err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)
						# discard contribution from regularization term because NaN errors.
				dgroup[i]=np.argmin(err_temp)+1;

			# Update regressor
			rms2=self.UpdateRegression(LeastSqauresLamda);
			# Update the centers of each group
			rmsUpdate2=self.UpdateCenters(ReGroup);	
			# Update rms if regroup in update center
			if(rmsUpdate2):
				rms2=self.UpdateRegression(LeastSqauresLamda);

			if self.print_in: print "Grouped state: "+ str(dgroup.T)
			if self.print_in: print "RMS now: "+str(rms2)
			print "RMS&Group: "+str(rms2)+ str(dgroup.T)

			if (rms2<minrms):
				if (minrms-rms2<0.0001):
					minrmscount+=1;
				else: minrmscount=1;
				minrms=rms2;
				minstartcount=0;
				print "New Min RMS: "+str(minrms)
			if (rms2-minrms<0.0001 and rms2>=minrms):
				minrmscount+=1;
				#if (minrmscount>10):
				#	print "Minrms: "+str(minrms)+" for "+str(minrmscount)
			if (minrmscount>10):
				break
			# repeat doing
			if (minstartcount>1000):
				oldminrms=minrms;
				minrms=rms2;
				print "Best Min RMS before changing MINRMS: "+str(oldminrms)
				minstartcount=0;
			minstartcount+=1;

		partition_groups=dgroup.copy()
		coef=wts.copy()
		finalgroup=list(set(dgroup.T[0].tolist()))
		finalgroup.sort()
		if self.print_in: print "Final group number: "+str(len(finalgroup))
		#print rms1,rms2
		#print coef
		return (rms2,coef,partition_groups);

if (__name__=="__main__"):
	try:
		# create Kplane object and initialize
		kp=KPlaneObj();
		kp.initializeByFile("data.txt",num_group=8,lamda=0.0,print_in=False,randomgroup=True);

		#Update the partitioning by k-plane regression
		rms, coef, partition_groups=kp.kplane(ReGroup=True);

	except:
		traceback.print_exc()

#raw_input()