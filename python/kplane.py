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
lambda is the regularization constant

Author: Bao Wang, wangbaonj@gmail.com, 09/30/2015 based on Matlab
Reviser: Zhixiong Zhao, 10/06/2015 based on python
########################################################################################
'''

import os,sys,traceback
import numpy as np

#The number of groups
num_group=8;

class KPlaneObj():
	def __init__(self,x,y,dgroup,num_group=3,Clamda=0.001,print_in=False):
		if (isinstance(x,np.ndarray)):
			self.x=x
		else:
			self.x=np.array(x)
		if (isinstance(y,np.ndarray)):
			self.y=y
		else:
			self.y=np.array(y)
		if (isinstance(dgroup,np.ndarray)):
			self.dgroup=dgroup
		else:
			self.dgroup=np.array(dgroup)
		self.num_group=num_group
		self.Clamda=Clamda
		self.print_in=print_in

def kplane(x,y,num_group,dgroup,Clamda,print_in=False):
	rms1=10.0; rms2=0.0;
	rms=0.0;
	coef=0.0;
	m, n=x.shape;
	# Initialize the centers of the groups
	centers=np.zeros((num_group, n));

	err=np.zeros((m, 1));
	err_temp=np.zeros((num_group, 1));
	
	# averge x value for each group
	for i in range(num_group):
		try:
			centers[i, :]=sum(x.take(np.where(dgroup.ravel()==i+1)[0],axis=0))/(len(np.where(dgroup.ravel()==i+1)[0])); 
		# no data in a group
		except ZeroDivisionError:
			centers[i, :]=nan
	#Calculate the weights for each group
	wts=np.zeros((num_group, n));
	for i in range(num_group):
		#500*eye(n) is the regularization in the regression
		dgindex=np.where(dgroup.ravel()==i+1)[0]
		tmpx=x.take(dgindex,axis=0)
		tmpy=y.take(dgindex,axis=0)
		# Regression
		wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;
		#Regression error for each data
		#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
		# Instead method:
		absnow=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T
		for j in range(len(absnow)):
			err[dgindex[j]]=absnow[j];
	rms1=np.sqrt((err.T.dot(err).ravel()[0])/m);

	if print_in: print "Starting RMS: " + str(rms1)

	# Update the partition and regressor via EM algorithm
	while(np.abs(rms1-rms2)>0.001):
		rms1=rms2;
		#Update the partition
		for i in range(m):
			for j in range(num_group):
			#lambda*norm2(x(i, :)-centers(j, :)) is the regularization term
				err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)

				# discard contribution from average x value
				#+Clambda*np.linalg.norm(x[i, :]-centers[j, :]);

			dgroup[i]=np.argmin(err_temp)+1;

		#Update regressor
		for i in range(num_group):
			#500*eye(n) is the regularization in the regression
			dgindex=np.where(dgroup.ravel()==i+1)[0]
			tmpx=x.take(dgindex,axis=0)
			tmpy=y.take(dgindex,axis=0)

			# update weigh until no data in the group
			if (len(tmpx)!=0):
				wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;

			#Updated Regression error for each data
			#	Matlab method fail in python
			#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
			# 	Instead method:
			absnow=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T
			for j in range(len(absnow)):
				err[dgindex[j]]=absnow[j];

		#Update the centers of each group
		for i in range(num_group):
			try:
				centers[i, :]=sum(x.take(np.where(dgroup.ravel()==i+1)[0],axis=0))/(len(np.where(dgroup.ravel()==i+1)[0]));
			# No data in a group
			except ZeroDivisionError:
				#print "ZeroDivisionError at "+ str(i)
				centers[i, :]=np.nan

		rms2=np.sqrt((err.T.dot(err).ravel()[0])/m);
		if print_in: print "Grouped state: "+ str(dgroup.T)
		if print_in: print "RMS now: "+str(rms2)

	partition_groups=dgroup
	coef=wts
	rms=rms2
	finalgroup=list(set(dgroup.T[0].tolist()))
	finalgroup.sort()
	if print_in: print "Final group number: "+str(len(finalgroup))
	#print rms1,rms2
	#print coef
	return (rms,coef,partition_groups);

def kplane2(x,y,num_group,dgroup,Clamda,print_in=False):
	rms1=10.0; rms2=0.0;
	rms=0.0;
	coef=0.0;
	m, n=x.shape;

	#Calculate the weights for each group
	err=np.zeros((m, 1));
	err_temp=np.zeros((num_group, 1));
	wts=np.zeros((num_group, n));
	for i in range(num_group):
		#500*eye(n) is the regularization in the regression
		dgindex=np.where(dgroup.ravel()==i+1)[0]
		tmpx=x.take(dgindex,axis=0)
		tmpy=y.take(dgindex,axis=0)
		# Regression
		wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;
		#Regression error for each data
		#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
		# Instead method:
		absnow=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T
		for j in range(len(absnow)):
			err[dgindex[j]]=absnow[j];

	# Initialize the centers of the groups
	centers=np.zeros((num_group, n));
	# averge x value for each group
	for i in range(num_group):
		try:
			centers[i, :]=sum(x.take(np.where(dgroup.ravel()==i+1)[0],axis=0))/(len(np.where(dgroup.ravel()==i+1)[0])); 
		# no data in a group
		except ZeroDivisionError:
			print "Group disappear: "+ str(i+1)
			maxerrn=np.argmax(err)
			err[maxerrn]=0.0
			maxerrn2=np.argmax(err)
			err[maxerrn2]=0.0
			centers[i, :]=sum(x.take((maxerrn,maxerrn2),axis=0))/2;
			dgroup[maxerrn]=i+1
			dgroup[maxerrn2]=i+1

	rms1=np.sqrt((err.T.dot(err).ravel()[0])/m);
	if print_in: print "Starting RMS: " + str(rms1)

	# Update the partition and regressor via EM algorithm
	while(np.abs(rms1-rms2)>0.00001 or rms1<rms2):
		rms1=rms2;
		#Update the partition
		for i in range(m):
			for j in range(num_group):
			#lambda*norm2(x(i, :)-centers(j, :)) is the regularization term
				err_temp[j]=abs(y[i]-(wts[j, :].dot(x[i, :].T)).T)+Clambda*np.linalg.norm(x[i, :]-centers[j, :]);

			dgroup[i]=np.argmin(err_temp)+1;

		#Update regressor
		for i in range(num_group):
			#500*eye(n) is the regularization in the regression
			dgindex=np.where(dgroup.ravel()==i+1)[0]
			tmpx=x.take(dgindex,axis=0)
			tmpy=y.take(dgindex,axis=0)

			# update weigh until no data in the group
			if (len(tmpx)!=0):
				wts[i, :]=(np.linalg.inv( tmpx.T.dot(tmpx)+500*np.eye(n) ).dot(tmpx.T.dot(tmpy)) ).T;

			#Updated Regression error for each data
			#	Matlab method fail in python
			#err.take(dgindex,axis=0)=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T;
			# 	Instead method:
			absnow=(np.abs(tmpy.T-(wts[i, :].dot(tmpx.T)).T)).T
			for j in range(len(absnow)):
				err[dgindex[j]]=absnow[j];

		#Update the centers of each group
		for i in range(num_group):
			try:
				centers[i, :]=sum(x.take(np.where(dgroup.ravel()==i+1)[0],axis=0))/(len(np.where(dgroup.ravel()==i+1)[0]));
			# No data in a group
			except ZeroDivisionError:
				print "Group disappear: "+ str(i+1)
				# Move max err point on a line
				maxerrn=np.argmax(err)
				err[maxerrn]=0.0
				maxerrn2=np.argmax(err)
				err[maxerrn2]=0.0
				centers[i, :]=sum(x.take((maxerrn,maxerrn2),axis=0))/2;
				dgroup[maxerrn]=i+1
				dgroup[maxerrn2]=i+1

		rms2=np.sqrt((err.T.dot(err).ravel()[0])/m);
		if print_in: print "Grouped state: "+ str(dgroup.T)
		if print_in: print "RMS now: "+str(rms2)

	partition_groups=dgroup
	coef=wts
	rms=rms2
	finalgroup=list(set(dgroup.T[0].tolist()))
	finalgroup.sort()
	if print_in: print "Final group number: "+str(len(finalgroup))
	#print rms1,rms2
	#print coef
	return (rms,coef,partition_groups);


if (__name__=="__main__"):
	try:

		#load data
		data=np.loadtxt('data.txt') #m data*dlen(n) para
		dlen=data.shape[0]
		plen=data.shape[1]

		ymat=data[:,0];
		ymat.shape=(dlen,1)
		xmat=data[:,1:]
		dgroup=np.zeros((dlen,1),dtype='int32')

		#regressors
		coef=np.zeros((num_group, plen-1)); 
		#The partitions of the groups
		partition_groups=np.zeros((dlen,1));
		#Regularization constant in the kplane regression
		Clambda=0.01;

		for k in range(0,num_group-1):
			dgroup[k*np.floor(dlen/num_group):(k+1)*np.floor(dlen/num_group)]=k+1;
		dgroup[(num_group-1)*np.floor(dlen/num_group):dlen]=num_group;

		#Update the partitioning by k-plane regression
		rms, coef, partition_groups=kplane2(xmat, ymat, num_group, dgroup, Clambda,True)

	except:
		traceback.print_exc()

#raw_input()