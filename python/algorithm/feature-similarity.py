#! /usr/bin/env python
import sys,heapq,string,traceback
import numpy as np
import Data2D as D2
from optparse import OptionParser


####### Class defined here ############




####### Function defined here ############

def readfile(self,filename,delim=None,convert=str):
	"""Read data from file. 
	delim is the delimiter for data.
	convert is the method to convert the data."""
	if (not filename):
		raise IOError("No file was found!");
	f=open(filename);
	datas=[]
	if (convert==str):
		for line in f:
			data=line.strip().split(delim)
			datas.append(data)
	else:
		for line in f:
			data=line.strip().split(delim)
			datas.append(map(convert,data))
	f.close()
	return datas

def similarity(vals1,vals2,maxs,mins,weights,coeff=2):
	"""Calculate the similarity between two data;
	sim=weight[i]*(1-|vals1[i]-vals2[i]|/(maxs[i]-mins[i]))^coeff 
	vals1 and vals2 are the features for two data; 
	maxs and mins are the max/min feature in the data group; 
	weights is the weight list for the features
	coeff is for the power function"""
	sv=0.0
	for i in range(len(vals1)):
		if ((maxs[i]-mins[i])==0): 
			#sv += weights[i];
			continue;
		sv+=weights[i]*pow(1-(abs(vals1[i]-vals2[i])/float(maxs[i]-mins[i])), coeff)
	return sv

def outputmatrix(filename,datas):
	"""Write out matrix"""

	# Check data
	if (len(datas)<1 or len(datas[0])<1 ):
		raise ValueError("Error matrix data..")
		exit(1)
	
	# Check data type
	if (isinstance(datas[0][0],float) and datas[0][0]<10 and datas[0][0]>=0):
		formatlambda=(lambda val: format(val, '6.4f'));
	elif (isinstance(datas[0][0],int) and datas[0][0]<1000000 and datas[0][0]>-100000):
		formatlambda=(lambda val: format(val, '6d'));
	else:
		formatlambda=(lambda val: str(val));

	# write out matrix
	fw=open(filename,'w')
	for i in range(len(datas)):
		fw.write(" ".join(map(formatlambda,(datas[i])))+"\n")
	fw.close()


def mostsimilarMOL(filename=None, weighs=None, 
	neednum=1, datas=None,outfilesimmatrix=None):
	""" Calculate similarity matrix and top #neednum similar mol index."""

	if (not filename and (not datas)):
		raise ValueError("Error blank file name or blank data")
		exit(1)

	# Read data
	elif (filename and (not datas)):
		fr=open(filename)
		datas=[]
		for line in fr:
			tmp=line.strip().split()
			datas.append(map(float,tmp))
		fr.close()

	# Read file as use datas as index
	elif(filename and datas):
		fr=open(filename)
		needindex=map(int,datas);
		datas=[]
		for line in fr:
			tmp=map(float, line.strip().split())
			tmp2=[]
			for cindex in needindex:
				tmp2.append(tmp[cindex-1]);
			datas.append(tmp2);
		fr.close()		
	# Use giving datas
	else:
		datas=map(float,datas)

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
	if (not weighs or len(weighs)<=0):
		print "No weights are given. Set all weights to 1."
		weighs=[1]*numfeature
	elif (numfeature != len(weighs)):
		raise ValueError("Length not the same for features ("+str(numfeature)+") and giving weights")
		exit(1)
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

	if (outfilesimmatrix): outputmatrix(outfilesimmatrix,simmatrix)

	# find out top similar id, start from 1
	dlen='<'+str(len(str(numdata))); # < used for left alignment
	formatlambda=(lambda val: format(val, dlen+'d'));
	simmols=[]
	for i in range(numdata):
		largestvalue=heapq.nlargest(neednum+1, simmatrix[i])
		minlim=largestvalue[-1]
		outnums=[];
		for j in range(numdata):
			if (simmatrix[i][j] >= minlim and i != j): outnums.append(j+1)
		simmols.append(outnums);
		#print " ".join(map(formatlambda,outnums))
	return (simmatrix, simmols)

def readtraindatas(filename):
	"""Read data for training. """
	if (not filename):
		raise IOError("No train file was found!");
	f=open(filename)
	datas=[]
	for line in f:
		data=line.strip().split()
		datas.append(data)
	f.close()
	return datas

def traindatas(datas):
	"""Train data to get weights.
	Input as: [[d11, d12, d13..], [d21, d22, d23..], .."""
	weights=[]

	return weights

def predictdata(data, weights):
	'''Predict value for giving data with given weights'''
	pvalue=0;
	return pvalue;

def errordata(pvals, rvals):
	"""Errors between giving predicted values (list) and real values (list)"""
	return 0;

if (__name__=="__main__"):
	helpdes='''Calculate features of molecules based on Pybel and Openbabel.
	# For one file, use -i option to assign the input file;
	# For many files, use -m option to assign a file containing file name without extension.
	# -f option can assign the file format. It must be given when using -m option. 
	# Without -f option and using -i option, the format will be deduced based on file extension.
	# -t option will print the title for features.'''

	parser = OptionParser(description=helpdes) 
	parser.add_option("-i", "--input", action="store", 
					dest="input", default="",
					help="Read input data from input file")
	parser.add_option("-n", "--neednum", action="store", 
					dest="neednum", default="",
					help="Output given amount of most similar molecules")
	parser.add_option("-w", "--weighs", action="store", 
					dest="weighs", default="",
					help="Giving all weights, as -w 1,2,3. Default blank to set all to 1.")
	parser.add_option("-d","--dindex", action="store", 
					dest="dataindex", default="",
					help="Giving data index for features, as -d 1,2,3. Default blank to use all features.")
	parser.add_option("-s", "--simout", action="store", 
					dest="simmatrix", default="",
					help="Output file for similarity matrix.")
	parser.add_option("-o", "--output", action="store", 
					dest="output", default="",
					help="The output file to save result")
	parser.add_option("-t", "--train", action="store", 
					dest="train", default="",
					help="Training set file")
	parser.add_option("-r", "--real", action="store", 
					dest="real", default="",
					help="File saving real value for data")
	(options, args) = parser.parse_args()
	
	if (len(sys.argv)<2):
		print "Please assign an input file or a file containing all file prefix!"
		parser.print_help()
		#parser.print_description()
		#parser.print_usage()
		exit(1)

	if (len(sys.argv)==2):
		options.input=sys.argv[1];

	neednum=1
	if (int(options.neednum)):
		neednum=int(options.neednum);

	weighs=None
	if (options.weighs):
		tmp=options.weighs.split(",")
		if (len(tmp)>1):
			weighs=[ int(string.strip(i)) for i in tmp]

	dataindex=None
	if options.dataindex:
		tmp=options.dataindex.split(",")
		if (len(tmp)>1):
			dataindex=[ int(string.strip(i)) for i in tmp]

	outfilesimmatrix=""	
	if (options.simmatrix):
		outfilesimmatrix=options.simmatrix

	trainfile=""	
	if (options.train):
		trainfile=options.train

	realfile=""	
	if (options.real):
		realfile=options.real

	simmatrix,simmols=mostsimilarMOL(options.input,weighs=weighs, neednum=neednum, 
		datas=dataindex,outfilesimmatrix=outfilesimmatrix);

	if (options.output):
		outputmatrix(options.output,simmols);
