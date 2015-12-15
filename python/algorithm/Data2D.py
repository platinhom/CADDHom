'''Module for 2D data based on numpy.ndarray'''

__author__="Zhixiong Zhao"
__date__="2015.12.08"

import string
#import traceback
import numpy as np

### Test Data
a=[[1,2,3],[2,3,4],[3,4,5],[4,5,6]]
b=Data2D(a)
c=Data2D([1,1.2,1.4,1.6]).T

class mlarray(np.ndarray):
	@property
	def I(self):
		''' Inverse of the array.'''
		return np.linalg.inv(self.data)

class Data2D(object):
	"""Class to deal with the 2D data in numpy array."""
	def InitData(self, data, dtype="float64", name='Data'):
		'''Initial data from list to array'''
		if (isinstance(data,list) or isinstance(data,tuple)):
			arr=np.array(data,dtype=dtype);
		elif (isinstance(data,np.ndarray)):
			arr=data
		if arr.ndim >2:
			raise ValueError("Data dimension > 2!");
			exit();
		else:
			# set the object property
			self.name=name;
			self.data=arr;
			# No data but show ndim=1
			if self.data.size==0:
				print "Warning: Blank data!"
				self.data.shape=(0,0);
			# Dimension corect, translate 1D to one row
			# [1,2,3..]
			if self.data.ndim==1:
				self.data.shape=(1,self.data.shape[0]);
			return self.data

	def __init__(self,data=None,
		filename=None, dtype="float64", delim=None,name=None):
		if (filename):
			self.readfile(filename,delim=delim,dtype=dtype);
			if (isinstance(name,str)):
				self.name=name;
		elif (data is not None):
			if (isinstance(data,list) or isinstance(data,tuple) or isinstance(data,np.ndarray)):
				self.InitData(data,name=name,dtype=dtype)
			else:
				#Todo: optimize
				raise TypeError("Error input type")
		else:
			self.data=None;
			self.name=""
		if (name):
			self.name=name
		self.dtype=dtype;

	def __len__(self):
		return self.data.__len__();

	def __str__(self):
		return self.data.__str__();

	def __repr__(self):
		return self.data.__repr__();

	# important for give Data2D to numpy
	def __getitem__(self,i):
		return self.data.__getitem__(i);

	def __delitem__(self,i):
		return self.data.__delitem__(i);

	def __setitem__(self,i,y):
		return self.data.__setitem__(i,y);

	# for slice
	def __getslice__(self,i,j):
		return self.data.__getslice__(i,j);

	def __delslice__(self,i,j):
		return self.data.__delslice__(i,j);

	def __setslice__(self,i,j,y):
		return self.data.__setslice__(i,j,y);

	def __iter__(self):
		return self.data.__iter__();

	def __add__(self,y):
		return self.data.__add__(y);

	def __sub__(self,y):
		return self.data.__sub__(y);


	# ndarray-based
	@property
	def ndim(self):
	    return 2

	@property
	def size(self):
		return self.data.size
	
	@property
	def shape(self):
		return self.data.shape

	@shape.setter
	def shape(self,value):
		self.resize(value);

	@property
	def row(self):
		return self.data.shape[0]
	
	@property
	def col(self):
		return self.data.shape[1]	

	@property
	def T(self):
		'''Transpose of the array.'''
		return Data2D(data=self.data.T)

	@property
	def I(self):
		''' Inverse of the array.'''
		return Data2D(np.linalg.inv(self.data))
	
	@property
	def ADT(self):
		'''array.dot(array.T)'''
		return Data2D(self.data.dot(self.data.T))
	
	@property
	def TDA(self):
		'''array.T.dot(array)'''
		return Data2D(self.data.T.dot(self.data))	


	def resize(self,row=1,col=1):
		if (isinstance(row,tuple) or isinstance(row,list)):
			if (len(row)==2):
				col=row[1];
				row=row[0];
			elif (len(row)==1):
				row=row[0]
			else:
				raise ValueError("Error shape tuple dimension!")
		self.data.resize((row,col));

	def reshape(self,row=1,col=1):
		if (isinstance(row,tuple) or isinstance(row,list)):
			if (len(row)==2):
				col=row[1];
				row=row[0];
			elif (len(row)==1):
				row=row[0]
			else:
				raise ValueError("Error shape tuple dimension!")
		return self.data.reshape((row,col));

	def dot(self, data):
		if (isinstance(data,np.ndarray)):
			return Data2D(self.data.dot(data))
		if (not isinstance(data,Data2D)):
			data=self.convert(data)
		return Data2D(self.data.dot(data.data))

	def tolist(self):
		'''Convert to a list'''
		return self.data.tolist();

	##### Custom Funtions #####
	
	## array related 
	def append(self,data,column=False):
		if (isinstance(data,list) or isinstance(data,tuple)):
			data=np.array(data)
		elif (isinstance(data,np.ndarray) or isinstance(data,Data2D)):
			pass
		else:
			raise TypeError("Type error for data:"+str(type(data)));

		# Perceive number of row/column of given data
		nrow=1;ncol=1;
		if data.ndim is 1:
			if column: 
				nrow=len(data)
				data.shape=(nrow,1)
			else:
				ncol=len(data)
				data.shape=(1,ncol)
		else:
			nrow=data.shape[0]
			ncol=data.shape[1]
		# Check length
		axis=0
		if (column):
			axis=1
			if (self.row!=nrow):
				raise ValueError("Length Data Row "+str(self.row)+"is not equal to given data lengt:"+str(nrow));
		else:
			if (self.col!=ncol):
				raise ValueError("Length Data Column "+str(self.col)+"is not equal to given data lengt:"+str(ncol));

		if data.ndim <= 2:
			self.data=np.append(self.data,data,axis)
			#olddata=self.data
			#if (not column):
			#	self.data=np.zeros((self.row+nrow,self.col),dtype=self.dtype)
			#	self.data[:self.row,:]=olddata
			#	self.data[self.row:,:]=data
			#	self.row+=nrow;			
			#else:
			#	self.data=np.zeros((self.row,self.col+ncol),dtype=self.dtype)
			#	self.data[:,:self.col]=olddata
			#	self.data[:,self.col:]=data	
			#	self.col+=ncol;
			#del olddata
		else:
			raise ValueError("Error dimension! "+str(data.ndim));

	def func4each(self,func,column=False, tolist=False):
		'''Return Data2D/list containing result based on func for each row(default)/column'''
		if (tolist):
			if (not column):
				return [ func(i) for i in self.data]
			else:
				return [ func(i) for i in self.data.T]
		else:
			if (not column):
				return Data2D( data=[ [func(i)] for i in self.data],name=func.__name__)
			else:
				return Data2D(data=[[ func(i) for i in self.data.T]],name=func.__name__)

	def getchild(self, childlist, column=False, name="Child", dtype="float64"):
		"""Creat a child Data2D object
		Given a list for index, starting from zero"""
		# Todo: Optimize algorithm
		if (isinstance(childlist,list) or isinstance(childlist,tuple) ):
			if (not column):
				childdatas=[ self.data[i].tolist() for i in childlist ];
				return Data2D(data=childdatas, name=name, dtype=dtype);
			else:
				tdata=self.data.T
				childdatas=[ tdata[i].tolist() for i in childlist ];
				return Data2D(data=childdatas, name=name, dtype=dtype).T;

		elif (isinstance(childlist,str)):
			clist=[]
			# parse the string, 0,2-4 mean [0], [2], [3],[4]
			# use , as delimiter
			tmp=[ string.strip(i) for i in childlist.strip().split(',') ]
			for i in tmp:
				if (i is ''):
					continue;
				if ('-' in i):
					tmp2= i.split('-')
					if tmp2[0]=='':tmp2[0]=0
					if tmp2[1]=='':tmp2[1]=self.row-1
					if int(tmp2[1])>=self.row:
						tmp2[1]=self.row-1
					## Don't support negative number
					#if int(tmp2[0])<0:
					#	tmp2[0]=self.row+int(tmp2[0])
					#if int(tmp2[1])<0:
					#	tmp2[1]=self.row+int(tmp2[1])
					if (int(tmp2[1])>=int(tmp2[0]) and int(tmp2[0])>=0 and int(tmp2[1])<self.row):
						for j in range(int(tmp2[0]),int(tmp2[1])+1):
							clist.append(j)
					else:
						print tmp2
				else:
					clist.append(int(i))
			clist.sort()
			if (not column):
				return Data2D(data=[ self.data[i].tolist() for i in clist], name=name, dtype=dtype);
			else:
				tdata=self.data.T
				return Data2D(data=[ tdata[i].tolist() for i in clist], name=name, dtype=dtype).T;
		else:
			raise TypeError("Error given child type, should be list/tuple or string")
			return None;

	## User define

	def convert(self,data, dtype="float64", name='Data'):
		if (isinstance(data,list) or isinstance(data,tuple) or isinstance(data,np.ndarray)):
			data=Data2D(data,dtype=dtype,name=name);
		if (isinstance(data,Data2D)):
			return data
		else: 
			raise TypeError("Can't convert to Data2D type!");
			return None

	def readfile(self,filename,delim=None,dtype="float64"):
		"""Read data from file. 
		delim is the delimiter for data.
		convert is the method to convert the data."""
		if (not filename):
			raise IOError("No file was found!");
		# read file
		f=open(filename);
		datas=[]
		for line in f:
			data=line.strip().split(delim)
			datas.append(data)
		f.close()
		# Create New data from list
		self.InitData(datas,dtype=dtype,name=filename);
		del datas[:]
		return self.data

	def savefile(self,filename,delim=" ",outformat=""):
		'''outformat is the format formula in format function'''
		f=open(filename,'w');
		for i in range(self.row):
			f.write(delim.join( \
				[ format(item, outformat) for item in self.data[i].tolist()]) \
			+"\n");
		f.close()

	## test function for filter
	#def testfunc(self,d):
	#	if d[0]>1:return True; 
	#	else: return False
	def filter(self,func, column=False):
		'''Use a filter function to filter the data.
		When func(data[i]) return True, save it. 
		The filter can one by one in row or in column.'''
		copydata=None
		if column:
			copydata=self.T
		else:
			copydata=self
		outdata=np.zeros(copydata.shape)	
		savenum=0;
		for d in copydata:
			if (func(d)):
				outdata[savenum]=d
				savenum+=1;
		finaldata=outdata[:savenum,:]
		if (column):
			return Data2D(finaldata.T, dtype=self.dtype)
		else:
			return Data2D(finaldata, dtype=self.dtype)

	def each4min(self,column=False, tolist=False):
		'''Min for each row/column'''
		out= np.amin(self.data,axis=0) if column else np.amin(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4max(self,column=False, tolist=False):
		'''Max for each row/column'''
		out= np.amax(self.data,axis=0) if column else np.amax(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4sum(self,column=False, tolist=False):
		'''Sum for each row/column'''
		out= np.sum(self.data,axis=0) if column else np.sum(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4median(self,column=False, tolist=False):
		'''Median for each row/column.
		Given a vector V of length N, the median of V is the middle value of a sorted copy of V, 
		V_sorted - i.e., V_sorted[(N-1)/2], when N is odd. 
		When N is even, it is the average of the two middle values of V_sorted.'''
		out= np.median(self.data,axis=0) if column else np.median(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4average(self,weights=None, column=False, tolist=False):
		'''Average for each row/column.
		Different to mean, it can be given a weight for each element!'''
		out= np.average(self.data,weights=weights,axis=0) if column else np.average(self.data,weights=weights,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4mean(self,column=False, tolist=False):
		'''Mean for each row/column'''
		out= np.mean(self.data,axis=0) if column else np.mean(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4std(self,column=False, tolist=False):
		'''Standard deviation for each row/column'''
		out= np.std(self.data,axis=0) if column else np.std(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4var(self,column=False, tolist=False):
		'''Variance for each row/column'''
		out= np.var(self.data,axis=0) if column else np.var(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4ptp(self,column=False, tolist=False):
		'''Max-Min for each row/column'''
		out= np.ptp(self.data,axis=0) if column else np.ptp(self.data,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def each4percentile(self, q, column=False, tolist=False):
		'''(Max-Min)*q/100 for each row/column'''
		out= np.percentile(self.data,q,axis=0) if column else np.percentile(self.data,q,axis=1).reshape(self.data.shape[0],1)
		return  out.ravel().tolist() if tolist else Data2D(out,dtype=self.dtype)

	def LeastSquares(self,dataY,rlamda=500):
		'''Least-Squares method to calculate weight W'''
		mx, nx=self.data.shape;
		# data: (my, 1)
		if (dataY.ndim==1):
			my=dataY.size;
			ny=1
		# data: (my, ny)
		elif (dataY.ndim==2):
			my, ny=dataY.shape;
		else:
			raise ValueError, "data for y dimension > 2!";

		wts=np.zeros((nx,1))
		# data: (mx, nx) * (my, 1)
		if (mx==my and ny==1):
			# (self.TDA+rlamda*np.eye(nx)).I.dot(self.T.dot(Y)).T
			wts=(np.linalg.inv(self.data.T.dot(self.data)+rlamda*np.eye(nx)).dot(self.data.T.dot(dataY)))
		# data: (mx, nx) * (1, my)
		elif( mx==ny and my==1):
			wts=(np.linalg.inv(self.data.T.dot(self.data)+rlamda*np.eye(nx)).dot(self.data.T.dot(dataY.T)))
		else:
			raise ValueError, "Not match array size for x and y!";
		return Data2D(data=wts,name="weight");


#########################################
class LT2D(object):
	"""Linear Transformation for ndarray/Data2D"""
	def __init__(self, X, Y, rlamda=500):
		X=self.normX(X);
		self.X=X;
		Y=self.normY(Y,checkXY=True);
		if (self.checkXY(X,Y)):
			self.Y=Y
			self.rlamda=rlamda
			# W for weights for input X/Y
			self.W=self.LeastSquares(self.X, self.Y, rlamda=self.rlamda)
			# train for training LT2D
			self._train=None
		else:
			raise ValueError("Error for matching X-Y size..")

	@property
	def row(self):
		'''row number of X'''
		return self.X.shape[0]

	@property
	def col(self):
		'''column number of X'''
		return self.X.shape[1]

	@property
	def error(self):
		'''Error array from current Y'=XW and real Y'''
		return self.calcError()

	@property
	def rms(self):
		'''Self RMS'''
		return self.calcRMS()

	@property
	def WT(self):
		'''WT for weights from training set for test set/validation'''
		if (self.train): return self.train.W

	@property
	def train(self):
		# Saving a training LT2D object
		if (not self._train):
			raise ValueError('Train model is not set!')
			return None
		return self._train

	@train.setter
	def train(self,value):
		if (not isinstance(value,LT2D)):
			raise ValueError("Error, argument should be LT2D class")
		self._train=value;

	def normX(self, X):
		'''Normalize X as ndarray'''
		if (not (isinstance(X,Data2D) or isinstance(X,np.ndarray))):
			X=np.array(X)
		return X

	def normY(self, Y, checkXY=False):
		'''Normalize Y as ndarray ~ (n,1) '''
		if (not (isinstance(Y,Data2D) or isinstance(Y,np.ndarray))):
			Y=np.array(Y)
		if (Y.shape[1]!=1 and Y.shape[0]==1):
			Y.resize(Y.shape[1],Y.shape[0]);
		if (Y.shape[1]!=1):
			raise ValueError("Y column is not 1!")
		# When set Y, check whether consistent to X. 
		# Need checkXY=True and set X before.
		if (checkXY and not self.checkXY(X=self.X,Y=Y)):
			raise ValueError("Y row is not equal to X row!")
		return Y

	def checkXY(self, X, Y):
		'''Check X-Y shape consistent'''
		return X.row==Y.row

	###### Reset some value and recalculate

	def resetX(self, X, rlamda=None):
		'''reset X and W'''
		if (rlamda is None): 
			rlamda=self.rlamda;
		self.X=self.normX(X);
		self.W=self.LeastSquares(self.X,self.Y,rlamda=rlamda)

	def resetY(self,Y, rlamda=None):
		'''reset Y and W'''
		if (not rlamda): rlamda=self.rlamda;
		self.Y=self.normY(Y)
		self.W=self.LeastSquares(self.X, self.Y, rlamda=rlamda)

	def resetR(self,rlamda):
		'''recalculate weight based on given rlamda'''
		self.W=self.LeastSquares(self.X, self.Y, rlamda=rlamda)

	##### Calculate Linear Tranformation methods
	def LeastSquares(self,dataX=None,dataY=None,rlamda=500):
			'''Least Sqaures for y data and x data.
			- dataX/Y is ndarray/Data2D object
			- Return the weight for x'''
		#try:
			if (dataX is None):
				dataX=self.X
			if (dataY is None):
				dataY=self.Y
			# data: (mx,1)
			if (dataX.ndim==1):
				mx=dataX.size;
				nx=1
			# data: (mx, nx)
			elif (dataX.ndim==2):
				mx, nx=dataX.shape;
			else:
				raise ValueError, "data for x dimension > 2!";
			# data: (my, 1)
			if (dataY.ndim==1):
				my=dataY.size;
				ny=1
			# data: (my, ny)
			elif (dataY.ndim==2):
				my, ny=dataY.shape;
			else:
				raise ValueError, "data for y dimension > 2!";

			# data: (mx, nx) * (my, 1)
			wts=None
			if (mx==my and ny==1):
				wts=np.zeros((nx,1))
				wts=(np.linalg.inv(dataX.T.dot(dataX)+rlamda*np.eye(nx)).dot(dataX.T.dot(dataY)))
				#wts=np.zeros((1,nx))
				#wts[:, 0]=(np.linalg.inv(dataX.T.dot(dataX)+rlamda*np.eye(nx)).dot(dataX.T.dot(dataY))).T
			# data: (mx, nx) * (1, my)
			elif (mx==ny and my==1):
				#wts=np.zeros((1,nx))
				#wts[:,0]=(np.linalg.inv(dataX.T.dot(dataX)+rlamda*np.eye(nx)).dot(dataX.T.dot(dataY.T))).T
				wts=np.zeros((nx,1))
				wts=(np.linalg.inv(dataX.T.dot(dataX)+rlamda*np.eye(nx)).dot(dataX.T.dot(dataY.T)))
			### Should never happen if give correct X and Y
			# data: (my, ny)*(mx, 1)
			elif ( mx==ny and nx==1):
				#wts=np.zeros((1,ny))
				#wts[:,0]=(np.linalg.inv(dataY.T.dot(dataY)+rlamda*np.eye(ny)).dot(y.T.dot(dataX))).T
				wts=np.zeros((ny,1))
				wts=(np.linalg.inv(dataY.T.dot(dataY)+rlamda*np.eye(ny)).dot(y.T.dot(dataX)))
			# data: (mx, nx) * (1, my)
			elif (mx==my and mx==1):
				#wts=np.zeros((1,ny))
				#wts[:,0]=(np.linalg.inv(dataY.T.dot(dataY)+rlamda*np.eye(ny)).dot(y.T.dot(dataX.T))).T
				wts=np.zeros((ny,1))
				wts=(np.linalg.inv(dataY.T.dot(dataY)+rlamda*np.eye(ny)).dot(y.T.dot(dataX.T)))
			else:
				raise ValueError, "Not match array size for x and y!";
			return Data2D(data=wts,name="weight");
		#except:
		#	traceback.print_exc()
		#	exit(1)

	def calcY(self,weight):
		'''Calculate the Y' for given Weight based on current X'''
		if (isinstance(weight,list) or isinstance(weight,tuple)):
			weight=Data2D(weight,name="weight")
		if (weight.shape[0] is 1): 
			weight.resize((weight.size,1));
		Y=self.X.dot(weight)
		return Y

	def calcWeight(self,Y, rlamda=None):
		'''Calculate the Weight for given Y based on current X'''
		if (isinstance(Y,list) or isinstance(Y,tuple)):
			Y=Data2D(Y,name="Y")
		if (rlamda is None):
			rlamda=self.rlamda;
		if (Y.shape[0] is 1): 
			Y.resize((Y.size,1));
		return self.LeastSquares(self.X, Y, rlamda=rlamda);

	def calcError(self, Y=None, Abs=False):
		'''Calculate the Error array for given Y and current Y
		- Abs: Calculate Absolute Error array.
		- if Y is not given, calculate current delta Y'''
		if (isinstance(Y,list) or isinstance(Y,tuple)):
			Y=Data2D(Y,name="Y")
		if (Y is None): 
			Y=self.calcY(self.W);
		if (Y.shape[0] is 1): 
			Y.resize((Y.size,1));
		if (Abs):
			return np.abs(Y-self.Y)
		else:		
			return Y-self.Y

	def calcRMS(self,Err=None, Y=None, weight=None):
		'''Calculate the RMS based on given data
		- If Err is given, calculate RMS directly.
		- If no Err and Y given, calculate Error based on given Y firstly
		- If no Err & Y, but weight given, calculate Y based on given weight firstly
		- If no Err & Y & weight, calculate RMS based on current Y'=XW '''

		if (isinstance(Err,list) or isinstance(Err,tuple)):
			Err=Data2D(Err,name="Error")
		if (isinstance(Y,list) or isinstance(Y,tuple)):
			Y=Data2D(Y,name="Y")
		if (isinstance(weight,list) or isinstance(weight,tuple)):
			weight=Data2D(weight,name="weight")

		if (Err is None and Y is None and weight is None):
			weight=self.W
			if (weight.shape[0] is 1): weight.resize((weight.size,1));			
			Err=self.X.dot(weight)-self.Y;
		if (Err is None and Y is None and weight is not None):
			if (weight.shape[0] is 1): weight.resize((weight.size,1));
			Y=self.calcY(weight)
		if (Err is None and Y is not None):
			if (Y.shape[0] is 1): Y.resize((Y.size,1));
			Err=self.calcError(Y)
		# Calculate the RMS for deltaY Error Array
		return np.sqrt((Err.T.dot(Err)).ravel()[0]/Err.shape[0])

	def getErrorRMS(self,weight=None,rlamda=None):
		'''Calculate Error and RMS based on current X and Y.
		- When weight is given, using the given weight;
		- When weight is not given and rlamda is given, recalculate the weight
		- When weight and rlamda is not given, using the current W''' 
		if ( weight is None and rlamda is None): 
			weight=self.W;
		elif ( weight is None and rlamda is not None):
			weight=self.calcWeight(self.Y,rlamda);
		return self.calcError(self.calcY(weight)),self.calcRMS(weight=weight)

	##### Methods for train-test set
	def testY(self, train=None):
		'''Calculate the Y' based on training set W'''
		if (train is None):
			WT=self.WT
		else:
			WT=train.W
		return self.calcY(WT)

	def testError(self, train=None):
		'''Calculate the Error array based on training set W'''
		return self.calcError(self.testY(train));

	def testRMS(self, train=None):
		'''Calculate the RMS based on training set W'''
		Err=self.testError(train);
		return np.sqrt((Err.T.dot(Err)).ravel()[0]/Err.shape[0])

