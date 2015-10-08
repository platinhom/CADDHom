import __init__
import os,sys,math
from copy import deepcopy

class FileFormator:
	filename=""
	mode='r'
	extension=[]
	handle=None
	
	def open(self,filename,mode='r'):
		self.filename=filename
		self.mode=mode
		self.handle=open(filename,mode)

	def close(self):
		self.handle.close();

	def write(self,line):
		self.handle.write(line);

	def readline(self):
		return self.handle.readline();

	def readlines(self):
		return self.handle.readlines();
