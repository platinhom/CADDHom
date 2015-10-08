import __init__
import copy

class HHBase():
	attrflags={}
	def has(self,attr):
		return self.__dict__.has_key(attr)
	def getflag(self,flag):
		return attrflags[flag];
	def __getattr__(self,attr):
		if self.defaults.has_key(attr):
			return copy.deepcopy(self.defaults[attr])
		else:
			raise AttributeError(attr)