
class HHBase():
	def has(self,attr):
		return self.__dict__.has_key(attr)
	def __getattr__(self,attr):
        if self.defaults.has_key(attr):
            return copy.deepcopy(self.defaults[attr])
        else:
            raise AttributeError(attr)