import __init__
from HHAtom import *

class Bond:
	bgn=None
	end=None
	idx_bgn=0
	idx_end=0
	order=0 #1,2,3,4(breaking),5(aromatic), 6(resonance,1.5)
	btype="" #1,2,3 AR,AM,DU
	index=0
	mol=None

	def __init__(self):
		pass

	def __str__(self):
		return "Bond "+str(self.index)+": Atoms "+str(self.idx_bgn)+" - "+str(self.idx_end)
				+"; Order: "+str(self.order)+"; BondType: "+self.btype;

	def CreateBond_Order(self,bgn,end,order=1,btype="1"):
		bond=Bond();
		bond.bgn=bgn
		bond.end=end
		bond.order=order
		bond.idx_bgn=bgn.index
		bond.idx_end=end.index

	def SetBgnFromIdx(self):
		if (mol):
			atom=self.mol.atoms[self.idx_bgn-1];
			if (atom.index==self.idx_bgn):
				self.bgn=atom;
			else:
				raise ValueError("Can't find atom at index: "+self.idx_bgn)
		else: raise ValueError("No molecule for atom at index: "+self.idx_bgn)

	def SetEndFromIdx(self):
		if (mol):
			atom=self.mol.atoms[self.idx_end-1];
			if (atom.index==self.idx_end):
				self.end=atom;
			else:
				raise ValueError("Can't find atom at index: "+self.idx_end)
		else: raise ValueError("No molecule for atom at index: "+self.idx_end)
	
	def SetAtomsFromIdx(self):
		self.SetBgnFromIdx();
		self.SetEndFromIdx();

