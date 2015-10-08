import __init__

from HHChain import * 
from HHResidue import *
from HHAtom import *

class Molecule:

	name=""
	resname=""
	atoms=[]
	reses=[]
	bonds=[]
	rings=[]
	nfrag=1
	mtype=""

	def __init__(self):
		pass

	def GetNumAtom(self):
		return len(self.atoms);

	def GetNumRes(self):
		return len(self.reses);

	def GetNumBond(self):
		return len(self.bonds);

	def GetNumRing(self):
		return len(self.rings);

	def GetNumFrag(self):
		return self.nfrag;

	def __str__(self):
		out=self.name+": "+str(len(self.atoms))+" atoms;"+str(len(self.reses))+" residues;"
		return out

	attrflags={
		"Perceive_Connect":False,
		"Perceive_Aromatic":False,
		"Perceive_Ring":False,
		"Perceive_AtomType":False,
		"Perceive_BondType":False,
		"Perceive_FuncGroup":False
	}