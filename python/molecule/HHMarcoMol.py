# -*- coding: utf8 -*-
'''
For macromolecule, protein, nucleac acid, polymer
Sequence
Protein
NucAcid
Complex
MacroStruture
Polymer
'''

__author__="Zhixiong Zhao"
__version__="0.1"
__date__="2015.10.10"

import __init__
from HHMolecule import *

class Sequence:
	pass;

class Protein(HHMolecule.Molecule):
    def __init__(self):
        self.Allatoms={}
        self.Allres={}

class NucAcid:
	def __init__(self):
		pass

class Complex:
	def __init__(self):
		pass

class MacroStruture:
	def __init__(self):
		pass