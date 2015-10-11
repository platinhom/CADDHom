# -*- coding: utf8 -*-
''' Ligand processor:
Ligand
Conformation
Force field
2D Feature
3D Feature
Fingerprint
SMILE
SMART
PH4 Feature, Pharmacophore.
Similarity
LigCluster
QSAR
Optimize
'''
__author__="Zhixiong Zhao"
__version__="0.1"
__date__="2015.10.10"

import __init__
from HHMolecule import *

class Ligand:        
    def __init__(self):
        self.Allatoms={}
        self.resname=''

class PH4Feature:
	'''Pharmacophore feature'''
	pass

class Pharmacophore:
	pass