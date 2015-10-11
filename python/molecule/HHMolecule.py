'''
Molecule and its components:
Element -> Atom
Bond, Ring
Kekulize, Aromatic
Chiral
FuncGroup
Residue, Chain
Molecule
'''

__author__="Zhixiong Zhao"
__version__="0.1"
__date__="2015.10.10"
import __init__

from basic.geometry import Point,Vector


# -------------------------- Element ---------------------------------

def GetMass(element_name):
	return atomic_mass[element_name.upper()];

def GetElementNum(element_name):
	return atomic_number[element_name.upper()];

class Element:
	def __init__(self,element_name=""):
		assert isinstance(element_name,str);
		self.element_name=element_name;
		self.element_valid=self.ValidateElement();

	def ValidateElement(self,name=None):
		if (name):
			return (name in atomic_mass)
		else:
			return (self.element_name.upper() in atomic_mass )

	def GetMass(self):
		return atomic_mass[self.element_name.upper()]

	def GetElementNumber(self):
		return atomic_number[self.element_name.upper()]

    # Deduce Element from name
    def DeduceElementFromName(self,name=""):
        name=name.strip();
        uname=name.upper();
        if (len(name)==1 and uname in atomic_name ):
            return uname;
        elif (len(name)>=2):
        	uname2=uname[:2]
        	if (uname2=="CA" and name[:2]=="CA" || 
        		uname2=="CO" and name[:2]=="CO"):
        		return "C"
        	if (uname2 in atomic_name):
        		return (atomic_name[uname2])
        	elif (uname2[0] in atomic_name):
        		return uname2[0];
        raise ValueError("Can't deduce element from name! : "+name);

#-------------------- Atom -----------------------
class Atom(Element):

    name=''
    atype=''
    index=0
    atomid=0
    resname=''
    resid=0
    res=None
    mol=None
    chain=None
    isring=False
    isaromatic=False
    coordinates = Point(99999, 99999, 99999)
    x=coordinates.coors()[0]
    y=coordinates.coors()[1]
    z=coordinates.coors()[2]
    fcharge=0.0
    pcharge=0.0

    #element_name
    #element_valid

    def __init__(self):
        Element.__init__(self,"")
        self.undo_coordinates = Point(99999, 99999, 99999)
        self.line=""
        self.PDBIndex = ""
        self.IndeciesOfAtomsConnecting=[]
    
    def __str__(self):
        return self.CreatePDBLine()
		
    def atom_radius(self):
        element=self.element.upper()
        if element=="H": return 1.2
        if element=="C": return 1.7
        if element=="N": return 1.55
        if element=="O": return 1.52
        if element=="F": return 1.47
        if element=="P": return 1.8
        if element=="S": return 1.8
        return 1.0 # the default
    
    # function to determine if the atom belongs to a protein
    # returns true or false
    def belongs_to_protein(self):
        protein_residues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        if self.resname.strip() in protein_residues: return True
        else: return False

    def CopyOf(self):
        newatom = atom()
        newatom.chain = self.chain
        newatom.resid = self.resid
        newatom.name = self.name
        newatom.resname = self.resname
        newatom.coordinates = self.coordinates.CopyOf()
        newatom.undo_coordinates = self.undo_coordinates.CopyOf()
        newatom.element = self.element
        newatom.PDBIndex = self.PDBIndex
        newatom.line = self.line
        for index in self.IndeciesOfAtomsConnecting:
            newatom.IndeciesOfAtomsConnecting.append(index)
        return newatom

    # Reads text (PDB format) into an atom object
    # Requires: A string containing the PDB line
    def ReadPDBLine(self, Line):
        self.line = Line
        self.name = Line[11:16].strip()
        self.chain = Line[21:22]
        if Line[22:26].strip() != "":
            self.resid = int(Line[22:26])
        else:
            self.resid = 0
        
        if len(self.name)==1: # redo using rjust
            self.name = self.name + "  "
        elif len(self.name)==2:
            self.name = self.name + " "
        elif len(self.name)==3:
            self.name = self.name + " " # This line is necessary for babel to work, though many PDBs in the PDB would have this line commented out
        
        self.coordinates = point(float(Line[30:38]), float(Line[38:46]), float(Line[46:54]))
        
        if len(Line) >= 79:
            self.element = Line[76:79].strip().upper() # element specified explicitly at end of life
        elif self.element == "": # try to guess at element from name
            two_letters = self.name[0:2].strip().upper()
            if two_letters=='BR':
                self.element='BR'
            elif two_letters=='CL':
                self.element='CL'
            elif two_letters=='BI':
                self.element='BI'
            elif two_letters=='AS':
                self.element='AS'
            elif two_letters=='AG':
                self.element='AG'
            elif two_letters=='LI':
                self.element='LI'
            elif two_letters=='HG':
                self.element='HG'
            elif two_letters=='MG':
                self.element='MG'
            elif two_letters=='RH':
                self.element='RH'
            elif two_letters=='ZN':
                self.element='ZN'
            else: #So, just assume it's the first letter.
                self.element = self.name[0:1].strip().upper()
                
        # Any number needs to be removed from the element name
        self.element = self.element.replace('0','')
        self.element = self.element.replace('1','')
        self.element = self.element.replace('2','')
        self.element = self.element.replace('3','')
        self.element = self.element.replace('4','')
        self.element = self.element.replace('5','')
        self.element = self.element.replace('6','')
        self.element = self.element.replace('7','')
        self.element = self.element.replace('8','')
        self.element = self.element.replace('9','')

        self.PDBIndex = Line[6:12].strip()
        self.resname = Line[16:20]
        if self.resname.strip() == "": self.resname = " MOL"
    # Creates a PDB line from the atom object
    # Returns: PDB String
    def CreatePDBLine(self):

        #if len(self.name) > 1: self.name = self.name[:1].upper() + self.name[1:].lower()

        output = "ATOM "
        #output = output + str(index).rjust(6) + self.name.rjust(5) + self.residue.rjust(4)
        output = output + str(self.index).rjust(6) + self.name.rjust(5) + self.resname.rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(18)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + self.element_name.rjust(24) # + "   " + str(uniqueID) #This last part must be removed
        return output

    # Sets the undo point for later undoing
    def SetUndoPoint(self):
        self.undo_coordinates = self.coordinates.CopyOf()
        
    # Resets coordinate values after translations or rotations ("Undo")
    def Undo(self):
        self.coordinates = self.undo_coordinates.CopyOf()

#---------------------------  Bond --------------------------------------
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

#---------------------------  Ring --------------------------------------
class Ring:
    def __init(self,atomindex_tuple):
        self.atomindex_tuple

# -------------------------- Aromatic -----------------------------------
class Aromatic:
	pass
# -------------------------- Kekulize -----------------------------------
class Kekulize:
	pass
# -------------------------- Chiral -------------------------------------
class Chiral:
	pass

#-------------------------- FuncGroup -----------------------------------
class FuncGroup:
	def __init__(self):
		pass

#--------------------------- Residue  -----------------------------------
class Residue:
    def __init__(self):
        self.resname=''
        self.resid=0
        self.chain=''
        self.CA=Atom()
        self.attr=''
        self.Allatoms={}


#---------------------------- Chain -------------------------------------
class Chain:
	def __init__(self):
		pass


#--------------------------- Molecule -----------------------------------
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


######################## DATA PART ################################

# hvd values obtained from http://www.webelements.com/ and recorded to their known accuracy.
# Use upper case to validate the element name
# Only normal elements Include X, Du
# No R and so on.

atomic_name = {
	'H':'H','C':'C','N':'N','O':'O','F':'F','P':'P','S':'S','B':'B','I':'I','K':'K','X':'X','DU':'Du','CL':'Cl','BR':'Br',
	'NA':'Na','SI':'Si','CA':'Ca','MG':'Mg','AL':'Al','TI':'Ti','FE':'Fe','CU':'Cu','NI':'Ni'
}

atomic_mass = {
	'H'	 :	 1.00794,
	'He' :	 4.002602,
	'HE' :	 4.002602,
	'Li' :	 6.941,
	'LI' :	 6.941,
	'Be' :	 9.012182,
	'BE' :	 9.012182,
	'B'	 :	10.811,
	'C'	 :	12.0107,
	'N'	 :	14.0067,
	'O'	 :	15.9994,
	'F'	 :	18.9984032,
	'Ne' :	20.1797,
	'NE' :	20.1797,
	'Na' :	22.989770,
	'NA' :	22.989770,
	'Mg' :	24.3050,
	'MG' :	24.3050,
	'Al' :	26.981538,
	'AL' :	26.981538,
	'Si' :	28.0855,
	'SI' :	28.0855,
	'P'	 :	30.973761,
	'S'	 :	32.065,
	'Cl' :	35.453,
	'CL' :	35.453,
	'Ar' :	39.948,
	'AR' :	39.948,
	'K'	 :	39.0983,
	'Ca' :	40.078,
	'CA' :	40.078,
	'Sc' :	44.955910,
	'SC' :	44.955910,
	'Ti' :	47.867,
	'TI' :	47.867,
	'V'	 :	50.9415,
	'Cr' :	51.9961,
	'CR' :	51.9961,
	'Mn' :	54.938049,
	'MN' :	54.938049,
	'Fe' :	55.845,
	'FE' :	55.845,
	'Co' :	58.933200,
	'CO' :	58.933200,
	'Ni' :	58.6934,
	'NI' :	58.6934,
	'Cu' :	63.546,
	'CU' :	63.546,
	'Zn' :	65.39,
	'ZN' :	65.39,
	'Ga' :	69.723,
	'GA' :	69.723,
	'Ge' :	72.64,
	'GE' :	72.64,
	'As' :	74.92160,
	'AS' :	74.92160,
	'Se' :	78.96,
	'SE' :	78.96,
	'Br' :	79.904,
	'BR' :	79.904,	  
	'Kr' :	83.80,
	'KR' :	83.80,
	'Rb' :	85.4678,
	'RB' :	85.4678,
	'Sr' :	87.62,
	'SR' :	87.62,
	'Y'	 :	88.90585,
	'Zr' :	91.224,
	'ZR' :	91.224,
	'Nb' :	92.90638,
	'NB' :	92.90638,
	'Mo' :	95.94,
	'MO' :	95.94,
	'Tc' :	98,
	'TC' :	98,
	'Ru' : 101.07,
	'RU' : 101.07,
	'Rh' : 102.90550,
	'RH' : 102.90550,
	'Pd' : 106.42,
	'PD' : 106.42,
	'Ag' : 107.8682,
	'AG' : 107.8682,
	'Cd' : 112.411,
	'CD' : 112.411,
	'In' : 114.818,
	'IN' : 114.818,
	'Sn' : 118.710,
	'SN' : 118.710,
	'Sb' : 121.760,
	'SB' : 121.760,
	'Te' : 127.60,
	'TE' : 127.60,
	'I'	 : 126.90447,
	'Xe' : 131.293,
	'XE' : 131.293,
	'Cs' : 132.90545,
	'CS' : 132.90545,
	'Ba' : 137.327,
	'BA' : 137.327,
	'La' : 138.9055,
	'LA' : 138.9055,
	'Ce' : 140.116,
	'CE' : 140.116,
	'Pr' : 140.90765,
	'PR' : 140.90765,
	'Nd' : 144.24,
	'ND' : 144.24,
	'Pm' : 145,
	'PM' : 145,
	'Sm' : 150.36,
	'SM' : 150.36,
	'Eu' : 151.964,
	'EU' : 151.964,
	'Gd' : 157.25,
	'GD' : 157.25,
	'Tb' : 158.92534,
	'TB' : 158.92534,
	'Dy' : 162.50,
	'DY' : 162.50,
	'Ho' : 164.93032,
	'HO' : 164.93032,
	'Er' : 167.259,
	'ER' : 167.259,
	'Tm' : 168.93421,
	'TM' : 168.93421,
	'Yb' : 173.04,
	'YB' : 173.04,
	'Lu' : 174.967,
	'LU' : 174.967,
	'Hf' : 178.49,
	'HF' : 178.49,
	'Ta' : 180.9479,
	'TA' : 180.9479,
	'W'	 : 183.84,
	'Re' : 186.207,
	'RE' : 186.207,
	'Os' : 190.23,
	'OS' : 190.23,
	'Ir' : 192.217,
	'IR' : 192.217,
	'Pt' : 195.078,
	'PT' : 195.078,
	'Au' : 196.96655,
	'AU' : 196.96655,
	'Hg' : 200.59,
	'HG' : 200.59,
	'Tl' : 204.3833,
	'TL' : 204.3833,
	'Pb' : 207.2,
	'PB' : 207.2,
	'Bi' : 208.98038,
	'BI' : 208.98038,
	'Po' : 208.98,
	'PO' : 208.98,
	'At' : 209.99,
	'AT' : 209.99,
	'Rn' : 222.02,
	'RN' : 222.02,
	'Fr' : 223.02,
	'FR' : 223.02,
	'Ra' : 226.03,
	'RA' : 226.03,
	'Ac' : 227.03,
	'AC' : 227.03,
	'Th' : 232.0381,
	'TH' : 232.0381,
	'Pa' : 231.03588,
	'PA' : 231.03588,
	'U'	 : 238.02891,
	'Np' : 237.05,
	'NP' : 237.05,
	'Pu' : 244.06,
	'PU' : 244.06,
	'Am' : 243.06,
	'AM' : 243.06,
	'Cm' : 247.07,
	'CM' : 247.07,
	'Bk' : 247.07,
	'BK' : 247.07,
	'Cf' : 251.08,
	'CF' : 251.08,
	'Es' : 252.08,
	'ES' : 252.08,
	'Fm' : 257.10,
	'FM' : 257.10,
	'Md' : 258.10,
	'MD' : 258.10,
	'No' : 259.10,
	'NO' : 259.10,
	'Lr' : 262.11,
	'LR' : 262.11,
	'Rf' : 261.11,
	'RF' : 261.11,
	'Db' : 262.11,
	'DB' : 262.11,
	'Sg' : 266.12,
	'SG' : 266.12,
	'Bh' : 264.12,
	'BH' : 264.12,
	'Hs' : 269.13,
	'HS' : 269.13,
	'Mt' : 268.14,
	'MT' : 268.14,
	}

atomic_number = {
	'H'	 :	 1,
	'He' :	 2,
	'HE' :	 2,
	'Li' :	 3,
	'LI' :	 3,
	'Be' :	 4,
	'BE' :	 4,
	'B'	 :	 5,
	'C'	 :	 6,
	'N'	 :	 7,
	'O'	 :	 8,
	'F'	 :	 9,
	'Ne' :	10,
	'NE' :	10,
	'Na' :	11,
	'NA' :	11,
	'Mg' :	12,
	'MG' :	12,
	'Al' :	13,
	'AL' :	13,
	'Si' :	14,
	'SI' :	14,
	'P'	 :	15,
	'S'	 :	16,
	'Cl' :	17,
	'CL' :	17,
	'Ar' :	18,
	'AR' :	18,
	'K'	 :	19,
	'Ca' :	20,
	'CA' :	20,
	'Sc' :	21,
	'SC' :	21,
	'Ti' :	22,
	'TI' :	22,
	'V'	 :	23,
	'Cr' :	24,
	'CR' :	24,
	'Mn' :	25,
	'MN' :	25,
	'Fe' :	26,
	'FE' :	26,
	'Co' :	27,
	'CO' :	27,
	'Ni' :	28,
	'NI' :	28,
	'Cu' :	29,
	'CU' :	29,
	'Zn' :	30,
	'ZN' :	30,
	'Ga' :	31,
	'GA' :	31,
	'Ge' :	32,
	'GE' :	32,
	'As' :	33,
	'AS' :	33,
	'Se' :	34,
	'SE' :	34,
	'Br' :	35,
	'BR' :	35,
	'Kr' :	36,
	'KR' :	36,
	'Rb' :	37,
	'RB' :	37,
	'Sr' :	38,
	'SR' :	38,
	'Y'	 :	39,
	'Zr' :	40,
	'ZR' :	40,
	'Nb' :	41,
	'NB' :	41,
	'Mo' :	42,
	'MO' :	42,
	'Tc' :	43,
	'TC' :	43,
	'Ru' :	44,
	'RU' :	44,
	'Rh' :	45,
	'RH' :	45,
	'Pd' :	46,
	'PD' :	46,
	'Ag' :	47,
	'AG' :	47,
	'Cd' :	48,
	'CD' :	48,
	'In' :	49,
	'IN' :	49,
	'Sn' :	50,
	'SN' :	50,
	'Sb' :	51,
	'SB' :	51,
	'Te' :	52,
	'TE' :	52,
	'I'	 :	53,
	'Xe' :	54,
	'XE' :	54,
	'Cs' :	55,
	'CS' :	55,
	'Ba' :	56,
	'BA' :	56,
	'La' :	57,
	'LA' :	57,
	'Ce' :	58,
	'CE' :	58,
	'Pr' :	59,
	'PR' :	59,
	'Nd' :	60,
	'ND' :	60,
	'Pm' :	61,
	'PM' :	61,
	'Sm' :	62,
	'SM' :	62,
	'Eu' :	63,
	'EU' :	63,
	'Gd' :	64,
	'GD' :	64,
	'Tb' :	65,
	'TB' :	65,
	'Dy' :	66,
	'DY' :	66,
	'Ho' :	67,
	'HO' :	67,
	'Er' :	68,
	'ER' :	68,
	'Tm' :	69,
	'TM' :	69,
	'Yb' :	70,
	'YB' :	70,
	'Lu' :	71,
	'LU' :	71,
	'Hf' :	72,
	'HF' :	72,
	'Ta' :	73,
	'TA' :	73,
	'W'	 :	74,
	'Re' :	75,
	'RE' :	75,
	'Os' :	76,
	'OS' :	76,
	'Ir' :	77,
	'IR' :	77,
	'Pt' :	78,
	'PT' :	78,
	'Au' :	79,
	'AU' :	79,
	'Hg' :	80,
	'HG' :	80,
	'Tl' :	81,
	'TL' :	81,
	'Pb' :	82,
	'PB' :	82,
	'Bi' :	83,
	'BI' :	83,
	'Po' :	84,
	'PO' :	84,
	'At' :	85,
	'AT' :	85,
	'Rn' :	86,
	'RN' :	86,
	'Fr' :	87,
	'FR' :	87,
	'Ra' :	88,
	'RA' :	88,
	'Ac' :	89,
	'AC' :	89,
	'Th' :	90,
	'TH' :	90,
	'Pa' :	91,
	'PA' :	91,
	'U'	 :	92,
	'Np' :	93,
	'NP' :	93,
	'Pu' :	94,
	'PU' :	94,
	'Am' :	95,
	'AM' :	95,
	'Cm' :	96,
	'CM' :	96,
	'Bk' :	97,
	'BK' :	97,
	'Cf' :	98,
	'CF' :	98,
	'Es' :	99,
	'ES' :	99,
	'Fm' : 100,
	'FM' : 100,
	'Md' : 101,
	'MD' : 101,
	'No' : 102,
	'NO' : 102,
	'Lr' : 103,
	'LR' : 103,
	'Rf' : 104,
	'RF' : 104,
	'Db' : 105,
	'DB' : 105,
	'Sg' : 106,
	'SG' : 106,
	'Bh' : 107,
	'BH' : 107,
	'Hs' : 108,
	'HS' : 108,
	'Mt' : 109,
	'MT' : 109
	}

implicit_valence = {
	'H'	 :	{0:1,1:0},
	'C'	 :	{0:4,1:3,2:2,3:1,4:0},
	'N'	 :	{0:3,1:2,2:1,3:0,4:0},
	'O'	 :	{0:2,1:1,2:0},
	'F'	 :	{0:1,1:0},
	'Cl' :	{0:1,1:0},
	'CL' :	{0:1,1:0},
	'Br' :	{0:1,1:0},
	'BR' :	{0:1,1:0},
	'I'	 :	{0:1,1:0},
	'S'	 :	{0:2,1:2,2:0,3:1,4:0,5:1,6:0}, # ambiguity?
	'K'	 :	{0:0,1:0},	   # as drawn
	'Cu' :	{0:0,1:0,2:0}, # as drawn
	'CU' :	{0:0,1:0,2:0}, # as drawn
	'Zn' :	{0:0,1:0,2:0,4:0}, # as drawn
	'ZN' :	{0:0,1:0,2:0,4:0}, # as drawn
	'Mg' :	{0:1,1:0},
	'MG' :	{0:1,1:0},
	'Ca' :	{0:1,1:0},
	'CA' :	{0:1,1:0},
	'P'	 :	{0:3,1:2,2:1,3:0,4:0,5:0,6:0},
	}