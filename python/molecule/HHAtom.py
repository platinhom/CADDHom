import __init__
from geometry.HHPoint import *
from HHElement import *

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


