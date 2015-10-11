# -*- coding: utf-8 -*-
'''
Molecular File Format Process
Support: mol2/pdb/sdf/pqr/gjf/xyz/xyzr
'''

__author__="Zhixiong Zhao"
__version__="0.1"
__date__="2015.10.10"


import __init__
import os,sys,math
from copy import deepcopy
import molecule.HHMolecule 
import molecule.HHAtom
import molecule.HHResidue
import molecule.HHBond
import geometry.HHPoint

Mol=molecule.HHMolecule.Molecule
Atom=molecule.HHAtom.Atom
Res=molecule.HHResidue.Residue
Bond=molecule.HHBond.Bond
Point=geometry.HHPoint.Point

class FileFormator(object):
	'''Base class for all format'''

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

#----------------------- mol2 ---------------------------

class MOL2(FileFormator):
    extension=['mol2'];

    def CreateAtomLine(self, atom, lenatom=4, lenres=3):
        output=atom.index.rjust(lenatom)+" "+atom.name.ljust(5)
        output+=("%.4f" % atom.coordinates.x).rjust(11) + ("%.4f" % atom.coordinates.y).rjust(11)+ ("%.4f" % atom.coordinates.z).rjust(11)+ ' '
        output+=atom.atype.ljust(6)+str(atom.resid).rjust(lenres)+ ' ' + atom.resname.ljust(6)+ atom.pcharge.rjust(9)+ os.linesep
        return output

    def CreateBondline(bond,lenbond=4):
        output=bond.index.rjust(lenbond)+" "+bond.idx_bgn.rjust(lenbond)+" "+\
                bond.idx_end.rjust(lenbond)+"  "+bond.btype.lower().ljust(lenbond)+ os.linesep
        return output

    def WriteObj(self,obj):
        if (isinstance(obj,Atom)):
            self.write(CreateAtomLine(obj))
        elif(isinstance(obj,Res) or isinstance(obj,Mol)):
            for atom in obj.atoms:
               self.write(CreateAtomLine(atom))
        elif(isinstance(obj,Bond)):
            self.write(CreateBondline(obj));
        else:
            self.write(str(obj));

    def ReadAtomLine(self, Line):
        items=Line.split()
        atom=Atom()
        atom.index = int(items[0])
        atom.atomid = int(items[0])
        atom.name = items[1]
        atom.coordinates = Point(float(items[2]), float(items[3]), float(items[4]))
        atom.atype=items[5]
        #sybyl type
        #atom.element_name=atom.atype[0:2].strip('.').strip()
        atom.element_name=atom.DeduceElementFromName(atom.name);
        if len(items)==9:
            atom.resid = int(items[6])
            atom.resname = items[7]
            atom.charge = items[8]
        return atom;

    def ReadBondLine(self, Line):
        items=Line.split()
        bond=Bond()
        bond.index = int(items[0])
        bond.idx_bgn = int(items[1])
        bond.idx_bgn = int(items[2])
        bond.btype = items[3]
        return bond;

    def WriteMolFile(self,mol,filename):
        self.open(filename,'w');
        self.write("@<TRIPOS>MOLECULE\n")
        self.write(mol.name+'\n')
        self.write("%5d %5d %5d %5d %5d \n", mol.GetNumAtom(), mol.GetNumBond(), mol.GetNumFrag(), 0, 0);

        self.write("@<TRIPOS>ATOM\n");
        self.WriteObj(mol);
        self.write("@<TRIPOS>BOND\n");

    def ReadMolFile(self, filename):
        self.open(filename,'r');
        findmol=False;
        findatom=False;
        findbond=False;
        nextmol=False;
        mols=[]
        mol=None
        for line in self.handle:
            if (line[:17] == "@<TRIPOS>MOLECULE"):
                findmol=True;
                findatom=False;
                findbond=False;
                if (nextmol):
                    mols.append(mol)
                    nextmol=False;
                mol=Mol()
                continue;
            if (line[:13] == "@<TRIPOS>ATOM"):
                findatom=True;
                findmol=False;
                nextmol=True;
                continue;
            if (line[:13] == "@<TRIPOS>BOND"):
                findatom=False;
                findbond=True;
                continue;
            if (findbond and line[:9]=="@<TRIPOS>"):
                findbond=False;
                continue;
            if (findatom):
                atom=self.ReadAtomLine(line);
                atom.mol=mol;
                mol.atoms.append();
            if (findbond):
                bond=self.ReadBondLine(line);
                bond.mol=mol;
                bond.SetAtomsFromIdx()
                mol.bonds.append(bond);
        mols.append(mol);
        self.close();
        if (len(mols)==1):return mols[0];
        elif (len(mols)>1):return mols;
        elif (len(mols)==0):return None;

#----------------------- pdb ---------------------------	
class PDB(FileFormator):
	# PDB class
    def __init__ (self):
        self.AllAtoms={}

    # Loads a PDB from a file
    # Requires: FileName, a string containing the filename
    def LoadPDB(self, FileName):

        autoindex = 1

        self.__init__()

        # Now load the file into a list
        file = open(FileName,"r")
        lines = file.readlines()
        file.close()

        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    TempAtom = atom()
                    TempAtom.ReadPDBLine(line)

                    self.AllAtoms[autoindex] = TempAtom # because points files have no indecies
                    autoindex = autoindex + 1

    # Saves the PDB object to a PDB file
    # Requires: filename to be saved
    def SavePDB(self, filename):

        if len(self.AllAtoms) > 0: # so the pdb is not empty (if it is empty, don't save)

            file = open(filename,"w")

            # write coordinates
            for atomindex in self.AllAtoms:
                file.write(self.AllAtoms[atomindex].CreatePDBLine() + "\n")

            file.close()

    # identifies the greatest distance between any two atoms of the PDB
    # Used for quickly comparing two PDBs so RMSD alignment not necessary if they're different
    # Returns float = distance
    def farthest_away_atoms_dist(self):
        dist_max = 0.0
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                if dist > dist_max: dist_max = dist
        self.max_inter_atom_distance = dist_max
        return dist_max

    # identifies the smallest distance between any two atoms of the PDB
    # Used for quickly comparing two PDBs so RMSD alignment not necessary if they're different
    # Returns float = distance
    def closest_atoms_dist(self):
        dist_min = 999999999.99
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                if dist < dist_min: dist_min = dist
        self.min_inter_atom_distance = dist_min
        return dist_min

    # Creates a "distance fingerprint" so PDB's can be quickly compared w/o RMSD alignment
    # Returns: a list containing sorted distances
    def distance_fingerprint(self):
        fingerprint = []
        keys = self.AllAtoms.keys()
        for key_index1 in range(0,len(keys)-1):
            key1 = keys[key_index1]
            atom1 = self.AllAtoms[key1]
            for key_index2 in range(key_index1+1,len(keys)):
                key2 = keys[key_index2]
                atom2 = self.AllAtoms[key2]
                dist = atom1.coordinates.dist_to(atom2.coordinates)
                fingerprint.append(dist)
        fingerprint.sort()
        self.distance_fingerprint = fingerprint
        return fingerprint

    # To detect if two distance fingerprints are sufficiently similar
    # Requires: Two fingerprints (lists) and a tolerance (float)
    # Returns: True or False
    def fingerprint_same_as(self, fingerprint, tol):
        if len(self.distance_fingerprint) <> len(fingerprint): return False

        for index in range(len(self.distance_fingerprint)):
                item1 = self.distance_fingerprint[index]
                item2 = fingerprint[index]
                if math.fabs(item1-item2) > tol: return False
        return True

    # Print out info about the PDB
    def print_out_info(self):
        for index in self.AllAtoms:
            print self.AllAtoms[index].CreatePDBLine()

    # Translate the molecule
    def TranslateMolecule(self, x, y, z):
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].coordinates.x = self.AllAtoms[atomindex].coordinates.x + x
            self.AllAtoms[atomindex].coordinates.y = self.AllAtoms[atomindex].coordinates.y + y
            self.AllAtoms[atomindex].coordinates.z = self.AllAtoms[atomindex].coordinates.z + z

    # Rotate the PDB around it's center.
    # Requires: degrees to rotate, as float
    def RotatePDB(self, thetax, thetay, thetaz):

        # first, identify the geometric center
        x = 0.0
        y = 0.0
        z = 0.0
        count = 0
        for index in self.AllAtoms:
                if self.AllAtoms[index].element != "H":
                        count = count + 1
                        x = x + self.AllAtoms[index].coordinates.x
                        y = y + self.AllAtoms[index].coordinates.y
                        z = z + self.AllAtoms[index].coordinates.z
        x = x / count
        y = y / count
        z = z / count

        # now, move the pdb to the origin
        self.TranslateMolecule(-x, -y, -z)

        # now rotate
        sinx = math.sin(thetax)
        siny = math.sin(thetay)
        sinz = math.sin(thetaz)
        cosx = math.cos(thetax)
        cosy = math.cos(thetay)
        cosz = math.cos(thetaz)

        cosy_cosz = cosy * cosz
        sinx_siny_cosz_plus_cosx_sinz = sinx * siny * cosz + cosx * sinz
        sinx_sinz_minus_cosx_siny_cosz = sinx * sinz - cosx * siny * cosz
        cosy_sinz = cosy * sinz
        cosx_cosz_minus_sinx_siny_sinz = cosx * cosz - sinx * siny * sinz
        cosx_siny_sinz_plus_sinx_cosz = cosx * siny * sinz + sinx * cosz
        sinx_cosy = sinx * cosy
        cosx_cosy = cosx * cosy


        for atomindex in self.AllAtoms:
            vector = self.AllAtoms[atomindex].coordinates

            new_x = vector.x * cosy_cosz + vector.y * sinx_siny_cosz_plus_cosx_sinz + vector.z * sinx_sinz_minus_cosx_siny_cosz
            new_y = -vector.x * cosy_sinz + vector.y * cosx_cosz_minus_sinx_siny_sinz + vector.z * cosx_siny_sinz_plus_sinx_cosz
            new_z = vector.x * siny - vector.y * sinx_cosy + vector.z * cosx_cosy

            self.AllAtoms[atomindex].coordinates = point(new_x, new_y, new_z)

        # now move it back from the origin
        self.TranslateMolecule(x, y, z)

    # Same as with the atom class
    def SetUndoPoint(self): # you can restore all atom positions to some undo point. This sets that point.
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].SetUndoPoint()

    def Undo(self):
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].Undo()

#----------------------- pqr ---------------------------

class PQR(FileFormator):
	# PQR format
	pass

#----------------------- sdf ---------------------------

class SDF(FileFormator):
	# sdf format
	pass

#----------------------- gjf ---------------------------

class GJF(FileFormator):
	# gjf format
	pass

#----------------------- xyzr --------------------------

class XYZR(FileFormator):
	# xyzr format
	pass

#----------------------- xyz ---------------------------

class XYZ(FileFormator):
	# xyz format
	pass

if __name__=="__main__":
    mr=MOL2()
    a=mr.ReadMolFile("test.mol2");
    print a
    print a.atoms[0]


