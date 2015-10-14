
# -*- coding: utf-8 -*-
"""
Created on 2015-10-05
@author: Zhixiong Zhao
"""

import __init__
from HHFormat import *
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

if __name__=="__main__":
    mr=MOL2()
    a=mr.ReadMolFile("test.mol2");
    print a
    print a.atoms[0]




