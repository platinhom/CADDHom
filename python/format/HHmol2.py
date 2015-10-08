
# -*- coding: utf-8 -*-
"""
Created on 2015-10-05
@author: Zhixiong Zhao
"""

from copy import deepcopy
import os
import math

class Point:
    x=99999.0
    y=99999.0
    z=99999.0    
    def __init__ (self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z
        
    def coors(self):
        coor=(self.x,self.y,self.z)
        return coor

class Atom:
    def __init__(self):
        self.name=''
        self.element=''        
        self.resname=''
        self.resid=0
        self.chain=''
        self.coordinates = Point(99999, 99999, 99999)
        self.x=self.coordinates.coors()[0]
        self.y=self.coordinates.coors()[1]
        self.z=self.coordinates.coors()[2]
        self.undo_coordinates = Point(99999, 99999, 99999)
        self.line=""
        self.PDBIndex = ""
        self.MOL2Index =""
        self.type=''
        self.charge=''
        self.IndeciesOfAtomsConnecting=[]
        self.in_ring=False

    def ReadMOL2Line(self, Line):
        self.line = Line
        items=Line.split()
        self.MOL2Index = items[0]
        self.atomname = items[1]
        self.chain = ''
        self.coordinates = Point(float(items[2]), float(items[3]), float(items[4]))
        self.type=items[5]
        self.element=self.type[0:2].strip('.')
        if len(items)==9:
            self.resid = int(items[6])
            self.resname = items[7]
            self.charge = items[8]

    def CreateMOL2Line(self, mol2=MOL2()):
        if mol2.len_atom==0:
            len_atom=4
            len_resid=3
            len_resname=6
        else:
            len_atom=mol2.len_atom
            len_resid=mol2.len_resid
            len_resname=mol2.len_resid
        output=self.MOL2Index.rjust(len_atom)+ ' ' +self.atomname.ljust(5)
        output=output+("%.3f" % self.coordinates.x).rjust(10) + ("%.4f" % self.coordinates.y).rjust(11)+ ("%.4f" % self.coordinates.z).rjust(11)+ ' '
        output=output+self.type.ljust(6)+str(self.resid).rjust(len_resid)+ ' ' + self.resname.ljust(len_resname)+ self.charge.rjust(9)+ os.linesep
        return output
        

class MOL2:
    def __init__(self):
        self.AllAtoms={}
        self.Conects=Conect()
        self.heavy_atoms_Conects=Conect()
        self.bonds={}
        len_atom=0
        len_resid=0
        len_resname=0
        
    def LoadMOL2(self, FileName):
        autoindex = 1
        self.__init__()
        # Now load the file into a list
        file = open(FileName,"r")
        lines = file.readlines()
        file.close()
        self.filename=FileName
        conects=[]
        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    TempAtom = Atom()
                    TempAtom.ReadPDBLine(line)
                    self.AllAtoms[autoindex] = TempAtom # because points files have no indecies
                    autoindex = autoindex + 1
            items=line.split()
            if items[0]=='CONECT':
                conects.append(items[1:])
        self.Conects.conects=conects[:]
        self.Conects.restart()
        self.Conects.find_ring_atoms()
        heavyconects=[]
        hydrogens={}
        for i in self.Conects.conects:
            for j in self.AllAtoms:
                if i[0]==self.AllAtoms[j].PDBIndex:
                    self.Conects.atoms_index[i[0]]=j
                    if self.AllAtoms[j].element!='H':
                        heavyconects.append(i)
                        self.heavy_atoms_Conects.atoms_index[i[0]]=j
                    elif self.AllAtoms[j].element=='H':
                        hydrogens[i[0]]=j
        copy_heavy=deepcopy(heavyconects)
        for i in copy_heavy:
            for j in i:
                if j in hydrogens:
                    heavyconects[copy_heavy.index(i)].remove(j)
        self.heavy_atoms_Conects.conects=heavyconects[:]
        self.heavy_atoms_Conects.restart()
        self.heavy_atoms_Conects.find_ring_atoms()
