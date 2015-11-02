#! /usr/bin/env python
import math
import pybel2
import openbabel as op
import numpy as np

def pqrbug(filename):
# Bug in Openbabel for PQR format reader
# Return the string for the file, which can be read by:
#   pybel.readstring('pqr',pqrbug(filename)).
#   ob.OBConversion().ReadString(obmol, string)
# BUG has been removed in Mac version.
	f=open(filename);
	lines=f.readlines();
	out=""
	for line in lines:
		if line[:6]=="ATOM  " or line[:6]=="HETATM":
			#print line[:-3];
			out+=line[:-3]+'\n';
		else: out+=line;
	f.close()
	return out;

Molecule=pybel2.Molecule
Atom=pybel2.Atom
Bond=pybel2.Bond

def calcdipoleAtoms(*atoms):
	# give many atom as input
	# Best method for list of atoms should calcdipoleAtoms(*Atomlist)
	if (len(atoms)<= 0):
		raise TypeError("Errors: No Input Atoms!")
		return 0.0
	# if giving a list of atoms. 
	if (isinstance(atoms[0],list)):
		atoms=atoms[0]
	if (not isinstance(atoms[0],Atom)):
		raise TypeError("Errors: Input should be Atom!")
		return 0.0
	dx=0.0;dy=0.0;dz=0.0
	for atom in atoms:
		coor=atom.coords
		charge=atom.partialcharge
		dx+=coor[0]*charge
		dy+=coor[1]*charge
		dz+=coor[2]*charge
	dipole=math.sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2))
	return dipole

def calcdipoleBond(bond):
	# Give a bond as input
	bgn=bond.bgn
	end=bond.end
	bgncoor=bgn.coords
	bgncharge=bgn.partialcharge
	endcoor=end.coords
	endcharge=end.partialcharge
	dipole=math.sqrt(pow((bgncharge*bgncoor[0]+endcharge*endcoor[0]),2)
			+pow((bgncharge*bgncoor[1]+endcharge*endcoor[1]),2)
			+pow((bgncharge*bgncoor[2]+endcharge*endcoor[2]),2));
	return dipole

def atomnumBondPair(bond):
	return tuple(sorted([bond.bgn.atomicnum,bond.end.atomicnum]))

def atomNumHyd(atom):
	return (atom.atomicnum,atom.hyb)

def distance(atom1,atom2):
	return atom1.OBAtom.GetDistance(atom2.OBAtom);

def MolInfo(mol):
	print "Mol Formula:",mol.formula
	print "Mol Weight:",mol.molwt
	smile=mol.write('smi').strip()
	print "Mol SMILE:",smile
	
	dipole=calcdipoleAtoms(*mol.atoms)
	TNatms=mol.OBMol.NumAtoms();
	HEatms=mol.OBMol.NumHvyAtoms();
	Hatms=mol.OBMol.NumAtoms()-mol.OBMol.NumHvyAtoms();
	TNbnds=mol.OBMol.NumBonds();

	moldesc=mol.calcdesc()
	sbnds=int(moldesc['sbonds']);
	dbnds=int(moldesc['dbonds']);
	tbnds=int(moldesc['tbonds']);
	abnds=int(moldesc['abonds']);

	print "Mol dipole:",dipole;
	print "Total Atoms number:", TNatms;
	print "Heavy Atom number:", HEatms;
	print "Hydrogen number:", Hatms;
	print "Bond number:", TNbnds;
	print "Single Bond number:",sbnds
	print "Double Bond number:",dbnds
	print "Triple Bond number:",tbnds
	print "Aromatic Bond number:",abnds

	return [mol.formula,mol.molwt,smile,TNatms,HEatms,Hatms,TNbnds,sbnds,dbnds,tbnds,abnds]

def descVar(*args):
	# Return [max, min, sum, average, std]
	if (isinstance(args[0],list)):
		args=args[0]
	mx=max(args)
	mi=min(args)
	sumall=math.fsum(args)
	aver=sumall/len(args)
	var= math.fsum((pow(x-aver,2) for x in args)) /(len(args)-1)
	std=math.sqrt(var)
	return (mx,mi,sumall,aver,std)

def CalcFeatures(mol):
	atoms=mol.atoms;
	atomshyb=[atomNumHyd(atom) for atom in atoms];
	bonds=mol.bonds;
	# H,C,N,O,F,P,S,Cl,Br,I
	elements=[1,6,7,8,9,15,16,17,35,53]

	# molecule 
	molinfo=MolInfo(mol);

	# partial charge
	acDict={};
	acDictAbs={};
	pcs=[atom.partialcharge for atom in mol]
	pcsAbs=[abs(pc) for pc in pcs]
	for atom in mol:
		an=atom.atomicnum
		acDict[an]=acDict.get(an,0.0)+atom.partialcharge;
		acDictAbs[an]=acDictAbs.get(an,0.0)+abs(atom.partialcharge);
	print "Element Partial Charge:",acDict
	print "Element Abs Partial Charge:",acDictAbs
	print "Partial Charge Max, Min, Sum, Average, Std:",descVar(pcs)
	print "Abs Partial Charge Max, Min, Sum, Average, Std:",descVar(pcsAbs)

	EleHyb={}
	# hybridization
	for atom in mol:
		ehyb=atomNumHyd(atom)
		if (ehyb[0] is 6):
			if (ehyb[1] is 1):
				EleHyb["C1"]=EleHyb.get("C1",0)+1
			elif (ehyb[1] is 2):
				EleHyb["C2"]=EleHyb.get("C2",0)+1
			elif (ehyb[1] is 3):
				EleHyb["C3"]=EleHyb.get("C3",0)+1
		if (ehyb[0] is 7):
			if (ehyb[1] is 1):
				EleHyb["N1"]=EleHyb.get("N1",0)+1
			elif (ehyb[1] is 2):
				EleHyb["N2"]=EleHyb.get("N2",0)+1
			elif (ehyb[1] is 3):
				EleHyb["N3"]=EleHyb.get("N3",0)+1
		if (ehyb[0] is 8):
			if (ehyb[1] is 1):
				EleHyb["O1"]=EleHyb.get("O1",0)+1
			elif (ehyb[1] is 2):
				EleHyb["O2"]=EleHyb.get("O2",0)+1
			elif (ehyb[1] is 3):
				EleHyb["O3"]=EleHyb.get("O3",0)+1
		if (ehyb[0] is 16):
			if (ehyb[1] is 1):
				EleHyb["S1"]=EleHyb.get("S1",0)+1
			elif (ehyb[1] is 2):
				EleHyb["S2"]=EleHyb.get("S2",0)+1
			elif (ehyb[1] is 3):
				EleHyb["S3"]=EleHyb.get("S3",0)+1
	print "Element Hybridization:",EleHyb

if __name__ =="__main__":
	mol=pybel2.readstring('pqr',pqrbug(filename));