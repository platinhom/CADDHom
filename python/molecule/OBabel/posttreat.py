#! /usr/bin/env python
import math
import pybel2

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

def calcdipole(bond):
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

def distance(atom1,atom2):
	return atom1.OBAtom.GetDistance(atom2.OBAtom);

if __name__ =="__main__":
	mol=pybel2.readstring('pqr',pqrbug(filename));