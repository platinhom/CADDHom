#! /usr/bin/env python
import math, sys, os , platform
from collections import Iterable
from optparse import OptionParser 
import pybel2
import openbabel as op
import numpy as np
# http://openbabel.org/dev-api/classOpenBabel_1_1OBMol.shtml
 
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
			newline=line
			# Bug in Openbabel 2.3.2 windows version
			if (platform.system()=="Windows"):
				newline=newline[:-3]+'\n';
			# Br ->B..
			if (line[12:15]==" Br"): newline=newline[:12]+"Br "+newline[15:] 
			if (line[12:15]==" Cl"): newline=newline[:12]+"Cl "+newline[15:] 
			if (line[12:15]==" Na"): newline=newline[:12]+"Na "+newline[15:] 
			if (line[12:15]==" Mg"): newline=newline[:12]+"Mg "+newline[15:] 
			out+=newline
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
	# Return bond atoms atomic number pair, such as C-O return (6,7)
	return tuple(sorted([bond.bgn.atomicnum,bond.end.atomicnum]))

def atomNumHyd(atom):
	# Return Atom's (atomic number, hydribazation) pair.
	return (atom.atomicnum,atom.hyb)

def distance(atom1,atom2):
	# Return distance between two atoms.
	return atom1.OBAtom.GetDistance(atom2.OBAtom);

def MolInfo(mol,printInfo=True):
	# Ruturn a list containing molecular information/features
	smile=mol.write('smi').strip().split()[0]
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

	if (printInfo):
		print "Mol Formula:",mol.formula
		print "Mol Weight:",mol.molwt
		print "Mol SMILE:",smile
		print "Mol dipole:",dipole;
		print "Total Atoms number:", TNatms;
		print "Heavy Atom number:", HEatms;
		print "Hydrogen number:", Hatms;
		print "Bond number:", TNbnds;
		print "Single Bond number:",sbnds
		print "Double Bond number:",dbnds
		print "Triple Bond number:",tbnds
		print "Aromatic Bond number:",abnds

	return [mol.formula,mol.molwt,smile,dipole,TNatms,HEatms,Hatms,TNbnds,sbnds,dbnds,tbnds,abnds]

def descVar(*args):
	# Return [max, min, sum, average, std] for given data
	if (isinstance(args[0],list)):
		args=args[0]
	mx=max(args)
	mi=min(args)
	sumall=math.fsum(args)
	aver=sumall/len(args)
	var= math.fsum((pow(x-aver,2) for x in args)) /(len(args))
	std=math.sqrt(var)
	return (mx,mi,sumall,aver,std)

def featureDict2List(ftype, fdict):
	features=[]
	for f in ftype:
		if (isinstance(fdict[f],Iterable)):
			features+=list(fdict[f])
		else:
			features.append(fdict[f])
	return features

def CalcFeatures(mol,printInfo=True):
	atoms=mol.atoms;
	atomshyb=[atomNumHyd(atom) for atom in atoms];
	bonds=mol.bonds;
	# H,C,N,O,F,P,S,Cl,Br,I
	elements=[1,6,7,8,9,15,16,17,35,53]

	# molecule 
	molinfo=MolInfo(mol,printInfo=printInfo);
	elecounts={}
	for ele in elements:
		elecounts[ele]=0
	for atom in mol:
		an=atom.atomicnum
		elecounts[an]=elecounts.get(an,0)+1

	# partial charge
	acDict={};
	acDictAbs={};
	pcs=[atom.partialcharge for atom in mol]
	pcsAbs=[abs(pc) for pc in pcs]
	pcdict={}
	pcdesc={}
	pcAbsdict={}
	pcAbsdesc={}
	for ele in elements:
		pcdict[ele]=[]
		pcAbsdict[ele]=[]
		acDict[ele]=0.0
		acDictAbs[ele]=0.0
	for atom in mol:
		an=atom.atomicnum
		pcdict[an].append(atom.partialcharge)
		pcAbsdict[an].append(abs(atom.partialcharge))
		acDict[an]=acDict.get(an,0.0)+atom.partialcharge;
		acDictAbs[an]=acDictAbs.get(an,0.0)+abs(atom.partialcharge);
	for ele in elements:
		if (len(pcdict[ele])>0):
			pcdesc[ele]=descVar(pcdict[ele])
			pcAbsdesc[ele]=descVar(pcAbsdict[ele])
		else:
			pcdesc[ele]=descVar([0.0])
			pcAbsdesc[ele]=descVar([0.0])
	elePCfeatures=featureDict2List(elements,acDict)+featureDict2List(elements,acDictAbs)+featureDict2List(elements,pcdesc)+featureDict2List(elements,pcAbsdesc)
	if (printInfo):
		# mol feature
		print "Partial Charge Max, Min, Sum, Average, Std:",descVar(pcs)
		print "Abs Partial Charge Max, Min, Sum, Average, Std:",descVar(pcsAbs)
		# element feature
		print "Element Partial Charge:",acDict
		print "Element Abs Partial Charge:",acDictAbs
		print "Element Partial Charge Max, Min, Sum, Average, Std:",pcdesc
		print "Abs Element Partial Charge Max, Min, Sum, Average, Std:",pcAbsdesc

	# hybridization
	EleHyb={}
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
	
	hybtypes=["C1","C2","C3","N1","N2","N3","O1","O2","O3","S1","S2","S3"];
	hybcount={}
	for t in hybtypes:
		hybcount[t]=EleHyb.get(t,0)
	hybfeatures=featureDict2List(hybtypes,hybcount)
	if (printInfo):
		#print "Element Hybridization:",EleHyb
		print "Element Hybridization count:",hybcount

	# dipole
	dipoles= [calcdipoleBond(bond) for bond in mol.bonds]
	bndpair= [atomnumBondPair(bond) for bond in mol.bonds ]
	bpd={}
	for i in range(len(dipoles)):
		dp=dipoles[i]
		bp=bndpair[i]
		if not bpd.has_key(bp):
			bpd[bp]=[]
		bpd[bp].append(dp)
	bpneed=[(1,6),(1,7),(1,8),(1,16),(6,6),(6,7),(6,8),(6,9),(6,15),(6,16),(6,17),(6,35),(6,53),
			(7,8),(8,15),(8,16),(15,16),(16,16)]
	bpddesc={}
	for bpn in bpneed:
		bpddesc[bpn]=descVar(bpd.get(bpn,[0.0]))
	dipolefeautures=featureDict2List(bpneed,bpddesc)
	if printInfo:
		# mol feature
		print "Bond Dipoles Max, Min, Sum, Average, Std:", descVar(dipoles)
		# element feature
		print "Bond Atom Pair Dipoles Max, Min, Sum, Average, Std:",bpddesc

	# Mol Formula, Mol Weight, Mol SMILE, Mol dipole, Total Atoms number, Heavy Atom number, Hydrogen number, 
	# Bond number, Single Bond number, Double Bond number, Triple Bond number, Aromatic Bond number
	# For H,C,N,O,F,P,S,Cl,Br,I
	#   Element number
	# Partial Charge, Abs Partial Charge, Bond Dipoles: Max, Min, Sum, Average, Std
	# 
	# For H,C,N,O,F,P,S,Cl,Br,I
	# 	Element Partial Charge 
	# 	Element Abs Partial Charge 
	# 	Element Partial Charge Max, Min, Sum, Average, Std 
	# 	Abs Element Partial Charge Max, Min, Sum, Average, Std
	#
	# For ["C1","C2","C3","N1","N2","N3","O1","O2","O3","S1","S2","S3"]
	# 	Element Hybridization count
	#
	# For [(1,6),(1,7),(1,8),(1,16),(6,6),(6,7),(6,8),(6,9),(6,15),(6,16),(6,17),(6,35),(6,53),(7,8),(8,15),(8,16),(15,16),(16,16)]
	# For HC,HN,HO,HS,CC,CN,CO,CF,CP,CS,CCl,CBr,CI,NO,OP,OS,PO,SS
	# 	Bond Atom Pair Dipoles Max, Min, Sum, Average, Std
	outlist=molinfo+featureDict2List(elements, elecounts)+list(descVar(pcs))+list(descVar(pcsAbs))+list(descVar(dipoles)) \
			+elePCfeatures+hybfeatures+dipolefeautures
	if printInfo:print outlist
	return [ str(f) for f in outlist ]

def featureString():
	# 12 mol feature
	fstr="Mol_Formula Mol_Weight Mol_SMILE Mol_dipole Total_Atoms_number Heavy_Atom_number Hydrogen_number "
	fstr+="Bond_number Single_Bond_number Double_Bond_number Triple_Bond_number Aromatic_Bond_number "
	# 10 element partial charge feature
	for i in ["H","C","N","O","F","P","S","Cl","Br","I"]:
		fstr+=(i+"_"+"num"+" ")
	# 15 mol feature
	for i in ["PartCharge","AbsPartCharge","Bond_Dipole"]:
		for j in ["Max","Min","Sum","Aver","Std"]:
			fstr+=(i+"_"+j+" ")
	# 120 element partial charge feature
	for i in ["H","C","N","O","F","P","S","Cl","Br","I"]:
		fstr+=(i+"_"+"PC"+" ")
	for i in ["H","C","N","O","F","P","S","Cl","Br","I"]:
		fstr+=(i+"_"+"APC"+" ")
	for i in ["H","C","N","O","F","P","S","Cl","Br","I"]:
		for j in ["Max","Min","Sum","Aver","Std"]:
			fstr+=(i+"_"+"PC"+"_"+j+" ")
	for i in ["H","C","N","O","F","P","S","Cl","Br","I"]:
		for j in ["Max","Min","Sum","Aver","Std"]:
			fstr+=(i+"_"+"APC"+"_"+j+" ")
	# 12 hybrid feature
	for i in ["C1","C2","C3","N1","N2","N3","O1","O2","O3","S1","S2","S3"]:
		fstr+=(i+"_"+"Hyb"+" ")
	# 90 atom pair bond dipole feature
	for i in [ "HC","HN","HO","HS","CC","CN","CO","CF","CP","CS","CCl","CBr","CI","NO","OP","OS","PO","SS"]:
		for j in ["Max","Min","Sum","Aver","Std"]:
			fstr+=(i+"_"+"DP"+"_"+j+" ")
	#print fstr
	return fstr	

if __name__ =="__main__":
	helpdes='''Calculate features of molecules based on Pybel and Openbabel.
	# For one file, use -i option to assign the input file;
	# For many files, use -m option to assign a file containing file name without extension.
	# -f option can assign the file format. It must be given when using -m option. 
	# Without -f option and using -i option, the format will be deduced based on file extension.
	# -t option will print the title for features.'''

	parser = OptionParser(description=helpdes) 
	parser.add_option("-i", "--input", action="store", 
                    dest="input", default="",
                    help="Read input data from input file")
	parser.add_option("-m", "--multi", action="store", 
					dest="multi", default="",
              		help="File containing file name without extension, format must be assigned!")
	parser.add_option("-f", "--format", action="store", 
					dest="format", default="",
              		help="Input file format")
	parser.add_option("-o", "--output", action="store", 
					dest="output", default="",
              		help="The output file to save result")
	parser.add_option("-t", "--title", action="store_true", 
					dest="title", default=False,
              		help="Print the feature title")
	(options, args) = parser.parse_args()
	stdout=sys.stdout
	if (options.output!=""):
		ftmp=open(options.output,'w');
		sys.stdout=ftmp
	if (options.title): print featureString()
	if (options.input != ""):
		filename=options.input
		fnamelist=os.path.splitext(filename)
		fformat=options.format
		if (fformat==""):
			fformat=fnamelist[1][1:].lower()
		if (fformat=="pqr"):
			mol=pybel2.readstring('pqr',pqrbug(filename));
		else:
			mol=pybel2.readfile(fformat,filename).next();
		features=CalcFeatures(mol,printInfo=False)
		print fnamelist[0]+" "+" ".join(features)
	elif (options.multi != "" and options.format != ""):
		fin=open(options.multi)
		flist=fin.readlines()
		fin.close()
		fformat=options.format
		for f in flist:
			try:
				filename=f.strip()+"."+fformat
				if (fformat=="pqr"):
					mol=pybel2.readstring('pqr',pqrbug(filename));
				else:
					mol=pybel2.readfile(fformat,filename).next();
				features=CalcFeatures(mol,printInfo=False)
				print f.strip()+" "+" ".join(features)	
			except IOError:
				print f.strip()		
	else:
		raise ValueError("No input file!")
		exit(1)
	sys.stdout=stdout;
	if (options.output!=""):ftmp.close()