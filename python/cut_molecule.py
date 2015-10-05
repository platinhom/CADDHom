from copy import deepcopy

def sortlen(f1,f2):
    if len(f1)<len(f2):return -1
    else:return 0

class Conect:
    ##The connect information of molecule.Nothing to initial.To copy this object,use deepcopy. Need:
    ##from copy import deepcopy
    def __init__(self):
        self.conects=[]
        self.atom_conects={}#atoms and their conects
        self.atoms=[]
        self.ends=[]#atoms at the end of the sidechain
        self.ring_atoms=[]#atoms in the rings
        self.ring_linkers=[]#atoms link the ring together
        self.ring_link_atoms=[]#the atoms that link to the ring atoms
        self.sidechain_atoms=[]#atoms in the sidechain
        self.rings=[]#search rings based on [0] information,non-chemical sense..
        self.longly_ring=[]#non-fused ring
        self.fused_atoms={}#fused atoms and information,non-chemical sense..
        self.fused_rings={}#fused rings and information,non-chemical sense..
        self.fragment=[]
        self.fragment_rings=[]#the big ring fragment.

    ##restart the atoms,end points,and the atom_conect attributes from conects.
    def restart(self):
        for i in self.conects:
            self.atoms.append(i[0])
            self.atom_conects[i[0]]=i[1:]
            if len(i)==2:
                self.ends.append(i[0])
        
    def read_conects_PDB(self,pdbfile):
        f=open(pdbfile,'r')
        lines=f.readlines()
        f.close()
        conects=[]
        for line in lines:
            items=line.split()
            if items[0]=='CONECT':
                conects.append(items[1:])
        self.conects=conects[:]
        self.restart()
        return conects

    ##cut and refit a bond, offer two atoms, even the conects.
    def cut(self,atom1,atom2,conects=''):
        if conects=='':
            self.atom_conects[atom1].remove(atom2)
            self.atom_conects[atom2].remove(atom1)
        else:
            conects[atom1].remove(atom2)
            conects[atom2].remove(atom1)

    def rebond(self,atom1,atom2,conects=''):
        if conects=='':
            self.atom_conects[atom1].append(atom2)
            self.atom_conects[atom2].append(atom1)
        else:
            conects[atom1].append(atom2)
            conects[atom2].append(atom1)

    ##remove the sidechain, conserve the rings and their linker
    def remove_sidechain(self):
        for end_atom in self.ends:
            startp=end_atom
            while True:
                if len(self.atom_conects[startp])==1:
                    reach=self.atom_conects[startp][0]
                    self.cut(reach,startp)
                    reach=self.move_bond(startp,0)##move to and return next atom
                    startp=reach
                else:break
        self.refresh()
        self.ends=[]
            
f='D:\Documents and Settings\zhaozx\Desktop/ligand_ring2.pdb'
a=Conect()
a.read_conects_PDB(f)
##a.remove_sidechain()
##b=len(a.atoms)
##c=len(a.atom_conects)
##a.find_ring_atoms()
##a.refresh()
##ringlink=a.link_to_ring()
##hi=a.fragmentation()
print 'OK!'
