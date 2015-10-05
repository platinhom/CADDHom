from copy import deepcopy


class Conect:
    ##The connect information of molecule.Nothing to initial.To copy this object,use deepcopy. Need:
    ##from copy import deepcopy
    def __init__(self):
        self.conects=[]
        self.atom_conects={}
        self.atoms=[]
        self.ends=[]
        self.ring_atoms=[]
        self.rings=[]
        
    def allatoms(self):
        allatoms=[]
        for i in self.conects:
            allatoms.append(i[0])
        self.atoms=allatoms
        return allatoms

    def atomconnects(self):
        atomsconects={}
        for i in self.conects:
            atomsconects[i[0]]=i[1:]
        self.atom_conects=deepcopy(atomsconects)
        return atomsconects
        
    def end_atoms(self):
        endatoms=[]
        for conect in self.conects:
            if len(conect)==2:
                endatoms.append(conect[0])
        self.ends=endatoms
        return endatoms

    ##refresh the bonded atoms
    def refresh(self):
        conects_copy=deepcopy(self.atom_conects)
        for i in conects_copy:
            if len(conects_copy[i])==0:
                del self.atom_conects[i]
                self.atoms.remove(i)
        return
    
    ##refit the atoms,end points,and the atom_conect attributes from conects.
    def refit(self):
        self.atoms=self.allatoms()
        self.ends=self.end_atoms()
        self.atom_conects=self.atomconnects()

    ##cut a bond, offer two atoms.    
    def cut(self,atom1,atom2):
        self.atom_conects[atom1].remove(atom2)
        self.atom_conects[atom2].remove(atom1)

    ##from start point to another point,remove the bond.
    ##Here,num is the index to decide the move path
    def move_bond(self,startp,num=0):
        reach=self.atom_conects[startp][num]
#        if reach not in self.atom_conects:return ('Ring',startp,reach)
        self.cut(reach,startp)#remove the bond
        return reach
    
    ##Further del atoms from start point(an end atom) until reach branch atom 
    def remove_line(self,startp):
        if len(self.atom_conects[startp])!=1:return
        while True:
            if len(self.atom_conects[startp])==1:
                reach=self.move_bond(startp,0)##move to and return next atom
                startp=reach
            else:break
        return
    
    ##remove the sidechain, conserve the rings and their linker
    def remove_sidechain(self):
        for end_atom in self.ends:
            self.remove_line(end_atom)
        self.refresh()
        self.ends=[]

    ##To find the ring atoms and return some rings.
    def find_ring_atoms(self,startpoint=''):
        if self.ends!=[]:self.remove_sidechain()##remove sidechain 
        if startpoint=='':startpoint=self.atoms[0]
        ring=[]
        ringstart=''
        ringend=''
        ring_switch=False
        save_ring=False
        restore=self.atom_conects
        connect=deepcopy(restore)
        startp=startpoint
        breakpoint=[startp,self.atom_conects[startp][0]]##(point_stand,point_goto)
        while len(restore[startpoint])!=0:
            reach=connect[startp][0]##move and delete start point. Key for judgement of ring.
            del connect[startp]
            if reach not in connect: reach=('Ring',startp,reach)
            else: connect[reach].remove(startp)
                        
            if startp==ringstart:##save the rings and ring atoms.
                if ring_switch==True:
                    save_ring=True
                    #print 'now save ring!'
            if save_ring==True:
                ring.append(startp)
            if reach==ringend:
                ring.append(reach)
                #print 'ring is',ring
                self.rings.append(ring)
                for i in ring:
                    if i not in self.ring_atoms:self.ring_atoms.append(i)
                ring=[]
                ring_switch=False
                save_ring=False
                
            #print 'startp',startp,'reach',reach
            #when find a ring.restart the start point and save the ring information.
            if 'Ring' in reach:
                self.cut(reach[1],reach[2])
                print restore[reach[1]],restore[reach[2]]
                ring_switch=True
                ring=[]
                ringstart=reach[2]
                ringend=reach[1]
                #print 'start save ring, ringstart is',ringstart,'ringend is',ringend
                breakpoint=[startpoint,restore[startpoint][0]]
                connect=deepcopy(restore)
                reach=startpoint#
                startp=reach#    
                continue

            ##If find a end point of the chain.Breakpoint is use to cut the sidechain.
            elif len(connect[reach])==0:
                restore[breakpoint[1]].remove(breakpoint[0])
                restore[breakpoint[0]].remove(breakpoint[1])
                if len(restore[startpoint])==0:break
                breakpoint=[startpoint,restore[startpoint][0]]
                connect=deepcopy(restore)
                reach=startpoint#
                startp=reach#
                continue
            
            if len(connect[reach])>=2:
                breakpoint=[reach,connect[reach][0]]       
            print restore[startp], restore[reach]
            startp=reach
        self.refit()

        
    ##read the CONECT information from pdb file. Need to give the file name.
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
        self.refit()
        return conects

    
f='C:/Users/Hom/Desktop/ligand_ring.pdb'
a=Conect()
a.read_conects_PDB(f)
a.remove_sidechain()
len(a.atoms)
len(a.atom_conects)
