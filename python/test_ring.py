from copy import deepcopy


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

    ##cut a bond, offer two atoms, even the conects.
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

    ##To find the fused ring fragment and fused atoms.
    def fused_find(self,rings=[]):
        if rings==[]:rings=self.rings[:]
        fused_atoms={}
        fragment_rings=[]
        fused_last=False
        for i in rings[:-1]:
            fuse=False#to save the longly ring
            for j in rings[rings.index(i)+1:]:
                fused=list(set(i)&set(j))#fused or not                
                if fused!=[]:
                    fuse=True
                    if j==rings[-1]:fuse_last=True#last ring is not included in i.
                    #Follow to save the fused fragment. Sometimes,some framents contains more than two rings.
                    if list(set(i)|set(j)) not in fragment_rings:fragment_rings.append(list(set(i)|set(j)))
                    #follow to save the fused ring information
                    if tuple(i) not in self.fused_rings:self.fused_rings[tuple(i)]=[]
                    self.fused_rings[tuple(i)].append(j)
                    if tuple(j) not in self.fused_rings:self.fused_rings[tuple(j)]=[]
                    self.fused_rings[tuple(j)].append(i)
                    for k in fused: #save the fused atoms informations
                        if k not in fused_atoms:fused_atoms[k]=[]
                        fused_atoms[k].append((i,j))
                elif tuple(i) in self.fused_rings:fuse=True #to avoid false judgement of fused ring to the longly rings
            #the follow to save the longly ring~
            if fuse==False:
                fragment_rings.append(i)
                self.longly_ring.append(i)
        if fuse_last==False:
            fragment_rings.append(rings[-1])
            self.longly_ring.append(ring[-1])
            
        #start the loop to merge the preliminary fragment to bigger fragment. Colud the loop judge from above?
        loop=True
        while loop:
            loop=False
            if len(fragment_rings)<=1:break
            further_rings=fragment_rings[:]
            for i in fragment_rings[:-1]:
                for j in fragment_rings[fragment_rings.index(i)+1:]:
                        if list(set(i)&set(j))!=[]:
                                loop=True #need to further loop
                                fused=list(set(i)|set(j))
                                if i in further_rings:further_rings.remove(i)
                                if j in further_rings:further_rings.remove(j)
                                if fused not in further_rings:further_rings.append(fused)
            fragment_rings=further_rings[:]      
        self.fused_atoms=deepcopy(fused_atoms)
        self.fragment_rings=fragment_rings[:]
        
    ##to find the atoms link to the ring
    def find_ring_link_atoms(self,ringatoms=[],conects={}):
        if ringatoms==[]:ringatoms=self.ring_atoms
        if conects=={}:conects=self.atom_conects
        for ra in ringatoms:
            for link in conects[ra]:
                if (link not in ringatoms) and (link not in self.ring_link_atoms):#if not for Conect.self, the ringlinker maybe error
                    self.ring_link_atoms.append(link)

    ##To find the ring atoms and return some rings.
    def find_ring_atoms(self,startpoint=''):
        if self.ends!=[]:self.remove_sidechain()##remove sidechain
        if len(self.atoms)==0:return ##It's a chain molecule!
        if startpoint=='':startpoint=self.atoms[0]
        ringcore=[]
        for i in self.atom_conects: ##save the main core atoms.
            ringcore.append(i)
        ring=[]
        ringstart=''
        ringend=''
        ring_switch=False
        save_ring=False
        restore=self.atom_conects
        connect=deepcopy(restore)
        startp=startpoint
        breakpoint=[startp,self.atom_conects[startp][0]]#(point_stand,point_goto)
        while len(restore[startpoint])!=0:
            reach=connect[startp][0]#move and delete start point. Key for judgement of ring.
            del connect[startp]
            if reach not in connect: reach=('Ring',startp,reach)
            else: connect[reach].remove(startp)

            if startp==ringstart:#save the rings and ring atoms.
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
                #print restore[reach[1]],restore[reach[2]]
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

            #If find a end point of the chain.Breakpoint is use to cut the sidechain.
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
            #print restore[startp], restore[reach]
            startp=reach
        self.ring_linkers=list(set(ringcore)-set(self.ring_atoms))
        self.refit()#refit the bond,atoms and ends.
        self.sidechain_atoms=list(set(self.atoms)-set(ringcore))
        self.find_ring_link_atoms()#find the atoms link to the rings 
        self.fused_find()#find the fused information

    def removelink(self,conects={}):
        if conects=={}:conects=self.atom_conects
        for i in conects:
            for j in range(0,conects[i].count('link')):
                conects[i].remove('link')

##    def fragmentation(self):
##        conects=deepcopy(self.atom_conects)
####        fragment_end_atoms={}
####        fragments=[]
##        fragment_temp=[]
##        fragment_lib=[]
##        processed_atoms=set([])
####        joint_temp={}
##        ##fragmentation of the sidechain
##        for i in self.ends:
####            fragment_end_atoms[i]=[i]
##            count=1 #include the start point
##            startp=i
##            while True:
##                print startp,count
##                reach=self.atom_conects[startp][0]
##                fragment_temp.append(startp)
##                if count==3:
####                    fragments.append(fragment_end_atoms[i])
##                    processed_atoms.update(fragment_temp)
##                    self.cut(startp,reach,conects)
##                    fragment_lib.append(fragment_temp)
##                    fragment_temp=[]
##                    count=0
##                if reach not in self.ring_atoms:
##                    self.cut(reach,startp)#remove the bond
####                    fragment_end_atoms[i].append(reach)#
##                    if len(self.atom_conects[reach])!=1:
####                        if reach not in joint_temp:joint_temp[reach]=set([])
####                        joint_temp[reach].update(fragment_end_atoms[i])
##                        if startp in conects[reach]:#the situation:count=3 before.
##                            self.cut(startp,reach,conects)
##                            fragment_lib.append(fragment_temp)
##                            fragment_temp=[]
##                        break
##                    startp=reach
##                    count=count+1
####                    print fragment_end_atoms
##                else:
##                    self.cut(startp,reach,conects)
####                    if reach not in joint_temp:joint_temp[reach]=set([reach])
####                    joint_temp[reach].update(fragment_end_atoms[i])
##                    fragment_lib.append(fragment_temp)
##                    fragment_temp=[]
##                    break
####            print 'temp is',joint_temp
####            self.fragment.append(fragment_end_atoms[i])
##            self.atom_conects=deepcopy(conects)
##        print processed_atoms
##        print 'fragment_lib is',fragment_lib
##        return conects

    def fragmentation(self):
        conects=deepcopy(self.atom_conects)
        conects_origin=deepcopy(self.atom_conects)#as reference
        reconects={}
        for i in self.atoms:reconects[i]=[]
        fragment_temp=[]
        fragment_lib=[]
        processed_atoms=set([])
        for i in self.ring_link_atoms:
            for j in conects_origin[i]:
                if j in self.ring_atoms:
                    conects[i].insert(conects[i].index(j),'rlink')
                    conects[j].insert(conects[j].index(i),'rlink')
                    self.cut(i,j,conects)
        self.atom_conects=deepcopy(conects)
        for i in self.ends:
            count=1 #include the start point
            startp=i
            while True:
                print startp,count
                reach=self.atom_conects[startp][0]
                if reach=='rlink':
                    fragment_temp=[]
                    break
                fragment_temp.append(startp)
                self.cut(startp,reach)
                if count==3:
                    processed_atoms.update(fragment_temp)
                    fragment_lib.append(fragment_temp)
                    self.rebond(fragment_temp[0],fragment_temp[1],reconects)
                    self.rebond(fragment_temp[2],fragment_temp[1],reconects)
                    conects[reach].insert(conects[reach].index(startp),'link')
                    self.cut(startp,reach,conects)
                    for i in fragment_temp:del conects[i]
                    fragment_temp=[]
                    count=0
                if len(self.atom_conects[reach])>1:
                    fragment_temp=[]
                    break
                else:
                    startp=reach
                    count+=1
            self.atom_conects=deepcopy(conects_origin)


                    
                
                    
                    


##                if reach not in self.ring_atoms:
##                    self.cut(reach,startp)#remove the bond
####                    fragment_end_atoms[i].append(reach)#
##                    if len(self.atom_conects[reach])!=1:
####                        if reach not in joint_temp:joint_temp[reach]=set([])
####                        joint_temp[reach].update(fragment_end_atoms[i])
##                        if startp in conects[reach]:#the situation:count=3 before.
##                            self.cut(startp,reach,conects)
##                            fragment_lib.append(fragment_temp)
##                            fragment_temp=[]
##                        break
##                    startp=reach
##                    count=count+1
####                    print fragment_end_atoms
##                else:
##                    self.cut(startp,reach,conects)
####                    if reach not in joint_temp:joint_temp[reach]=set([reach])
####                    joint_temp[reach].update(fragment_end_atoms[i])
##                    fragment_lib.append(fragment_temp)
##                    fragment_temp=[]
##                    break
####            print 'temp is',joint_temp
####            self.fragment.append(fragment_end_atoms[i])
##            self.atom_conects=deepcopy(conects)
##        print processed_atoms
##        print 'fragment_lib is',fragment_lib
        return (conects,reconects)

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


f='D:\Documents and Settings\zhaozx\Desktop/ligand_ring2.pdb'
a=Conect()
a.read_conects_PDB(f)
a.remove_sidechain()
b=len(a.atoms)
c=len(a.atom_conects)
a.find_ring_atoms()
a.refresh()
hi=a.fragmentation()
print 'OK!'
