from copy import deepcopy
ligand='''CONECT    1   11    5    2
CONECT    2   38   39    1    3
CONECT    3    2    4
CONECT    4   25   15    3
CONECT    5   40    1    6
CONECT    6    9    5    7
CONECT    7   42    6    8
CONECT    8    7   41
CONECT    9   43   10    6
CONECT   10   44   11    9
CONECT   11   45   10    1
CONECT   12   21   17   13
CONECT   13   18   14   12
CONECT   14   46   13   15
CONECT   15   16   14    4
CONECT   16   47   17   15
CONECT   17   16   35   12
CONECT   18   48   19   13
CONECT   19   49   20   18
CONECT   20   50   19   21
CONECT   21   51   20   12
CONECT   22   52   53   30   23
CONECT   23   54   22   29   24
CONECT   24   55   56   23   25
CONECT   25   57   24   26    4
CONECT   26   58   59   27   25
CONECT   27   60   28   30   26
CONECT   28   61   62   29   27
CONECT   29   63   64   28   23
CONECT   30   65   66   22   27
CONECT   31   33   37   32
CONECT   32   31
CONECT   33   67   34   31
CONECT   34   68   33   35
CONECT   35   34   36   17
CONECT   36   69   37   35
CONECT   37   70   36   31
CONECT   38    2
CONECT   39    2
CONECT   40    5
CONECT   41    8
CONECT   42    7
CONECT   43    9
CONECT   44   10
CONECT   45   11
CONECT   46   14
CONECT   47   16
CONECT   48   18
CONECT   49   19
CONECT   50   20
CONECT   51   21
CONECT   52   22
CONECT   53   22
CONECT   54   23
CONECT   55   24
CONECT   56   24
CONECT   57   25
CONECT   58   26
CONECT   59   26
CONECT   60   27
CONECT   61   28
CONECT   62   28
CONECT   63   29
CONECT   64   29
CONECT   65   30
CONECT   66   30
CONECT   67   33
CONECT   68   34
CONECT   69   36
CONECT   70   37'''
class Conect:
    def __init__(self,connects):
        self.conect=deepcopy(connects)
    def cut(self,atom1,atom2):
        self.conect[atom1].remove(atom2)
        self.conect[atom2].remove(atom1)

class Ring:
    def __init__(self,ringatoms):
        self.atoms=tuple(ringatoms)
        self.fused=False
        self.fused_atoms=[]

class Fragment:
    pass

conects=ligand.replace('CONECT','').split('\n')
atoms=[]
atom_connect={}
end_atoms=[]
for i in conects:
    conect=i.split()
    atom=conect[0]
    atoms.append(atom)
    atom_connect[atom]=conect[1:]
    if len(conect)==2:end_atoms.append(atom)

connect_copy=deepcopy(atom_connect)
atoms_nosc=atoms[:]
remove_atoms=[] #collect the non-ring atoms 

def next_atom(startp,connects):
    reach=connects[startp][0]
    del connects[startp]
    if reach not in connects:return ('ring',startp,reach)
    connects[reach].remove(startp)
    return reach

def remove_line(startp,connects,atoms): ##remove the sidechain and change the ligand.
    if len(connects[startp])!=1:return
    while True:
        if len(connects[startp])==1:
            reach=next_atom(startp,connects)##move to next atom
            atoms.remove(startp)
            remove_atoms.append(startp)
            startp=reach
        else:break
    return
        
for end_atom in end_atoms:remove_line(end_atom,connect_copy,atoms_nosc)
print 'atoms_nosc has',len(atoms_nosc)
print 'connect_copy has',len(connect_copy)

class branch:
    def __init__(self,point,come,go):
        self.point=point
        self.come=come
        self.go=go

def restart_ring_loop():
    breakpoint=branch(startpoint,'None',connect_restore[startpoint][0])
    connect=deepcopy(connect_restore)
    reach=startpoint#
    startp=reach#


startpoint=atoms_nosc[2]
print 'startpoint is '+startpoint
rings=[]
ring_atoms=[]
path=[]
ringstart=''
ringend=''
ring_switch=False
save_ring=False
connect_restore=deepcopy(connect_copy)
connect=deepcopy(connect_copy)
startp=startpoint
breakpoint=branch(startp,'None',connect[startp][0])

while len(connect_restore[startpoint])!=0:
    reach=next_atom(startp,connect)
    if reach==ringstart:
        if ring_switch==True:
            save_ring=True
##            print 'now save ring!'
    if save_ring==True:
        path.append(reach)
    if reach==ringend:
##        print 'ring is',path
        rings.append(path)
        for i in path:
            if i not in ring_atoms:ring_atoms.append(i)
        path=[]
        ring_switch=False
        save_ring=False
            
##    print 'startp',startp,'reach',reach#
    if 'ring' in reach:
        connect_restore[reach[1]].remove(reach[2])
        connect_restore[reach[2]].remove(reach[1])
##        print connect_restore[reach[1]],connect_restore[reach[2]]
        ring_switch=True
        path=[]
        ringstart=reach[2]
        ringend=reach[1]
##        print 'start save ring, ringstart is',ringstart,'ringend is',ringend
        breakpoint=branch(startpoint,'None',connect_restore[startpoint][0])
        connect=deepcopy(connect_restore)
        reach=startpoint#
        startp=reach#    
        continue
    
    elif len(connect[reach])==0:
        connect_restore[breakpoint.go].remove(breakpoint.point)
        connect_restore[breakpoint.point].remove(breakpoint.go)
        if len(connect_restore[startpoint])==0:break
        breakpoint=branch(startpoint,'None',connect_restore[startpoint][0])
        connect=deepcopy(connect_restore)
        reach=startpoint#
        startp=reach#
        continue
    
    if len(connect[reach])>=2:
        breakpoint=branch(reach,startp,connect[reach][0])
        
##    print connect_restore[startp], connect_restore[reach]
    startp=reach

individual_ring=[]
ringstemp=rings[:]
ringcontact={}
rings_contact={}
rings_array={}
ringlargest=0
for i in rings:
    if len(i)>ringlargest:ringlargest=len(i)

for i in rings:
    countring=0
    count=0
    for j in i:
        for k in rings:
            if k!=i:
                for kn in k:
                    if kn==j:count+=1
                if count>0:
                    if tuple(i) not in rings_contact:rings_contact[tuple(i)]=[]
                    if k not in rings_contact[tuple(i)]:rings_contact[tuple(i)].append(k)
                    countring+=1
            count=0
    if countring not in ringcontact:ringcontact[countring]=[]
    ringcontact[countring].append(i)
rings_array=deepcopy(rings_contact)
for i in rings_contact:
    for j in rings_contact[i]:
        for k in rings_contact[tuple(j)]:
            if k not in rings_contact[i]:rings_array[i].append(k)

rings_array2={}
rings_array3=[]
rings_array4={}
ringcount=0
for i in rings_array:
    if list(i) not in rings_array3:
        rings_array4[ringcount]=rings_array[i]
        rings_array2[i]=rings_array[i]
        ringcount+=1
    for j in rings_array[i]:
        if j not in rings_array3:
            rings_array3.append(j)

fuse_atoms={}
for i in rings_array4:
    fuse_atoms[i]=[]
    for j in rings_array4[i]:
        indexj=rings_array4[i].index(j)
        for jn in rings_array4[i][indexj]:
            for k in rings_array4[i]:
                indexk=rings_array4[i].index(k)
                if k!=j:
                    for kn in rings_array4[i][indexk]:
                        if kn==jn:
                            if jn not in fuse_atoms[i]:fuse_atoms[i].append(jn)

tempconnect={}
for i in rings[2]:
    tempconnect[i]=connect_copy[i]
    for j in connect_copy[i]:
        if j not in rings[2]:tempconnect[i].remove(j)
