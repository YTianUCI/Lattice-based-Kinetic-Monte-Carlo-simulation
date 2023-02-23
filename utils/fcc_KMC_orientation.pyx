#for fcc crystal
import cython

import numpy as np
import random as rd

# cdef class lattice_site:
#     cdef public int occupy
#     cdef public int number
#     cdef public float position[3]
#     cdef public int neighbor
#     def __init__(self,number,position):
#         self.occupy=0
#         self.number=number
#         self.position=position
#         self.neighbor=0
#         self.neighbor_list=[]
#     def neighborlist_init(self,lattice,rmax):
#         for i in lattice:
#             m=np.array(i.position)-np.array(self.position)
#             if i.number is self.number:continue
#             if abs(m[0])<rmax and abs(m[1])<rmax and abs(m[2])<rmax: 
#                 if (m[0]**2+m[1]**2+m[2]**2)<rmax**2:
#                     self.neighbor_list.append(i.number)

class lattice_site:
    def __init__(self,number,position):
        self.full_neighbor=0
        self.occupy=0
        self.number=number
        self.position=position
        self.neighbor=0
        self.neighbor_list=[]
        self.hopping_time=0
    def neighborlist_init(self,lattice,rmax):
        for i in lattice:
            m=np.array(i.position)-np.array(self.position)
            if i.number==self.number:continue
            if abs(m[0])<rmax and abs(m[1])<rmax and abs(m[2])<rmax: 
                if m[0]**2+m[1]**2+m[2]**2<rmax**2:
                    self.neighbor_list.append(i.number)

cpdef int neighbor_atom_renew(lattice,int number):
    cdef int x=0
    for j in lattice[number].neighbor_list:
        if lattice[j].occupy==1:x=x+1
    lattice[number].neighbor=x
cdef float rate_list[15];
bond_energy=0.1;Kb=8.617e-5;T=400;v=10**13/500*T
def rate(bond):
    return v*np.exp(-bond_energy*bond/Kb/T)
for i in range(15):
    rate_list[i]=rate(i)
cpdef float sums(float a,float b):
    return a+b
cpdef float subs(float a,float b):
    return a-b
cpdef bint judge_surface(bint occupy, int neighbor):
    if occupy and neighbor!=12:return True;
    else: return False;
# cpdef slct_(dict surface, float time, float a=rd.random()):
#     cdef float rate_sum=0
#     for j in surface:
#         rate_sum=sums(rate_sum,surface[j])
#     cdef float Ri=a*rate_sum
#     time=time-np.log(rd.random())/rate_sum
#     cdef int k=-1
#     for j in surface:
#         Ri=subs(Ri,surface[j]);
#         if Ri<=0:return j,time;

cpdef slct(dict surface,float time, float a=rd.random()):
    cdef float rate_sum=0
    cdef int atom_list[10000]
    cdef float rate_list[10000]
    cdef int m=0
    for j in surface:
        rate_sum=sums(rate_sum,surface[j])
        rate_list[m]=rate_sum
        atom_list[m]=j
        m+=1
    cdef float Ri=a*rate_sum
    time=time-np.log(rd.random())/rate_sum
    cdef int i=0
    for i from 0<=i<m:
        if Ri<rate_list[i]:return atom_list[i],time; 


cpdef restricted_hopping(lattice,dict surface,int hopping_atom):
    possible_site=[]
    cdef int hopped_position=-1
    for j in lattice[hopping_atom].neighbor_list:
        if lattice[j].occupy==0 and lattice[j].neighbor!=1 and lattice[j].neighbor!=2:possible_site.append(j)
    cdef int x=len(possible_site)
    if x==0: 
        for j in lattice[hopping_atom].neighbor_list:
            if lattice[j].occupy==0 and lattice[j].neighbor!=1:possible_site.append(j)
    x=len(possible_site)
    if x==0: raise Exception("isolated atom") 
    cdef int neighbor=0    
    for j in possible_site:
        neighbor_=lattice[j].neighbor  
        if neighbor<neighbor_:neighbor=neighbor_;hopped_position=j
    changed_atoms=set(lattice[hopping_atom].neighbor_list+lattice[hopped_position].neighbor_list)
    lattice[hopping_atom].occupy=0
    lattice[hopped_position].occupy=1
    for j in lattice[hopping_atom].neighbor_list:
        lattice[j].neighbor-=1
    for j in lattice[hopped_position].neighbor_list:
        lattice[j].neighbor+=1
    for j in changed_atoms:
        if surface.__contains__(j): surface.pop(j)
        if lattice[j].occupy==1 and lattice[j].neighbor!=lattice[j].full_neighbor:
            surface[j]=rate_list[lattice[j].neighbor]
    if lattice[hopped_position].neighbor==0:print(hopping_atom,hopped_position);raise Exception("hopped_position") 

cpdef free_hopping(lattice,surface,int hopping_atom):
    possible_site=[]
    cdef int hopped_position=-1
    for j in lattice[hopping_atom].neighbor_list:
        if lattice[j].occupy==0 and lattice[j].neighbor!=1 and lattice[j].neighbor!=2:possible_site.append(j)
    cdef int x=len(possible_site)
    if x==0: raise Exception("isolated atom") 
    hopped_position=rd.choice(possible_site);
    changed_atoms=set(lattice[hopping_atom].neighbor_list+lattice[hopped_position].neighbor_list)
    lattice[hopping_atom].occupy=0
    lattice[hopped_position].occupy=1
    for j in lattice[hopping_atom].neighbor_list:
        lattice[j].neighbor-=1
    for j in lattice[hopped_position].neighbor_list:
        lattice[j].neighbor+=1
    for j in changed_atoms:
        if surface.__contains__(j): surface.pop(j)
        if lattice[j].occupy==1 and lattice[j].neighbor!=lattice[j].full_neighbor:
            surface[j]=rate_list[lattice[j].neighbor]
    if lattice[hopped_position].neighbor==0:print(hopping_atom,hopped_position);raise Exception("hopped_position") 

