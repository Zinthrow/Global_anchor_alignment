"""
Created on Sat Nov  4 17:01:01 2017

@author: alex
"""
import os
import numpy as np

def S(i,j):
    xi = m[i-1] # i element of m-mer
    yj = n[j-1] # j element of n-mer
    if xi == yj:
        return p1
    else:
        return p2

def Ds(i,j):  #score at (i,j) coordinates for diagonal movement
    return max(D[j-1][i-1] + S(i,j),
               Ix[j-1][i-1] + S(i,j),
               Iy[j-1][i-1] + S(i,j))
    
def Iys(i,j): #score at (i,j) coordinates for vertical movement
    return max(D[j-1][i] + g + s,
               Iy[j-1][i] + s)

def Ixs(i,j): #score at (i,j) coordinates for horizonatal movement
    return max(D[j][i-1] + g + s,
               Ix[j][i-1] + s)
                        
def initialize():
    for indy, y in enumerate(n):
        for indx, x in enumerate(m):
            if indy > 0 and indx > 0:
                D[indy][indx] = Ds(indx,indy)
                Ix[indy][indx] = Ixs(indx,indy)
                Iy[indy][indx] = Iys(indx,indy)
            elif indx is 0 and indy is 0:
                D[0][0] = 0
                Ix[0][0] = g
                Iy[0][0] = g
            elif indy is 0 and indx > 0:
                Ix[indy][indx] = Ixs(indx,indy)
            elif indx is 0 and indy > 0:
                Iy[indy][indx] = Iys(indx,indy)
               
            if indx in anchor_pairs:
                if anchor_pairs[indx] == indy:
                    D[indy][indx] = D[indy][indx] + p1
                    Ix[indy][indx] = Ix[indy][indx] + p1
                    Iy[indy][indx] = Iy[indy][indx] + p1
                    print ("tripped")
def box_max(i,j):
    return  max(D[j][i],
                Ix[j][i],
                Iy[j][i])
               
def align():
    m_align = ""
    n_align = ""
    indx = m_len-1
    indy = n_len-1
    #current = 0
    if box_max(indx,indy) == D[indy][indx]:
        m_align = m[indx] + m_align
        n_align = n[indy] + n_align
    elif box_max(indx,indy) == Ix[indy][indx]:
        m_align = m[indx] + m_align
        n_align = "-" + n_align
    elif box_max(indx,indy) == Iy[indy][indx]:
        m_align = "-" + m_align
        n_align = n[indy] + n_align
    while indx !=0 and indy != 0:
        #if box_max(indx,indy) == D
        max_dir = max(box_max(indx-1,indy-1),
                      box_max(indx-1,indy),
                      box_max(indx,indy-1))
        if max_dir == box_max(indx-1,indy-1):
            m_align = m[indx-1] + m_align
            n_align = n[indy-1] + n_align
            indx = indx-1
            indy = indy-1 
        elif max_dir == box_max(indx-1,indy):
            m_align = m[indx-1] + m_align
            n_align = "-" + n_align
            indx = indx-1
        elif max_dir == box_max(indx,indy-1):
            m_align = "-" + m_align
            n_align = n[indy-1] + n_align
            indy = indy-1
        
    print (m_align)
    print (n_align)       
filename = "input2.txt"
if os.path.exists(filename):
    lines = open(filename, 'r')
    lines = list(lines)
    
parameters = []
for char in lines[0]:
    if char.isnumeric() == True:
        parameters.append(int(char))
        
k = parameters[0] # number of anchor pairs
p1 = parameters[1] #match score
p2 = -parameters[2] #mismatch score
g = -parameters[3] #gap penalty
s = -parameters[4] #extending gap penalty
m = lines[1][:-1] # m-mer sequence of length m : used as the X-axis 
m_len = len(m)          
n = lines[2][:-1] # n-mer sequence of length n : used as the Y-axis
n_len = len(n)
inf = float('inf')

anchor_pairs = {}
if k == 0:
    g = 0
elif k >= 1:
    s = s*k 
    for anchor in range(k):
        pair = lines[anchor+3][:-1].split()
        anchor_pairs[int(pair[0])] = int(pair[1])
    

D = np.array([[-inf]*m_len]*n_len) #diagonal movement
Ix = np.array([[-inf]*m_len]*n_len) #horizontal movement
Iy = np.array([[-inf]*m_len]*n_len) #vertical movement

initialize()

print (D)
print (Ix)
print (Iy)
align()

