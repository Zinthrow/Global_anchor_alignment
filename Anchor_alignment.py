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
    
def Ixs(i,j): #score at (i,j) coordinates for horizonatal movement
    return max(D[j-1][i] + g + s,
               Ix[j-1][i] + s)

def Iys(i,j): #score at (i,j) coordinates for horizonatal movement
    return max(D[j][i-1] + g + s,
               Iy[j][i-1] + s)
                        
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
                Iy[indy][indx] = Iys(indx,indy)
            elif indx is 0 and indy > 0:
                Ix[indy][indx] = Ixs(indx, indy)

filename = "input5.txt"
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
m = lines[1][:-1] # m-mer sequence of length m             
n = lines[2][:-1] # n-mer sequence of length n        
inf = float('inf')

D = np.array([[-inf]*len(m)]*len(n)) #diagonal movement
Ix = np.array([[-inf]*len(m)]*len(n)) #horizontal movement
Iy = np.array([[-inf]*len(m)]*len(n)) #vertical movement

initialize()

#np.set_printoptions(threshold=np.inf)
print (Ix)           
       


