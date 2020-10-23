#IN THE NAME OF ALLAH
#new algorithm
#read line 

import numpy as np
import time


def linalg_norm(a,b):
    #print("from function")
    return np.linalg.norm(a - b)

# get starting time
start = time.time()

Atom=np.loadtxt('atom.data')
Bond=np.loadtxt('bond.data')
dist_array_TVB=[]

new_Bond_TVB=np.zeros(shape=(3000,4))
for i in range(len(Atom)):
    #print("############atom################## ",i+1)
    for j in range(len(Atom)):
        if Atom[i,-1]==Atom[j,-1] and j>i:
            #print("compute bond between",i+1,j+1)
            dist=linalg_norm(Atom[i,4:],Atom[j,4:])
            #print("the i TVB-bond is",i+1,"dist",dist)
            dist_array_TVB.append(float(dist))
        
            if 3.33<dist<3.35 :
                #print ("the",k,"bond")
                #print("the TVB-bond between",i+1,j+1,"dist",dist)
                new_Bond_TVB[i,0]=float(dist)
                new_Bond_TVB[i,1]=3
                new_Bond_TVB[i,2]=Atom[i,0]
                new_Bond_TVB[i,3]=Atom[j,0]
                    
np.savetxt("new_Bond_TVB.txt", new_Bond_TVB,fmt="%3d")
np.savetxt("dist_array_TVB.txt", dist_array_TVB,fmt="%3d")

elapsed_time = (time.time() - start)
print("execution time",elapsed_time)
