#IN THE NAME OF ALLAH
#new algorithm
#read line 

import numpy as np

def linalg_norm(a,b):
    #print("from function")
    return np.linalg.norm(a - b)



Atom=np.loadtxt('atom.data')
Bond=np.loadtxt('bond.data')
Angle=np.loadtxt('angle.data')
dist_array=[]

new_Bond=Bond.copy()
for i in range(len(Bond)):
    #print("Bond",i,"is",Bond[i])
    #print("##############################################################################")
    a=Bond[i,2]
    b=Bond[i,3]
    #print(a)
    #print(b)
    #print(Atom[int(a)-1])
    #print(Atom[int(b)-1])
    dist=linalg_norm(Atom[int(a)-1,4:],Atom[int(b)-1,4:])
    print("the i bond is",i,"dist",dist)
    dist_array.append(float(dist))
    
    if dist>2.24:
        new_Bond[i,1]=2
np.savetxt("new_Bond.txt", new_Bond,fmt="%3d")
np.savetxt("dist_array.txt", dist_array,fmt="%3d")


