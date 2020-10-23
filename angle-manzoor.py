import numpy as np

Atom=np.loadtxt('atom.data')
Bond=np.loadtxt('bond.data')
Angle=np.loadtxt('angle.data')
Dihedral=np.loadtxt('dihedral.data')

new_Angle=Angle.copy()
    
angle_array=[]#np.array()
for i in range((len(Angle))):
    a=Atom[int(Angle[i,2])-1,4:]
    b=Atom[int(Angle[i,3])-1,4:]
    c=Atom[int(Angle[i,4])-1,4:]

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle) 

    b=np.degrees(angle)
    #print("the angle",i+1,"is",b)
    angle_array.append(b)

    if 90<b<93 :
        new_Angle[i,1]=1
    else:
        new_Angle[i,1]=2

np.savetxt("new_Angle.txt", new_Angle,fmt="%3d")
np.savetxt("Angle_array.txt",angle_array,fmt="%3d")
