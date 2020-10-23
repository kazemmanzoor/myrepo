# Program to calculate dihedral angle 
# Equation from http://www.stem2.org/je/proteina.pdf
#
# Author: Dr. Walter Filgueira de Azevedo Jr.
# Date: April, 9th 2016.

def initial_vectors():
    """Function to set up initial vectors"""
    import numpy as np
    
    # Set initial values for arrays
    #p1 = np.zeros(3)
    #p2 = np.zeros(3)
    #p3 = np.zeros(3)
    #p4 = np.zeros(3)


    Atoms=np.loadtxt('atom.data')
    Bond=np.loadtxt('bond.data')
    Angle=np.loadtxt('angle.data')
    Dihedral=np.loadtxt('dihedral.data')
    
    # Set initial coordinates (http://www.stem2.org/je/proteina.pdf)
    p1 = Atoms[20,4:]
    p2 = Atoms[0,4:]
    p3 = Atoms[1,4:]
    p4 = Atoms[3,4:]

    print("p1=",p1)
    print("p2=",p2)
    print("p3=",p3)
    print("p4=",p4)
    



    return p1,p2,p3,p4

def calc_q_vectors(p1,p2,p3,p4):
    """Function to calculate q vectors"""
    import numpy as np
    
    # Calculate coordinates for vectors q1, q2 and q3
    q1 = np.subtract(p2,p1) # b - a
    q2 = np.subtract(p3,p2) # c - b
    q3 = np.subtract(p4,p3) # d - c

    #print("q1=",q1)
    #print("q2=",q2)
    #print("q3=",q3)
    
    return q1,q2,q3

def calc_cross_vectors(q1,q2,q3):
    """Function to calculate cross vectors"""    
    import numpy as np
    
    # Calculate cross vectors
    q1_x_q2 = np.cross(q1,q2)
    q2_x_q3 = np.cross(q2,q3)

    #print("q1_x_q2=",q1)
    #print("q2_x_q3=",q2)
    
    
    return q1_x_q2, q2_x_q3 

def calc_normals(q1_x_q2,q2_x_q3):  
    """Function to calculate normal vectors to planes"""
    import numpy as np
    
    # Calculate normal vectors
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))
    
    return n1,n2

def calc_orthogonal_unit_vectors(n2,q2):
    """Function to calculate orthogonal unit vectors"""
    import numpy as np
    
    # Calculate unit vectors
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2,q2)))
    u2 = np.cross(u3,u1)
    
    return u1,u2,u3
    
def calc_dihedral_angle(n1,u1,u2,u3):
    """Function to calculate dihedral angle"""
    import numpy as np
    import math
    
    # Calculate cosine and sine
    cos_theta = np.dot(n1,u1)
    sin_theta = np.dot(n1,u2)
    
    # Calculate theta
    theta = -math.atan2(sin_theta,cos_theta)    # it is different from atan2 from fortran math.atan2(y,x)
    theta_deg = np.degrees(theta)
    
    # Show results
    #print("theta (rad) = %8.3f"%theta)
    #print("theta (deg) = %8.3f"%theta_deg)
    return theta_deg

##################################################################################################

"""Function to set up initial vectors"""
import numpy as np
    
# Set initial values for arrays

Atom=np.loadtxt('atom.data')
Bond=np.loadtxt('bond.data')
Angle=np.loadtxt('angle.data')
Dihedral=np.loadtxt('dihedral.data')

new_Dihedral=Dihedral.copy()
    
# Set initial coordinates (http://www.stem2.org/je/proteina.pdf)

dihedral_array=[]#np.array()
for i in range((len(Dihedral))):
    p1=Atom[int(Dihedral[i,2])-1,4:]
    p2=Atom[int(Dihedral[i,3])-1,4:]
    p3=Atom[int(Dihedral[i,4])-1,4:]
    p4=Atom[int(Dihedral[i,5])-1,4:]
    


    # Call calc_q_vectors(p1,p2,p3,p4) function
    q1,q2,q3 = calc_q_vectors(p1,p2,p3,p4)
    
    # Call calc_cross_vectors(q1,q2,q3) function
    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1,q2,q3)
    
    # Call calc_normalss(q1_x_q2,q2_x_q3) function
    n1, n2 = calc_normals(q1_x_q2,q2_x_q3)
    
    # Call calc_orthogonal_unit_vectors(n2,q2) function
    u1,u2,u3 = calc_orthogonal_unit_vectors(n2,q2)
    
    # Call calc_dihedral_angle(u1,u2,u3) function
    b=calc_dihedral_angle(n1,u1,u2,u3)
    print("the dihedral",i+1,"is",b)
    dihedral_array.append(b)

    if 80<np.abs(b)<81:
        new_Dihedral[i,1]=3
    elif 27<np.abs(b)<28:
        new_Dihedral[i,1]=1
    elif 76<np.abs(b)<77:
        new_Dihedral[i,1]=2
    else:
        new_Dihedral[i,1]=4
        
np.savetxt("new_Dihedral.txt", new_Dihedral,fmt="%3d")
np.savetxt("dihedral_array.txt", dihedral_array,fmt="%3d")
    
    


