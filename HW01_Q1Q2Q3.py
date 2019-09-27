import math
import re
from math import pow
import numpy as np
from Bio.PDB import *
from numpy import cross, eye, dot
from scipy.linalg import expm, norm

print("Please check output files Q1.txt, Q2.txt and Q3.txt")
#1. Read the pdb file and retrieve the coordinates of the backbone atoms (specifically those of Ni, CAi, and Ci.)
fh1 = open("2gb1.pdb", "r")
fh2 = open("Q1.txt", "w")
lstN = []
lstCA = []
lstC = []
i=0
for line in fh1:
	if (('\t'.join(line.split())).split("\t")[0] == "ATOM"):
		if (('\t'.join(line.split())).split("\t")[2] == "N"):
			lstN.append((('\t'.join(line.split())).split("\t")[6:9]))
			fh2.write("N"+str(i+1)+str(lstN[i])+"\n")
		if (('\t'.join(line.split())).split("\t")[2] == "CA"):
			lstCA.append((('\t'.join(line.split())).split("\t")[6:9]))
			fh2.write("CA"+str(i+1)+str(lstCA[i])+"\n")
		if (('\t'.join(line.split())).split("\t")[2] == "C"):
			lstC.append((('\t'.join(line.split())).split("\t")[6:9]))
			fh2.write("C"+str(i+1)+str(lstC[i])+"\n")
			i += 1

fh2.close()
fh1.close()

#2. Calculate and output the mean and std (standard deviation) for each of the following variables
fh2 = open("Q2.txt", "w")

#2a.bond lengths: Ni-CAi, CAi-Ci, Ci-Ni+1
#2c.distance between adjacent CAs, i.e., between CAi and CAi+1
distNCA = []
distCAC = []
distCN2 = []
distCACA = []
for i in range(0,len(lstN)):
	distNCA.append(math.sqrt((float(lstN[i][0]) - float(lstCA[i][0]))**2 + (float(lstN[i][1]) - float(lstCA[i][1]))**2 + (float(lstN[i][2]) - float(lstCA[i][2]))**2))
	#print(str(i)+"\t"+str(distNCA[i]))
	
	distCAC.append(math.sqrt((float(lstC[i][0]) - float(lstCA[i][0]))**2 + (float(lstC[i][1]) - float(lstCA[i][1]))**2 + (float(lstC[i][2]) - float(lstCA[i][2]))**2))
	#print(str(i)+"\t"+str(distCAC[i]))

for i in range(0,len(lstN)-1):
	distCN2.append(math.sqrt((float(lstN[i+1][0]) - float(lstC[i][0]))**2 + (float(lstN[i+1][1]) - float(lstC[i][1]))**2 + (float(lstN[i+1][2]) - float(lstC[i][2]))**2))
	#print(str(i)+"\t"+str(distCN2[i]))
	
	distCACA.append(math.sqrt((float(lstCA[i][0]) - float(lstCA[i+1][0]))**2 + (float(lstCA[i][1]) - float(lstCA[i+1][1]))**2 + (float(lstCA[i][2]) - float(lstCA[i+1][2]))**2))
	#print(str(i)+"\t"+str(distCACA[i]))

#2a.bond lengths: Ni-CAi, CAi-Ci, Ci-Ni+1
#2c.distance between adjacent CAs, i.e., between CAi and CAi+1

fh2.write("Mean Ni-CAi bond length Angstrom: "+str(np.mean(distNCA))+"\n")
fh2.write("Standard deviation Ni-CAi bond length Angstrom: "+str(np.std(distNCA))+"\n")
fh2.write("Mean CAi-C bond length Angstrom: "+str(np.mean(distCAC))+"\n")
fh2.write("Standard deviation CAi-C bond length Angstrom: "+str(np.std(distCAC))+"\n")
fh2.write("Mean Ci-Ni+1 bond length Angstrom: "+str(np.mean(distCN2))+"\n")
fh2.write("Standard deviation Ci-Ni+1 bond length Angstrom: "+str(np.std(distCN2))+"\n")
fh2.write("Mean distance between CAi-CAi+1 Angstrom: "+str(np.mean(distCACA))+"\n")
fh2.write("Standard deviation CAi-CAi+1 Angstrom: "+str(np.std(distCACA))+"\n")

#2b.bond angles: Ni-CAi-Ci, CAi-Ci-Ni+1, Ci-Ni+1-CAi+1
parser = PDBParser()
structure = parser.get_structure('PHA-L', '2gb1.pdb')
coordN = []
coordCA = []
coordC = []
i = 0
for atom in structure.get_atoms():
    if(atom.get_name() == "N"):
    	#print(atom.get_vector())
    	coordN.append(atom.get_vector())
    elif(atom.get_name() == "CA"):
    	#print(atom.get_vector())
    	coordCA.append(atom.get_vector())
    elif(atom.get_name() == "C"):
    	#print(atom.get_vector())
    	coordC.append(atom.get_vector())

angNCAC = []
angCACN2 = []
angCN2CA2 = []

for i in range(0,len(coordN)):
	angNCAC.append(57.2958*calc_angle(coordN[i], coordCA[i], coordC[i]))

for i in range(0,len(coordN)-1):
	angCACN2.append(57.2958*calc_angle(coordCA[i], coordC[i], coordN[i+1]))
	angCN2CA2.append(57.2958*calc_angle(coordC[i], coordN[i+1], coordCA[i+1]))

fh2.write("\n\nMean Ni-CAi-Ci bond angle degrees: "+str(np.mean(angNCAC))+"\n")
fh2.write("Standard deviation Ni-CAi-Ci bond angle degrees: "+str(np.std(angNCAC))+"\n")
fh2.write("Mean CAi-Ci-Ni+1 bond angle degrees: "+str(np.mean(angCACN2))+"\n")
fh2.write("Standard deviation CAi-Ci-Ni+1 bond angle degrees: "+str(np.std(angCACN2))+"\n")
fh2.write("Mean Ci-Ni+1-CAi+1 bond angle degrees: "+str(np.mean(angCN2CA2))+"\n")
fh2.write("Standard deviation Ci-Ni+1-CAi+1 bond angle degrees: "+str(np.std(angCN2CA2))+"\n")

fh2.close()

fh2 = open("Q3.txt", "w")
#3. For residue 30, which is a Phenylalanine, compute and output its torsional angles: phi, psi, and omega
vector1 = coordC[28]
vector2 = coordN[29]
vector3 = coordCA[29]
vector4 = coordC[29]
fh2.write("Phi: "+str(57.2958*calc_dihedral(vector1, vector2, vector3, vector4)))

vector1 = coordN[29]
vector2 = coordCA[29]
vector3 = coordC[29]
vector4 = coordN[30]
fh2.write("\nPsi: "+str(57.2958*calc_dihedral(vector1, vector2, vector3, vector4)))

vector1 = coordCA[29]
vector2 = coordC[29]
vector3 = coordN[30]
vector4 = coordCA[30]
fh2.write("\nOmega: "+str(57.2958*calc_dihedral(vector1, vector2, vector3, vector4)))
fh2.close()

'''
#4. 
#print(coordC[28])
#print(coordN[29])
b1 = coordCA[29]-coordN[29]
b2 = coordC[29]-coordCA[29]

b1vec = []
b1vec.append(b1[0])
b1vec.append(b1[1])
b1vec.append(b1[2])

b2vec = []
b2vec.append(b2[0])
b2vec.append(b2[1])
b2vec.append(b2[2])

#print(np.cross(b1vec, b2vec))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

v = b2vec
axis = b1vec
#v = b1vec
#axis = b2vec
theta = 64.62665048440938/57.2958

newCVec = np.dot(rotation_matrix(axis, theta), v)
arrCA = []
arrCA.append(coordCA[29][0])
arrCA.append(coordCA[29][1])
arrCA.append(coordCA[29][2])

newC = newCVec + arrCA

#b2dash = [1.769,4.919,0.973]-arrCA
b2dash = []
b2dash.append(1.769-arrCA[0])
b2dash.append(4.919-arrCA[1])
b2dash.append(0.973-arrCA[2])

b2dashvec = []
b2dashvec.append(b2dash[0])
b2dashvec.append(b2dash[1])
b2dashvec.append(b2dash[2])
v = b2dashvec
newOVec = np.dot(rotation_matrix(axis, theta), v)
newO = newOVec + arrCA

b3vec = []
b3vec.append(coordN[30][0] - newC[0])
b3vec.append(coordN[30][1] - newC[1])
b3vec.append(coordN[30][2] - newC[2])

b2vec = []
b2vec.append(newC[0] - coordCA[29][0])
b2vec.append(newC[1] - coordCA[29][1])
b2vec.append(newC[2] - coordCA[29][2])

v = b3vec
axis = b2vec
theta = 44.65345268888827/57.2958
newNVec = np.dot(rotation_matrix(axis, theta), v)
newN = newNVec + newC

print("New co-ordinates of C after setting Phi to 0: ",newC)
print("New co-ordinates of O after setting Phi to 0: ",newO)
'''