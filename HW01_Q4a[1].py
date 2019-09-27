import math
import re
from math import pow
import numpy as np
from Bio.PDB import *
from numpy import cross, eye, dot
from scipy.linalg import expm, norm

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

parser = PDBParser()
structure = parser.get_structure('PHA-L', '2gb1.pdb')

#Setting phi to 0
n30 = (structure[0]['A'][30]['N']).get_coord()
ca30 = (structure[0]['A'][30]['CA']).get_coord()
axis = ca30-n30
phitheta = 64.62665048440938/57.2958

vectorList = []
atomNameList = []
c = 1
for res in structure.get_residues():
    if(c >= 30):
    	for atom in res.get_atoms():
    		atomNameList.append(atom.get_name())
    		vectorList.append(atom.get_coord())
    c += 1

phi0Coord = []
phi0Coord.append(n30)
phi0Coord.append(ca30)
for i in range(2,len(vectorList)):
	v = vectorList[i] - vectorList[1]
	newVec = np.dot(rotation_matrix(axis, phitheta), v)
	newCoord = newVec + ca30
	phi0Coord.append(newCoord)
	#print("New co-ordinates of "+atomNameList[i]+" after setting Phi to 0: ",newCoord)

#Setting psi to 0
axis = phi0Coord[2]-ca30
psitheta = 44.65345268888827/57.2958
psi0Coord = [] 
#0 - newC
#1 - newO
nextNForPsi = atomNameList[1:].index("N")
#print(phi0Coord[nextNForPsi-1])
#print(atomNameList[nextNForPsi+1])
for i in range(nextNForPsi+1,len(phi0Coord)):
	v = phi0Coord[i] - phi0Coord[2]
	newVec = np.dot(rotation_matrix(axis, psitheta), v)
	newCoord = newVec + phi0Coord[2]
	psi0Coord.append(newCoord)
	#print("New co-ordinates of "+atomNameList[i+2]+" after setting Psi to 0: ",newCoord)

fh1 = open("2gb1.pdb", "r")
fh2 = open("2GB1-new.pdb","w")
#round(2.675, 2)
i = 0
allNewCoord = []
for i in range(0,nextNForPsi+1):
	allNewCoord.append(phi0Coord[i])
	#print(phi0Coord[i])
for i in range(0,len(psi0Coord)):
	allNewCoord.append(psi0Coord[i])
	#print(psi0Coord[i])
i = 0
for line in fh1:
	if(('\t'.join(line.split())).split("\t")[0] == "ATOM"):
		if(int(('\t'.join(line.split())).split("\t")[5]) > 29):
			fh2.write(('\t'.join(line.split())).split("\t")[0].ljust(6)+
			('\t'.join(line.split())).split("\t")[1].rjust(5)+
			('\t'.join(line.split())).split("\t")[2].center(6)+
			('\t'.join(line.split())).split("\t")[3].ljust(4)+
			('\t'.join(line.split())).split("\t")[4].rjust(1)+
			('\t'.join(line.split())).split("\t")[5].rjust(4)+
			str(round(allNewCoord[i][0],3)).rjust(12)+
			str(round(allNewCoord[i][1],3)).rjust(8)+
			str(round(allNewCoord[i][2],3)).rjust(8)+
			('\t'.join(line.split())).split("\t")[9].rjust(6)+
			('\t'.join(line.split())).split("\t")[10].rjust(6)+
			('\t'.join(line.split())).split("\t")[11].rjust(12)+"  \n")
			i += 1
		else:
			fh2.write(line)
	else:
		fh2.write(line)	
fh2.close()
fh1.close()

print("After setting phi and psi angles to zero, new co-ordinates are computed and saved into the required format as 2GB1-new.pdb")