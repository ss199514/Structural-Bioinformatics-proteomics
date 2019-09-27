import math
import re
#import pandas as pd
from math import pow
import numpy as np
from Bio.PDB import *
from numpy import cross, eye, dot
from scipy.linalg import expm, norm

fh2 = open("Q5.txt", "w")
parser = PDBParser()
structure = parser.get_structure('PHA-L', '2gb1.pdb')
chi1 = ["CYS", "SER", "THR", "VAL"]
chi2 = ["ASN", "ASP", "HIS", "ILE", "LEU", "PHE", "PRO", "TRP", "TYR"]
chi3 = ["GLN", "GLU", "MET"]
chi4 = ["LYS"]
chi5 = ["ARG"]
residueList = []
residueNumber = []
atomList = []
vectorList = []
resnum = 1
for residue in structure.get_residues():
	if residue.get_resname() in chi1:
		#print(residue.get_resname(), resnum)
		n = 0
		for atom in residue:
			if n < 6:
				#print(n)
				if(atom.get_name() != "O" and atom.get_name() != "C"):
					#print(atom.get_name(), atom.get_vector())
					residueList.append(residue.get_resname())
					residueNumber.append(resnum)
					atomList.append(atom.get_name())
					vectorList.append(atom.get_vector())
			else:
				next
			n += 1
	elif residue.get_resname() in chi2:
		#print(residue.get_resname(), resnum)
		n = 0
		for atom in residue:
			if n < 7:
				#print(n)
				if(atom.get_name() != "O" and atom.get_name() != "C"):
					#print(atom.get_name(), atom.get_vector())
					residueList.append(residue.get_resname())
					residueNumber.append(resnum)
					atomList.append(atom.get_name())
					vectorList.append(atom.get_vector())
			else:
				next
			n += 1
	elif residue.get_resname() in chi3:
		#print(residue.get_resname(), resnum)
		n = 0
		for atom in residue:
			if n < 8:
				#print(n)
				if(atom.get_name() != "O" and atom.get_name() != "C"):
					#print(atom.get_name(), atom.get_vector())
					residueList.append(residue.get_resname())
					residueNumber.append(resnum)
					atomList.append(atom.get_name())
					vectorList.append(atom.get_vector())
			else:
				next
			n += 1
	elif residue.get_resname() in chi4:
		#print(residue.get_resname(), resnum)
		n = 0
		for atom in residue:
			if n < 9:
				#print(n)
				if(atom.get_name() != "O" and atom.get_name() != "C"):
					#print(atom.get_name(), atom.get_vector())
					residueList.append(residue.get_resname())
					residueNumber.append(resnum)
					atomList.append(atom.get_name())
					vectorList.append(atom.get_vector())
			else:
				next
			n += 1
	elif residue.get_resname() in chi5:
		#print(residue.get_resname(), resnum)
		n = 0
		for atom in residue:
			if n < 10:
				#print(n)
				if(atom.get_name() != "O" and atom.get_name() != "C"):
					#print(atom.get_name(), atom.get_vector())
					residueList.append(residue.get_resname())
					residueNumber.append(resnum)
					atomList.append(atom.get_name())
					vectorList.append(atom.get_vector())
			else:
				next
			n += 1	
	resnum += 1
fh2.write("Residue\tBackboneTorsion\tAtomsList\tValue\n")
chinum = 1
for i in range(0,len(residueList)-3):
	if(residueList[i] == residueList[i+1] and residueList[i] == residueList[i+2] and residueList[i] == residueList[i+3]):
		#print(residueList[i], residueNumber[i], atomList[i])
		#print(residueList[i+1], residueNumber[i+1], atomList[i+1])
		#print(residueList[i+2], residueNumber[i+2], atomList[i+2])
		#print(residueList[i+3], residueNumber[i+3], atomList[i+3])
		#print("chi"+str(chinum))
		dih = 57.2958*calc_dihedral(vectorList[i], vectorList[i+1], vectorList[i+2], vectorList[i+3])
		fh2.write(str(residueList[i])+str(residueNumber[i])+"\tchi"+str(chinum)+"\t")
		fh2.write(str(atomList[i])+"-"+str(atomList[i+1])+"-"+str(atomList[i+2])+"-"+str(atomList[i+3]))
		fh2.write(":\t"+str(dih) +"\n")
		#print(str(57.2958*calc_dihedral(vectorList[i], vectorList[i+1], vectorList[i+2], vectorList[i+3])))
		chinum += 1
	else:
		chinum = 1
		next

print("List of all backbone torsions is in the output file Q5.txt")
	
fh2.close()
