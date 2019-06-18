#!/usr/bin/python

#######################################################################################################################
#
#  NIO  v. 1.0
#  By Hassan Harb
#  Program that reads in SCF Densities of ground and detached states and returns the eigenvalues and eigenvectors of
#  natural ionization orbitals
#
#  usage: python nio.py ground.fchk detached.fchk
#
#  For more details, check: J. Chem. Phys., 144, 204117 (2016)
#
#  Last edited by Hassan Harb, June 18, 2019
#
#######################################################################################################################

from __future__ import division
import sys
import math
import cmath
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal 
from BEATLES import *

#Function: Symmetrize: Transform packed array into NBasis x NBasis 
def symmetrize(a):
  Nbas = int((np.sqrt(8*len(a)+1)-1)/2)
  b = np.zeros((Nbas,Nbas))
  n = 0
  for i in range(0,Nbas):
    for j in range(0,i+1):
      b[i,j]=a[n]
      b[j,i]=a[n]
      n=n+1
  return b

# Function: convert output to scientific notation
def sci_notation(n):
    a = '%.8E' % n
    return '%.8E' % Decimal(n.real)

# Function: Calculate Ncol, the column vector that shows the MO contributions of NIOs
def CalNCol(T,k,NBasis,S):
   Ti = np.zeros(NBasis)
   for i in range(0,NBasis):
      Ti[i] = T[i,k].real

   T_in = np.dot(Ti,Ti)
   N = np.outer(Ti,Ti)
   I = np.identity(NBasis)
   N  = np.multiply(N,I)
   N_cont = np.sum(N)
   print "Inner product =", T_in
   print "Contraction = ", N_cont
   print "trace N = ", np.trace(N)

   NCol = np.zeros(NBasis)
   for i in range(0,NBasis):
       for j in range(0,NBasis):
          NCol[i] = NCol[i] + N[j,i]
       if (NCol[i] < 0.0000000001 ):
          NCol[i] = 0
   NIOPops(NCol,NBasis)

#Function: Calculate the percentage contributions of MOs in each NIO
def NIOPops(NCol,NBasis):
    for i in range(0,NBasis):
       if (np.absolute(NCol[i]) >= 0.01):
           print "The Contribution from MO ", i+1 , " is ", np.around(NCol[i],decimals=2)*100, "percent \n"
    print "---------------\n"
#Part 1: Read in the matrix files from two checkpoint files

NBasis = 0
filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = "NIO-"+filename1
acc = 8 #accuracy to the acc's decimal place
minEigVal = 0.2

print "NIO: Calculates the density difference between two checkpoint files.\n"
print "Checkpoint 1:", filename1
print "Checkpoint 2:", filename2
 
with open(filename1,'r') as origin:
    for line in origin:
        if  "Number of basis functions" in line:
          words = line.split()
          for i in words:
              for letter in i:
                  if(letter.isdigit()):
                      NBasis = NBasis*10 + int(letter)   

print "Number of Basis Functions = ", NBasis, "\n"
MOElements = NBasis * NBasis
PElements = int(NBasis*(NBasis+1)/2)
NBasis2=0

with open(filename2,'r') as origin:
    for line in origin:
        if  "Number of basis functions" in line:
          words = line.split()
          for i in words:
              for letter in i:
                  if(letter.isdigit()):
                      NBasis2 = NBasis2*10 + int(letter)

if (NBasis != NBasis2):
    print "ERROR: Files have different numbers of basis functions"
    exit()

MOlines = int(MOElements/5) + 1
Plines = int(PElements/5) + 1

if (MOElements % 5 ==0):
   MOlines = int(MOElements/5)
if (PElements % 5 ==0):
   Plines = int(PElements/5)

MOrawa = np.zeros(MOElements)
MOrawb = np.zeros(MOElements)
TotalPraw1 = np.zeros(PElements) 
SpinPraw1 = np.zeros(PElements)

TotalPraw2 = np.zeros(PElements)
SpinPraw2 = np.zeros(PElements)

p = 0
r = 0
AOE = 0
with open(filename1,'r') as origin:
    for i, line  in enumerate(origin):
        if "Alpha Orbital Energies" in line:
              AOE = i
        if  "Alpha MO coefficients" in line:              
              i=i+1
              AMO=i
              j=i+MOlines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   MOrawa[r] = nextline[p]
                   r = r+1
                 p = 0
        if "Beta MO coefficients" in line:
              r = 0
              i=i+1
              BMO = i
              j=i+MOlines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   MOrawb[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Total SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              j=i+Plines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   TotalPraw1[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Spin SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              j=i+Plines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   SpinPraw1[r] = nextline[p]
                   r = r+1
                 p = 0

with open(filename2,'r') as origin:
    for i, line  in enumerate(origin):
        if  "Total SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              j=i+Plines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   TotalPraw2[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Spin SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              j=i+Plines-1
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   SpinPraw2[r] = nextline[p]
                   r = r+1
                 p = 0

MOrawa = np.around(MOrawa,decimals=acc)
MOrawb = np.around(MOrawb,decimals=acc)
TotalPzraw1 = np.around(TotalPraw1,decimals=acc)
SpinPraw1 = np.around(SpinPraw1,decimals=acc)
TotalPraw2 = np.around(TotalPraw2,decimals=acc)
SpinPraw2 = np.around(SpinPraw2,decimals=acc)

PalphaRaw1 = (np.add(TotalPraw1,SpinPraw1)) * 0.5
PbetaRaw1 = (np.subtract(TotalPraw1,SpinPraw1)) * 0.5 

PalphaRaw2 = (np.add(TotalPraw2,SpinPraw2)) * 0.5
PbetaRaw2 = (np.subtract(TotalPraw2,SpinPraw2)) * 0.5

#print "Ground state: Alpha Density matrix =\n", PalphaRaw1
#print "Ground state: Beta Density matrix =\n", PbetaRaw1

#print "Detached state: Alpha Density matrix =\n", PalphaRaw2
#print "Detached state: Beta Density matrix =\n", PbetaRaw2

#Part 2: The NIO Code

P1a = np.zeros((NBasis,NBasis))
P2a = np.zeros((NBasis,NBasis))
P1b = np.zeros((NBasis,NBasis))
P2b = np.zeros((NBasis,NBasis))

#Density Matrices

P1a = symmetrize(PalphaRaw1)
P2a = symmetrize(PalphaRaw2)
P1b = symmetrize(PbetaRaw1)
P2b = symmetrize(PbetaRaw2)

#print "Alpha Density 1 = \n", P1a
#print "Beta Density 1 = \n", P1b
#print "Alpha Density 2 = \n", P2a
#print "Beta Density 2 = \n", P2b

Pa1, Pb1 = MatGrab(filename1,NBasis,2)
Pa2, Pb2 = MatGrab(filename2,NBasis,2)

#Difference Density Matrices
Dpa = Pa2 - Pa1
Dpb = Pb2 - Pb1 

print "Alpha Difference Density = \n", Dpa
print "Beta Difference Density = \n", Dpb

#MO Coefficient Matrix
C = np.zeros((NBasis,NBasis))
Cb = np.zeros((NBasis,NBasis))
t=0
for i in range(0,NBasis):
   for j in range(0,NBasis):
     C[j,i]=MOrawa[t]
     Cb[j,i]=MOrawb[t]
     t=t+1

#C = np.transpose(C)
#Cb = np.transpose(Cb)

#Inverse of C:
CInv = np.linalg.inv(C)
CbInv = np.linalg.inv(Cb)
#Overlap Matrix:
# S = (C**(-1)t*(C**(-1))
S_old = np.dot(np.transpose(CInv),CInv)
Sb = np.dot(np.transpose(CbInv),CbInv)
#Checkpoint: <Dp.S> = -1

C_test = MatGrab(filename1,NBasis,1)
Cb_test = MatGrab(filename1,NBasis,-1)
S = GetOverlap(C_test,NBasis)
C_test = column2square(C_test,NBasis)
Cb_test = column2square(Cb_test,NBasis)

DpSa = np.dot(Dpa,S)
DpSb = np.dot(Dpb,S)

print "Alpha <DeltaP S> = ", np.trace(DpSa)
print "Beta <DeltaP S> = ", np.trace(DpSb)

C = C_test
Cb = Cb_test

#Calculate S**(0.5)
Svals, Svecs = np.linalg.eig(S)
Sval_minhalf = (np.diag(Svals**(0.5)))
Shalf = np.dot(Svecs,np.dot(Sval_minhalf,np.transpose(Svecs)))
Shalf = Shalf.real

S_try = np.dot(Shalf,Shalf)
diff = S_try - S

#Calculae S**(0.5).Dpa.S**(0.5) and S**(0.5).Dpb.S**(0.5)
#Let these two matrices be Z1 and Z2, respectively

Sleft, Sval, Sright = np.linalg.svd(S)

I = np.eye(NBasis)

Sval_I = np.multiply(Sval,I)

for i in range(0,NBasis):
    Sval_I[i,i] = Sval_I[i,i]**(0.5)


S_test = np.dot(Sleft,np.dot(Sval_I,Sright))

Shalf = S_test

Z1 = np.dot(Shalf,np.dot(Dpa,Shalf))
Z2 = np.dot(Shalf,np.dot(Dpb,Shalf))


#Diagonaize Z1 and Z2, obtain eigenvalues (e1 and e2) and eigenvectors (U1 and U2)

e1, U1 = np.linalg.eig(np.dot(Shalf,np.dot(Dpa,Shalf)))
e2, U2 = np.linalg.eigh(np.dot(Shalf,np.dot(Dpb,Shalf)))

e1 = e1.real
e2 = e2.real
U1 = U1.real
U2 = U2.real

print "Alpha NIO Eigenvectors = \n", U1
print "Beta NIO Eigenvectors = \n", U2
print "Alpha NIO Eigenvalues = \n", e1
print "Beta NIO Eigenvalues = \n", e2

#Backtransform U1 and U2 into the AO Basis

V1 = np.dot(np.linalg.inv(Shalf),U1)
V2 = np.dot(np.linalg.inv(Shalf),U2)

V1 = V1.real
V2 = V2.real

print "Alpha NIO Coefficients =\n", V1, "\n"
print "Beta NIO Coefficients =\n", V2, "\n" 

########## HH +

# Screen the NIO eigenvalues, pick up the ones we need.


# Calculate T=CTranspose.S.U

Ta = np.dot(np.transpose(C),np.dot(S,V1))
Tb = np.dot(np.transpose(Cb),np.dot(S,V2))


for i in range(0,NBasis):
    if (np.absolute(e1[i].real) >= minEigVal):
        print "The ", i,"th element of e1 has an eigenvalue of ", np.around(e1[i].real,decimals=3) ," perform population analysis\n"
        print "Population analysis on Alpha NIO ", i+1, "\n"
        print "----------------------------\n"
        CalNCol(Ta,i,NBasis,S)        

    if (np.absolute(e2[i].real) >= minEigVal):
        print "The ", i,"th element of e2 has an eigenvalue of ", np.around(e2[i].real,decimals=3) ," perform population analysis\n"
        print "Population analysis on Beta NIO ", i+1, "\n"
        print "----------------------------\n"
        CalNCol(Tb,i,NBasis,S)

# Here we need to define a function that takes in C, S, and either of V1 or V2 (Depending on whether e1 or e2 has a non zero value)
# We then extract the ith row of T, named Ti, and calculate N which is the outer product of Ti
# After that we calculate NCol which is a column vector that sums over the rows of N 

# Variables needed for the function: 
#    1. Matrix T (NBasis,NBasis): inputted as either Ta or Tb (T = C(t)*S*V)
#    2. Scalar i: indicates the position of the desired NIO
#    3. Number of Basis Functions: NBasis

# Generate a Matrix (Call it N) from the outer product of T(i) columns
#Ti = np.zeros(NBasis)
#for i in range(0,NBasis):
#   Ti[i] = Tb[1,i]
#N = np.outer(Ti,Ti)
#NCol = np.zeros(NBasis)
#for i in range(0,NBasis):
#    for j in range(0,NBasis):
#       NCol[i] = NCol[i] + N[i,j]
#print "T1 =", T1, "\n"
#print "Ti =", Ti, "\n"
#print "N = ", N, "\n"
#print "NCol = ", NCol, "\n"
# Sum over the rows of the newly-formed Matrix N to get a column vector 

#print "Ut.U = ", np.dot(np.transpose(U1),U1)
#print "Ut.U = ", np.dot(np.transpose(U2),U2)
#print "Ct.S.C = ", np.dot(np.transpose(C),np.dot(S,C))
#print "Snippet test: complete\n"
#print "Tb.Tbt = ", np.dot(Tb,np.transpose(Tb))
#print "Ta.Tat = ", np.dot(Ta,np.transpose(Ta))


##########  HH -

# Part 3: write everything to a new fchk file:
# We need to replace Alpha Orbital Energies with e1 & Beta orbital energies with e2
# We also need to replace Alpha MO Coefficients with V1 & Beta MO Coefficients with V2

# print "U.Ut = ", np.dot(U1,np.transpose(U1))

pointer=0
counter=1
with open(filename1,'r') as origin:
  data = origin.readlines()
  if "Alpha Orbital Energies" in line:
    AOE = i  
    BOE = AOE + int(NBasis/5) + 1
  with open(filename3,'w') as f2:

      print "Writing results to new output file: ", filename3, " ... "

#Read in data from old fchk (unchanged)
      while (pointer < AOE+1):
         f2.write(data[pointer])
         pointer = pointer+1
#Write the Alpha eigenvalues, overwrite the ALpha Orbital energies
      for j in range(0,NBasis):
          f2.write(" ")
          if (e1[j] >= 0):
             f2.write(" ")
          f2.write(str(sci_notation(e1[j].real)))
          if (counter%5 == 0):
              f2.write("\n")
              counter=0
          counter=counter+1
      counter =1      
#Write the Beta eigenvalues, overwrite the Beta Orbital energies
      BOE = AOE + (int(NBasis/5)+2)
 #     print "BOE at line =", BOE
      if (NBasis%5 != 0):
          f2.write("\n")
      if (NBasis%5 == 0):
          BOE = BOE - 1 
      f2.write(data[BOE])
      for j in range(0,NBasis):
          f2.write(" ")
          if (e2[j] >= 0):
             f2.write(" ")
          f2.write(str(sci_notation(e2[j].real)))
          if (counter%5 ==0):
              f2.write("\n")
              counter=0
          counter = counter+1
      counter =1
#      f2.write("\n")
#Write the Alpha NIO eigenvectors, overwrite the Alpha MO Coefficients
      AMO = BOE + (int(NBasis/5)+2)
      if (NBasis%5 != 0):
          f2.write("\n")
      if (NBasis%5 == 0):
          AMO = AMO - 1
      f2.write(data[AMO])
      for i in range(0,NBasis):
          for j in range(0,NBasis):
               f2.write(" ")
               if (V1[j,i] >= 0):
                  f2.write(" ")
               f2.write(str(sci_notation(V1[j,i].real)))
               if (counter%5 ==0):
                   f2.write("\n")
                   counter=0
               counter = counter + 1
         # counter =1
      counter = 1
#      f2.write("\n")
#Write the Beta NIO eigenvectors, overwrite the Beta MO coefficients
      BMO = AMO + (int(NBasis*NBasis/5))+2
      if (NBasis%5 != 0):
          f2.write("\n")
      if (NBasis%5 == 0):
          BMO = BMO - 1
      f2.write(data[BMO])
      for i in range(0,NBasis):
             for j in range(0,NBasis):
                  f2.write(" ")
                  if (V2[j,i] >= 0):
                     f2.write(" ")
                  f2.write(str(sci_notation(V2[j,i].real)))
                  if (counter%5 ==0):
                      f2.write("\n")
                      counter=0
                  counter = counter + 1
           #  counter = 1
      counter = 1
      if (NBasis%5 != 0):
         f2.write("\n")
#Copy the rest of the checkpoint file (unchanged)      
      pointer = BMO + (int(NBasis*NBasis/5))+2
      while (pointer < len(data)):
         f2.write(data[pointer])
         pointer = pointer+1
print "NBasis = ", NBasis
print "Copying data to new checkpoint file: DONE"
print "NIO Calculation done!"

