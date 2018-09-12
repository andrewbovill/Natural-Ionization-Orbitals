#!/usr/bin/python

#######################################################################################################################
#  
#  NIOcal v. 1.0
#  By Hassan Harb
#  University of California Merced
#  Program that reads in SCF Densities of ground and detached states and returns the eigenvalues and eigenvectors of 
#  natural ionization orbitals
#  
#  For more details, check: J. Chem. Phys., 144, 204117 (2016) 
#
#  Last edited by Hassan Harb, September 6, 2018
#
####################################################################################################################### 

from __future__ import division
import sys
import math
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal

# PRaw1 = P of ground state
# PRaw2 = P of detached state

# Need to pull out Palpha and Pbeta, C(MO)alpha and C(MO)beta
# FCHK file contains TOTAL SCF Density and Spin SCF Density
# Need to calculate Alpha SCF Density and Beta SCF Density from total & spin SCF Density


# Function: Transform packed array into NBasis x NBasis matrix
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

def sci_notation(n):
    a = '%.8E' % n
    return '%.8E' % Decimal(n)
#    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

# Read in arrays 
PRaw1 = genfromtxt('./ground_P.dat',dtype=None)
PRaw2 = genfromtxt('./detached_P.dat',dtype=None)

# Calculate dimensions of matrices
NBasis  = int((np.sqrt(8*len(PRaw1)+1)-1)/2)

print "Number of basis functions = ", NBasis
print "\n"
# print "PRaw1\n", PRaw1
# print "PRaw2\n", PRaw2

# Transform arrays into matrices, call function symmetrize

# P1 = NBasis x NBasis Matrix, density matrix of gorund state
# P2 = NBasis x NBasis Matrix, density matrix of detached state

P1 = np.zeros((NBasis,NBasis))
P2 = np.zeros((NBasis,NBasis))

P1 = symmetrize(PRaw1)
P2 = symmetrize(PRaw2)

print "Pg = \n", P1, "\n"
print "Pi = \n", P2, "\n"

# Calculate Delta P, Dp = P2 - P1
Dp = P2 - P1
print "Delta P = \n", Dp, "\n"

# Read MO Coefficients
Craw = genfromtxt('./ground_C.dat',dtype=None)
C = np.zeros((NBasis,NBasis))
n=0
for i in range(0,NBasis):
  for j in range(0,NBasis):
#    Nc=len(Craw)
    C[j,i]=Craw[n]
    n=n+1
# Calculate Overlap Matrix
# Note to self: Once S is calculated, next checkpoint should be <PS>

CInv = np.linalg.inv(C)
#Overlap Matrix S
# S = (C**(-1))t*(C**(-1))
S = np.dot(np.transpose(CInv,CInv))
print "Overlap Matrix = \n", S

#Checkpoint: <Dp.S> = -1

DpSa = np.dot(Dpa,S)
DpSb = np.dot(Dpb,S)

print "Alpha <DeltaP S> = ", np.trace(DpSa)
print "Beta <DeltaP S> = ", np.trace(DpSb)

#Calculate S**(0.5)
Svals, Svecs = np.linalg.eig(S)
Sval_minhalf = (np.diag(Svals**(0.5)))
Shalf = np.dot(Svecs,np.dot(Sval_minhalf,np.transpose(Svecs)))
print "S**(0.5) = \n", Shalf

#Calculae S**(0.5).Dpa.S**(0.5) and S**(0.5).Dpb.S**(0.5)
#Let these two matrices be Z1 and Z2, respectively

Z1 = np.dot(Shalf,np.dot(Dpa,Shalf))
Z2 = np.dot(Shalf,np.dot(Dpb,Shalf))

print "Z1 = \n", Z1
print "Z2 = \n", Z2

#Diagonaize Z1 and Z2, obtain eigenvalues (e1 and e2) and eigenvectors (U1 and U2)

e1, U1 = np.linalg.eig(Z1)
e2, U2 = np.linalg.eig(Z2)

e1 = np.sort(e1)
e2 = np.sort(e2)

print "Alpha NIO Eigenvalues = \n", e1
print "Beta NIO Eigenvalues = \n", e2
print "Alpha NIO Eigenvectors = \n", U1
print "Beta NIO Eigenvectors = \n", U2

#Backtransform U1 and U2 into the AO Basis

V1 = np.dot(np.linalg.inv(Shalf),U1)
V2 = np.dot(np.linalg.inv(Shalf),U2)

print "V1 =\n", V1, "\n"
print "V2 =\n", V2, "\n"



