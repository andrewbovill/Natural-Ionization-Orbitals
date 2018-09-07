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
    for j in range(i,Nbas):
      b[i,j]=a[n]
      b[j,i]=a[n]
      n=n+1
  return b

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


