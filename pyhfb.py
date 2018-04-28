# -----------------------------------------------------------------
# python code for solving HFB equation for axially deformed nuclei
# -----------------------------------------------------------------
# April, 2018 Jiangming Yao
# -----------------------------------------------------------------
# ...... import libraries
import numpy as np
from numpy import linalg as LA
from numpy import array, dot, diag, reshape, transpose
from scipy.linalg import eigvalsh
from scipy.integrate import odeint, ode
from sys import argv
import pandas as pd
import sys as sys

# import constants which can be used directly and their values are not changeable
from constant import *
# import common variables which are changeable: cm.xxx
import common as cm


# Headline
print("----------------------------------------------\n")
print("This is a HFB code for axially deformed nuclei\n")
print("----------------------------------------------\n")

class HO_Parameters:
    def __init__(self,emax,hwHO):
        self.emax=int(emax)
        self.nmax=int(emax/2)   # emax=2*n+l
        self.lmax=int(emax)
        self.hwHO=int(hwHO)
        self.tmax=1             # 0,1
        self.ljmax=int(2*emax+1)
        if(emax == 4):
        	self.nljmax=15
        elif(emax == 6):
        	self.nljmax=28
        elif(emax == 8):
        	self.nljmax=45
        else:
        	print(f" The HO basis with emax={emax} is not implemented yet\n")
        	sys.exit(0) 

class TBPC:
    def __init__(self,emax):
    	self.bmax=emax
    	self.a12max=emax



# Initiliztion of Global variables

emax=int(input("Please enter a number for emax:"))
hwHO=int(input("Please enter a number for the hwHO:"))
HO=HO_Parameters(emax,hwHO)


def Error_MSN(msn1,msn2):
	print("Error occurs in"+msn1+":"+msn2+"\n")
 
#------------------------------- 
# Read one-body matrix elements
#-------------------------------
def Read_ME1B(File1b): 
    check = 0
    me1b=pd.read_csv(File1b,delimiter=r'\s+',header=None,skiprows=1)
    me1b.columns=['t','lj','n1','n2','vme1b']  

    #---------------------------------------------
    # In Python: index in an array starts from 0
    # e.g. A[5] contains actually "a[0],...,a[4]""
    H1B=np.zeros((HO.tmax+1,HO.ljmax+1,HO.nmax+1,HO.nmax+1))  # e.g. n=0,1,..,nmax

    print(f"itmax={HO.tmax}, ljmax={HO.ljmax}, nmax={HO.nmax}\n")
    for item in me1b.index.values:
        it = me1b['t'][item]   # it=0 for n and =1 for proton
        lj = me1b['lj'][item]
        n1 = me1b['n1'][item]
        n2 = me1b['n2'][item] 

		# check if the index is outside of the model space
        check = [1 if it > 1 else 0]
        check = [1 if lj > HO.ljmax else 0]
        check = [1 if n1 > HO.nmax else 0] 
        check = [1 if n2 > HO.nmax else 0]  
        if check == 1:
            Error_MSN("Read_ME1B","Model Space are Not Proper for the Matrix Elements")
            sys.exit(0) 
            # give the value to H1B
        H1B[it,lj,n1,n2] = me1b['vme1b'][item] 
    return H1B

def Read_ME2B(File2b): 
	me2b=pd.read_csv(File2b,delimiter=r'\s+',header=None,skiprows=1)
	me2b.columns=['t1','t2','t3','t4', 'n1','n2','n3','n4','lj1','lj2','lj3','lj4','J12','vme2b']  
	H2B=np.zeros((bmax,a12max,a12max))

	for item in me1b.index.values:
		t1 = me2b['t1'][item]
		t2 = me2b['t2'][item]
		n1 = me2b['n1'][item]
		n2 = me2b['n2'][item]

		lj1 = me2b['lj1'][item]
		lj2 = me2b['lj2'][item]

		J12 = me2b['J12'][item] 

		b12 = b12_index(T12,P12,J12)
		a12 = a12_index(t1,t2,n1,n2,lj1,lj2)
		a34 = a12_index(t3,t4,n3,n4,lj3,lj4)

		H2B[b12,a12,a34] = me2b['vme2b'][item] 
	return H1B

def b12_index(T12,P12,J12):
	return T12+P12+J12




#----------------------------------------------------
# The main program starts here
#----------------------------------------------------
# initialization of model space
#  HO=HO_Parameters(6,20)
def main():
    print(f"Pi is {pi}.")
    print(f"The emax of HO is {HO.emax} \n"
  	    f"The nmax of HO is {HO.nmax} \n" 
  	    f"The lmax of HO is {HO.lmax} \n" 
  	    f"The frequency of HO is {HO.hwHO}\n"
  	    f"The number of levels in j-scheme is {HO.nljmax}\n") 

#   Read Hamiltonian
    H1B = Read_ME1B("me1b.dat")  
#   print(H1B[:,:,:,:])
#   H2B = Read_ME2B("me2b.dat")

# Diagonalize the H1B using the libarary subroutine 
# dot(H1B[:,:], v[:,i]) = w[i] * v[:,i]
    w, v = LA.eigh(H1B)
  # w: The eigenvalues, each repeated according to its multiplicity.
  # v: The normalized (unit “length”) eigenvectors, 
  #    such that the column v[:,i] is the eigenvector corresponding to the eigenvalue w[i].
    print("The eigen values:")
    print(w) 


#------------------------------------------------------------------------------
# make executable
#------------------------------------------------------------------------------
if __name__ == "__main__":
  main()
