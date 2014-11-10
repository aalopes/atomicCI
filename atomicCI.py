#!/usr/bin/python

"""
    atomicCI.py
    
    Program for calculating the one-particle reduced density matrix (1-RDM) of
    small atomic systems. It makes use of the PyQuante library by Rick Muller
    (which you should download separately).
    To execute it simply run the associated Bash script, or execute the program
    with two input parameters: "atomic number" and "name of basis set".
    e.g. ./atomicCI.py 2 6-31g**.
      
----------------------------------------------------------------------------
This file is part of atomicCI.

Copyright (c) 2012, Alexandre Lopes
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the <organization> nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    
"""

from PyQuante.Ints      import getbasis
from PyQuante.Molecule  import Molecule
from PyQuante.CGBF      import coulomb
from numpy              import *
from functions	import *
import math
import itertools
import sys

# Read command line arguments
# First arg is atomic number, second basis name
K = int(sys.argv[1])
basisName = sys.argv[2]

molec = Molecule('mol',[(K,(0,0,0))])

# Create a basis for this molecule
# Access the CGBF atributes of the CGBF number n by doing 
# basis.bfs[CGBF_n].atribute (or, since it has some operator overloading defined, 
# basis[n] = basis.bfs[n])
# E.g. basis.bfs[0].powers (basis[0].powers)
# One can also access the contraction coefficients and the exponents of the 
# PGBF making up the CGBF by doing
# basis[CGBF_n].pcoefs[PGBF_i]
# basis[CGBF_n].pexps[PGBF_i]


# Use "6-31g**" for "best" results
# Use "STO-3G" for debbuging
bas = getbasis(molec,basis_data="6-31g**")


# Grab length of basis and number of electrons and print - some of this 
# obviously will not make sense if there's more than one type of atom
lenBasis        = 2*len(bas)
nEle            = molec.get_nel()
nucNum          = nEle # Nuclear number - one of the reasons this only works 
                       # for atoms. 
                       # If we want molecules, we must be able to change this 
                       # according to the molecule
lenNPartBasis   = math.factorial(lenBasis)/(math.factorial(nEle) * \
                  math.factorial(lenBasis-nEle) )
print "\n"
print "Using a single particle spatial basis consisting of", lenBasis/2, \
      "contracted gaussians."
print "Using a single particle basis length of", lenBasis, "."
print "The N particle basis has a length", lenNPartBasis, "."

# Let us now construct the spin-orbital basis. It will be a list consisting of 
# tuples. Each element of the list is a spin-orbital and each tuple contains two 
# values. A spin number (+1 or -1) and a CGBF object: 
# [(spin_i,CGBF_i), ..., (...)]

soBasis = []
for bfc in bas.bfs:
    soBasis.append((1,bfc))
    soBasis.append((-1,bfc))
    

# Uncomment or comment this block as necessary
# List spin and the exponents, coef, powers and pos (which indicates to which 
# atom they belong) making up the orbital basis functions
i = 0
for basis in soBasis:
    print "\n"
    print "Spin-orbital basis function number", i
    print "Spin"                              , soBasis[i][0]
    print "Powers"                            , soBasis[i][1].powers
    print "Contraction coeff."                , soBasis[i][1].pcoefs
    print "Primitive exponents"               , soBasis[i][1].pexps
    print "Center"                            , soBasis[i][1].origin
    print "___________________________________ "
    i = i+1

# Gram matrix
S=[[i.overlap(j) for i in bas.bfs] for j in bas.bfs]

# Choleski matrix (for orthonormalizing the basis)
R=matrix(linalg.cholesky(S)).T.I

# now can recover the Hamiltonian w.r.t. an ONB by computing
# R.T*h*R
# where h is the Hamiltonian expressed in the skew bas

# Let us first calculate h1 and h2 (acting only on a single particle and a 
# double particle space).
# We then use these results when working with the N particle basis.

h1 = matrix([[bra.kinetic(ket)+K*bra.nuclear(ket,(0,0,0)) for bra in bas.bfs] 
            for ket in bas.bfs])
h1Ortho = R.T*h1*R
h1Ortho = array(h1Ortho)

# h2Ortho 4-index tensor
# note that one must have [ij|kl] = [ji|kl] = [kl|ij] = ...
h2 = array([[[[coulomb(bra1,ket1,bra2,ket2) for ket2 in bas.bfs] for bra2 in 
             bas.bfs] for ket1 in bas.bfs] for bra1 in bas.bfs])
                
# Orthonormalizing h2
lenBasisR = lenBasis/2
H2m = kron(R.T,R.T)*matrix(reshape(h2,(lenBasisR**2,lenBasisR**2)))*kron(R,R)
H2Ortho = reshape(array(H2m),(lenBasisR,lenBasisR,lenBasisR,lenBasisR))
# nPartBasis = [[chi_i,chi_j, ...],[...],...,[...]] where each [...] is a 
# representation of an N particle slater determinant (|chi_i, chi_j, ... >) 
# and chi_i = i is a number corresponding to the 1 particle basis function in 
# soBasis - so this is just a list of "slater determinants". One must refer to 
# "soBasis" to get the CGBF objects

nPartBasis = []
i = 0
for basis in itertools.combinations(range(lenBasis),nEle):
    nPartBasis.append(list(basis))
# Note that len(nPartBasis) has to be = to lenNPartBasis!


# Now we need to get the Hamiltonian matrix in the orthonormal Slater basis. 

# Let H1 be the single particle part of the Hamiltonian and H2 the two particle 
# part of the Hamiltonian so that H = H1 + H2
# braNpart and ketNpart will be lists of the form [SO_i, SO_j, ... SO_k] 
# (representing a single slater determinant |SO_i, SO_j ,...> ), where SO_n 
# is a non negative integer integer indicating the basis spin-orbital.
# Note that we will be changing the basis to a orthonormal one.

H1 = zeros( (lenNPartBasis,lenNPartBasis) )
H2 = zeros( (lenNPartBasis,lenNPartBasis) )

# We shall fill only the lower triangular part of H since it is Hermitian
# i.e. j<=i

for i in range(lenNPartBasis):
    for j in range(i+1):
        #print "Building the Hamiltonian matrix. Calculating element ", i+1,j+1
        braNpart = nPartBasis[i]
        ketNpart = nPartBasis[j]
        H1[i,j] = matrixElementH1(braNpart,ketNpart,soBasis,nucNum,h1Ortho)
        H2[i,j] = matrixElementH2(braNpart,ketNpart,soBasis,H2Ortho)
        #print braNpart, ketNpart

print "\n Single particle Hamiltonian matrix \n", H1,
print "\n \n Electron interaction Hamiltonian matrix \n", H2
H = H1 + H2
print "\n Hamiltonian matrix \n", H

# Now diagonalize the matrix (since it's symmetric we can use that to reduce 
# the complexity of the problem - use eigh instead of eig)

energy, vect = linalg.eigh(H)

print "\n Energy \n", energy

# We now proceed to calculate the 1PDM. Note that above the GS vector was 
# calculated in a Slater type basis, where we orthonormalized the single 
# particle basis functions. Therefore, this Slater basis is an orthonormal one.
# We will do it like \rho = <GS| a_i^\dagger a_j |GS>, i,j = 1, ..., lenBasis.
# So, let us operate with a_i^\dagger a_j in every basis element (this will give 
# us a basis element - up to a sign of which we will keep track afterwards). 
# We then decompose  a_i^\dagger a_j|GS> into basis vectors (up to a reordering 
# of the creation operations/SO in the Slater determinant which may give rise to 
# a sign). We then compute the inner product of the Bra and the Ket, by using 
# their decomposition into the Slater basis.

# Calculating the result of operating with the creation and annihilation 
# operators on the Slater basis
# For annihilation:
# Find the SO we wish to destroy, get its position, p, (where p = 0, ..., nEle) 
# and the sign should be (-1)^p, p = 0, ..., nEle-1
# For creation:
# Let us say we want to create SO n. Then find the SO of smallest index i, s.t.
# i>n. 
# We should then create the SO n in position i and the sign should be (-1)^(i).

# Save the sign as well psi[basis,i,j]=[basis,sgn] (or 0 if the state is null)

# Initializing psi. Psi will be a 3 indices list of size lenNPartBasis x 
# lenBasis x lenBasis.
psi = [ [ [0 for i in range(lenBasis)  ] for j in range(lenBasis) ] for 
           k in range(lenNPartBasis) ]

# Note that one could make a function out of this and call it when calculating 
# rho_1
# but that seems to make the run time longer!
# We shall fill only the lower triangular part of rho since it is Hermitian
# i.e. j<=i
for i in range(lenBasis):
    for j in range(i+1):
        for nbasis in range(lenNPartBasis): 
            basis = nPartBasis[nbasis]
            # Check if the SO we want to annihilate is in the basis state
            if j in basis:
                auxSgn   = (-1)**basis.index(j)
                auxState = basis[:]             # force pass by value
                auxState.remove(j)              # if yes, then annihilate that SO
                                                # but before doing that compute 
                                                # the sign change of destroying it.
            else:
                psi[nbasis][i][j] = [0,0]       # otherwise give 0
                continue                        # and go to the next iteration 
            # Check if the SO we want to create is not in the basis state
            if i not in auxState:
                pos = getPos(auxState,i,nEle)   # Get the position of the first element higher than i
                sgn = auxSgn*(-1)**pos
                auxState.insert(pos,i)          # Include the SO in the right position
                psi[nbasis][i][j]=[auxState,sgn]
            else:
                psi[nbasis][i][j] = [0,0]


# Now all that is left is to iterate over row and columns of the density matrix, 
# and all the basis elements making up the GS and calculate the innerProducts

# initialize rho
rho1 = zeros( (lenBasis,lenBasis) )

# We shall fill only the lower triangular part of rho since it is Hermitian
# i.e. j<=i
for i in range(lenBasis):
    for j in range(i+1):
        for nBasisBra in range(lenNPartBasis):
            for nBasisKet in range(lenNPartBasis):
                basisBra  = nPartBasis[nBasisBra]
                basisKet  = psi[nBasisKet][i][j][0]
                sgn       = psi[nBasisKet][i][j][1]
                # Grab the right coeff from the GS
                coeffBra  = vect[0][nBasisBra]
                coeffKet  = vect[0][nBasisKet]
                if innerP(basisBra,basisKet): # if the inner product is not 0 
                    rho1[i,j] = rho1[i,j] + coeffBra*coeffKet*sgn


lambd, rhovec= linalg.eigh(rho1)

print "\n rho_1 eigenvalues \n", lambd

print " tr(rho_1) = ", lambd.sum()
