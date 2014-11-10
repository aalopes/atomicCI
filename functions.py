#!/usr/bin/python

"""
functions.py

Functions used for the calculation of the one-particle reduced density 
matrix (1-RDM) of small atomic systems. It makes use of the PyQuante 
library by Rick Muller.

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
import math
import itertools

def matrixElementH1(braNpart,ketNpart,soBasis,nucNum,h1):
    """
        Calculates matrix elements for the single particle Hamiltonian in an N 
        particle basis from a matrix with all the one particle integrals.
    """

# Calculating the integrals - see Szabo and Ostlund table 2.3 and 2.4 for this
# Maybe checking if the integrals are 0 should be the first check to speed this 
# up (since H is most probably sparse)

    ndiffSO, diffBraSO, diffKetSO, identElems, sgn = maxCoinc(braNpart,ketNpart)
    # If identical - case 1
    if ndiffSO == 0:
    # Note that we;re going to do this over al spin orbitals, but one could 
    # improve this by checking which spatial orbitals are doubly filled and 
    # calculating only one integral in those situations instead of 2
        sum = 0
        for n in braNpart:
            # Convert from orbital x spin to orbital only indices
            n = floor(float(n)/2)
            # nucNum must adapt to the atom in question if want wants to use 
            # molecules
            sum = sum + h1[n][n]
        return sum
      
    # If they differ on only one spin-orbital (one has to check this maximum
    # coincidence requirement - keep track of the permutation parity) - case 2
    if ndiffSO == 1:
        mp       = diffBraSO[0]
        pp       = diffKetSO[0]

        spinm   = soBasis[mp][0]
        spinp   = soBasis[pp][0]
                        
        # Convert from orbital x spin to orbital only indices
        m = int(floor(float(mp)/2))
        p = int(floor(float(pp)/2))

        # Check if the spin is different
        if spinInner(spinm,spinp) == 0:
            return 0
        else:
            return sgn*( h1[m][p] )
    
    # Else it's 0 - case 3
    else:
        return 0
        
# 
def matrixElementH2(braNpart,ketNpart,soBasis,h2):
    """
        Calculates matrix elements for the two particle Hamiltonian in the N 
        particles basis from a matrix with all the two particle integrals
    
    """
# Calculating the integrals - see Szabo and Ostlund table 2.3 and 2.4 for this
# Maybe checking if the integrals are 0 should be the first check to speed this 
# up (since H is most probably sparse).

# Here, specially for the spin inner product, note that [ij|kl] = <ik|jl>. Then 
# [ij|kl] = <i|j>.<k|l>


    ndiffSO, diffBraSO, diffKetSO, identElems, sgn = maxCoinc(braNpart,ketNpart)
    # If identical - case 1
    if ndiffSO == 0:
        sum = 0
        # Here one does not need to go through evey single orbital we only need 
        # to go n < m
        # Also, note that most probably m is different from n 
        # (since |K> = |...mn...>)
        for np in braNpart:
            for mp in braNpart:
                spinn   = soBasis[np][0]
                spinm   = soBasis[mp][0]
                # Convert from orbital x spin to orbital only indices
                n = int(floor(float(np)/2))
                m = int(floor(float(mp)/2))
                sum = sum+0.5*(h2[m][m][n][n] - spinInner(spinn,spinm)*h2[m][n][n][m])
        return sum
    # If they differ on only one spin-orbital (one has to check this maximum 
    # coincidence requirement - keep track of the permutation parity)
    # Then check if the two SO have different spin. If yes, this is immediately 
    # 0, if not, calculate the matrix element involving the orbital parts of the 
    # wavef. - case 2
    if ndiffSO == 1:
        mp       = diffBraSO[0]
        pp       = diffKetSO[0]
        
        spinm   = soBasis[mp][0]
        spinp   = soBasis[pp][0]
        
        # Convert from orbital x spin to orbital only indices
        m = int(floor(float(mp)/2))
        p = int(floor(float(pp)/2))

        sum = 0
        for np in identElems:
            spinn = soBasis[np][0]
            
            # Convert from orbital x spin to orbital only indices
            n = int(floor(float(np)/2))
            
            sum = sum + (sgn*(spinInner(spinm,spinp)* h2[m][p][n][n] - 
                  spinInner(spinm,spinn)*spinInner(spinn,spinp)* h2[m][n][n][p]))
                  # The sign can be introduced at the end of the loop avoiding 
                  # unnecessary multiplications
        return sum
        
    # If they differ on two spin-orbitals
    # Again the same thing, check first the spin part, if they match then 
    # proceed to the orbital part - case 3
    if ndiffSO == 2:
        mp      = diffBraSO[0]
        np      = diffBraSO[1]
        pp      = diffKetSO[0]
        qp      = diffKetSO[1]
        
        spinm   = soBasis[mp][0]
        spinn   = soBasis[np][0]
        spinp   = soBasis[pp][0]
        spinq   = soBasis[qp][0]

        # Convert from orbital x spin to orbital only indices
        m = int(floor(float(mp)/2))
        n = int(floor(float(np)/2))
        p = int(floor(float(pp)/2))
        q = int(floor(float(qp)/2))

        sum = sgn*( spinInner(spinm,spinp)*spinInner(spinn,spinq)*  
                    h2[m][p][n][q] - spinInner(spinm,spinq)* 
                    spinInner(spinn,spinp)* h2[m][q][n][p])
        return sum
    # If they differ on more than two spin-orbitals - case 4
    else:
        return 0
def maxCoinc(braNpart,ketNpart):
    """
        Checks in how many spin-orbitals they differ and return that value. 
        Return as well these spin-orbitals (the ones from the bra and the ones
        from the ket) and the sign that comes out of the transposition 
        operations necessary to put bra and ket into maximal coincidence
    """

    # Need to write what every variable represents
    ndiffSO     = 0
    diffBraSO   = []
    diffKetSO   = []
    identElems  = []   
    auxBra      = braNpart[:]
    auxKet      = ketNpart[:]
    nTransp     = 0
    sgn         = 1
    
    # If they are the same state, then skip the calculation
    if braNpart == ketNpart:
        return ndiffSO, diffBraSO, diffKetSO, identElems, sgn
    
    # Checking how many different SO there are between the bra and the ket - 
    # note that if they differ in more than 2 then the integrals we are 
    # interested in will be 0. Therefore we can skip the transposition.
    for soBra in braNpart:
        if soBra not in ketNpart:
            ndiffSO = ndiffSO + 1
            if ndiffSO == 3: # Note that if they differ in 3 or more the 
            # integrals are immediately 0, so we let ndiffSO = 3 even if 
            # they differ in more than 3 SO.
                return ndiffSO, diffBraSO, diffKetSO, identElems, sgn
        
    # This is how I'll do it: I'll keep the Bra fixed and transpose the Ket 
    # elements (while keeping track of the number of transpositions) until it is 
    # into maximal coincidence with the Bra.
    # One or two spin-orbitals will be different between the Bra and the Ket. 
    # One must then extract them and if there are two, we must keep them in the 
    # right order.
    
    # Extract the identical elements - one can do the extraction of the 
    # identical elements and the calculation of ndiffSO in one loop
    for soBra in braNpart:
        if soBra in ketNpart:
            identElems.append(soBra)
    # Now order the auxKet so it is in maximal coincidence with the Bra
    for element in identElems:
        if braNpart.index(element) == auxKet.index(element): # right position, skip
            continue
        else:
            idxBra          = braNpart.index(element)
            idxKet          = auxKet.index(element)
            # transpose elements
            auxKet[idxKet]  = auxKet[idxBra]
            auxKet[idxBra]  = element 
            sgn = -1*sgn        
    # now we only need to extract the different SO (in the right order)
    for i in range(len(braNpart)):
        soBra = braNpart[i]
        if soBra not in auxKet:
            diffBraSO.append(braNpart[i])
            diffKetSO.append(auxKet[i])
    return ndiffSO, diffBraSO, diffKetSO, identElems, sgn
        
# Calculating spin inner product - it should check the input and throw an error 
# if it's not correct (this obviously slows it down - but the integrals are the bottleneck)
def spinInner(spini,spinj):
    if spini == spinj:
        return 1
    else:
        return 0
        
def maxCoincSgn(braNpart,ketNpart):
    """
        Checks the parity of the transpositions necessary for the Slater Bra and 
        Ket to coincide
    
    """

    # Need to write what every variable represents
    ndiffSO     = 0
    diffBraSO   = []
    diffKetSO   = []
    identElems  = []   
    auxBra      = braNpart[:]
    auxKet      = ketNpart[:]
    nTransp     = 0
    sgn         = 1
    
    # If they are the same state, then skip the calculation
    if braNpart == ketNpart:
        return sgn
    
    # Cheking how many different SO there are between the bra and the ket - note 
    # that if they differ in more than 2 then the integrals we are interested in 
    # will be 0. Therefore we can skip the transposition.
    for soBra in braNpart:
        if soBra not in ketNpart:
            ndiffSO = ndiffSO + 1

        
    # Do it
    # This is how I'll do it: I'll keep the Bra fixed and transpose the Ket 
    # elements (while keeping track of the number of transpositions) until it
    # is into maximal coincidence with the Bra.

   
    # Extract the identical elements - one can do the extraction of the 
    # identical elements and the calculation of ndiffSO in one loop
    for soBra in braNpart:
        if soBra in ketNpart:
            identElems.append(soBra)
    # Now order the auxKet so it is in maximal coincidence with the Bra
    for element in identElems:
        if braNpart.index(element) == auxKet.index(element): # right position skip
            continue
        else:
            idxBra          = braNpart.index(element)
            idxKet          = auxKet.index(element)
            # transpose elements
            auxKet[idxKet]  = auxKet[idxBra]
            auxKet[idxBra]  = element 
            sgn = -1*sgn        
        return sgn


def getPos(auxState,i,nEle):
    """
        Gets the smallest j in auxState s.t j>i
        Note that auxState contains the SO indices after we removed one of them!
    """
    pos = 0
    for j in range(nEle-1):
        if auxState[j] > i:
            return j
            break
    # If we didnt find j>i, then i should go to the end of the list,
    # i.e. j+1
    return j+1

def innerP(basisBra,basisKet):
    """
        Calculates inner product between bra and ket
    """
    if basisBra == basisKet:
        return 1
    else:
        return 0
