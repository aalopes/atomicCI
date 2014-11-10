#!/usr/bin/python

"""
Script for parsing Gaussian input files (obtained from the EMSL website for 
example) and converting them into PyQuante's (http://pyquante.sourceforge.net/) 
format.
It reads a file from the stdin and outputs to the stdout (so you should redirect 
it to a file). This file should then be put in the appropriate folder. Read the 
README file for more instructions.

This script works for all "normal" orbitals, S, P, D, ... and SP hybrid orbitals.

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

import fileinput # for reading from the stdin
import re        # regular expressions package

# Problems here: detecting if it's an element name or an orbital!!
# Must solve this thing above!! - need a flag
# initialization
basis     = {}
comment   = re.compile('[*!]{1,}')     # so we can search for * and ! characters 
atmLine   = re.compile('[A-Za-z]{1,}') # so we can search for something containing one letter
orbitLine = re.compile('[A-Za-z]{1,}')

def name2no(element):
    # takes an element name (string) and spits out the atomic number
    table = {
    "H" : 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B" : 5,
    "C" : 6,
    "N" : 7,
    "O" : 8,
    "F" : 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
	"P":  15,
	"S":  16,
	"Cl": 17,
	"Ar": 18,
	"K":  19,
	"Ca": 20,
	"Sc": 21,
	"Ti": 22,
	"V":  23,
	"Cr": 24,
	"Mn": 25,
	"Fe": 26,
	"Co": 27,
	"Ni": 28,
	"Cu": 29,
	"Zn": 30
    }
    return table.get(element)

for line in fileinput.input():
    # split line string where there are blank spaces
    sline = line.split()
    if not sline:                           # in case it's a blank line
        continue
    if comment.search(line):                # in case the line has some * or ! characters
        continue
    if atmLine.search(line) and len(line.split()) == 2:            
                                            # if there is at least one letter and the line has 2 elements
                                            # then we found the element line
        element = line.split()[0]
        atno    = name2no(element)           # convert from element name to atomic number
        basis[atno] = []                     # add key and value to dictionary
    # if not, we should have found the orbital letter
    elif orbitLine.search(line) and len(line.split()) == 3:
                                             # if there's at least one letter and the line has 3 elements
                                             # then we found the orbital line 
        sline       = line.split()
        symbol      = sline[0]               # the orbital symbol  
        if symbol == "SP":                   # then we have found an SP orbital and these are a little bit different...
            spflag = 1
            nprim       = sline[1]            # number of primitives - not used
            basis[atno].append(("S",[]))      # for the S orbital
            basis[atno].append(("P",[]))      # for the P orbital

        else:
            spflag = 0
            nprim       = sline[1]               # number of primitives - not used
            basis[atno].append((symbol,[]))
    # if not, we have a line of exponents and coeff, if spflag = 0 we have a normal orbital
    elif spflag == 0:
        sline       = line.split()
        exp         = float(sline[0])
        coeff       = float(sline[1])
        basis[atno][-1][1].append((exp,coeff))
    # if not, we have a line of exponents and coeff, if spflag = 1 we have an SP orbital
    elif spflag == 1:
    # this is definitively not finished!
        sline       = line.split()
        exp         = float(sline[0])
        coeffS       = float(sline[1])
        coeffP       = float(sline[2])
        basis[atno][-2][1].append((exp,coeffS))
        basis[atno][-1][1].append((exp,coeffP))

print("basis_data = \\")
print(basis)
fileinput.close()