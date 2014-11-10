Use this script to parse Gaussian basis set files so they can be used in PyQuante. 

One should redirect the output to a file *.py, put it in PyQuante's Basis/ folder and edit Tools.py (in the same folder) and add the name of the output file.

An example would be:

 ./parseGaussian basis.inp > myset.py

Then, find out PyQuante's folder by typing in the Python Shell:

 import sys
 sys.path

an example output would be 

 /usr/local/lib/python2.7/dist-packages/PyQuante-1.6.4-py2.7-linux-x86_64.egg

Copy myset.py to the /Basis folder of PyQuante, and edit the script Tools.py (in the same folder). You should add a new pair of key-value to the dictionary basis_map. Something like

 'NameOfMySet':'myset'
    
You can then use this set in PyQuante by choosing NameOfMySet as the basis set's name.

A. A. Lopes
Universitaet Freiburg
2012
