#!/bin/bash

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

echo -e "This software is intended to be used for the calculation of the one particle reduced density matrix (1-RDM) of small atomic systems.\n
It makes use of the Python library PyQuante by Rick Muller.\n"

echo -e "Please enter the atomic number of the element for which you wish to calculate the 1-RDM (< 5).\n"

while [[ 1 ]]; do
    read element
    if  [[ "$element" == "1" ]] || [[ "$element" == "2" ]] || [[ "$element" == "3" ]] || [[ "$element" == "4" ]] ; then
        echo -e "Type the name of the basis you want to use.\n"
        echo -e  "Some valid examples are: 6-31g** or STO-3G. You can also use your custom basis if you have inserted it into PyQuante"
        read basis
        ./atomicCI.py $element $basis
        break
    else
        echo "Invalid atomic number!"
        break
    fi   

done