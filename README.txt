atomicCI v0.1 (22.11.2012)

Copyright 2012 Alexandre Lopes

-------------------------------------------------------------------------------
= Description =

atomicCI is a software for the calculation of one particle reduced density
matrices (1-RDMs) for small atomic systems (Z <= 4) using Full CI.

This software makes use of the PyQuante-1.6.4 Python library by Rick Muller
<http://pyquante.sourceforge.net/>.

-------------------------------------------------------------------------------
= How to use =

Run the atomic.sh Bash script and follow the on-screen instructions in order
to choose the atomic element and the basis set used for the calculation.

In case one needs to add more basis functions, one can make use of the
associated Python script parseGauss.py on <https://github.com/aalopes>.

-------------------------------------------------------------------------------
= Note =

This software is slow and not really efficient. It was produced as part of my
doctoral research and served its purpose. I have no reasons to optimize it, and
will not do so in the future.

-------------------------------------------------------------------------------

= License =

All code is licensed under the New BSD Licence.

See the COPYING file for the license text.

This software is provided "as is", without any warranty.
-------------------------------------------------------------------------------