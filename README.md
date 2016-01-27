pofadinr
====================

Distance methods from SNP data
------------------------------

This is the repository for the pofadinr R package

Source and binaries 
-------------------

Installation
------------

In R, type the following:

$> install.packages("pofadinr", repos = "http://www.omegahat.org/R", type="source")

Alternatively, you can download the binaries for your platform in install them in R.

Notes
-----

This package was written to make available in R some SNP based distance methods (Joly et al. 2015). However, more are available from the POFAD software.

This implementation is however much faster than the original program, which is possible by using the bit level operations for DNA sequences developped by Emmanuel Paradis.