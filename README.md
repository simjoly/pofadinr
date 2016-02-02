# pofadinr

**Distance methods from SNP data**

This is the repository for the pofadinr R package

## Installation

### Manual installation

Download either the source files or the binaries for the [latest release](https://github.com/simjoly/pofadinr/releases) and install them in R.

### Direct installation in R

You first need to install the devtools package.

$> install.packages("devtools")

And then, intall the pofadinr package from the github repository.

$> library(devtools)

$> install_github("simjoly/pofadinr")

### Source and binaries 

Check the [latest release](https://github.com/simjoly/pofadinr/releases).

## Notes

This package was written to make available in R some SNP based distance methods (Joly et al. 2015). However, more are available from the POFAD software.

This implementation is however much faster than the original program, which is possible by using the bit level operations for DNA sequences developped by Emmanuel Paradis.

## References

Joly, S., Bryant, D. and Lockhart, P.J. 2015. Flexible methods for estimating genetic distances from nucleotide data. *Methods in Ecology and Evolution*, 6, 938–948.