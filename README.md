# pofadinr

**Distance methods from SNP data**

This is the repository for the pofadinr R package

## About

The **pofadinr** R package implements the distances **genpofad** and **matchstates** based on single nucleotide polymorphisms (oly et al. 2015). They were developped to take into account intra-individual polymorphisms when estimating genetic distances between individuals. **genpofad**, in particular, provides an accurate measure of genetic distance and is useful to deal with hybrids and individuals of different ploidy levels.

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

This package was written to make available in R some SNP based distance methods (Joly et al. 2015). However, more methods are available in the [POFAD](http://github.com/simjoly/pofad) software.

Note that this implementation is much faster than that of the POFAD program. This is because pofadinr uses bit-wise operations for comparing nucleotides, following the [bit-level coding scheme for nucleotides](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html) developped by Emmanuel Paradis.

## References

Joly, S., Bryant, D. and Lockhart, P.J. 2015. Flexible methods for estimating genetic distances from nucleotide data. *Methods in Ecology and Evolution*, 6, 938–948.