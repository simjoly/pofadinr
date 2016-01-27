\name{consensus.org}
\alias{consensus.org}
\title{Consensus sequences from DNA alignments for several organisms}
\description{
  This function estimates the sequence consensus of several
  organisms from a DNA alignment using the IUPAC ambiguitiy
  nucleotides.
}
\usage{
consensus.org(x, mat)
}
\arguments{
  \item{x}{a nucleotide alignment of class \code{"DNAbin"}.}
  \item{mat}{a \code{data.frame} that contains the organisms
    for which consensus sequences will be estimated in the
    first column, and the alleles names that correspond 
    \emph{exactly} to the sequence names in the alignment 
    \code{x}. Note that several alleles can be associated to
    the same organism.}
}
\details{
  Computes consensus sequences for a list of organisms (individuals).
  All the alleles of the organisms listed in the data.frame \code{mat}
  need to be present in the matrix.

  Note that gaps are presently treated as missing characters
  only and as such they are ignored. Also, sites
  identified as polymorphic for a nucleotides are flagged as
  'not known' according to the DNAbin standards, even if both of these
  nucleotides can be present in a given organism. 
}

\value{
  A nucleotide alignment of class \link{DNAbin} where each row
  contains the consensus sequence of an organism listed in the input
  data frame \code{mat}.
}

\author{Simon Joly}

\seealso{
  \code{\link{consensus.dna}} 
}

\examples{
# Example data alignment
cat("4 11",
    "A     TRTACGTTGAC",
    "B     AYTAGTTAGWR",
    "C     AGTAATATGSC",
    "D     ARWGGTADGWT",
    file = "exdna.txt", sep = "\n")
# Read the alignment
require(ape)
ex.dna <- read.dna("exdna.txt", format = "sequential")
# Matrix assigning alleles to organisms
mat <- data.frame(organisms=rep(c("org1","org2"),each=2),
                  alleles=c("A","B","C","D"))
mat
# Compute the consensus sequences for the two organisms
conseq <- consensus.org(ex.dna,mat)
as.character(conseq)
}
