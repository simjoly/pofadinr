\name{consensus.dna}
\alias{consensus.dna}
\title{Consensus DNA sequences from DNA alignments}
\description{
  This function estimates the sequence consensus of an alignment
  using the IUPAC ambiguitiy nucleotides.
}
\usage{
consensus.dna(x)
}
\arguments{
  \item{x}{a nucleotide alignment of class \code{"DNAbin"}.}
}
\details{
  Computes consensus sequences. Note that gaps are presently treated as
  missing characters only and as such they are ignored. Also, sites
  identified as polymorphic for a nucleotides are flagged as
  'not known' according to the DNAbin standards, even if both of these
  nucleotides can be present in a given organism. 
}

\value{
  An object of class \link{DNAbin} with a single sequence.
}

\author{Simon Joly}

\seealso{
  \code{\link{dist.snp}}
}

\examples{
# Example data alignment
cat("2 11",
    "A     TRTACGTTGAC",
    "B     AYTAGTTAGWR",
    file = "exdna.txt", sep = "\n")
# Read the alignment
require(ape)
ex.dna <- read.dna("exdna.txt", format = "sequential")
con <- consensus.dna(ex.dna)
as.character(con)
}
