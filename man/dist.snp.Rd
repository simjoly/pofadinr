\name{dist.snp}
\alias{dist.snp}
\title{SNP distances from nucleotide data}
\description{
  This function implements different distance-based methods for 
  estimating genetic distances from nucleotide data.

  Its main interest resides in using the information of the IUPAC
  ambiguitiy nucleotides to get accurate distance between individuals
  that could show polymorphisms at certain site positions in a DNA
  sequence because of allelic variation. These polymorphisms need 
  to be coded using the ambiguity nucleotides codes (see 
  \code{consensus.dna}). For instance, the presence of the ambiguity
  nucleotide 'R' at a position indicates that this individual has 
  both an 'A' and a 'G' at this position in its genome, possibly 
  because of the presence of two alleles at this locus. These methods 
  are particularly useful for comparing heterozygous individuals and 
  comparing individuals of different ploidy levels.
}
\usage{
dist.snp(x, model = "GENPOFAD", pairwise.deletion = TRUE, 
    as.matrix = FALSE)
}
\arguments{
  \item{x}{a nucleotide alignment of class \code{"DNAbin"}.}
  \item{model}{a character string specifying the distance method used.
    Two choices are available: \code{"GENPOFAD"}, \code{"MATCHSTATES"}.}
  \item{pairwise.deletion}{a logical indicating whether to delete the
    sites with no nucleotide data in a pairwise way (DEFAULT). The 
    alternative is to delete the alignment sites that do not have 
    any nucleotide information (A, T, G, C; or a combination of these: 
    Y, R, W, ..., N). In practice, setting to FALSE will delete sites
    with at least one '-' or one '?'.}
  \item{as.matrix}{a logical indicating whether to return the results as
    a matrix. The default is to return an object of class
    \link[stats]{dist}.}
}
\details{
  Several methods exist for estimating distances from DNA sequences, 
  several of which are available in the \code{"dist.dna"} function of the 
  \pkg{ape} package. The methods available here are different by
  using information on all nucleotides present at an alignment site in
  an individual to estimate distances. To do this, the 
  polymorphisms at the DNA level need to be coded using the IUPAC 
  ambiguity nucleotides. For instance, the presence of the ambiguity
  nucleotide 'R' at a position indicates that this individual has 
  both an 'A' and a 'G' at this position in its genome, possibly 
  because of the presence of two alleles at this locus. 

  Two methods are available, \code{GENPOFAD} and \code{MATCHSTATES}, 
  which are described by Joly et al. (2015). These methods are estimated 
  at the nucleotide level and not at the gene level. The overall distance 
  is the mean distance of all alignment sites.

  Following Joly et al. (2015), let \eqn{A_X^i}{A_X_i} be the complete
  set of nucleotides for individual \eqn{X} at site \eqn{i} and let 
  \eqn{|A_X^i|}{|A_X_i|} be the number of nucleotide states observed for 
  individual \eqn{X} at site \eqn{i}.

\itemize{
  \item{\code{GENPOFAD}: }{
  The \code{genpofad} distance can be defined as one minus the ratio of
  the number of nucleotides shared between two individuals divided by
  the maximum number of nucleotides observed in either of the individuals
  at a given sequence site. Then, the \code{genpofad} distance between
  individual \eqn{X} and individual \eqn{Y} at site \eqn{i} is

    \deqn{\code{genpofad}_{XY}^i:=1-\frac{|A_X^i \cap A_Y^i |}{max(| A_X^i |,| A_Y^i |)}}{%
          GENPOFAD_XY_i := 1 - (|A_X_i UNION A_Y_i| / max(|A_X_i|,|A_Y_i|))}
  }

  \item{\code{MATCHSTATES}: }{
    Let \eqn{A_X^i \Delta A_Y^i}{|A_X_i| EXCLUSIVE_OR |A_Y_i|} denote
    the set of nucleotides that belong to either \eqn{A_X^i}{A_X_i}
    or \eqn{A_Y^i}{A_Y_i}, but not both. The \code{matchstates} distance
    between individual \eqn{X} and individual \eqn{Y} at site \eqn{i} is

    \deqn{\code{matchstates}_{XY}^i:=\frac{|A_X^i \Delta A_Y^i |}{|A_X^i |+|A_Y^i|}}{%
        MATCHSTATES_XY_i := (|A_X_i| EXCLUSIVE_OR |A_Y_i|) / |A_X_i |+|A_Y_i|}
  }
}

The overall distance is then:

\deqn{\frac{1}{s} \sum_{i=1}^{s}d_{XY}^{i}}{ 1/s \sum_s d_XY_i}

Please see Joly et al. (2015) for a comparison of these methods with other 
methods that exist to use nuclotide polymorphisms. Although \code{genpofad} 
and \code{matchstates} are probably the most useful methods (Joly et al. 2015),
additionnal methods and more options are available in the POFAD software
(Joly et al. 2006; Joly 2014).
}

\note{
  This function was developped using code from the \code{dist.dna} function
  of the \pkg{ape} package by Emmanuel Paradis. 
}

\value{
  an object of class \link[stats]{dist} (by default), or a numeric
  matrix if \code{as.matrix = TRUE}.
}

\references{
  Joly, S. and Bruneau, A. 2006. Incorporating allelic variation for 
  reconstructing the evolutionary history of organisms from multiple 
  genes: an example from Rosa in North America. \emph{Systematic 
  Biology}, \bold{55}, 623--636.

  Joly, S. 2014. \emph{POFAD: Phylogenetic Analysis From Allelic Data}.
  Distributed by the author. https://github.com/simjoly/pofad/releases. 

  Joly, S., Bryant, D. and Lockhart, P.J. 2015. Flexible methods for
  estimating genetic distances from nucleotide data. \emph{Methods in
  Ecology and Evolution}, \bold{6}, 938--948.
}

\author{Simon Joly, Emmanuel Paradis}

\examples{
# Example data alignment
cat("2 11",
    "A     TRTACGTTGAC",
    "B     AYTAGTTAGWR",
    file = "exdna.txt", sep = "\n")
# Read the alignment
require(ape)
ex.dna <- read.dna("exdna.txt", format = "sequential")
# Compute genpofad distance
d <- dist.snp(ex.dna, model="GENPOFAD", pairwise.deletion = TRUE)
d
}
