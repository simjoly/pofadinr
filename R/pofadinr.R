
dist.snp <- function(x, model = "GENPOFAD", pairwise.deletion = TRUE, 
                      as.matrix = FALSE)
{
  
    MODELS <- c("GENPOFAD", "MATCHSTATES")
    imod <- pmatch(toupper(model), MODELS)
    if (is.na(imod))
        stop(paste("'model' must be one of:",
                   paste("\"", MODELS, "\"", sep = "", collapse = " ")))

    if (is.list(x)) x <- as.matrix(x)
    nms <- dimnames(x)[[1]]
    n <- dim(x)
    s <- n[2]
    n <- n[1]

    if (!pairwise.deletion) {
        keep <- .C("GlobalDeletionDNA", x, n, s, rep(1L, s))[[4]]
        x <- x[, as.logical(keep)]
        s <- dim(x)[2]
    }
    Ndist <- n*(n - 1)/2
    d <- .C("pofad", x, as.integer(n), as.integer(s), imod,
            double(Ndist), as.integer(pairwise.deletion))
    d <- d[[5]]
    attr(d, "Size") <- n
    attr(d, "Labels") <- nms
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    attr(d, "method") <- model
    class(d) <- "dist"
    if (as.matrix) d <- as.matrix(d)
    d
}

consensus.dna <- function(x)
{
    # Compute the consensus sequence of the alignment given as a matrix
    # given in DNAbin format

    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)
    s <- n[2]
    n <- n[1]

    d <- .C("consensusDNA", x, raw(s), as.integer(n), as.integer(s))

    con <- matrix(d[[2]],ncol=s)
    class(con) <- "DNAbin"
    con
}

consensus.org <- function(x, mat)
{
  # Given a DNA alignment in DNAbin format and 
  # a matrix assigning alleles to organisms,
  # the function computes the consensus sequence of all
  # alleles for each organisms
  
  if (is.list(x)) x <- as.matrix(x)
  nms <- dimnames(x)[[1]]
  s <- dim(x)[2]
  
  if(!is.data.frame(mat)) mat <- as.data.frame(mat)
  
  if (!is.factor(mat[,1]) || !is.factor(mat[,2]))
    stop('Only dataframe with factors are accepted (no numeric columns)')
  
  alleles.missing <- mat[!(mat[,2] %in% nms),2]
  if (length(alleles.missing) > 0)
    stop("The following alleles are not found in the DNA matrix\n",
         paste(alleles.missing,collapse="\n"))
  
  #order alleles in data.frame as in alignment
  mat <- mat[match(nms,mat[,2]),]
  
  org <- levels(mat[,1])
  con <- matrix(raw(),nrow=length(org),ncol=s,dimnames=list(org,NULL))
  class(con) <- "DNAbin"
  
  for (i in 1: length(org)) {
    newx <- x[mat[,1]==org[i],]
    d <- .C("consensusDNA", newx, raw(s), as.integer(nrow(newx)), as.integer(s))
    con[i,] <- d[[2]]
  }
  con
}
