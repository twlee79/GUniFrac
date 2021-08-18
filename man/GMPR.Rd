\name{GMPR}
\alias{GMPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geometric Mean of Pairwise Ratios (GMPR) Normalization for Zero-inflated Count Data
}
\description{
 A robust normalization method for zero-inflated count data such as microbiome sequencing data. 
}
\usage{
GMPR(OTUmatrix, min_ct = 2, intersect_no = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OTUmatrix}{
An OTU count table, where OTUs are arranged in rows and samples in columns.
}
  \item{min_ct}{
The minimal number of OTU counts. Only those OTU pairs with at least \code{min_ct} counts are 
considered in the ratio calculation.   The default is 2.
}
  \item{intersect_no}{
The minimal number of shared OTUs between samples. Only those sample pairs sharing at least \code{intersect_no} OTUs 
are considered in geometric mean calculation. The default is 4.
}
}
\details{
Normalization is a critical step in microbiome sequencing data analysis to account for variable library sizes.
Microbiome data contains a vast number of zeros, which makes the traditional RNA-Seq normalization methods unstable. 
The proposed  GMPR normalization remedies this problem by switching the two steps in DESeq2 normalization:


First, to calculate rij, the median count ratio of nonzero counts between samples:
  rij=median(cki/ckj) (k in 1:OTU_number and cki, ckj is the non-zero count of the kth OTU)

Second, to calculate the size factor si for a given sample i:
  si=geometric_mean(rij)
  
}
\value{
A vector of GMPR size factor for each sample.
}
\references{
Li Chen, James Reeve, Lujun Zhang, Shenbing Huang, and Jun Chen. 2018. GMPR: A robust normalization method for zero-inflated count
data with application to microbiome sequencing data. PeerJ, 6, e4600.
}
\author{
Jun Chen and Lujun Zhang
}


\examples{
data(throat.otu.tab)
size.factor <- GMPR(t(throat.otu.tab))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalization}% use one of  RShowDoc("KEYWORDS")
\keyword{ microbiome }% __ONLY ONE__ keyword per line
