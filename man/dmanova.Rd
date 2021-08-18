
\name{dmanova}
\alias{dmanova}


\title{Distance-based Multivariate Analysis of Variance (Analytical P-value Calculation)}

\description{Analysis of variance using distance matrices --- for
  partitioning distance matrices among sources of variation and fitting
  linear models (e.g., factors, polynomial regression) to distance 
  matrices; calculate the analytical p-value based on pseudo-\eqn{F} statistic 
  without permutation.}

\usage{
dmanova(formula, data = NULL, positify = FALSE,
		contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		returnG = FALSE)
}


		
\arguments{

  \item{formula}{model formula. The LHS must be a dissimilarity matrix (either class \code{matrix} or class
  \code{dist}, e.g., from  \code{\link{vegdist}} or \code{\link{dist}}.  The RHS defines the independent 
  variables. These can be continuous variables or factors, they can be transformed
    within the formula, and they can have interactions as in a typical
    \code{\link{formula}}. }

  \item{data}{the data frame for the independent variables.}
  \item{positify}{a logical value indicating whether to make the Gower's matrix  positive definite
   using the \code{\link{nearPD}} function in \code{Matrix} package. This is equivalent to modifying the distance 
   matrix so that it has an Euclidean embedding. }
    \item{contr.unordered, contr.ordered}{contrasts used for the design
    matrix (default in R is dummy or treatment contrasts for unordered
    factors). }
   \item{returnG}{a logical value indicating whether the Gower's matrix should be returned.}
}

\details{

\code{dmanova} is a permutation-free method for approximating the p-value from  distance-based permutational multivariate
analysis of variance (PERMANOVA). PERMANOVA is  slow when the sample size is large. In contrast,
\code{dmanova} provides an analytical solution, which is several orders of magnitude faster for large sample sizes. The covariate of interest
should be put as the last term in \code{formula} while the variables to be adjusted are put before the covariate of interest.
}

\value{


  Function \code{dmanova} returns a list with the following components:

  \item{aov.tab}{typical AOV table showing sources of variation,
    degrees of freedom,  sums of squares, mean squares,
    \eqn{F} statistics, partial \eqn{R^2}{R-squared} and \eqn{P}
    values. }
  \item{df}{degree of freedom for the Chisquared distribution.} 
  \item{G}{ The Gower's matrix if \code{returnG} is true.}   
  \item{call}{ the call made}

}


\references{
Chen, J. & Zhang, X. 2021.  D-MANOVA: fast distance-based multivariate analysis of variance for large-scale microbiome association studies.
 Bioinformatics. https://doi.org/10.1093/bioinformatics/btab498

}

\author{Jun Chen and Xianyang Zhang}


\examples{

data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

groups <- throat.meta$SmokingStatus

# Rarefaction
otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff

# Calculate the UniFrac distance
unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs

# Test the smoking effect based on unweighted UniFrac distance, adjusting sex
# 'Sex' should be put before 'SmokingStatus'
dmanova(as.dist(unifracs[, , 'd_UW']) ~ Sex + SmokingStatus, data = throat.meta)


}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{adonis3}}
}

\keyword{distance}
\keyword{multivariate}


