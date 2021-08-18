
\name{adonis3}
\alias{adonis3}


\title{Permutational Multivariate Analysis of Variance Using Distance Matrices (Freedman-Lane permutation)}

\description{Analysis of variance using distance matrices --- for
  partitioning distance matrices among sources of variation and fitting
  linear models (e.g., factors, polynomial regression) to distance 
  matrices; uses a permutation test (Freedman-Lane permutation) with pseudo-\eqn{F} ratios.}

\usage{
adonis3(formula, data, permutations = 999, method = "bray",
    strata = NULL, contr.unordered = "contr.sum",
    contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)
}

\arguments{

  \item{formula}{model formula. The LHS must be either a community
    data matrix or a dissimilarity matrix, e.g., from
    \code{\link{vegdist}} or \code{\link{dist}}.  If the LHS is a data
    matrix, function \code{\link{vegdist}} will be used to find the
    dissimilarities. The RHS defines the independent variables. These
    can be continuous variables or factors, they can be transformed
    within the formula, and they can have interactions as in a typical
    \code{\link{formula}}. }

  \item{data}{ the data frame for the independent variables.}
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{method}{ the name of any method used in \code{\link{vegdist}} to
    calculate pairwise distances if the left hand side of the
    \code{formula} was a data frame or a matrix. }
  \item{strata}{ groups (strata) within which to constrain permutations.  }
  \item{contr.unordered, contr.ordered}{contrasts used for the design
    matrix (default in R is dummy or treatment contrasts for unordered
    factors). }
  \item{parallel}{number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{\dots}{Other arguments passed to \code{vegdist}.}
}

\details{

\code{adonis3} is the re-implementation of the famous \code{adonis} function in
the \code{vegan} package based on the Freedman-Lane permutation scheme.
(Freedman & Lane (1983), Hu & Satten (2020)).  \code{adonis} is the function 
for the analysis and partitioning sums of squares using dissimilarities. The original
implementation in the \code{vegan} package is directly based on the algorithm of Anderson (2001) and
performs a sequential test of terms. Statistical significance is calculated
based on permuting the distance matrix.  As shown in Chen & Zhang (2020+), such permutation
will lead to power loss in testing the effect of a covariate of interest while adjusting
for other covariates (confounders). The power loss is more evident  when the confounders' effects
are strong, the correlation between the covariate of interest and the confounders is high, and
the sample size is small. When the sample size is large than 100, the difference is usually small.
The new implementation is revised on the \code{adonis} function with the same interface.
}

\value{


  Function \code{adonis3} returns an object of class \code{"adonis"} with
  following components:

  \item{aov.tab}{typical AOV table showing sources of variation,
    degrees of freedom, sequential sums of squares, mean squares,
    \eqn{F} statistics, partial \eqn{R^2}{R-squared} and \eqn{P}
    values, based on \eqn{N} permutations. }
  \item{coefficients}{ matrix of coefficients of the linear model, with
    rows representing sources of variation and columns representing
    species; each column represents a fit of a species abundance to the
    linear model. These are what you get when you fit one species to
    your predictors. These are NOT available if you supply the distance
    matrix in the formula, rather than the site x species matrix} 
  \item{coef.sites}{ matrix of coefficients of the linear model, with
    rows representing sources of variation and columns representing
    sites; each column represents a fit of a sites distances (from all
    other sites) to the  linear model. These are what you get when you
    fit distances of one site to
    your predictors. }   
  \item{f.perms}{ an \eqn{N} by \eqn{m} matrix of the null \eqn{F}
    statistics for each source of variation based on \eqn{N}
    permutations of the data. The permutations can be inspected with
    \code{\link{permustats}} and its support functions.}
  \item{model.matrix}{the \code{\link{model.matrix}} for the right hand
    side of the formula.}
  \item{terms}{the \code{\link{terms}} component of the model.}
}



\references{
Anderson, M.J. 2001. A new method for non-parametric multivariate
analysis of variance. \emph{Austral Ecology}, \strong{26}: 32--46.

Freedman D. & Lane D. 1983. A nonstochastic interpretation of reported significance
levels. \emph{Journal of Business and Economic Statistics}, \strong{1}292--298.

Hu, Y. J. & Satten, G. A. 2020. Testing hypotheses about the microbiome using the 
linear decomposition model (LDM). Bioinformatics.


}

\author{Martin Henry H. Stevens (\code{adonis}) and Jun Chen
  (\code{adonis3}). }


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
adonis3(as.dist(unifracs[, , 'd_UW']) ~ Sex + SmokingStatus, data = throat.meta)


}

\keyword{distance}
\keyword{multivariate}


