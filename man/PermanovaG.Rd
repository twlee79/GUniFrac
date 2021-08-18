\name{PermanovaG}
\alias{PermanovaG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Permutational Multivariate Analysis of Variance Using Multiple Distance Matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
In practice, we do not know a priori which type of change happens in 
the microbiome. Each distance measure is most powerful in detecting
only a certain scenario. When multiple distance matrices are available,
separate tests using each distance matrix will lead to loss of power 
due to multiple testing correction. Combing the distance matrices
in a single test will improve power. PermanovaG combines multiple distance
matrices  by taking the minimum of the P values for individual distance
matrices. Significance is assessed by permutation.
}
\usage{
PermanovaG(formula, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
a formula, left side of the formula (\code{Y ~ X}) is a three dimensional ARRAY
containing the supplied distance matrices as produced by \code{\link{GUniFrac}}
 function. Or it could be a list of distance matrices.
}
  \item{data}{
%%     ~~Describe \code{dat} here~~
a data frame containing the covariates
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
parameter passing to \code{adonis} function
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a list containing:
\item{p.tab}{a data frame, columns: p-values for individual distance matrices and the omnibus test, 
rows: covariates. (Note: they are sequential p-values, put the variable of interest in the end)}
\item{aov.tab.list}{a list of \code{adonis} AOV tables for individual distance matrices}

}

\references{
%% ~put references to the literature/web site here ~
Chen, J., Bittinger, K., Charlson, E.S., Hoffmann, C., Lewis, J., Wu, G.D., Collman, R.G., Bushman, F.D. and Li, H.(2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
}
\author{
%%  ~~who you are~~
Jun Chen <chen.jun2@mayo.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{Rarefy}}, \code{\link{GUniFrac}}
}
\examples{
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

groups <- throat.meta$SmokingStatus

# Rarefaction
otu.tab.rff <- Rarefy(throat.otu.tab)$otu.tab.rff

# Calculate the UniFracs
unifracs <- GUniFrac(otu.tab.rff, throat.tree, alpha=c(0, 0.5, 1))$unifracs

#	Combine unweighted and weighted UniFrac for testing
PermanovaG(unifracs[, , c("d_1", "d_UW")] ~ groups)
# Combine d(0), d(0.5), d(1) for testing
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ regression }% __ONLY ONE__ keyword per line
