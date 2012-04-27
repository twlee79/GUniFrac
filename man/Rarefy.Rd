\name{Rarefy}
\alias{Rarefy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rarefy the OTU table to an equal sequencing depth
}
\description{
GUniFrac is also sensitive to different sequencing depth. To compare
microbiomes on an equal basis, rarefaction might be used.
}
\usage{
Rarefy(otu.tab, depth = min(rowSums(otu.tab)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu.tab}{
OTU count table, row - n sample, column - q OTU
}
  \item{depth}{
Required sequencing depth; If not specified, the lowest sequencing
depth is used.
}
}

\value{
 Return a LIST containing:
 \item{otu.tab.rff}{Rarefied OTU table}
 \item{discard}{IDs of samples that do not reach the specified sequencing depth}
%% ...
}
\references{
Jun Chen and Hongzhe Li(2012). Associating microbiome composition with 
environmental covariates using generalized UniFrac distances. (Submitted)
}
\author{
Jun Chen <chenjun@mail.med.upenn.edu>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GUniFrac}}, \code{\link{PermanovaG}}
}
\examples{
data(throat.otu.tab)
# Rarefaction
otu.tab.rff <- Rarefy(throat.otu.tab, 1024)$otu.tab.rff
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Rarefaction }
\keyword{ Ecology }% __ONLY ONE__ keyword per line
