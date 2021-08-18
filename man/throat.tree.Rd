\name{throat.tree}
\alias{throat.tree}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
UPGMA Tree of Throat Microbiome OTUs
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The OTU tree is constructed using UPGMA on the K80 distance matrice of
the OTUs. It is a rooted tree of class "phylo".
}
\usage{data(throat.tree)}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The OTUs are produced by the QIIME software. Singleton OTUs have
been discarded.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010) 
Disordered Microbial Communities in the Upper Respiratory Tract of 
Cigarette Smokers. PLoS ONE 5(12): e15216. 
}

\examples{
data(throat.tree)
}
\keyword{datasets}
