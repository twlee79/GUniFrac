\name{vaginal.otu.tab}
\alias{vaginal.otu.tab}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Vaginal Microbiome OTU Count Table
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
OTU count table from 16S V3-V5 targeted sequencing of the vaginal microbiome samples from the HMP project. A total of 780 OTUs from 381 samples.
}
\usage{data(vaginal.otu.tab)}

\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The OTU table was taken from R bioconductor "HMP16SData" package. OTUs with prevalence less than 10\% and maximum proportion less than 0.2\% were removed. 
This OTU table can be used for simulating vaginal microbiome sequencing data.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Schiffer L, Azhar R, Shepherd L, Ramos M, Geistlinger L, Huttenhower C, Dowd JB, Segata N, Waldron L (2019). 
“HMP16SData: Efficient Access to the Human Microbiome Project through Bioconductor.” American Journal of Epidemiology. doi: 10.1093/aje/kwz006.
}

\examples{
data(stool.otu.tab)
}
\keyword{datasets}
