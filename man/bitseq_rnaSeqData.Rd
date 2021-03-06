% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bitseq_rnaSeqData.R
\name{bitseq_rnaSeqData}
\alias{bitseq_rnaSeqData}
\title{Obtaining data in the format of example RNAseqDATA by using BitSeq output files}
\usage{
bitseq_rnaSeqData(t, trFileName)
}
\arguments{
\item{t}{Vector which contains the input values, i.e., time points. The file names 
for the corresponding time point is specified as the names of this vector.}

\item{trFileName}{Name of the transcriptome file.}
}
\value{
List of GP data within the right structure.
}
\description{
Function for obtaining the means and standard deviations at the available time points.
}
\examples{
# t=log(c(0,5,10,20,40,80,160,320,640,1280)+5)
# names(t)=c("t0000.rpkm","t0005.rpkm","t0010.rpkm","t0020.rpkm","t0040.rpkm",
# "t0080.rpkm","t0160.rpkm","t0320.rpkm","t0640.rpkm","t1280.rpkm")
# trFileName="example_tr"
# bitseq_rnaSeqData(t,trFileName)

}
\author{
Hande Topa, \email{hande.topa@helsinki.fi}
}
\keyword{GP}
\keyword{deviation}
\keyword{mean}
\keyword{standard}
