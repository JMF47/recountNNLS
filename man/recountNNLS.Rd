% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recountNNLS.R
\name{recountNNLS}
\alias{recountNNLS}
\title{NNLS for Tx Abundance Calculation}
\usage{
recountNNLS(pheno, jx_file = NULL, cores = 1)
}
\arguments{
\item{pheno}{The table of phenotype information from processPheno().}

\item{jx_file}{The path to the Rail-RNA junction coverage file if processing sample not already in recount2.}

\item{cores}{The number of processing cores to use.}
}
\value{
Returns an rse object of the estimated number of reads and the associated standard errors.
Each row represents a protein-coding gene, and each column represents a sample in the phenotype matrix.
}
\description{
This function applies NNLS to calculate the transcript abundance
from coverage of exonic and junction coverage statistics for samples
annotated in a phenotype matrix created using processPheno().
}
\examples{
project = 'SRP063581'
pheno = processPheno(project)
}
\keyword{recountNNLS}
