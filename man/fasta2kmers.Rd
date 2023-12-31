% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fasta2kmers}
\alias{fasta2kmers}
\alias{fasta2alpha}
\title{Create random initial PWM guess from fasta.}
\usage{
fasta2alpha(fasta, k)
}
\arguments{
\item{fasta}{Dataset of sequences.}

\item{k}{Size of kmers.}
}
\value{
The best kmers from dataset with respect to alpha.
}
\description{
Create random initial PWM guess from fasta.
}
