% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateInputData.R
\name{CreateMatchInput}
\alias{CreateMatchInput}
\title{Merge genotype data with metadata to create match file}
\usage{
CreateMatchInput(
  genotypeResults,
  filt = "No",
  level = 4,
  metaData = "rovbasemetadata.xlsx",
  confirmedDead = "No",
  colnamesWithSEP = "StrekkodeProve"
)
}
\arguments{
\item{genotypeResults}{dataframe with genotype results and unique
names of each sample}

\item{filt}{parameter that turn on and off filtering of data. If
set to No/NO/no the returned data is not filtered if set to
Yes/YES/yes samples with more than the missing data option level}

\item{level}{the maximum number of markers with missing data}

\item{metaData}{path and filename to the file with
metadata. Assumes that the file can be imported using LoadData
function and that the file contains sample name (Strekkode, SEP),
DNAid,  date of sampling and GPS coordinates.}

\item{confirmedDead}{is the sample confirmed dead or not. Allowed values
Yes/YES/yes or No/NO/no}

\item{colnamesWithSEP}{name of the column in metaData that contains
SEP number eg. SEP##### or M####}
}
\value{
dataframe with genetic and metadata merged
}
\description{
Merges genotype results with metadata creating a
dataframe suitable for matching genotype results with available
genetic database.
}
