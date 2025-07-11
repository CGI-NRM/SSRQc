% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{WriteIndFile}
\alias{WriteIndFile}
\title{Creates/Writes excel workbook with pre-specified header read from file
combined with data from an excel file with sample data in specific format}
\usage{
WriteIndFile(
  headerInfoFile,
  importSampleFile,
  sheet = 1,
  startRow = 6,
  outputFile = "toimportIndfinal.xlsx"
)
}
\arguments{
\item{headerInfoFile}{excel file with headers to be part of output}

\item{importSampleFile}{full path to excel file with sample data that
will be appended to the headerInfoFile}

\item{sheet}{which sheet in the excel file to write to}

\item{startRow}{integer corresponding to row in excel that data will be written to}

\item{outputFile}{file name to write data to. Should have the file ending .xlsx}
}
\value{
Summary of the data written to file and an excel file is
written to disk as a side effect.
}
\description{
Creates/Writes excel workbook with pre-specified header read from file
combined with data from an excel file with sample data in specific format
}
