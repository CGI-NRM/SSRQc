#' @title Automate tabular data import from different file formats
#'
#' @description Import of tabular data from different file formats
#'     requires use of different packages and syntax. Here we use file
#'     endings to guess the file format and select the most suitable
#'     import tool to generate a dataframe with nyntactically valid column
#'     names. Currently supports excel files (.xlsx, .xls), tab
#'     delimited files (.tsv, .tdf, .txt), open document format (.ods)
#'     and comma separated files (.csv). Column names are altered so the
#'     dataframe will only contain latin-ascii characters even if the read
#'     file contains special characters
#' @param filePath name of file to import data from (with path)
#' @param sheet which sheet to read data from in excel and ods file
#'     formats
#' @param naStrings entries should be converted to NA values on import
#' @param notAllowed entries that is not valid input and if found will
#'     inhibit import
#'     __NB! The function uses file endings to guess
#'     file formats. Will not work if file ending__
#'     __is not consistent with actual file type__
#' @return a dataframe of the imported data with column names that only
#'     have latin-ascii characters.
#' @importFrom readxl read_excel
#' @importFrom readODS read_ods
#' @importFrom utils read.table
LoadData <- function(filePath, sheet = 1, naStrings = "", notAllowed = "") {
    stopifnot("Can not import the file. Please give full path to the file you want to import" =
                  file.exists(filePath))
    if (endsWith(filePath, ".xls") || endsWith(filePath, ".xlsx")) {
        raw_data <- read_excel(path = filePath,
                                       col_names = TRUE,
                                       na = naStrings,
                                       sheet = sheet)
    } else if (endsWith(filePath, ".ods")) {
        raw_data <- read_ods(path = filePath,
                                      col_names = TRUE,
                                      na = naStrings,
                                      sheet = sheet)
    } else if (endsWith(filePath, ".tsv") ||
               endsWith(filePath, ".tdf") ||
               endsWith(filePath, ".txt")) {
        raw_data <- read.table(file = filePath,
                               header = TRUE,
                               na.strings = naStrings,
                               sep = "\t",
                               stringsAsFactors = FALSE)
    } else {
        raw_data <- read.table(file = filePath,
                               header = TRUE,
                               na.strings = naStrings,
                               sep = ",",
                               stringsAsFactors = FALSE)
    }
    names(raw_data) <- GenerateValidNames(raw_data)
    notA <- notAllowed %in% unlist(raw_data)
    if (any(notA)) {
        warning("Check file there are still:",
            notAllowed[notA], "in the file\n")
        return(data.frame(raw_data))
    } else {
        cat("Result file succesfully imported\n")
        return(data.frame(raw_data))
    }
}

#' @title Extract basename of string
#'
#' @description Seperates vector of character strings on @param sep
#'     and returns a vector where every entry contains all characters
#'     occuring prior to @sep in the input vector.
#' @param x a vector containing entries with names that where
#'     the names contain the character @param sep that will be used to
#'     split the names. Only characters or number prior to the first
#'     instance of @param sep will be retained.
#' @param sep character to split input name on.
#' @return a vector of names that correspond to leading part of names.
ExtractNames <- function(x, sep = "_") {
    unlist(lapply(strsplit(x = x,
                          split = sep,
                          fixed = TRUE),
                  "[[",
                  1)
           )
}

#' @title Generate valid dataframe with "proper" names
#'
#' @description Will generate a dataframe where column names have been cleaned up.
#' It will remove any white spaces and other illegal characters for names.
#' Will also remove any language specific characters that is not found in the
#' "latin-ascii" character set.
#'
#' @param x a dataframe with names containing the names
#' @return a vector of names that are valid and contain no strange characters.
#' @importFrom stringi stri_trans_general
GenerateValidNames <- function(x) {
    orgNames <- names(x)
    validNames <- gsub("\\.", replacement = "", x = make.names(orgNames))
    validNames <- stri_trans_general(validNames, "latin-ascii")
    return(validNames)
}

#' @title Quality control of  allele sizes from replicate runs.
#'
#' @description Quality control that compare replicate results of
#'  a sample and return 'OK' if all values are identical excluding
#'  missing data.  If values are inconsistent it will return
#'  'Check'. If there is only missing data on a given sample it
#'  will return 'No data'
#' @param data with genotyping results
#' @param locus of the marker to be analysed
QCRep <- function(data = data, locus = "") {
  data <- data[,grepl(locus, names(data))]
  QC <- data[,1] == data[,3] & data[,2] == data[,4]
  QC <- ifelse(QC, yes = "OK", no = "Check")
  QC[is.na(QC)] <- "No Data"
  return(QC)
}

#' @title Merges data from multiple dataframes to single dataframe
#'
#' @description accepts any number of dataframes as arguments and if
#'     row names are the same they will be merged a single dataframe.
#' @param ... list of data ot merged
#' @return datafram with merged data
#' @export
MergePlateS  <- function(...) {
    args <- list(...)
    mergeData <- lapply(args, row.names)
    mergeData <- Reduce(intersect,mergeData)
    if (length(args[[1]][,1]) != length(mergeData)) {
        cat("Can not merge object with different names\n")
    } else {
        cat("Data merged\n")
        return(cbind(...))
  }
}

#' @title Extract the max value from a range of values.
#'
#' @description Returns maximum value from a vector of values and
#'     returns NA if all data is NA.
#' @param x vector with values
#' @return the maximum value from a vector or NA if only missing data is found
MaxMod <- function(x) {
    ifelse(!all(is.na(x)),
           max(x, na.rm = TRUE), NA)
}

#' Extract the min value from a range of values. Returns NA without a
#' warning if all data is NA.
#' @param x vector with values
#' @return the minimum value of a vector or NA if only missing data is found
MinMod <- function(x) {
    ifelse(!all(is.na(x)),
            min(x, na.rm = TRUE), NA)
}

#' @title Creates/Writes excel workbook with pre-specified header read from file
#'  combined with data from selected dataframe.
#'
#' @param headerInfoFile excel file with headers to be part of output
#' @param data dataframe with data to be appended to the headerInfoFile
#' @param sheet which sheet in the excel file to write to
#' @param startRow integer corresponding to row in excel that dataframe will be written to
#' @param outputFile file name to write data to. Should have the file ending .xlsx
#'
#' @return Nothing, but will write and overwrite excel file on disk as a side effect.
#'
#' @importFrom openxlsx loadWorkbook writeData saveWorkbook
#' @export
WriteImportFile <- function(headerInfoFile, data, sheet = 1, startRow = 24, outputFile = "toimportfinal.xlsx") {
    header <- loadWorkbook(file = headerInfoFile)
    writeData(wb = header, sheet = sheet, x = data, colNames = FALSE, startRow = startRow)
    saveWorkbook(header, outputFile, overwrite = TRUE)
    cat(nrow(data) + startRow, "rows was written to", outputFile)
}

#' @title Creates/Writes excel workbook with pre-specified header read from file
#'  combined with data from an excel file with sample data in specific format
#'
#' @param headerInfoFile excel file with headers to be part of output
#' @param importSampleFile full path to excel file with sample data that
#'  will be appended to the headerInfoFile
#' @param sheet which sheet in the excel file to write to
#' @param startRow integer corresponding to row in excel that data will be written to
#' @param outputFile file name to write data to. Should have the file ending .xlsx
#'
#' @return Summary of the data written to file and an excel file is
#'  written to disk as a side effect.
#'
#' @importFrom openxlsx loadWorkbook writeData saveWorkbook
#' @export
WriteIndFile <- function(headerInfoFile, importSampleFile, sheet = 1, startRow = 6, outputFile = "toimportIndfinal.xlsx") {
    header <- loadWorkbook(file = headerInfoFile)
    toind <- read_excel(importSampleFile, skip = 22)
    toind <- as.data.frame(toind)
    names(toind) <- GenerateValidNames(toind)
    toind <- toind[,c("Individnavn", "KjonnID")]
    toind <- toind[!duplicated(toind$Individnavn),]
    toind <- toind[order(toind$Individnavn),]
    writeData(wb = header, sheet = sheet, x = toind, colNames = FALSE, startRow = startRow)
    saveWorkbook(header, outputFile, overwrite = TRUE)
    cat(nrow(toind) + startRow, "rows was written to", outputFile)
}

#' Generate trivial named based on the year that the individual was found and the origin.
#' NB! Only scat samples should be named with this and dead bears should instead have the SVA-ID
#' as trivial name.
#'
#' @return dataframe with two columns. One with individual name and the second containing the trivialname.
#'  the trivial names will look like this from a sample from Norrbotten in 2024: BD24-001.
#'
#' @param importData dataframe with samples to be imported
#' @param individualColumn name of the column with individualdata
#' @param countyColumn name of the column with county data. NB! this should be as the names are given in
#'  rovbase.
#' @param year what year are samples from.
#' @importFrom utils data
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate row_number ungroup
generateTrivialID <- function(importData, individualColumn = "individ", countyColumn = "Fylke", year) {
    # Load the dataset with county shortname
    counties <- SSRqc::counties
    # Ensure counties is available
    if (!exists("counties")) {
        stop("The 'counties' dataset could not be loaded.")
    }
    noDuplicatedInds <- !duplicated(importData[,individualColumn])
    importDataUnique <- importData[noDuplicatedInds,]
    importDataUnique[,countyColumn] <- counties[match(importDataUnique[,"Fylke"], names(counties))]

    importDataUnique <- importDataUnique[order(importDataUnique[,"Fylke"]),]
    importDataUnique <- importDataUnique %>%
        group_by(.data[[countyColumn]]) %>%
        #group_by(Fylke) %>%
        mutate(TrivialID = paste0(.data[[countyColumn]], year, "-", sprintf("%03d", row_number()))) %>%
        ungroup()

    return(importDataUnique)

    }
