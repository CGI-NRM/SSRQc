#' Helper function read facilitate data import of different file
#' formats
#'
#' @return a dataframe of the imported data
#' @param file_path name of file to import data from (with path)
#' @param sheet which sheet to read data from in excel and ods file formats
#' @param naStrings entries should be converted to NA values on import
#' __NB! The function uses file endings to guess file formats. Will not work if file ending__
#' __is not consistent with actual file type__
#' @export
LoadData <- function(filePath, sheet = 1, naStrings = c("NA", "-99", "0", "000", "No peaks in locus", "No peaks")) {
    if (endsWith(filePath, ".xls") | endsWith(filePath, ".xlsx")) {
        raw_data <- readxl::read_excel(path = filePath,
                                       col_names = TRUE,
                                       na = naStrings,
                                       sheet = sheet)
    } else if (endsWith(filePath, ".ods")) {
        raw_data <- readODS::read_ods(path = filePath,
                                      col_names = TRUE,
                                      na = naStrings,
                                      sheet = sheet)
    } else if (endsWith(filePath, ".tsv") |
               endsWith(filePath, ".tdf") |
               endsWith(filePath, ".txt")) {
        raw_data <- utils::read.table(file = filePath,
                               header = TRUE,
                               na.strings = naStrings,
                               sep = "\t",
                               stringsAsFactors = FALSE)
    } else {
        raw_data <- utils::read.table(file = filePath,
                               header = TRUE,
                               na.strings = naStrings,
                               sep = ",",
                               stringsAsFactors = FALSE)
    }
    raw_datanona <- stats::na.omit(raw_data)
    temp <- any(raw_datanona == "Unbinned peaks in locus" |
                raw_datanona == "Too many alleles")
    if (temp) {
        cat("Check Geneious file there are still 'Unbinned peaks' and/or 'Too
            many alleles' at some markers\n")
        return(NULL)
    } else {
    cat("Result file succesfully imported\n")    
    return(data.frame(raw_data))
    }
}

#' Function that extract individual names from index (SEP) or M numbers
#'
#' @return a vector of names that 10 long if SEP or 7 long in all
#' other cases
#' @param nameVector a vector containing entries with names that have
#' replicate name and SSR mix attached. The latter will be stripped in
#' and only the _actual_ names of samples will be retained.
#' 
ExtractNames <- function(nameVector) {
    ifelse(startsWith(nameVector, "SEP"),
           yes = substring(nameVector, 1, 10),
           no = substring(nameVector, 1, 7))
}

#' Check if allele sizes from replicate runs are the same. Returns 'OK'
#' if all looks identical and returns 'Check' if replicates are
#' inconsistent. Returns "No data" if there is no data available.
#' @param data with genotyping results
#' @param locus of the marker to analysed
QCRep <- function(data = data, locus = "") {
  data <- data[,grepl(locus, names(data))]
  QC <- data[,1] == data[,3] & data[,2] == data[,4]
  QC <- ifelse(QC, yes = "OK", no = "Check")
  QC[is.na(QC)] <- "No Data"
  QC
}


#' Add genotype results from SSR markers that are analysed on
#' different PCR mixes to a single dataframe
#' @param data1 First dataframe with results
#' @param data2 Second dataframe with results
#' @export
MergePlate <- function(data1, data2) {
    if (!identical(row.names(data1$Genotypes), row.names(data2$Genotypes))) {
        cat("Can not merge objects with different names\n")
    } else {
        cat("Data merged\n")
        return(cbind(data1, data2))
    }
}

#' Add genotype results different marker types to 
#' a single dataframe
#' @param data1 First dataframe with results
#' @param data2 Second dataframe with results
#' @export
MergePlateSum <- function(data1, data2) {
  if (!identical(row.names(data1), row.names(data2))) {
    cat("Can not merge objects with different names\n")
  } else {
    cat("Data merged\n")
    return(cbind(data1, data2))
  }
}

#' Extract the max value from a range of values. Returns NA without a
#' warning if all data is NA.
#' @param x vector with values
#' @return the maximum value from a vector or NA if only missing data
MaxMod <- function(x) {
    ifelse(!all(is.na(x)),
           max(x, na.rm=T), NA)
}

#' Extract the min value from a range of values. Returns NA without a
#' warning if all data is NA.
#' @param x vector with values
#' @return the minimum value of a vector or NA if only missing data
MinMod <- function(x) {
    ifelse(!all(is.na(x)),
           min(x, na.rm = T), NA)
}

#' Merge genotype data with metadata
#'
#' 

CreateMatchInput <- function(GenotypeResults, filt = "No") {
  header <- c("index", "date", "north", "east",
              "gender", "confirmed_dead", "G10L_1",
              "G10L_2", "MU05_1", "MU05_2", "MU09_1",
              "MU09_2", "MU10_1", "MU10_2", "MU23_1",
              "MU23_2", "MU50_1", "MU50_2", "MU51_1",
              "MU51_2", "MU59_1", "MU59_2", "Sex_1", "Sex_2")
  # metadata <- read_excel("~/ownCloud/Spillning2020/QCResults/rovbasedata.xlsx")
  metadata <- read_excel("~/ownCloud/Spillning2021/rovbasemeta.xlsx")
  metadata <- metadata[,c(1,4:6)]
  names(metadata) <- c("index", "north", "east", "date")
  tt <- merge(GenotypeResults[,c("G10L.1min", "G10L.2max",
                                 "MU05.1min", "MU05.2max",
                                 "MU09.1min", "MU09.2max",
                                 "MU10.1min", "MU10.2max",
                                 "MU23.1min", "MU23.2max",
                                 "MU50.1min", "MU50.2max",
                                 "MU51.1min", "MU51.2max",
                                 "MU59.1min", "MU59.2max")], 
              metadata, 
              by.x = "row.names", 
              by.y = "index")
  tt$confirmed_dead <- "No"
  tt <- tt[,c(1,20,18,19,21,2:17)]
  names(tt) <- header
  tt$gender[is.na(tt$gender)] <- "Okänt"
  if(filt == "No") {
    tt
  } else {
    tt[rowSums(is.na(tt[7:22]))<5,]
  }
}

CreateMatchInputDead <- function(GenotypeResults, filt = "No") {
  header <- c("index", "date", "north", "east",
              "gender", "confirmed_dead", "G10L_1",
              "G10L_2", "MU05_1", "MU05_2", "MU09_1",
              "MU09_2", "MU10_1", "MU10_2", "MU23_1",
              "MU23_2", "MU50_1", "MU50_2", "MU51_1",
              "MU51_2", "MU59_1", "MU59_2")
  metadata <- read_excel("/home/thomkall/Develop/SSRqc/inst/extdata/DeadbearsRovbase.xlsx")
  metadata <- metadata[,c(1,3:5)]
  names(metadata) <- c("index", "north", "east", "date")
  tt <- merge(GenotypeResults, metadata, by.x = "row.names", by.y = "index")
  tt$confirmed_dead <- "Yes"
  tt <- tt[,c("Row.names", "date", "north", "east", "consensus", "confirmed_dead", "G10L.1min",
        "G10L.2max", "MU05.1min", "MU05.2max", "MU09.1min", "MU09.2max",
        "MU10.1min", "MU10.2max", "MU23.1min", "MU23.2max", "MU50.1min",
        "MU50.2max", "MU51.1min", "MU51.2max", "MU59.1min", "MU59.2max")]
  names(tt) <- header
  tt$gender[is.na(tt$gender)] <- "Okänt"
  if(filt == "No") {
    tt
  } else {
    tt[rowSums(is.na(tt[7:22]))<5,]
  }
}

# CreateImportData <- function(fromMatchDB, RBData, fileImport, fileInd) {
#    cols <- c("DNAID","Strekkode","RovbaseID",
#              "ArtsID", "Provestatus", "KjonnID", "MetodeID",
#              "VurderingID", "Kontrollstatus", "AnalysertAvID",
#              "Merknad", "IndividID", "Individnavn", "ReferenseID",
#              "Referensetekst")
    
    
    
