#' Helper function read facilitate data import of different file
#' formats
#'
#' @return a dataframe of the imported data
#' @param file_path name of file to import data from (with path)
#' @param sheet which sheet to read data from in excel and ods file formats
#' @param na_strings entries should be converted to NA values on import
#' __NB! The function uses file endings to guess file formats. Will not work if file ending__
#' __is not consistent with actual file type__
#' @export
load_data <- function(file_path, sheet = 1, na_strings = c("NA", "-99", "0", "000", "No peaks in locus", "No peaks")) {
    if (endsWith(file_path, ".xls") | endsWith(file_path, ".xlsx")) {
        raw_data <- readxl::read_excel(path = file_path,
                                       col_names = TRUE,
                                       na = na_strings,
                                       sheet = sheet)
    } else if (endsWith(file_path, ".ods")) {
        raw_data <- readODS::read_ods(path = file_path,
                                      col_names = TRUE,
                                      na = na_strings,
                                      sheet = sheet)
    } else if (endsWith(file_path, ".tsv") |
               endsWith(file_path, ".tdf") |
               endsWith(file_path, ".txt")) {
        raw_data <- read.table(file = file_path,
                               header = TRUE,
                               na.strings = na_strings,
                               sep = "\t",
                               stringsAsFactors = FALSE)
    } else {
        raw_data <- read.table(file = file_path,
                               header = TRUE,
                               na.strings = na_strings,
                               sep = ",",
                               stringsAsFactors = FALSE)
    }
    raw_datanona <- na.omit(raw_data)
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
#' @param name_vector a vector containing entries with names that have
#' replicate name and SSR mix attached. The latter will be stripped in
#' and only the _actual_ names of samples will be retained.
#' 
extract_names <- function(name_vector) {
    ifelse(startsWith(name_vector, "SEP"),
           yes = substring(name_vector, 1, 10),
           no = substring(name_vector, 1, 7))
}

#' Check if allele sizes from replicate runs are the same. Returns 'OK'
#' if all looks identical and returns 'Check' if replicates are
#' inconsistent. Returns "No data" if there is no data available.
#' @param dataframe with genotyping results
#' @locus Name of the marker to analysed
qc_rep <- function(data, locus = "") {
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
mergePlate <- function(GenoMix1, GenoMix2) {
    if (!identical(row.names(GenoMix1$Genotypes), row.names(GenoMix2$Genotypes))) {
        cat("Can not merge objects with different names\n")
    } else {
        cat("Data merged\n")
        return(cbind(GenoMix1$Genotypes, GenoMix2$Genotypes))
    }
}

#' Extract the max value from a range of values. Returns NA if all
#' data is NA without generating a warning.
#' @param x vector with values
max_mod <- function(x) {
    ifelse(!all(is.na(x)),
           max(x, na.rm=T), NA)
}

#' Extract the min value from a range of values. Returns NA if all
#' data is NA without generating a warning.
#' @param x vector with values
min_mod <- function(x) {
    ifelse(!all(is.na(x))
          ,
           min(x, na.rm = T), NA)
}


