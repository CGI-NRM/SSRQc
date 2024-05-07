#' @title Automate tabular data import from different file formats
#'
#' @description Import of tabular data from different file formats
#'     requires use of different packages and syntax. Here we use file
#'     endings to quess the file format and select the most suitable
#'     import tool to generate a dataframe with column
#'     headers. Currently supports excel files (.xlsx, .xls), tab
#'     delimeted files (.tsv, .tdf, .txt), open document format (.ods)
#'     and comma separated files (.csv)
#' @param filePath name of file to import data from (with path)
#' @param sheet which sheet to read data from in excel and ods file
#'     formats
#' @param naStrings entries should be converted to NA values on import
#' @param notAllowed entries that is not valid input and if found will
#'     inhibit import __NB! The function uses file endings to guess
#'     file formats. Will not work if file ending__ __is not
#'     consistent with actual file type__
#' @return a dataframe of the imported data
#' @export
LoadData <- function(filePath, sheet = 1,
                     naStrings,
                     notAllowed) {
    if (endsWith(filePath, ".xls") || endsWith(filePath, ".xlsx")) {
        raw_data <- readxl::read_excel(path = filePath,
                                       col_names = TRUE,
                                       na = naStrings,
                                       sheet = sheet)
    } else if (endsWith(filePath, ".ods")) {
        raw_data <- readODS::read_ods(path = filePath,
                                      col_names = TRUE,
                                      na = naStrings,
                                      sheet = sheet)
    } else if (endsWith(filePath, ".tsv") ||
               endsWith(filePath, ".tdf") ||
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
    notA <- notAllowed %in% unlist(raw_data)
    if (any(notA)) {
        cat("Check export from Geneious there are still:",
            notAllowed[notA], "in the export\n")
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

#' @title Quality control of  allele sizes from replicate runs.
#'
#' @description Quality controls will compare all replicate run of
#'     sample and return 'OK' if all values are identical excluding
#'     missing data.  If values are inconsistent it will return
#'     'Check'. If there is only missing data on a given sample it
#'     will return 'No data'
#' @param data with genotyping results
#' @param locus of the marker to be analysed
QCRep <- function(data = data, locus = "") {
  data <- data[,grepl(locus, names(data))]
  QC <- data[,1] == data[,3] & data[,2] == data[,4]
  QC <- ifelse(QC, yes = "OK", no = "Check")
  QC[is.na(QC)] <- "No Data"
  return(QC)
}

#' @title Merges data from multiple dataframes to dataframe
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
#' @description Returns maximum value from a a vector of values and
#'     returns NA if all data is NA.
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
                min(x, na.rm = T), NA
        )
}

#' @title Merge genotype data with metadata to create match file
#'
#' @description Merges genotype results with metadata creating a
#'     dataframe suitable for matching genotype results with available
#'     genetic database.
#' @param genotypeResults dataframe with genotype results and unique
#'     names of each sample
#' @param filt parameter that turn on and off filtering of data. If
#'     set to No/NO/no the returned data is not filtered if set to
#'     Yes/YES/yes samples with more than the missing data option level
#' @param level the maximum number of markers with missing data
#' @param metaData path and filename to the file with
#'     metadata. Assumes that the file is an excel file and that the file
#'     contains sample name (Strekkode, SEP), DNAid,  date of sampling
#'     and GPS coordinates. 
#' @param confirmedDead is the sample confirmed dead or not. Allowed values
#'     Yes/YES/yes or No/NO/no
#' @param colnamesWithSEP name of the column in metaData that contains
#'     SEP number eg. SEP##### or M####
#' @return dataframe with genetic and metadata merged
#' @export
CreateMatchInput <- function(GenotypeResults, filt = "No",
                             level = 4,
                             metaData = "rovbasemetadata.xlsx",
                             confirmedDead = "No",
                             colnamesWithSEP = "Strekkode pröve") {
    cDead <- toupper(confirmedDead)
    stopifnot("confirmedDead can only be yes or no" =
                  cDead %in% c("YES","NO"))
    filt <- toupper(filt)
    stopifnot("filt can only be yes or no"= filt %in% c("YES","NO"))
    stopifnot("Can not import the metadatafile. Please give full path to file"=
              file.exists(metaData))
    metadata <- suppressWarnings(read_excel(metaData))
    matchData <- merge(GenotypeResults, metadata, by.x = "row.names", by.y = colnamesWithSEP)
    matchData$confirmed_dead  <- confirmedDead
    matchExport <- data.frame(index =  matchData$Row.names,
                              DNAid = matchData$`DNAID (Prøve)`,
                              Individ = matchData$Individ,
                              RbSex   = matchData$`Kjønn (Analyse)`,
                              date  = matchData$Funnetdato,
                              north = matchData$`Nord (UTM33/SWEREF99 TM)`,
                              east  = matchData$`Øst (UTM33/SWEREF99 TM)`,
                              gender = matchData$consensusSex,
                              confirmed_dead = matchData$confirmed_dead,
                              G1A_1  = matchData$G1A.1min,
                              G1A_2  = matchData$G1A.2max,
                              G1D_1  = matchData$G1D.1min,
                              G1D_2  = matchData$G1D.2max,
                              G10B_1  = matchData$G10B.1min,
                              G10B_2  = matchData$G10B.2max,
                              G10L_1  =  matchData$G10L.1min,
                              G10L_2  = matchData$G10L.2max,
                              MU05_1  = matchData$MU05.1min,
                              MU05_2  = matchData$MU05.2max,
                              MU09_1  = matchData$MU09.1min,
                              MU09_2  = matchData$MU09.2max,
                              MU10_1  = matchData$MU10.1min,
                              MU10_2  = matchData$MU10.2max,
                              MU15_1  = matchData$MU15.1min,
                              MU15_2  = matchData$MU15.2max,
                              MU23_1  = matchData$MU23.1min,
                              MU23_2  = matchData$MU23.2max,
                              MU50_1  = matchData$MU50.1min,
                              MU50_2  = matchData$MU50.2max,
                              MU51_1  = matchData$MU51.1min,
                              MU51_2  = matchData$MU51.2max,
                              MU59_1  = matchData$MU59.1min,
                              MU59_2  = matchData$MU59.2max)
    matchExport$gender[is.na(matchExport$gender)] <- "Okänt"
    NumberOfSamplesWithGenotype <- nrow(GenotypeResults)
    NumberOfSamplesInRovbaseMeta <- nrow(matchExport)
    if(filt == "NO") {
        cat(NumberOfSamplesInRovbaseMeta, "out of",
            NumberOfSamplesWithGenotype,
            "samples have information in rovbase", "\n")
        cat("Saved:", nrow(matchExport), "sample(s)", "\n")
        return(matchExport)
    } else {
        matchExport <- matchExport[rowSums(is.na(matchExport[,c("G10L_1",
                                                                    "G10L_2",
                                                                    "MU05_1",
                                                                    "MU05_2",
                                                                    "MU09_1",
                                                                    "MU09_2",
                                                                    "MU10_1",
                                                                    "MU10_2",
                                                                    "MU23_1",
                                                                    "MU23_2",
                                                                    "MU50_1",
                                                                    "MU50_2",
                                                                    "MU51_1",
                                                                    "MU51_2",
                                                                    "MU59_1",
                                                                    "MU59_2")]
                                                 ))<5,]
        cat(NumberOfSamplesInRovbaseMeta, "out of",
            NumberOfSamplesWithGenotype,
            "samples have information in rovbase", "\n")
        cat("Saved:", nrow(matchExport), "sample(s)", "\n")
        return(matchExport)
    }
}

#' Function to from the match database create import file according
#' the standard given by rovbase
#'
#' @return a dataframe with the neccessary content for import to
#' rovbase.
#'
#' @param fromMatchDB Dataframe with target data
#'
#'
#' @param RBData Dataframe with metadata from rovbase for the samples
#' in fromMatchDB

CreateImport <- function(fromMatchDB, RBData) {
        miss <- fromMatchDB$index[!fromMatchDB$index %in% RBData$Strekkode]
        completeFile <- merge(fromMatchDB, RBData,
                by.x = "index.x", by.y = "Strekkode"
        )
        nr <- nrow(completeFile)
        res <- data.frame(
                DNAID = completeFile$DNAID,
                Strekkode = completeFile$index.x,
                RovbaseID = completeFile$RBID.x,
                ArtsID = rep(3, nr),
                Prøvestatus = rep(1, nr),
                KjønnID = completeFile$gender.x,
                MetodeID = rep(18, nr),
                VurderingID = rep(1, nr),
                Kontrollstatus = rep(6, nr),
                AnalysertAvID = rep(10, nr),
                Merknad = rep("", nr),
                Individ = completeFile$individ,
                rbIND   = completeFile$Individ.x.x,
                Lan     = completeFile$Fylke.x.x,
                TrivialID = completeFile$TrivialID
        )
        sum <- list(
                miss = miss,
                res = res
        )
        return(sum)
}

WriteImportFile <- function(data = results, file = "toimportfinal.xlsx") {
    openxlsx::write.xlsx(x = data, file = file, overwrite = FALSE)
    }


createIndFile <- function(fromMatchDB) {
    tempData  <- data.frame(IndividNavn = c(NA, NA, NA,
                                            "IndividNavn"),
                            Kjønn = c("1 = Hane",
                                      "2 = Hona",
                                      "3 = Okänt",
                                      "Kjønn"))
    nyInds <- fromMatchDB[startsWith(fromMatchDB$Individ,"NRM_"),
                                        c("Individ", "Gender", "SVAID")]
                    names(nyInds)  <- c("IndividNavn", "Kjønn")
                    res <- rbind(tempData, nyInds)
                    return(res)
}

