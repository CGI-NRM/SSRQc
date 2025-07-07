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
#'     metadata. Assumes that the file can be imported using LoadData
#'     function and that the file contains sample name (Strekkode, SEP),
#'     DNAid,  date of sampling and GPS coordinates.
#' @param confirmedDead is the sample confirmed dead or not. Allowed values
#'     Yes/YES/yes or No/NO/no
#' @param colnamesWithSEP name of the column in metaData that contains
#'     SEP number eg. SEP##### or M####
#' @return dataframe with genetic and metadata merged
#' @importFrom lubridate dmy
#'
#' @export
#'
CreateMatchInput <- function(genotypeResults, filt = "No",
                             level = 4,
                             metaData = "rovbasemetadata.xlsx",
                             confirmedDead = "No",
                             colnamesWithSEP = "StrekkodeProve") {
    cDead <- toupper(confirmedDead)
    stopifnot("confirmedDead can only be yes or no" =
                  cDead %in% c("YES","NO"))
    filt <- toupper(filt)
    stopifnot("filt can only be yes or no" = filt %in% c("YES","NO"))
    stopifnot("Can not import the metadatafile. Please give full path to file"=
                  file.exists(metaData))
    metadata <- LoadData(metaData)
    matchData <- merge(genotypeResults, metadata, by.x = "row.names", by.y = colnamesWithSEP)
    matchData$confirmed_dead  <- confirmedDead
    matchExport <- data.frame(index =  matchData$Row.names,
                              DNAid = matchData[,"DNAIDProve"],
                              Individ = matchData$Individ,
                              RbSex   = matchData[,"KjonnAnalyse"],
                              date  = matchData$Funnetdato,
                              north = matchData[,"NordUTM33SWEREF99TM"],
                              east  = matchData[,"OstUTM33SWEREF99TM"],
                              gender = matchData$consensusSex,
                              confirmed_dead = matchData$confirmed_dead,
                              G1A_1  = matchData$G1A1.min,
                              G1A_2  = matchData$G1A2.max,
                              G1D_1  = matchData$G1D1.min,
                              G1D_2  = matchData$G1D2.max,
                              G10B_1  = matchData$G10B1.min,
                              G10B_2  = matchData$G10B2.max,
                              G10L_1  =  matchData$G10L1.min,
                              G10L_2  = matchData$G10L2.max,
                              MU05_1  = matchData$MU051.min,
                              MU05_2  = matchData$MU052.max,
                              MU09_1  = matchData$MU091.min,
                              MU09_2  = matchData$MU092.max,
                              MU10_1  = matchData$MU101.min,
                              MU10_2  = matchData$MU102.max,
                              MU15_1  = matchData$MU151.min,
                              MU15_2  = matchData$MU152.max,
                              MU23_1  = matchData$MU231.min,
                              MU23_2  = matchData$MU232.max,
                              MU50_1  = matchData$MU501.min,
                              MU50_2  = matchData$MU502.max,
                              MU51_1  = matchData$MU511.min,
                              MU51_2  = matchData$MU512.max,
                              MU59_1  = matchData$MU591.min,
                              MU59_2  = matchData$MU592.max)
    matchExport$gender[is.na(matchExport$gender)] <- "Ok\u00e4nt"
    NumberOfSamplesWithGenotype <- nrow(genotypeResults)
    NumberOfSamplesInRovbaseMeta <- nrow(matchExport)
    if(filt == "NO") {
        cat(NumberOfSamplesInRovbaseMeta, "out of",
            NumberOfSamplesWithGenotype,
            "samples have information in rovbase", "\n")
        cat("Saved:", nrow(matchExport), "sample(s) in dataframe", "\n")
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
        ))<level + 1,]
        cat(NumberOfSamplesInRovbaseMeta, "out of",
            NumberOfSamplesWithGenotype,
            "samples have information in rovbase", "\n")
        cat("Saved:", nrow(matchExport), "sample(s) in dataframe", "\n")
        return(matchExport)
    }
}

#' @title Function to from the match database and rovbasefiles create an import file according
#' the standard given by rovbase
#'
#' @return a dataframe with the neccessary content for import to
#' rovbase.
#'
#' @param matchDBFile file (.xlsx, .csv, tsv) with data to be imported. This are
#'  samples matched to local genetic database and that are to be imported to
#'  rovbase. The file should have the following columns (but can contain others):
#'  - index
#'  - date
#'  - north
#'  - east
#'  - gender
#'  - date_changed
#'  - confirmed_dead
#'  - individ
#'
#' @param metaDataFile file with metadata from rovbase. This data will be merged
#'  with samples informtation from the matchDBFile found in the fromMatchDB dataframe.
#'  __NB! samples found in the dataframe matchDBFile that is not found in rovbase
#'  will be saved in output of this function.__
#'  This dataframe needs to contain the following columns (but can contain others):
#'  - Strekkodeprove - SEP or M number. Original name in rb = "Strekkode (Prøve)"
#'  - DNAIDProve - DNA ID. Original name in rb = "DNAID (Prøve)"
#'  - RovbaseIDProver - Rovbase ID. Original name in rb = "RovbaseID (Prøve)"
#'  - Individ - Rovbase individ from rb.
#'  - KjonnIndivid - Rovbase sex. Original name in rb = "Kjønn (Individ)"
#'  - Provetype - Sample type. Original name in rb = "Prøvetype"
#'
#' @param year which years were the samples collected in NB! accepts the last two
#'  digits so 2024 should be 24.
#'
#' @param onlyNew Since new individuals also need to have a trivial id created this
#'  new inds should be run with this set to TRUE and samples that matched already
#'  known individuals will be handled correctly if this is set to FALSE.
#'
#' @param sampleType string setting the samples types that will be kept in the
#'  output. NB! We mostly do scat samples and this should then be left as "Spillning".
#'
#' @param methodID Number corresponding to the methods used for analyzing the samples.
#'  8 SSR markers that were originally used at NRM for sample identification goes under
#'  the methodID 18.
#'
#' @param finalVersion Boolean value if set to TRUE only the columns needed for import
#'  file to rovbase are saved. If FALSE a version with columns that can be manually checked
#'  are retained. NB! Use FALSE and manually clean the file in excel unless you are 100% sure
#'  that the data is consistent in the matching data and rovbase.
#'
#' @return List with two dataframes. The one named miss will contain the samples
#'  found in the matching database, but than are not found in rovbase. The second one
#'  contains as dataframe with data almost ready for import to rovbase.
#'  NB! The sex returned as "rbKjonnID" is the sex of the __individual__ in rovbase
#'  and hence not actual extracted from the specific sample eg. matching with
#'  "StrekkodeProve".
#'
#' @export

 # matchDBFile <- "~/Develop/SSRQc/inst/extdata/matchDB.csv"
 # metaDataFile <- "~/Develop/SSRQc/inst/extdata/rovbasemetadata.xlsx"
 # year = 24
 # onlyNew = FALSE
 # sampleType = "Spillning"
 # methodID = 18
 # finalVersion = FALSE
CreateSampleImport <- function(matchDBFile = "matchdb.csv",
                               metaDataFile = "rovbasemetadata.xlsx",
                               year = 24, onlyNew = TRUE, sampleType = "Spillning",
                               methodID = 18, finalVersion = FALSE) {
    stopifnot("onlyNew can only be a boolean value (TRUE or FALSE)" =
                  is.logical(onlyNew))
    stopifnot("finalVersion can only be a boolean value (TRUE or FALSE)" =
                  is.logical(onlyNew))
    stopifnot("Can not import the matchdatabase file. Please give full path to file"=
                  file.exists(matchDBFile))
    stopifnot("Can not import the metadata file. Please give full path to file"=
                  file.exists(metaDataFile))
    fullyear <- as.numeric(paste0(20, year))
    fromMatchDB <- LoadData(matchDBFile)
    metaData <- LoadData(metaDataFile)
    fromMatchDB$date <- dmy(fromMatchDB$date)
    fromMatchDB <- fromMatchDB[year(fromMatchDB$date) == fullyear,]
    if(onlyNew) {
        fromMatchDB <- fromMatchDB[grepl("NRM", fromMatchDB$individ),]
    } else {
        fromMatchDB <- fromMatchDB[grepl("^BI", fromMatchDB$individ),]
    }
    metaData$IndividShort <- unlist(lapply(strsplit(metaData$Individ, " "), `[`, 1))
    miss <- fromMatchDB$index[!fromMatchDB$index %in% metaData$StrekkodeProve]
    completeFile <- merge(fromMatchDB, metaData,
                          by.x = "index", by.y = "StrekkodeProve"
    )
    completeFile <- as.data.frame(completeFile)
    completeFile <- completeFile[completeFile[,"Provetype"] == sampleType,]
    completeFile <- completeFile[order(completeFile$date),]
    nr <- nrow(completeFile)
    if(onlyNew) {
        TrivialID <- generateTrivialID(importData = completeFile, individualColumn = "individ", countyColumn = "Fylke", year = year)
        TrivialID <- TrivialID[, c("individ", "TrivialID")]
        completeFile <- merge(completeFile, TrivialID, by = "individ")
    } else {
        TrivialID <- rep("",nr)
        completeFile <- cbind(completeFile, TrivialID)
        completeFile[,"KjonnIndivid"] <- metaData[,"KjonnIndivid"][match(completeFile$individ, metaData$IndividShort)]
    }
    if(finalVersion) {
        sex <- ifelse(completeFile$gender == "Hane",
                      yes = 1, no = ifelse(completeFile$gender == "Hona",
                                           yes = 2, no = 3))
        res <- data.frame(
            DNAID = completeFile[,"DNAIDProve"],
            Strekkode = completeFile$index,
            RovbaseID = completeFile[,"RovbaseIDProve"],
            ArtsID = rep(3, nr),
            ProvestatusID = rep(1, nr),
            KjonnID = sex,
            MetodeID = rep(methodID, nr),
            VurderingID = rep(1, nr),
            Kontrollstatus = rep(6, nr),
            AnalysertAvID = rep(10, nr),
            Merknad = rep("", nr),
            IndividID = completeFile[,"individ"],
            IndividNavn = completeFile[,"TrivialID"],
            ReferenseID = rep(8, nr))
    } else {
        res <- data.frame(
            DNAID = completeFile[,"DNAIDProve"],
            Strekkode = completeFile$index,
            RovbaseID = completeFile[,"RovbaseIDProve"],
            ArtsID = rep(3, nr),
            Provestatus = rep(1, nr),
            KjannID = completeFile[,"gender"],
            MetodeID = rep(methodID, nr),
            VurderingID = rep(1, nr),
            Kontrollstatus = rep(6, nr),
            AnalysertAvID = rep(10, nr),
            Merknad = rep("", nr),
            IndividID = completeFile[,"individ"],
            IndividNavn = completeFile[,"TrivialID"],
            ReferanseID = rep(8, nr),
            ReferenseTekst = rep("", nr),
            NRMIndivid = completeFile[,"individ"],
            rbIndivid   = completeFile[,"Individ"],
            rbKjonnID = completeFile[,"KjonnIndivid"],
            Lan     = completeFile[,"Fylke"])
    }
    sum <- list(
        miss = miss,
        res = res
    )
    return(sum)
}
