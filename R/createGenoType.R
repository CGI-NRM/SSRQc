#' @title
#' Extract most likely genotype data from replicates
#'
#' @description
#' Function to extract the most likely genotype scores from multiple
#' replicates of genotype data. The structure of the input data is a
#' table (excel etc) where the first column contains sample names and
#' the genotypes of interest can be in any other column in the input
#' data. Note that the function assumes that replicates will have the
#' same name after removing the tail part of the original name. Eg. if
#' replicates are named SEP1000001_mix1_rep1, SEP1000001_mix1_rep2 etc
#' the sample name will be SEP1000001.
#' NB! The naming used for loci in this implementation assumes that
#' there are two entries for every locus and the only allowed special
#' character is "_". The two alleles should have the same name and be followed
#' by a 1 for allele one and a 2 for allele two. See input examples in
#' vignette.
#'
#' @return a list of three dataframes: Genotypes, QC and number of
#' missing data points
#' @param file name of the file with replicate data
#' @param locus names of the loci to be analysed
#' @param naStrings string treated as NA in import
#' @param notAllowed strings that will stop import
#' @export
#'
CreateGenotype <- function(file = file,
                           locus = c("Locus1", "Locus2"),
                           naStrings = c("NA", "-99", "0", "000",
                           "No peaks in locus", "No peaks"),
                           notAllowed = c("Unbinned peaks",
                                          "Too many alleles",
                                          "Unbinned peaks in locus",
                                          "Peaks outside loci")) {
    dataPlate <- LoadData(file, naStrings = naStrings, notAllowed = notAllowed)
    if(!is.null(dataPlate) & any(grepl(paste0(locus, collapse = "|"), names(dataPlate)))) {
        dataPlate$names <- ExtractNames(dataPlate[,1])
        togrep <- paste(locus, collapse="|")
        dataPlate1 <- dataPlate[,grepl(togrep, names(dataPlate))]
        #colnames(dataPlate1) <- GenerateValidNames(dataPlate1)
        sampleDataMax <- stats::aggregate(dataPlate1, by = list(dataPlate$names),
                                MaxMod)
        row.names(sampleDataMax) <- sampleDataMax$Group.1
        sampleDataMax <- sampleDataMax[, - which(names(sampleDataMax) == "Group.1")]
        colnames(sampleDataMax) <- paste0(colnames(sampleDataMax), ".max")
        sampleDataMin <- stats::aggregate(dataPlate1, by = list(dataPlate$names),
                                          MinMod)
        row.names(sampleDataMin) <- sampleDataMin$Group.1
        sampleDataMin <- sampleDataMin[, - which(names(sampleDataMin) == "Group.1")]
        colnames(sampleDataMin) <- paste0(colnames(sampleDataMin), ".min")
        sampleData <- cbind(sampleDataMin, sampleDataMax)
        QC <- lapply(locus, QCRep, data = sampleData)
        QC <- as.data.frame(do.call(cbind, QC))
        colnames(QC) <- paste(locus, ".QC", sep = "")
        sampleData <- sampleData[,grepl("1.min|2.max", names(sampleData))]
        QCGenotypes <- cbind(sampleData, QC)
        QCGenotypes <- QCGenotypes[,order(names(QCGenotypes))]
        genotypes <- QCGenotypes[,!grepl("QC", names(QCGenotypes))]
        result <- list(QC = QCGenotypes,
                    Genotypes = genotypes,
                    genotype_missing = rowSums(is.na(genotypes)))
        return(result)
    } else {
        cat("There is no data from specified loci in the given file\n")
        return(NULL)
    }
}

#' Determine sex from genotype data
#'
#' This function will based on the selected number of loci determine the
#' most likely sex of a sample from genotype data.
#' NB! The function assumes that the marker with name Male have the following
#' properties:
#' 1. 152 152 == Female
#' 2.  98 152 == Male
#' 3.  98  98 == Male
#' 4.  NA 152 == Female
#' 5.  NA  98 == Male
#' 6.  NA  NA == NA
#'
#' The other markers that are y-chromosome markers have the following properties:
#' 1. NA NA == Female
#' 2. Any other value == Male
#'
#' The third marker is found in both male and females and can not
#' separate between sexes, but is used as a control that the sample can
#' be analysed for sex chromosomes. So strictly speaking a working sample irrespective of
#' sex should have data on this marker.
#'
#' Consensus for this is the extracted as follows:
#' 1. All marker the same sex -> this is the true sex of the sample
#' 2. All markers have NA -> sex is not determined NA is reported
#' 3. In case of mismatch ->  most common sex i extracted.
#' 4. The Male marker is better as the y-chromosome markers can not separate failed from female.
#'
#' @param data a dataframe with genotype results including loci that
#'     are linked to sex chromosomes
#' @param locus a vector of loci that contains genotypes for sex determination.
#' @return list with input data plus QC score and a vector with most
#' likely sex
#' @importFrom stats na.omit
#' @export
SexDetermination <- function(data = data, locus = c("Male",
                                                    "UaY15020",
                                                    "UarY369.4",
                                                    "Y318.2",
                                                    "SMCY",
                                                    "ZFX")) {
    sexMarkerMale  <- rep(NA, nrow(data))
    sexMarkerUaY15020  <- rep(NA, nrow(data))
    sexMarkerUarY369.4  <- rep(NA, nrow(data))
    sexMarkerY318.2    <- rep(NA, nrow(data))
    sexMarkerSMCY <- rep(NA, nrow(data))
    sexMarkerZFX <- rep(NA, nrow(data))

    if(any(grepl("Male", colnames(data)))) {
        resSexMarkerMale <- data[,grepl("Male", names(data))]
        resSexMarkerMale$merge <- paste(resSexMarkerMale[,1], resSexMarkerMale[,2])
        sexMarkerMale <- ifelse(resSexMarkerMale$merge == "NA 152"  |
                                resSexMarkerMale$merge == "152 152" |
                                resSexMarkerMale$merge == "152 NA",
                            yes = "Hona",
                            no = ifelse(resSexMarkerMale$merge == "98 152" |
                                            resSexMarkerMale$merge == "98 NA"  |
                                            resSexMarkerMale$merge == "NA 98"  |
                                            resSexMarkerMale$merge == "152 98",
                                        yes = "Hane",
                                        no = NA))

    }
    if(any(grepl("UaY15020", colnames(data)))) {
        resSexMarkerUaY15020    <- data[,grepl("UaY", names(data))]
        sexMarkerUaY15020 <- ifelse(is.na(resSexMarkerUaY15020[,1]),
                            yes = "Hona",
                            no  = "Hane")
    }

    if(any(grepl("UarY369.4", colnames(data)))) {
        resSexMarkerUarY369.4 <- data[,grepl("UarY369", names(data))]
        sexMarkerUarY369.4 <- ifelse(is.na(resSexMarkerUarY369.4[,1]),
                            yes = "Hona",
                            no  = "Hane")
    }
    if(any(grepl("ZFX", colnames(data)))) {
        resSexMarkerZFX <- data[,grepl("ZFX", names(data))]
        sexMarkerZFX <- ifelse(is.na(resSexMarkerZFX[,1]),
                            yes = NA,
                            no  = "Hona/Hane")
    }
    if(any(grepl("Y318.2", colnames(data)))) {
        resSexMarkerY318.2 <- data[, grepl("Y318", names(data))]
        sexMarkerY318.2 <- ifelse(is.na(resSexMarkerY318.2[,1]),
                                yes = "Hona",
                                no  = "Hane")
    }
    if(any(grepl("SMCY", colnames(data)))) {
        resSexMarkerSMCY    <- data[,grepl("SMCY", names(data))]
        sexMarkerSMCY <- ifelse(is.na(resSexMarkerSMCY[,1]),
                            yes = "Hona",
                            no  = "Hane")
    }
    geneticSex <- data.frame(sexMarkerZFX, sexMarkerMale,
        sexMarkerUarY369.4, sexMarkerUaY15020, sexMarkerY318.2,
        sexMarkerSMCY)
    SexCon <- function(sexVector) {
        sexVno2  <- factor(sexVector[-1])
        if(is.na(sexVector[1])) {
            NA
        } else if (names(sort(table(sexVno2),
                              decreasing = TRUE)[1]) == "Hane") {
            "Hane"
        } else {
            "Hona"
        }
    }
    SexQC  <- function(sexVector) {
        sexV2  <- sexVector[-1]
        ifelse(length(unique(na.omit(sexV2))) == 1,
               yes = "OK",
               no  = "Check!")
    }
    QC <- unlist(apply(geneticSex, MARGIN = 1, SexQC))
    consensusSex <- unlist(apply(geneticSex, MARGIN = 1, SexCon))
    geneticSex <- cbind(geneticSex, QC, consensusSex)
    colnames(geneticSex) <- c(names(geneticSex[-ncol(geneticSex)]), "Sex")
    rownames(geneticSex) <- rownames(data)
    consensusSex <- as.data.frame(consensusSex)
    rownames(consensusSex) <- rownames(data)
    return(list(QC = geneticSex,
                Genotypes = consensusSex))
}
