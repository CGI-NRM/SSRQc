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
#' there are spaces in the names. Hence the need to modify the names
#' with the replacement of "..." with ".".
#'
#' 
#' @return a list of three dataframes: Genotypes, QC and number of
#' missing data points
#' @param file name of the file with replicate data
#' @param locus names of the loci to be analysed
#' @export

createGT <- function(file = file, locus = c("Locus1", "Locus2")) {
    platta <- LoadData(file)
    if(!is.null(platta) {
        platta$names <- ExtractNames(platta[,1])
        togrep <- paste(locus, collapse="|")
        platta1 <- platta[,grepl(togrep, names(platta))]
        
        indDataMax <- stats::aggregate(platta1, by = list(platta$names),
                                MaxMod)
        row.names(indDataMax) <- indDataMax$Group.1
        indDataMax <- indDataMax[, - which(names(indDataMax) == "Group.1")]
        colnames(indDataMax) <- paste0(gsub(pattern = "...",
                                            x = colnames(indDataMax),
                                            replacement = ".",
                                            fixed = TRUE),
                                       "max")
        indDataMin <- stats::aggregate(platta1, by = list(platta$names),
                                MinMod)
        row.names(indDataMin) <- indDataMin$Group.1
        indDataMin <- indDataMin[, - which(names(indDataMin) == "Group.1")]
        colnames(indDataMin) <- paste0(gsub(pattern = "...",
                                            x = colnames(indDataMin),
                                            replacement = ".",
                                            fixed = TRUE),
                                       "min")
        indData <- cbind(indDataMin, indDataMax)
        QC <- lapply(locus, QCRep, data = indData)
        QC <- as.data.frame(do.call(cbind, QC))
        colnames(QC) <- paste(locus, ".QC", sep = "")
        indData <- indData[,grepl(".1min|.2max", names(indData))]
        QC_Gt <- cbind(indData, QC)
        QC_Gt <- QC_Gt[,order(names(QC_Gt))]
        gt <- QC_Gt[,!grepl("QC", names(QC_Gt))]
        res <- list(QC = QC_Gt,
                    Genotypes = gt,
                    gt.missing = rowSums(is.na(gt)))
        return(res)
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
#' Consensus for this is the extracted as follows:
#' 1. All marker same sex this is the true sex of the sample
#' 2. All markers have NA, the sex is not determined NA is reported
#' 3. In case of mismatch the most common sex i extracted.
#' 4. The Male marker is better as the y-chromosome markers can not separate failed from female.
#' @param data a dataframe with genotype results
#' @return list with input data plus QC score and a vector with most likely sex
sexdetermination <- function(data = data, locus = c("Male",
                                                    "UaY15020",
                                                    "UarY369.4",
                                                    "Y318.2",
                                                    "SMCY")) {
    genSexOrg  <- rep(NA, nrow(data))
    genSexUay  <- rep(NA, nrow(data))
    genSexUar  <- rep(NA, nrow(data))
    genSexY    <- rep(NA, nrow(data))
    genSexSMCY <- rep(NA, nrow(data))

    if(any(grepl("Male", colnames(data)))) {
        ResOrg <- data[,grepl("Male", names(data))]
        ResOrg$merge <- paste(ResOrg[,1], ResOrg[,2])
        genSexOrg <- ifelse(ResOrg$merge == "NA 152"  |
                                ResOrg$merge == "152 152" |
                                ResOrg$merge == "152 NA",
                            yes = "Hona",
                            no = ifelse(ResOrg$merge == "98 152" |
                                            ResOrg$merge == "98 NA"  |
                                            ResOrg$merge == "NA 98"  |
                                            ResOrg$merge == "152 98",
                                        yes = "Hane",
                                        no = NA))
        
    }
    if(any(grepl("UaY15020", colnames(data)))) {
        ResUay    <- data[,grepl("UaY", names(data))]
        genSexUay <- ifelse(is.na(ResUay[,1]),
                            yes = NA,
                            no  = "Hane")
    }
    
    if(any(grepl("UarY369.4", colnames(data)))) {
        ResUar    <- data[,grepl("Uar", names(data))]
        genSexUar <- ifelse(is.na(ResUar[,1]),
                            yes = NA,
                            no  = "Hane")
    }
    if(any(grepl("Y318.2", colnames(data)))) {
        ResY <- data[, grepl("ZFX", names(data))]
        ResY <- cbind(ResY, data[, grepl("318.2", names(data))])  
        ResY$merge <- paste(ResY[,1], ResY[,3])
        genSexY <- ifelse(ResY$merge == "NA 160"  |
                          ResY$merge == "160 NA",
                          yes = "Hona",
                          no  = ifelse(ResY$merge == "NA NA",
                                       yes = NA,
                                       no  = "Hane"))
    }
    if(any(grepl("SMCY", colnames(data)))) {
        ResSMCY    <- data[,grepl("SMCY", names(data))]
        genSexSMCY <- ifelse(is.na(ResSMCY[,1]),
                            yes = NA,
                            no  = "Hane")
    }
    genSex <- data.frame(genSexOrg, genSexUar, genSexUay, genSexY, genSexSMCY)
    sexQC <- function(vec) {
        if(all(is.na(vec))) {
            NA
        } else {        
        vec2 <- na.omit(vec)
        tt <- all(vec2 == vec2[1], na.rm = TRUE)
        ifelse(tt, yes = "OK", no = "Check!")
        }
        }
    QC <- unlist(apply(genSex, MARGIN = 1, sexQC))
    consensusExtract <- function(vec) {
        names(sort(table(factor(vec)), decreasing = TRUE))[1]
    }
    consensus <- unlist(apply(genSex, MARGIN = 1, consensusExtract))
    
    genSex <- cbind(genSex, QC, consensus)
    colnames(genSex) <- c(names(genSex[-ncol(genSex)]), "Sex")
    rownames(genSex) <- rownames(data)
    consensus <- as.data.frame(consensus)
    rownames(consensus) <- rownames(data)
    return(list(QC = genSex,
                Genotypes = consensus))
} 
