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
#' @return a list of three dataframes: Genotypes, QC and number of
#' missing data points
#' @param file name of the file with replicate data
#' @param locus names of the loci to be analysed

createGT <- function(file = file, locus = c("Locus1", "Locus2")) {
    platta <- load_data(file)
    if(!is.null(platta)) {
        platta$names <- extract_names(platta[,1])
        togrep <- paste(locus,collapse="|")
        platta1 <- platta[,grepl(togrep, names(platta))]
        
        indDataMax <- aggregate(platta1, by = list(platta$names),
                                max_mod)
        row.names(indDataMax) <- indDataMax$Group.1
        colnames(indDataMax) <- paste0(gsub(pattern = "...",
                                            x = colnames(indDataMax),
                                            replacement = ".",
                                            fixed = TRUE),
                                       "max")
        indDataMin <- aggregate(platta1, by = list(platta$names),
                                min_mod)
        row.names(indDataMin) <- indDataMin$Group.1
        colnames(indDataMin) <- paste0(gsub(pattern = "...",
                                            x = colnames(indDataMin),
                                            replacement = ".",
                                            fixed = TRUE),
                                       "min")
        indData <- cbind(indDataMin, indDataMax)
        QC <- lapply(locus, qc_rep, data = indData)
        QC <- as.data.frame(do.call(cbind, QC))
        colnames(QC) <- paste(locus, ".QC", sep = "")
        indData <- indData[,grepl(".1min|.2max", names(indData))]
        QC_Gt <- cbind(indData, QC)
        QC_Gt <- QC_Gt[,order(names(QC_Gt))]
        gt <- QC_Gt[,!grepl("QC", names(QC_Gt))]
        res <- list(QC = QC_Gt,
                    Genotypes = gt,
                    gt.sum = rowSums(is.na(gt)))
        return(res)
    } else {
        cat("")
    }
    
}

    
    
