createGenotypeMix1v2 <- function(file = file, na.strings = c(0, "NA",
                                                           "No peaks", "No peaks in locus", "No peaks in bins")) {
    platta <- load_data(file, na_strings = na.strings)
    platta <- platta[,!grepl("NED", names(platta))]
    platta$names <- extract_names(platta[,1])
    indDataMax <- aggregate(platta[,-1], by = list(platta$names), max, na.rm = TRUE)
    indDataMax[indDataMax == -Inf] <- NA
    row.names(indDataMax)  <- indDataMax$Group.1
    colnames(indDataMax) <- paste0(gsub(pattern = "...",
                                        x = colnames(indDataMax),
                                        replacement = ".",
                                        fixed = TRUE),
                                   "max")
    indDataMin <- aggregate(platta[,-1], by = list(platta$names), min, na.rm = TRUE)
    indDataMin[indDataMin == Inf] <- NA
    row.names(indDataMin)  <- indDataMin$Group.1
    colnames(indDataMin) <- paste0(gsub(pattern = "...",
                                        x = colnames(indDataMin),
                                        replacement = ".",
                                        fixed = TRUE),
                                   "min")
    indData <- cbind(indDataMin, indDataMax)
    indData <- indData[,grepl("MU" ,names(indData)) | grepl("UaY", names(indData))]
    indData <- indData[,order(colnames(indData))]
    MU09.QC <- qc_rep(indData, locus = "MU09")
    MU10.QC <- qc_rep(indData, locus = "MU10")
    MU05.QC <- qc_rep(indData, locus = "MU05")
    MU23.QC <- qc_rep(indData, locus = "MU23")
    UaY15020.QC <- qc_rep(indData, locus = "UaY15020")
    QC <- cbind(indData[,c("MU09.1min", "MU09.2max")], MU09.QC,
                indData[,c("MU10.1min", "MU10.2max")], MU10.QC,
                indData[,c("MU05.1min", "MU05.2max")], MU05.QC,
                indData[,c("MU23.1min", "MU23.2max")], MU23.QC,
                indData[,c("UaY15020.1min", "UaY15020.2max")], UaY15020.QC)
    gt <- QC[,!grepl("QC", names(QC))]
    res <- list(QC = QC,
                Genotypes = gt,
                gt.sum = rowSums(is.na(gt)))
    return(res)
}


createGenotypeMix2v2 <- function(file = file,
                               na.strings = c(0, "NA",
                                              "No peaks",
                                              "No peaks in locus",
                                              "No peaks in bins")) {
    # Create QC and genotype results from file with genotype results.
    # Most commonly the file is an export from Geneious or similar software.
    # This assumes that the structure of the file corresponds to the mix2 of
    # the Swedish bearSSR project.
    platta <- load_data(file, na_strings = na.strings)
    platta <- platta[,!grepl("NED", names(platta))]
    platta$names <- extract_names(platta[,1])
    indDataMax <- aggregate(platta[,-1], by = list(platta$names), max, na.rm = TRUE)
    indDataMax[indDataMax == -Inf] <- NA
    row.names(indDataMax)  <- indDataMax$Group.1
    colnames(indDataMax) <- paste0(gsub(pattern = "...",
                                        x = colnames(indDataMax),
                                        replacement = ".",
                                        fixed = TRUE
    ),
    "max")
    indDataMin <- aggregate(platta[,-1], by = list(platta$names), min, na.rm = TRUE)
    indDataMin[indDataMin == Inf] <- NA
    row.names(indDataMin)  <- indDataMin$Group.1
    colnames(indDataMin) <- paste0(gsub(pattern = "...",
                                        x = colnames(indDataMin),
                                        replacement = ".",
                                        fixed = TRUE),
                                   "min")
    indData <- cbind(indDataMin, indDataMax)
    indData <- indData[,grepl("MU", names(indData)) | grepl("G10", names(indData)) | grepl("UarY", names(indData))]
    indData <- indData[,order(colnames(indData))]
    #indData
    G10.QC <- qc_rep(indData, locus = "G10L")
    MU50.QC <- qc_rep(indData, locus = "MU50")
    MU51.QC <- qc_rep(indData, locus = "MU51")
    MU59.QC <- qc_rep(indData, locus = "MU59")
    UarY369.4.QC <- qc_rep(indData, locus = "UarY369.4")
    #indData <- indData[,grepl("max", names(indData))]
    QC <- cbind(indData[,c("G10L.1min", "G10L.2max")], G10.QC,
                indData[,c("MU50.1min", "MU50.2max")], MU50.QC,
                indData[,c("MU51.1min", "MU51.2max")], MU51.QC,
                indData[,c("MU59.1min", "MU59.2max")], MU59.QC,
                indData[,c("UarY369.4.1min", "UarY369.4.2max")], UarY369.4.QC)
    gt <- QC[,!grepl("QC", names(QC))]
    res <- list(QC = QC,
                Genotypes = gt,
                gt.sum = rowSums(is.na(gt)))
    return(res)
}
