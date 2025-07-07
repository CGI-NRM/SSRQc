test_that("CreateGenotype function works correctly", {
    testfile <- system.file("extdata", "Mix_UarA.csv", package="SSRqc")
    test_that("CreateGenotype returns correct output with valid input", {
        file <- testfile
        locus <- c("MU59","G10L","MU05","ZFX","Y318_2")
        result <- CreateGenotype(file, locus)
        expect_type(result, "list")
        expect_length(result, 3)
        expect_s3_class(result$Genotypes, "data.frame")
        expect_s3_class(result$QC, "data.frame")
        expect_type(result$genotype_missing, "double")
    })
    test_that("CreateGenotype handles empty data correctly", {
        file <- testfile
        locus <- c("Locus1", "Locus2")
        result <- CreateGenotype(file, locus)
        expect_null(result)
    })
    test_that("CreateGenotype handles not allowed strings correctly", {
        testfileunbinned <- system.file("extdata", "Mix_UarA_unbinned.csv", package="SSRqc")
        file <- testfileunbinned
        locus <- c("MU59", "Y318_2")
        notAllowed <- c("Unbinned peaks", "Too many alleles", "Unbinned peaks in locus", "Peaks outside loci")
        expect_warning(CreateGenotype(file, locus, notAllowed = notAllowed))
    })
})
