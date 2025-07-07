test_that("GenerateValidNames works", {
    testfile <- system.file("extdata", "Mix_UarA.csv", package="SSRqc")
    test <- CreateGenotype(testfile, locus = c("G10L", "MU05", "MU59", "Y318_2", "ZFX"))
    expect_equal(GenerateValidNames(test$Genotypes), c("G10L1min", "G10L2max",
                                                 "MU051min", "MU052max",
                                                 "MU591min", "MU592max",
                                                 "Y318_21min", "Y318_22max",
                                                 "ZFX1min", "ZFX2max"))
    expect_equal(GenerateValidNames(data.frame(1:4, 5:8)), c("X14", "X58"))
    expect_equal(GenerateValidNames(data.frame(bäö = 1:2)), "bao")
})
