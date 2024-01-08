test_that("ExtractNames works", {
    expect_equal(ExtractNames("test_test"), "test")
    expect_equal(ExtractNames("test.test", sep = "."), "test")
})
