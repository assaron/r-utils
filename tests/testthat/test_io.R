test_that("read.tsv detects row names", {
    t <- read.tsv(system.file("tests/data/mm.tcells.de.tsv", package="rUtils"))
    expect_equal(head(rownames(t), n=3), c("170942", "80876", "15937"))
})