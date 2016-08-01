context("exprs")

test_that("normalizeDE works on limma", {
    t <- read.tsv(system.file("tests/data/mm.tcells.de.tsv", package="rUtils"))
    rownames(t) <- as.character(rownames(t))
    
    t.norm <- normalizeGeneDE(t)
    expect_equal(head(t.norm$ID, n=3), c("170942", "80876", "15937"))
    expect_true(abs(t.norm$pval[1] / 4.86e-11 - 1) < 0.01)
})



test_that("collapseBy works with NAs", {
    t <- read.gct(system.file("tests/data/GSE63040.gct", package="rUtils"))
    fData(t)$symbol <- gsub("na", NA, fData(t)$symbol)
    t1 <- collapseBy(t, fData(t)$symbol)
    expect_equal(length(rownames(t1)), 8)
})

test_that("normalize.rows works with NAs", {
    t <- read.gct(system.file("tests/data/GSE63040.gct", package="rUtils"))
    exprs(t)[2,2] <- NA
    expect_equal(sum(is.na(normalize.rows(exprs(t))[2, ])), 1)
})

test_that("zScore works with NAs", {
    t <- read.gct(system.file("tests/data/GSE63040.gct", package="rUtils"))
    exprs(t)[2,2] <- NA
    expect_equal(sum(is.na(zScore(exprs(t))[2, ])), 1)
})
