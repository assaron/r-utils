context("exprs")

test_that("normalizeDE works on limma", {
    t <- read.tsv(system.file("tests/data/mm.tcells.de.tsv", package="rUtils"))
    
    t.norm <- normalizeGeneDE(t)
    expect_equal(head(t.norm$ID, n=3), c("170942", "80876", "15937"))
    expect_true(abs(t.norm$pval[1] / 4.86e-11 - 1) < 0.01)
})