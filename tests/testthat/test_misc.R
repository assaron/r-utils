context("misc")

test_that("%o% works", {
    a <- matrix(-12:11, nrow=6)
    b <- apply(a, 1, sum %o% abs)
    c <- apply(a, 1, function(x) sum(abs(x)))
    expect_equal(b, c)
})

test_that("ulength works", {
    expect_equal(ulength(c("a", "b", "a")), 2)
})

test_that("replaceNA works", {
    expect_equal(replaceNA(c(1, NA, 3), 0), c(1, 0, 3))
})

test_that("%f% works", {
    expect_equal("x %s y" %f% "z", "x z y")
})
