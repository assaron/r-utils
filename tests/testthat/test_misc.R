context("misc")

test_that("%o% works", {
    a <- matrix(-12:11, nrow=6)
    b <- apply(a, 1, sum %o% abs)
    c <- apply(a, 1, function(x) sum(abs(x)))
    expect_equal(b, c)
})
