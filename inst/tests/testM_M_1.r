unstable <- function () {
    expect_error(M_M_1(6,3))
}

correct <- function () {
    sol <- M_M_1(3, 6)
    
    expect_equal(sol$out$l, 1)
    expect_equal(sol$out$lq, 0.5)
    expect_equal(sol$out$w, 1/3)
    expect_equal(sol$out$wq, 1/6)
    expect_equal(sol$out$eff, 2)
    expect_equal(sol$out$rho, 0.5)
}

inputerrors <- function() {
    expect_error(M_M_1(1, 1))
    expect_error(M_M_1("a", 3))
    expect_error(M_M_1(3, "a"))
}

unstable()
correct()
inputerrors()