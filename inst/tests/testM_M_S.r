unstable <- function () {
  expect_error(M_M_S(8, 3, 2))
}

correct <- function () {
  sol <- M_M_S(3, 2, 4)
  
  expect_equal(sol$out$l, 1398/905)
  expect_equal(sol$out$lq, 81/1810)
  expect_equal(sol$out$w, 466/905)
  expect_equal(sol$out$wq, 27/1810)
  expect_equal(sol$out$eff, 932/905)
  expect_equal(sol$out$rho, 3/8)
}

inputerrors <- function() {
  expect_error(M_M_S("a", 2, 4))
  expect_error(M_M_S(3, "a", 4))
  expect_error(M_M_S(3, 2, "a"))
}

unstable()
correct()
inputerrors()