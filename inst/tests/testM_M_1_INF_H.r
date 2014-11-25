unstable <- function () {
}

correct <- function () {
  sol <- M_M_1_INF_H(1/40, 1/10, 4)
  
  expect_equal(1-Pn(sol, 0), 71/103)
  expect_equal(Pn(sol, 4), 3/103)
  expect_equal(sol$out$l, 128/103)
}

inputerrors <- function() {
  expect_error(M_M_1_INF_H("a", 2, 4))
  expect_error(M_M_1_INF_H(3, "a", 4))
  expect_error(M_M_1_INF_H(3, 2, "a"))
}

unstable()
correct()
inputerrors()