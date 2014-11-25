unstable <- function () {
}

correct <- function () {
  sol <- M_M_1_K(4, 1, 3)
  
  expect_equal(Pn(sol, 0), 1/341)
  expect_equal(sol$out$barrho, 340/341)
  expect_equal(sol$out$w, 1252/340)
  expect_equal(sol$out$lq, 912/341)
}

inputerrors <- function() {
  expect_error(M_M_1_K("a", 2, 4))
  expect_error(M_M_1_K(3, "a", 4))
  expect_error(M_M_1_K(3, 2, "a"))
}

unstable()
correct()
inputerrors()