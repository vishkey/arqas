vals <- runif(10000)

foo <- function() {
  for (i in 1:10000) {
    indices_validos <- which(vals != -1)
    indice_min <- which.min(vals[indices_validos])
    vals[indices_validos] <- vals[indices_validos] - vals[indice_min]
    vals[indice_min] <- -1
  }
}

foo2 <- function() {
  for (i in 1:10000) {
    indice_min <- which.min(vals)
    vals <- vals -vals[indice_min]
    vals[indice_min] <- NA
  }
}

print(system.time(foo()))
print(system.time(foo2()))