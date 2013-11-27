Ciclic <- function(n) {
  diag(n)[c(2:n,1),]
}

Tamdem <- function(n) {
  diag(n+1)[2:(n+1), 1:n]
}