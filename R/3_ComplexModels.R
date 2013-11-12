#' Obtains the main characteristics of a M/M/1/\eqn{\infty}/H queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param h Population size
#' @return
#' Returns the next information of a M/M/1/\eqn{\infty}/H model:
#' \item{rho}{Constant \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{barlambda}{Mean service rate \eqn{\bar{\lambda}}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export 
#' @family AnaliticalModels
M_M_1_INF_H <- function(lambda=3, mu=6, h=5) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (h < 1) stop("Argument 'h' must be equal or greater than one")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  
  obj$servers <- 1
  
  
  obj$h <- h
  
  rho <- lambda/mu
  cn <- c(1, rho*(h-(1:h)+1))
  cn <- cumprod(cn)
  
  p0 <- 1/sum(cn)
  
  l <- sum((1:h)*cn[-1]*p0)
  barlambda <- lambda*(h-l)
  w <- l/barlambda
  wq <- w - (1/mu)
  lq <- barlambda * wq
  eff <- mu * w
  barrho <- rho*(h-l)
  
  obj$out <- list(rho = rho, barrho = barrho, barlambda = barlambda, l=l, lq=lq, wq=wq, w=w, eff=eff, p0=p0, cn=cn)
  oldClass(obj) <- c("M_M_1_INF_H", "M_M_S_INF_H", oldClass(obj))
  return(obj)
}

exportToUI(M_M_1_INF_H, "M/M/1/INF/H", c("numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_1_INF_H
#' @details
#' \code{Pn.M_M_1_INF_H} implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @export
Pn.M_M_1_INF_H <- function(qm, n) {
  Pn.M_M_S_INF_H(qm, n)
}

#' @rdname maxCustomers
#' @method maxCustomers M_M_1_INF_H
#' @details
#' \code{maxCustomers.M_M_1_INF_H} implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @export
maxCustomers.M_M_1_INF_H <- function(qm) {
  maxCustomers.M_M_S_INF_H(qm)
}

#' @rdname Qn
#' @method Qn M_M_1_INF_H
#' @details
#' \code{Qn.M_M_1_INF_H} implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @export
Qn.M_M_1_INF_H <- function(qm, n) {
  Qn.M_M_S_INF_H(qm, n)
}

#' @rdname FWq
#' @method FWq M_M_1_INF_H
#' @details
#' \code{FWq.M_M_1_INF_H} implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @export
FWq.M_M_1_INF_H <- function(qm, x) {
  FWq.M_M_S_INF_H(qm, x)
}

#' @rdname FW
#' @method FW M_M_1_INF_H
#' @details
#' \code{FW.M_M_1_INF_H} implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @export
FW.M_M_1_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limits: 0:Inf\n")}
  
  mu <- rate(qm$serv.distr)
  A <- S <- rep(1, length(x))
  B <- rep(Qn(qm, 0), length(x))
  
  if (qm$h > 1) {
    for(n in 1:(qm$h-1)) {
      A <- A*((mu*x)/n)
      S <- S + A
      B <- B + Qn(qm, n)*S
    }
  }
  return(1-B*exp(-mu*x))
}

#' Obtains the main characteristics of a M/M/s/\eqn{\infty}/H queueing model
#'  
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param s Number of servers
#' @param h Population size
#' @return 
#' Returns the next information of a M/M/s/\eqn{\infty}/H model:
#' \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{barlambda}{Mean effective arrival rate \eqn{\bar{\lambda}}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export
#' @family AnaliticalModels
M_M_S_INF_H <- function(lambda=3, mu=6, s=3, h=5) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (s <= 0) stop ("Argument 's' must be greather than zero")
  if (s > h) stop("Argument 'h' must be equal or greater than argument 's'")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  
  obj$servers <- s
  
  obj$h <- h
  
  rho <- lambda/(s*mu)
  
  cn <- c(1, (lambda/(1:s*mu))*(h-(1:s)+1), rho*(h-((s+1):h)+1))
  cn <- cumprod(cn)
  
  p0 <- 1/sum(cn)
  
  l <- sum((1:h)*(cn[-1]*p0))
  barlambda <- lambda*(h-l)
  w <- l/barlambda
  wq <- w-(1/mu)
  lq <- barlambda*wq
  eff <- mu*w
  barrho <- rho*(h-l)
  
  obj$out <- list(rho = rho, barrho = barrho, barlambda = barlambda, l=l, lq=lq, wq=wq, w=w, eff=eff, p0=p0, cn=cn)
  oldClass(obj) <- c("M_M_S_INF_H", oldClass(obj))
  return(obj)       
}

exportToUI(M_M_S_INF_H, "M/M/s/INF/H", c("numeric", "numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_S_INF_H
#' @details
#' \code{Pn.M_M_S_INF_H} implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @export
Pn.M_M_S_INF_H <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("P(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > qm$h, 0, {
    pn <- c(qm$out$p0, qm$out$cn[-1]*qm$out$p0)
    pn[n+1]
  })
}

#' @rdname maxCustomers
#' @method maxCustomers M_M_S_INF_H
#' @details
#' \code{maxCustomers.M_M_S_INF_H} implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @export
maxCustomers.M_M_S_INF_H <- function(qm) {
  return(qm$h)
}


#' @rdname Qn
#' @method Qn M_M_S_INF_H
#' @details
#' \code{Qn.M_M_S_INF_H} implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @export
Qn.M_M_S_INF_H <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(", n, "): Index out of limits: 0:Inf\n", sep=""))}
  
  return(ifelse(n > (qm$h-1), 0, ((qm$h-n)*Pn(qm, n))/(qm$h-qm$out$l)))
}

#' @rdname FWq
#' @method FWq M_M_S_INF_H
#' @details
#' \code{FWq.M_M_S_INF_H} implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @export
FWq.M_M_S_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("Wq(t): Index out of limits: 0:Inf\n")}
  
  mu <- rate(qm$serv.distr)
  A <- S <- rep(1, length(x))
  B <- rep(Qn(qm, qm$servers), length(x))
  
  if(qm$servers >=qm$h+2) {
    for(n in (qm$servers+1):(qm$h-1)) {
      A <- A*((qm$servers*mu*x)/(n-qm$servers))
      S <- S + A
      B <- B + Qn(qm, n)*S
    }
  }
  return(1-B*exp(-qm$servers*mu*x))  
}

#' @rdname FW
#' @method FW M_M_S_INF_H
#' @details
#' \code{FW.M_M_S_INF_H} implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @export
FW.M_M_S_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limits: 0:Inf")}
  
  mu <- rate(qm$serv.distr)
  integrateaux <- function(t) {
    fwaux <- function(x) {FWq(qm, t-x)*mu*exp(mu*x*-1)}
    integrate(fwaux, lower=0, upper=t)
  }
  #Calculamos W(t) a partir de la integral
  return(unlist(sapply(x, integrateaux)[1,], use.names=FALSE))    
}

#' Obtains the main characteristics of a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param s Number of servers
#' @param h Population size
#' @param y Number of replacements
#' @return
#' Returns the next information of a M/M/s/\eqn{\infty}/H/Y model: 
#' \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{barlambda}{Effective arrival rate \eqn{\bar{\lambda}}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)} \eqn{C_{n}}}
#' \item{p0}{Probability of 0 customers in the system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export
#' @family AnaliticalModels
M_M_S_INF_H_Y <- function(lambda=3, mu=6, s=3, h=5, y=3) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (s <= 0) stop ("Argument 's' must be greather than zero")
  if (h <= 0) stop("Argument 'h' must be greather than zero")
  if (y <= 0) stop("Argument 'y' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  
  if (s > y+h) {
    stop("The sum of arguments 'y' and 'h' must be equal or greater than s\n")
    return()
  }
  else {
    obj$servers <- s
  }
  obj$h <- h
  obj$y <- y
  
  rho <- lambda/(s*mu)
  if (s <= y) {
    cn <- c(1, (h*lambda)/(1:s*mu), rep(rho*h, y-(s+1)+1), rho*(h+y-((y+1):(y+h))+1))      
  }
  else {
    cn <- c(1, (h*lambda)/(1:y*mu), ((h+y-((y+1):s)+1)*lambda)/(((y+1):s)*mu), (h+y-((s+1):(y+h))+1)*rho)
  }
  cn <- cumprod(cn)
  p0 <- 1/sum(cn)
  
  n <- y:(y+h)
  barlambda <- lambda*(h - sum((n-y)*p0*cn[n+1]))
  barrho <- barlambda/(s*mu)
  n <- 1:(y+h)
  l <- sum(n*p0*cn[n+1])
  w <- l/barlambda
  wq <- w-(1/mu)
  lq <- barlambda*wq
  eff <- mu*w
  
  obj$out <- list(rho = rho, barrho = barrho, barlambda = barlambda, l=l, lq=lq, wq=wq, w=w, eff=eff, p0=p0, cn=cn)
  oldClass(obj) <- c("M_M_S_INF_H_Y", "M_M_S_INF_H", oldClass(obj)) 
  return(obj)
}

exportToUI(M_M_S_INF_H_Y, "M/M/s/INF/H with Y replacements", c("numeric", "numeric", "numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_S_INF_H_Y
#' @details
#' \code{Pn.M_M_S_INF_H_Y} implements the method for a M/M/s/\eqn{\infty}/H/Y queueing model
#' @export
Pn.M_M_S_INF_H_Y <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("P(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$y + qm$h), 0, {
    pn <- c(qm$out$p0, qm$out$cn[-1]*qm$out$p0)
    pn[n+1] 
  })          
}

#' @rdname maxCustomers
#' @method maxCustomers M_M_S_INF_H_Y
#' @details
#' \code{maxCustomers.M_M_S_INF_H_Y} implements the method for a M/M/s/\eqn{\infty}/H/Y queueing model
#' @export
maxCustomers.M_M_S_INF_H_Y<- function(qm) {
  return(qm$y + qm$h)
}


#' @rdname Qn
#' @method Qn M_M_S_INF_H_Y
#' @details
#' \code{Qn.M_M_S_INF_H_Y} implements the method for a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' @export
Qn.M_M_S_INF_H_Y <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$y+qm$h-1), 0, {
    emes <- qm$y:(qm$y+qm$h)
    sumemes <- sum((emes-qm$y)*Pn(qm,emes))       
    ifelse(n <= (qm$y-1), qm$h*Pn(qm, n)/(qm$h-sumemes), (qm$h+qm$y-n)*Pn(qm, n)/(qm$h-sumemes))
  }) 
}

#' @rdname FWq
#' @method FWq M_M_S_INF_H_Y
#' @details
#' \code{FWq.M_M_S_INF_H_Y} implements the method for a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' @export
FWq.M_M_S_INF_H_Y <- function(qm, x) {
  return(FWq.M_M_S_INF_H(qm,x))
}

#' @rdname FW
#' @method FW M_M_S_INF_H_Y
#' @details
#' \code{FW.M_M_S_INF_H_Y} implements the method for a M/M/s/\eqn{\infty}/H/ with Y replacements queueing model
#' @export
FW.M_M_S_INF_H_Y <- function(qm, x) {
  return(FW.M_M_S_INF_H(qm, x))
}

#' Obtains the main characteristics of a M/M/\eqn{\infty} queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @return
#' Returns the next information of a M/M/\eqn{\infty} model: 
#' \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}} (\eqn{L_{q}=0} in this model)}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}} (\eqn{W_{q}=0} in this model)}
#' \item{eff}{Efficiency of the system  \eqn{Eff = W/(W-W_q)}}
#' @export 
#' @family AnaliticalModels
M_M_INF <- function(lambda=3, mu=6) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  
  barrho <- 0
  p0 <- exp(-lambda/mu)        
  l <- lambda/mu
  lq <- 0
  w <- 1/mu
  wq <- 0
  eff <- 1
  obj$out <- list(rho = l, barrho = barrho, l=l, lq=lq, wq=wq, w=w, eff=eff, p0=p0)
  oldClass(obj) <- c("M_M_INF", oldClass(obj))
  return(obj)
}

exportToUI(M_M_INF, "M/M/INF", c("numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_INF
#' @details
#' \code{Pn.M_M_INF} implements the method for a M_M_INF queueing model
#' @export
Pn.M_M_INF <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop("P(n): Index out of limits: 0:Inf\n")}
  
  lambda <- rate(qm$arr.distr)
  mu <- rate(qm$serv.distr)        
  cn <- c(1, lambda/((1:maxval)*mu))
  cn <- cumprod(cn)
  
  pn <- c(qm$out$p0, cn[-1]*qm$out$p0)
  return(pn[n+1])
}

#' @rdname FW
#' @method FW M_M_INF
#' @details
#' \code{FW.M_M_INF} implements the method for a M/M/\eqn{\infty} queueing model
#' @export
FW.M_M_INF <- function(qm, x) {
  rep(0, length(x))
}

#' @rdname FWq
#' @method FWq M_M_INF
#' @details
#' \code{FWq.M_M_INF} implements the method for a M/M/\eqn{\infty} queueing model
#' @export
FWq.M_M_INF <- function(qm, x) {
  rep(0, length(x))
}
#' Print the main characteristics of a M/M/S/\eqn{\infty}/H model
#' @param x a M_M_S_INF_H object
#' @param ... Further arguments passed to or from other methods.
#' @method print M_M_S_INF_H
#' @keywords internal
#' @export
print.M_M_S_INF_H <- function(x, ...) {
  cat("Model: ", class(x)[1])
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$barrho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
}

#' Print the main characteristics of a M/M/\eqn{\infty} model
#' @param x a M_M_INF object
#' @param ... Further arguments passed to or from other methods.
#' @method print M_M_INF
#' @keywords internal
#' @export
print.M_M_INF <- function(x, ...) {
  cat("Model: ", class(x)[1])
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$barrho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
}