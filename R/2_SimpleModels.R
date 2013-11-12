#' Obtains the main characteristics of a M/M/1 queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @return
#' Returns the next information of a M/M/1 model:
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export
#' @family AnaliticalModels 
M_M_1 <- function(lambda=3, mu=6) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  obj$servers <- 1
  
  rho <- lambda/mu
  if (rho >= 1) {stop("non-stationary model\n")}
  
  l <- lambda/(mu-lambda)
  wq <- lambda/(mu*(mu-lambda)) 
  lq <- lambda*wq
  w <- l/lambda
  eff <- mu*w
  
  obj$out <- list(rho = rho, barrho=rho, l=l, lq=lq, wq=wq, w=w, eff=eff)
  oldClass(obj) <- c("M_M_1", "M_M_S", oldClass(obj))
  return(obj)        
}

exportToUI(M_M_1, "M/M/1", c("numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_1
#' @details
#' \code{Pn.M_M_1} implements the method for a M/M/1 queueing model
#' @export
Pn.M_M_1 <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  rho <- qm$out$rho
  n <- floor(n)
  if (min(n) < 0) {stop("P(n): Index out of limits: 0:Inf\n")}
  
  return((1-rho)*rho^n)
}

#' @rdname FW
#' @method FW M_M_1
#' @details
#' \code{FW.M_M_1} implements the method for a M/M/1 queueing model
#' @export
FW.M_M_1 <- function(qm, x) {
  lambda <- rate(qm$arr.distr)
  mu <- rate(qm$serv.distr)
  return(ifelse(x <0, stop("W(t): Index out of limits: 0:inf\n"), 1-exp((lambda-mu)*x)))
}

#' @rdname FWq
#' @method FWq M_M_1
#' @details
#' \code{FWq.M_M_1} implements the method for a M/M/1 queueing model
#' @export
FWq.M_M_1 <- function(qm, x) {
  lambda <- rate(qm$arr.distr)
  mu <- rate(qm$serv.distr)
  return(ifelse(x <0, stop("Wq(t): Index out of limits: 0:inf\n"), 1-(lambda/mu)*exp((lambda-mu)*x))) 
}

#' Obtains the main characteristics of a M/M/s queueing model
#'
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param s Number of servers
#' @return
#' Returns the next information of a M/M/s model:
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)} \eqn{C_{n}}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system  \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W - W_q)}}
#' @export       
#' @family AnaliticalModels 
M_M_S <- function (lambda=3, mu=6, s=2) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (s <= 0) stop ("Argument 's' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  obj$servers <- s
  
  rho <- lambda/(mu*s)
  
  if (rho >= 1) {stop("non-stationary model")}
  
  if (s < 2)
    cn <- c(1, lambda/(s*mu-lambda))
  else
    cn <- c(1, lambda/((1:(s-1))*mu), lambda/(s*mu-lambda))
  
  cn <- cumprod(cn)
  p0 <- 1/sum(cn)
  
  lq <- (cn[s+1]*lambda*p0)/(s*mu-lambda)
  wq <- lq/lambda
  w <- wq + 1/mu
  l <- lambda * w
  eff <- mu * w
  
  oldClass(obj) <- c("M_M_S", oldClass(obj))
  obj$out <- list(rho = rho, barrho=rho, cn = cn, p0 = p0, l=l, lq=lq, wq=wq, w=w, eff=eff)
  return(obj)        
}

exportToUI(M_M_S, "M/M/s", c("numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_S
#' @details
#' \code{Pn.M_M_S} implements the method for a M/M/S queueing model
#' @export
Pn.M_M_S <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  #comprobamos que el indice menor sea correcto
  minval <- min(n)
  if (minval < 0) stop("P(n): Parameter 'n' must be positive")
  
  maxval <- max(n)
  pn <- c(qm$out$p0, qm$out$p0 * qm$out$cn[-1][-qm$servers])
  if (maxval > (qm$servers-1)) {
    times <- maxval - (qm$servers-1)
    pnadd <- c(pn[qm$servers], rep(qm$out$rho, times))
    pnadd <- cumprod(pnadd)
    pn <- c(pn, pnadd[-1])
  }
  return(pn[n+1])
}

#' @rdname FW
#' @method FW M_M_S
#' @details
#' \code{FW.M_M_S} implements the method for a M/M/S queueing model
#' @export
FW.M_M_S <- function(qm, x) {
  lambda <- rate(qm$arr.distr)
  mu <- rate(qm$serv.distr)
  if (lambda/mu == (qm$servers-1)) {
    return(ifelse(x < 0, 0, 1-(1+qm$out$cn[qm$servers+1]*qm$out$p0*x*mu)*exp(-mu*x)))
  }
  else {
    return(ifelse(x < 0, 0, 1 + ((lambda - qm$servers*mu+mu*FWq(qm, 0))/(qm$servers*mu-lambda-mu))*exp(-mu*x) + ((qm$out$cn[qm$servers+1]*mu*qm$out$p0)/(qm$servers*mu-lambda-mu))*exp(-(qm$servers*mu-lambda)*x)))
  }
  
}

#' @rdname FWq
#' @method FWq M_M_S
#' @details
#' \code{FWq.M_M_S} implements the method for a M/M/S queueing model
#' @export
FWq.M_M_S <- function(qm, x) {
  lambda <- rate(qm$arr.distr)
  mu <- rate(qm$serv.distr)
  return(ifelse(x < 0, 0, 1-qm$out$cn[qm$servers+1]*qm$out$p0*exp(-(qm$servers*mu-lambda)*x)))
}

#' Obtains the main characteristics of a M/M/1/K queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param k Maximun size of the queue
#' @return
#' Returns the next information of a M/M/1/K model:
#' \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{barlambda}{Effective arrival rate \eqn{\bar{\lambda}}}
#' \item{l}{Expected mean number of customers in the system \eqn{L}}
#' \item{lq}{Expected mean number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export 
#' @family AnaliticalModels 
M_M_1_K <- function(lambda=3, mu=6, k=2) {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (k <= 0) stop ("Argument 'k' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  obj$servers <- 1
  obj$k <- k
  
  rho <- lambda/mu
  
  if (rho != 1) {
    l <- (rho/(1-rho)) - ((k+2)*rho^(k+2)/(1-rho^(k+2)))
    barlambda <- lambda*(rho^(k+1)-1)/(rho^(k+2)-1)
  }
  else {
    l <- (k+1)/2
    barlambda <- lambda*(k+1)/(k+2)
  }
  w <- l/barlambda
  wq <- w-1/mu
  lq <- barlambda*wq
  eff <- mu*w
  barrho <- barlambda/mu
  oldClass(obj) <- c("M_M_1_K", "M_M_S_K", oldClass(obj))
  
  obj$out <- list(rho = rho, barrho = barrho, barlambda = barlambda, l=l, lq=lq, wq=wq, w=w, eff=eff)
  return(obj)
}

exportToUI(M_M_1_K, "M/M/1/K", c("numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_1_K
#' @details
#' \code{Pn.M_M_1_K} implements the method for a M/M/1/K queueing model
#' @export
Pn.M_M_1_K <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("P(n): Index out of limits: 0:Inf\n", sep=""))}

  ifelse(n > (qm$k+1), 0, {
      rho <- qm$out$rho
      if (rho == 1) {
        rep(1/(qm$k+2), length(n))
      } else {
        ((rho-1)/(rho^(qm$k+2)-1))*rho^n
      } 
  })
}

#' @rdname maxCustomers
#' @method maxCustomers M_M_1_K
#' @details
#' \code{maxCustomers.M_M_1_K} implements the method for a M/M/1/K queueing model
#' @export
maxCustomers.M_M_1_K <- function(qm) {
  return(qm$k+1)
}

#' @rdname Qn
#' @method Qn M_M_1_K
#' @details
#' \code{Qn.M_M_1_K} implements the method for a M/M/1/K queueing model
#' @export
Qn.M_M_1_K <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > qm$k, 0, (Pn(qm, n)/(1-Pn(qm, qm$k+1))))     
}

#Para el caso FWq usamos la version del M_M_S_K que nos hace el calculo de manera eficiente
#' @rdname FWq
#' @method FWq M_M_1_K
#' @details
#' \code{FWq.M_M_1_K} implements the method for a M/M/1/K queueing model
#' @export
FWq.M_M_1_K <- function(qm, x) {
  FWq.M_M_S_K(qm, x)
}

#La FW si que la implementamos ya que para el M_M_1_K podemos evitar la integral
#' @rdname FW
#' @method FW M_M_1_K
#' @details
#' \code{FW.M_M_1_K} implements the method for a M/M/1/K queueing model
#' @export
FW.M_M_1_K <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limites: 0:Inf\n")}
  
  mu <- rate(qm$serv.distr)
  A <- S <- rep(1, length(x))
  B <- rep(Qn(qm, 0), length(x))
  for(n in 1:qm$k) {
    A <- A*((mu*x)/n)
    S <- S + A
    B <- B + Qn(qm, n)*S
  }
  return(1 - B*exp(-mu*x))
}

#' Obtains the main characteristics of a M/M/S/k queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param s Number of servers
#' @param k Maximun size of the queue
#' @return
#' Returns the next information of a M/M/S/K model:
#' \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
#' \item{barlambda}{Effective arrival rate \eqn{\bar{\lambda}}}
#' \item{cn}{Constant coefficients used in the computation of \eqn{P(n)}}
#' \item{pks}{Probability of \eqn{K+s} customers in the system \eqn{P_{K+s}}}
#' \item{p0}{Probability of empty system \eqn{P_{0}}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue\eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-W_q)}}
#' @export   
#' @family AnaliticalModels 
M_M_S_K <- function(lambda=3, mu=6, s=2, k=3)  {
  if (lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (s <= 0) stop ("Argument 's' must be greather than zero")
  if (k <= 0) stop ("Argument 'k' must be greather than zero")
  
  obj <- MarkovianModel(Exp(lambda), Exp(mu))
  obj$servers <- s
  obj$k <- k
  
  rho <- lambda/(s*mu)
  if (s < 2)
    cn <- c(1)
  else
    cn <- c(1, lambda/((1:(s-1))*mu))
  
  if (rho == 1) {
    cn <- c(cn, k+1)
    cn <- cumprod(cn)
    
    p0 <- 1/sum(cn)
    lq <- ((k*(k+1))/2)*(cn[s]*p0)
  }
  else {
    cn <- c(cn, (rho-rho^(k+2))/(1-rho))
    cn <- cumprod(cn)
    
    p0 <- 1/sum(cn)
    lq <- (((1+k*rho^(k+1)-(k+1)*rho^k)*rho^2)/((1-rho)^2))*(cn[s]*p0)
  }
  pks <- cn[s]*p0*rho^(k+1)
  barlambda <- lambda*(1-pks)
  barrho <- rho*(1-pks)
  wq <- lq/barlambda
  w <- wq + 1/mu
  l <- barlambda*w
  eff <- mu*w
  oldClass(obj) <- c("M_M_S_K", oldClass(obj))
  
  obj$out <- list(rho = rho, barrho = barrho, barlambda = barlambda, cn = cn, p0 = p0, pks = pks, l=l, lq=lq, wq=wq, w=w, eff=eff)
  return(obj)
}

exportToUI(M_M_S_K, "M/M/s/K", c("numeric", "numeric", "numeric", "numeric"), "markovian")

#' @rdname Pn
#' @method Pn M_M_S_K
#' @details
#' \code{Pn.M_M_S_K} implements the method for a M/M/S/K queueing model
#' @export
Pn.M_M_S_K <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  n <- floor(n)
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("P(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$k+qm$servers), 0, {
      pn <- c(qm$out$p0, qm$out$cn[-1][-qm$servers]*qm$out$p0)
      pnadd <- c(pn[qm$servers], rep(qm$out$rho, qm$servers+qm$k))
      pnadd <- cumprod(pnadd)
      pn <- c(pn, pnadd[-1])
      pn[n+1]
  })
}

#' @rdname maxCustomers
#' @method maxCustomers M_M_S_K
#' @details
#' \code{maxCustomers.M_M_S_K} implements the method for a M/M/S/K queueing model
#' @export
maxCustomers.M_M_S_K <- function(qm) {
  return(qm$k+qm$servers)
}

#' @rdname Qn
#' @method Qn M_M_S_K
#' @details
#' \code{Qn.M_M_S_K} implements the method for a M/M/S/K queueing model
#' @export
Qn.M_M_S_K <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  n <- floor(n)
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$k+qm$servers-1), 0, (Pn(qm, n)/(1-qm$out$pks)))
}

#' @rdname FWq
#' @method FWq M_M_S_K
#' @details
#' \code{FWq.M_M_S_K} implements the method for a M/M/S/K queueing model
#' @export
FWq.M_M_S_K <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("Wq(t): Index out of limits: 0:Inf")}
  
  mu <- rate(qm$serv.distr)
  A <- S <- rep(1, length(x))
  B <- rep(Qn(qm, qm$servers), length(x))
  for(n in (qm$servers+1):(qm$k+qm$servers-1)) {
    A <- A*((qm$servers*mu*x)/(n-qm$servers))
    S <- S + A
    B <- B + Qn(qm, n)*S
  }
  return(1-B*exp(-qm$servers*mu*x)) 
}

#' @rdname FW
#' @method FW M_M_S_K
#' @details
#' \code{FW.M_M_S_K} implements the method for a M/M/S/K queueing model
#' @export
FW.M_M_S_K <- function(qm, x) {
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

#' Print the main characteristics of a M_M_S_K model
#' @param x M_M_S_K class object
#' @param ... Further arguments passed to or from other methods.
#' @method print M_M_S_K
#' @keywords internal
#' @export
print.M_M_S_K <- function(x, ...) {
  cat("Model: ", class(x)[1])
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$barrho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
}