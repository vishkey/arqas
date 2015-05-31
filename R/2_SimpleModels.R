#' Obtains the main characteristics of a M/M/1 queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @return
#' Returns the next information of a M/M/1 model:
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_n}: \eqn{C_n}}{\out{<i>P<sub>n</sub>: C<sub>n</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_q}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_q}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub>)</i>}}}
#' @examples
#' #A workstation with a single processor
#' #runs programs with CPU time following
#' #an exponential distribution with mean 3 minutes.
#' #The programs arrives to the workstation following
#' #a Poisson process with an intensity of 15
#' #programs for hour.
#' 
#' M_M_1(lambda=15, mu=60/3)
#' @export
#' @family AnaliticalModels 
M_M_1 <- function(lambda=3, mu=6) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  
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


#' @describeIn Pn Implements the method for a M/M/1 queueing model
#' @method Pn M_M_1
#' @usage NULL
#' @export
Pn.M_M_1 <- function(qm, n) {
  #n must be integer
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  rho <- qm$out$rho
  n <- floor(n)
  if (min(n) < 0) {stop("P(n): Index out of limits: 0:Inf\n")}
  
  return((1-rho)*rho^n)
}

#' @describeIn FW Implements the method for a M/M/1 queueing model
#' @method FW M_M_1
#' @usage NULL
#' @export
FW.M_M_1 <- function(qm, x) {
  lambda <- rate(qm$arrivalDistribution)
  mu <- rate(qm$serviceDistribution)
  return(ifelse(x <0, stop("W(t): Index out of limits: 0:inf\n"), 1-exp((lambda-mu)*x)))
}

#' @describeIn FWq Implements the method for a M/M/1 queueing model
#' @method FWq M_M_1
#' @usage NULL
#' @export
FWq.M_M_1 <- function(qm, x) {
  lambda <- rate(qm$arrivalDistribution)
  mu <- rate(qm$serviceDistribution)
  return(ifelse(x <0, stop("Wq(t): Index out of limits: 0:inf\n"), 1-(lambda/mu)*exp((lambda-mu)*x))) 
}

#' Obtains the main characteristics of a M/M/s queueing model
#'
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param s Number of servers
#' @return
#' Returns the next information of a M/M/s model:
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_n}: \eqn{C_{n}}}{\out{<i>P<sub>n</sub></i>: <i>C<sub>n</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system:  \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>0</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>0</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W - W_q)}}{\out{<i>Eff = W/(W - W<sub>q</sub></i>)}}} 
#' @examples
#' #A workstation with three processors
#' #runs programs with CPU time following
#' #an exponential distribution with mean 3 minutes.
#' #The programs arrives to the workstation following
#' #a Poisson process with an intensity of 15
#' #programs for hour.
#' 
#' M_M_S(lambda=15, mu=60/3, s=3)      
#' @export
#' @family AnaliticalModels 
M_M_S <- function (lambda=3, mu=6, s=2) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(s)      | s <= 0) stop ("Argument 's' must be greather than zero")
  
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


#' @describeIn Pn Implements the method for a M/M/s queueing model
#' @method Pn M_M_S
#' @usage NULL
#' @export
Pn.M_M_S <- function(qm, n) {
  #n must be integer
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  #Check the min interval
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

#' @describeIn FW Implements the method for a M/M/S queueing model
#' @method FW M_M_S
#' @usage NULL
#' @export
FW.M_M_S <- function(qm, x) {
  lambda <- rate(qm$arrivalDistribution)
  mu <- rate(qm$serviceDistribution)
  if (lambda/mu == (qm$servers-1)) {
    return(ifelse(x < 0, 0, 1-(1+qm$out$cn[qm$servers+1]*qm$out$p0*x*mu)*exp(-mu*x)))
  }
  else {
    return(ifelse(x < 0, 0, 1 + ((lambda - qm$servers*mu+mu*FWq(qm, 0))/(qm$servers*mu-lambda-mu))*exp(-mu*x) + ((qm$out$cn[qm$servers+1]*mu*qm$out$p0)/(qm$servers*mu-lambda-mu))*exp(-(qm$servers*mu-lambda)*x)))
  }
  
}

#' @describeIn FWq Implements the method for a M/M/S queueing model
#' @method FWq M_M_S
#' @usage NULL
#' @export
FWq.M_M_S <- function(qm, x) {
  lambda <- rate(qm$arrivalDistribution)
  mu <- rate(qm$serviceDistribution)
  return(ifelse(x < 0, 0, 1-qm$out$cn[qm$servers+1]*qm$out$p0*exp(-(qm$servers*mu-lambda)*x)))
}

#' Obtains the main characteristics of a M/M/1/K queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param k Maximun size of the queue
#' @return
#' Returns the next information of a M/M/1/K model:
#' \item{rho}{Constant coefficient: \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&rho;</i>}}}
#' \item{barlambda}{Effective arrival rate: \ifelse{latex}{\eqn{\bar{\lambda}}}{\out{<i>&#862;&lambda;</i>}}}
#' \item{l}{Mean number of customers in the system: \eqn{L}}
#' \item{lq}{Mean number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @examples
#' #A workstation with a single processor
#' #runs programs with CPU time following
#' #an exponential distribution with mean 3 minutes.
#' #The programs arrive to the workstation following
#' #a Poisson process with an intensity of 15
#' #programs for hour.
#' #The workstation has a limited memory and only
#' #one program is allowed to wait if the processor
#' #is busy.
#' 
#' M_M_1_K(lambda=15, mu=60/3, k=1)
#' @export 
#' @family AnaliticalModels 
M_M_1_K <- function(lambda=3, mu=6, k=2) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(k)      | k < 0) stop ("Argument 'k' must be greather than zero")
  
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


#' @describeIn Pn Implements the method for a M/M/1/K queueing model
#' @method Pn M_M_1_K
#' @usage NULL
#' @export
Pn.M_M_1_K <- function(qm, n) {
  #n must be an integer
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
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

#' @describeIn maxCustomers Implements the method for a M/M/1/K queueing model
#' @method maxCustomers M_M_1_K
#' @usage NULL
#' @export
maxCustomers.M_M_1_K <- function(qm) {
  return(qm$k+1)
}

#' @describeIn Qn Implements the method for a M/M/1/K queueing model
#' @method Qn M_M_1_K
#' @usage NULL
#' @export
Qn.M_M_1_K <- function(qm, n) {
  #n must be an integer
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be an integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > qm$k, 0, (Pn(qm, n)/(1-Pn(qm, qm$k+1))))     
}


#' @describeIn FWq Implements the method for a M/M/1/K queueing model
#' @method FWq M_M_1_K
#' @usage NULL
#' @export
FWq.M_M_1_K <- function(qm, x) {
  FWq.M_M_S_K(qm, x)
}

#' @describeIn FW Implements the method for a M/M/1/K queueing model
#' @method FW M_M_1_K
#' @usage NULL
#' @export
FW.M_M_1_K <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limites: 0:Inf\n")}
  
  mu <- rate(qm$serviceDistribution)
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
#' \item{rho}{Constant coefficient: \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&rho;</i>}}}
#' \item{barlambda}{Effective arrival rate: \ifelse{latex}{\eqn{\bar{\lambda}}}{\out{<i>&#862;&lambda;</i>}}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_n}: \eqn{C_n}}{\out{<i>P<sub>n</sub></i>: <i>C<sub>n</sub></i>}}}
#' \item{pks}{Probability of having \eqn{K+s} customers in the system: \ifelse{latex}{\eqn{P_{K+s}}}{\out{<i>P<sub>K+s</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @examples
#' #A workstation with three processors
#' #runs programs with CPU time following
#' #an exponential distribution with mean 3 minutes.
#' #The programs arrive to the workstation following
#' #a Poisson process with an intensity of 15
#' #programs for hour.
#' #The workstation has a limited memory and only
#' #one program is allowed to wait if the processor
#' #is busy.
#' 
#' M_M_S_K(lambda=15, mu=60/3, s=3, k=1)
#' @export   
#' @family AnaliticalModels 
M_M_S_K <- function(lambda=3, mu=6, s=2, k=3)  {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(s)      | s <= 0) stop ("Argument 's' must be greather than zero")
  if (!is.numeric(k)      | k < 0) stop ("Argument 'k' must be greather than zero")
  
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


#' @describeIn Pn Implements the method for a M/M/s/K queueing model
#' @method Pn M_M_S_K
#' @usage NULL
#' @export
Pn.M_M_S_K <- function(qm, n) {
  #n must be an integer
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
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

#' @describeIn maxCustomers Implements the method for a M/M/S/K queueing model
#' @method maxCustomers M_M_S_K
#' @usage NULL
#' @export
maxCustomers.M_M_S_K <- function(qm) {
  return(qm$k+qm$servers)
}

#' @describeIn Qn Implements the method for a M/M/S/K queueing model
#' @method Qn M_M_S_K
#' @usage NULL
#' @export
Qn.M_M_S_K <- function(qm, n) {
  #n must be an integer
  if (!all.equal(n, as.integer(n))) stop("Q(n): Argument 'n' must be integer")
  if (any(is.na(n))) stop("Q(n): Argument 'n' invalid")
  
  n <- floor(n)
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("Q(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$k+qm$servers-1), 0, (Pn(qm, n)/(1-qm$out$pks)))
}

#' @describeIn FWq Implements the method for a M/M/S/K queueing model
#' @method FWq M_M_S_K
#' @usage NULL
#' @export
FWq.M_M_S_K <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("Wq(t): Index out of limits: 0:Inf")}
  
  mu <- rate(qm$serviceDistribution)
  A <- S <- rep(1, length(x))
  B <- rep(Qn(qm, qm$servers), length(x))
  if ((qm$servers+1) < (qm$k+qm$servers-1)) {
    for(n in (qm$servers+1):(qm$k+qm$servers-1)) {
      A <- A*((qm$servers*mu*x)/(n-qm$servers))
      S <- S + A
      B <- B + Qn(qm, n)*S
      print(n)
    }
  }
  return(1-B*exp(-qm$servers*mu*x)) 
}

#' @describeIn FW Implements the method for a M/M/S/K queueing model
#' @method FW M_M_S_K
#' @usage NULL
#' @export
FW.M_M_S_K <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limits: 0:Inf")}
  mu <- rate(qm$serviceDistribution)
  integrateaux <- function(t) {
    fwaux <- function(x) {FWq(qm, t-x)*mu*exp(mu*x*-1)}
    integrate(fwaux, lower=0, upper=t)
  }
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