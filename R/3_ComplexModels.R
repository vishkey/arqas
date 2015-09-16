#' Obtains the main characteristics of a M/M/1/\eqn{\infty}/H queueing model
#' 
#' @param lambda Mean arrival rate
#' @param mu Mean service rate 
#' @param h Population size
#' @return
#' Returns the next information of a M/M/1/\eqn{\infty}/H model :
#' \item{rho}{Constant: \eqn{\lambda/\mu}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{{\out{<i>&#862;&rho;</i>}}}}
#' \item{barlambda}{Mean effective arrival rate: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&lambda;</i>}}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_n}: \eqn{C_n}}{\out{P<sub>n</sub>: <i>C<sub>n</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_q}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_q}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System Efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub>)</i>}}}
#' @examples
#' # A computer system with five workstations must
#' # make a back- up from time to time, preventing
#' # users from using the system.
#' # The time between ending a back -up until the
#' # next begins is random and follows an
#' # exponential distribution with mean 2 hours. 
#' # The duration of the backup is also random and
#' # follows an exponential distribution with mean
#' # 5 minutes.
#' # There exists a single tape drive to perform the
#' # back-up process and the station will wait if it is
#' # busy.
#' 
#' M_M_1_INF_H(lambda =1/2, mu=60/5, h=5)
#' 
#' 
#' @export 
#' @family AnaliticalModels
M_M_1_INF_H <- function(lambda=1/2, mu=60/5, h=5) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(h)      | h < 1) stop("Argument 'h' must be equal or greater than one")
  
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

#' @describeIn Pn implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @method Pn M_M_1_INF_H
#' @usage NULL
#' @export
Pn.M_M_1_INF_H <- function(qm, n) {
  Pn.M_M_S_INF_H(qm, n)
}

#' @describeIn maxCustomers Implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @method maxCustomers M_M_1_INF_H
#' @usage NULL
#' @export
maxCustomers.M_M_1_INF_H <- function(qm) {
  maxCustomers.M_M_S_INF_H(qm)
}

#' @describeIn Qn Implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @method Qn M_M_1_INF_H
#' @usage NULL
#' @export
Qn.M_M_1_INF_H <- function(qm, n) {
  Qn.M_M_S_INF_H(qm, n)
}

#' @describeIn FWq Implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @method FWq M_M_1_INF_H
#' @usage NULL
#' @export
FWq.M_M_1_INF_H <- function(qm, x) {
  FWq.M_M_S_INF_H(qm, x)
}

#' @describeIn FW Implements the method for a M/M/1/\eqn{\infty}/H queueing model
#' @method FW M_M_1_INF_H
#' @usage NULL
#' @export
FW.M_M_1_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limits: 0:Inf\n")}
  
  mu <- rate(qm$serviceDistribution)
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
#' \item{rho}{Constant coefficient: \eqn{\lambda/\mu}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&rho;</i>}}}
#' \item{barlambda}{Mean effective arrival rate: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&lambda;</i>}}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_{n}}: \eqn{C_n}}{\out{P<sub>n</sub>: <i>C<sub>n</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_q}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_q}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @export
#' @examples
#' # A computer system with five workstations must
#' # make a back- up from time to time, preventing
#' # users from using the system.
#' # The time between ending a back -up until the
#' # next begins is random and follows an
#' # exponential distribution with mean 2 hours. 
#' # The duration of the backup is also random and
#' # follows an exponential distribution with mean
#' # 5 minutes.
#' # There are 2 tape drives to perform the
#' # back-up process, the station remained on hold
#' # if both were occupied.
#' 
#' M_M_S_INF_H(lambda=1/2, mu=60/5, s=2, h=5)
#' @family AnaliticalModels
M_M_S_INF_H <- function(lambda=1/2, mu=60/5, s=2, h=5) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(s)      | s <= 0) stop ("Argument 's' must be greather than zero")
  if (!is.numeric(h)      | (s > h)) stop("Argument 'h' must be equal or greater than argument 's'")
  
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

#' @describeIn Pn Implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @usage NULL 
#' @method Pn M_M_S_INF_H
#' @export
Pn.M_M_S_INF_H <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
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
#' @examples
#' maxCustomers(M_M_S_INF_H(lambda=1/2, mu=60/5, s=2, h=5))
#' @export
maxCustomers.M_M_S_INF_H <- function(qm) {
  return(qm$h)
}


#' @describeIn Qn Implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @method Qn M_M_S_INF_H
#' @usage NULL
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

#' @describeIn FWq Implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @method FWq M_M_S_INF_H
#' @usage NULL
#' @export
FWq.M_M_S_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("Wq(t): Index out of limits: 0:Inf\n")}
  
  mu <- rate(qm$serviceDistribution)
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

#' @describeIn FW Implements the method for a M/M/s/\eqn{\infty}/H queueing model
#' @method FW M_M_S_INF_H
#' @usage NULL
#' @export
FW.M_M_S_INF_H <- function(qm, x) {
  minval <- min(x)
  if (minval < 0) {stop("W(t): Index out of limits: 0:Inf")}
  
  mu <- rate(qm$serviceDistribution)
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
#' \item{rho}{Constant coefficient: \eqn{\lambda/\rho}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&rho;</i>}}}
#' \item{barlambda}{Effective arrival rate: \ifelse{latex}{\eqn{\bar{\lambda}}}{\out{<i>&#862;&lambda;</i>}}}
#' \item{cn}{Coefficients used in the computation of \ifelse{latex}{\eqn{P_{n}}: \eqn{C_n}}{\out{P<sub>n</sub>: <i>C<sub>n</sub></i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_q}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_q}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @examples
#' #A bank has 5 ATMs. Occasionally one ot them is 
#' #damaged until one of the two hired technicians
#' #repairs it. It is known that the mean time to repair
#' #each ATM follows an exponential distribution with mean
#' #10 minutes, while the distribution of time an ATM
#' #works is also exponential
#' #with mean 2 hours. The bank has an ATM extra to
#' #replace a damaged one.
#' 
#' M_M_S_INF_H_Y(lambda=1/2, mu=60/10, s=2, h=5, y=1)
#' @export
#' @family AnaliticalModels
M_M_S_INF_H_Y <- function(lambda=3, mu=6, s=3, h=5, y=3) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  if (!is.numeric(s)      | s <= 0) stop ("Argument 's' must be greather than zero")
  if (!is.numeric(h)      | h <= 0) stop("Argument 'h' must be greather than zero")
  if (!is.numeric(y)      | y <= 0) stop("Argument 'y' must be greather than zero")
  
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

#' @describeIn Pn Implements the method for a M/M/s/\eqn{\infty}/H/Y queueing model
#' @method Pn M_M_S_INF_H_Y
#' @usage NULL
#' @export
Pn.M_M_S_INF_H_Y <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop(paste("P(n): Index out of limits: 0:Inf\n", sep=""))}
  
  ifelse(n > (qm$y + qm$h), 0, {
    pn <- c(qm$out$p0, qm$out$cn[-1]*qm$out$p0)
    pn[n+1] 
  })          
}

#' @describeIn maxCustomers Implements the method for a M/M/s/\eqn{\infty}/H/Y queueing model
#' @method maxCustomers M_M_S_INF_H_Y
#' @usage NULL
#' @export
maxCustomers.M_M_S_INF_H_Y<- function(qm) {
  return(qm$y + qm$h)
}


#' @describeIn Qn Implements the method for a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' @method Qn M_M_S_INF_H_Y
#' @usage NULL
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

#' @describeIn FWq Implements the method for a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' @method FWq M_M_S_INF_H_Y
#' @usage NULL
#' @export
FWq.M_M_S_INF_H_Y <- function(qm, x) {
  return(FWq.M_M_S_INF_H(qm,x))
}

#' @describeIn FW Implements the method for a M/M/s/\eqn{\infty}/H with Y replacements queueing model
#' @method FW M_M_S_INF_H_Y
#' @usage NULL
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
#' \item{rho}{Constant coefficient: \eqn{\lambda/\mu}}
#' \item{barrho}{Traffic intensity: \ifelse{latex}{\eqn{\bar{\rho}}}{\out{<i>&#862;&rho;</i>}}}
#' \item{p0}{Probability of empty system: \ifelse{latex}{\eqn{P_{0}}}{\out{<i>P<sub>0</sub></i>}}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_q} (\eqn{L_q} = 0 in this model)}{\out{<i>L<sub>q</sub> (L<sub>q</sub> = 0 in this model)</i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_q} (\eqn{W_q} = 0 in this model)}{\out{<i>W<sub>q</sub> (W<sub>q</sub> = 0 in this model)</i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @examples
#' #The number of people turning on their television sets
#' #on Saturday evening during prime time can be described
#' #rather well by a Poisson distribution with a mean of
#' #100000/hr. 
#' #There are five major TV stations, and a given person
#' #chooses among them essentially at random.
#' #Surveys have also shown that the average person tunes
#' #in for 90 min and that viewing times are approximately
#' #exponentially distributed.
#' M_M_INF(lambda=100000/5, mu=60/90)
#' @export 
#' @family AnaliticalModels
M_M_INF <- function(lambda=3, mu=6) {
  if (!is.numeric(lambda) | lambda <= 0) stop("Argument 'lambda' must be greather than zero")
  if (!is.numeric(mu)     | mu <= 0) stop("Argument 'mu' must be greather than zero")
  
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


#' @describeIn Pn Implements the method for a M_M_INF queueing model
#' @method Pn M_M_INF
#' @usage NULL
#' @export
Pn.M_M_INF <- function(qm, n) {
  #Comprobamos que n sea entero
  if (!all.equal(n, as.integer(n))) stop("P(n): Argument 'n' must be an integer")
  if (any(is.na(n))) stop("P(n): Argument 'n' invalid")
  
  minval <- min(n)
  maxval <- max(n)
  if (minval < 0) {stop("P(n): Index out of limits: 0:Inf\n")}
  
  lambda <- rate(qm$arrivalDistribution)
  mu <- rate(qm$serviceDistribution)   
  if (maxval > 0) {
    cn <- c(1, lambda/((1:maxval)*mu))
    cn <- cumprod(cn)
  } else
    cn <- c(1)
  if (any(is.infinite(cn[n+1]))) warning("Possible overflow. Infinite detected. Try to change the parameters.")
  return(cn[n+1]*qm$out$p0)
}

#' @describeIn FW Implements the method for a M/M/\eqn{\infty} queueing model
#' @method FW M_M_INF
#' @usage NULL
#' @export
FW.M_M_INF <- function(qm, x) {
  return(p(qm$serviceDistribution)(x))
}

#' @describeIn FWq Implements the method for a M/M/\eqn{\infty} queueing model
#' @method FWq M_M_INF
#' @usage NULL
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