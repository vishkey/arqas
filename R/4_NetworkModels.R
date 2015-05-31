#' Ciclic Routing Matrix
#' 
#' @param n Number of nodes
#' @keywords internal
Ciclic <- function(n) {
  diag(n)[c(2:n,1),]
}

#' Serial Routing Matrix
#' 
#' @param n Number of nodes
#' @keywords internal
Tamdem <- function(n) {
  diag(n+1)[2:(n+1), 1:n]
}

#' Obtains the main characteristics of an Open Jackson network model
#' 
#' @param lambda Vector of arrival rates at each node
#' @param mu Vector of mean service rates
#' @param s Vector with the number of servers at each node
#' @param p Routing matrix, where \ifelse{latex}{\eqn{p_{ij}}}{\out{p<sub>ij</sub>}} is the routing probability from node i to node j
#' @return
#' Returns the next information of an Open Jackson network model:
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{l}{Vector with the number of customers in the nodes: \eqn{L}}
#' \item{lq}{Vector with the number of customers in the queue at each node: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Vector with the waiting time in each node: \eqn{W}}
#' \item{wq}{Vector with the waiting time in the queue at each node: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{lt}{Number of customers in the network: \ifelse{latex}{\eqn{L_{Total}}}{\out{<i>L<sub>Total</sub></i>}}}
#' \item{lqt}{Number of customers in all the queues: \ifelse{latex}{\eqn{L_{qTotal}}}{\out{<i>L<sub>qTotal</sub><i>}}}
#' \item{wt}{Total waiting time in the network: \ifelse{latex}{\eqn{W_{Total}}}{\out{<i>W<sub>Total</sub></i>}}}
#' \item{wqt}{Total waiting time in all the queues: \ifelse{latex}{\eqn{W_{qTotal}}}{\out{<i>W<sub>qTotal</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @examples
#' #Two servers recieve 20 tasks for minute the first one,
#' #and 30 tasks for minute the second one.
#' #The unique processsor in the first server can manage
#' #100 tasks for minute, while the two processors in the
#' #second server only can manage 25 task for minute. 
#' #When a task is near to finish in the server 2, it creates
#' #a new task in the server 1 with a probability of 25%,
#' #the task ends in the other case.
#' #The tasks that ends in the server 1 creates a new one
#' #in the same server the 20% of the times and creates
#' #a new one in the server 2 the 10% of the times, ending
#' #in other case.
#' 
#' OpenJacksonNetwork(lambda=c(20, 30), 
#'                    mu=c(100, 25), 
#'                    s=c(1,2), 
#'                    p=matrix(c(0.2,0.1,
#'                               0.25,0), 2, byrow = TRUE)) 
#' @export
#' @family AnaliticalModels
OpenJacksonNetwork <- function(lambda=c(20, 30), mu=c(100, 25), s=c(1,2), p=matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2)) {
  if (!is.numeric(lambda)) stop("Argument 'lambda' must be numeric.")
  if (!is.numeric(mu)) stop("Argument 'mu' must be numeric.")
  if (!is.numeric(s)) stop("Argument 's' must be numeric.")
  if (!is.numeric(p)) stop("Argument 'p' must be numeric.")
  
  #Comprobar parametros de entrada tienen la misma longitud, y P cumple las condiciones  
  sizelambda <- length(lambda)
  sizemu <- length(mu)
  sizes <- length(s)
  sizep <- nrow(p)
  
  if (sizelambda != sizemu || sizelambda != sizes || sizemu != sizes || sizelambda != sizep || sizemu != sizep || sizes != sizep)
    stop("Arguments 'lambda', 'mu','s' and 'p' must have the same length")
  
  if (nrow(p) != sizes || ncol(p) != sizes)
    stop(simpleError(paste("Argument 'p' must have ", sizes, " rows and columns.")))
  
  if (any(p < 0 | p > 1) || any(rowSums(p) > 1)) {
    stop (simpleError("Argument 'p' must have values between 0 and 1 and each row must sum 1 or less."))
  }
  obj <- list(lambda=lambda, mu=mu, servers=s, prob=p)
  
  id <- diag(length(lambda))
  A <- id-t(p)
  
  if (abs(det(A)) <= (.Machine$double.eps ^ 0.5)) {
    stop(simpleError("Not stationary system"))
  }
  # V <- c(solve(A)%*%lambda)
  barlambda <- solve(A, lambda)
  if (any(barlambda > mu*s)) {
    stop(simpleError("Not stationary system"))
  }
  obj$nodes <- mapply(M_M_S, barlambda, mu, s, SIMPLIFY=FALSE)
  lq <- sapply(obj$nodes, function(x) x$out$lq)
  l <- lq+(barlambda/mu)
  w <- lq/barlambda + 1/mu
  wq <- lq/barlambda
  
  lt <- sum(l)
  lqt <- sum(lq)
  wt <- lt/sum(lambda)
  wqt <- lqt/sum(lambda)
  
  obj$A <- A
  obj$out <- list(lq = lq, l=l, w=w, wq=wq, lt=lt, lqt=lqt, wt=wt, wqt=wqt) 
  oldClass(obj) <- c("OpenJackson", "Network", "MarkovianModel")
  return(obj)    
}

#' Data for a Open Network Example
#' 
#' @return The solution of the example
#' @keywords internal
SN_Example <- function () {
  p <- matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2)
  s <- c(1, 2)
  lambda <- c(20 ,30)
  mu <- c(100, 25)
  OpenJacksonNetwork(lambda, mu, s, p)
}

#' Gets the model of the selected node
#' 
#' @param net network
#' @param i node
#' @return object model
#' @method node OpenJackson
#' @keywords internal
#' @export
node.OpenJackson <- function(net, i) {
  if (i <= 0 || i > length(net$mu)) {stop(paste("node: Index out of limits: 1:", length(net$mu), "\n", sep=""))}
  
  return(net$nodes[[i]])
}

#' @describeIn Pn Implements the method for a Open Jackson Network model
#' @method Pn OpenJackson
#' @usage NULL
#' @export
Pn.OpenJackson <- function(qm, n) {
  if (length(n) != length(qm$nodes)) {stop("P(n): Length of indexes must be equal to the number of nodes in the network\n")}
  
  probs <- mapply(Pn, qm$nodes, n)
  return(prod(probs)) 
}

#' @rdname P0i
#' @method P0i OpenJackson
#' @details
#' \code{P0i.OpenJackson} implements the method for an Open Jackson Network model
#' @export
P0i.OpenJackson <- function(net, i) {
  if (i <= 0 || i > length(net$lambda)) {
    stop(paste("P0i: Index out of limits: 1:", length(net$lambda), "\n", sep=""))
  }
  return(net$lambda[i]/sum(net$lambda))
}

#' @rdname Pi0
#' @method Pi0 OpenJackson
#' @details
#' \code{Pi0.OpenJackson} implements the method for an Open Jackson Network model
#' @export
Pi0.OpenJackson <- function(net, i) {
  if (i <= 0 || i > nrow(net$prob)) {
    stop(paste("Pi0: Index out of limits: 1:", nrow(net$prob), "\n", sep=""))
  }
  return(1-sum(net$prob[i,]))
}

#' Computes \ifelse{latex}{\eqn{f_{i}(n)}}{\out{<i>f<sub>i</sub>(n)</i>}}
#' 
#' @param ps partial solution of closed network
#' @param i node
#' @param n clients
#' @return \ifelse{latex}{\eqn{f_{i}(n)}}{\out{<i>f<sub<i</sub>(n)</i>}}
#' @keywords internal
f_close <- function(ps, i, n) {
  rho <- ps$out$rho[i]
  s <- ps$servers[i]
  
  f <- function(n) {
    if (n == 0) {return(1)}
    if (n <= s) {
      return((rho^n)/prod(1:n))
    }
    else {
      return((rho^n)/(prod(1:s)* s^(n-s)))
    }
  }
  return(mapply(f, n))
}

#' Computes matrix G
#' 
#' @param ps partial solution of closed network
#' @return matrix G
#' @keywords internal
calculateG <- function(ps) { 
  res <- matrix(c(1), nrow=ps$k, ncol=(ps$n+1))
  for(i in 2:(ps$n+1)) {
    res[1,i] <- f_close(ps, 1, i-1)
    if (ps$k > 1)
      for(j in 2:ps$k)
        res[j, i] <- sum(res[j-1, i:1] * f_close(ps, j, 0:(i-1)))
  }
  return(res)
  
}

#' Probability of i clients in the last node of the network
#' 
#' @param ps partial solution of closed network
#' @param i clients
#' @return P(i) in the last node
#' @keywords internal
pnlast <- function(ps, i) {
  if (ps$k == 1) return(1)
  if (min(i) < 0 || max(i)>ps$n) stop(paste("Pnlast: Argument 'i' must be between 0 and ", ps$n, sep=""))
  return((f_close(ps, ps$k, i) * ps$g[ps$k-1, ps$n-i+1])/ps$g[ps$k, ps$n+1])
}

#' Data for a Closed Network Example
#' 
#' @return The solution of the example
#' @keywords internal
CN_example <- function() {
  mu <- c(5,5,10,15)
  s <- c(2,2,1,1)
  p <- array(c(0.25,0.15,0.5,0.4,0.15,0.35,0.25,0.3,0.2,0.2,0.15,0.25,0.4,0.30,0.1,0.05), dim=c(4,4))
  nClients <- 3
  ClosedJacksonNetwork(mu, s, p, 3)
}

#' Obtains the main characteristics of a Closed Jackson Network model
#' 
#' @param mu Vector of mean service rates 
#' @param s Vector of servers at each node
#' @param p Routing matrix, where \ifelse{latex}{\eqn{p_{ij}}}{\out{<i>p<sub>ij</sub></i>}} is the routing probability from node i to node j
#' @param n Number of customers in the network
#' @return Returns the next information of a Closed Jackson Network model:
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{l}{Number of customers in the system: \eqn{L}}
#' \item{lq}{Number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Waiting time in the system: \eqn{W}}
#' \item{wq}{Waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_q)}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' @export 
#' @examples
#' #An system have 4 workstations interconnected.
#' #For the control of the system there is three
#' #tasks in continuous execution in some of the
#' #workstations. Once the task ends, this creates
#' #a copy of itself and sends it tu execute in
#' #some of the other three, following the next
#' #probabilities table
#' 
#' # Origin-destiny     1      2     3     4
#' #      1           0.25   0.15  0.20  0.40
#' #      2           0.15   0.35  0.20  0.30
#' #      3           0.50   0.25  0.15  0.10
#' #      4           0.40   0.30  0.25  0.05
#' 
#' #The servers 1 and 2 have two processors and
#' #each of one have a process time with exponential
#' #distribution and capacitiy of 5 tasks for
#' #minute.
#' #The servers 3 and 4 have a single processor
#' #and they can serve 10 and 15 task for minute
#' #respectively.
#' 
#' ClosedJacksonNetwork(mu=c(5,5,10,15),
#'                      s=c(2,2,1,1),
#'                      p=matrix(c(0.25, 0.15, 0.20, 0.40,
#'                                 0.15, 0.35, 0.20, 0.30,
#'                                 0.50, 0.25, 0.15, 0.10,
#'                                 0.40, 0.30, 0.25, 0.05), 4, byrow = TRUE),
#'                      n = 3)
#'                                
#'@family AnaliticalModels
ClosedJacksonNetwork <- function(mu=c(5,5,10,15), s=c(2,2,1,1), p=matrix(c(0.25, 0.15, 0.20, 0.40, 0.15, 0.35, 0.20, 0.30,0.50, 0.25, 0.15, 0.10,0.40, 0.30, 0.25, 0.05), 4, byrow = TRUE), n=3) {
  if (!is.numeric(mu)) stop("Argument 'mu' must be numeric.")
  if (!is.numeric(s)) stop("Argument 's' must be numeric.")
  if (!is.numeric(p)) stop("Argument 'p' must be numeric.")
  if (!is.numeric(n)) stop("Argument 'n' must be numeric.")
  
  sizemu <- length(mu)
  sizes <- length(s)
  sizep <- if(is.null(output <- nrow(p))) length(p) else output
  
  if (sizemu != sizes || sizemu != sizep || sizes != sizep)
    stop("Arguments 'mu', 's' and 'p' must have the same length")
  if (n < 0)
    stop("Argument 'n' must be greather than 0")
  
  if(!is.null(nrow(p)))
      if (any(p < 0 | p > 1) || any(rowSums(p) != 1)) {
        stop (simpleError("Argument 'p' must have values between 0 and 1 and each row must sum 1."))
      }
  obj <- list(mu=mu, servers=s, prob=p, n=n)
  obj$k <- k <- length(mu)
  if (!is.null(numCol <- ncol(p)) && length(p) > 2) {
     id <- diag(k)
     trasp <- t(p)
     id <- id[-nrow(id),]
     aux <- id - trasp[-nrow(trasp),]
     A <- aux[,-1]
     B <- -aux[,1]
     obj$lambda <- lambda <- c(1, solve(A, B))
     print(lambda)
  } else if (length(p) == 2) {
     obj$lambda <- lambda <- c(1, trasp[2,1]/trasp[2,2])  
  } else {
     obj$lambda <- lambda <- c(1)
  }
  
  rho <- lambda/mu
  obj$out$rho <- rho
  nodes <- matrix(c(1), nrow=k, ncol=4, dimnames=list(c(1:k), c("L", "Lq", "W","Wq")))

  shiftdown <- c(k, 1:(k-1))
  obj$out$gkn <- NULL
  for (node in k:1) { 
    #Calculamos G
    obj$g <- calculateG(obj)
    if (is.null(obj$out$gkn)) obj$out$gkn <- obj$g[k, n+1]
    l <- sum(1:n * pnlast(obj, 1:n))
    if (node == k) {
      barlambda <- 0
      for (i in 1:n) {
        if (i <= obj$servers[k]) {
          barlambda <- barlambda + obj$mu[k]*i*pnlast(obj, i)
        } else {
          barlambda <- barlambda + obj$mu[k]*obj$servers[k]*pnlast(obj, i)
        }
      }
      c <- barlambda/obj$lambda[k]
    } else {
      barlambda <- c*obj$lambda[node]
    }
    w <- l/barlambda
    wq <- w - (1/obj$mu[k])
    lq <- l/w*wq
    nodes[node,] <- c(l, lq, w, wq)
    
    #Desplazamos mu, rho, s Y P
    if (obj$k > 1) {
      obj$mu <- obj$mu[shiftdown]
      obj$out$rho <- obj$out$rho[shiftdown]
      obj$servers <- obj$servers[shiftdown]
      obj$prob <- obj$prob[shiftdown, shiftdown]
    }
  }
  #obj$out$nodes <- nodes
  obj$out$l <- as.numeric(nodes[,"L"])
  obj$out$lq <- as.numeric(nodes[, "Lq"])
  obj$out$w <- as.numeric(nodes[, "W"])
  obj$out$wq <- as.numeric(nodes[, "Wq"])
  oldClass(obj) <- c("ClosedJackson", "Network", "MarkovianModel")
  return(obj)
}

#' @describeIn Pn Implements the method for a Closed Jackson Network model
#' @method Pn ClosedJackson
#' @usage NULL
#' @export
Pn.ClosedJackson <- function(qm, n) {
  if (length(n) != qm$k) stop("P(n): Argument 'n' must be a list of the same length than the number of nodes in the network")
  if (sum(n) != qm$n) stop("P(n): the sum of all values of n must be equal to the clients in the network")
  acum <- 1
  for(i in 1:(qm$k)) {
    acum <- acum * f_close(qm, i, n[i])
  }
  probs <- (1/qm$out$gkn)* acum
  return(probs)
}

#' @rdname Pi
#' @method Pi ClosedJackson
#' @details
#' \code{Pi.ClosedJackson} implements the method for a Closed Jackson Network model
#' @export
Pi.ClosedJackson <- function(net, n, node) {
  if ((node < 1) || (node > net$k)) stop(paste("Pi: node must be in the range 1:", net$k, sep=""))
  
  #Calculammos el desplazamiento
  if ((node+1) <= net$k) {
    shift <- c((node+1):net$k ,1:node)
    #Desplazamos los valores necesarios de la solucion parcial
    net$out$rho <- net$out$rho[shift]
    net$servers <- net$servers[shift]
  } 
  g <- calculateG(net) 
  return(c((f_close(net, net$k, n[n<=net$n]) * g[net$k-1, net$n-n[n<=net$n]+1])/g[net$k, net$n+1], rep(0, length(n[n>net$n]))))      
}

#' Print the main characteristics of a Open Jackson Network model
#' @param x a OpenJackson object
#' @param ... Further arguments passed to or from other methods.
#' @method print OpenJackson
#' @keywords internal
#' @export
print.OpenJackson <- function(x, ...) {
  cat("Model: ", class(x)[1], "\n")
  aux <- matrix(c(x$out$l, x$out$lt, x$out$lq, x$out$lqt, x$out$w, x$out$wt, x$out$wq, x$out$wqt), ncol=4, dimnames=list(c(as.character(1:length(x$s)), "Total"), c("L", "Lq", "W", "Wq")))
  print(aux)
}

#' Print the main characteristics of a Closed Jackson Network model
#' @param x a ClosedJackson object
#' @param ... Further arguments passed to or from other methods.
#' @method print ClosedJackson
#' @keywords internal
#' @export
print.ClosedJackson <- function(x, ...) {
  cat("Model: ", class(x)[1], "\n")
  aux <- matrix(c(x$out$l, x$out$lq, x$out$w, x$out$wq), ncol=4, dimnames=list(as.character(1:x$k), c("L", "Lq", "W", "Wq")))
  print(aux)
}