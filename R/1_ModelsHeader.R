#' Define a queueing model
#'
#' Constructor for \code{MarkovianModel} class.
#'
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @return 
#' An object of class \code{MarkovianModel}, a list with the following components:
#' \item{arrivalDistribution}{Arrival distribution}
#' \item{serviceDistribution}{Service distribution}
#' @export
MarkovianModel <- function (arrivalDistribution = Exp(1), serviceDistribution = Exp(1)) {
  packagearr <- attr(class(arrivalDistribution), "package")
  packageserv <- attr(class(serviceDistribution), "package")
  if (is.null(packagearr) || packagearr != "distr") stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
  if (is.null(packageserv)|| packageserv!= "distr") stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
  
  obj <- list(arrivalDistribution=arrivalDistribution, serviceDistribution=serviceDistribution)
  oldClass(obj) <- "MarkovianModel"
  return(obj)
}

#' Distribution function of the waiting time in the system
#' 
#' Returns the value of the cumulative distribution function of the waiting time in the system
#' for a queueing model
#' \deqn{W(x) = P(W \le x)}
#' 
#' @param qm Queueing model
#' @param x Time
#' @return \deqn{W(x)}
#' @export
FW <- function(qm, x) {UseMethod("FW", qm)}

#'  Distribution function of the waiting time in the queue
#'  
#'  Returns the value of the cumulative distribution function of waiting time in the queue
#'  \deqn{W_{q} = P(W_{q} \le x)}
#'  
#'  @param qm Queueing model
#'  @param x Time
#'  @return \deqn{W_{q}(x)}
#'  @export
FWq <- function(qm, x) {UseMethod("FWq", qm)}

#' Steady-state probability of having n customers in the system
#' 
#' Returns the probability of having n customers in the given queueing model
#' 
#' @param qm Queueing model
#' @param n Customers
#' @return \deqn{P(n)}
#' @export
Pn <- function(qm, n) {UseMethod("Pn", qm)}

#' Steady-state probability of finding n customers in the system when a new customer arrives
#' 
#' Returns the probability of n customers in the system in the moment of the
#' arrival of a customer.
#' 
#' @param qm Queueing model
#' @param n Customers
#' @return \deqn{Q(n)}
#' @export
Qn <- function(qm, n) {UseMethod("Qn", qm)}

#' Steady-state probability of 0 customers in the system on the node i of an Open Jackson Network.
#' 
#' Returns the value of the probability of 0 customers in node i of an Open Jackson Network.
#' 
#' @param net Network
#' @param i Node
#' @return \deqn{P_{0,i}()}
#' @export
P0i <- function(net, i) {UseMethod("P0i", net)}

#' Steady-state probability of i customers in the system on the node 0 of an Open Jackson Network.
#' 
#' Returns the value of the probability of i customers in node 0 of an Open Jackson Network.
#' 
#' @param net Network
#' @param i Customers
#' @return \eqn{P_{i,0}()}
#' @keywords internal
Pi0 <- function(net, i) {UseMethod("Pi0", net)}

#' Returns the queueing model which corresponds to the node i of the network
#' 
#' @param net Network
#' @param i Node
#' @return \code{MarkovianModel} object
#' @export
node <- function(net, i) {UseMethod("node", net)}

#' Steady-state probability of n customers in the node i of a network.
#' 
#' Returns the value \eqn{P_{i}(n)} in the node i of a network
#' 
#' @param net Network
#' @param n Customers
#' @param node Node
#' @return P(n) in the selected node
#' @export
Pi <- function(net, n, node) {UseMethod("Pi", net)}


#' @rdname FW
#' @method FW MarkovianModel
#' @details
#' \code{FW.MarkovianModel} implements the default method (generates a message)
#' @export
FW.MarkovianModel <- function(qm, x) {stop(simpleError("W(t): Model not defined"))}

#' @rdname FWq
#' @method FWq MarkovianModel
#' @details
#' \code{FWq.MarkovianModel} implements the default method (generates a message)
#' @export
FWq.MarkovianModel <- function(qm, x) {stop(simpleError("Wq(t): Model not defined"))}

#' @rdname Pn
#' @method Pn MarkovianModel
#' @details
#' \code{Pn.MarkovianModel} implements the default method (generates a message)
#' @export
Pn.MarkovianModel <- function(qm, n) {stop(simpleError("Pn(t): Model not defined"))}

#' @rdname Qn
#' @method Qn MarkovianModel
#' @details
#' \code{Qn.MarkovianModel} implements the default method (generates a message).
#' @export
Qn.MarkovianModel <- function(qm, n) {stop(simpleError("Qn(t): Model not defined"))}

#' Print the main characteristics of a queueing model
#' @param x MarkovianModel object
#' @param ... Further arguments passed to or from other methods.
#' @method print MarkovianModel
#' @keywords internal
#' @export
print.MarkovianModel <- function(x, ...) {
  cat("Model: ", class(x)[1])
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$rho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
}

#' Shows the main graphics of the parameters of a Markovian Model
#' 
#' @param object Markovian Model
#' @param t Range of t
#' @param n Range of n
#' @param ... Further arguments passed to or from other methods.
#' @method summary MarkovianModel
#' @export
summary.MarkovianModel <- function(object, t=list(range=seq(object$out$w, object$out$w*3, length.out=100)), n=c(0:5), ...) {
  if (!is.null(t) && !is.null(n)) {
      par(mfrow=c(2,1)) 
      summaryWtWqt(object, t, "graphics")
      summaryPnQn(object, n, "graphics")
  } else {
    if (!is.null(t)) {
        summaryWtWqt(object, t)
    }
    if (!is.null(n))  {
        summaryPnQn(object, n)
    }
  }
}

#' Shows a plot of W(t) and Wq(t) values of a Markovian Model
#' 
#' @param object Markovian Model
#' @param t Range of t
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @export
summaryWtWqt <- function(object, t, graphics="ggplot2") {
  try({
    epsilon <- 0.001
    
    if (is.list(t)) {
      searchvalues <- FW(object, t$range)
      closeone <- t$range[which.min(1-searchvalues)]
      t <- seq(0, closeone, 0.01)
    }
    data <- data.frame(t, "W"=FW(object, t), "Wq"=FWq(object, t))
    switch(graphics,
           "graphics" = {                  
                  plot(data$t, data$W, col="red", type="l", ylim=c(0,1),  xlab="t", ylab="W(t) & Wq(t)", ann=FALSE)
                  lines(data$t, data$Wq, col="blue")
                  legend("bottomright", c("W", "Wq"), lty =c(1,1), col = c("red", "blue"), bty="n")
                  title(main=paste("Distribution functions of waiting times (t from ", data$t[1], " to ", data$t[length(data$t)], ")", sep=""))
           },
           "ggplot2" = {
                  data <- melt(data, id.var="t") 
                  qplot(t, value, data=data, geom="line", colour=variable, 
                        main=paste("Distribution functions of waiting times (t from ", data$t[1], " to ", data$t[length(data$t)], ")", sep=""),
                        ylab="Cumulative Probability") + scale_colour_discrete(name="")
           })
  })
}

#' Shows a plot of P(n) and Q(n) values of a Markovian Model
#' 
#' @param object Markovian Model
#' @param n Range of n
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @export
summaryPnQn <- function(object, n, graphics="ggplot2") {
  switch(graphics,
         "graphics" = {
              tryCatch({
                barplot(rbind(Qn(object, n), Pn(object, n)), names.arg=n, col=c("blue", "red"),  legend.text=c("Qn", "Pn"), beside=TRUE)
                legend(0, 0, legend=c("Qn", "Pn"))
              }, error=function(e) {
                barplot(Pn(object, n), col="red", legend.text="Pn")
              })
              title(main=paste("Probability of n customers in the system (n from ", n[1], " to ", n[length(n)], ")", sep=""))
         },
         "ggplot2" = {
             tryCatch({
               data <- melt(data.frame(n, "Pn"=Pn(object, n), "Qn"=Qn(object, n)), id.var="n")
              qplot(n, value, data=data, geom="bar", stat="identity", fill=variable, position="dodge",
                    main=paste("Probability of n customers (n from ", n[1], " to ", n[length(n)], ")", sep=""),
                    ylab="Probabilities") + scale_fill_discrete(name="")
             }, error= function(e) {
                 data <- melt(data.frame(n, "Pn"=Pn(object, n)), id.var="n")
                 qplot(n, value, data=data, geom="bar", stat="identity", fill=variable,
                       main=paste("Probability of n customers (n from ", n[1], " to ", n[length(n)], ")", sep=""),
                       ylab="Probability") + scale_fill_discrete(name="")
             })
         })
}

#' Returns the maximun value for n that satisfies P(n) > 0
#' 
#' @param qm object MarkovianModel
#' @export
maxCustomers <- function(qm) {UseMethod("maxCustomers", qm)}

#' @rdname maxCustomers
#' @method maxCustomers MarkovianModel
#' @details
#' \code{maxCustomers.MarkovianModel} implements the default method. Returns infinite.
#' @export
maxCustomers.MarkovianModel <- function(qm) {
      return(Inf)
}