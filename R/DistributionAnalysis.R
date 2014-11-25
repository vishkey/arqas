#' Computes the estimated parameters of the giving distributions for the input data
#' 
#' @param data data to estimate parameters
#' @param ldistr A list of distributions names for estimate its parameters
#' @return A list of estimate parameters for each distribution
#' @export
#' @family DistributionAnalysis
fitData <- function (data, ldistr= c("exp", "norm", "weibull", "unif", "lnorm", "gamma")) {
  if (is.null(data)) stop("Argument 'data' must be an array of numeric values")
  options(warn=-1)
  res <-lapply(ldistr, function(x) {tryCatch(fitdist(data, x, method="mle", start=NULL),
                                             error = function(e) {
                                               print(e)
                                               return(NA)
                                             })})
  names(res) <- ldistr
  res <- res[!is.na(res)]
  class(res) <- c("FitList", class(res))
  return(res)
}


#' Computes the p-value of the chi-square test and Kolmogorov-Smirnov test and each statistic value
#' 
#' @param lfitdata a list of fitted data
#' @return the p-values and the statistics values of the chi-square test and Kolmogorov-Smirnov test
#' @export
#' @family DistributionAnalysis
goodnessFit <- function(lfitdata) {
  res <- data.frame(distrnames=NULL, chisq= NULL, chisq.pvalue=NULL, ks= NULL, ks.pvalue=NULL)
  
  for(i in lfitdata){
      callfun <- paste("ks.test(i$data, 'p", i$distname, "', ", sep="")
      argnames <- names(i$estimate)
      j <- 0
      if (length(i$estimate) > 1){
        for(j in 1:(length(i$estimate)-1)) {
          callfun <- paste(callfun, argnames[j], "=", i$estimate[j], ", ", sep="")
        }
      }
      callfun <- paste(callfun, argnames[j+1], "=", i$estimate[j+1], ")", sep="")
      
      testks <- eval(parse(text=callfun))
      argschisq <- with(i, c(list(x=data, distribution = distname), estimate) )
      dataChisq <- do.call("chisq.test.cont", args = argschisq)
      res <- rbind(res, data.frame(distrnames=i$distname, chisq=dataChisq$statistic, chisq.pvalue=format(round(dataChisq$p.value, 3), nsmall=3), ks=testks$statistic, ks.pvalue=format(round(testks$p.value, 3), nsmall=3)))
  }
  class(res) <- c("GoodnessFit", class(res))
  return(res)
}

#' Shows three plots:
#'      The histogram and theoretical densities
#'      The empirical and thoeretical CDF's
#'      The Q-Q plot
#' 
#' @param lfitdata a list of fitted data 
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @export
#' @family DistributionAnalysis
summaryFit <- function(lfitdata, graphics="ggplot2") {
  ltext <- names(lfitdata)
  if (ltext[1] == "estimate")
    ltext <- lfitdata$distname
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  switch (graphics, 
          "graphics" = {denscomp(lfitdata, legendtext=ltext)
                        cdfcomp(lfitdata, legendtext=ltext)
                        qqcomp(lfitdata,  legendtext=ltext)},
          "ggplot2" = {grid.arrange(denscompggplot2(lfitdata),
                                 cdfcompggplot2(lfitdata),
                                 qqcompggplot2(lfitdata), ncol=2
                       )}
          )
}

#' Density histogram Plot using the package ggplot2
#' 
#' @param lfitdata a list of data fitted
#' @export
#' @family DistributionAnalysis
denscompggplot2 <- function(lfitdata){
   ltext <- names(lfitdata)
   if (ltext[1] == "estimate") {
     n <- lfitdata$n
     y <- lfitdata$data
   }
   else {
     n <- lfitdata[[1]]$n
     y <- lfitdata[[1]]$data
   }
   aux <- hist(y, plot=FALSE)
   data <- data.frame(data=aux$mids, Density=aux$density)
   resplot <- ggplot(data, aes(x=data, y=Density)) + geom_histogram(stat="identity", binwidh = 1)
   x <- seq(aux$breaks[1], aux$breaks[length(aux$breaks)], length.out=n)
   for(fit in lfitdata) {
      eldata <- do.call(paste("d", fit$distname, sep=""), c(list(x), as.list(fit$estimate)))
      aux2 <- data.frame(data=x, Density=eldata, Distributions=fit$distname)
      Density <- Distributions <- NULL
      resplot <- resplot + geom_line(data=aux2, size=0.25, aes(x=data, y=Density, colour=Distributions))
   }
  resplot
}

#' Cumulative Density plot using the package ggplot2
#' 
#' @param lfitdata a list of fitted data 
#' @export
#' @family DistributionAnalysis
cdfcompggplot2 <- function(lfitdata){
  ltext <- names(lfitdata)
  if (ltext[1] == "estimate") {
    n <- lfitdata$n
    y <- lfitdata$data
  }
  else {
    n <- lfitdata[[1]]$n
    y <- lfitdata[[1]]$data
  }
  cumulative <- cumsum(rep(1/n, n))
  y <- sort(y)
  data <- data.frame(data=y, Density=cumulative)
  resplot <- ggplot(data, aes(x=data, y=Density)) + geom_line()
  x <- seq(y[1], y[n], length.out=n)
  for(fit in lfitdata) {
    eldata <- do.call(paste("r", fit$distname, sep=""), c(list(x), as.list(fit$estimate)))
    eldata <- sort(eldata)
    aux2 <- data.frame(data=eldata, Density=cumulative, Distributions=fit$distname)
    Density <- Distributions <- NULL
    resplot <- resplot + geom_line(data=aux2, size=0.25, aes(x=data, y=Density, colour=Distributions))
  }
  resplot
}

#' Q-Q Plot using the package ggplot2
#' 
#' @param lfitdata a list of fitted data 
#' @export
#' @family DistributionAnalysis
qqcompggplot2 <- function(lfitdata) {
  ltext <- names(lfitdata)
  if (ltext[1] == "estimate") {
    n <- lfitdata$n
    y <- lfitdata$data
  }
  else {
    n <- lfitdata[[1]]$n
    y <- lfitdata[[1]]$data
  }
  data <- data.frame(data=y, Density=y)
  resplot <- ggplot(data, aes(x=data, y=Density)) + geom_line()
  y <- sort(y)
  x <- cumsum(rep(1/n, n))
  for(fit in lfitdata) {
    eldata <- do.call(paste("q", fit$distname, sep=""), c(list(x), as.list(fit$estimate)))
    aux2 <- data.frame(data=y, Density=eldata, Distributions=fit$distname)
    Density <- Distributions <- NULL
    resplot <- resplot + geom_line(data=aux2, size=0.25, aes(x=data, y=Density, colour=Distributions))
  }
  resplot
}


#-------------------------------------------------------------------------------
# chisq.test.cont(x, distribution, nclasses, output, nestpar,...)
#-------------------------------------------------------------------------------
# Realiza el test ji-cuadrado de bondad de ajuste para una distribución continua
# discretizando en intervalos equiprobables.
# Parámetros:
#   distribution = "norm","unif",etc
#   nclasses = floor(length(x)/5)
#   output = TRUE
#   nestpar = 0= nº de parámetros estimados
#   ... = parámetros distribución
# Ejemplo:
#   chisq.test.cont(x, distribution="norm", nestpar=2, mean=mean(x), sd=sqrt((nx-1)/nx)*sd(x))
#-------------------------------------------------------------------------------
chisq.test.cont <- function(x, distribution = "norm", nclasses = min(100, max(3, floor(length(x)/5))), 
                            output = FALSE, nestpar = 0, ...) {
  # Funciones distribución
  q.distrib <- eval(parse(text = paste("q", distribution, sep = "")))
  d.distrib <- eval(parse(text = paste("d", distribution, sep = "")))
  
  # Puntos de corte
  q <- q.distrib((1:(nclasses - 1))/nclasses, ...)
  tol <- sqrt(.Machine$double.eps)
  xbreaks <- c(min(x) - tol, q, max(x) + tol)
  
  # Gráficos y frecuencias
  if (output) {
    xhist <- hist(x, breaks = xbreaks, freq = FALSE, lty = 2, border = "grey50")
    curve(d.distrib(x, ...), add = TRUE)
  } else {
    xhist <- hist(x, breaks = xbreaks, plot = FALSE)
  }
  
  # Cálculo estadístico y p-valor
  O <- xhist$counts  # Equivalente a table(cut(x, xbreaks)) pero más eficiente
  E <- length(x)/nclasses
  DNAME <- deparse(substitute(x))
  METHOD <- "Pearson's Chi-squared test"
  STATISTIC <- sum((O - E)^2/E)
  names(STATISTIC) <- "X-squared"
  PARAMETER <- nclasses - nestpar - 1
  names(PARAMETER) <- "df"
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  
  # Preparar resultados
  classes <- format(xbreaks)
  classes <- paste("(", classes[-(nclasses + 1)], ",", classes[-1], "]", 
                   sep = "")
  RESULTS <- list(classes = classes, observed = O, expected = E, residuals = (O - 
                                                                                E)/sqrt(E))
  if (output) {
    cat("\nPearson's Chi-squared test table\n")
    print(as.data.frame(RESULTS))
  }
  if (any(E < 5)) 
    warning("Chi-squared approximation may be incorrect")
  structure(c(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, 
                   method = METHOD, data.name = DNAME), RESULTS), class = "htest")
}
