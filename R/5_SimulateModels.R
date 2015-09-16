#' Checks if a object belongs to a determinate package
#' 
#' @param obj Object to check
#' @param packagename Name of the package
#' @return TRUE if obj belongs to the package, FALSE if not
#' @keywords internal
belong <- function(obj, packagename) {
  if (length(obj) > 1)
    packageobj <- sapply(lapply(obj, class), attr, "package")
  else
    packageobj <- attr(class(obj), "package")
  
  return(all(packageobj==packagename))
}

#' Combines a list of independent simulations of a queueing model, computing the mean and variance of each characteristic of interest
#' 
#' @param listsims A list of independent simulations
#' @return an object with the mean and estimated precision of estimated parameters L, Lq, W, Wq, Rho and Eff.
#' @examples
#' combineSimulations(G_G_1(nsim=5))
#' @export
combineSimulations <- function(listsims) {
  if (class(listsims)[1] != "list") 
    if (length(intersect(class(listsims), "SimulatedModel")) > 0)
      return(listsims)
    else
      stop(simpleError("The argument 'listsims' must be a list of queue models"))
  else if (length(listsims) == 1) return(listsims[[1]])
  
  getAttribute <- function(qm, attr) {return(getElement(getElement(qm, "out"), attr))}
    
  l <- sapply(listsims, getAttribute, "l")
  lq <- sapply(listsims, getAttribute, "lq")
  w <- sapply(listsims, getAttribute, "w")
  wq <- sapply(listsims, getAttribute, "wq")
  rho <- sapply(listsims, getAttribute, "rho")
  eff <- sapply(listsims, getAttribute, "eff")
  
  res <- listsims[[1]]
  res$out$l <- if(is.null(nrow(l)))
                 list(mean=mean(l), sd=sd(l), summary = summary(l))
               else
                 list(mean=sapply(1:nrow(l), function(i){mean(l[i,])}), sd=sapply(1:nrow(l), function(i){sd(l[i,])}), summary=sapply(1:nrow(l), function(i){summary(l[i,])}))
  res$out$lq <- if(is.null(nrow(lq)))
                   list(mean=mean(lq), sd=sd(lq), summary = summary(lq))
                else
                 list(mean=sapply(1:nrow(lq), function(i){mean(lq[i,])}), sd=sapply(1:nrow(lq), function(i){sd(lq[i,])}), summary=sapply(1:nrow(l), function(i){summary(lq[i,])}))
  res$out$w <- if(is.null(nrow(w)))
                 list(mean=mean(w), sd=sd(w), summary = summary(w))
               else
                 list(mean=sapply(1:nrow(w), function(i){mean(w[i,])}), sd=sapply(1:nrow(w), function(i){sd(w[i,])}), summary=sapply(1:nrow(l), function(i){summary(w[i,])}))
  res$out$wq <- if(is.null(nrow(wq)))
                  list(mean=mean(wq), sd=sd(wq), summary = summary(wq))
                else
                  list(mean=sapply(1:nrow(wq), function(i){mean(wq[i,])}), sd=sapply(1:nrow(wq), function(i){sd(wq[i,])}), summary=sapply(1:nrow(l), function(i){summary(wq[i,])}))
  res$out$rho <- list(mean=mean(rho), sd=sd(rho), summary = summary(rho))
  res$out$eff <- list(mean=mean(eff), sd=sd(eff), summary = summary(eff))
  
  if (is.null(nrow(wq))) {
    sizepn <- function(qm) {return(length(getAttribute(qm, "pn")))}
    maxsizepn <- max(sapply(listsims, sizepn))
    
    getPn <- function(qm) {aux <- getAttribute(qm, "pn")
                           if ((dif <- maxsizepn - length(aux)) > 0)
                             aux <- c(aux, rep(NA, dif))
                           return(aux)}
    pn <- sapply(listsims, getPn)
    res$out$pn <- apply(pn, 1, mean, na.rm = TRUE)
  } else {
    sizepn <- function(qm) {return(dim(getAttribute(qm, "pn"))[1])}
    maxsizepn <- max(sapply(listsims, sizepn))
    
    pn <- array(NA, dim=c(maxsizepn, length(listsims[[1]]$s), length(listsims)))
    for (i in 1:length(listsims)) {
        aux <- getAttribute(listsims[[i]], "pn")
        dims <- dim(aux)
        if ((dif <- maxsizepn - dims[1]) > 0)
          aux <- rbind(aux, matrix(NA, nrow=dif, ncol=dims[2]))
        pn[,,i] <- aux
    }
    res$out$pn <- apply(pn, c(1,2), mean, na.rm=TRUE) 
  }
  return(res)
}

#' Parallelize an queue model function between multiple processors
#' 
#' @param modelfunction Function that implements a simulated queue model
#' @param parameters List of parameters of the simulated queue model
#' @param nsim Number of simulations
#' @param nproc Processors to use in the simulation.
#' @return The main characteristcs of the given model
#' @keywords internal
#' @export
ParallelizeSimulations <- function(modelfunction, parameters, nsim, nproc) {
  if (nproc > 1) {
    cl <- parallel::makeCluster(nproc)
    registerDoParallel(cl, cores=nproc)
    simulations <- NULL
    parallelRes <- foreach(simulations=iterators::idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(do.call(modelfunction, parameters)))
      }
      return(res)
    }
    do.call(parallel::stopCluster, list(cl))
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  } else {
    res <-list()
    for(i in 1:nsim)
      res <- c(res, list(do.call(modelfunction, parameters)))
    return(if (length(res)==1) res[[1]] else res)
  }
}

#' Obtains the main characteristics of a G/G/1 model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1 model:
#' \item{pn}{Stores all the empirical steady-state probabilities of having n customers, with n from 0 to staClients+nClients: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{<i>L, L<sub>q</sub>, W, W<sub>q</sub></i>}},\emph{Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' 
#' G_G_1(Norm(10, 0.5), Unif(5,6), staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_1 <- function(arrivalDistribution = Exp(3), serviceDistribution = Exp(6), staClients = 100, nClients = 1000, historic = FALSE, nsim=10, nproc=1) {
      if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
      if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
      G_G_1_secuential <- function(arrivalDistribution, serviceDistribution, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (!is.numeric(staClients) | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")

        tArr <- r(arrivalDistribution) (staClients+nClients)
        tServ <- r(serviceDistribution) (staClients+nClients)
        if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
        
        iArr <- iServ <- 1
        sysClients <- simClients<- 0
        a <- 0
        b <- -1
        cron <- d <- c <- 0
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, staclients=staClients, nclients=nClients)
        if (historic) hist <- matrix(nrow=(staClients+nClients), ncol=7, dimnames=list(1:(staClients+nClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
        
        while (simClients < staClients) {
          if (sysClients > 0)
            tMin <- min(a, b)
          else
            tMin <- a
  
          cron <- cron + tMin
          
          if (tMin == a) {
            simClients <- simClients + 1     
            if (sysClients == 0) {
              b <- tServ[iServ]
              iServ <- iServ+1
            } else {
              c <- c + tMin*sysClients
              d <- d + tMin*(sysClients-1)
              b <- b - tMin
            }
            a <- tArr[iArr]
            iArr <- iArr + 1
            sysClients <- sysClients+1
          } else {
            c <- c + tMin*sysClients
            d <- d + tMin*(sysClients-1)
            sysClients <- sysClients - 1
            if (sysClients == 0) {
              b <- -1
            } else {
              b <- tServ[iServ]
              iServ <- iServ + 1
            }
            a <- a - tMin
          }
          if (historic) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients, ] <- c(l,lq,w,wq, sysClients, l-lq, cron)
          }
        }
        acumsta <- cron
        simClients <- 0
        cron <- d <- c <- 0
        tnClients <- numeric(nClients)
        while (simClients < nClients) {
          if (sysClients > 0) {
            tMin <- min(a, b)
          }
          else {
            tMin <- a
          }
          cron <- cron + tMin
          tnClients[sysClients+1] <- tnClients[sysClients+1] + tMin
  
          
          if (tMin == a) {
              simClients <- simClients + 1
              if (sysClients == 0) {
                  b <- tServ[iServ]
                  iServ <- iServ+1
              }
              else {
                  c <- c + tMin*sysClients
                  d <- d + tMin*(sysClients-1)
                  b <- b - tMin
              }
              sysClients <- sysClients+1
              a <- tArr[iArr]
              iArr <- iArr + 1
          }
          else {
              c <- c + tMin*sysClients
              d <- d + tMin*(sysClients-1)
              
              sysClients <- sysClients - 1
              if (sysClients == 0) {
                b <- -1
              } else {
                b <- tServ[iServ]
                iServ <- iServ + 1
              }
              a <- a - tMin
          }
          #In each iteration, stores the evolution of the values
          if (historic && simClients > 0) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
    
            hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, l-lq, acumsta+cron)
          }
        }
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        eff <- w/(w-wq)
        rho <- l-lq
        
        minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
        if (length(minprob) > 0)
          pn <- tnClients[1:max(minprob)]/cron
        else
          pn <- c(0, 0)
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_1", "SimulatedModel")
        
        return(obj)
    }  
    
    ParallelizeSimulations(G_G_1_secuential, list(arrivalDistribution, serviceDistribution, staClients, nClients, historic), nsim, nproc)
}

#' Obtains the main characteristics of a G/G/s model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Number of servers
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S model:
#' \item{pn}{vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_S(Norm(10, 0.5), Unif(5,6), 2, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_S <- function (arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=2, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
   if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
   if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
    G_G_S_secuential <- function(arrivalDistribution, serviceDistribution, s, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (!is.numeric(s) | s <= 0) stop("Argument 's' must be greather than 0.")
        if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (!is.numeric(nClients) | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) (nClients+staClients)
        tServ <- r(serviceDistribution) (nClients+staClients)
        if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
        
        iArr <-1
        iServ <- 0
        bussyservs <- numeric()
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
       
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, s=s, staclients=staClients, nclients=nClients)
        tnClients <- numeric(nClients)
        if (historic) hist <- matrix(nrow=(staClients+nClients), ncol=7, dimnames=list(1:(staClients+nClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
        
        while(simClients < staClients) {
          if (sysClients > 0) {
            if (length(bussyservs) > 0)
              tmin <- min(tArr[iArr], min(tServ[bussyservs]))
            else
              tmin <- tArr[iArr]
          }
          else tmin <- tArr[iArr]
          
          if (tmin == tArr[iArr]) {
            simClients <- simClients + 1
            cron <- cron + tmin
            iArr <- iArr + 1    
            c <- c + tmin*sysClients  
            #if exists any server free
            if (sysClients < s) {
              sysClients <- sysClients + 1
              #update the time of the servers working
              if (length(bussyservs) > 0)
                tServ[bussyservs] <- tServ[bussyservs] - tmin
              #Add a new server to the array
              iServ <- iServ+1
              bussyservs <- c(bussyservs, iServ)
            } else {
              d <- d + tmin*(sysClients-s)
              sysClients <- sysClients + 1
              if (length(bussyservs) > 0)
                tServ[bussyservs] <- tServ[bussyservs] - tmin
            }
          }else {
            finishserver <- which.min(tServ[bussyservs])
            cron <- cron + tmin
            c <- c + tmin*sysClients
            
            if (sysClients > s)
              d <- d + tmin*(sysClients-s)
            
            sysClients <- sysClients-1
            if (length(bussyservs) > 0)
              tServ[bussyservs] <- tServ[bussyservs] - tmin
            
            bussyservs <- bussyservs[-finishserver]
            if (sysClients >= s) {
              iServ <- iServ+1
              bussyservs <- c(bussyservs, iServ)
            }
            tArr[iArr] <- tArr[iArr] - tmin
          }
          
          if (historic) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, cron)
          }
        }
        acumsta <- cron
        simClients <- 0
        cron <- d <- c <- 0
        while(simClients < nClients) {
          #print(bussyservs)
          #Check what type of event happen
          if (sysClients > 0) {
            if (length(bussyservs) > 0)
              tmin <- min(tArr[iArr], min(tServ[bussyservs]))
            else
              tmin <- tArr[iArr]
          }
          else tmin <- tArr[iArr]
          
          if (tmin == tArr[iArr]) {
            simClients <- simClients + 1
            cron <- cron + tmin
            tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
            iArr <- iArr + 1
            c <- c + tmin*sysClients     
            #if exists any server free
            if (sysClients < s) {
              sysClients <- sysClients + 1
              #update the time of the servers working
              if (length(bussyservs) > 0)
                tServ[bussyservs] <- tServ[bussyservs] - tmin
              #Add a new server to the array
              iServ <- iServ+1
              bussyservs <- c(bussyservs, iServ)
            } else {
              d <- d + tmin*(sysClients-s)
              sysClients <- sysClients + 1
              if (length(bussyservs) > 0)
                tServ[bussyservs] <- tServ[bussyservs] - tmin
            }
          }else {
             finishserver <- which.min(tServ[bussyservs])
             cron <- cron + tmin
             c <- c + tmin*sysClients
             
             if (sysClients > s)
               d <- d + tmin*(sysClients-s)
             
             tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
             sysClients <- sysClients-1
             if (length(bussyservs) > 0)
                tServ[bussyservs] <- tServ[bussyservs] - tmin
             
             bussyservs <- bussyservs[-finishserver]
             if (sysClients >= s) {
               iServ <- iServ+1
               bussyservs <- c(bussyservs, iServ)
             }
             tArr[iArr] <- tArr[iArr] - tmin
          }
          #In each iteration, stores the evolution of the values
          if (historic && simClients>0) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients+staClients, ] <- c(l,lq,w,wq,sysClients, (l-lq)/s, acumsta+cron)
          }
        }
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        eff <- w/(w-wq)
        rho <- (l-lq)/s
        
        minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
        if (length(minprob) > 0)
          pn <- tnClients[1:max(minprob)]/cron
        else
          pn <- c(0, 0)
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_S", "SimulatedModel")
        
        return(obj)
    }
    ParallelizeSimulations(G_G_S_secuential, list(arrivalDistribution, serviceDistribution, s, staClients, nClients, historic), nsim, nproc)
}

#'Obtains the main characteristics of a G/G/1/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1/K model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @examples
#' G_G_1_K(Norm(10, 0.5), Unif(5,6), 5, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_1_K <- function (arrivalDistribution=Exp(3), serviceDistribution=Exp(6), K=2, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
      if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
      if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
      G_G_1_K_secuential <- function(arrivalDistribution, serviceDistribution, K, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (!is.numeric(K) | K <= 0) stop("Argument 'K' must be greather than 0.")
        if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (!is.numeric(nClients) | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
        
        iArr <- iServ <- 1
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
        tnClients <- numeric(nClients)
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, k=K, staclients=staClients, nclients=nClients)
        if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(nClients+staClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
        
        while(simClients < staClients) {
          
          if (sysClients > 0) {
            tmin <- min(tArr[iArr], tServ[iServ])
          } else {
            tmin <- tArr[iArr]
          }
          
          if(tmin == tArr[iArr]) {
            cron <- cron + tmin
            iArr <- iArr + 1
            
            if (sysClients == (K+1)) {
              tServ[iServ] <- tServ[iServ] - tmin
              c <- c + tmin*sysClients
              d <- d + tmin*(sysClients-1)
              
            } else {
              simClients <- simClients + 1
              if (sysClients == 0) iServ <- iServ + 1
              else {
                c <- c+tmin*sysClients
                d <- d+tmin*(sysClients-1)
                tServ[iServ] <- tServ[iServ] - tmin
              }
    
              sysClients <- sysClients + 1
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            d <- d + tmin*(sysClients-1)
            sysClients <- sysClients - 1
            
            if(sysClients != 0) {
              iServ <- iServ + 1
            }
            tArr[iArr] <- tArr[iArr] - tmin
          }
          
          if (historic) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            hist[simClients, ] <- c(l,lq,w,wq, sysClients, l-lq, cron)
          }
        }
        acumsta <- cron
        simClients <- 0
        cron <- c <- d <- 0
        while(simClients < nClients) {
          
          if (sysClients > 0) {
            tmin <- min(tArr[iArr], tServ[iServ])
          } else {
            tmin <- tArr[iArr]
          }
          
          
          if(tmin == tArr[iArr]) { 
            cron <- cron + tmin
            tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
            iArr <- iArr + 1
            
            if (sysClients == (K+1)) {
              tServ[iServ] <- tServ[iServ] - tmin
              c <- c + tmin*sysClients
              d <- d + tmin*(sysClients-1)
             
            } else {
              simClients <- simClients + 1
              if (sysClients == 0) {
                iServ <- iServ + 1
              } else {
                c <- c+tmin*sysClients
                d <- d+tmin*(sysClients-1)
                tServ[iServ] <- tServ[iServ] - tmin
              }
              sysClients <- sysClients + 1
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            d <- d + tmin*(sysClients-1)
            tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
            sysClients <- sysClients - 1
            
            if(sysClients != 0) {
              iServ <- iServ + 1
            }
            tArr[iArr] <- tArr[iArr] - tmin
          }
          #In each iteration, stores the evolution of the values
          if (historic && simClients > 0) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            hist[(simClients+staClients), ] <- c(l,lq,w,wq, sysClients, l-lq, acumsta+cron)
          }
        }
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        rho <- l-lq
        eff <- w/(w-wq)
        minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
        if (length(minprob) > 0)
          pn <- tnClients[1:max(minprob)]/cron
        else
          pn <- c(0, 0)
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_1_K", "SimulatedModel")
        
        return(obj)
    }
    
    ParallelizeSimulations(G_G_1_K_secuential, list(arrivalDistribution, serviceDistribution, K, staClients, nClients, historic), nsim, nproc)
}

#' Obtains the main characteristics of a G/G/s/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Number of servers
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S/K model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_S_K(Norm(10, 0.5), Unif(5,6), 2, 5, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_S_K <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=2, K=3, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
    G_G_S_K_secuential <- function(arrivalDistribution, serviceDistribution, s, K, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (!is.numeric(s) | s <= 0) stop("Argument 's' must be greather than 0.")
        if (!is.numeric(K) | K <= 0) stop("Argument 'K' must be greather than 0.")
        if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (!is.numeric(nClients)   | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
        
        iArr <- iServ <- 1
        bussyservers <- rep(NA, s)
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, s=s, k=K, staclients=staClients,  nclients=nClients)
        tnClients <- numeric(nClients)
        if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(staClients+nClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
        
        while(simClients < staClients) {
          indice_min_b <- which.min(bussyservers)
          
          if (sysClients > 0) 
            tmin <- min(tArr[iArr], bussyservers[indice_min_b])
          else
            tmin <- tArr[iArr]
          
          if (tmin == tArr[iArr]) {
            cron <- cron + tmin
            iArr <- iArr + 1
            c <- c + tmin*sysClients
            if(sysClients == (K+s)){
              d <- d + tmin*(sysClients-s)
              bussyservers <- bussyservers-tmin
            } else {
              simClients <- simClients + 1
              if (sysClients < s) {
                sysClients <- sysClients + 1
                iArr <- iArr + 1
                bussyservers <- bussyservers - tmin
                bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
                iServ <- iServ + 1
              } else {
                d <- d + tmin*(sysClients-s)
                sysClients <- sysClients + 1
                bussyservers <- bussyservers-tmin
              }
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            if (sysClients > s)
              d <- d + tmin*(sysClients-s)
            
            sysClients <- sysClients - 1
          
            bussyservers <- bussyservers-tmin
            if (sysClients < s)
              bussyservers[indice_min_b] <- NA
            else {
              bussyservers[indice_min_b] <- tServ[iServ]
              iServ <- iServ + 1
            }
            tArr[iArr] <- tArr[iArr] - tmin
          }
          if (historic) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, cron)
          }
        }
        acumsta <- cron
        simClients <- 0
        cron <- c <- d <- 0
        while(simClients < nClients) {
          indice_min_b <- which.min(bussyservers)
          
          if (sysClients > 0) 
            tmin <- min(tArr[iArr], bussyservers[indice_min_b])
          else
            tmin <- tArr[iArr]
          
          if (tmin == tArr[iArr]) {
            cron <- cron + tmin
            tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
            iArr <- iArr + 1
            c <- c + tmin*sysClients
            if(sysClients == (K+s)){
              d <- d + tmin*(sysClients-s)
              bussyservers <- bussyservers-tmin
            } else {
              simClients <- simClients + 1
              if (sysClients < s) {
                sysClients <- sysClients + 1
                iArr <- iArr + 1
                bussyservers <- bussyservers - tmin
                bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
                iServ <- iServ + 1
              } else {
                d <- d + tmin*(sysClients-s)
                sysClients <- sysClients + 1
                bussyservers <- bussyservers-tmin
              }
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            if (sysClients > s)
              d <- d + tmin*(sysClients-s)
            
            tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
            sysClients <- sysClients - 1
            
            bussyservers <- bussyservers-tmin
            if (sysClients < s)
              bussyservers[indice_min_b] <- NA
            else {
              bussyservers[indice_min_b] <- tServ[iServ]
              iServ <- iServ + 1
            }
            tArr[iArr] <- tArr[iArr] - tmin
          }
          #In each iteration, stores the evolution of the values
          if (historic && simClients > 0) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, acumsta+cron)
          }
        }
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        rho <- (l-lq)/s
        eff <- w/(w-wq)
        minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
        if (length(minprob) > 0)
          pn <- tnClients[1:max(minprob)]/cron
        else {
          pn <- c(0, 0)
        }
        
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_S_K", "SimulatedModel")
        
        return(obj)
    }
    ParallelizeSimulations(G_G_S_K_secuential, list(arrivalDistribution, serviceDistribution, s, K, staClients, nClients, historic), nsim, nproc)
}


#' Obtains the main characteristics of a G/G/1/\eqn{\infty}/H model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1/\eqn{\infty}/H model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_1_INF_H(Norm(10, 0.5), Unif(5,6), 10, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_1_INF_H <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), H=5, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
  G_G_1_INF_H_secuential <- function(arrivalDistribution, serviceDistribution, H, staClients, nClients, historic) {
    if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
    if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
    if (!is.numeric(H) | H <= 0) stop("Argument 'H' must be greather than 0.")
    if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
    if (!is.numeric(nClients)   | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
    
    tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
    tServ <- r(serviceDistribution) ((nClients+staClients)*2)
    if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
    if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
    
    possibleClients <- tArr[1:H]
    iServ <- 1
    iArr <- H+1
    sysClients <- 0
    simClients <- cron <- c <- d <- 0
    
    obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, h=H, staclients=staClients, nclients=nClients)
    tnClients <- numeric(nClients)
    if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(nClients+staClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
    
    while(simClients < staClients) {
      
      indice_minimo_a <- which.min(possibleClients)
      if(sysClients>0)
        minimo<-min(possibleClients[indice_minimo_a], tServ[iServ], na.rm=TRUE)
      else
        minimo<-possibleClients[indice_minimo_a]
      if (minimo  == tServ[iServ]) {
        cron <- cron + minimo
        c <- c + minimo*sysClients
        d <- d + minimo*(sysClients-1)
        sysClients <- sysClients-1   
        if (sysClients!=0)
          iServ <- iServ + 1      
        possibleClients<-possibleClients-minimo
        possibleClients[which(is.na(possibleClients))[1]]<- tArr[iArr]
        iArr <- iArr+1
      } else {
        simClients<-simClients+1
        cron <- cron + minimo
        if(sysClients==0) {              
          sysClients<-sysClients+1
          possibleClients <-possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
        } else{
          c <- c + minimo*sysClients
          d <- d + minimo*(sysClients-1)
          sysClients<-sysClients+1
          tServ[iServ]<-tServ[iServ]-minimo
          possibleClients<-possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
        }
      }
      if (historic) {
        l <- c/cron
        lq <- d/cron
        w <- c/simClients
        wq <- d/simClients
        
        hist[simClients, ] <- c(l,lq,w,wq, sysClients, l-lq, cron)
      }
    }
    acumsta <- cron
    simClients <- cron <- c <- d <- 0
    
    while(simClients < nClients) {
      indice_minimo_a <- which.min(possibleClients)
      
      if (sysClients>0)
        minimo <- min(possibleClients[indice_minimo_a],tServ[iServ], na.rm=TRUE)
      else
        minimo <- possibleClients[indice_minimo_a]
      
      if (minimo == tServ[iServ]) {
        cron <- cron + minimo
        c <- c + minimo*sysClients
        d <- d + minimo*(sysClients-1)
        tnClients[sysClients+1] <- tnClients[sysClients+1]+minimo
        sysClients<-sysClients-1
        if (sysClients!=0)
          iServ <- iServ + 1
        possibleClients<-possibleClients-minimo
        possibleClients[which(is.na(possibleClients))[1]]<- tArr[iArr]
        iArr <- iArr + 1
      } else {
        simClients <- simClients+1
        
        if (sysClients==0) {
          cron <- cron + minimo
          tnClients[sysClients+1] <- tnClients[sysClients+1]+minimo
          sysClients <- sysClients+1
          possibleClients <- possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
          iServ <- iServ + 1
        } else {
          cron <- cron+minimo
          c <- c + minimo*sysClients
          d <- d + minimo*(sysClients-1)
          tnClients[sysClients+1] <- tnClients[sysClients+1]+minimo
          sysClients <- sysClients+1
          tServ[iServ] <- tServ[iServ]-minimo
          possibleClients<-possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
        }
      }
      #In each iteration, stores the evolution of the values
      if (historic && simClients > 0) {
        l <- c/cron
        lq <- d/cron
        w <- c/simClients
        wq <- d/simClients
        
        hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, l-lq, acumsta+cron)
      }
    }
    l <- c/cron
    lq <- d/cron
    w <- c/nClients
    wq <- d/nClients
    rho <- l - lq
    eff <- w/(w-wq)
    minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
    if (length(minprob) > 0)
      pn <- tnClients[1:max(minprob)]/cron
    else
      pn <- c(0, 0)
    
    if (historic)
      obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
    else
      obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
    
    oldClass(obj) <-  c("G_G_1_INF_H", "SimulatedModel")
    return(obj)
  }
  ParallelizeSimulations(G_G_1_INF_H_secuential, list(arrivalDistribution, serviceDistribution, H, staClients, nClients, historic), nsim, nproc)
}

#' Obtains the main characteristics of a G/G/S/\eqn{\infty}/H  model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Number of servers
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S/\eqn{\infty}/H model
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system:\eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_S_INF_H(Norm(10, 0.5), Unif(5,6), 3, 10, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels

G_G_S_INF_H <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=3, H=5, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
    G_G_S_INF_H_secuential <- function(arrivalDistribution, serviceDistribution, s, H, staClients, nClients, historic) {
          if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
          if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
          if (!is.numeric(s) | s <= 0) stop("Argument 's' must be greather than 0.")
          if (!is.numeric(H) | H <= 0) stop("Argument 'H' must be greather than 0.")
          if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
          if (!is.numeric(nClients)   | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
          tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
          tServ <- r(serviceDistribution) ((nClients+staClients)*2)
          if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
          if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
          
          possibleClients <- tArr[1:H]
          bussyservers <- rep(NA, s)
          iServ <- 1
          iArr <- H+1
          sysClients <- 0
          simClients <- cron <- c <- d <- 0
          nClicien <- 100/nClients
          obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, s=s, h=H, staclients=staClients, nclients=nClients)
          tnClients <- numeric(nClients)
          if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(nClients+staClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
          
          while(simClients < staClients) {        
            if(sysClients>0) {
              indice_minimo_a <- which.min(possibleClients)
              indice_minimo_b <- which.min(bussyservers)
              tmin <- min(possibleClients[indice_minimo_a], bussyservers[indice_minimo_b])
            }
            else {
              indice_minimo_a <- which.min(possibleClients)
              tmin <- possibleClients[indice_minimo_a]
            }
            
            if (length(indice_minimo_a > 0) && (tmin == possibleClients[indice_minimo_a])) {
              simClients <- simClients + 1
              if (sysClients < s) {
                cron <- cron + tmin
                c <- c + tmin*sysClients
                sysClients <- sysClients + 1
                possibleClients <- possibleClients - tmin
                possibleClients[indice_minimo_a] <- NA
                bussyservers <- bussyservers - tmin
                bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
                iServ <- iServ + 1
              } else {
                cron <- cron + tmin
                c <- c + tmin*sysClients
                d <- d + tmin*(sysClients-s)
                sysClients <- sysClients + 1
                bussyservers <- bussyservers - tmin
                possibleClients <- possibleClients - tmin
                possibleClients[indice_minimo_a] <- NA
              }
            } else {
              cron <- cron + tmin
              c <- c + tmin*sysClients
              if (sysClients > s)
                d <- d + tmin*(sysClients-s)
              sysClients <- sysClients - 1 
              bussyservers <- bussyservers-tmin
              if (sysClients < s)
                bussyservers[indice_minimo_b] <- NA
              else {
                bussyservers[indice_minimo_b] <- tServ[iServ]
                iServ <- iServ + 1
              }
              possibleClients <- possibleClients - tmin
              possibleClients[which(is.na(possibleClients))[1]] <- tArr[iArr]
              iArr <- iArr+1
            }
            if (historic) {
              l <- c/cron
              lq <- d/cron
              w <- c/simClients
              wq <- d/simClients
              
              hist[simClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, cron)
            }
          }
          acumsta <- cron
          simClients <- cron <- c <- d <- 0
#           progressbar <- txtProgressBar(min=0, max=nClients, style=3)
#           auxProg <- !(((0:nClients)/nClients*100) %% 5)
#           jump <- 10
          while (simClients < nClients) {
#             if (auxProg[simClients+1]) setTxtProgressBar(progressbar, simClients)
            if(sysClients>0) {
              indice_minimo_a <- which.min(possibleClients)
              indice_minimo_b <- which.min(bussyservers)
              tmin <- min(possibleClients[indice_minimo_a], bussyservers[indice_minimo_b])
            }
            else {
              indice_minimo_a <- which.min(possibleClients)
              tmin <- possibleClients[indice_minimo_a]
            }
            
            if (length(indice_minimo_a > 0) && (tmin == possibleClients[indice_minimo_a])) {
              simClients <- simClients + 1
              if (sysClients < s) {
                cron <- cron + tmin
                tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
                c <- c + tmin*sysClients
                sysClients <- sysClients + 1
                possibleClients <- possibleClients - tmin
                possibleClients[indice_minimo_a] <- NA
                bussyservers <- bussyservers - tmin
                bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
                iServ <- iServ + 1
              } else {
                cron <- cron + tmin
                tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
                c <- c + tmin*sysClients
                d <- d + tmin*(sysClients-s)
                sysClients <- sysClients + 1
                bussyservers <- bussyservers - tmin
                possibleClients <- possibleClients - tmin
                possibleClients[indice_minimo_a] <- NA
              }
            } else {
              cron <- cron + tmin
              c <- c + tmin*sysClients
              if (sysClients > s)
                d <- d + tmin*(sysClients-s)
              tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
              sysClients <- sysClients - 1 
              bussyservers <- bussyservers - tmin
              if (sysClients < s)
                bussyservers[indice_minimo_b] <- NA
              else {
                bussyservers[indice_minimo_b ] <- tServ[iServ]
                iServ <- iServ + 1
              }
              possibleClients <- possibleClients - tmin
              possibleClients[which(is.na(possibleClients))[1]] <- tArr[iArr]
              iArr <- iArr+1
            } 
            #In each iteration, stores the evolution of the values
            if (historic && simClients > 0) {
              l <- c/cron
              lq <- d/cron
              w <- c/simClients
              wq <- d/simClients
              
              hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, acumsta+cron)
            }
          }
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
          rho <- (l - lq)/s
          eff <- w/(w-wq)
          minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
          if (length(minprob) > 0)
            pn <- tnClients[1:max(minprob)]/cron
          else
            pn <- c(0, 0)
          if (historic)
            obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
          else
            obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
          oldClass(obj) <-  c("G_G_S_INF_H", "SimulatedModel")
          return(obj)
    }
  ParallelizeSimulations(G_G_S_INF_H_secuential, list(arrivalDistribution, serviceDistribution, s, H, staClients, nClients, historic), nsim, nproc)
}

#' Obtains the main characteristics of a G/G/S/\eqn{\infty}/H with Y replacements model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Number of servers
#' @param H Population size
#' @param Y Number of replacements
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1/S/\eqn{\infty}/H/Y model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_S_INF_H_Y(Norm(10, 0.5), Unif(5,6), 3, 10, 2, staClients=10, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels
G_G_S_INF_H_Y <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=3, H=5, Y=3, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
  G_G_S_INF_H_Y_secuential <- function(arrivalDistribution, serviceDistribution, s, H, Y, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution' must be a valid Class of the Distr package")
        if (!is.numeric(s) | s <= 0) stop("Argument 's' must be greather than 0.")
        if (!is.numeric(H) | H <= 0) stop("Argument 'H' must be greather than 0.")
        if (!is.numeric(Y) | Y <= 0) stop("Argument 'Y' must be greather than 0.")
        if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (!is.numeric(nClients)   | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
        
        possibleClients <- tArr[1:H]
        #   bussyservers <- rep(-1, s)
        bussyservers <- rep(NA, s)
        iServ <- 1
        iArr <- H+1
        sysClients <- 0
        simClients <- cron <- c <- d <- 0
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, s=s, h=H, y=Y, staclients=staClients, nclients=nClients)
        tnClients <- numeric(nClients)
        if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(nClients+staClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
        
        while(simClients < staClients) {    
          if(sysClients>0) {
            indice_minimo_b <- which.min(bussyservers)
            indice_minimo_a <- which.min(possibleClients)
            tmin <- min(possibleClients[indice_minimo_a], bussyservers[indice_minimo_b])
          }else {
            indice_minimo_a <- which.min(possibleClients)
            tmin <- possibleClients[indice_minimo_a]
          }
          
          if (length(indice_minimo_a) > 0 && (tmin == possibleClients[indice_minimo_a])) {
            simClients <- simClients + 1
            if (sysClients < s) {
              cron <- cron + tmin
              c <- c + tmin*sysClients
              sysClients <- sysClients + 1
              possibleClients <- possibleClients - tmin
              possibleClients[indice_minimo_a] <- NA
        
              bussyservers <- bussyservers - tmin
              bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
              iServ <- iServ + 1
              if (sysClients <= Y) {
                possibleClients[indice_minimo_a] <- tArr[iArr]
                iArr <- iArr + 1
              }
            } else {
              cron <- cron + tmin
              c <- c + tmin*sysClients
              d <- d + tmin*(sysClients-s)
              sysClients <- sysClients + 1
              bussyservers <- bussyservers - tmin
              possibleClients <- possibleClients - tmin
              possibleClients[indice_minimo_a] <- NA
              if (sysClients <= Y) {
                possibleClients[indice_minimo_a] <- tArr[iArr]
                iArr <- iArr + 1
              }
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            if (sysClients > s)
              d <- d + tmin*(sysClients-s)
            
            sysClients <- sysClients - 1 
            bussyservers <- bussyservers - tmin
            if (sysClients < s)
              bussyservers[indice_minimo_b] <- NA
            else {
              bussyservers[indice_minimo_b ] <- tServ[iServ]
              iServ <- iServ + 1
            }
            possibleClients <- possibleClients - tmin
            if (sysClients >= Y) {
              possibleClients[which(is.na(possibleClients))[1]] <- tArr[iArr]
              iArr <- iArr+1
            }
          }
          if (historic) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, cron)
          }
        }
        acumsta <- cron
        simClients <- cron <- c <- d <- 0
        while (simClients < nClients) {
          
          if(sysClients>0) {
            indice_minimo_b <- which.min(bussyservers)
            indice_minimo_a <- which.min(possibleClients)
            tmin <- min(possibleClients[indice_minimo_a], bussyservers[indice_minimo_b])
          } else {
            indice_minimo_a <- which.min(possibleClients)
            tmin <- possibleClients[indice_minimo_a]
          }
        
          
          if (length(indice_minimo_a) > 0 && (tmin == possibleClients[indice_minimo_a])) {
            simClients <- simClients + 1
            if (sysClients < s) {
              cron <- cron + tmin
              tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
              c <- c + tmin*sysClients
              sysClients <- sysClients + 1
              possibleClients <- possibleClients - tmin
              possibleClients[indice_minimo_a] <- NA
              if (sysClients <= Y) {
                possibleClients[indice_minimo_a] <- tArr[iArr]
                iArr <- iArr + 1
              }
              bussyservers <- bussyservers - tmin
              bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
              iServ <- iServ + 1
            } else {
              cron <- cron + tmin
              tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
              c <- c + tmin*sysClients
              d <- d + tmin*(sysClients-s)
              sysClients <- sysClients + 1
              bussyservers <- bussyservers - tmin
              possibleClients <- possibleClients - tmin
              possibleClients[indice_minimo_a] <- NA
              if (sysClients <= Y) {
                possibleClients[indice_minimo_a] <- tArr[iArr]
                iArr <- iArr + 1
              }
            }
          } else {
            cron <- cron + tmin
            c <- c + tmin*sysClients
            if (sysClients > s)
              d <- d + tmin*(sysClients-s)
            tnClients[sysClients+1] <- tnClients[sysClients+1]+tmin
            sysClients <- sysClients - 1 
            bussyservers <- bussyservers - tmin
            if (sysClients < s)
              bussyservers[indice_minimo_b] <- NA
            else {
              bussyservers[indice_minimo_b ] <- tServ[iServ]
              iServ <- iServ + 1
            }
            possibleClients <- possibleClients - tmin
            if (sysClients >= Y) {
              possibleClients[which(is.na(possibleClients))[1]] <- tArr[iArr]
              iArr <- iArr+1
            }
          } 
          #In each iteration, stores the evolution of the values
          if (historic && simClients > 0) {
            l <- c/cron
            lq <- d/cron
            w <- c/simClients
            wq <- d/simClients
            
            hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, (l-lq)/s, acumsta+cron)
          }
        }
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        rho <- (l - lq)/s
        eff <- w/(w-wq)
        minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
        if (length(minprob) > 0)
          pn <- tnClients[1:max(minprob)]/cron
        else
          pn <- c(0, 0)
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        
        oldClass(obj) <-  c("G_G_S_INF_H_Y", "SimulatedModel")
        return(obj)
  }
  ParallelizeSimulations(G_G_S_INF_H_Y_secuential, list(arrivalDistribution, serviceDistribution, s, H, Y, staClients, nClients, historic), nsim, nproc)
}

#' Obtains the main characteristics of a G/G/\eqn{\infty} model by simulation
#' 
#' @param arrivalDistribution Arrival distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param serviceDistribution Service distribution (object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param staClients Number of customers used in stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/\eqn{\infty} model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Empirical number of customers in the system: \eqn{L}}
#' \item{lq}{Empirical number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Empirical waiting time in the system: \eqn{W}}
#' \item{wq}{Empirical waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @examples
#' G_G_INF(Norm(10, 0.5), Unif(2,4), staClients=50, nClients=100, nsim=10)
#' @export
#' @family SimulatedModels
G_G_INF <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim) || nsim <= 0) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("Argument 'nproc' must be an integer greather than 0")
  
  G_G_INF_secuential <- function(arrivalDistribution, serviceDistribution, staClients, nClients, historic) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (!is.numeric(nClients) | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
      tServ <- r(serviceDistribution) ((nClients+staClients)*2)
      if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
      if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
      
      iServ <- iArr <- 1
      sysClients <- 0
      simClients <- cron <- c <- d <- 0
      a <- 0
      b <- -1
      
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, staclients=staClients, nclients=nClients)
      tnClients <- numeric(nClients)
      if (historic) hist <- matrix(nrow=(nClients+staClients), ncol=7, dimnames=list(1:(nClients+staClients), c("L", "Lq", "W","Wq", "Clients", "Intensity", "tClient")))
      
      while(simClients < staClients) {
        indice_j <- which(b != -1)
        indice_j_menos_uno <- which(b==-1)
        
        if (length(indice_j_menos_uno) == 0) {
          b <- c(b, -1)
          indice_j_menos_uno <- length(b)
        }
        
        if (sysClients > 0) {
          posibles_b <- b[indice_j]
          tmin <- min(a, posibles_b)
        } else
          tmin <- a
        
        cron <- cron + tmin
        if (tmin == a) {
          simClients <- simClients + 1
          c <- c + tmin*sysClients
          sysClients <- sysClients + 1
          a <- tArr[iArr]
          iArr <- iArr + 1
          b[indice_j] <- b[indice_j] - tmin
          b[indice_j_menos_uno[1]] <- tServ[iServ]
          iServ <- iServ + 1
        } else {
          indice_minimo <- which(b==tmin)
          c <- c + tmin*sysClients
          sysClients <- sysClients - 1
          b[indice_j] <- b[indice_j] - tmin
          b[indice_minimo] <- -1
          a <- a - tmin
        }
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/simClients
          wq <- d/simClients
          
          hist[simClients, ] <- c(l,lq,w,wq, sysClients, l-lq, cron)
        }
      }
      acumsta <- cron
      simClients <- cron <- c <- d <- 0
      
      while(simClients < nClients) {
        #print(simClients)
        indice_j <- which(b != -1)
        indice_j_menos_uno <- which(b==-1)
        
        if (length(indice_j_menos_uno) == 0) {
          b <- c(b, -1)
          indice_j_menos_uno <- length(b)
        }
        
        if (sysClients > 0) {
          posibles_b <- b[indice_j]
          tmin <- min(a, posibles_b)
        } else
          tmin <- a
        
        cron <- cron + tmin
        tnClients[sysClients+1] <- tnClients[sysClients+1] + tmin
        if (tmin == a) {
          simClients <- simClients + 1
          c <- c + tmin*sysClients
          sysClients <- sysClients + 1
          a <- tArr[iArr]
          iArr <- iArr + 1
          b[indice_j] <- b[indice_j] - tmin
          b[indice_j_menos_uno[1]] <- tServ[iServ]
          iServ <- iServ + 1
        } else {
          indice_minimo <- which(b==tmin)
          c <- c + tmin*sysClients
          sysClients <- sysClients - 1
          b[indice_j] <- b[indice_j] - tmin
          b[indice_minimo] <- -1
          a <- a - tmin
        }
        #In each iteration, stores the evolution of the values
        if (historic && simClients > 0) {
          l <- c/cron
          lq <- d/cron
          w <- c/simClients
          wq <- d/simClients
          
          hist[simClients+staClients, ] <- c(l,lq,w,wq, sysClients, l-lq, acumsta+cron)
        }
      }
      l <- c/cron
      lq <- d/cron
      w <- c/nClients
      wq <- d/nClients
      rho <- 0
      eff <- w/(w-wq)
      minprob <- which(tnClients>(.Machine$double.eps ^ 0.5))
      if (length(minprob) > 0)
        pn <- tnClients[1:max(minprob)]/cron
      else
        pn <- c(0, 0)
      if (historic)
        obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      else
        obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      oldClass(obj) <-  c("G_G_INF", "SimulatedModel")
      
      return(obj)
  }
  ParallelizeSimulations(G_G_INF_secuential, list(arrivalDistribution, serviceDistribution, staClients, nClients, historic), nsim, nproc)
}

SCN_example <- function (sta, trans) {
  serviceDistribution <- c(Exp(5), Exp(5), Exp(10), Exp(15))
  s <- c(2,2,1,1)
  p <- array(c(0.25,0.15,0.5,0.4,0.15,0.35,0.25,0.3,0.2,0.2,0.15,0.25,0.4,0.30,0.1,0.05), dim=c(4,4))
  nClients <- 3
  ClosedNetwork(serviceDistribution, s, p, sta, nClients, trans)
}

#' Obtains the main characteristics of a Closed Network model by simulation
#' 
#' @param serviceDistribution Service distributions for the nodes of the network (Each element must be an object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Vector of servers at each node
#' @param p Routing matrix, where \ifelse{latex}{\eqn{p_{ij}}}{\out{<i>p<sub>ij</sub></i>}} is the routing probability from node i to node j
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers in the system
#' @param transitions Number of transitions between nodes used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a Closed Network model:
#' \item{pn}{Vector of empirical steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_n}}{\out{<i>P<sub>n</sub></i>}} (Only the probabilities bigger than 0 are included)}
#' \item{l}{Vector of empirical number of customers in the nodes: \eqn{L}}
#' \item{lq}{Vector of empirical number of customers in the queues of the nodes: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{lqt}{Empirical number of customers in the all the queues: \ifelse{latex}{\eqn{L_{qTotal}}}{\out{<i>L_<sub>qTotal</sub></i>}}}
#' \item{w}{Vector of empirical waiting times in the nodes: \eqn{W}}
#' \item{wq}{Vector of empirical waiting times in the queues of the nodes: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{Empirical system efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Empirical traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @examples
#' ClosedNetwork(serviceDistribution = c(Exp(5), Exp(5), Exp(10), Exp(15)),
#'               s                   = c(2,2,1,1),
#'               p                   = matrix(c(0.25, 0.15, 0.2,  0.4,
#'                                              0.15, 0.35, 0.2,  0.3,
#'                                              0.5,  0.25, 0.15, 0.1,
#'                                              0.4,  0.3,  0.25, 0.05), 4, byrow=TRUE),
#'               nClient             = 3,
#'               staClients          = 10,
#'               transitions         = 100,
#'               nsim                = 10)
#' @export
#' @family SimulatedModels

ClosedNetwork <- function(serviceDistribution=c(Exp(5), Exp(5), Exp(10), Exp(15)), s=c(2,2,1,1), p=array(c(0.25,0.15,0.5,0.4,0.15,0.35,0.25,0.3,0.2,0.2,0.15,0.25,0.4,0.30,0.1,0.05), dim=c(4,4)), staClients=100, nClients=3, transitions=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim)) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc)) stop("Argument 'nproc' must be an integer greather than 0")
  
  ClosedNetwork_secuential <- function(serviceDistribution, s, p, staClients, nClients, transitions, historic) {
      if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution' must be a valid Class of the Distr package")
      if (!is.numeric(p))  stop("All elements in argument 'p' must be numerical.")
      if (!is.numeric(s) | any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
      if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (!is.numeric(nClients) | nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      if (!is.numeric(transitions) | transitions < 0) stop("Argument 'transitions' must be equal or greather than 0.")
      
      nodesServ <- length(serviceDistribution)
      nodesS <- length(s)
      nodesP <- nrow(p)
    
      #print(paste(nodesServ, ", ", nodesS, ", ", nodesP))
      if (nodesServ != nodesS || nodesServ != nodesP || nodesS != nodesP) stop("The arguments 'serviceDistribution', 's' and 'p' must have the same dimensions.")
      
      nodes <- length(s)
      maxserv <- max(s)
      map <- function(f) {f(staClients+transitions)}
      tServ <- matrix(sapply(lapply(serviceDistribution, r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=nodes )
      if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
      
      iServ <- rep(1, nodes)
      rand <- r(Unif()) (staClients+transitions)
      simClients <- 0
      c <- d <- cron <- in_node <- rep(0, nodes)
      
      obj <- list(serviceDistribution=serviceDistribution, s= s, prob=p, staclients=staClients, transitions=transitions, nclients=nClients)
      if (historic) hist <- array(dim=c(nodes, 5, staClients+transitions), dimnames=list(1:nodes, c("L", "Lq", "W","Wq", "tClient"), 1:(staClients+transitions)))
      
      sysClients <- floor((nClients)*(s^(-1))/sum(s^-1))
      res <- nClients-sum(sysClients)
      sysClients[1:res] <- sysClients[1:res] + 1
      
      b <- matrix(nrow=maxserv, ncol=nodes)
      for(i in 1:nodes) {
        if (s[i] > sysClients[i]) {
          if (sysClients[i] == 0) next
          b[1:sysClients[i], i] <- tServ[1:sysClients[i], i]
          iServ[i] <- iServ[i] + sysClients[i]
        } else {
          b[1:s[i],i] <- tServ[1:s[i],i]
          iServ[i] <- iServ[i] + s[i]
        }
      }
      while (simClients < staClients) {
        indice_minimo_b <- which.min(b) #al ser los libre NA, devuelve el indice menor del que no es NA
        tmin <- b[indice_minimo_b]
        nodex <- ceiling(indice_minimo_b/maxserv)
        
        simClients <- simClients + 1
        sysClients[nodex] <- sysClients[nodex] - 1
        b <- b - tmin  #los NA siguen siendo NA al restar
        
        if (sysClients[nodex] < s[nodex]) {
          b[indice_minimo_b] <- NA
        } else {
          b[indice_minimo_b] <- tServ[iServ[nodex], nodex]
          iServ[nodex] <- iServ[nodex] + 1
        }
        n_s <- rand[simClients]
        acum <- 0
        node_dest <- 1
        while (acum < n_s) {
          acum <- acum + p[nodex, node_dest]
          node_dest <- node_dest+1
        }
        node_dest <- node_dest-1
        sysClients[node_dest] <- sysClients[node_dest] + 1
        in_node[node_dest] <- in_node[node_dest]+ 1
        
        if (sysClients[node_dest] <= s[node_dest]) {
          indice_minimo_nodo_minimo <- which(is.na(b[,node_dest]))
          b[indice_minimo_nodo_minimo[1], node_dest] <- tServ[iServ[node_dest], node_dest]
          iServ[node_dest] <- iServ[node_dest] + 1
        }
        cron[nodex] <- cron[nodex] + tmin
        for (j in 1:nodes) {
          if ((j == nodex) && (j != node_dest)) {
            c[j] <- c[j]+tmin*(sysClients[j]+1)
            if ((sysClients[j]+1) > s[j])
              d[j] <- d[j]+tmin*(sysClients[j]+1-s[j])
          } else {
            if ((j == node_dest) && (j != nodex)) {
              c[j] <- c[j] + tmin*(sysClients[j]-1)
              if ((sysClients[j]-1) > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
            } else {
              c[j] <- c[j] + tmin*sysClients[j]
              if (sysClients[j] > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-s[j])
            }
          }
        }
        if (historic && all(in_node>0)) {
          l <- c/sum(cron)
          lq <- d/sum(cron)
          w <- c/in_node
          wq <- d/in_node
          hist[,,simClients] <- array(c(l, lq, w, wq, cron), dim= c(nodes, 5))
        }
      }
      simClients <- 0
      cron <- c <- d <- in_node <- rep(0, nodes)
      prob <- matrix(c(0), nrow=(staClients + transitions), ncol=nodes)
      while(simClients < transitions) {    
        indice_minimo_b <- which.min(b) #al ser los libre NA, devuelve el indice menor del que no es NA
        tmin <- b[indice_minimo_b]
        nodex <- ceiling(indice_minimo_b/maxserv)
        simClients <- simClients + 1
        sysClients[nodex] <- sysClients[nodex] - 1
        b <- b - tmin
        
        if (sysClients[nodex] < s[nodex]) {
          b[indice_minimo_b] <- NA
        } else {
          b[indice_minimo_b] <- tServ[iServ[nodex], nodex]
          iServ[nodex] <- iServ[nodex] + 1
        }
        
        n_s <- rand[simClients+staClients]
        acum <- 0
        node_dest <- 1
        while (acum < n_s) {
          acum <- acum + p[nodex, node_dest]
          node_dest <- node_dest+1
        }
        node_dest <- node_dest-1
        sysClients[node_dest] <- sysClients[node_dest] + 1
        in_node[node_dest] <- in_node[node_dest]+ 1
        
        if (sysClients[node_dest] <= s[node_dest]) {
          indice_minimo_nodo_minimo <- which(is.na(b[,node_dest]))
          b[indice_minimo_nodo_minimo[1], node_dest] <- tServ[iServ[node_dest], node_dest]
          iServ[node_dest] <- iServ[node_dest] + 1
        }
        
        cron[nodex] <- cron[nodex] + tmin
        for (j in 1:nodes) {
          if ((j == nodex) && (j != node_dest)) {
            c[j] <- c[j]+tmin*(sysClients[j]+1)
            if ((sysClients[j]+1) > s[j])
              d[j] <- d[j]+tmin*(sysClients[j]+1-s[j])
            prob[sysClients[j]+2, j] <- prob[sysClients[j]+2,j] + tmin
          } else {
            if ((j == node_dest) && (j != nodex)) {
              c[j] <- c[j] + tmin*(sysClients[j]-1)
              if ((sysClients[j]-1) > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
              prob[sysClients[j], j] <- prob[sysClients[j], j] + tmin
            } else {
              c[j] <- c[j] + tmin*sysClients[j]
              if (sysClients[j] > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-s[j])
              prob[sysClients[j]+1, j] <- prob[sysClients[j]+1, j] + tmin
            }
          }
        }
        if (historic && all(in_node>0)) {
          l <- c/sum(cron)
          lq <- d/sum(cron)
          w <- c/in_node
          wq <- d/in_node
          hist[,,simClients+staClients] <- array(c(l, lq, w, wq, cron), dim= c(nodes, 5))
        }
      }
      l <- c/sum(cron)
      lq <- d/sum(cron)
      w <- c/in_node
      wq <- d/in_node
      lqt <- sum(lq)
      rho <- l - lq
      eff <- w/(w-wq)
      nozero <- apply(prob, 1, function(v){any(v>(.Machine$double.eps ^ 0.5))})
      pn <- prob[1:max(which(nozero)),]/cron 
      obj$out <- list(l=l, lq=lq, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      #obj$out$data <- array(c(l, lq, w, wq), dim= c(nodes, 4), dimnames=list(1:nodes, c("L", "Lq", "W", "Wq")))
      if (historic) {
        obj$out$historic <- hist
      }
      oldClass(obj) <-  c("Closed", "SimulatedNetwork", "SimulatedModel")
      
      return(obj)
  }
  ParallelizeSimulations(ClosedNetwork_secuential, list(serviceDistribution, s, p, staClients, nClients, transitions, historic), nsim, nproc)
}

SON_Example <- function (sta, trans) {
  p <- matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2)
  s <- c(1, 2)
  arrivalDistribution <- c(Exp(20), Exp(30))
  serviceDistribution <- c(Exp(100), Exp(25))
  OpenNetwork(arrivalDistribution, serviceDistribution, s, p, sta, trans)
}

#' Obtains the main characteristics of an Open Network model by simulation
#'  
#' @param arrivalDistribution Vector indicating the arrival distribution at each node (Each element must be an object of S4-class \code{distr} 
#' defined in \pkg{distr} package or the no_distr() object)
#' @param serviceDistribution Vector indicating the service distribution at each node (Each element must be an object of S4-class \code{distr} 
#' defined in \pkg{distr} package)
#' @param s Vector of servers in each node
#' @param p Routing matrix, where \eqn{p_{ij}} is the routing probability from node i to node j
#' @param staClients Number of customers used in the stabilization stage
#' @param transitions Number of transitions between nodes used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of an Open network model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Vector of expected number of customers in the nodes: \eqn{L}}
#' \item{lq}{Vector of expected number of customers in the queues of the nodes: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{lqt}{Expected number of customers in all the queues: \ifelse{latex}{\eqn{L_{qTotal}}}{\out{<i>L<sub>qTotal</sub></i>}}}
#' \item{w}{Vector of expected waiting times in the nodes: \eqn{W}}
#' \item{wq}{Vector of expected waiting time in the queues of the nodes: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @examples
#' OpenNetwork(arrivalDistribution = c(Exp(20), no_distr()), 
#'             serviceDistribution = c(Exp(100), Exp(25)),
#'             s                   = c(1,2),
#'             p                   = matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2),
#'             staClients          = 10,
#'             transitions         = 100,
#'             nsim                = 10)
#' @export
#' @family SimulatedModels

         
OpenNetwork <- function(arrivalDistribution=c(Exp(20), Exp(30)), serviceDistribution=c(Exp(100), Exp(25)), s=c(1, 2), p=matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2), staClients=100, transitions=1000, historic=FALSE, nsim=10, nproc=1) {
  if (!is.numeric(nsim)) stop("Argument 'nsim' must be an integer greather than 0")
  if (!is.numeric(nproc)) stop("Argument 'nproc' must be an integer greather than 0")

  OpenNetwork_secuential <- function(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic) {
      nodesArrival <- length(arrivalDistribution)
      nodesService <- length(serviceDistribution)
      nodesS <- length(s)
      nodesP <- nrow(p)

      if (nodesArrival != nodesService || nodesArrival != nodesS || nodesArrival != nodesP || nodesService != nodesS || nodesService != nodesP || nodesS != nodesP)
        stop("The arguments 'arrivalDistribution', 'serviceDistribution', 's' and 'p' must have the same dimensions.")
      obj <- list(arrivalDistribution=arrivalDistribution, serviceDistribution=serviceDistribution, s= s, prob=p, staclients=staClients, transitions=transitions)
      
      nodeswitharrivals <- sapply(arrivalDistribution, class) != "no_distr"
      if (!any(nodeswitharrivals)) stop("At least one node must have an Arrival Distribution.")
      arrivalDistribution <- pairlist(arrivalDistribution[nodeswitharrivals], which(nodeswitharrivals))
       
      if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution' must be a valid Class of the Distr package")
      if (!belong(arrivalDistribution[[1]], "distr")) stop("All elements in argument 'ArrivalDistribution' must be a valid Class of the Distr package")
      if (!is.numeric(s) | any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
      if (!is.numeric(staClients) | staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (!is.numeric(transitions) | transitions < 0) stop("Argument 'transitions' must be equal or greather than 0.")
      if (!is.numeric(p))  stop("All elements in argument 'p' must be numerical.")
      
     
      nodes <- length(s)
      maxserv <- max(s)
      #Adaptamos la matriz de probabilidades para tener encuenta las salidas
      p <- matrix(c(1-rowSums(p), as.vector(p)), nrow=nrow(p), ncol=ncol(p)+1)
      #raux <- function(v) {if (is.na(v)) return(NA) else return(r(v))}
      map <- function(f) {f((staClients+transitions)*2)}
      tServ <- matrix(sapply(lapply(serviceDistribution, r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=nodes )
      tArr <- matrix(sapply(lapply(arrivalDistribution[[1]], r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=length(arrivalDistribution[[2]]), dimnames=list(1:((staClients+transitions)*2), arrivalDistribution[[2]]))
      if (any(tArr < 0)) stop("There's a problem with the Arrival Distribution, please check the parameters are correct.")
      if (any(tServ < 0))stop("There's a problem with the Service Distribution, please check the parameters are correct.")
      
      iServ <- rep(1, nodes)
      iArr <- rep(1, nodes)
      rand <- r(Unif()) ((staClients+transitions)*2)
      irand <- 1
      simClients <- 0
      
      if (historic) hist <- array(dim=c(nodes, 5, staClients+transitions), dimnames=list(1:nodes, c("L", "Lq", "W","Wq", "tClient"), 1:(staClients+transitions)))
      
      a <- array(NA, dim=nodes)
      for(i in arrivalDistribution[[2]]) {
        a[i] <- tArr[iArr[i], as.character(i)]
        iArr[i] <- iArr[i] + 1
      }
     b <- matrix(NA, nrow=maxserv, ncol=nodes)
     sysClients <- rep(0, nodes)
     cron <- 0
     c <- d <- entradas_nodo <- rep(0, nodes)
      while (simClients < staClients) {
        if (sum(sysClients) > 0) {
          indice_minimo_b <- which.min(b) #al ser los libre NA, devuelve el indice menor del que no es NA
          indice_minimo_a <- which.min(a)
          columna_minimo_b <- ceiling(indice_minimo_b/maxserv)
          tmin <- min(a[indice_minimo_a], b[indice_minimo_b])
        } else {
          indice_minimo_a <- which.min(a)
          tmin <- a[indice_minimo_a]
        }
        
        if (tmin == a[indice_minimo_a]) {
          simClients <- simClients + 1
          entradas_nodo[indice_minimo_a] <- entradas_nodo[indice_minimo_a] + 1
          if (sysClients[indice_minimo_a] < s[indice_minimo_a]) {
            sysClients[indice_minimo_a] <- sysClients[indice_minimo_a] + 1
            a <- a - tmin
            a[indice_minimo_a] <- tArr[iArr[indice_minimo_a], as.character(indice_minimo_a)]
            iArr[indice_minimo_a] <- iArr[indice_minimo_a] + 1
            
            b <- b- tmin
            indice_minimo_nodo_minimo <- which(is.na(b[, indice_minimo_a]))
            b[indice_minimo_nodo_minimo[1], indice_minimo_a] <- tServ[iServ[indice_minimo_a], indice_minimo_a]
            iServ[indice_minimo_a] <- iServ[indice_minimo_a] + 1  
          } else {
            sysClients[indice_minimo_a] <- sysClients[indice_minimo_a] + 1
            a <- a - tmin
            a[indice_minimo_a] <- tArr[iArr[indice_minimo_a], as.character(indice_minimo_a)]
            iArr[indice_minimo_a] <- iArr[indice_minimo_a] + 1
            b <- b - tmin
          }
          
          cron <- cron + tmin
          for (j in 1:nodes) {
            if (j == indice_minimo_a) {
              c[j] <- c[j]+ tmin*(sysClients[j]-1)
              if (sysClients[j]-1 > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
            } else {
              c[j] <- c[j] + tmin*sysClients[j]
              if (sysClients[j]>s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-s[j])
            }
          }
        } else {
          sysClients[columna_minimo_b] <- sysClients[columna_minimo_b] - 1
          b <- b - tmin
          
          if (sysClients[columna_minimo_b] < s[columna_minimo_b]) {
            b[indice_minimo_b] <- NA
          } else {
            b[indice_minimo_b] <- tServ[iServ[columna_minimo_b], columna_minimo_b]
            iServ[columna_minimo_b] <- iServ[columna_minimo_b] + 1
          }
          a <- a - tmin
          nsim <- rand[irand]
          irand <- irand+1
          acumulador <- 0
          j <- 0
          while (acumulador < nsim) {
            acumulador <- acumulador+ p[columna_minimo_b, j+1]
            j <- j + 1
          }
          j <- j-1
          nodo_transicion <- j
          
          if (j > 0) {
            sysClients[j] <- sysClients[j] + 1
            if (sysClients[j] <= s[j]) {
              indice_minimo_nodo_minimo <- which(is.na(b[,j]))
              b[indice_minimo_nodo_minimo[1], j] <- tServ[iServ[j], j]
              iServ[j] <- iServ[j] + 1
            }
          }
          
          for (j in 1:nodes) {
            if ((j==columna_minimo_b) && (j!=nodo_transicion)) {
              c[j] <- c[j] + tmin*(sysClients[j]+1)
              if (sysClients[j]+1 > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]+1-s[j])
            } else {
              if ((j == nodo_transicion) && (j!=columna_minimo_b)) {
                c[j] <- c[j] + tmin*(sysClients[j]-1)
                if (sysClients[j]-1 > s[j])
                  d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
              } else {
                c[j] <- c[j]+tmin*sysClients[j]
                if (sysClients[j] > s[j])
                  d[j] <- d[j] + tmin*(sysClients[j]-s[j])
              }         
            }
          }
        }
        if (historic && all(entradas_nodo>0)) {
          l <- c/cron
          lq <- d/cron
          w <- c/entradas_nodo
          wq <- d/entradas_nodo
          hist[,,simClients] <- array(c(l, lq, w, wq, cron), dim= c(nodes, 5))
        }
      }
      simClients <- 0
      cron <- 0
      c <- d <- entradas_nodo <- rep(0, nodes)
      prob <- matrix(c(0), nrow=(staClients+transitions), ncol=nodes)
      while (simClients < transitions) {
        if (sum(sysClients) > 0) {
          indice_minimo_b <- which.min(b) #al ser los libre NA, devuelve el indice menor del que no es NA
          indice_minimo_a <- which.min(a)
          if (length(indice_minimo_b) == 0)
            tmin <- a[indice_minimo_a]
          else {
            tmin <- min(a[indice_minimo_a], b[indice_minimo_b])
            columna_minimo_b <- ceiling(indice_minimo_b/maxserv)
          }
        } else {
          indice_minimo_a <- which.min(a)
          tmin <- a[indice_minimo_a]
        }
        if (tmin == a[indice_minimo_a]) {
          simClients <- simClients + 1
          entradas_nodo[indice_minimo_a] <- entradas_nodo[indice_minimo_a] + 1
          if (sysClients[indice_minimo_a] < s[indice_minimo_a]) {
            sysClients[indice_minimo_a] <- sysClients[indice_minimo_a] + 1
            a <- a - tmin
            a[indice_minimo_a] <- tArr[iArr[indice_minimo_a], as.character(indice_minimo_a)]
            iArr[indice_minimo_a] <- iArr[indice_minimo_a] + 1
            
            b <- b - tmin
            indice_minimo_nodo_minimo <- which(is.na(b[, indice_minimo_a]))
            b[indice_minimo_nodo_minimo[1], indice_minimo_a] <- tServ[iServ[indice_minimo_a], indice_minimo_a]
            iServ[indice_minimo_a] <- iServ[indice_minimo_a] + 1  
          } else {
            sysClients[indice_minimo_a] <- sysClients[indice_minimo_a] + 1
            a <- a - tmin
            a[indice_minimo_a] <- tArr[iArr[indice_minimo_a], as.character(indice_minimo_a)]
            iArr[indice_minimo_a] <- iArr[indice_minimo_a] + 1
            b <- b - tmin
          }
          
          cron <- cron + tmin
          for (j in 1:nodes) {
            if (j == indice_minimo_a) {
              c[j] <- c[j]+ tmin*(sysClients[j]-1)
              if (sysClients[j]-1 > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
              prob[sysClients[j], j] <- prob[sysClients[j], j] + tmin
            } else {
              c[j] <- c[j] + tmin*sysClients[j]
              if (sysClients[j]>s[j])
                d[j] <- d[j] + tmin*(sysClients[j]-s[j])
              prob[sysClients[j]+1, j] <- prob[sysClients[j]+1, j] + tmin
            }
          }
        } else {
          sysClients[columna_minimo_b] <- sysClients[columna_minimo_b] - 1
          b <- b - tmin
          
          if (sysClients[columna_minimo_b] < s[columna_minimo_b]) {
            b[indice_minimo_b] <- NA
          } else {
            b[indice_minimo_b] <- tServ[iServ[columna_minimo_b], columna_minimo_b]
            iServ[columna_minimo_b] <- iServ[columna_minimo_b] + 1
          }
          a <- a - tmin
          nsim <- rand[irand]
          irand <- irand+1
          acumulador <- 0
          j <- 0
          while (acumulador < nsim) {
            acumulador <- acumulador+ p[columna_minimo_b, j+1]
            j <- j + 1
          }
          j <- j-1
          nodo_transicion <- j
          
          if (j > 0) {
            sysClients[j] <- sysClients[j] + 1
            entradas_nodo[j] <- entradas_nodo[j]+1
            if (sysClients[j] <= s[j]) {
              indice_minimo_nodo_minimo <- which(is.na(b[,j]))
              b[indice_minimo_nodo_minimo[1], j] <- tServ[iServ[j], j]
              iServ[j] <- iServ[j] + 1
            }
          }
          
          cron <- cron + tmin
          for (j in 1:nodes) {
            if ((j==columna_minimo_b) && (j!=nodo_transicion)) {
              c[j] <- c[j] + tmin*(sysClients[j]+1)
              if (sysClients[j]+1 > s[j])
                d[j] <- d[j] + tmin*(sysClients[j]+1-s[j])
              prob[sysClients[j]+2, j] <- prob[sysClients[j]+2, j]+ tmin
            } else {
              if ((j == nodo_transicion) && (j!=columna_minimo_b)) {
                c[j] <- c[j] + tmin*(sysClients[j]-1)
                if (sysClients[j]-1 > s[j])
                  d[j] <- d[j] + tmin*(sysClients[j]-1-s[j])
                prob[sysClients[j], j] <- prob[sysClients[j], j] + tmin
              } else {
                c[j] <- c[j]+tmin*sysClients[j]
                if (sysClients[j] > s[j])
                  d[j] <- d[j] + tmin*(sysClients[j]-s[j])
                prob[sysClients[j]+1, j] <- prob[sysClients[j]+1, j] + tmin
              }         
            }
          }
        }
        if (historic && all(entradas_nodo>0)) {
          l <- c/cron
          lq <- d/cron
          w <- c/entradas_nodo
          wq <- d/entradas_nodo
          hist[,,simClients+staClients] <- array(c(l, lq, w, wq, cron), dim= c(nodes, 5))
        }
      }
      l <- c/cron
      lq <- d/cron
      w <- c/entradas_nodo
      wq <- d/entradas_nodo
      lqt <- sum(lq)
      lt <- sum(l)
      rho <- l - lq
      eff <- w/(w-wq)
      nozero <- apply(prob, 1, function(v){any(v>(.Machine$double.eps ^ 0.5))})
      pn <- prob[1:max(which(nozero)),]/cron
      obj$out <- list(l=l, lq=lq, lt=lt, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      if (historic) {
        obj$out$historic <- hist
      }
      oldClass(obj) <-  c("Open", "SimulatedNetwork", "SimulatedModel")
      return(obj)
  }
  ParallelizeSimulations(OpenNetwork_secuential, list(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic), nsim, nproc)
}

#' @describeIn Pn Implements the method for a Simulated model
#' @method Pn SimulatedModel
#' @usage NULL
#' @export
Pn.SimulatedModel <- function(qm, n) {
  if (any(oldClass(qm) == "SimulatedNetwork")) stop ("Use the method Pi instead.")
  
  ifelse(length(qm$out$pn) <= n, 0, qm$out$pn[n+1])
}

#' @rdname Pi
#' @method Pi SimulatedNetwork
#' @details
#' \code{Pi.SimulatedNetwork} implements the method for a SimulatedNetwork model
#' @export
Pi.SimulatedNetwork <- function(net, n, node) {
  ifelse(length(net$out$pn) <= n, 0, net$out$pn[n+1, node])
}


#' Shows basic statistical information of the variable
#' 
#' @param object SimulatedModel
#' @param var Main characteristic of an queue model
#' @param ... Further argumets
#' @method summary SimulatedModel
#' @details
#' \code{summary.SimulatedModel} implements the function for an object of class SimulatedModel.
#' @export
summary.SimulatedModel <- function(object, var="l", ...) {
  if (!any(var == c("l", "lq", "w", "wq", "rho", "eff"))) stop("Var must be one of the following values: 'l', 'lq', 'w', 'wq', 'rho' or 'eff'")
  getElement(object$out, var)$summary
}

#' Print the main characteristics of a SimulatedModel object
#' @param x SimulatedModel object
#' @param ... Further arguments passed to or from other methods.
#' @method print SimulatedModel
#' @keywords internal
#' @export
print.SimulatedModel <- function(x, ...) {
  cat("Model: ", class(x)[1], "\n")
  if (is.list(x$out$l)) {
    res <- data.frame("\u03BC" = c(x$out$l$mean, x$out$lq$mean, x$out$w$mean, x$out$wq$mean, x$out$rho$mean, x$out$eff$mean),
               "sd" = c(x$out$l$sd, x$out$lq$sd, x$out$w$sd, x$out$wq$sd, x$out$rho$sd, x$out$eff$sd),
                row.names = c("L", "Lq", "W", "Wq", "Intensity", "Efficiency"))
    print(res)
  }
  else {
    cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensity =\t", x$out$rho , "\n")
    cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEfficiency =\t", x$out$eff, "\n\n")
  }
}

#' Print the main characteristics of a SimulatedNetwork object
#' @param x SimulatedNetwork object
#' @param ... Further arguments passed to or from other methods.
#' @method print SimulatedNetwork
#' @keywords internal
#' @export
print.SimulatedNetwork <- function(x, ...) {
  cat("Model: ", class(x)[1], "\n")
  if (is.list(x$out$l)) {
    print(data.frame("L"=x$out$l[c("mean", "sd")], "Lq"=x$out$lq[c("mean", "sd")], "W"=x$out$w[c("mean", "sd")], "Wq"=x$out$wq[c("mean", "sd")]))
  }else{
    print(data.frame("L" = x$out$l, "Lq" = x$out$lq, "W" = x$out$w, "Wq"= x$out$wq))
    
  }
}


# #' Shows a plot of the evolution of a variable during the simulation
# #' 
# #' @param x Simulated Model
# #' @param minrange Number of customers needed to establish the start of the plot
# #' @param maxrange Number of customers needed to establish the end of the plot
# #' @param var This variable indicates the parameter of the queue to show in graphic (L, Lq, W, Wq, Clients, Intensity)
# #' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
# #' @param depth Number of points printed in the plot
# #' @param nSimulation Only used when the var param is equal to "Clients". Selects one of the multiple simulations to show the evolution of the Clients.

# #' @param ... Further arguments passed to or from other methods.
# #' @export
# summarySimple <- function(object, minrange, maxrange, var, graphics, depth, ...) {UseMethod("summarySimple", object)}


#' @rdname plot
#' @method plot SimulatedModel
#' @details
#' \code{plot.SimulatedModel} implements the function for an object of class SimulatedModel.
#' @export
plot.SimulatedModel <- function(x, minrange=1, maxrange, var="L", graphics="ggplot2", depth=maxrange-minrange, ...) {
  if (is.null(x$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  if (length(intersect(var, c("L", "Lq", "W", "Wq", "Clients", "Intensity"))) <= 0) stop("Argument 'var' must be any of the following values: 'L', 'Lq', 'W', 'Wq', 'Clients', 'Intensity'.")
  truerange <- seq(minrange, maxrange, length.out=depth)
  #eval(parse(text=paste("data <- data.frame(n=x$out$historic[truerange, 'tClient'], '", var, "'=x$out$historic[truerange, '", var, "'])", sep="")))
  eval(parse(text=paste("data <- data.frame(n=truerange, '", var, "'=x$out$historic[truerange, '", var, "'])", sep="")))
  switch(graphics,
        "graphics" =  {eval(parse(text=paste("plot(data$n, data$", var, ", col='red', type='l')\n
                                              abline(v=x$Staclients, untf = FALSE, col='black')\n
                                              title(main='Evolution of ", var, "')", sep="")))},
         "ggplot2" = {data <- reshape::melt(data, id.var="n")
                      eval(parse(text=paste("ggplot2::qplot(n, value, data=data, geom='line', colour=variable, 
                      main='Evolution of ", var, ".', ylab='", var , "') + geom_vline(xintercept=x$out$historic[x$Staclients, 'tClient'], linetype='dotdash') + scale_colour_discrete(name='') + theme(legend.position='none')", sep="")))})
}

#' @rdname plot
#' @method plot SimulatedNetwork
#' @details
#' \code{plot.SimulatedNetwork} implements the function for an object of class SimulatedNetwork.
#' @export
plot.SimulatedNetwork <- function(x, minrange=1, maxrange, var="L", graphics="ggplot2", depth=maxrange-minrange, nSimulation=NULL, ...) {
  if (is.null(x$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  if (length(intersect(var, c("L", "Lq", "W", "Wq", "Clients", "Intensity"))) <= 0) stop("Argument 'var' must be any of the following values: 'L', 'Lq', 'W', 'Wq', 'Clients', 'Intensity'.")
  truerange <- seq(minrange, maxrange, length.out=depth)
  nnodes <- dim(x$out$historic)[1]
  eval(parse(text=paste("data <- data.frame(n=NULL, var = NULL, node=NULL)\n
                         for(i in 1:nnodes)\n 
                            data <- rbind(data, data.frame(n=1:length(truerange), var=x$out$historic[i, '", var, "', truerange], node=rep(i, length(truerange))))\n",sep="")))
  switch(graphics,
         "graphics" =  {plot(data$n[data$node==1], data$var[data$node==1], col=2, type='l', xlab="n", ylab=var)
                        if (nnodes > 1){
                          for(j in 2:nnodes){
                            lines(data$n[data$node==j], data$var[data$node==j], col=j+1, type='l')
                          }
                        }
                        abline(v=x$Staclients, untf = FALSE, col='black')
                        legend('topright', paste('Node', 1:nnodes, sep=' '), lty =c(1), col = (1:nnodes)+1, bty='n')
                        title(main=paste("Evolution of ", var, sep=""))},
         "ggplot2" = {
                      ggplot2::ggplot(data, aes(x=n, y=var, colour=factor(node), order=factor(node))) +
                               geom_line(na.rm=TRUE) +  xlab("n") + ylab(var) + 
                               ggtitle(paste("Evolution of ", var, sep="")) + scale_colour_discrete(name="Nodes")
                     }
  )
}

#' @rdname plot
#' @method plot list
#' @param showMean Shows the mean of all the simulations
#' @param showValues Shows the values of all the simulations
#' @details
#' \code{plot.list} implements the function for an object of class list
#' @export
plot.list <- function(x, minrange=1, maxrange, var="L", graphics="ggplot2", depth=maxrange-minrange+1, nSimulation=1, showMean=TRUE, showValues=TRUE, ...) {
      if (length(intersect(var, c("L", "Lq", "W", "Wq", "Clients", "Intensity"))) <= 0) stop("Argument 'var' must be any of the following values: 'L', 'Lq', 'W', 'Wq', 'Clients', 'Intensity'.")
      if (!is.numeric(nSimulation)) nSimulation <- 1
      firstHistoric <- x[[1]]$out$historic
      if (length(dim(firstHistoric)) > 2)
        plotdata <- data.frame(x=NULL, val=NULL, sim=NULL)
      else
        plotdata <- data.frame(x=NULL, val=NULL, sim=NULL, node=NULL)
      index <- 1
      dimensions <- NULL
      lapply(x, function(qm){
          historic <- getElement(getElement(qm, "out"), "historic")
          dimensions <<- dim(historic)
          truerange <- seq(minrange, maxrange, length.out=depth)
          if (is.null(historic)) stop(simpleError("The argument 'historic' must be TRUE to show the plot."))
          if (length(dim(historic)) == 2)
            plotdata <<- rbind(plotdata, data.frame(x=truerange, val=historic[truerange, var], sim=rep(index, length(truerange))))
          else {
            if (var == "Clients" || var=="Intensity") stop("The argument 'var' can't be 'Clients' for a Network.")
            for (i in 1:dimensions[1])
              plotdata <<- rbind(plotdata, data.frame(x=truerange, val=historic[i, var, truerange], sim=rep(index, length(truerange)), node=rep(i, length(truerange))))
          }
          index <<- index+1
      })
      if (showMean) {
        if (length(dimensions) == 2) {
          meanvalue <- sapply(x, function(foo){getElement(getElement(foo, "out"), "historic")[,var]})[seq(minrange, maxrange, length.out=depth), ]
          meanvalue <- apply(meanvalue, 1, mean, na.rm=TRUE)
          meanvalue <- data.frame(x = seq(minrange, maxrange, length.out=depth), val=meanvalue)
        } else {
          meanvalue <- data.frame(x=NULL, val=NULL, node=NULL)
          for (i in 1:dimensions[1]) {
            meanaux <- sapply(x, function(foo){getElement(getElement(foo, "out"), "historic")[i,var,]})[seq(minrange, maxrange, length.out=depth), ]
            meanaux <- apply(meanaux, 1, mean, na.rm=TRUE)
            meanvalue <- rbind(meanvalue, data.frame(x=seq(minrange, maxrange, length.out=depth), val=meanaux, node=rep(i, depth)))
          }
        }
      }
      sim <- x <- val <- NULL
      switch(graphics, 
      "graphics"= {if (length(dimensions) == 2) 
                     if (var == "Clients")
                        barplot(plotdata$val[plotdata$sim==nSimulation], main = paste("Evolution of Clients in simulation ", nSimulation, sep=""), xlab="n", ylab="Number of clients")
                     else {
                        plot(plotdata$x[plotdata$sim==1], plotdata$val[plotdata$sim==1], col=1, type='l', xlab="n", ylab=var)
                        if (dimensions[1] > 1)
                            for(j in 2:dimensions[1])
                              lines(plotdata$x[plotdata$sim==j], plotdata$val[plotdata$sim==j], col=1, type='l')
                        lines(meanvalue$x, meanvalue$val, col=2, type='l')
                        legend('topright', legend = c("Mean value"), col = 2, cex = 0.8, pch = 1)
                     }
                   else {
                     
                   } 
      },
      "ggplot2" = {if (length(dimensions) == 2) 
                      if (var == "Clients")
                        ggplot2::ggplot(subset(plotdata, sim == nSimulation), aes(x=x, y=val, order=factor(sim))) + geom_histogram(stat="identity", na.rm=TRUE) + xlab("n") + ylab("Number of clients") + ggtitle(paste("Evolution of ", var, " in simulation ", nSimulation, sep="")) + theme(legend.position="none")
                      else {
                        plotbase <- ggplot2::ggplot()
                        if (showValues)
                          plotbase <- plotbase + ggplot2::geom_line(data=plotdata, aes(x=x, y=val, order=factor(sim)), na.rm=TRUE)
                        if (showMean)
                          plotbase <- plotbase + ggplot2::geom_line(data=meanvalue, aes(x=x, y=val, colour='red'))
                        
                        plotbase + ggplot2::xlab("n") + ggplot2::ylab(var) +
                                   ggplot2::ggtitle(paste("Evolution of ", var, sep="")) +
                                   ggplot2::scale_colour_discrete(name  =element_blank(), breaks=c("red"), labels=c("Mean value")) +
                                   ggplot2::theme(legend.position="top")
                      }
                  else {
                    plotbase <- ggplot2::ggplot()
                    if (showValues) 
                      plotbase <- plotbase + geom_line(data=plotdata, aes(x=x, y=val, colour=factor(node), alpha=0.95, order=factor(sim)), na.rm=TRUE) + scale_colour_discrete(name="Nodes") + scale_alpha_continuous(name="", breaks=NULL, labels=NULL)
                    if (showMean)
                      plotbase <- plotbase + geom_line(data=meanvalue, aes(x=x, y=val, colour=factor(node)), size=1) + scale_colour_discrete(name="Nodes")
                    
                    plotbase + xlab("n") + ylab(var) +
                               ggtitle(paste("Evolution of ", var, sep="")) 
                  }
                }
      )
}
