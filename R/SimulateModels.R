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

#' Combines a list of independent simulations of a queueing model, calculating the mean and variance of each characteristic of interest
#' 
#' @param listsims A list of independent simulations
#' @return an object with the mean and variance of parameters L, Lq, W, Wq, Rho and Eff.
#' @export
combineSimulations <- function(listsims) {
  if (class(listsims)[1] != "list") 
    if (length(intersect(class(listsims), "SimulatedModel")) > 0)
      return(listsims)
    else
      stop(simpleError("The argument listsims must be a list of queue models"))
  else if (length(listsims) == 1) return(listsims[[1]])
  
  getAttribute <- function(qm, attr) {return(getElement(getElement(qm, "out"), attr))}

  l <- sapply(listsims, getAttribute, "l")
  lq <- sapply(listsims, getAttribute, "lq")
  w <- sapply(listsims, getAttribute, "w")
  wq <- sapply(listsims, getAttribute, "wq")
  rho <- sapply(listsims, getAttribute, "rho")
  eff <- sapply(listsims, getAttribute, "eff")
  
  res <- listsims[[1]]
  k <- 2/sqrt(length(args))
  res$out$l <- if(is.null(nrow(l)))
                 list(mean=mean(l), error=k*sd(l))
               else
                 list(mean=sapply(1:nrow(l), function(i){mean(l[i,])}), error=sapply(1:nrow(l), function(i){k*sd(l[i,])}))
  res$out$lq <- if(is.null(nrow(lq)))
                   list(mean=mean(lq), error=k*sd(lq))
                else
                 list(mean=sapply(1:nrow(lq), function(i){mean(lq[i,])}), error=sapply(1:nrow(lq), function(i){k*sd(lq[i,])}))
  res$out$w <- if(is.null(nrow(w)))
                 list(mean=mean(w), error=k*sd(w))
               else
                 list(mean=sapply(1:nrow(w), function(i){mean(w[i,])}), error=sapply(1:nrow(w), function(i){k*sd(w[i,])}))
  res$out$wq <- if(is.null(nrow(wq)))
                  list(mean=mean(wq), error=k*sd(wq))
                else
                  list(mean=sapply(1:nrow(wq), function(i){mean(wq[i,])}), error=sapply(1:nrow(wq), function(i){k*sd(wq[i,])}))
  res$out$rho <- list(mean=mean(rho), error=k*sd(rho))
  res$out$eff <- list(mean=mean(eff), error=k*sd(eff))
  return(res)
}

#' Obtains the main characteristics of a G/G/1 model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1 model:
#' \item{pn}{Stores all the positives steady-state probabilities of having n customers, with n from 0 to staClients+nClients: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{<i>L, L<sub>q</sub>, W, W<sub>q</sub></i>}},\emph{Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels

G_G_1 <- function(arrivalDistribution = Exp(3), serviceDistribution = Exp(6), staClients = 100, nClients = 1000, historic = FALSE, nsim=10, nproc=1) {
      G_G_1_secuential <- function(arrivalDistribution, serviceDistribution, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")

        tArr <- r(arrivalDistribution) (staClients+nClients)
        tServ <- r(serviceDistribution) (staClients+nClients)
        if (any(tArr < 0)) stop("There's a problem with de Arrival Distribution, please check the parameters are correct.")
        if (any(tServ < 0))stop("There's a problem with de Service Distribution, please check the parameters are correct.")
        
        iArr <- iServ <- 1
        sysClients <- simClients<- 0
        a <- 0
        b <- -1
        cron <- d <- c <- 0
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Staclients=staClients, Nclients=nClients)
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
          #En cada iteración almacenamos la evolucion de los valores
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
          pn <- c()
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_1", "SimulatedModel")
        
        return(obj)
    }      
    if (nproc > 1) {
      cl <- makeCluster(nproc)
      registerDoParallel(cl, cores=nproc)
      parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
        belong <- function(obj, packagename) {
          if (length(obj) > 1)
            packageobj <- sapply(lapply(obj, class), attr, "package")
          else
            packageobj <- attr(class(obj), "package")
          
          return(all(packageobj==packagename))
        }
        res <- list()
        for(i in 1:simulations) {
          res <- c(res, list(G_G_1_secuential(arrivalDistribution, serviceDistribution, staClients, nClients, historic)))
        }
        return(res)
      }
      stopCluster(cl)
      return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
    } else {
      res <-list()
      for(i in 1:nsim)
        res <- c(res, list(G_G_1_secuential(arrivalDistribution, serviceDistribution, staClients, nClients, historic)))
      return(if(length(res)==1) res[[1]] else res)
    }
}

#' Obtains the main characteristics of a G/G/s model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S model:
#' \item{pn}{vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S <- function (arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=2, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
      G_G_S_secuential <- function(arrivalDistribution, serviceDistribution, s, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (s <= 0) stop("Argument 's' must be greather than 0.")
        if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) (nClients+staClients)
        tServ <- r(serviceDistribution) (nClients+staClients)
        iArr <-1
        iServ <- 0
        bussyservs <- numeric()
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
       
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, Staclients=staClients, Nclients=nClients)
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
          #En cada iteración almacenamos la evolucion de los valores
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
        pn <- tnClients[1:max(which(tnClients>0))]/cron
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_S", "SimulatedModel")
        
        return(obj)
    }
    if (nproc > 1) {
      cl <- makeCluster(nproc)
      registerDoParallel(cl, cores=nproc)
      parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
        belong <- function(obj, packagename) {
          if (length(obj) > 1)
            packageobj <- sapply(lapply(obj, class), attr, "package")
          else
            packageobj <- attr(class(obj), "package")
          
          return(all(packageobj==packagename))
        }
        res <- list()
        for(i in 1:simulations) {
          res <- c(res, list(G_G_S_secuential(arrivalDistribution, serviceDistribution, s, staClients, nClients, historic)))
        }
        return(res)
      }
      stopCluster(cl)
      return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
    } else {
      res <-list()
      for(i in 1:nsim)
        res <- c(res, list(G_G_S_secuential(arrivalDistribution, serviceDistribution, s, staClients, nClients, historic)))
      return(if (length(res)==1) res[[1]] else res)
    }
}

#'Obtains the main characteristics of a G/G/1/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1/K model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @export
#' @family SimulatedModels

G_G_1_K <- function (arrivalDistribution=Exp(3), serviceDistribution=Exp(6), K=2, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
      G_G_1_K_secuential <- function(arrivalDistribution, serviceDistribution, K, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (K <= 0) stop("Argument 'K' must be greather than 0.")
        if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        iArr <- iServ <- 1
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
        tnClients <- numeric(nClients)
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Staclients=staClients, Nclients=nClients)
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
          #En cada iteración almacenamos la evolucion de los valores
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
        pn <- tnClients[1:max(which(tnClients>0))]/cron
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_1_K", "SimulatedModel")
        
        return(obj)
    }
    if (nproc > 1) {
      cl <- makeCluster(nproc)
      registerDoParallel(cl, cores=nproc)
      parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
        belong <- function(obj, packagename) {
          if (length(obj) > 1)
            packageobj <- sapply(lapply(obj, class), attr, "package")
          else
            packageobj <- attr(class(obj), "package")
          
          return(all(packageobj==packagename))
        }
        res <- list()
        for(i in 1:simulations) {
          res <- c(res, list(G_G_1_K_secuential(arrivalDistribution, serviceDistribution, K, staClients, nClients, historic)))
        }
        return(res)
      }
      stopCluster(cl)
      return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
    } else {
      res <-list()
      for(i in 1:nsim)
        res <- c(res, list(G_G_1_K_secuential(arrivalDistribution, serviceDistribution, K, staClients, nClients, historic)))
      return(if (length(res)==1) res[[1]] else res)
    }
}

#' Obtains the main characteristics of a G/G/s/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S/K model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S_K <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=2, K=3, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
      G_G_S_K_secuential <- function(arrivalDistribution, serviceDistribution, s, K, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
        if (s <= 0) stop("Argument 's' must be greather than 0.")
        if (K <= 0) stop("Argument 'K' must be greather than 0.")
        if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        iArr <- iServ <- 1
        bussyservers <- rep(NA, s)
        sysClients <- simClients<- 0
        cron <- d <- c <- 0
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, Staclients=staClients,  Nclients=nClients)
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
          #En cada iteración almacenamos la evolucion de los valores
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
        pn <- tnClients[1:max(which(tnClients>0))]/cron
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        oldClass(obj) <-  c("G_G_S_K", "SimulatedModel")
        
        return(obj)
    }
    if (nproc > 1) {
      cl <- makeCluster(nproc)
      registerDoParallel(cl, cores=nproc)
      parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
        belong <- function(obj, packagename) {
          if (length(obj) > 1)
            packageobj <- sapply(lapply(obj, class), attr, "package")
          else
            packageobj <- attr(class(obj), "package")
          
          return(all(packageobj==packagename))
        }
        res <- list()
        for(i in 1:simulations) {
          res <- c(res, list(G_G_S_K_secuential(arrivalDistribution, serviceDistribution, s, K, staClients, nClients, historic)))
        }
        return(res)
      }
      stopCluster(cl)
      return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
    } else {
      res <-list()
      for(i in 1:nsim)
        res <- c(res, list(G_G_S_K_secuential(arrivalDistribution, serviceDistribution, s, K, staClients, nClients, historic)))
      return(if (length(res)==1) res[[1]] else res)
    }
}


#' Obtains the main characteristics of a G/G/1/\eqn{\infty}/H model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/1/\eqn{\infty}/H model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels

G_G_1_INF_H <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), H=5, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  G_G_1_INF_H_secuential <- function(arrivalDistribution, serviceDistribution, H, staClients, nClients, historic) {
    if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
    if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
    if (H <= 0) stop("Argument 'H' must be greather than 0.")
    if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
    if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
    
    tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
    tServ <- r(serviceDistribution) ((nClients+staClients)*2)
    possibleClients <- tArr[1:H]
    iServ <- 1
    iArr <- H+1
    sysClients <- 0
    simClients <- cron <- c <- d <- 0
    
    obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, H=H, Staclients=staClients, Nclients=nClients)
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
      #En cada iteración almacenamos la evolucion de los valores
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
    #pn <- tnClients[1:max(which(tnClients>0))]/cron
    pn <- NA
    if (historic)
      obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
    else
      obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
    
    oldClass(obj) <-  c("G_G_1_INF_H", "SimulatedModel")
    return(obj)
  }
  if (nproc > 1) {
    cl <- makeCluster(nproc)
    registerDoParallel(cl, cores=nproc)
    parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(G_G_1_INF_H_secuential(arrivalDistribution, serviceDistribution, H, staClients, nClients, historic)))
      }
      return(res)
    }
    stopCluster(cl)
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  } else {
    res <-list()
    for(i in 1:nsim)
      res <- c(res, list(G_G_1_INF_H_secuential(arrivalDistribution, serviceDistribution, H, staClients, nClients, historic)))
    return(if (length(res)==1) res[[1]] else res)
  }
}

#' Obtains the main characteristics of a G/G/S/\eqn{\infty}/H  model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/S/\eqn{\infty}/H model
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system:\eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S_INF_H <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=3, H=5, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
      G_G_S_INF_H_secuential <- function(arrivalDistribution, serviceDistribution, s, H, staClients, nClients, historic) {
          if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
          if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
          if (s <= 0) stop("Argument 's' must be greather than 0.")
          if (H <= 0) stop("Argument 'H' must be greather than 0.")
          if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
          if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
          tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
          tServ <- r(serviceDistribution) ((nClients+staClients)*2)
          possibleClients <- tArr[1:H]
          bussyservers <- rep(NA, s)
          iServ <- 1
          iArr <- H+1
          sysClients <- 0
          simClients <- cron <- c <- d <- 0
          nClicien <- 100/nClients
          obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, H=H, Staclients=staClients, Nclients=nClients)
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
            #En cada iteración almacenamos la evolucion de los valores
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
          pn <- tnClients[1:max(which(tnClients>0))]/cron
          if (historic)
            obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
          else
            obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
          oldClass(obj) <-  c("G_G_S_INF_H", "SimulatedModel")
          return(obj)
    }
if (nproc > 1) {
  cl <- makeCluster(nproc)
  registerDoParallel(cl, cores=nproc)
  parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
    belong <- function(obj, packagename) {
      if (length(obj) > 1)
        packageobj <- sapply(lapply(obj, class), attr, "package")
      else
        packageobj <- attr(class(obj), "package")
      
      return(all(packageobj==packagename))
    }
    res <- list()
    for(i in 1:simulations) {
      res <- c(res, list(G_G_S_INF_H_secuential(arrivalDistribution, serviceDistribution, s, H, staClients, nClients, historic)))
    }
    return(res)
  }
  stopCluster(cl)
  return(ifelse(length(parallelRes) == 1, parallelRes[[1]], parallelRes))
} else {
  res <-list()
  for(i in 1:nsim)
    res <- c(res, list(G_G_S_INF_H_secuential(arrivalDistribution, serviceDistribution, s, H, staClients, nClients, historic)))
  return(if (length(res)==1) res[[1]] else res)
}
}

#' Obtains the main characteristics of a G/G/S/\eqn{\infty}/H with Y replacements model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
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
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels
G_G_S_INF_H_Y <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), s=3, H=5, Y=3, staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  G_G_S_INF_H_Y_secuential <- function(arrivalDistribution, serviceDistribution, s, H, Y, staClients, nClients, historic) {
        if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
        if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution' must be a valid Class of the Distr package")
        if (s <= 0) stop("Argument 's' must be greather than 0.")
        if (H <= 0) stop("Argument 'H' must be greather than 0.")
        if (Y <= 0) stop("Argument 'Y' must be greather than 0.")
        if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
        if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
        
        tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
        tServ <- r(serviceDistribution) ((nClients+staClients)*2)
        possibleClients <- tArr[1:H]
        #   bussyservers <- rep(-1, s)
        bussyservers <- rep(NA, s)
        iServ <- 1
        iArr <- H+1
        sysClients <- 0
        simClients <- cron <- c <- d <- 0
        
        obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, H=H, Y=Y, Staclients=staClients, Nclients=nClients)
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
          #En cada iteración almacenamos la evolucion de los valores
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
        pn <- tnClients[1:max(which(tnClients>0))]/cron
        if (historic)
          obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        else
          obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
        
        oldClass(obj) <-  c("G_G_S_INF_H_Y", "SimulatedModel")
        return(obj)
  }
  if (nproc > 1) {
    cl <- makeCluster(nproc)
    registerDoParallel(cl, cores=nproc)
    parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(G_G_S_INF_H_Y_secuential(arrivalDistribution, serviceDistribution, s, H, Y, staClients, nClients, historic)))
      }
      return(res)
    }
    stopCluster(cl)
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  } else {
    res <-list()
    for(i in 1:nsim)
      res <- c(res, list(G_G_S_INF_H_Y_secuential(arrivalDistribution, serviceDistribution, s, H, Y, staClients, nClients, historic)))
    return(if (length(res)==1) res[[1]] else res)
  }
}

#' Obtains the main characteristics of a G/G/\eqn{\infty} model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param staClients Number of customers used in stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @param nsim Number of simulations
#' @param nproc Processors used in the simulation.
#' @return
#' Returns the next information of a G/G/\eqn{\infty} model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_{n}}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Expected number of customers in the system: \eqn{L}}
#' \item{lq}{Expected number of customers in the queue: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{w}{Expected waiting time in the system: \eqn{W}}
#' \item{wq}{Expected waiting time in the queue: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation}
#' @export
#' @family SimulatedModels
G_G_INF <- function(arrivalDistribution=Exp(3), serviceDistribution=Exp(6), staClients=100, nClients=1000, historic=FALSE, nsim=10, nproc=1) {
  G_G_INF_secuential <- function(arrivalDistribution, serviceDistribution, staClients, nClients, historic) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
      tServ <- r(serviceDistribution) ((nClients+staClients)*2)
      iServ <- iArr <- 1
      sysClients <- 0
      simClients <- cron <- c <- d <- 0
      a <- 0
      b <- -1
      
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Staclients=staClients, Nclients=nClients)
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
        #En cada iteración almacenamos la evolucion de los valores
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
      pn <- tnClients[1:max(which(tnClients>0))]/cron
      if (historic)
        obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      else
        obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      oldClass(obj) <-  c("G_G_INF", "SimulatedModel")
      
      return(obj)
  }
  if (nproc > 1) {
    cl <- makeCluster(nproc)
    registerDoParallel(cl, cores=nproc)
    parallelRes <- foreach(simulations=idiv(nsim, chunks=nproc), .combine='c', .packages=c("distr", "doParallel")) %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(G_G_INF_secuential(arrivalDistribution, serviceDistribution, staClients, nClients, historic)))
      }
      return(res)
    }
    stopCluster(cl)
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  } else {
    res <-list()
    for(i in 1:nsim)
      res <- c(res, list(G_G_INF_secuential(arrivalDistribution, serviceDistribution, staClients, nClients, historic)))
    return(if (length(res)==1) res[[1]] else res)
  }
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
#' @param serviceDistribution Service distributions for the nodes of the network
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
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system: \ifelse{latex}{\eqn{P_n}}{\out{<i>P<sub>n</sub></i>}}}
#' \item{l}{Vector of expected number of customers in the nodes: \eqn{L}}
#' \item{lq}{Vector of expected number of customers in the queues of the nodes: \ifelse{latex}{\eqn{L_{q}}}{\out{<i>L<sub>q</sub></i>}}}
#' \item{lqt}{Expected number of customers in the all queues: \ifelse{latex}{\eqn{L_{qTotal}}}{\out{<i>L_<sub>qTotal</sub></i>}}}
#' \item{w}{Vector of expected waiting times in the nodes: \eqn{W}}
#' \item{wq}{Vector of expected waiting times in the queues of the nodes: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @export
#' @family SimulatedModels

ClosedNetwork <- function(serviceDistribution=c(Exp(5), Exp(5), Exp(10), Exp(15)), s=c(2,2,1,1), p=array(c(0.25,0.15,0.5,0.4,0.15,0.35,0.25,0.3,0.2,0.2,0.15,0.25,0.4,0.30,0.1,0.05), dim=c(4,4)), staClients=100, nClients=3, transitions=1000, historic=FALSE, nsim=10, nproc=1) {
  ClosedNetwork_secuential <- function(serviceDistribution, s, p, staClients, nClients, transitions, historic) {
      if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution' must be a valid Class of the Distr package")
      if (any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
      if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      if (transitions < 0) stop("Argument 'transitions' must be equal or greather than 0.")
      
      nodesServ <- length(serviceDistribution)
      nodesS <- length(s)
      nodesP <- nrow(p)
    
      #print(paste(nodesServ, ", ", nodesS, ", ", nodesP))
      if (nodesServ != nodesS || nodesServ != nodesP || nodesS != nodesP) stop("The arguments 'serviceDistribution', 's' and 'p' must have the same dimensions.")
      
      nodes <- length(s)
      maxserv <- max(s)
      map <- function(f) {f(staClients+transitions)}
      tServ <- matrix(sapply(lapply(serviceDistribution, r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=nodes )
      iServ <- rep(1, nodes)
      rand <- r(Unif()) (staClients+transitions)
      simClients <- 0
      c <- d <- cron <- in_node <- rep(0, nodes)
      
      obj <- list(serviceDistribution=serviceDistribution, Servs= s, Prob=p, Staclients=staClients, Transitions=transitions, Nclients=nClients)
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
      nozero <- apply(prob, 1, function(v){any(v>0)})
      pn <- prob[1:max(which(nozero)),]/cron 
      obj$out <- list(l=l, lq=lq, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      #obj$out$data <- array(c(l, lq, w, wq), dim= c(nodes, 4), dimnames=list(1:nodes, c("L", "Lq", "W", "Wq")))
      if (historic) {
        obj$out$historic <- hist
      }
      oldClass(obj) <-  c("Closed", "SimulatedNetwork", "SimulatedModel")
      
      return(obj)
  }
  if (nproc > 1) {
    cluster <- makeCluster(nproc)
    registerDoParallel(cluster, cores=nproc)
    parallelRes <- foreach(simulations=idiv(nsim, chunks = nproc), .combine=c, .multicombine=TRUE, .packages="distr") %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(ClosedNetwork_secuential(serviceDistribution, s, p, staClients, nClients, transitions, historic)))
      }
      return(res)
    }
    stopCluster(cluster)
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  }
  res <-list()
  for(i in 1:nsim) {
    res <- c(res, list(ClosedNetwork_secuential(serviceDistribution, s, p, staClients, nClients, transitions, historic)))
  }
  return(if (length(res)==1) res[[1]] else res)
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
#' @param arrivalDistribution PairList indicating the arrival distribution and the node that uses it.
#' @param serviceDistribution Vector of service distribution in each node
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
#' \item{lqt}{Expected number of customers in all queues: \ifelse{latex}{\eqn{L_{qTotal}}}{\out{<i>L<sub>qTotal</sub></i>}}}
#' \item{w}{Vector of expected waiting times in the nodes: \eqn{W}}
#' \item{wq}{Vector of expected waiting time in the queues of the nodes: \ifelse{latex}{\eqn{W_{q}}}{\out{<i>W<sub>q</sub></i>}}}
#' \item{eff}{System efficiency: \ifelse{latex}{\eqn{Eff = W/(W-W_{q})}}{\out{<i>Eff = W/(W-W<sub>q</sub></i>)}}}
#' \item{rho}{Traffic intensity: \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \ifelse{latex}{\eqn{L}, \eqn{L_q}, \eqn{W} and  \eqn{W_q}}{\out{L, L<sub>q</sub>, W, W<sub>q</sub>}}\emph{, Customers in the system, Rho and Elapsed time} during the simulation.}
#' @export
#' @family SimulatedModels

         
OpenNetwork <- function(arrivalDistribution=c(Exp(20), Exp(30)), serviceDistribution=c(Exp(100), Exp(25)), s=c(1, 2), p=matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2), staClients=100, transitions=1000, historic=FALSE, nsim=10, nproc=1) {
  OpenNetwork_secuential <- function(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic) {
      nodesArrival <- length(arrivalDistribution)
      nodesService <- length(serviceDistribution)
      nodesS <- length(s)
      nodesP <- nrow(p)
      
      obj <- list(arrivalDistribution=arrivalDistribution, serviceDistribution=serviceDistribution, Servs= s, Prob=p, Staclients=staClients, Transitions=transitions)
      if (nodesArrival != nodesService || nodesArrival != nodesS || nodesArrival != nodesP || nodesService != nodesS || nodesService != nodesP || nodesS != nodesP)
        stop("The arguments 'arrivalDistribution', 'serviceDistribution', 's' and 'p' must have the same dimensions.")
    
      nodeswitharrivals <- sapply(arrivalDistribution, class) != "no_distr"
      if (!any(nodeswitharrivals)) stop("At least one node must have an Arrival Distribution.")
      arrivalDistribution <- pairlist(arrivalDistribution[nodeswitharrivals], which(nodeswitharrivals))
       
      if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution' must be a valid Class of the Distr package")
      if (!belong(arrivalDistribution[[1]], "distr")) stop("All elements in argument 'ArrivalDistribution' must be a valid Class of the Distr package")
      if (any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
      if (staClients < 0) stop("Argument 'staClients' must be equal or greather than 0.")
      if (transitions < 0) stop("Argument 'transitions' must be equal or greather than 0.")
      
      nodes <- length(s)
      maxserv <- max(s)
      #Adaptamos la matriz de probabilidades para tener encuenta las salidas
      p <- matrix(c(1-rowSums(p), as.vector(p)), nrow=nrow(p), ncol=ncol(p)+1)
      #raux <- function(v) {if (is.na(v)) return(NA) else return(r(v))}
      map <- function(f) {f((staClients+transitions)*2)}
      tServ <- matrix(sapply(lapply(serviceDistribution, r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=nodes )
      tArr <- matrix(sapply(lapply(arrivalDistribution[[1]], r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=length(arrivalDistribution[[2]]), dimnames=list(1:((staClients+transitions)*2), arrivalDistribution[[2]]))
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
      rho <- l - lq
      eff <- w/(w-wq)
     
      nozero <- apply(prob, 1, function(v){any(v>0)})
      pn <- prob[1:max(which(nozero)),]/cron
      obj$out <- list(l=l, lq=lq, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      #obj$out$data <- array(c(l, lq, w, wq), dim= c(nodes, 4), dimnames=list(1:nodes, c("L", "Lq", "W", "Wq")))
      #obj$out$Lqt <- sum(obj$out$data[,1])
      if (historic) {
        obj$out$historic <- hist
      }
      oldClass(obj) <-  c("Open", "SimulatedNetwork", "SimulatedModel")
      
      return(obj)
  }
  if (nproc > 1) {
    cluster <- makeCluster(nproc)
    registerDoParallel(cluster, cores=nproc)
    parallelRes <- foreach(transitions=idiv(nsim, chunks = nproc), .combine=c, .multicombine=TRUE, .packages="distr") %dopar% {
      belong <- function(obj, packagename) {
        if (length(obj) > 1)
          packageobj <- sapply(lapply(obj, class), attr, "package")
        else
          packageobj <- attr(class(obj), "package")
        
        return(all(packageobj==packagename))
      }
      res <- list()
      for(i in 1:simulations) {
        res <- c(res, list(OpenNetwork_secuential(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic)))
      }
      return(res);
    }
    stopCluster(cluster)
    return(if (length(parallelRes)==1) parallelRes[[1]] else  parallelRes)
  }
  res <-list()
  for(i in 1:nsim)
    res <- c(res, list(OpenNetwork_secuential(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic)))
  return(if (length(res)==1) res[[1]] else res)
}

#' @rdname Pn
#' @method Pn SimulatedModel
#' @details
#' \code{Pn.SimulatedModel} implements the method for a SimulatedModel model
#' @export
Pn.SimulatedModel <- function(qm, n) {
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


#' Print the main characteristics of a SimulatedModel object
#' @param x SimulatedModel object
#' @param ... Further arguments passed to or from other methods.
#' @method print SimulatedModel
#' @keywords internal
#' @export
print.SimulatedModel <- function(x, ...) {
  cat("Model: ", class(x)[1])
  if (is.list(x$out$l)) {
    cat("\nL =\t", x$out$l$mean, "\u00B1", x$out$l$error, "\tW =\t", x$out$w$mean, "\u00B1", x$out$w$error, "\t\tIntensity =\t", x$out$rho$mean , "\u00B1", x$out$rho$error, "\n")
    cat("Lq =\t", x$out$lq$mean, "\u00B1", x$out$lq$error, "\tWq =\t", x$out$wq$mean, "\u00B1", x$out$wq$error, "\tEfficiency =\t", x$out$eff$mean, "\u00B1", x$out$eff$error, "\n\n")
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
  print(data.frame("L"=x$out$l, "Lq"=x$out$lq, "W"=x$out$w, "Wq"=x$out$wq))
}

# summary.SimulatedModel <- function(object, range=NULL, ...) {
#   if (is.null(object$out$historic)) stop("Argument 'historic' should be TRUE to show the plots")
#   if (is.null(range))  {
#     maxrange <- length(object$out$historic[,"L"])
#     minrange <- 1
#   }
#   else {
#     maxrange <- max(range)
#     minrange <- min(range)
#   }
#   par(mfrow=c(3,2))
#   barplot(object$out$historic[minrange:maxrange,"Clients"])
#   plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Intensity"], type="l")
#   plot(minrange:maxrange, object$out$historic[minrange:maxrange,"L"], type="l")
#   plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Lq"], type="l")
#   plot(minrange:maxrange, object$out$historic[minrange:maxrange,"W"], type="l")
#   plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Wq"], type="l")
# }

#' Shows a plot of the evolution of a variable during the simulation
#' 
#' @param object Simulated Model
#' @param minrange First client choosen in the simulation
#' @param maxrange Last client choosen in the simulation
#' @param var Variable to show graphic (L, Lq, W, Wq)
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @param depth Number of point printed in the plot
#' @export
summarySimple <- function(object, minrange, maxrange, var, graphics, depth) {UseMethod("summarySimple", object)}


#' @rdname summarySimple
#' @method summarySimple SimulatedModel
#' @details
#' \code{summarySimple.SimulatedModel} implements the function for an object of class SimulatedModel.
#' @export
summarySimple.SimulatedModel <- function(object, minrange, maxrange, var, graphics="ggplot2", depth=maxrange-minrange) {
  if (is.null(object$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  if (length(intersect(var, c("L", "Lq", "W", "Wq", "Clients", "Intensity"))) <= 0) stop("Argument 'var' must be any of the following values: 'L', 'Lq', 'W', 'Wq', 'Clients', 'Intensity'.")
  eval(parse(text=paste("data <- data.frame(t=object$out$historic[minrange:maxrange, 'tClient'], '", var, "'=object$out$historic[minrange:maxrange, '", var, "'])", sep="")))
  switch(graphics,
        "graphics" =  {eval(parse(text=paste("plot(data$t, data$", var, ", col='red', type='l')\n
                                              abline(v=object$Staclients, untf = FALSE, col='black')\n
                                              legend('bottomright', c('", var, "'), lty =c(1), col = c('red'), bty='t')\n
                                              title(main='Evolution of ", var, "')", sep="")))},
         "ggplot2" = {data <- melt(data, id.var="t")
                      eval(parse(text=paste("qplot(t, value, data=data, geom='line', colour=variable, 
                      main='Evolution of ", var, ".', ylab='", var , "') + geom_vline(xintercept=object$out$historic[object$Staclients, 'tClient'], linetype='dotdash') + scale_colour_discrete(name='') + theme(legend.position='none')", sep="")))})
}

#' @rdname summarySimple
#' @method summarySimple SimulatedNetwork
#' @details
#' \code{summarySimple.SimulatedNetwork} implements the function for an object of class SimulatedNetwork.
#' @export
summarySimple.SimulatedNetwork <- function(object, minrange, maxrange, var, graphics="ggplot2", depth=maxrange-minrange) {
  if (is.null(object$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  if (length(intersect(var, c("L", "Lq", "W", "Wq", "Clients", "Intensity"))) <= 0) stop("Argument 'var' must be any of the following values: 'L', 'Lq', 'W', 'Wq', 'Clients', 'Intensity'.")
  switch(graphics,
         "graphics" =  {eval(parse(text=paste("plot(data$t, data$", var, ", col='red', type='l')\n
                                              abline(v=object$Staclients, untf = FALSE, col='black')\n
                                              legend('bottomright', c('", var, "'), lty =c(1), col = c('red'), bty='t')\n
                                              title(main='Evolution of ", var, "')", sep="")))},
         "ggplot2" = {
                      cumPlot <- qplot(x=ifelse(is.na(aux<-object$out$historic[1, "tClient", minrange:maxrange]), 0, aux), y=ifelse(is.na(aux<-object$out$historic[1, var, minrange:maxrange]), 0, aux), geom="line", colour="red")
                      for(node in 2:(dim(object$out$historic)[1])) {
                        cumPlot <- cumPlot + geom_line(aes(x=ifelse(is.na(aux<-object$out$historic[node, "tClient", minrange:maxrange]), 0, aux), y=ifelse(is.na(aux<-object$out$historic[node, var, minrange:maxrange]), 0, aux)), colour=colours(distinct=TRUE)[node*2+1]) + theme(element_rect(colour = node*2+1))
                        #data <- melt(data[node,,], id.vars=c("t"))
                        #eval(parse(text=paste("qplot(t, value, data=data, geom='line', colour=node, 
                        #main='Evolution of ", var, ".', ylab='", var , "') + geom_vline(xintercept=object$out$historic[1, 'tClient', object$Staclients], linetype='dotdash') + theme(legend.position='none')", sep="")))
                      } 
                      cumPlot
                     }
  )
}

#' @rdname summarySimple
#' @method summarySimple list
#' @details
#' \code{summarySimple.list} implements the function for an object of class list
#' @export
summarySimple.list <- function(object, minrange, maxrange, var, graphics="ggplot2", depth=100) {
  index <- 1
  rawdata <- lapply(object, function(qm){
                        array <- getElement(getElement(qm, "out"), "historic")
                        if (is.null(array)) stop(simpleError("The variable 'historic' must be active to show the plot"))
                        res <- NULL
                        if (length(dim(array)) == 2) {
                          res <- data.frame(var=array[seq(minrange, maxrange, length.out=depth),var], t=array[seq(minrange, maxrange, length.out=depth), "tClient"], sta=rep(array[qm$Staclient, "tClient"], depth), sim = rep(index, depth))
                        } else {
                          index <- seq(minrange, maxrange, length.out=depth)
                          res <- data.frame(var=matrix(array[,var, index ], ncol=1), t=matrix(array[,"tClient", index], ncol=1), node=matrix(sapply(1:dim(array)[1], function(i){rep(i, length(index))}), ncol=1), sta=rep(array[1, "tClient", qm$Staclient], depth))
                        }
                        index <<- index + 1
                        return(res)
                        })
  data <- rawdata[[1]]
  for(i in rawdata[-1]) {
    data <- rbind(data, i)
  }
  if(is.null(data$node)) {
     ggplot(data, aes(x=t, y=var, group=sim)) + geom_point(size=0.001, colour=276, na.rm=TRUE) + geom_vline(xintercept=data$sta, linetype='dotdash', colour="red") + ggtitle(paste("Evolution of ", var, sep="")) + theme(legend.position="none")
  }
   else {
     ggplot(data, aes(x=t, y=var, colour=node, group=node)) + geom_point(size=0.001, aes(alpha=0.98), na.rm=TRUE) + stat_smooth(se=FALSE, na.rm=TRUE) + geom_vline(xintercept=data$sta, linetype='dotdash', colour="red") + ggtitle(paste("Evolution of ", var, sep="")) + theme(legend.position="none")
  }
}

# summarySimple.list <- function(object, minrange, maxrange, var, graphics="ggplot2") {
#   lcolors <- colors(distinct=TRUE)
#   switch(graphics,
#          "ggplot2" = {
#             cumplot <- summarySimple(object[[1]], minrange, maxrange, var, graphics)
#             aux<-10
#             for (i in object[-1]){
#                 actualplot <- summarySimple(i, minrange, maxrange, var, graphics)
#                 cumplot <- cumplot + geom_line() + 
#                            geom_line(mapping=actualplot$mapping, colour=lcolors[aux]) 
#                 aux <- aux + 10
#             }
#             cumplot
#          }
#   )
# }

# summarywwq <- function(object, minrange, maxrange, graphics="ggplot2") {
#   if (is.null(object$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
#   data <- data.frame(n=minrange:maxrange, "W"=object$out$historic[minrange:maxrange, "W"], "Wq"=object$out$historic[minrange:maxrange, "Wq"])
#   switch(graphics,
#          "graphics" =  {plot(data$n, data$W, col="red", type="l")
#                         lines(data$n, data$Wq, col="blue")
#                         legend("bottomright", c("W", "Wq"), lty =c(1,1), col = c("red", "blue"), bty="n")
#                         title(main=paste("Evolution of W and Wq", sep=""))},
#          "ggplot2" = {data <- melt(data, id.var="n")
#                       qplot(n, value, data=data, geom="line", colour=variable, 
#                             main="Evolution of W and Wq.", ylab="Mean time")})
# }