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

#' Obtains the main characteristics of a G/G/1 model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/1 model:
#' \item{pn}{Stores all the positives steady-state probabilities of having n customers, with n from 0 to staClients+nClients}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and  \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels

G_G_1 <- function(arrivalDistribution = Exp(1), serviceDistribution = Exp(1), staClients = 100, nClients = 1000, historic = FALSE) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) (staClients+nClients)
      tServ <- r(serviceDistribution) (staClients+nClients)
      iArr <- iServ <- 1
      sysClients <- simClients<- 0
      a <- 0
      b <- -1
      
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Nclients=nClients)
      if (historic) hist <- matrix(nrow=nClients, ncol=6, dimnames=list(1:nClients, c("L", "Lq", "W","Wq", "Clients", "Intensity")))
      
      while (simClients < staClients) {
        if (sysClients > 0)
          tMin <- min(a, b)
        else
          tMin <- a

        if (tMin == a) {
          simClients <- simClients + 1     
          if (sysClients == 0) {
            b <- tServ[iServ]
            iServ <- iServ+1
          } else
            b <- b - tMin
          a <- tArr[iArr]
          iArr <- iArr + 1
          sysClients <- sysClients+1
        } else {
          sysClients <- sysClients - 1
          if (sysClients == 0) {
            b <- -1
          } else {
            b <- tServ[iServ]
            iServ <- iServ + 1
          }
          a <- a - tMin
        }
      }
      simClients <- 0
      cron <- d <- c <- 0
      tnClients <- numeric(nClients)
      actualper <- seq(0, 100, 10)
      percentages <- (actualper*nClients)/100
      iper <- 1 
      while (simClients < nClients) {
        if (simClients > percentages[iper]) {
          iper <- iper+1
          message(paste(actualper[iper], "%", sep=""))
        }
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
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
  
          hist[simClients, ] <- c(l,lq,w,wq, sysClients, l-lq)
        }
      }
      l <- c/cron
      lq <- d/cron
      w <- c/nClients
      wq <- d/nClients
      eff <- w/(w-wq)
      rho <- l-lq
      pn <- tnClients[1:max(which(tnClients>0))]/cron
      if (historic)
        obj$out <- list(historic=hist, l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      else
        obj$out <- list(l=l, lq=lq, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
      oldClass(obj) <-  c("G_G_1", "SimulatedModel")
      
      return(obj)
}

exportToUI(G_G_1, "G/G/1",  c("distr", "distr", "numeric", "numeric", "boolean"),  c("G_G_1", "SimulatedModel"))

#' Obtains the main characteristics of a G/G/s model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter used to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/S model:
#' \item{pn}{vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S <- function (arrivalDistribution=Exp(1), serviceDistribution=Exp(1), s=2, staClients=100, nClients=1000, historic=FALSE) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (s <= 0) stop("Argument 's' must be greather than 0.")
      if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) (nClients+staClients)
      tServ <- r(serviceDistribution) (nClients+staClients)
      iArr <-1
      iServ <- 0
      bussyservs <- numeric()
      sysClients <- simClients<- 0
      cron <- d <- c <- 0
     
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, Nclients=nClients)
      tnClients <- numeric(nClients)
      if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
      
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
          iArr <- iArr + 1    
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
            sysClients <- sysClients + 1
            if (length(bussyservs) > 0)
              tServ[bussyservs] <- tServ[bussyservs] - tmin
          }
        }else {
          finishserver <- which.min(tServ[bussyservs])
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
      }
      simClients <- 0
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
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
          
          hist[simClients, ] <- c(l,lq,w,wq)
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

exportToUI(G_G_S, "G/G/s", c("distr", "distr", "numeric", "numeric", "numeric", "boolean"),  c("G_G_S", "SimulatedModel"))

#'Obtains the main characteristics of a G/G/1/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/1/K model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation.}
#' @export
#' @family SimulatedModels

G_G_1_K <- function (arrivalDistribution=Exp(1), serviceDistribution=Exp(1), K=2, staClients=100, nClients=1000, historic=FALSE) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (K <= 0) stop("Argument 'K' must be greather than 0.")
      if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
      tServ <- r(serviceDistribution) ((nClients+staClients)*2)
      iArr <- iServ <- 1
      sysClients <- simClients<- 0
      cron <- d <- c <- 0
      tnClients <- numeric(nClients)
      
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Nclients=nClients)
      if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
      
      while(simClients < staClients) {
        
        if (sysClients > 0) {
          tmin <- min(tArr[iArr], tServ[iServ])
        } else {
          tmin <- tArr[iArr]
        }
        
        
        if(tmin == tArr[iArr]) {
          iArr <- iArr + 1
          
          if (sysClients == (K+1)) {
            tServ[iServ] <- tServ[iServ] - tmin
            
          } else {
            simClients <- simClients + 1
            if (sysClients == 0)ciServ <- iServ + 1
            else tServ[iServ] <- tServ[iServ] - tmin
  
            sysClients <- sysClients + 1
          }
        } else {
          sysClients <- sysClients - 1
          
          if(sysClients != 0) {
            iServ <- iServ + 1
          }
          tArr[iArr] <- tArr[iArr] - tmin
        }
      }
      simClients <- 0
      
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
              sysClients <- sysClients + 1
              iServ <- iServ + 1
            } else {
              c <- c+tmin*sysClients
              d <- d+tmin*(sysClients-1)
              sysClients <- sysClients + 1
              tServ[iServ] <- tServ[iServ] - tmin
            }
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
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
          hist[simClients, ] <- c(l,lq,w,wq)
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

exportToUI(G_G_1_K, "G/G/1/K", c("distr", "distr", "numeric", "numeric", "numeric", "boolean"), c("G_G_1_K", "SimulatedModel"))

#' Obtains the main characteristics of a G/G/s/K model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param K Maximun size of the queue
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/S/K model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W}and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S_K <- function(arrivalDistribution=Exp(1), serviceDistribution=Exp(1), s=2, K=3, staClients=100, nClients=1000, historic=FALSE) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (s <= 0) stop("Argument 's' must be greather than 0.")
      if (K <= 0) stop("Argument 'K' must be greather than 0.")
      if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
      
      tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
      tServ <- r(serviceDistribution) ((nClients+staClients)*2)
      iArr <- iServ <- 1
      bussyservers <- rep(NA, s)
      sysClients <- simClients<- 0
      cron <- d <- c <- 0
      
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, Nclients=nClients)
      tnClients <- numeric(nClients)
      if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
      
      while(simClients < staClients) {
        indice_min_b <- which.min(bussyservers)
        
        if (sysClients > 0) 
          tmin <- min(tArr[iArr], bussyservers[indice_min_b])
        else
          tmin <- tArr[iArr]
        
        if (tmin == tArr[iArr]) {
          if(sysClients == (K+s)){
            iArr <- iArr + 1
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
              sysClients <- sysClients + 1
              iArr <- iArr + 1
              bussyservers <- bussyservers-tmin
            }
          }
        } else {
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
      }
      simClients <- 0
      
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
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
          
          hist[simClients, ] <- c(l,lq,w,wq)
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

exportToUI(G_G_S_K, "G/G/s/K",  c("distr", "distr", "numeric", "numeric", "numeric", "numeric", "boolean"),  c("G_G_S_K", "SimulatedModel"))

#' Obtains the main characteristics of a G/G/1/\eqn{\infty}/H model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/1/\eqn{\infty}/H model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels

G_G_1_INF_H <- function(arrivalDistribution=Exp(1), serviceDistribution=Exp(1), H=2, staClients=100, nClients=1000, historic=FALSE) {
    if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
    if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
    if (H <= 0) stop("Argument 'H' must be greather than 0.")
    if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
    if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
    
    tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
    tServ <- r(serviceDistribution) ((nClients+staClients)*2)
    possibleClients <- tArr[1:H]
    iServ <- 1
    iArr <- H+1
    sysClients <- 0
    simClients <- 0
    
    obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, H=H, Nclients=nClients)
    tnClients <- numeric(nClients)
    if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
    
    while(simClients < staClients) {
      indice_minimo_a <- which.min(possibleClients)
      
      if(sysClients>0)
        minimo<-min(possibleClients[indice_minimo_a], tServ[iServ])
      else
        minimo<-possibleClients[indice_minimo_a]
      
      if (length(indice_minimo_a) == 0 || minimo==tServ[iServ]) {
        sysClients <- sysClients-1
        if (sysClients==0)
          tServ[iServ] <- -1
        else
          iServ <- iServ + 1
        possibleClients<-possibleClients-minimo
        possibleClients[which(is.na(possibleClients))]<- tArr[iArr]
        iArr <- iArr+1
      } else {
        simClients<-simClients+1
        
        if(sysClients==0) {
          sysClients<-sysClients+1
          possibleClients <-possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
          iServ <- iServ + 1
        } else{
          sysClients<-sysClients+1
          tServ[iServ]<-tServ[iServ]-minimo
          possibleClients<-possibleClients-minimo
          possibleClients[indice_minimo_a]<- NA
        }
      }
    }
    
    simClients <- cron <- c <- d <- 0

    while(simClients < nClients) {
      indice_minimo_a <- which.min(possibleClients)
      
      if (sysClients>0)
        minimo <- min(possibleClients[indice_minimo_a],tServ[iServ])
      else
        minimo <- possibleClients[indice_minimo_a]
      
      if (length(indice_minimo_a) == 0 || minimo==tServ[iServ]) {
        cron <- cron + minimo
        c <- c + minimo*sysClients
        d <- d + minimo*(sysClients-1)
        tnClients[sysClients+1] <- tnClients[sysClients+1]+minimo
        sysClients<-sysClients-1
        if (sysClients==0)
          tServ[iServ] <- -1
        else
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
      if (historic) {
        l <- c/cron
        lq <- d/cron
        w <- c/nClients
        wq <- d/nClients
        
        hist[simClients, ] <- c(l,lq,w,wq)
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
    
    oldClass(obj) <-  c("G_G_1_INF_H", "SimulatedModel")
    return(obj)
}

exportToUI(G_G_1_INF_H, "G/G/1/INF/H", c("distr", "distr", "numeric", "numeric", "numeric", "boolean"), c("G_G_1_INF_H", "SimulatedModel"))

#' Obtains the main characteristics of a G/G/S/\eqn{\infty}/H  model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param s Number of servers
#' @param H Population size
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/S/\eqn{\infty}/H model
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels

G_G_S_INF_H <- function(arrivalDistribution=Exp(1), serviceDistribution=Exp(1), s=2, H=2, staClients=100, nClients=1000, historic=FALSE) {
      if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
      if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
      if (s <= 0) stop("Argument 's' must be greather than 0.")
      if (H <= 0) stop("Argument 'H' must be greather than 0.")
      if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
      if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
  
      tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
      tServ <- r(serviceDistribution) ((nClients+staClients)*2)
      possibleClients <- tArr[1:H]
      bussyservers <- rep(NA, s)
      iServ <- 1
      iArr <- H+1
      sysClients <- 0
      simClients <- 0
      nClicien <- 100/nClients
      obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, H=H, Nclients=nClients)
      tnClients <- numeric(nClients)
      if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
      
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
            sysClients <- sysClients + 1
            possibleClients <- possibleClients - tmin
            possibleClients[indice_minimo_a] <- NA
            bussyservers <- bussyservers - tmin
            bussyservers[which(is.na(bussyservers))[1]] <- tServ[iServ]
            iServ <- iServ + 1
          } else {
            sysClients <- sysClients + 1
            bussyservers <- bussyservers - tmin
            possibleClients <- possibleClients - tmin
            possibleClients[indice_minimo_a] <- NA
          }
        } else {
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
      }
      simClients <- cron <- c <- d <- 0
      progressbar <- txtProgressBar(min=0, max=nClients, style=3)
      auxProg <- !(((0:nClients)/nClients*100) %% 5)
      jump <- 10
      while (simClients < nClients) {
        if (auxProg[simClients+1]) setTxtProgressBar(progressbar, simClients)
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
        if (historic) {
          l <- c/cron
          lq <- d/cron
          w <- c/nClients
          wq <- d/nClients
          
          hist[simClients, ] <- c(l,lq,w,wq)
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
      oldClass(obj) <-  c("G_G_S_INF_H", "SimulatedModel")
      close(progressbar)
      return(obj)
}

exportToUI(G_G_S_INF_H, "G/G/s/INF/H", c("distr", "distr", "numeric", "numeric", "numeric", "numeric", "boolean"), c("G_G_S_INF_H", "SimulatedModel"))

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
#' @return
#' Returns the next information of a G/G/1/S/\eqn{\infty}/H/Y model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels
G_G_S_INF_H_Y <- function(arrivalDistribution=Exp(1), serviceDistribution=Exp(1), s=2, H=2, Y=3, staClients=100, nClients=1000, historic=FALSE) {
  if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
  if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
  if (s <= 0) stop("Argument 's' must be greather than 0.")
  if (H <= 0) stop("Argument 'H' must be greather than 0.")
  if (Y <= 0) stop("Argument 'Y' must be greather than 0.")
  if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
  if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
  
  tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
  tServ <- r(serviceDistribution) ((nClients+staClients)*2)
  possibleClients <- tArr[1:H]
#   bussyservers <- rep(-1, s)
  bussyservers <- rep(NA, s)
  iServ <- 1
  iArr <- H+1
  sysClients <- 0
  simClients <- 0
  
  obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Servs=s, H=H, Y=Y, Nclients=nClients)
  tnClients <- numeric(nClients)
  if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
  
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
  }
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
    if (historic) {
      l <- c/cron
      lq <- d/cron
      w <- c/nClients
      wq <- d/nClients
      
      hist[simClients, ] <- c(l,lq,w,wq)
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

exportToUI(G_G_S_INF_H_Y, "G/G/s/INF/H with Y replacements", c("distr", "distr", "numeric", "numeric", "numeric", "numeric", "numeric", "boolean"),  c("G_G_S_INF_H_Y", "SimulatedModel"))

#' Obtains the main characteristics of a G/G/\eqn{\infty} model by simulation
#' 
#' @param arrivalDistribution Arrival distribution
#' @param serviceDistribution Service distribution
#' @param staClients Number of customers used in stabilization stage
#' @param nClients Number of customers used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a G/G/\eqn{\infty} model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Expected number of customers in the system \eqn{L}}
#' \item{lq}{Expected number of customers in the queue \eqn{L_{q}}}
#' \item{w}{Expected waiting time in the system \eqn{W}}
#' \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation}
#' @export
#' @family SimulatedModels
G_G_INF <- function(arrivalDistribution=Exp(1), serviceDistribution=Exp(1), staClients=100, nClients=1000, historic=FALSE) {
  if (!belong(arrivalDistribution, "distr")) stop("Argument 'arrivalDistribution' must be a valid Class of the Distr package")
  if (!belong(serviceDistribution, "distr")) stop("Argument 'serviceDistribution'must be a valid Class of the Distr package")
  if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
  if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
  
  tArr <- r(arrivalDistribution) ((nClients+staClients)*2)
  tServ <- r(serviceDistribution) ((nClients+staClients)*2)
  iServ <- iArr <- 1
  sysClients <- 0
  simClients <- 0
  a <- 0
  b <- -1
  
  obj <- list(arrivalDistribution = arrivalDistribution, serviceDistribution=serviceDistribution, Nclients=nClients)
  tnClients <- numeric(nClients)
  if (historic) hist <- matrix(nrow=nClients, ncol=4, dimnames=list(1:nClients, c("L", "Lq", "W","Wq")))
  
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
    
    if (tmin == a) {
      simClients <- simClients + 1
      sysClients <- sysClients + 1
      a <- tArr[iArr]
      iArr <- iArr + 1
      b[indice_j] <- b[indice_j] - tmin
      b[indice_j_menos_uno[1]] <- tServ[iServ]
      iServ <- iServ + 1
    } else {
      indice_minimo <- which(b==tmin)
      sysClients <- sysClients - 1
      b[indice_j] <- b[indice_j] - tmin
      b[indice_minimo] <- -1
      a <- a - tmin
    }
  }
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
    if (historic) {
      l <- c/cron
      lq <- d/cron
      w <- c/nClients
      wq <- d/nClients
      
      hist[simClients, ] <- c(l,lq,w,wq)
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

exportToUI(G_G_INF, "G/G/INF", c("distr", "distr", "numeric","numeric", "boolean"),  c("G_G_INF", "SimulatedModel"))

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
#' @param p Routing matrix, where \eqn{p_{ij}} is the routing probability from node i to node j
#' @param staClients Number of customers used in the stabilization stage
#' @param nClients Number of customers in the system
#' @param transitions Number of transitions between nodes used in the simulation stage
#' @param historic Parameter to activate/deactivate the historic information
#' @return
#' Returns the next information of a Closed Network model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Vector of expected number of customers in the nodes \eqn{L}}
#' \item{lq}{Vector of expected number of customers in the queues of the nodes \eqn{L_{q}}}
#' \item{lqt}{Expected number of customers in the all queues}
#' \item{w}{Vector of expected waiting times in the nodes \eqn{W}}
#' \item{wq}{Vector of expected waiting times in the queues of the nodes \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq} during the simulation.}
#' @export
#' @family SimulatedModels

ClosedNetwork <- function(serviceDistribution, s, p, staClients, nClients, transitions, historic=FALSE) {
  if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution'must be a valid Class of the Distr package")
  if (any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
  if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
  if (nClients <= 0) stop("Argument 'nClients' must be greather than 0.")
  if (transitions <= 0) stop("Argument 'transitions' must be greather than 0.")
  
  ini <- proc.time()[3]
  nodes <- length(s)
  maxserv <- max(s)
  map <- function(f) {f(staClients+transitions)}
  tServ <- matrix(sapply(lapply(serviceDistribution, r), map, simplify=TRUE), nrow=((staClients+transitions)*2), ncol=nodes )
  iServ <- rep(1, nodes)
  rand <- r(Unif()) (staClients+transitions)
  simClients <- 0
  
  obj <- list(serviceDistribution=serviceDistribution, Servs= s, Prob=p, Nclients=nClients)
  if (historic) hist <- array(dim=c(nodes, 4, transitions), dimnames=list(1:nodes, c("L", "Lq", "W","Wq"), 1:transitions))
  
  sysClients <- floor((nClients-1)*(s^(-1))/sum(s^-1))
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
  ini <- proc.time()[3] - ini
  estab <- proc.time()[3]
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
    j <- 1
    while (acum < n_s) {
      acum <- acum + p[nodex, j]
      j <- j+1
    }
    j <- j-1
    
    sysClients[j] <- sysClients[j] + 1
    if (sysClients[j] <= s[j]) {
      indice_minimo_nodo_minimo <- which(is.na(b[,j]))
      b[indice_minimo_nodo_minimo[1], j] <- tServ[iServ[j], j]
      iServ[j] <- iServ[j] + 1
    }
  }
  simClients <- cron <- 0
  c <- d <- in_node <- rep(0, nodes)
  prob <- matrix(c(0), nrow=transitions, ncol=nodes)
  estab <- proc.time()[3] - estab
  simula <- proc.time()[3]
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
    
    cron <- cron + tmin
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
    if (historic) {
      l <- c/cron
      lq <- d/cron
      w <- c/in_node
      wq <- d/in_node
      hist[,,simClients] <- array(c(l, lq, w, wq), dim= c(nodes, 4))
    }
  }
  simula <- proc.time()[3] - simula
  l <- c/cron
  lq <- d/cron
  w <- c/in_node
  wq <- d/in_node
  lqt <- sum(lq)
  rho <- l - lq
  eff <- w/(w-wq)
  pn <- prob/cron
  obj$out <- list(l=l, lq=lq, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
  #obj$out$data <- array(c(l, lq, w, wq), dim= c(nodes, 4), dimnames=list(1:nodes, c("L", "Lq", "W", "Wq")))
  if (historic) {
    obj$out$historic <- hist
  }
  oldClass(obj) <-  c("Closed", "SimulatedNetwork", "SimulatedModel")
  
  return(obj)
}

exportToUI(ClosedNetwork, "Closed Network", c("Vdistr", "vector", "matrix", "numeric", "numeric", "numeric", "boolean"), c("Closed", "SimulatedNetwork", "SimulatedModel"))

SON_Example <- function (sta, trans) {
  p <- matrix(c(0.2, 0.25, 0.1, 0), nrow=2, ncol=2)
  s <- c(1, 2)
  arrivalDistribution <- pairlist(c(Exp(20), Exp(30)), c(1,2))
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
#' @return
#' Returns the next information of an Open network model:
#' \item{pn}{Vector of steady-state probabilities of having n customers in the system \eqn{P(n)}}
#' \item{l}{Vector of expected number of customers in the nodes \eqn{L}}
#' \item{lq}{Vector of expected number of customers in the queues of the nodes \eqn{L_{q}}}
#' \item{lqt}{Expected number of customers in all queues}
#' \item{w}{Vector of expected waiting times in the nodes \eqn{W}}
#' \item{wq}{Vector of expected waiting time in the queues of the nodes \eqn{W_{q}}}
#' \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
#' \item{rho}{Traffic intensity \eqn{\rho}}
#' \item{historic}{Optional parameter that stores the evolution of L, Lq, W and Wq during the simulation.}
#' @export
#' @family SimulatedModels

         
OpenNetwork <- function(arrivalDistribution, serviceDistribution, s, p, staClients, transitions, historic=FALSE) {
  if (!belong(serviceDistribution, "distr")) stop("All elements in argument 'serviceDistribution'must be a valid Class of the Distr package")
  if (!belong(arrivalDistribution[[1]], "distr")) stop("All elements in first position of argument 'serviceDistribution'must be a valid Class of the Distr package")
  if (any(s <= 0)) stop("All elements in argument 's' must be greather than 0.")
  if (staClients <= 0) stop("Argument 'staClients' must be greather than 0.")
  if (transitions <= 0) stop("Argument 'transitions' must be greather than 0.")
  
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
  
  obj <- list(serviceDistribution=serviceDistribution, Servs= s, Prob=p)
  if (historic) hist <- array(dim=c(nodes, 4, transitions), dimnames=list(1:nodes, c("L", "Lq", "W","Wq"), 1:transitions))
  
  a <- array(NA, dim=nodes)
  for(i in arrivalDistribution[[2]]) {
    a[i] <- tArr[iArr[i], as.character(i)]
    iArr[i] <- iArr[i] + 1
  }
  
 b <- matrix(NA, nrow=maxserv, ncol=nodes)
 sysClients <- rep(0, nodes)
  
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
      
      if (j > 0) {
        sysClients[j] <- sysClients[j] + 1
        if (sysClients[j] <= s[j]) {
          indice_minimo_nodo_minimo <- which(is.na(b[,j]))
          b[indice_minimo_nodo_minimo[1], j] <- tServ[iServ[j], j]
          iServ[j] <- iServ[j] + 1
        }
      }
    }
  }
  simClients <- 0
  cron <- 0
  c <- d <- entradas_nodo <- c(0, nodes)
  prob <- matrix(c(0), nrow=transitions, ncol=nodes)
  
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
  }
  l <- c/cron
  lq <- d/cron
  w <- c/entradas_nodo
  wq <- d/entradas_nodo
  lqt <- sum(lq)
  rho <- l - lq
  eff <- w/(w-wq)
  pn <- prob/cron
  obj$out <- list(l=l, lq=lq, lqt=lqt, w=w, wq=wq, pn=pn, rho=rho, eff=eff)
  #obj$out$data <- array(c(l, lq, w, wq), dim= c(nodes, 4), dimnames=list(1:nodes, c("L", "Lq", "W", "Wq")))
  #obj$out$Lqt <- sum(obj$out$data[,1])
  if (historic) {
    obj$out$Historic <- hist
  }
  oldClass(obj) <-  c("Open", "SimulatedNetwork", "SimulatedModel")
  
  return(obj)
}

exportToUI(OpenNetwork, "Open Network", c("Vdistr", "distr", "vector", "matrix", "numeric", "numeric", "boolean"),  c("Open", "SimulatedNetwork", "SimulatedModel"))

#' Print the main characteristics of a SimulatedModel object
#' @param x SimulatedModel object
#' @param ... Further arguments passed to or from other methods.
#' @method print SimulatedModel
#' @keywords internal
#' @export
print.SimulatedModel <- function(x, ...) {
  cat("Model: ", class(x)[1])
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$rho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
}

#' Print the main characteristics of a SimulatedNetwork object
#' @param x SimulatedNetwork object
#' @param ... Further arguments passed to or from other methods.
#' @method print SimulatedNetwork
#' @keywords internal
#' @export
print.SimulatedNetwork <- function(x, ...) {
  cat("Model: ", class(x)[1], "\n")
  cat("\nL =\t", x$out$l, "\tW =\t", x$out$w, "\t\tIntensidad =\t", x$out$rho , "\n")
  cat("Lq =\t", x$out$lq, "\tWq =\t", x$out$wq, "\tEficiencia =\t", x$out$eff, "\n\n")
  cat("LqT: ", x$out$lqt, "\n")
}

#' Shows the main graphics of the parameters of a Simulated Model
#' 
#' @param object Simulated Model
#' @param range Range of the graphics
#' @method summary SimulatedModel
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.SimulatedModel <- function(object, range=NULL, ...) {
  if (is.null(object$out$historic)) stop("Argument 'historic' should be TRUE to show the plots")
  if (is.null(range))  {
    maxrange <- length(object$out$historic[,"L"])
    minrange <- 1
  }
  else {
    maxrange <- max(range)
    minrange <- min(range)
  }
  par(mfrow=c(3,2))
  barplot(object$out$historic[minrange:maxrange,"Clients"])
  plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Intensity"], type="l")
  plot(minrange:maxrange, object$out$historic[minrange:maxrange,"L"], type="l")
  plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Lq"], type="l")
  plot(minrange:maxrange, object$out$historic[minrange:maxrange,"W"], type="l")
  plot(minrange:maxrange, object$out$historic[minrange:maxrange,"Wq"], type="l")
}

#' Shows a plot of the evolution of L and Lq during the simulation
#' 
#' @param object Simulated Model
#' @param minrange Start client of simulation
#' @param maxrange Last client of simulation
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @export
summaryllq <- function(object, minrange, maxrange, graphics="ggplot2") {
  if (is.null(object$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  data <- data.frame(n=minrange:maxrange, "L"=object$out$historic[minrange:maxrange, "L"], "Lq"=object$out$historic[minrange:maxrange, "Lq"])
  switch(graphics,
        "graphics" =  {plot(data$n, data$L, col="red", type="l")
                       lines(data$n, data$Lq, col="blue")
                       legend("bottomright", c("L", "Lq"), lty =c(1,1), col = c("red", "blue"), bty="n")
                       title(main=paste("Evolution of L and Lq", sep=""))},
         "ggplot2" = {data <- melt(data, id.var="n")
                      qplot(n, value, data=data, geom="line", colour=variable, 
                      main="Evolution of L and Lq.", ylab="Mean customers")})
}

#' Shows a plot of the evolution of W and Wq during the simulation
#' 
#' @param object Simulated Model
#' @param minrange Start client of simulation
#' @param maxrange Last client of simulation
#' @param graphics Type of graphics: "graphics" use the basic R plot and "ggplot2" the library ggplot2
#' @export
summarywwq <- function(object, minrange, maxrange, graphics="ggplot2") {
  if (is.null(object$out$historic)) stop("Argument 'historic' must be TRUE to show the plots")
  data <- data.frame(n=minrange:maxrange, "W"=object$out$historic[minrange:maxrange, "W"], "Wq"=object$out$historic[minrange:maxrange, "Wq"])
  switch(graphics,
         "graphics" =  {plot(data$n, data$W, col="red", type="l")
                        lines(data$n, data$Wq, col="blue")
                        legend("bottomright", c("W", "Wq"), lty =c(1,1), col = c("red", "blue"), bty="n")
                        title(main=paste("Evolution of W and Wq", sep=""))},
         "ggplot2" = {data <- melt(data, id.var="n")
                      qplot(n, value, data=data, geom="line", colour=variable, 
                            main="Evolution of W and Wq.", ylab="Mean time")})
}