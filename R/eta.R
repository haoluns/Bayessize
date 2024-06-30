#' Calculating Squared Drift Parameter
#'
#' This function computes the squared drift parameter given the type I and II error rates,
#' the number of tests, and the alpha-spending function.
#' The code is revised from the algorithm/code by Lai (2013).
#'
#' @param alpha The desired type I error rate.
#' @param beta The desired type II error rate.
#' @param ntest The number of tests.
#' @param type The type of the alpha-spending function, 1 for O'Briend-Fleming's design and 2 for Pocock's design.
#'
#' @import mvtnorm
#' @return The value of the drift parameter.
#' @export
#'
#' @examples
#' powervec=c(0.519,0.532,0.559,0.581,0.601,0.622,0.606,0.646,0.661,0.683,0.723,0.721,0.709,
#'     0.758,0.753,0.764,0.827,0.782,0.805,0.806,0.836,0.85,0.85,0.854,0.877,0.888,0.8835,
#'     0.884,0.898,0.9175,0.9025,0.911,0.922,0.9265,0.936,0.935,0.9465,0.935,0.943,0.948)
#' sizevec = seq(from=300,to=1080,by=20)
#' alpha = 0.1
#' ntest = 3
#' type =  2
#'
#' eta(alpha,1-powervec[1],ntest,type)
#'
#' @references Reference: Lai D (2013) Sample size determination for group sequential test under fractional Brownian motion. Electronic J of Stat 7, 1957–1967.
#'
eta<-function(alpha,beta,ntest,type){

  k = ntest
  alpha1 <- rep(0,k+1)
  alpha2 <- alpha1

  h  = 0.5

  covmatrix <- matrix(rep(0,(k+1)*(k+1)),ncol=k+1,byrow=T)
  tij <- (seq(0:k)-1)/k

  # the covariance matric of B(t)/sqrt(t)
  for (i in c(1:(k+1))) {
    for (j in c(1:(k+1))) {
      covmatrix[i,j] <- (1/2)*(tij[i]^(2*h)+tij[j]^(2*h)-(abs(tij[i]-tij[j]))^(2*h))/sqrt(tij[i]*tij[j])
    }
  }
  # print(covmatrix)
  covbm <- covmatrix
  # one sided
  for (i in c(1:k)) {
    alpha1[i+1] <- 2-2*pnorm(qnorm(1-alpha/2)/sqrt(i/k))
    alpha2[i+1] <- alpha*log(1+(exp(1)-1)*(i/k))
  }
  # print(alpha1)
  # print(alpha2)
  # compute the boundaries given the alpha value
  a1b <- rep(0,k)
  a2b <- a1b
  # one sided
  a1b[1] <- qnorm(1-alpha1[2])
  a2b[1] <- qnorm(1-alpha2[2])
  # fx1, one sided
  fx1 <- function(x,ub,covm,tprob) {
    # kn number of (upper) boundary already known
    kn <- length(ub)
    lb <- rep(-Inf,kn)
    pmv <- pmvnorm(lower=c(lb,x),upper=c(ub,Inf),sigma=covm)[1]
    tprob-pmv
  }
  a1ub <- NULL
  a2ub <- a1ub
  if (k > 1 ) {
    for (i in c(2:k)) {
      a1ub <- c(a1ub,a1b[i-1])
      a2ub <- c(a2ub,a2b[i-1])
      # print(a1ub)
      a1b[i] <- uniroot(fx1,interval=c(0.01,10),lower=0.01,upper=10,ub=a1ub,
                        covm=covmatrix[2:(i+1),2:(i+1)],tprob=alpha1[i+1]-alpha1[i])$root
      a2b[i] <- uniroot(fx1,interval=c(0.01,10),lower=0.01,upper=10,ub=a2ub,
                        covm=covmatrix[2:(i+1),2:(i+1)],tprob=alpha2[i+1]-alpha2[i])$root
    }
  }
  #print(a1b)
  #print(a2b)
  # search for the mean, hence the drift parameter theta=mean/sqrt(t)
  fxlai <- function(x,ub,covm,tprob) {
    # x is the drift parameter
    # kn number of (upper) boundary already known
    kn <- length(ub)
    lb <- rep(-Inf,kn)
    lmean <- x*sqrt(seq(1:kn)/kn)
    pmv <- pmvnorm(lower=lb,upper=ub,mean=lmean,sigma=covm)[1]
    tprob-pmv
  }
  # compute theta
  if (type==1){
    theta <- uniroot(fxlai,interval=c(0.01,10),lower=0.01,upper=10,ub=a1b,
                     covm=covmatrix[2:(k+1),2:(k+1)],tprob=beta)$root
  }

  if(type==2){
    theta <- uniroot(fxlai,interval=c(0.01,10),lower=0.01,upper=10,ub=a2b,
                     covm=covmatrix[2:(k+1),2:(k+1)],tprob=beta)$root
  }

  return(theta)
}

