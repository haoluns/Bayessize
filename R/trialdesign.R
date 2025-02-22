
#' Trialdesign fitted linear models
#'
#' This function creates an object that represents the planning of a clinical trial, 
#' which requires vector input of a series of powers and sample sizes,
#' as well as the type I rate, the number of tests, and the alpha-spending function.
#'
#' @param powervec The vector of simluated power values.
#' @param sizevec The vector of simluated sample size values.
#' @param alpha The desired type I error rate.
#' @param ntest The number of tests.
#' @param type The type of the alpha-spending function, 1 for O'Briend-Fleming's design aand 2 for Pocock's design.
#'
#' @return An S3 object of class "trialdesign" that contains the fitted linear model between the sample size and squared drift paramters, as well as all the function inputs.
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
#' newdesign = trialdesign(powervec, sizevec, alpha,ntest,type)
#'
trialdesign <- function(powervec, sizevec, alpha,ntest,type){

  tt  =sapply(1:length(powervec),function(i){eta(alpha,1-powervec[i],ntest,type)})
  tt2 = tt^2
  model = lm(sizevec ~ 0 + tt2)

  out = list(model = model, powervec = powervec, sizevec = sizevec, sqdrift = tt2, alpha = alpha, ntest = ntest, type = type)
  class(out) = "trialdesign"
  return(out)
}

