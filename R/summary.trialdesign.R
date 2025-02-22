#' Summary Method for the "Trialdesign" Object
#'
#' Summary method for class "trialdesign".
#'
#' @param object An object of class "trialdesign".
#'
#' @return Prints the number of tests, the type I error rate and the type of alpha-spending function of the trial design.
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
#' summary(newdesign)
#'
summary.trialdesign <- function(object){

  base::print(paste0("A Bayesian sequential design with ", ntest, "analyses and type I error rate of ", alpha))
  if (type == 1){
    spending = "O'Brien-Fleming's"
  }
  if (type == 2){
    spending = "Pocock's"
  }
  base::print(paste0(spending, "alpha-spending function"))

}

