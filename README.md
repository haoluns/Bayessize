# Bayessize
 An R package for computing the sample size in Bayesian sequential trials based on fitting linear model between the sample size and squared drift paramters.

Yueyang Han, Haolun Shi, Ruitao Lin, Thomas Murray

# Example

To start calibrating the sample size of a Bayesian sequential design, we begin with a series of sample sizes and their 
corresponding simulated power values. Assuming that under the type I error rate of 0.1, number of interim tests of 3, and 
under the Pocock's design, the power values under sample sizes 300,320, ..., 1080 are as follows.

``` {.r language="R"}
 powervec=c(0.519,0.532,0.559,0.581,0.601,0.622,0.606,0.646,0.661,0.683,0.723,0.721,0.709,
    0.758,0.753,0.764,0.827,0.782,0.805,0.806,0.836,0.85,0.85,0.854,0.877,0.888,0.8835,
    0.884,0.898,0.9175,0.9025,0.911,0.922,0.9265,0.936,0.935,0.9465,0.935,0.943,0.948)
sizevec = seq(from=300,to=1080,by=20)
alpha = 0.1
ntest = 3
type =  2

```

We use function trialdesign() to create an object of class "trialdesign", which that represents the planning of a clinical trial, 
which requires vector input of a series of powers and sample sizes, as well as the type I rate, the number of tests, and the alpha-spending function. We may then use the plot(), predict() and summary() methods on the created object.

``` {.r language="R"}
newdesign = trialdesign(powervec, sizevec, alpha,ntest,type)
```

The summary() method prints the number of tests, the type I error rate and the type of alpha-spending function of the trial design. The plot() method creates a scatterplot of the relationship between the sample size and the squared drift parameter. The predict() method computes the sample size needed given the target power (in this case 0.8).
``` {.r language="R"}
summary(newdesign)
plot(newdesign)
predict(newdesign,0.8)
```