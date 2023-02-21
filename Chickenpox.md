Odds ratio, Relative risk ratio, Risk difference and generalized linear
model
================

## Estimating Risk Ratios and Risk Differences Using Regression (Naimi & Whitcomb, 2020)

Recently, I discover a paper wrote by Naimi and Whitcomb (2020), in
which these authors provided the mean to estimate the odds ratio, the relative
risk ratio and the risk difference using the `glm` function and by
selecting the appropriate link function. They provided R codes to allow
reproduction of their results. Here, we will explore how to use those
different links with another concrete example.

## Odds ratio, Relative risk ratio and risk difference ratio

Traditionally estimated in epidemiological studies, we will discover
three estimates, that are, the odds ratio, the relative risk ratio and the risk
difference, and how to obtain them from a Generalized linear models in
R.

The odds ratio expresses the ratio between the odd of a binomial
distributed outcome (e.g., a particular disease) given the exposition to
a given factor (or given a specific factor) over the odd of the outcome
given no-exposure. Before continuing, what’s a odd. A odd expresses the
ratio between the probability of the outcome occurring $P(Y = 1)$ over
the probability of the outcome non-occurring $1 - P(Y = 1)$, then, if
the odd is greater than 1, it means that the outcome is more likely to
occur, whereas if less than 1, the outcome is less likely to occur.
Therefore, if the odds ratio is greater than 1, the outcome is more
likely to occur given the exposition to a specific factor, whereas if
less than one the outcome is less likely to occur given the exposure.
The odds ratio can be computed as follow(s):

$$\ OR = \frac{p_1/1-p_1}{p_2/1-p_2}$$

It can be easily computed too from a contingency table, with the
following equation:

$$\ OR = \frac{a/c}{b/d}$$ where a and d denote the main-diagonal cells,
and b and c denote the off-diagonal cells.

A second estimate is the relative risk ratio, which is defined as the
ratio between two probabilities: the probability of the outcome given
the exposition (or in the exposed group) and the probability of the
outcome given a non-exposition. This estimate is depicted as follows:

$$\ RR = \frac{p_1}{p_2}$$

One can obtain this estimate from a 2 by 2 contingency table via the
following formula:

$$\ RR = \frac{a/(a+b)}{c/(c+d)}$$

For instance, we wish to estimate the RR ratio of a particular disease
in an exposed group vs in a non-exposed group. If the estimated RR is
greater than 1, it means that the outcome is more likely to occur given
the factor, whereas if less than 1, the outcome is less likely to occur
given the factor. Indeed, a RR = 3 means that the outcome is 3 times
more likely to occur given the exposition relative to no exposition.

The last estimate that we will describe is the Risk
difference/Attributable risk. “RD or AR is defined as the difference in
risk of a condition such as a disease between an exposed group and an
unexposed group” (Kim, 2017). Simply. Contrary to the previous effect
size which express the $ratio$ between probabilities, in the current
case it expresses the $difference$ as depicted by the following formula:

$$\ RR = p_1 - p_2$$

One can obtain this estimate from a 2 by 2 contingency table via the
following formula:

$$\ RR = a/(a+b) - c/(c+d)$$

## Generalized linear models

Generalized linear models are special cases of model that are often use
in R in order to explore how a predictor variable (or more) predicts a
to-be explained variable. Then, we have to select the appropriate
distribution (e.g., are the data normally distributed ?), and a
corresponding link function. In the case of logistic regression,
that models the relationship between a binary outcome and predictor
variables, the appropriate distribution is the binomial distribution and
the corresponding link function will dependent about the estimate of
interest for a main purpose. As we’re speaking about three different
estimates, we will focus on three link functions.

`logit`: $\ log\frac{P(Y = 1)}{1 - P(Y = 1)}$

`log`: $log(P)$

`identity`: $P$

## Create the data: Estimating Odds ratios, Risk Ratios and Risk Differences Using Regression (Wrensch et al., 1997)

For the main purpose, we will borrow the data from Wrensch and
colleagues. (1997). In this case-control study, those authors estimated
the strength of association between an episode of chickenpox (disease
caused by the Varicella Zoster virus) and the risk to be affected by
glioma (a brain tumor). Overall they found that there is an inverse link
between a prior infection with Chickenpox and the risk to develop a glioma. Let's create the
data in a contingency table to appreciate the data.

``` r
DataM <- matrix(c(267, 348, 114, 66),
                nrow = 2, byrow = T); DataM
```

    ##      [,1] [,2]
    ## [1,]  267  348
    ## [2,]  114   66

Before perform the analysis, we need some functions that I’ve programmed
to keep the analysis easier. The first one allows to take the output
of the logistic regression model, and estimated a corresponding odds
ratio or a risk ratio along with the corresponding 95% confidence
interval.

``` r
effORRR <- function(model) {
  CoefSum <- coef(summary(model))
  data <- matrix(CoefSum, nrow = 2, ncol = 4)
  Ratio <- round(exp(data[2,1]), 2)
  Ll <- round(exp(sum(data[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*data[2,2])), 2)
  Ul <- round(exp(sum(data[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*data[2,2])), 2)
  cat("OR/RR =", Ratio, "95% CI [", Ll, ",", Ul, "]", "\n")
}
```

The next function allows simply to extract the risk difference and the
corresponding 95%CI.

``` r
effRD <- function(model) {
    coefModel <- as.matrix(summary(model)$coefficients)
  RD <- round(coefModel[2,1], 2)
  Ll <- round(sum(coefModel[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*coefModel[2,2]), 2)
  Ul <- round(sum(coefModel[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*coefModel[2,2]), 2)
  cat("Risk difference = ", RD, "95%CI [", Ll, ",", Ul, "]\n")
}
```

The last function will allow to check whether the extracted estimates
correspond well to what we are looking for, by computed them from the
contingency table (please note that we are working with only one
categorical predictor variable, hence, it returns raw estimates (i.e.,
un-adjusted estimates) that we are looking for.

``` r
ratio <- function(data) {
  OR <- ((data[1,1]/sum(data[1,1], data[1,2]))/(1-(data[1,1]/sum(data[1,1], data[1,2]))))/((data[2,1]/sum(data[2,1], data[2,2]))/(1-(data[2,1]/sum(data[2,1], data[2,2]))))
  RR <- (data[1,1]/sum(data[1,1], data[1,2]))/ (data[2,1]/sum(data[2,1], data[2,2]))
  data <- addmargins(data)
  data <- as.matrix(data)
  RD <- (data[1,1]/data[1,3]) - (data[2,1]/data[2,3]) 
  cat("OR = ", OR, "RR = ", RR, "RD = ", RD)
}
```

Well, now, we have to create a data frame in order to analyze it with a logistic regression model. One can see that this not so tough to create it. With `id` we just defined the total number of participants. The next line allows to create a vector, in which we repeated the number `1` 267 times (corresponding to a number of persons that were exposed and that were cases), `0` 114 times (corresponding to a number of persons that were unexposed and cases) and etc... So why separate them like this ? Because, when collapsing this vector with the following vector (i.e., the 267 + 114 participnats that were both exposed and unexposed and that were cases and the 348 + 66 participants that were both exposed and unexposed and that were non-cases) will create a matrix in which the number of persons that were exposed or not and that were cases or not match well (i.e., capturing well the 4 combinations). 


``` r
if(!exists("data_wide")) {
  id <- 1: sum(267, 114, 348, 66)
  exposition <- c(rep(1, 267), rep(0, 114), rep(1, 348), rep(0, 66))
  Case <- c(rep(1, sum(267, 114)), rep(0, sum(348, 66)))
  data_wide <- cbind(id, exposition, Case)
  data_wide <- as.data.frame(data_wide)
}
```

For the two following models we are searching whether there is an
association between the predictor variable and the to-be explained
variable. More particularly, we want to assess whether the predictor
variable predicts significantly the outcome, that we will assess using the `glm`.

HO: There is no association between a history of Chickenpox and the risk
to develop a Glioma 

H1: There is an association between a history of
Chickenpox and the risk to develop a Glioma

H0: The history of Chickenpox does not predict the risk to develop the risj to develop a glioma

H1: The history of Chickenpox predicts the risk to develop a glioma

In term of our model, the paramater $\beta_1$ quantify the prediction
and the aim is to discover whether it is significant.

HO: $\beta_1 = 0$

H1: $\beta_1 =/= 0$

The first model that we are building for use a binomial distribution
with a `logit` link function.

``` r
model1 <- glm(Case~exposition, family = binomial(logit), data = data_wide, control = list(trace = TRUE)); summary(model1)
```

    ## Deviance = 1079.341 Iterations - 1
    ## Deviance = 1078.449 Iterations - 2
    ## Deviance = 1078.448 Iterations - 3
    ## Deviance = 1078.448 Iterations - 4

    ## 
    ## Call:
    ## glm(formula = Case ~ exposition, family = binomial(logit), data = data_wide, 
    ##     control = list(trace = TRUE))
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.417  -1.067  -1.067   1.292   1.292  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.5465     0.1547   3.534  0.00041 ***
    ## exposition   -0.8115     0.1748  -4.643 3.43e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1100.7  on 794  degrees of freedom
    ## Residual deviance: 1078.4  on 793  degrees of freedom
    ## AIC: 1082.4
    ## 
    ## Number of Fisher Scoring iterations: 4

At this step, we want to interpret our results of course. Let takes the
column `Estimate` and the row `expoisition`, we see a $\beta_1$
parameter that equals -0.8115. We see that this parameter is
significant, thus rejecting the null hypothesis stating that the exposition does not 
predict the outcome. In fact, this parameter reflects the
change in $log-odds$ for one-unit change in the predictor variable.
Then, given the fact that the predictor variable is binary, 0(no)/1(yes), a shift from 0 to
1, results in a change in the log-odds of -0.8115. Therefore, the
negative estimate reveals that the non-exposition results as an outcome
that is more likely to occur compared to the exposition. Put another
way, a shift in the exposition results in an outcome that is less likely to
occur, thus an exposition decreases the likelihood to
develop the disease, thus, offering a protection against the development
of the disease. In the context of our example, it means that a prior
Chickenpox infection protects against the likelihood to develop a
glioma. Before continuing what is the relationship between the $\beta$
parameter and the odds ratio ? Some lines before, we have just said this
parameter reflect the $log-odds$ changes per 1-unit changes, then this
estimate is the natural logarithm of the odds ratio (remember the
formula of the link function), as a consequence exponentiate the $\beta$
parameter gives you the odds ratio:

$$e^{\beta} = Odds ratio$$

At this step you probably ask yourself how to obtain mathematically the
corresponding confidence interval:

$$Ll = e^{\beta - 1.96*se}$$

where -1.96 denotes the 2.5% percentile (hence 1.96 denotes the 97.5%
percentile) of the normal distribution, that is the standard deviation
from the mean (0) where the area under the curve equals $\alpha$/2 and
$se$ denotes the standard error the $\beta$ parameter.

The same logic is applied to obtain the upper limit:

$$Ul = e^{\beta + 1.96*se}$$

all that our first function allows from the output of our model.

``` r
effORRR(model1)
```

    ## OR/RR = 0.44 95% CI [ 0.32 , 0.63 ]

The returned odds ratio equals 0.44, which indicates that the odds of
developing a glioma given a prior infection to chickenpox is 0.44 higher
than the odd of developing a glioma given no prior infection, then an
exposition reduces the risk to develop the disease (i.e., because the OR
is below 1). This odds ratio expresses the strengh between the two variables.
We obtain as well a corresponding 95%CI. A confidence
interval is quite tough to interpret, often interpreted as the
population parameter lies between both the lower and upper limits, which
is false, rather, it offers well a reliable estimate of the study
precision and the evidential weight in favor of an effect of interest
(see Matthews, 2001). Moreover, in order to reject the null hypothesis, 
one has to check whether the confidence interval does not include 1. 
It is however problematic to focus on whether the CI includes or not 1 regarding our statistical 
hypothesis (e.g., Wasserstein et al., 2020). I'm pretty sure that it is because... personally,
I think that it is potentially inetresting to check that whenether both the lower and the upper limits are close to the point estimate.

Now, let us carry out another `glm` model, rather than use
`link = logit`, we will use `link = log`.

``` r
model2 <- glm(Case~exposition, family = binomial("log"), data = data_wide); summary(model2)
```

    ## 
    ## Call:
    ## glm(formula = Case ~ exposition, family = binomial("log"), data = data_wide)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.417  -1.067  -1.067   1.292   1.292  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.45676    0.05671  -8.054 8.02e-16 ***
    ## exposition  -0.37762    0.07304  -5.170 2.35e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1100.7  on 794  degrees of freedom
    ## Residual deviance: 1078.4  on 793  degrees of freedom
    ## AIC: 1082.4
    ## 
    ## Number of Fisher Scoring iterations: 5

Unsurprisingly, the estimate is again significant and negative. The
returned value is different from the previous one. Indeed, we have
estimated a log risk. In the context of our exemple, it highlights that
having no episode of chickenpox results in an increased risk to develop
a glioma. Exponentiate this estimate returns thus the Relative risk
ratio.

``` r
effORRR(model2)
```

    ## OR/RR = 0.69 95% CI [ 0.59 , 0.79 ]

Thus, have had a chickenpox episode increased by 0.69 the risk to
develop a glioma compared to have had nothing. Thus, a chickenpox
infection protects against the risk, hence, we want to quantify its
percentage of protection:

``` r
(1-0.69)*100
```

    ## [1] 31

So, according to the relative risk ratio, the chickenpox infection
protects by 31% against the risk to develop a glioma. Quite simple no ?

Let’s now perform a new model by which we want to estimate the Risk
difference:

``` r
model3 <- glm(Case ~ exposition, family = binomial("identity"), data = data_wide); summary(model3)
```

    ## 
    ## Call:
    ## glm(formula = Case ~ exposition, family = binomial("identity"), 
    ##     data = data_wide)
    ## 
    ## Deviance Residuals: 
    ##    Min      1Q  Median      3Q     Max  
    ## -1.417  -1.067  -1.067   1.292   1.292  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  0.63333    0.03592  17.633  < 2e-16 ***
    ## exposition  -0.19919    0.04110  -4.846 1.26e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 1100.7  on 794  degrees of freedom
    ## Residual deviance: 1078.4  on 793  degrees of freedom
    ## AIC: 1082.4
    ## 
    ## Number of Fisher Scoring iterations: 2

and let’s display the results with 95%CI:

``` r
effRD(model3)
```

    ## Risk difference =  -0.2 95%CI [ -0.28 , -0.12 ]

Now, we come back to the RR ratio, as revealed in Naimi and Whitcomb
(2020), it appears that use the `log` link function with a binomial
distribution (for estimating log relative risk ratio) can be accompanied
by an error. In that case it is necessary to use a poisson distribution
and the `log` link function.

``` r
library(sandwich)
model4 <- glm(Case~exposition, data = data_wide, family = poisson("log")); summary(model4)
```

    ## 
    ## Call:
    ## glm(formula = Case ~ exposition, family = poisson("log"), data = data_wide)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.1255  -0.9318  -0.9318   0.7328   0.7328  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -0.45676    0.09366  -4.877 1.08e-06 ***
    ## exposition  -0.37762    0.11188  -3.375 0.000738 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 560.48  on 794  degrees of freedom
    ## Residual deviance: 549.70  on 793  degrees of freedom
    ## AIC: 1315.7
    ## 
    ## Number of Fisher Scoring iterations: 5

However, because the distribution is not what one has to use, it
returned incorrect standard errors. Therefore, one has to use the robust
sandwich variance estiamtor using the `sandwich` package as follows:

``` r
se <- sqrt(sandwich(model4)[2,2]); se
```

    ## [1] 0.07304571

``` r
coefM4 <- as.matrix(summary(model4)$coefficients)
RR <- round(exp(coefM4[2,1]), 2)
Ll <- round(exp(sum(coefM4[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*se)), 2)
Ul <- round(exp(sum(coefM4[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*se)), 2)
cat("RR = ", RR, "95%CI [", Ll, ",", Ul, "]\n")
```

    ## RR =  0.69 95%CI [ 0.59 , 0.79 ]

Another mean to obtain the correct standard error:

``` r
vcov_sandwich <- vcovHC(model4, type = "HC1"); vcov_sandwich
```

    ##              (Intercept)   exposition
    ## (Intercept)  0.003224486 -0.003224486
    ## exposition  -0.003224486  0.005349133

``` r
se <- sqrt(vcov_sandwich[2,2]); se
```

    ## [1] 0.07313777

We have just highlighted a main issue with model 4, as a consequence the
$z$ value and the corresponding $p$-value are both incorrect, therefore one
can correct it as follows. Note that I use the `integrate` function since the
$p$-value reflects the area under the curve.

``` r
zvalue <- coefM4[2,1]/se
pvalue <- integrate(dnorm, lower = abs(zvalue), upper = Inf)$value
pvalue <- 2*pvalue
cat("Z-value =", zvalue, "p-value = ", pvalue, "\n")
```

    ## Z-value = -5.163067 p-value =  2.429361e-07

We have previously seen how to obtain multiple estimates as a function
of the link function used with `glm`, to be sure that it offers the real
estimate that we wanted, we will run a function and use it to estimate
those effect sizes from the 2 by 2 contingency table previously coded.

``` r
ratio <- function(data) {
  OR <- ((data[1,1]/sum(data[1,1], DataM[1,2]))/(1-(data[1,1]/sum(data[1,1], data[1,2]))))/((data[2,1]/sum(data[2,1], data[2,2]))/(1-(data[2,1]/sum(data[2,1], data[2,2]))))
  RR <- (data[1,1]/sum(data[1,1], DataM[1,2]))/ (data[2,1]/sum(data[2,1], data[2,2]))
  data <- addmargins(data)
  data <- as.matrix(data)
  ARD <- (data[2,1]/data[2,3]) - (data[1,1]/data[1,3])
  cat("OR = ", OR, "RR = ", RR, "ARD = ", ARD)
}
```

``` r
ratio(DataM)
```

    ## OR =  0.4441924 RR =  0.6854942 ARD =  0.199187

We see that our models estimated well all our effect sizes.

In other tutorials I will explore how to perform such logistic
regression models with a continuous predictor variable and with multiple
predictor variable, hence allowing to obtain adjusted odds/relative risk
ratio (as much as Data from Wrensch et al., 1997 were adjusted on age).
As well I will explore how to perform different tests to assess the
association between the variables (i.e., different tests that allow to
assess whether $\beta =/= 0$).
