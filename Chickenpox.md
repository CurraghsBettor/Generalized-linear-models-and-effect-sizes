The glm function and the multiple
================

## Estimating Risk Ratios and Risk Differences Using Regression (Naimi & Whitcomb, 2020)

Recently, I discover a paper wrote by Naimi and Whitcomb (2020), in
which the authors provided the mean to estimate odds ratio, relative
risk ratio and risk difference ratio with the $glm()$ function by
selecting the appropriate link function. They provided R codes to allow
reproduction of their results. Here, we will explore how to use those
different links with another concrete example.

## Odds ratio, Relative risk ratio and risk difference ratio

Traditionally estimated in epidemiological studies, we will discover
three estimates, that are, odds ratio, relative risk ratio and risk
difference, and how to obtain them using a Generalized linear models in
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

$\ OR = \frac{p_1/1-p_1}{p_2/1-p_2}$

It can be easily computed too from a contingency table, with the
following equation:

$\ OR = \frac{a/c}{b/d}$, where a and d denote the main-diagonal cells,
and b and c denote the off-diagonal cells.

A second estimate is the relative risk ratio, which is defined as the
ratio between two probabilities: the probability of the outcome given
the exposition (or in the exposed group) and the probability of the
outcome given a non-exposition. This estimate is depicted as follows:

$\ RR = \frac{p_1}{p_2}$

If we wish to estimate the RR of a particular disease in an exposed
group vs a non-exposed group. If the estimated RR is greater than 1, it
means that the outcome is more likely to occur given the factor, whereas
if less than 1, the outcome is less likely to occur given the factor. if
the estimated RR is 6, it means that the exposed group is 6 time more
likely to develop the disease compared the other group. This estimate is
relatively easy to understand and be interpreted compared to the odds
ratio. Indeed, a RR = 3 means that the outcome is 3 times more likely to
occur given the exposition relative to no exposition.

Risk difference

## Generalized linear models

Generalized linear models are special cases of model that are often use
in R in order to explore how a predictor variable (or more) predict the
to-be explained variable (). Then, we have to select the appropriate
distribution (e.g., are the data normally distributed ?), and a
corresponding link function. In the case of logistic regression model,
that models the relationship between a binary outcome and predictor
variables, the appropriate distribution is the binomial distribution and
the corresponding link function will dependent about the estimate of
interest for a main purpose. As we’re speaking about three different
estimates, we will focus on three link functions.

logit link: $\ log\frac{P(Y = 1)}{1 - P(Y = 1)}$

log link: $log(P)$

identity: $P$

## Create the data: Estimating Odds ratios, Risk Ratios and Risk Differences Using Regression (Wrensch et al., 1997)

For the main purpose, we will borrow the data from Wrensch and
colleagues. (1997). In this paper the authors estimated the strength of
association between an episode of chickenpox and the risk to be affected
by glioma that is a brain tumor. Overall they found that there is an
inverse link between the exposition to the Chickenpox and the outcome.
Let create the data in a contingency table, to appreciate the data.

``` r
DataM <- matrix(c(267, 348, 114, 66),
                nrow = 2, byrow = T); DataM
```

    ##      [,1] [,2]
    ## [1,]  267  348
    ## [2,]  114   66

Before perform the analysis, we need some functions that I’ve programmed
to render the analysis easier. The first one allows to take the output
of the logistic regression model, and estimated a corresponding odds
ratio or risk ratio along with the corresponding 95% confidence
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

The second one is dependent from the sandwich package. It has to be used
when perform a glm with a poisson distribution and a log link, as we
will see later and discribe why does not use it might results as
problematic.

``` r
effRR <- function(model) {
  library(sandwich)
  se <- sqrt(sandwich(model)[2,2])
  RR <- round(exp(model[2,1]), 2)
  Ll <- round(exp(sum(model[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*model[2,2])), 2)
  Ul <- round(exp(sum(model[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*model[2,2])), 2)
  cat("RR =", RR, "95%CI [", Ll, ",", Ul, "]")
}
```

The next function allows simply to extract the risk difference.

``` r
effRD <- function(model) {
  RD <- round(model[2,1], 2)
  Ll <- round(sum(model[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*model[2,2]), 2)
  Ul <- round(sum(model[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*model[2,2]), 2)  
}
```

The last function will allow to check whether the extracted estimates
correspond well to what we are looking for, by computed them from the
contingency table.

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

Well, now, we have to create a data frame in order to analyze it with a
logistic regression model, as follows. One can see that this not so
tough to create it. \[Developper\]

``` r
if(!exists("data_wide")) {
  id <- 1: sum(267, 114, 348, 66)
  exposition <- c(rep(1, 267), rep(0, 114), rep(1, 348), rep(0, 66))
  Case <- c(rep(1, sum(267, 114)), rep(0, sum(348, 66)))
  data_wide <- cbind(id, exposition, Case)
  data_wide <- as.data.frame(data_wide)
}
```

The first model use a binomial distribution with a logit function.

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
column $Estimate$ and the row $expoisition$, we see a $\beta_1$
parameter that equals -0.8115. We see that this parameter is
significant, thus rejecting the null hypothesis of null association
between the outcome and the exposition. This parameter reflects the
change in $log-odds$ for one-unit change in the predictor variable.
Then, the predictor variable is binary, yes/no, a shift from 0 to 1,
results in a change in the outcome of -0.8115. Therefore, the negative
estimate reveals that the non-exposition results as an outcome that is
more likely to occur whereas a shift to exposition results in an outcome
that is less likely to occur, thus an exposition results in a decrease
of the likelihood to develop the disease then, offering a protection
against the development of the disease. In the context of our example,
it means that a prior Chickenpox infection protects against the
likelihood to develop a glioma. To obtain the corresponding odds ratio
one has just to exponentiate the $\beta$ parameter Let now use the first
function that we’ve programmed to find the corresponding odds ratio
along with the corresponding confidence interval.

Before continuing what is the relationship between the $\beta$ parameter
and the odds ratio ? Some lines before, we have just said this parameter
reflect the $log-odds$ changes per 1-unit changes, then it is the
natural logarithm of the odds ratio (remember the formula of the link
function), as a consequence exponentiate the $\beta$ parameter give the
odds ratio:

$e^{\beta} = Odds ratio$

At this step you probably ask you how to obtain mathematically the
corresponding confidence interval:

$Ll = e^{\beta - 1.96*se}$, where 1.96 denotes the quantile of the
normal distribution and se denotes the standard error the $\beta$
parameter. The same logic is applied to obtain the upper limit:
$Ul = e^{\beta + 1.96*se}$ all that our function allows from the output
of our model.

``` r
effORRR(model1)
```

    ## OR/RR = 0.44 95% CI [ 0.32 , 0.63 ]

The returned odds ratio equals 0.44. We obtain as well the corresponding
95%CI. A confidence interval is quite tough to interpret,often
interpreted as the population parameter lies between both the lower and
upper limits, which is false. Rather, it offers well a reliable estimate
of the study precision and the evidential weight in favor the effect
(see Matthews, 2001).

``` r
library(finalfit)
explanatory = c("Case")
dependent = c("exposition")
or_plot(data_wide, dependent, explanatory)
```

    ## Attente de la réalisation du profilage...
    ## Attente de la réalisation du profilage...
    ## Attente de la réalisation du profilage...

    ## Warning: Removed 2 rows containing missing values (`geom_errorbarh()`).

![](Chickenpox_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

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

Unsurprisingly, the estimate is again significant and negative,
highlighting that a non-exposition results in an increased risk to
develop the disease. Exponentiate this estimate returns thus the
Relative risk ratio.

``` r
effORRR(model2)
```

    ## OR/RR = 0.69 95% CI [ 0.59 , 0.79 ]

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

We have previously seen how to obtain multiple estimates as a function
of the link function used with glm, to be sure that it offers reliable
estimate, we will run a function and use it to estimate those effect
sizes from the 2 by 2 contingency table previously coded.

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

We see that our models estimated well all our effect sizes, note that
I’ve programmed the mean to obtain the absolute risk difference, then
the returned value correspond well to the absolute value of the
estimated $\beta$ parameter with the third model.

Latter they investigate \[… \] Adjusted odds ratio. Until we do not have
any information about.

In this later section, I want to discuss about which estimate is the
best.

In some cases the explanatory variable is continuous rather than
dichotomous. For instance we want to investigate to which extent the
factor age is related then the odds ratio is computed.
