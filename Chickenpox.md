Odds ratio, Relative risk ratio, Risk difference and generalized linear
model
================

## Estimating Risk Ratios and Risk Differences Using Regression Model (Naimi & Whitcomb, 2020)

Recently, I have discovered a paper written by Naimi and Whitcomb (2020), in
which, these authors provided the mean to estimate the odds ratio (OR), the relative
risk ratio (RR), and the risk difference (RD) using the `glm` function and
by selecting the appropriate link function. They provided R codes that allow for the
reproduction of their results. Here, we will explore how to use these
different models with another concrete example.

## Odds ratio, Relative risk ratio and risk difference ratio

Traditionally estimated in epidemiological studies, let's consider 
three estimates: 1) the odds ratio, 2) the relative risk ratio, and 3) the risk
difference. This is possible to obtain these estimates using Generalized linear models in
R with the same function.

The OR expresses the ratio between the odd of a binomial
distributed outcome (e.g., a particular disease) given the exposition to
a given factor (or given a specific factor) over the odd of the outcome
given no-exposure. Wait, what’s an odd ? An odd expresses the
ratio between the probability of the outcome occurring $P(Y = 1)$ over
the probability of the outcome non-occurring $1 - P(Y = 1)$. As such, if
the odd is greater than 1, it means that the outcome is more likely to
occur, whereas if less than 1, the outcome is less likely to occur.
Therefore, if the odds ratio is greater than 1, the outcome is more
likely to occur given the exposition to a specific factor, whereas if
less than 1 the outcome is less likely to occur given the exposure.
The odds ratio can be computed as follows:

$$\ OR = \frac{p_1/1-p_1}{p_2/1-p_2}$$

It can be also computed from a 2*2 contingency table, with the
following equation:

$$\ OR = \frac{a/c}{b/d}$$ where $a$ and $d$ denote the main-diagonal cells,
and $b$ and $c$ denote the off-diagonal cells.

The second estimate is the RR, which is defined as the
ratio between two probabilities: the probability of the outcome given
the exposition (or in the exposed group) and the probability of the
outcome given the absence of exposition. This estimate is depicted as follows:

$$\ RR = \frac{p_1}{p_2}$$

We can obtain this estimate from a 2*2 contingency table as follows:

$$\ RR = \frac{a/(a+b)}{c/(c+d)}$$

For instance, we wish to estimate the RR ratio of a particular disease
in an exposed group vs in a non-exposed group. If the estimated RR is
greater than 1, it means that the outcome is more likely to occur given
the factor or the exposition, whereas if less than 1, the outcome is less likely to occur
given the factor or the exposition. For instance, if the RR = 3, that means that the outcome is 3 times
more likely to occur given the exposition relative to the absence of the exposition.

The last estimate that we're going to describe is the Risk
difference/Attributable risk. “RD or AR is defined as the difference in
risk of a condition such as a disease between an exposed group and an
unexposed group” (Kim, 2017). Simply! Contrary to the previous effect
size which express the $ratio$ between probabilities, it expresses 
the $difference$ as depicted with the following formula:

$$\ RD = p_1 - p_2$$

This estimate can be found from a 2*2 contingency table as follows:

$$\ RD = a/(a+b) - c/(c+d)$$

## Generalized linear models

Generalized linear models are special cases of models, often used
in R to explore how a predictor variable (or more) predicts a
to-be explained variable. Therefore, we have to select the appropriate
distribution (e.g., we ask ourselves whether the data are normally distributed, or binomially distributed, etc...), and a
corresponding link function. For the main purpose, we're going to model
the relationship between a binary outcome and a binary predictor
variable. The appropriate distribution is the binomial distribution and
the corresponding link function will dependent about the estimate of
interest. As we’re speaking about three different
estimates, we will focus on three link functions. (give defintions of the link function)

`logit`: $\ log\frac{P(Y = 1)}{1 - P(Y = 1)}$

`log`: $log(P)$

`identity`: $P$

## Program the data (Wrensch et al., 1997): Estimate Odds ratios, Risk Ratios, and Risk Differences Using Regression 

We're going to borrow the data published in Wrensch and
colleagues. (1997). In this case-control study, these authors estimated
the strength of association between an episode of chickenpox (disease
caused by the Varicella Zoster virus) and the risk to be affected by
glioma (a brain tumor). Overall they found an inverse link
between a prior infection with Chickenpox and the risk to develop a glioma (i.e., a prior infection decreases the probablity to develop a glioma). Let us create the data in a contingency table to better appreciate the data.

``` r
DataM <- matrix(c(267, 348, 114, 66),
                nrow = 2, byrow = T); DataM
```

    ##      [,1] [,2]
    ## [1,]  267  348
    ## [2,]  114   66

Before performing the analysis, we need some functions that I’ve programmed
to keep the analysis easier. The first one allows for taking the output
of the logistic or the log-binomial regression model, and estimate the OR or the RR along with the corresponding 95% confidence
intervals.

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

The next function allows simply for extracting the RD and the
corresponding 95% CI.

``` r
effRD <- function(model) {
    coefModel <- as.matrix(summary(model)$coefficients)
  RD <- round(coefModel[2,1], 2)
  Ll <- round(sum(coefModel[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*coefModel[2,2]), 2)
  Ul <- round(sum(coefModel[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*coefModel[2,2]), 2)
  cat("Risk difference = ", RD, "95%CI [", Ll, ",", Ul, "]\n")
}
```
Please note that we are working with only one categorical predictor variable, hence, all models will return raw estimates (i.e.,
un-adjusted estimates).

Well, now, we have to create a data frame in order to analyze the data with our models of interest. To briefly explain, with `id` we define the total number of observations. The next line allows for create a vector, in which, the number `1` is repeated 267 times (corresponding to the number of individuals that were exposed and that were cases), and `0` 114 times (corresponding to the number of individuals that were unexposed and cases) and, etc... So why separate them like this ? Because, when merging this vector with the following vector (i.e., the 267 + 114 participants that were both exposed and unexposed and that were cases and the 348 + 66 participants that were both exposed and unexposed and that were non-cases) will create a matrix in which the number of individuals that were exposed or not and that were cases or not will match (i.e., capturing the 4 combinations). 


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
variable predicts significantly the outcome, that we will assess using the `glm` function.

HO: There is no association between a history of Chickenpox and the risk
to develop a Glioma 

H1: There is an association between a history of
Chickenpox and the risk to develop a Glioma

H0: The history of Chickenpox does not predict the risk to develop the risk to develop a glioma

H1: The history of Chickenpox predicts the risk to develop a glioma

In term of our models, the paramater $\beta_1$ quantify the prediction
and the aim of modelling is to discover whether the prediction is significant.

HO: $\beta_1 = 0$

H1: $\beta_1 =/= 0$

The first model that we are building use a binomial distribution, 
with a `logit` link function, and is known as logistic regression.

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

Of course, at this step, we want to interpret our results. Let us take the
column `Estimate` and the row `exposition`, we see that the $\beta_1$
parameter equals -0.8115. We see that this parameter is
significant, rejecting the null hypothesis stating the null association between the exposition and the outcome. 
This parameter reflects the difference in the $log-odd$ of the outcome per one-unit change in the predictor variable.
Given that the predictor variable is binary, 0(no)/1(yes), a shift from 0 to
1, results in a change in the log-odd of -0.8115. Therefore, the
negative estimate reveals that outcome that is more likely to occur given the absence of exposition compared to the presence of the exposition. Put another way, a shift in the exposition means that the outcome is less likely to
occur. Thus, an exposition decreases the likelihood to
develop the disease (offering a protection against the development
of the disease). In the context of our example, it means that a prior
Chickenpox infection protects against the likelihood to develop a
glioma.  
What is the relationship between the $\beta$
parameter and the odds ratio ? Earlier, we have said that this
parameter reflects the difference in $log-odds$ per 1-unit change, and this estimate is the natural logarithm of the odds ratio (remember the
formula of the link function), as a consequence exponentiate the $\beta$
parameter gives you the odds ratio:

$$e^{\beta} = Odds ratio$$

At this step you could ask yourself how to mathematically obtain the
corresponding confidence interval. The follwoing line allows for etimate the lower limit:

$$Ll = e^{\beta - 1.96*se}$$

where -1.96 denotes the 2.5% percentile (therefore, 1.96 denotes the 97.5%
percentile) of the reduced centered Normal distribution, that is, the standard deviation
from the mean (0) where the area under the curve equals $\alpha$/2. 
$se$ denotes the standard error of the $\beta$ parameter.

The same logic is applied to obtain the upper limit:

$$Ul = e^{\beta + 1.96*se}$$

all that our first function enable to do.

``` r
effORRR(model1)
```

    ## OR/RR = 0.44 95% CI [ 0.32 , 0.63 ]

The returned odds ratio equals 0.44, which indicates that the odds of
developing a glioma given a prior infection to chickenpox is 0.44 higher
than the odd of developing a glioma given no prior infection. Thjerefore, an
exposition reduces the risk to develop the disease (i.e., because the OR
is below 1). This odds ratio expresses the strengh between the two variables (given the main issue and diffculty to interpret the odd, though confusing, I prefer to come back to the notion of risk without posit the rare-disease assumption).
Also, our function returns corresponding 95% confidence intervam. A confidence
interval is quite tough to interpret. Often interpreted as the
population parameter lies between both the lower and upper limits, this interpretation is false. Rather, it offers well a reliable estimate of the precision of a given study (e.g., the evidential weight in favor of an effect of interest; see Matthews, 2001). Moreover, in order to reject the null hypothesis, one has to check whether the confidence interval does not include 1. 
Howver, by focusing on whether the CI includes or not 1 in order or not to reject H0 seems problematic (e.g., Wasserstein et al., 2020). In my opinion, we could perform strong conclusions by checking whether the CI includes 1 if and only if
the width of the CI is small and quite close the point-estimate, because the narrower the CI, the more the evidential weight (Matthews, 2001).

Now, let us carry out another `glm` model. Rather than use
`link = logit`, we will use `link = log`, this model is known as the log-binomial regression.

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

Unsurprisingly, the estimate is again significant and negative. Contrary to the previous 
estimate which reveals a change in the $log-odds$ of the outcome for one-unit difference in the predictor variable, 
the current estimate reveals the change in the $log-risk$ per one-unit change. 
In the context of our exemple, it means that a shift from no
prior infection to a prior infection results in a difference in the $log-risk$ of the outcome of -0.37762, highlighting an increased risk to develop a glioma for patients without any episode of chickenpox. Exponentiate this estimate returns the Relative risk ratio.

``` r
effORRR(model2)
```

    ## OR/RR = 0.69 95% CI [ 0.59 , 0.79 ]

To interpret the exponentiate estimate, we may say that have had a chickenpox episode increased by 0.69 the risk to
develop a glioma compared to have had nothing. According to such a results, a prior chickenpox
infection protects against the risk to develop a glioma. A last step would be to quantify the
percentage of protection, provided as follows:

``` r
(1-0.69)*100
```

    ## [1] 31

So, the chickenpox infection protects by 31% against the risk to develop a glioma. 
Quite simple no ?

Let us now perform a new model in which we want to estimate the Risk
difference. In this case, the $\beta_1$ coeffciient returns this estimate.

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

Now let’s display the results with 95%CI:

``` r
effRD(model3)
```

    ## Risk difference =  -0.2 95%CI [ -0.28 , -0.12 ]

Now, let us come back to the RR ratio. As revealed in Naimi and Whitcomb
(2020), it appears that use the `log` link function with a binomial
distribution (for estimating log relative risk ratio) returns a error, explained by a fail a convergence (see Williamson et al., 2013; Huang, 2021). In this case, we should use a poisson regression model with a `log` link function.

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

To estimate RR, poisson regression models return correct $\beta$ estimates (though the distribution does not correspodn to what we have to use a priori), however, standard errors are not. Then, this method provide the mean to inscrease the risk to make type II errors (Huang, 2021) Therefore, as adviced by Naimi and Whitcomb (2020), we have to use the so-called robust sandwich variance estimator. The `sandwich` package allows for obtain correct standard errors (see also Huang, 2021 when data are analysed with generalized linear mixed effect models). Thus, using the following line, we obtain the correct standard error for $\beta_1$:

``` r
se <- sqrt(sandwich(model4)[2,2]); se
```

    ## [1] 0.07304571

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


``` r
coefM4 <- as.matrix(summary(model4)$coefficients)
RR <- round(exp(coefM4[2,1]), 2)
Ll <- round(exp(sum(coefM4[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = T)*se)), 2)
Ul <- round(exp(sum(coefM4[2,1], qnorm(0.025, mean = 0, sd = 1, lower.tail = F)*se)), 2)
cat("RR = ", RR, "95%CI [", Ll, ",", Ul, "]\n")
```

    ## RR =  0.69 95%CI [ 0.59 , 0.79 ]

The $z$ value correspoding to the estimate and therefore the $p$-value are incorrect, we can correct that as follows:

``` r
zvalue <- coefM4[2,1]/se
pvalue <- integrate(dnorm, lower = abs(zvalue), upper = Inf)$value
pvalue <- 2*pvalue
cat("Z-value =", zvalue, "p-value = ", pvalue, "\n")
```

    ## Z-value = -5.163067 p-value =  2.429361e-07


In other tutorials, I will, probably, explore how to perform such logistic
regressions with multiple predictor variable allowing for obtain adjusted odds/relative risk
ratio along different issues related to the current topic. Also, probably, I will explore how to perform similar analyses by using mixed effect models.
