---
title: "Sequential Monte Carlo & Predictive Process (HW 2)"
author: Arthur Lui
date: "4 March 2016"
geometry: margin=1in
fontsize: 12pt
header-includes: 
    - \usepackage{bm}
    - \pagestyle{empty}
    - \newcommand{\norm}[1]{\left\lVert#1\right\rVert}
    - \newcommand{\p}[1]{\left(#1\right)}
    - \newcommand{\bk}[1]{\left[#1\right]}
    - \newcommand{\bc}[1]{ \{#1\} }
    - \newcommand{\abs}[1]{ \left|#1\right| }
    - \newcommand{\mat}{ \begin{pmatrix} }
    - \newcommand{\tam}{ \end{pmatrix} }
---

# Probit Model
Note that for the probit model with 400 covariates, the $95\%$ credible intervals for the truly zero-valued coefficients contain zero $89\%$ of the time. And the truly zero-valued coefficients are all captured. Figure 1 shows the credible intervals and true values of 10 of the 400 covariates. The first six variables are truly non-zero, while the other 394 parameters are truly 0. Only the first 20 covariates are plotted, but the trend is similar for the remaining variables. 

| |$\beta\ne0$|$\beta=0$|
|:---:|:---:|:---:|
|$\hat\beta\ne 0$|6|42|
|$\hat\beta=   0$|0|352|

Table: This table summarizes how often the 95% credible intervals of the probit model coefficients contain 0. The credible intervals for the coefficients that are truly zero always exlude 0. But they never contain the true value of the coefficients, as seen in Figure 1. 

![Probit model 95% credible intervals. Plotted are the credible intervals and true values of 10 of the 400 covariates. The first six variables are truly non-zero, while the other 394 parameters are truly 0. ](../code/R/output/probit.pdf)

The coverage for the non-zero coefficients is 0. That is, they never contain the true value of the parameter. and the average length of the 95% credible intervals is 1.52 overall. The average 95% credible interval for the non-zero coefficients is 1.44. 

# Sequential Monte Carlo (SMC)
The SMC performs poorly at recovering the true model with $p=400$. Only one sample is produced in the sampler. The average length of the credible intervals is 0, and the coverage is 0.

# Predictive Process

The predictive process was fit to the simulated data with 1000 observations and 50 knots chosen at random in the three-dimensional location space. The priors used for the model were $\sigma^2\sim$ Inverse-Gamma$(2,.5)$, $\phi\sim$ Unif$(.1,3)$, and $\tau^2\sim$ Inverse-Gamma$(2,5)$. That is only $\phi$ was given an informative prior, and can be considered as fixed. In fact, with the prior, $\phi$ converged to the lower bound of the prior (.12). The priors for $\sigma^2$ and $\tau^2$ were chosen such that they have infinite priors, and mean at their scale parameters. The predictive process obtains for $\sigma^2$'s 95% credible interval (.54,.64), which nearly contains the true value of the parameter (.5). (See Figure 2)

![Posterior for hyperparameters](../code/R/output/gpPost.pdf)

After ordering the data by the magnitude of the true mean function, we see that the estimated mean function follows the mean function closely (Figure 2). The residuals sorted by the mean function, are plotted in Figure 3. When the mean function returns extreme values, the residuals are greater. This is probably a result of the estimated mean function interpolating towards the mean outside of the fitted grid-points.

![Mean function evaluated at the 1000 grid locations, sorted by the value of the mean function.](../code/R/output/gpOrderedData.pdf)

![Residuals sorted by the mean function.](../code/R/output/gpResid.pdf)

Finally, the mean function is plotted in Figure 4. The estimated mean function follows the true mean function closely. The variability of the mean function is small (close .01). This suggests that the predictive model is good at estimating the mean function with little uncertainty for 1000 observations evaluated at 50 knots.

![An awesome 3D plot. Top left: estimated mean function. Top right: true mean function. Bottom left: posterior mean-squared loss. Bottom right: posterior standard deviation.](../code/R/output/plot3d.pdf)

