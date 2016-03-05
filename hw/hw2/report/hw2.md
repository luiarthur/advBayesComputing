---
title: "Sequential Monte Carlo & Predictive Process (HW 2)"
author: Arthur Lui
date: "4 March 2016"
geometry: margin=1in
fontsize: 12pt
header-includes: 
    - \usepackage{bm}
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
