---
title: "AMS 268 Final Project Description"
author: Arthur Lui
date: "7 March 2016"
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

Explore the Gaussian Process model fused with a Generalized Double Pareto prior.

The GP model to be explored is 
$$
\begin{aligned}
  y|f(X) &= f(X) + \epsilon,~\epsilon \sim N(0,\sigma^2I)\\
  f(X)|d_1,...,d_p &\sim GP(0,\kappa(\cdot,\cdot))\\
  d_1,...,d_p &\sim GDP(a,b)\\
  \\
  \Rightarrow y|d_1,...,d_p &\sim N(0,\sigma^2I + K) \\
  d_1,...,d_p &\sim GDP(a,b)\\
\end{aligned}
$$
with $K_{ij} = \kappa(x_i,x_j)$, covariance function $\kappa(x_i,x_j) = \tau \exp\p{-\phi\sqrt{\sum_{k=1}^p (x_{ik}-x_{jk})^2 d_k^2}}$.


The density for the GDP($a,b$) is $f(x;a,b) = \frac{1}{2b}\p{1 + \frac{\abs{x}}{ab}}^{-a-1}$

