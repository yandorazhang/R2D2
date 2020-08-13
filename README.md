R2D2: Bayesian linear regression using the R2-D2 shrinkage prior
====
  
## Overview
This package is designed for fitting linear regression model using a Bayesian global and local shrinkage prior, i.e., R2-D2 shrinkage prior. The R2-D2 prior contains two types: conditional R2-D2 and marginal R2-D2. The conditional R2-D2 is a prior conditioned on the the design matrix, while the marginal R2-D2 is free of the design matrix. In general, the marginal R2-D2 shrinkage prior is recommended. 


## R2D2 Installation

R2D2 software can be installed via Github. To install the latest version of R2D2 package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("yandorazhang/R2D2",force=T)
```


## Readme
To use the conditional R2-D2 prior, type `r2d2cond()`.

To use the marginal R2-D2 prior, type `r2d2marg()`.


Too use the Dirichlet-Laplace prior proposed in [Bhattacharya et al. 2015](https://www.tandfonline.com/doi/abs/10.1080/01621459.2014.960967), type `dl()`. 



## Citation

Please cite the following paper when you use R2D2:


[Zhang, Yan Dora, et al. "Bayesian Regression Using a Prior on the Model Fit: The R2-D2 Shrinkage Prior." arXiv:1609.00046v3 (2020)](https://arxiv.org/abs/1609.00046v3)



## Contact the Author
Author: Yan Dora Zhang, Brian P. Naughton, Howard D. Bondell, Brian J. Reich

Software Developer: Yuxi Cai, Yan Dora Zhang

Software Maintainer: Yan Dora Zhang (doraz@hku.hk or  yandorazhang@gmail.com)

