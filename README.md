## varlasso package
This repo includes an R package for fitting Bayesian Cector Autoregressive (VAR) models with optional shrinkage priors on coefficient estimates

Documentation for the package is here: [https://atsa-es.github.io/varlasso/](https://atsa-es.github.io/varlasso/)

This work extends the MARSS package (which provides similar models in a maximum likelihood setting). 
MARSS: [https://cran.r-project.org/web/packages/MARSS/index.html](https://cran.r-project.org/web/packages/MARSS/index.html)

And more material can be found on our Applied Time Series website, [https://atsa-es.github.io/](https://atsa-es.github.io/)

## Installation

To install the package, use 

``` r
# install.packages("remotes")
remotes::install_github("atsa-es/varlasso")
```

## Citations

A paper describing applications of these methods to VAR problems is in PeerJ,

Ward, E.J., K.N. Marshall, and M.D. Scheuerell. 2022. Regularizing priors for Bayesian VAR applications to large ecological datasets, *PeerJ*  

And the citation for this repository is

[![DOI](https://zenodo.org/badge/360700576.svg)](https://zenodo.org/badge/latestdoi/360700576)
