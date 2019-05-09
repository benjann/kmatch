# kmatch
Multivariate-distance and propensity-score matching, including entropy balancing, 
inverse probability weighting, (coarsened) exact matching, and regression adjustment

`kmatch` matches treated and untreated observations with respect to covariates
and, if outcome variables are provided, estimates treatment effects based on
the matched observations, optionally including regression adjustment
bias-correction. Multivariate (Mahalanobis) distance matching as well as
propensity score matching is supported, either using kernel matching, ridge
matching, or nearest-neighbor matching. For kernel and ridge matching, several
methods for data-driven bandwidth selection such as cross-validation are
offered. In addition, several alternative matching and reweighting methods are
supported (coarsened exact matching, inverse probability weighting, entropy
balancing). The package also includes various commands for evaluating balancing
and common-support violations.

Stata version 11 or newer as well as `moremata` and `kmatch` is required. 

To install `kmatch` in Stata, type

    . ssc install kmatch, replace
    . ssc install moremata, replace
    . ssc install kdens, replace

