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


To install `kmatch` from the SSC Archive, type

    . ssc install kmatch, replace

in Stata. Stata version 11 or newer is required. Furthermore, the the `moremata` and `kdens` 
packages are required. To install these packages from the SSC Archive, type

    . ssc install moremata, replace
    . ssc install kdens, replace

---

Installation from GitHub:

    . net install kmatch, replace from(https://raw.githubusercontent.com/benjann/kmatch/master/)

---

Main changes:

    09mar2020
    - survey estimation is now supported through options -svy- and -subpop()-
    - ifgen() now stores the IFs even if nose is specified
    - options that generate variables are no longer allowed with vce(bootstrap) 
      or vce(jackknife)

    19jan2020
    - parsing of variable list failed if parentheses were used in factor variable
      specifications; this is fixed
    - options were allowed within outcome equations, but these options were 
      ignored; error is now returned if options are specified within outcome
      equations
    - vce() was only allowed if outcome variables have been specified; vce(analytic)
      and vce(cluster ...) are now also allowed without outcome variables
    - the strata variable stored by generate() (_KM_strata) was also filled in
      for observation outside of the estimation sample; this is fixed, i.e. the 
      variable is now set to missing for these observations

    30jul2019
    - in case of weighted data, balancing weights returned by -kmatch eb- were scaled
      in terms of sample size instead of sum of sampling weights; this did not
      affect treatment effect estimation, but lead to erroneous balancing 
      diagnostics in case of ATE; this is fixed

    29may2019
    - in case of multiple outcome equations, only the first equation was displayed
      in the output header; this is fixed
    - outcome equations in the header are now numbered
    - in case of duplicate outcomes only the duplicates were prefixed by a number 
      in the coefficient vector; this has been changed; now all outcomes receive a 
      prefix if there are duplicates
    - improved documentation (more examples)

    08may2019:
    - kmatch now computes approximate standard errors based on influence functions
      (assuming the matching weights to be fixed); corresponding vce() is -analytic-
      (default) or -cluster clustvar-; option -nose- suppresses SE estimation; option
      ifgenerate() stores the influence functions
    - option -comsup- without arguments now restricts obs to minimum PS range;
      returns by comsup changed
    - option -wor- added (nearest-neighbor matching without replacement)
    - option -keepall- added
    - coarsened exact matching now supported
    - new ebalance option to apply entropy balancing after matching
    - new -kmatch eb- command for entropy balancing
    - new -kmatch ipw- command for Inverse Probability Matching
    - new -kmatch em- command for Exact Matching (just for convenience)
    - new -kmatch ra- command for Regression adjustment (just for convenience)
    - dy() now supported by all subcommands
    - generate() now stores an additional variable containing strata ID (matching 
      subcommands only)
    - caliper() is now allowed as a synonym for bwidth()
    - new -idgenerate()- option to store IDs of matched controls
    - new -bwadjust()- option to adjust bandwidth by specified factor
    - new -maxiter()- option to restrict the maximum number of iterations for propensity
      score estimation; maxiter() calls -set maxiter-; the original value is restored
      after running the PS command; default is maxiter() = min(50,c(maxiter))
    - outcome variables no longer have to be unique
    - comlumn "bandwidth" no longer displayed in matching table in cases where 
      there is no bandwidth; title of column is "Caliper" in case of nn-matching;
      formating of matching table now takes linesize into account; revides alsp some 
      other aspect of results display
    - kmatch summarize/csummarize now also support skewness
    - kmatch csummarize used iweights for computation of variances/standard 
      deviations; this was appropriate for frequency weighted data but not in other
      cases; this is fixed
    - kmatch density and kmatch box could crash in case of negative matching weights
      (which are possible with ridge matching); the commands now treat negative
      weight as zero
    - exact matching returned error if there were no matches; this is fixed
    - -kmatch md, nn()- could crash under some data constellations if -bwidth()- 
      was specified; this is fixed
      [explanation: this was due to select() returning 0x0 if the input is 1x1 and 
      no elements are selected; subsequent code expected 0x1, as is returned in all
      other cases, i.e. if input is rx1 with r!=1]
    - kmatch could crash if the treatment variable had no variance (i.e. if one 
      of the groups was empty); this should now be fixed
    - kmatch now returns error if the variable names requested by generate() or dy()
      are not unique
    - kmatch md was not running in Stata 11 and 12 because it made use of Mata 
      function selectindex() that was not available prior to Stata 13; this is fixed

    22jun2017
    - kmatch returned error if no covariates and no ematch() variables were 
      specified; this is fixed
    - csummarize used the standard deviation of the matched sample (instead of the 
      standard deviation of the total sample) to compute the standardized 
      differences; this is fixed

    13jun2017
    - bw(cv over, weighted) had a bug so that wrong weights were used; this is fixed
    - matching algorithms now handle ties in the data more efficiently
    - matching weights with fweights now the same as in expanded data
    - kmatch csummmarize did not account for weights when computing statistics 
      for the unmatched; this is fixed

    09jun2017
    - kmatch csummarize had a wrong label in the rightmost column of the variance
      table; this is fixed
    - kmatch density, kmatch cdensity, and kmatch cbox did not always include labels
      for the variables; this is fixed

    08jun2017
    - penalty added for large bandwidths in bwidth(cv); suboption -nolimit- 
      deactivates the penalty

    07jun2017
    - results from -kmatch ps, nn(1)- were not always equal to results from
      -teffects psmatch, nn(1)-; this is fixed
    - results from -kmatch md, nn(#)- are not always equal to results from
      -teffects nnmatch, nn(#)-; this has to do with the fact that 
      -teffects nnmatch- treats controls as tied if their distance to the treatment
      observation does not differ by more than -dtolerance()-; setting 
      -dtolerance()- to a very small value, e.g. to -smallestdouble()-, should
      make results from -kmatch md- and -teffects nnmatch- equal
    - results from -kmatch md, nn(#)- with bias adjustment are not always equal to
      the results of -teffects nnmatch, nn(#)- with bias adjustment; this is because
      collinear variables in the treatment or control group are handled differently;
      results sould be equal if only non-collinear variables are included in the
      bias adjustment
  
    02jun2017
    - PM: now using 90% quantile of nonzero differences

    30may2017
    - there was a bug with how ties were handled in the cv-outcome algorithm so that
      results were wrong (and unstable); this is fixed  
    - changed some of the labeling/naming in output and returns
    - changed how CV results are returned
    - options noatt and noatc in cvplot are now called notreated and nountreated
    - -kmatch md- crashed if there were no covariate; this is fixed
    - bandwidth selection is now skipped if there are no covariates
    - notes about over category when computing BW displayed counter instead of over
      value; this is fixed

    23may2017
    - -cvplot, sort range()- did not connect all displayed points; this is fixed

    22may2017
    - fweights: results from CV are now the same as in the expanded data
    - CV with respect to outcome in -kmatch md- did not work with weights; this 
      is fixed

    20may2017
    - PM algorithm now takes account of weights when computing the minimum distance quantile

    19may2017
    - MD: added epsilon(h2) to h2 to compensate for possible roundoff error

    18may2017
    - option -sharedbw- added
    - -kmatch md- used the original X instead of the normalized X for 
       no-outcome-CV (unless mdmethod(1) was specified); this is fixed
