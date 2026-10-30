# MedPower 0.1.0

## VIF corrections and design inputs

* Add user-specified `vif_m` and `vif_y` arguments to `MedPower_CC()`,
  `MedPower_CB()`, `MedPower_BC()`, and `MedPower_BB()`. The factors multiply
  the corresponding coefficient variance factors and default to 1.
* Define `vif_m` from the auxiliary linear regression of exposure on adjustment
  covariates. Define `vif_y` as the incremental residual-sum-of-squares ratio for
  adding covariates to the auxiliary regression for the exposed-group mediator
  term, accounting for exposure-mediator interaction.
* Add `sigmaM_A1` to CC and CB to distinguish the mediator SD among exposed
  individuals from the adjusted mediator residual SD. Add `meanM_A1` to CB
  to specify the corresponding exposed-group mediator mean.
* Add optional `beta0` and `pM_A1` inputs to BC and BB for specifying a conditional
  mediator intercept and an exposed-group mediator probability in an adjusted
  design. Retain `pM` calibration for a no-covariate reference model.
* Return the VIFs and relevant design assumptions with the calculation results.
  Preserve the existing leading positional arguments and primary result fields.

## Input validation and numerical behavior

* Validate probabilities, standard deviations, sample sizes, and VIF inputs.
  VIFs must be finite numeric scalars of at least 1.
* Retain the original upper-tail power approximation and explicitly require
  `beta1 > 0` and `theta2 + theta3 > 0`. Unsupported effect directions now produce
  an error; the calculation is not the full two-sided normal power expression.
* Report an error when the target power is not reached within the sample-size
  search limit of 100,000, rather than returning the search limit as a solution.
* Use numerical root finding for the reference mediator intercept in BC and BB,
  and log-scale calculations to improve numerical stability in CB, BC, and BB.
* Reject simultaneous explicit `pM` and `beta0` inputs in BC and BB. An explicit
  `pM_A1` requires an explicit `beta0`.

## Documentation and tests

* Document how researchers can specify VIFs from pilot or external data, or
  from an assumed joint distribution at the design stage. VIFs remain
  user-specified; the package does not calculate them automatically.
* Clarify that coefficient variance corrections are exact for correctly
  specified homoskedastic linear models with compatible design and residual
  variance inputs, and approximate for logistic models. Power calculations
  remain asymptotic approximations.
* Describe the model-specific assumptions, including the rare-outcome and
  normal-moment approximations used by CB.
* Add a README with examples for all four model types and a reproducible
  auxiliary-regression example for specifying VIFs.
* Add regression tests for default VIFs, separate and joint VIF effects,
  effective sample-size scaling, sample-size minimality, and invalid VIF inputs.
* Update the package description and declare the use of the standard R package
  `stats`.
