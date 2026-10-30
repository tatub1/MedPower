# MedPower

MedPower calculates approximate power and sample size for joint significant
tests of natural indirect effects in mediation models. It supports a binary
exposure coded 0/1, continuous or binary mediators and outcomes, and an
exposure-mediator interaction in the outcome model.

Version 0.1.0 adds user-specified variance inflation factors (VIFs) for covariate
adjustment. MedPower requires these factors as inputs; it does not estimate them
automatically from data or from the effect sizes. A reproducible example of how
to specify VIFs at the design stage is provided below.

| Function | Mediator | Outcome |
| --- | --- | --- |
| `MedPower_CC()` | Continuous | Continuous |
| `MedPower_CB()` | Continuous | Binary |
| `MedPower_BC()` | Binary | Continuous |
| `MedPower_BB()` | Binary | Binary |

## Installation

Install the development version from [GitHub](https://github.com/tatub1/MedPower):

```r
# Run once if remotes is not installed:
# install.packages("remotes")
remotes::install_github("tatub1/MedPower")
```

Load the package:

```r
library(MedPower)
```

For reproducible research, record the package version and GitHub commit or
release tag used for the analysis.

## Quick start: continuous mediator and continuous outcome

Specify exactly one of `power` (to calculate sample size) or `n` (to calculate
power). The following values illustrate the interface and should be replaced
with study-specific assumptions.

```r
# No additional variance inflation is assumed.
cc_reference <- MedPower_CC(
  beta1 = 0.5,
  theta2 = 0.29,
  power = 0.8,
  vif_m = 1,
  vif_y = 1
)
cc_reference$required_n
# [1] 209

# Hold the other inputs fixed and increase both variance factors.
cc_inflated <- MedPower_CC(
  beta1 = 0.5,
  theta2 = 0.29,
  power = 0.8,
  vif_m = 1.5,
  vif_y = 1.5
)
cc_inflated$required_n
# [1] 313
```

At a fixed sample size, the same inflation reduces power:

```r
cc_power <- MedPower_CC(
  beta1 = 0.5,
  theta2 = 0.29,
  n = 209,
  vif_m = 1.5,
  vif_y = 1.5
)
cc_power$computed_power
# [1] 0.568487
```

The sample-size result contains `required_n` and `achieved_power`. The power
result contains `input_n` and `computed_power`. Both also report the VIFs and
other design inputs used in the calculation.

## Specifying VIFs at the design stage

Let `A` be the exposure, `M` the mediator, and `C` the set of planned adjustment
covariates. The outcome model includes `A`, `M`, `A:M`, and `C`. With
`M0 = (1 - A) * M` and `M1 = A * M`, it can be reparameterized using `A`, `M0`,
`M1`, and `C`. The coefficient of `M1` is `theta2 + theta3`, the mediator effect
among exposed individuals.

- `vif_m` adjusts the variance of `beta1`, the exposure coefficient in the
  mediator model. Set it to `1 / (1 - R_squared)`, using the ordinary linear
  auxiliary regression of `A` on `C`.
- `vif_y` adjusts the variance of `theta2 + theta3`. Fit the auxiliary linear
  regressions `M1 ~ A + M0` and `M1 ~ A + M0 + C` to the same observations, then
  divide the residual sum of squares of the first model by that of the second.
  This is the incremental inflation for adding `C`; an ordinary VIF for `M`
  alone does not represent this quantity when the outcome includes interaction.

The auxiliary regressions are linear even when the mediator or outcome model
is logistic. Logistic pseudo-R-squared values should not be substituted.

VIFs can be obtained from suitable pilot or external individual-level data, or
from a simulated design sample based on an assumed joint distribution of
`A`, `M`, and `C`. An exposure proportion and effect coefficients alone do not
identify these VIFs. Published covariance information is useful only if it
supplies the joint moments needed for the planned design, including the
interaction terms.

### Reproducible CC design example

This example generates a synthetic design sample to estimate design quantities.
The 10,000 simulated records are for estimating these quantities; they are not
the required study sample size. The example assumes linear mediator and outcome
models, with mediator residual SD 0.8 and outcome residual SD 1.

```r
set.seed(20261030)
n_design <- 10000L

design <- data.frame(
  C1 = rnorm(n_design),
  C2 = rbinom(n_design, size = 1, prob = 0.5)
)
design$A <- rbinom(
  n_design, size = 1,
  prob = plogis(-0.2 + 0.5 * design$C1 + 0.4 * design$C2)
)
design$M <- 0.2 + 0.5 * design$A +
  0.4 * design$C1 + 0.3 * design$C2 +
  rnorm(n_design, sd = 0.8)
```

For pilot data, replace the generation code above with a data frame containing
`A`, `M`, and the actual planned covariates. Adapt the formulas to those covariates.
Use the same analysis rows for the auxiliary regressions:

```r
variables <- c("A", "M", "C1", "C2")
d <- design[complete.cases(design[, variables]), variables]
stopifnot(all(d$A %in% c(0, 1)), length(unique(d$A)) == 2L)

d$M0 <- (1 - d$A) * d$M
d$M1 <- d$A * d$M

fit_A <- lm(A ~ C1 + C2, data = d)
fit_base <- lm(M1 ~ A + M0, data = d)
fit_full <- lm(M1 ~ A + M0 + C1 + C2, data = d)

vif_m_design <- 1 / (1 - summary(fit_A)$r.squared)
vif_y_design <- deviance(fit_base) / deviance(fit_full)

c(vif_m = vif_m_design, vif_y = vif_y_design)
```

Use a sufficiently large, representative design sample, ensure the regression
designs have full rank and positive residual sums of squares, and align the
missing-data handling with the planned study. If only a small pilot is available,
allow for uncertainty in the design inputs through sensitivity analyses.

Pass the estimated VIFs and the other compatible design inputs explicitly:

```r
cc_design <- MedPower_CC(
  pA = mean(d$A),
  beta1 = 0.5,
  theta2 = 0.29,
  theta3 = 0,
  sigmaM = 0.8,
  sigmaM_A1 = sd(d$M[d$A == 1]),
  sigmaY = 1,
  vif_m = vif_m_design,
  vif_y = vif_y_design,
  power = 0.8
)
cc_design$required_n
```

Here `sigmaM` is the residual SD in `M ~ A + C`, while `sigmaM_A1` is the SD of
`M` among individuals with `A = 1`, averaging over covariates in that group.
They generally differ when covariates predict the mediator. The outcome model
assumed for this calculation has mediator coefficient 0.29, interaction
coefficient 0, and homoskedastic residual SD 1. No observed outcome data are
needed to obtain the auxiliary VIFs.

## Other mediator and outcome types

These examples illustrate the additional design inputs required with covariates.
Their numerical values are hypothetical; they are not estimates from the CC
design sample above.

```r
# Continuous mediator, binary outcome.
cb_example <- MedPower_CB(
  beta0 = 0.14, beta1 = 0.26,
  theta0 = -4, theta1 = log(2),
  theta2 = log(1.5) * 4 / 5, theta3 = log(1.5) / 5,
  sigmaM = 0.8, sigmaM_A1 = 1, meanM_A1 = 0.5,
  vif_m = 1.2, vif_y = 1.25,
  power = 0.8
)
cb_example$required_n
# [1] 2451

# Binary mediator, continuous outcome.
bc_example <- MedPower_BC(
  theta0 = 0, theta1 = 0.1,
  beta0 = -0.6, beta1 = 0.39, pM_A1 = 0.46,
  theta2 = 0.2, theta3 = 0.05, sigmaY = 1,
  vif_m = 1.2, vif_y = 1.25,
  power = 0.8
)
bc_example$required_n
# [1] 1521

# Binary mediator, binary outcome.
bb_example <- MedPower_BB(
  theta0 = -4, theta1 = log(2),
  beta0 = -0.6, beta1 = 0.39, pM_A1 = 0.46,
  theta2 = log(3) * 4 / 5, theta3 = log(3) / 5,
  vif_m = 1.2, vif_y = 1.25,
  power = 0.8
)
bb_example$required_n
# [1] 1532
```

### Interpretation of design inputs

| Input | Meaning |
| --- | --- |
| `pA` | Probability that the exposure equals 1. |
| `sigmaM` | Residual SD in the adjusted linear mediator model (CC and CB). |
| `sigmaY` | Residual SD in the adjusted linear outcome model (CC and BC). |
| `sigmaM_A1` | SD of the mediator among the exposed, before conditioning on covariates (CC and CB). |
| `meanM_A1` | Mean of the mediator among the exposed, before conditioning on covariates (CB). |
| `beta0` | Mediator-model intercept at `A = 0` and `C = 0` (CB, BC and BB). |
| `pM_A1` | Probability of mediator value 1 among the exposed, averaging over their covariate distribution (BC and BB). |
| `pM` | Marginal mediator probability used to calibrate `beta0` in a no-covariate reference model (BC and BB). |

For CC and CB, the default `sigmaM_A1 = sigmaM` retains the original
no-covariate calculation. Specify the group SD separately for an adjusted design.
For CB, `meanM_A1` defaults to `beta0 + beta1`, which is generally insufficient
when covariates shift the mediator mean among the exposed.

For BC and BB, an adjusted design should supply both a calibrated `beta0` and
`pM_A1` from the same joint distribution. Do not explicitly supply `pM` together
with `beta0`. A marginal `pM` does not identify a conditional logistic intercept
when covariates are present. In logistic outcome models, `theta0` is a conditional
log-odds intercept and is not generally the logit of marginal outcome prevalence.

## Assumptions and scope

- All inputs must refer to the same study population, variable coding, and
  planned adjustment covariates. Specify the additional mediator-distribution
  inputs as well as the VIFs where needed.
- With compatible residual variances and design inputs, VIF corrections to
  coefficient variances are exact under correctly specified homoskedastic linear
  models. The joint-significance power calculation still uses an asymptotic
  approximation.
- For logistic models, multiplication by these unweighted VIFs is a working
  approximation. The full adjusted logistic information generally depends on
  the joint covariate distribution and outcome or mediator probabilities. The
  scalar VIFs and supplied group summaries do not determine it exactly.
- CB retains a rare-outcome approximation and a normal-moment expression for
  the mediator distribution among exposed individuals. Normal mediator
  residuals alone do not guarantee that this group distribution is normal.
- BB uses the saturated no-covariate logistic information as its reference
  calculation; its VIF adjustment with covariates remains approximate. This
  does not remove the causal or rare-outcome assumptions needed for the
  manuscript's natural-indirect-effect interpretation.
- The retained power expression uses only the upper rejection tail with the
  critical value `qnorm(1 - alpha / 2)`. It is not the full two-sided normal power
  expression. The functions require `beta1 > 0` and `theta2 + theta3 > 0`.
- VIFs must be finite numeric scalars of at least 1. A VIF of 1 means that no
  additional variance inflation is assumed; it does not establish absence of
  confounding.
- The required-sample-size search is limited to 100,000 individuals. An
  unattainable target within that range produces an error.

When logistic approximations may be inaccurate, assess the design using
simulations of the planned adjusted analysis. The auxiliary design example
above estimates input quantities; it is not a validation of power accuracy.

## Function documentation

```r
help("MedPower_CC", package = "MedPower")
help("MedPower_CB", package = "MedPower")
help("MedPower_BC", package = "MedPower")
help("MedPower_BB", package = "MedPower")
```

## License

MIT license. See [LICENSE](LICENSE) and [LICENSE.md](LICENSE.md).
