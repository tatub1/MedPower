#' Sample size or power calculation for mediation analysis (BC)
#'
#' Calculate power or sample size for the joint significance test when the
#' mediator is binary, the outcome is continuous, and exposure A is coded 0/1.
#'
#' @param alpha Significance level. Default is 0.05.
#' @param pA Proportion with A = 1. Default is 0.5.
#' @param pM Marginal probability Pr(M = 1) for the no-covariate reference model,
#'   used only when beta0 is NULL. Default is 0.5. It cannot identify the
#'   conditional mediator intercept in a covariate-adjusted model.
#' @param sigmaY Residual standard deviation in the adjusted linear outcome
#'   model Y ~ A + M + A:M + C. Default is 1.
#' @param theta0 Outcome-model intercept, retained for compatibility. It does
#'   not enter the homoskedastic outcome slope-variance formula.
#' @param theta1 Exposure coefficient in the outcome model, retained for
#'   compatibility. It does not enter the outcome slope-variance formula.
#' @param beta1 Exposure coefficient in the logistic mediator model (log-odds).
#' @param theta2 Mediator coefficient in the linear outcome model at A = 0.
#' @param theta3 Exposure-mediator interaction coefficient. Default is 0.
#' @param n Sample size, an integer of at least 2. Supply exactly one of n or power.
#' @param power Target power between 0 and 1.
#' @param vif_m User-specified additional variance inflation factor for beta1.
#'   The manuscript uses 1 / (1 - R_M^2), where R_M^2 is from the auxiliary
#'   linear regression A ~ C, not a logistic pseudo-R-squared. The correction
#'   is approximate for the logistic mediator model. Must be finite and >= 1.
#'   Default is 1, meaning no additional variance inflation is assumed.
#' @param vif_y User-specified additional variance inflation factor for
#'   delta1 = theta2 + theta3. Define M0 = (1 - A) * M and M1 = A * M.
#'   The manuscript's incremental factor for adding C is RSS_base / RSS_full
#'   from the auxiliary linear regressions M1 ~ A + M0 and M1 ~ A + M0 + C,
#'   fitted to the same observations. It is not the ordinary VIF of M1 relative
#'   to an intercept-only model. Must be finite and >= 1. Default is 1.
#' @param beta0 Optional intercept of the logistic mediator model at A = 0
#'   and C = 0. If NULL, it is recovered from pA, pM, and beta1 under a
#'   no-covariate reference model. Supply a calibrated beta0 for a covariate-
#'   adjusted design; marginal pM and beta1 do not identify that intercept.
#'   Do not explicitly supply both beta0 and pM. The default pM is ignored
#'   when beta0 is supplied.
#' @param pM_A1 Probability Pr(M = 1 | A = 1), averaging over the covariates
#'   in the exposed group, used for the outcome-model variance. Defaults to
#'   plogis(beta0 + beta1) in the no-covariate reference model. For a covariate-
#'   adjusted design, supply both beta0 and pM_A1 from the same assumed joint
#'   distribution or suitable pilot data. A supplied pM_A1 requires beta0.
#'   Do not substitute the probability at C = 0 for the group probability.
#'
#' @details
#' VIFs and adjusted-design inputs are user-specified. This function does not
#'   infer the joint distribution of A, M, and C or calibrate a conditional
#'   intercept from marginal prevalence in a covariate-adjusted model.
#'   The no-covariate reference calibration solves
#'   pM = (1 - pA) * plogis(beta0) + pA * plogis(beta0 + beta1).
#'   This equality is not generally valid with covariates in the mediator model.
#'
#' The mediator variance factor uses the reference probabilities
#'   r0 = plogis(beta0) and r1 = plogis(beta0 + beta1):
#'   vif_m * (1 / ((1 - pA) * r0 * (1 - r0)) +
#'            1 / (pA * r1 * (1 - r1))).
#'   With covariates, these are probabilities at C = 0, not group averages;
#'   the formula with vif_m is a working approximation for the conditional
#'   logistic coefficient variance. Group probabilities alone also do not
#'   determine its full logistic Fisher information. Use coherent inputs and
#'   check design accuracy by simulation when covariate effects are substantial.
#'   A VIF of 1 does not assert absence of confounding.
#'
#' The original paper's power approximation is retained:
#'   Phi(sqrt(n) * effect / sqrt(variance_factor) - qnorm(1 - alpha / 2)).
#'   It retains only the upper rejection tail and requires positive beta1 and
#'   positive theta2 + theta3. It is not the full two-sided normal power formula.
#'   The sample-size search is limited to n <= 100000; an unattainable target
#'   within that range produces an error.
#'
#' For the correctly specified homoskedastic linear outcome model, the variance
#'   factor is sigmaY^2 * vif_y / (pA * pM_A1 * (1 - pM_A1)). This design-based
#'   correction is exact with the true exposed-group probability, appropriate
#'   incremental vif_y, and residual variance. Using a probability at C = 0
#'   in place of the group probability does not preserve that exactness.
#'   The logistic mediator correction and overall power remain approximate.
#'
#' @returns A list with input_n and computed_power when n is supplied, or
#'   required_n and achieved_power when power is supplied. Both forms also
#'   report vif_m, vif_y, beta0, pM_A1, sigmaY, and reference probabilities
#'   pM_A0_reference, pM_A1_reference, pM_reference. The reference probabilities
#'   are not adjusted-design group or marginal probabilities in general.
#'   mediator_input and pM_A1_source identify how inputs were obtained.
#'
#' @export
#' @examples
#' # Original example.
#' MedPower_BC(theta0 = -4, theta1 = log(2), beta1 = 0.39,
#'             theta2 = log(3) * 4 / 5, theta3 = log(3) / 5, power = 0.8)
#'
#' # Illustrative adjusted-design inputs: supply beta0 instead of pM.
#' MedPower_BC(theta0 = 0, theta1 = 0.1, beta0 = -0.6, beta1 = 0.39,
#'             pM_A1 = 0.46, sigmaY = 1, theta2 = 0.2, theta3 = 0.05,
#'             vif_m = 1.2, vif_y = 1.25, power = 0.8)
MedPower_BC <- function(
    alpha = 0.05,
    pA = 0.5,
    pM = 0.5,
    sigmaY = 1,
    theta0,
    theta1,
    beta1,
    theta2,
    theta3 = 0,
    n = NULL,
    power = NULL,
    vif_m = 1,
    vif_y = 1,
    beta0 = NULL,
    pM_A1 = NULL
) {
  if (is.null(n) && is.null(power)) stop("Either 'n' or 'power' must be specified.")
  if (!is.null(n) && !is.null(power)) stop("Specify only one of 'n' or 'power'.")

  scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
    }
  }
  inputs <- list(alpha = alpha, pA = pA, theta0 = theta0, theta1 = theta1,
                 beta1 = beta1, theta2 = theta2, theta3 = theta3,
                 vif_m = vif_m, vif_y = vif_y)
  for (nm in names(inputs)) scalar(inputs[[nm]], nm)
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  if (pA <= 0 || pA >= 1) stop("'pA' must be between 0 and 1.")
  if (vif_m < 1 || vif_y < 1) stop("'vif_m' and 'vif_y' must be at least 1.")
  delta1 <- theta2 + theta3
  if (!is.finite(delta1)) stop("'theta2 + theta3' must be finite.")
  if (beta1 <= 0 || delta1 <= 0) {
    stop(paste("The retained upper-tail power approximation requires positive",
               "'beta1' and positive 'theta2 + theta3'."))
  }
  scalar(sigmaY, "sigmaY")
  if (sigmaY <= 0) stop("'sigmaY' must be positive.")
  beta0_supplied <- !is.null(beta0)
  pM_A1_supplied <- !is.null(pM_A1)
  if (beta0_supplied) {
    scalar(beta0, "beta0")
    if (!missing(pM)) {
      stop("Supply 'beta0' or 'pM', not both; 'pM' only calibrates the no-covariate reference model.")
    }
    mediator_input <- "beta0"
  } else {
    if (pM_A1_supplied) {
      stop("Supply a calibrated 'beta0' when supplying 'pM_A1'.")
    }
    scalar(pM, "pM")
    if (pM <= 0 || pM >= 1) stop("'pM' must be between 0 and 1.")
    # A monotone root avoids overflow/cancellation in the quadratic formula.
    upper <- stats::qlogis(pM)
    lower <- upper - beta1
    marginal_difference <- function(b0) {
      (1 - pA) * stats::plogis(b0) + pA * stats::plogis(b0 + beta1) - pM
    }
    if (lower == upper || marginal_difference(upper) == 0) {
      beta0 <- upper
    } else if (marginal_difference(lower) == 0) {
      beta0 <- lower
    } else {
      beta0 <- stats::uniroot(marginal_difference, c(lower, upper),
                              tol = 1e-12, maxiter = 1000L)$root
    }
    mediator_input <- "pM_reference"
  }

  eta_M1 <- beta0 + beta1
  if (!is.finite(eta_M1)) stop("'beta0 + beta1' must be finite.")
  pM_A0_reference <- stats::plogis(beta0)
  pM_A1_reference <- stats::plogis(eta_M1)
  if (!pM_A1_supplied) pM_A1 <- pM_A1_reference
  scalar(pM_A1, "pM_A1")
  if (pM_A1 <= 0 || pM_A1 >= 1) {
    stop("The supplied or derived 'pM_A1' must be numerically between 0 and 1.")
  }

  # Work on the log scale so probabilities near a boundary do not suffer
  # cancellation from 1 - plogis(eta) in the reference-model information.
  log_inv_bernoulli_var <- function(eta) {
    -stats::plogis(eta, log.p = TRUE) -
      stats::plogis(eta, lower.tail = FALSE, log.p = TRUE)
  }
  log_add <- function(x, y) {
    m <- max(x, y)
    m + log(exp(x - m) + exp(y - m))
  }
  log_E <- log_inv_bernoulli_var(beta0) - log1p(-pA)
  log_G <- log_inv_bernoulli_var(eta_M1) - log(pA)
  VarM <- exp(log(vif_m) + log_add(log_E, log_G))
  VarY <- exp(2 * log(sigmaY) + log(vif_y) - log(pA) -
                log(pM_A1) - log1p(-pM_A1))
  outcome_assumptions <- list(sigmaY = sigmaY)
  if (!is.finite(VarM) || !is.finite(VarY) || VarM <= 0 || VarY <= 0) {
    stop("Variance calculation is outside the supported numeric range.")
  }
  Test_beta <- beta1 / sqrt(VarM)
  Test_theta <- delta1 / sqrt(VarY)
  if (!is.finite(Test_beta) || !is.finite(Test_theta)) {
    stop("Effect calculation is outside the supported numeric range.")
  }
  compute_power <- function(n_val) {
    power1 <- stats::pnorm(Test_beta * sqrt(n_val) - stats::qnorm(1 - alpha / 2))
    power2 <- stats::pnorm(Test_theta * sqrt(n_val) - stats::qnorm(1 - alpha / 2))
    power1 * power2
  }
  assumptions <- list(
    vif_m = vif_m, vif_y = vif_y, beta0 = beta0, pM_A1 = pM_A1,
    pM_A0_reference = pM_A0_reference, pM_A1_reference = pM_A1_reference,
    pM_reference = (1 - pA) * pM_A0_reference + pA * pM_A1_reference,
    mediator_input = mediator_input,
    pM_A1_source = if (pM_A1_supplied) "user" else "reference_model"
  )
  assumptions <- c(assumptions, outcome_assumptions)

  if (!is.null(n)) {
    scalar(n, "n")
    if (n < 2 || n != floor(n)) stop("'n' must be an integer of at least 2.")
    return(c(list(input_n = n, computed_power = compute_power(n)), assumptions))
  }
  scalar(power, "power")
  if (power <= 0 || power >= 1) stop("'power' must be between 0 and 1.")
  min_n <- 2
  max_n <- 100000
  if (compute_power(max_n) < power) {
    stop("Target power is not reached within the search limit of n = 100000.")
  }
  while (min_n < max_n) {
    mid <- floor((min_n + max_n) / 2)
    if (compute_power(mid) < power) min_n <- mid + 1 else max_n <- mid
  }
  c(list(required_n = min_n, achieved_power = compute_power(min_n)), assumptions)
}
