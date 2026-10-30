#' Sample size or power calculation for mediation analysis (CB)
#'
#' Calculate power or sample size for the joint significance test when the
#' mediator is continuous, the outcome is binary, and exposure A is coded 0/1.
#'
#' @param alpha Significance level. Default is 0.05.
#' @param pA Proportion with A = 1. Default is 0.5.
#' @param sigmaM Residual standard deviation in the adjusted linear mediator
#'   model M ~ A + C. This is not the marginal standard deviation of M.
#'   Default is 1.
#' @param beta0 Intercept in the adjusted mediator model, at A = 0 and C = 0.
#'   Used with beta1 to form the default meanM_A1.
#' @param beta1 Exposure coefficient in the mediator model.
#' @param theta0 Log-odds intercept in the adjusted outcome model, at
#'   A = 0, M = 0, and C = 0. This is not the logit of marginal prevalence.
#' @param theta1 Exposure coefficient in the logistic outcome model, at M = 0.
#' @param theta2 Mediator coefficient in the logistic outcome model, at A = 0.
#' @param theta3 Exposure-mediator interaction coefficient in the logistic
#'   outcome model. Default is 0.
#' @param n Sample size, an integer of at least 2. Supply exactly one of n or power.
#' @param power Target power between 0 and 1.
#' @param vif_m User-specified additional variance inflation factor for beta1
#'   due to adjustment for C. Equal to 1 / (1 - R_M^2), where R_M^2 is the
#'   ordinary linear-regression R-squared from A ~ C. Must be finite and >= 1.
#'   Default is 1, meaning no additional variance inflation is assumed.
#' @param vif_y User-specified additional variance inflation factor used as a
#'   working approximation for delta1 = theta2 + theta3 in the logistic model.
#'   Define M0 = (1 - A) * M and M1 = A * M. Estimate the manuscript's factor
#'   as RSS_base / RSS_full from the auxiliary linear regressions M1 ~ A + M0
#'   and M1 ~ A + M0 + C, fitted to the same observations. This is an incremental
#'   factor for adding C, not the ordinary VIF of M1 relative to an intercept-only
#'   model. Must be finite and >= 1. Default is 1.
#' @param sigmaM_A1 Standard deviation of M among subjects with A = 1,
#'   sqrt(Var(M | A = 1)), before conditioning on C. Defaults to sigmaM for
#'   compatibility with the original no-covariate calculation. Supply this
#'   separately when covariates affect the mediator distribution.
#' @param meanM_A1 Mean of M among subjects with A = 1, E(M | A = 1), before
#'   conditioning on C. Defaults to beta0 + beta1 for compatibility with the
#'   original no-covariate calculation. In a linear mediator model with
#'   covariates, it equals beta0 + beta1 + beta_C' E(C | A = 1); specify the
#'   group mean explicitly when that last contribution is nonzero. Centering
#'   covariates in the whole population does not ensure E(C | A = 1) = 0.
#'
#' @details
#' VIFs and mediator distribution summaries are user-specified inputs. They
#' may be estimated from suitable pilot data or simulated design data based
#' on an assumed joint distribution of A, M, and C. All inputs must refer to
#' the same design population, coding, and planned adjustment covariates.
#' This function does not estimate VIFs, infer the joint distribution, or
#' calibrate outcome-model intercepts from marginal outcome prevalence.
#'
#' For the correctly specified homoskedastic linear mediator model, the VIF
#' correction to the coefficient variance is exact with consistent design
#' and residual-variance inputs. For the logistic outcome model it is a
#' working approximation, not an exact logistic information calculation.
#' A VIF of 1 does not assert absence of confounding.
#'
#' Writing delta1 = theta2 + theta3, the outcome variance factor is
#' vif_y / (sigmaM_A1^2 * J), where
#' J = pA * exp(theta0 + theta1 + meanM_A1 * delta1 +
#'              sigmaM_A1^2 * delta1^2 / 2).
#' This retains the manuscript's Whittemore rare-outcome approximation and
#' normal-moment expression for the distribution of M given A = 1. The latter
#' moment expression is exact if that distribution is normal; normal mediator
#' residuals alone do not ensure it is normal after averaging over covariates.
#' These summaries and an unweighted VIF do not in general recover the full
#' logistic Fisher information with covariates. When outcome risks are not
#' small, the mediator distribution is nonnormal, or covariate effects are
#' substantial, assess accuracy using the planned analysis in design simulations.
#'
#' The original paper's power approximation is retained:
#' Phi(sqrt(n) * effect / sqrt(variance_factor) - qnorm(1 - alpha / 2)).
#' It retains only the upper rejection tail and requires positive beta1 and
#' positive delta1. It is not the full two-sided normal power formula.
#' The sample-size search is limited to n <= 100000. Inputs outside these
#' supported conditions produce an error rather than a sample-size result.
#'
#' @returns A list with input_n and computed_power when n is supplied, or
#'   required_n and achieved_power when power is supplied. Both forms also
#'   report vif_m, vif_y, sigmaM, sigmaM_A1, and meanM_A1 used in the calculation.
#'
#' @export
#' @examples
#' # Original example: no additional variance inflation.
#' MedPower_CB(beta0 = 0.14, beta1 = 0.26, theta0 = -2.001,
#'             theta1 = log(2), theta2 = log(1.5) * 4 / 5,
#'             theta3 = log(1.5) / 5, power = 0.8)
#'
#' # Illustrative design inputs, to be replaced by study-specific assumptions.
#' MedPower_CB(beta0 = 0.14, beta1 = 0.26, theta0 = -4,
#'             theta1 = log(2), theta2 = log(1.5) * 4 / 5,
#'             theta3 = log(1.5) / 5, sigmaM = 0.8,
#'             sigmaM_A1 = 1, meanM_A1 = 0.5,
#'             vif_m = 1.2, vif_y = 1.25, power = 0.8)
MedPower_CB <- function(
    alpha = 0.05,
    pA = 0.5,
    sigmaM = 1,
    beta0,
    beta1,
    theta0,
    theta1,
    theta2,
    theta3 = 0,
    n = NULL,
    power = NULL,
    vif_m = 1,
    vif_y = 1,
    sigmaM_A1 = sigmaM,
    meanM_A1 = beta0 + beta1
) {
  if (is.null(n) && is.null(power)) stop("Either 'n' or 'power' must be specified.")
  if (!is.null(n) && !is.null(power)) stop("Specify only one of 'n' or 'power'.")

  scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
    }
  }
  inputs <- list(alpha = alpha, pA = pA, sigmaM = sigmaM,
                 beta0 = beta0, beta1 = beta1, theta0 = theta0,
                 theta1 = theta1, theta2 = theta2, theta3 = theta3,
                 vif_m = vif_m, vif_y = vif_y,
                 sigmaM_A1 = sigmaM_A1, meanM_A1 = meanM_A1)
  for (nm in names(inputs)) scalar(inputs[[nm]], nm)
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  if (pA <= 0 || pA >= 1) stop("'pA' must be between 0 and 1.")
  if (sigmaM <= 0 || sigmaM_A1 <= 0) {
    stop("All standard deviations must be positive.")
  }
  if (vif_m < 1 || vif_y < 1) stop("'vif_m' and 'vif_y' must be at least 1.")

  delta1 <- theta2 + theta3
  if (!is.finite(delta1)) stop("'theta2 + theta3' must be finite.")
  if (beta1 <= 0 || delta1 <= 0) {
    stop(paste("The retained upper-tail power approximation requires positive",
               "'beta1' and positive 'theta2 + theta3'."))
  }

  VarM <- sigmaM^2 / (pA * (1 - pA)) * vif_m

  # Log arithmetic avoids separate overflow/underflow of the exponential terms.
  log_J <- log(pA) + theta0 + theta1 + meanM_A1 * delta1 +
    0.5 * (sigmaM_A1 * delta1)^2
  log_VarY <- log(vif_y) - 2 * log(sigmaM_A1) - log_J
  VarY <- exp(log_VarY)
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
  assumptions <- list(vif_m = vif_m, vif_y = vif_y, sigmaM = sigmaM,
                      sigmaM_A1 = sigmaM_A1, meanM_A1 = meanM_A1)

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
    if (compute_power(mid) < power) {
      min_n <- mid + 1
    } else {
      max_n <- mid
    }
  }
  c(list(required_n = min_n, achieved_power = compute_power(min_n)), assumptions)
}
