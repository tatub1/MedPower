#' Sample size or power calculation for mediation analysis (CC)
#'
#' Calculate power or sample size for the joint significance test when the
#' mediator and outcome are continuous and exposure A is coded 0/1.
#'
#' @param alpha Significance level. Default is 0.05.
#' @param pA Proportion with A = 1. Default is 0.5.
#' @param sigmaM Residual standard deviation in the adjusted mediator model
#'   M ~ A + C. This is not the marginal standard deviation of M. Default is 1.
#' @param sigmaY Residual standard deviation in the adjusted outcome model
#'   Y ~ A + M + A:M + C. Default is 1.
#' @param beta1 Exposure coefficient in the mediator model.
#' @param theta2 Mediator coefficient in the outcome model.
#' @param theta3 Exposure-mediator interaction coefficient. Default is 0.
#' @param n Sample size, an integer of at least 2. Supply exactly one of n or power.
#' @param power Target power between 0 and 1.
#' @param vif_m User-specified additional variance inflation factor for beta1
#'   due to adjustment for C. Equal to 1 / (1 - R_M^2), where R_M^2 is the
#'   ordinary linear-regression R-squared from A ~ C. Must be finite and >= 1.
#'   Default is 1, meaning no additional variance inflation is assumed.
#' @param vif_y User-specified additional variance inflation factor for
#'   delta1 = theta2 + theta3. Define M0 = (1 - A) * M and M1 = A * M.
#'   Compute RSS_base / RSS_full from the auxiliary linear regressions
#'   M1 ~ A + M0 and M1 ~ A + M0 + C, fitted to the same observations.
#'   This is the incremental factor for adding C, not the ordinary VIF of M1
#'   relative to an intercept-only model. Must be finite and >= 1. Default is 1.
#' @param sigmaM_A1 Standard deviation of M among subjects with A = 1,
#'   sqrt(Var(M | A = 1)), before conditioning on C. Defaults to sigmaM for
#'   compatibility with the original no-covariate calculation. Supply this
#'   separately when covariates affect the mediator distribution.
#'
#' @details
#' The function requires user-specified VIFs; it does not infer them from effect
#' sizes. They may be estimated using suitable pilot data or simulated design
#' data based on an assumed joint distribution of A, M, and C. Input quantities
#' must refer to the same design population and planned adjustment covariates.
#'
#' With consistent residual-variance and design inputs, the VIF correction is
#' exact for the coefficient variance under correctly specified homoskedastic
#' linear models. The power calculation remains an asymptotic approximation.
#' A VIF of 1 does not assert absence of confounding.
#'
#' This revision retains the original paper's power approximation:
#' Phi(sqrt(n) * effect / sqrt(variance_factor) - qnorm(1 - alpha / 2)).
#' It retains only the upper rejection tail and requires positive beta1 and
#' positive theta2 + theta3. It is not the full two-sided normal power formula.
#' The sample-size search is limited to n <= 100000. Inputs outside these
#' supported conditions produce an error rather than a sample-size result.
#'
#' @returns A list with input_n and computed_power when n is supplied, or
#'   required_n and achieved_power when power is supplied. Both forms also
#'   report vif_m, vif_y, sigmaM, sigmaM_A1, and sigmaY used in the calculation.
#'
#' @export
#' @examples
#' # Original example: no additional variance inflation.
#' MedPower_CC(beta1 = 0.5, theta2 = 0.29, power = 0.8)
#'
#' # Illustrative design inputs, to be replaced by study-specific assumptions.
#' MedPower_CC(beta1 = 0.5, theta2 = 0.29, power = 0.8,
#'             sigmaM = 0.8, sigmaM_A1 = 1, sigmaY = 1,
#'             vif_m = 1.2, vif_y = 1.25)
MedPower_CC <- function(
    alpha = 0.05,
    pA = 0.5,
    sigmaM = 1,
    sigmaY = 1,
    beta1,
    theta2,
    theta3 = 0,
    n = NULL,
    power = NULL,
    vif_m = 1,
    vif_y = 1,
    sigmaM_A1 = sigmaM
) {
  if (is.null(n) && is.null(power)) stop("Either 'n' or 'power' must be specified.")
  if (!is.null(n) && !is.null(power)) stop("Specify only one of 'n' or 'power'.")

  scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop(sprintf("'%s' must be a finite numeric scalar.", name), call. = FALSE)
    }
  }
  inputs <- list(alpha = alpha, pA = pA, sigmaM = sigmaM, sigmaY = sigmaY,
                 beta1 = beta1, theta2 = theta2, theta3 = theta3,
                 vif_m = vif_m, vif_y = vif_y, sigmaM_A1 = sigmaM_A1)
  for (nm in names(inputs)) scalar(inputs[[nm]], nm)
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must be between 0 and 1.")
  if (pA <= 0 || pA >= 1) stop("'pA' must be between 0 and 1.")
  if (sigmaM <= 0 || sigmaY <= 0 || sigmaM_A1 <= 0) {
    stop("All standard deviations must be positive.")
  }
  if (vif_m < 1 || vif_y < 1) stop("'vif_m' and 'vif_y' must be at least 1.")

  VarM <- sigmaM^2 / (pA * (1 - pA)) * vif_m
  VarY <- sigmaY^2 / (pA * sigmaM_A1^2) * vif_y
  if (!is.finite(VarM) || !is.finite(VarY) || VarM <= 0 || VarY <= 0) {
    stop("Variance calculation is outside the supported numeric range.")
  }
  Test_beta <- beta1 / sqrt(VarM)
  Test_theta <- (theta2 + theta3) / sqrt(VarY)
  if (!is.finite(Test_beta) || !is.finite(Test_theta)) {
    stop("Effect calculation is outside the supported numeric range.")
  }
  if (Test_beta <= 0 || Test_theta <= 0) {
    stop(paste("The retained upper-tail power approximation requires positive",
               "'beta1' and positive 'theta2 + theta3'."))
  }

  compute_power <- function(n_val) {
    power1 <- stats::pnorm(Test_beta * sqrt(n_val) - stats::qnorm(1 - alpha / 2))
    power2 <- stats::pnorm(Test_theta * sqrt(n_val) - stats::qnorm(1 - alpha / 2))
    power1 * power2
  }
  assumptions <- list(vif_m = vif_m, vif_y = vif_y, sigmaM = sigmaM,
                      sigmaM_A1 = sigmaM_A1, sigmaY = sigmaY)

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
