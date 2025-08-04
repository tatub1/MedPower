#' Sample size or power calculation for mediation analysis (variable type BC)
#'
#' This function computes either the statistical power or the required sample size
#' for detecting a natural indirect effect in a mediation model where the mediator
#' is binary and the outcome is continuous (variable type BC).
#'
#' @param alpha Significance level. Default is 0.05.
#' @param pA Proportion of treated or exposed individuals. Default is 0.5.
#' @param pM Marginal probability of the binary mediator being 1. Default is 0.5.
#' @param sigmaY Standard deviation of the continuous outcome. Default is 1.
#' @param theta0 Intercept in the outcome model.
#' @param theta1 Effect of exposure on outcome.
#' @param beta1 Effect of exposure on binary mediator (log-odds scale).
#' @param theta2 Effect of mediator on outcome.
#' @param theta3 Interaction effect between exposure and mediator on outcome.
#' @param n Sample size. If specified, the function returns the corresponding power.
#' @param power Target power. If specified, the function returns the required sample size.
#'
#' @returns A list containing either:
#' \describe{
#'   \item{required_n}{Sample size required to achieve the target power (if \code{power} is specified).}
#'   \item{achieved_power}{Power achieved with the given sample size (if \code{n} is specified).}
#' }
#'
#' @export
#'
#' @examples
#' MedPower_BC(
#'   theta0 = -4,
#'   theta1 = log(2),
#'   beta1 = 0.39,
#'   theta2 = log(3) * 4 / 5,
#'   theta3 = log(3) * 1 / 5,
#'   power = 0.8
#' )
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
    power = NULL
) {
  if (is.null(n) && is.null(power)) stop("Either 'n' or 'power' must be specified.")
  if (!is.null(n) && !is.null(power)) stop("Specify only one of 'n' or 'power'.")

  q <- pM * (1 + exp(beta1)) + pA * (1 - exp(beta1)) - 1
  beta0 <- log((q + sqrt(q^2 + 4 * pM * (1 - pM) * exp(beta1))) / (2 * (1 - pM) * exp(beta1)))
  exp_beta0 <- (q + sqrt(q^2 + 4 * pM * (1 - pM) * exp(beta1))) / (2 * (1 - pM) * exp(beta1))

  E <- (1 + exp(beta0))^2 / (exp(beta0) * (1 - pA))
  G <- (1 + exp(beta0 + beta1))^2 / (exp(beta0 + beta1) * pA)
  VarM <- E + G

  VarY <- sigmaY^2 * (1 + exp(beta0 + beta1))^2 / (exp(beta0 + beta1) * pA)

  Test_beta <- beta1 / sqrt(VarM)
  Test_theta <- (theta2 + theta3) / sqrt(VarY)

  compute_power <- function(n_val) {
    power1 <- pnorm(Test_beta * sqrt(n_val) - qnorm(1 - alpha / 2))
    power2 <- pnorm(Test_theta * sqrt(n_val) - qnorm(1 - alpha / 2))
    return(power1 * power2)
  }

  if (!is.null(n)) {
    pow <- compute_power(n)
    return(list(
      input_n = n,
      computed_power = pow
    ))
  }

  find_min_n <- function(target_power, min_n = 2, max_n = 100000) {
    while (min_n < max_n) {
      mid <- floor((min_n + max_n) / 2)
      if (compute_power(mid) < target_power) {
        min_n <- mid + 1
      } else {
        max_n <- mid
      }
    }
    return(min_n)
  }

  n_found <- find_min_n(power)
  power_at_n <- compute_power(n_found)

  return(list(
    required_n = n_found,
    achieved_power = power_at_n
  ))
}
