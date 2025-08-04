#' Sample size or power calculation for mediation analysis (variable type CB)
#'
#' This function computes either the statistical power or the required sample size
#' for detecting a natural indirect effect in a mediation model where the mediator
#' is continuous and the outcome is binary (variable type CB).
#'
#' @param alpha Significance level. Default is 0.05.
#' @param pA Proportion of treated or exposed individuals. Default is 0.5.
#' @param sigmaM Standard deviation of the continuous mediator. Default is 1.
#' @param beta0 Intercept in the mediator model.
#' @param beta1 Effect of exposure on mediator.
#' @param theta0 Intercept in the outcome model.
#' @param theta1 Effect of exposure on outcome.
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
#' MedPower_CB(
#' beta0 = 0.14,
#' beta1 = 0.26,
#' theta0 = -2.001,
#' theta1 = log(2),
#' theta2 = log(1.5) * 4 / 5,
#' theta3 = log(1.5) * 1 / 5,
#' power = 0.8
#' )
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
    power = NULL
) {
  if (is.null(n) && is.null(power)) stop("Either 'n' or 'power' must be specified.")
  if (!is.null(n) && !is.null(power)) stop("Specify only one of 'n' or 'power'.")

  VarM <- sigmaM^2 / (pA * (1 - pA))

  mom1 <- exp((beta0 + beta1) * (theta2 + theta3))
  mom2 <- exp((sigmaM^2 * (theta2 + theta3)^2) / 2)
  mom <- mom1 * mom2
  J <- pA * exp(theta0 + theta1) * mom
  VarY <- 1 / (sigmaM^2 * J)

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
