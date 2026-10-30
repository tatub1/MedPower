# Regression tests for user-specified VIF corrections in all four models.
# These tests check software behavior, not the scientific accuracy of the
# logistic VIF approximation for arbitrary covariate distributions.

local({
  cases <- list(
    CC = list(
      fun = MedPower_CC,
      args = list(beta1 = 0.5, theta2 = 0.29)
    ),
    CB = list(
      fun = MedPower_CB,
      args = list(
        beta0 = 0.14, beta1 = 0.26, theta0 = -2.001,
        theta1 = log(2), theta2 = log(1.5) * 4 / 5,
        theta3 = log(1.5) / 5
      )
    ),
    BC = list(
      fun = MedPower_BC,
      args = list(
        theta0 = -4, theta1 = log(2), beta1 = 0.39,
        theta2 = log(3) * 4 / 5, theta3 = log(3) / 5
      )
    ),
    BB = list(
      fun = MedPower_BB,
      args = list(
        theta0 = -4, theta1 = log(2), beta1 = 0.39,
        theta2 = log(3) * 4 / 5, theta3 = log(3) / 5
      )
    )
  )

  for (model in names(cases)) {
    local({
      label <- model
      case <- cases[[label]]
      calculate <- function(...) {
        do.call(case$fun, c(case$args, list(...)))
      }

      test_that(paste(label, "defaults are equivalent to VIFs of one"), {
        expect_equal(
          calculate(power = 0.8),
          calculate(power = 0.8, vif_m = 1, vif_y = 1)
        )
        expect_equal(
          calculate(n = 500),
          calculate(n = 500, vif_m = 1, vif_y = 1)
        )
      })

      test_that(paste(label, "each VIF affects sample size and power"), {
        baseline <- calculate(power = 0.8, vif_m = 1, vif_y = 1)
        fixed_n <- baseline$required_n
        baseline_power <- calculate(n = fixed_n)$computed_power
        settings <- list(
          mediator_only = c(1.5, 1),
          outcome_only = c(1, 1.5),
          both = c(1.5, 1.5)
        )

        for (setting in settings) {
          increased_n <- calculate(
            power = 0.8, vif_m = setting[1], vif_y = setting[2]
          )
          decreased_power <- calculate(
            n = fixed_n, vif_m = setting[1], vif_y = setting[2]
          )
          # Integer sample sizes can stay unchanged after a small power decrease.
          expect_gte(increased_n$required_n, fixed_n)
          expect_lt(decreased_power$computed_power, baseline_power)
        }
        expect_gt(
          calculate(power = 0.8, vif_m = 1.5, vif_y = 1.5)$required_n,
          fixed_n
        )
      })

      test_that(paste(label, "common VIF inflation rescales effective sample size"), {
        # Doubling both variance factors and n preserves power.
        expect_equal(
          calculate(n = 1000, vif_m = 2, vif_y = 2)$computed_power,
          calculate(n = 500, vif_m = 1, vif_y = 1)$computed_power,
          tolerance = 1e-12
        )
        expect_equal(
          calculate(n = 1000, vif_m = 2.4, vif_y = 3)$computed_power,
          calculate(n = 500, vif_m = 1.2, vif_y = 1.5)$computed_power,
          tolerance = 1e-12
        )
      })

      test_that(paste(label, "VIF-adjusted sample size is sufficient and minimal"), {
        result <- calculate(power = 0.8, vif_m = 1.2, vif_y = 1.5)
        at_n <- calculate(n = result$required_n, vif_m = 1.2, vif_y = 1.5)
        before_n <- calculate(
          n = result$required_n - 1, vif_m = 1.2, vif_y = 1.5
        )

        expect_gte(at_n$computed_power, 0.8)
        expect_lt(before_n$computed_power, 0.8)
        expect_equal(result$achieved_power, at_n$computed_power)
        expect_equal(result$vif_m, 1.2)
        expect_equal(result$vif_y, 1.5)
      })

      test_that(paste(label, "invalid VIF inputs are rejected"), {
        invalid_values <- list(
          0.5, 0, -1, Inf, NA_real_, NaN, "1.5", c(1, 2), numeric(0)
        )
        for (bad_vif in invalid_values) {
          expect_error(calculate(n = 500, vif_m = bad_vif), "vif_m")
          expect_error(calculate(n = 500, vif_y = bad_vif), "vif_y")
        }
      })
    })
  }
})
