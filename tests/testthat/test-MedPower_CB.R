test_that("MedPower_CB returns expected sample size for given power", {
  result <- MedPower_CB(
    beta0 = 0.14,
    beta1 = 0.26,
    theta0 = -2.001,
    theta1 = log(2),
    theta2 = log(1.5) * 4 / 5,
    theta3 = log(1.5) * 1 / 5,
    power = 0.8
  )
  expect_equal(result$required_n, 502)
  expect_equal(result$achieved_power, 0.8, tolerance = 1e-2)
})
