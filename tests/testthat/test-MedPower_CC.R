test_that("MedPower_CC returns expected sample size for given power", {
  result <- MedPower_CC(
    beta1 = sqrt(0.25),
    theta2 = 0.29,
    theta3 = 0,
    power = 0.8
  )
  expect_equal(result$required_n, 209)
  expect_equal(result$achieved_power, 0.8, tolerance = 1e-2)
})
