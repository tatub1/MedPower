test_that("MedPower_BC returns expected sample size for given power", {
  result <- MedPower_BC(
    pA = 0.5,
    pM = 0.5,
    sigmaY = 1,
    theta0 = 0.14,
    theta1 = 0.14,
    beta1 = log(3),
    theta2 = 0.14 * 4/5,
    theta3 = 0.14 * 0/5,
    power = 0.8
  )
  expect_equal(result$required_n, 5393)
  expect_equal(result$achieved_power, 0.8, tolerance = 1e-2)
})
