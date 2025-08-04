test_that("MedPower_BB returns expected sample size for given power", {
  result <- MedPower_BB(
    pA = 0.5,
    pM = 0.5,
    theta0 = -4,
    theta1 = log(2),
    beta1 = 0.39,
    theta2 = log(3) * 4/5,
    theta3 = log(3) * 1/5,
    power = 0.8
  )
  expect_equal(result$required_n, 1288)
  expect_equal(result$achieved_power, 0.8, tolerance = 1e-2)
})
