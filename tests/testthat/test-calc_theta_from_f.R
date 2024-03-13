
# Testing across various pA's and F's
test_that("calc_theta_from_f and calc_f_from_theta are consistent", {
  pAs <- seq(0.01, 0.99, by = 0.1) # Range of pA values
  fs <- seq(0.01, 0.99, by = 0.1)  # Range of f values
  
  for (pA in pAs) {
    for (f in fs) {
      # Calculate theta from f
      theta <- calc_theta_from_f(pA, f)
      
      # Calculate f from the obtained theta
      # Note: calc_f_from_theta may return multiple roots, so we need to find the one closest to the original f
      fs_calculated <- calc_f_from_theta(pA, theta)
      
      # Check if the original f is close to any of the calculated f values
      expect_true(any(abs(fs_calculated - f) < 1e-5))
    }
  }
})


# Testing across various pA's and Thetas
test_that("calc_f_from_theta and calc_theta_from_f are consistent", {
  pAs <- seq(0.01, 0.99, by = 0.1) # Range of pA values
  thetas <- seq(1, 10, by = 0.5)   # Range of theta values for which the calculation is valid

  for (pA in pAs) {
    for (theta in thetas) {
      # Calculate f from theta
      fs <- calc_f_from_theta(pA, theta)

      # For each f calculated, check if calculating theta gives back the original theta
      # This loop is necessary because calc_f_from_theta can return multiple values
      for (f in fs) {
        theta_calculated <- calc_theta_from_f(pA, f)

        # Check if the calculated theta is close to the original theta
        expect_true(abs(theta_calculated - theta) < 1e-5)
      }
    }
  }
})
