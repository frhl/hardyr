
# Test for get_allele_counts function
test_that("get_allele_counts returns correct values", {
  expect_equal(get_allele_counts(N = 100, K = 40, M = 10), c(10, 20, 70))
  expect_equal(get_allele_counts(N = 50, K = 10, M = 5), c(5, 0, 45))
})

test_that("get_allele_counts handles errors correctly", {
  expect_error(get_allele_counts(N = 100, K = 250, M = 10))
  expect_error(get_allele_counts(N = 100, K = 40, M = 30))
})

# Test for get_mac function
test_that("get_mac returns correct minimum allele count", {
  expect_equal(get_mac(N = 100, K = 40, M = 10), min(40, 160)) # 40 is the count for minor allele (hom_alt+hets), 160 is for major allele (hom_ref+hets)
  expect_equal(get_mac(N = 50, K = 20, M = 10), min(20, 80))
})

test_that("get_mac handles errors correctly", {
  expect_error(get_mac(N = 100, K = 250, M = 10))
  expect_error(get_mac(N = 100, K = 40, M = 30))
})

