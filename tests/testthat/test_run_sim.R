library(testthat)
library(BMsimulation)




test_that("incorrect input m1 (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 1.1, m2 = 2, n.pred = 100, ran = 5, nm = 1, nhours = 2,
                     ton = 100, burnin = 50, degree = 1, nbeta = 3,
                     ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                     num_for_cores = 1))
})

test_that("incorrect input m2 (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 20, m2 = 2.2, n.pred = 100, ran = 5, nm = 1, nhours = 2,
                        ton = 100, burnin = 50, degree = 1, nbeta = 3,
                        ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                        num_for_cores = 1))
})

test_that("incorrect input n.pred (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 20, m2 = 2, n.pred = 1.5, ran = 5, nm = 1, nhours = 2,
                        ton = 100, burnin = 50, degree = 1, nbeta = 3,
                        ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                        num_for_cores = 1))})

  test_that("incorrect input nm (positive integer)", {
    # Call run_func with some test values
    expect_error(run_func(m1 = 20, m2 = 2, n.pred = 100, ran = 5, nm = 1.1, nhours = 2,
                          ton = 100, burnin = 50, degree = 1, nbeta = 3,
                          ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                          num_for_cores = 1))})

test_that("incorrect input nhours (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 20, m2 = 2, n.pred = 100, ran = 5, nm = 1, nhours = 2.2,
                     ton = 100, burnin = 50, degree = 1, nbeta = 3,
                     ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                     num_for_cores = 1))})

test_that("incorrect input ton (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 20, m2 = 2, n.pred = 100, ran = 5, nm = 1, nhours = 2,
                     ton = 100.1, burnin = 50, degree = 1, nbeta = 3,
                     ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                     num_for_cores = 1))})

test_that("incorrect input burnin (positive integer)", {
  # Call run_func with some test values
  expect_error(run_func(m1 = 20, m2 = 2, n.pred = 100, ran = 5, nm = 1, nhours = 2,
                     ton = 100, burnin = 50.1, degree = 1, nbeta = 3,
                     ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                     num_for_cores = 1))})


test_that("run_func returns the correct output format", {
  # Call run_func with some test values
  result <- run_func(m1 = 20, m2 = 2, n.pred = 100, ran = 5, nm = 1, nhours = 2,
                     ton = 100, burnin = 50, degree = 1, nbeta = 3,
                     ntheta = 2, nab = 2, nerror = 2, cov.model = "exponential", seed = 42,
                     num_for_cores = 1)

  # Check if the output is a list with two elements
  expect_is(result, "list")
  expect_length(result, 2)

  # Check the first element of the output list (allresult)
  expect_is(result[[1]], "list")
  expect_length(result[[1]], 2)})  # Assuming nhours=2
  # Add more tests for other components of 'simulated_data' if necessary.

