if (interactive()) options(mc.cores = parallel::detectCores())

tol_level <- 1e-6

iter <- 1000
warmup <- 500
thin <- 1
chains <- 1

n_ts <- 4
n_time <- 70
n_burn <- 50
xx <- matrix(0, n_time, n_ts)
set.seed(123)
xx[1, ] <- rnorm(n_ts)
B <- matrix(rnorm(n = n_ts * n_ts, mean = 0, sd = 0.1), n_ts, n_ts)
diag(B) <- runif(n_ts, min = 0.7, max = 0.9)

for (i in 2:n_time) {
  xx[i, ] <- B %*% xx[i - 1, ] + rnorm(n_ts, mean = 0, sd = 0.04)
}
xx <- xx[-c(1:n_burn), ]

yy <- xx + rnorm(n = nrow(xx) * ncol(xx), mean = 0, sd = 0.02)
df <- data.frame(
  "time" = rep(1:nrow(yy), ncol(yy)),
  "species" = sort(rep(1:ncol(yy), nrow(yy))),
  "y" = c(yy)
)


# ------------------------------------------------------
test_that("model with normal priors works ", {
  set.seed(123)
  m <- fit(
    data = df, chains = chains, iter = iter,
    warmup = warmup, thin = thin
  )
  expect_equal(class(m$fit)[1], "stanfit")

  pars <- rstan::extract(m$fit)

  expect_equal(mean(pars$sigma_obs), 0.0324634, tolerance = tol_level)
  expect_equal(mean(pars$sigma_proc), 0.01769837, tolerance = tol_level)

  b <- c(0.698, -0.376, 0.314, -0.561, -0.319)
  expect_equal(round(c(apply(pars$Bmat, c(2, 3), mean)), 3)[1:5], b,
    tolerance = tol_level
  )
})

# ------------------------------------------------------
test_that("model with student-t priors works ", {
  set.seed(123)
  m <- fit(
    data = df, chains = chains, iter = iter,
    warmup = warmup, thin = thin,
    off_diag_prior = "student-t"
  )
  expect_equal(class(m$fit)[1], "stanfit")

  pars <- rstan::extract(m$fit)
  expect_equal(mean(pars$nu), 20.85633, tolerance = tol_level)

  # try model with nu known
  set.seed(123)
  m <- fit(
    data = df, chains = chains, iter = iter,
    warmup = warmup, thin = thin, nu = 4,
    off_diag_prior = "student-t"
  )
  expect_equal(class(m$fit)[1], "stanfit")
})

# ------------------------------------------------------
test_that("model with laplace priors works ", {
  set.seed(123)
  m <- fit(
    data = df, chains = chains, iter = iter,
    warmup = warmup, thin = thin,
    off_diag_prior = "laplace"
  )
  expect_equal(class(m$fit)[1], "stanfit")

  pars <- rstan::extract(m$fit)
  #
  expect_equal(mean(pars$sigma_obs), 0.0204, tolerance = 0.01)
  expect_equal(mean(pars$sigma_proc), 0.0362, tolerance = 0.01)

  b <- c(0.809, 0.034, 0.047, -0.128, -0.066)
  expect_equal(round(c(apply(pars$Bmat, c(2, 3), mean)), 3)[1:5], b,
  tolerance = 0.01
  )
})

# ------------------------------------------------------
test_that("model with horseshoe priors works ", {
  set.seed(123)
  m <- fit(
    data = df, chains = chains, iter = iter,
    warmup = warmup, thin = thin,
    off_diag_prior = "hs"
  )
  expect_equal(class(m$fit)[1], "stanfit")

  # pars <- rstan::extract(m$fit)
  #
  # # expect_equal(mean(pars$sigma_obs), 0.017672, tolerance = tol_level)
  # # expect_equal(mean(pars$sigma_proc), 0.03836926, tolerance = tol_level)
  #
  # # b <- c(0.796,0.020,0.009,-0.038,-0.033)

  # expect_equal(round(c(apply(pars$Bmat, c(2, 3), mean)), 3)[1:5], b,
  #  tolerance = tol_level
  # )
})
