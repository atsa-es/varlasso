#' Fit a varlasso model to multivariate time series data
#'
#' \code{fit} is the primary function for fitting models with varlasso.
#'
#' @param data The data, as a long format dataframe
#' @param shared_q Optional vector with
#'   integers indicating which process variance parameters are shared; defaults
#'   to constant process variances shared across species.
#' @param shared_r Optional vector with
#'   integers indicating which observation variance parameters are shared;
#'   defaults to shared observation variances across species.
#' @param fixed_r Optional scalar or vector of values to be passed in for a model
#' where observation error variances are not estimated. Defaults to NULL, and them being
#' estimated. Cannot be set to zero.
#' @param est_trend Whether or not to estimate a species-specific trend, defaults to FALSE
#' @param varss Boolean, whether to fit a VARSS model (defaults to true) or VAR model
#' @param b_diag A 2 element list specifying the mean and sd of a normal prior on the diagonal elements.
#' These elements are `mu` (default= 0.7) and `sd` (default = 1)
#' @param b_diag_min The minimum value of the B diagonal, defaults to 0. Setting this lower,
#' e.g. -1, allows for oscillations
#' @param b_offdiag A 2 element list specifying the mean and sd of a normal prior on the off-diagonal elements.
#' These elements are `mu` (default= 0) and `sd` (default = 1)
#' @param off_diag_prior Character string denoting which prior to use for
#' off-diagonal elements. Can be "normal" (independent priors estimated for each)
#' parameter, "normal_pp" (normal prior with partial pooling, where model is hierarchical with an
#' estimated standard deviation), "student-t" (hierarchical student-t model used with estimated scale
#' and df), "laplace" (double-exponential model used with estimated scale) or
#' "hs" (regularized horseshoe prior used, following parameterization in rstanarm)
#' @param sigma_pro A 2 element list specifying the mean and sd of a half Student-t prior on the process
#' variance standard deviation. Elements of the list are `mu` (default= 0) and `sd` (default = 1), and the
#' prior df is fixed at 3
#' @param sigma_obs A 2 element list specifying the mean and sd of a half Student-t prior on the process
#' variance standard deviation. Elements of the list are `mu` (default= 0) and `sd` (default = 1), and the
#' prior df is fixed at 3
#' @param nu Optional known parameter (df) for Student-t priors
#' @param sigma_scale A 2 element list specifying the `df` and `sd` parameters for a
#' half Student-t prior on the scale parameter. This is only used for the Laplace and Student-t priors;
#' default values are 3 for `df` and 2 for `sd`
#' @param hs A list containing the hyperparameters of the horseshoe prior. Default values
#' are `df` (df of the Student-t prior for local shrinkage, defaults to 1),
#' `df_global` (Student-t df of the global shrinkage, defaults to 1),
#' `df_slab` (df of the Student-t prior on the slab regularization, defaults to 4),
#' `scale_global` (Scale for Student-t prior on global shrinkage, defaults to 1),
#' `scale_slab` (Scale of the Student-t prior on shrikage parameter, defaults to 2).
#' @param unique_reg Boolean, whether to estimate unique regularization parameters for the Student-t or Laplace
#' models (defaults to false)
#' @param pars_to_monitor Optional string vector of which parameters to monitor. By default this is
#' NULL and only most important parameters are monitored (variances, B matrix, etc). This can be
#' the word "all", in which case all parameters are saved (this may be good for
#' loo::loo_moment_match() for example).
#' @param iter Number of MCMC iterations, defaults to 1000
#' @param warmup Warmup / burn in phase, defaults to 500
#' @param thin MCMC thin, defaults to 1
#' @param chains MCMC chains, defaults to 3
#' @param save_log_lik Whether to save log_lik, defaults to FALSE for speed
#' @param estimation Whether to do MCMC sampling ("sampling"; default), maximum posterior estimation
#'  ("optimizing") or Stan's variational algrotim ("vb")
#' @param est_hessian Whether to estimate hessian matrix if doing MAP estimation, defaults to FALSE
#' @param prior_only Whether to only generate samples from the prior distribution, defaults to FALSE
#' @param ... Extra arguments to pass to sampling
#'
#' @return an object of class 'stanfit'
#'
#' @import rstan
#' @import methods
#' @import Rcpp
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' iter <- 50
#' warmup <- 20
#' thin <- 1
#' chains <- 1
#'
#' n_ts <- 4
#' n_time <- 70
#' n_burn <- 50
#' xx <- matrix(0, n_time, n_ts)
#' set.seed(123)
#' xx[1, ] <- rnorm(n_ts)
#' B <- matrix(rnorm(n = n_ts * n_ts, mean = 0, sd = 0.1), n_ts, n_ts)
#' diag(B) <- runif(n_ts, min = 0.7, max = 0.9)
#'
#' for (i in 2:n_time) {
#'   xx[i, ] <- B %*% xx[i - 1, ] + rnorm(n_ts, mean = 0, sd = 0.04)
#' }
#' xx <- xx[-c(1:n_burn), ]
#'
#' yy <- xx + rnorm(n = nrow(xx) * ncol(xx), mean = 0, sd = 0.02)
#' df <- data.frame(
#'   "time" = rep(1:nrow(yy), ncol(yy)),
#'   "species" = sort(rep(1:ncol(yy), nrow(yy))),
#'   "y" = c(yy)
#' )
#'
#' # fit basic model with normal independent priors
#' m <- fit(data = df, chains = chains, iter = iter, warmup = warmup, thin = thin)
#' # fit same model with shrinkage using student-t priors
#' m_t <- fit(
#'   data = df, chains = chains, iter = iter, warmup = warmup, thin = thin,
#'   off_diag_prior = "student-t"
#' )
#' # fit same model using student-t priors and known df
#' m_t_known <- fit(
#'   data = df, chains = chains, iter = iter, warmup = warmup, thin = thin,
#'   off_diag_prior = "student-t", nu = 4
#' )
#' # fit same model using laplace priors
#' m_lp <- fit(
#'   data = df, chains = chains, iter = iter, warmup = warmup, thin = thin,
#'   off_diag_prior = "laplace"
#' )
#' # fit same model using horseshoe priors
#' m_hs <- fit(
#'   data = df, chains = chains, iter = iter, warmup = warmup, thin = thin,
#'   off_diag_prior = "hs"
#' )
#'
#' # include example for just sampling from priors
#' m <- fit(
#'   data = df, chains = chains, iter = iter,
#'   warmup = warmup, thin = thin, prior_only = TRUE
#' )
#' }
fit <- function(data,
                shared_q = NULL,
                shared_r = NULL,
                fixed_r = NULL,
                est_trend = FALSE,
                varss = TRUE,
                b_diag = list(mu = 0.7, sd =1),
                b_diag_min = 0,
                b_offdiag = list(mu = 0, sd = 1),
                off_diag_prior = c("normal", "student-t", "laplace", "hs","normal_pp"),
                sigma_pro = list(mu = 0, sd = 1),
                sigma_obs = list(mu = 0, sd = 1),
                nu = NULL,
                sigma_scale = list(df = 3, sd = 2),
                hs = list(df = 1,
                                df_global = 1,
                                df_slab = 4,
                                scale_global = 1,
                                scale_slab = 2),
                unique_reg = FALSE,
                pars_to_monitor = NULL,
                iter = 1000,
                warmup = floor(iter / 2),
                thin = 1,
                chains = 3,
                save_log_lik = FALSE,
                estimation = c("sampling","optimizing","vb"),
                est_hessian = FALSE,
                prior_only = FALSE,
                ...) {

  # check data
  if ("species" %in% names(data) == FALSE) {
    stop("Error: make sure one of the names of the dataframe is 'species'")
  }
  if ("time" %in% names(data) == FALSE) {
    stop("Error: make sure one of the names of the dataframe is 'time'")
  }
  if ("y" %in% names(data) == FALSE) {
    stop("Error: make sure one of the names of the dataframe is 'y'")
  }
  # filter out data that include NA observations
  data <- data[which(!is.na(data$y)), ]

  n_pos <- nrow(data)
  n_time <- length(seq(min(data$time, na.rm = T), max(data$time, na.rm = T)))
  n_spp <- length(unique(data$species))
  n_off <- n_spp * (n_spp - 1)

  n_q <- length(unique(shared_q))
  n_r <- length(unique(shared_r))
  if (n_q == 0) {
    id_q <- rep(1, n_spp)
    n_q <- 1
  } else {
    id_q <- shared_q
  }
  if (n_r == 0) {
    id_r <- rep(1, n_spp)
    n_r <- 1
  } else {
    id_r <- shared_r
  }

  # construct indices for off diagonal elements
  rc_off <- matrix(0, n_off, 2)
  indx <- 1
  for (i in 1:n_spp) {
    for (j in 1:n_spp) {
      if (i != j) {
        rc_off[indx, ] <- c(i, j)
        indx <- indx + 1
      }
    }
  }

  # make sure priors are in right form
  if (length(b_diag$mu) %in% c(1, n_spp) == FALSE) {
    stop("Error: mean on diagonal for B (b_mu_diag) needs to be a scalar or vector whose length is the number of species")
  } else {
    if (length(b_diag$mu) == 1) {
      b_diag$mu <- rep(b_diag$mu, n_spp)
    }
  }
  if (length(b_diag$sd) %in% c(1, n_spp) == FALSE) {
    stop("Error: sd on diagonal for B (b_sd_diag) needs to be a scalar or vector whose length is the number of species")
  } else {
    if (length(b_diag$sd) == 1) {
      b_diag$sd <- rep(b_diag$sd, n_spp)
    }
  }
  if (length(b_offdiag$mu) %in% c(1, n_off) == FALSE) {
    stop("Error: mean on off diagonal for B (b_mu) needs to be a scalar or vector whose length is n_species * (n_species - 1)")
  } else {
    if (length(b_offdiag$mu) == 1) {
      b_offdiag$mu <- rep(b_offdiag$mu, n_off)
    }
  }
  if (length(b_offdiag$sd) %in% c(1, n_off) == FALSE) {
    stop("Error: sd on off diagonal for B (b_sd) needs to be a scalar or vector whose length is n_species * (n_species - 1)")
  } else {
    if (length(b_offdiag$sd) == 1) {
      b_offdiag$sd <- rep(b_offdiag$sd, n_off)
    }
  }

  if (length(off_diag_prior) > 1) off_diag_prior <- off_diag_prior[1]
  if (off_diag_prior %in% c("normal", "student-t", "laplace", "hs","normal_pp") == FALSE) {
    stop("Error: off_diag_prior must be one of 'normal','student-t','laplace','hs'")
  }
  off_diag_priors <- match(off_diag_prior, c("normal", "student-t", "laplace", "hs","normal_pp"))
  off_diag_priors <- off_diag_priors - 1 # index at 0

  nu_known <- 0
  if (off_diag_priors == 1) {
    if (!is.null(nu)) {
      if (nu < 2) stop("Error: if nu is entered, it must be >= 2")
      nu_known <- nu
    }
  }

  if(!is.null(fixed_r)) {
    if(length(fixed_r)==1) {
      fixed_r = rep(fixed_r, n_spp)
    } else {
      if(length(fixed_r) == n_spp) {
        # do nothing, right length
      } else {
        stop("Error: if fixed_r is used, it must be length = 1 or number of species")
      }
    }

  } else {
    fixed_r = rep(0, n_spp)
  }
  data_list <- list(
    n_pos = n_pos,
    n_time = n_time,
    n_spp = n_spp,
    n_off = n_spp * (n_spp - 1),
    n_q = n_q,
    n_r = n_r,
    fixed_r = fixed_r,
    est_trend = as.numeric(est_trend),
    id_q = id_q,
    id_r = id_r,
    time_index = data$time,
    species_index = data$species,
    yy = data$y,
    rc_off = rc_off,
    b_mu = b_offdiag$mu,
    b_sd = b_offdiag$sd,
    b_mu_diag = b_diag$mu,
    b_sd_diag = b_diag$sd,
    b_diag_min = b_diag_min,
    nu_known = nu_known,
    sigma_proc_sd = sigma_pro$sd,
    sigma_proc_mu = sigma_pro$mu,
    sigma_obs_sd = sigma_obs$sd,
    sigma_obs_mu = sigma_obs$mu,
    off_diag_priors = off_diag_priors,
    sigma_scale_df = sigma_scale$df,
    sigma_scale_sd = sigma_scale$sd,
    hs_df = hs$df,
    hs_df_global=hs$df_global,
    hs_df_slab=hs$df_slab,
    hs_scale_global=hs$scale_global,
    hs_scale_slab=hs$scale_slab,
    priors_only = as.numeric(prior_only),
    process_only = 1 - as.numeric(varss),
    est_unique_reg = as.numeric(unique_reg)
  )

  pars <- c("sigma_proc", "Bmat", "x0")
  if (sum(fixed_r) != 0) pars <- c(pars,"sigma_obs")
  if (save_log_lik == TRUE) pars <- c(pars, "log_lik")
  if (est_trend == TRUE) pars <- c(pars, "U")
  if (off_diag_priors > 0) pars <- c(pars, "sigma_scale")
  if (off_diag_priors == 1) {
    pars <- c(pars, "lambda2")
    if (is.null(nu)) {
      # nu is estimated
      pars <- c(pars, "nu")
    }
  }
  if (off_diag_priors == 2) {
    pars <- c(pars, "lambda2")
  }
  if (off_diag_priors == 3) {
    pars <- c(pars, "hs_global", "hs_local", "hs_slab")
  }

  if (!is.null(pars_to_monitor)) {
    pars <- pars_to_monitor
    # if all then monitor everything
    if (pars_to_monitor == "all") {
      pars <- NULL
    }
  }

  if (estimation[1] == "sampling") {
    out <- rstan::sampling(
      object = stanmodels$varlasso,
      data = data_list,
      pars = pars,
      warmup = warmup,
      iter = iter,
      thin = thin,
      chains = chains, ...
    )
  }
  if (estimation[1] == "optimizing") {
    sampling_args <- list(
      object = stanmodels$varlasso,
      data = data_list,
      # verbose = verbose,
      hessian = est_hessian,
      ...
    )
    out <- do.call(optimizing, sampling_args)
  }
  if (estimation[1] == "vb") {
    sampling_args <- list(
      object = stanmodels$varlasso,
      data = data_list,
      iter = iter,
      pars = pars,
      ...
    )
    out <- do.call(vb, sampling_args)
  }

  return(list(fit = out, data = data_list, pars = pars))
}
