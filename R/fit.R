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
#'   defaults to shared observation variances across species
#' @param b_mu_diag Prior mean for diagonal elements. Can be vector or scalar
#' @param b_sd_diag Prior sd for diagonal elements. Can be vector or scalar
#' @param b_mu Prior mean for off diagonal elements, when Normal priors used.
#' Can be vector or scalar, defaults to 0
#' @param b_sd Prior sd for off diagonal elements, when Normal priors used.
#' Can be vector or scalar, defaults to 1
#' @param off_diag_prior Character string denoting which prior to use for
#' off-diagonal elements. Can be "normal" (independent priors estimated for each)
#' parameter, "student-t" (hierarchical student-t model used with estimated scale
#' and df), "laplace" (double-exponential model used with estimated scale) or
#' "hs" (regularized horseshoe prior used, following parameterization in rstanarm)
#' @param prior_sigma_pro_mu Optional, prior mean for for process variance scale. The
#' prior is a Student-t, with default mu = 0
#' @param prior_sigma_pro_sd Optional, prior standard deviation for process variance scale. The
#' prior is a Student-t, with default sd = 2
#' @param prior_sigma_obs_mu Optional, prior mean for for observation variance scale. The
#' prior is a Student-t, with default mu = 0
#' @param prior_sigma_obs_sd Optional, prior standard deviation for observation variance scale. The
#' prior is a Student-t, with default sd = 2
#' @param nu Optional known parameter (df) for Student-t priors
#' @param sigma_scale_df Degrees of freedom for Student-t prior on scaling variance for Student-t and Laplace models, defaults to 3
#' @param sigma_scale_sd Standard deviation for Student-t prior on scaling variance for Student-t and Laplace models, defaults to 2
#' @param global_scale scale parameter for regularized horseshoe prior, defaults to 0.1
#' @param slab_scale scale parameter for regularized horseshoe prior, defaults to 2.5
#' @param slab_df df parameter for regularized horseshoe prior, defaults to 4
#' @param pars_to_monitor Optional string vector of which parameters to monitor. By default this is
#' NULL and only most important parameters are monitored (variances, B matrix, etc). This can be
#' the word "all", in which case all parameters are saved (this may be good for
#' loo::loo_moment_match() for example).
#' @param iter Number of MCMC iterations, defaults to 1000
#' @param warmup Warmup / burn in phase, defaults to 500
#' @param thin MCMC thin, defaults to 1
#' @param chains MCMC chains, defaults to 3
#' @param save_log_lik Whether to save log_lik, defaults to FALSE for speed
#' @param map_estimation Whether to do MAP / penalized maximum likelihood rather than full Bayesian, defaults to FALSE
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
#' m <- fit(data = df, chains = chains, iter = iter,
#' warmup = warmup, thin = thin,prior_only = TRUE)
#' }
fit <- function(data,
                shared_q = NULL,
                shared_r = NULL,
                b_mu_diag = 0.7,
                b_sd_diag = 1,
                b_mu = 0,
                b_sd = 1,
                off_diag_prior = c("normal", "student-t", "laplace", "hs"),
                prior_sigma_pro_mu = 0,
                prior_sigma_pro_sd = 2,
                prior_sigma_obs_mu = 0,
                prior_sigma_obs_sd = 2,
                nu = NULL,
                sigma_scale_df = 3,
                sigma_scale_sd = 2,
                global_scale = 0.1,
                slab_scale = 2.5,
                slab_df = 4,
                pars_to_monitor = NULL,
                iter = 1000,
                warmup = floor(iter / 2),
                thin = 1,
                chains = 3,
                save_log_lik=FALSE,
                map_estimation=FALSE,
                est_hessian=FALSE,
                prior_only=FALSE,
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
  if (length(b_mu_diag) %in% c(1, n_spp) == FALSE) {
    stop("Error: mean on diagonal for B (b_mu_diag) needs to be a scalar or vector whose length is the number of species")
  } else {
    if (length(b_mu_diag) == 1) {
      b_mu_diag <- rep(b_mu_diag, n_spp)
    }
  }
  if (length(b_sd_diag) %in% c(1, n_spp) == FALSE) {
    stop("Error: sd on diagonal for B (b_sd_diag) needs to be a scalar or vector whose length is the number of species")
  } else {
    if (length(b_sd_diag) == 1) {
      b_sd_diag <- rep(b_sd_diag, n_spp)
    }
  }
  if (length(b_mu) %in% c(1, n_off) == FALSE) {
    stop("Error: mean on off diagonal for B (b_mu) needs to be a scalar or vector whose length is n_species * (n_species - 1)")
  } else {
    if (length(b_mu) == 1) {
      b_mu <- rep(b_mu, n_off)
    }
  }
  if (length(b_sd) %in% c(1, n_off) == FALSE) {
    stop("Error: sd on off diagonal for B (b_sd) needs to be a scalar or vector whose length is n_species * (n_species - 1)")
  } else {
    if (length(b_sd) == 1) {
      b_sd <- rep(b_sd, n_off)
    }
  }

  if (length(off_diag_prior) > 1) off_diag_prior <- off_diag_prior[1]
  if (off_diag_prior %in% c("normal", "student-t", "laplace", "hs") == FALSE) {
    stop("Error: off_diag_prior must be one of 'normal','student-t','laplace','hs'")
  }
  off_diag_priors <- match(off_diag_prior, c("normal", "student-t", "laplace", "hs"))
  off_diag_priors <- off_diag_priors - 1 # index at 0

  nu_known <- 0
  if (off_diag_priors == 1) {
    if (!is.null(nu)) {
      if (nu < 2) stop("Error: if nu is entered, it must be >= 2")
      nu_known <- nu
    }
  }

  data_list <- list(
    n_pos = n_pos,
    n_time = n_time,
    n_spp = n_spp,
    n_off = n_spp * (n_spp - 1),
    n_q = n_q,
    n_r = n_r,
    id_q = id_q,
    id_r = id_r,
    time_index = data$time,
    species_index = data$species,
    yy = data$y,
    rc_off = rc_off,
    b_mu = b_mu,
    b_sd = b_sd,
    b_mu_diag = b_mu_diag,
    b_sd_diag = b_sd_diag,
    nu_known = nu_known,
    sigma_proc_sd = prior_sigma_pro_sd,
    sigma_proc_mu = prior_sigma_pro_mu,
    sigma_obs_sd = prior_sigma_obs_sd,
    sigma_obs_mu = prior_sigma_obs_mu,
    off_diag_priors = off_diag_priors,
    sigma_scale_df = sigma_scale_df,
    sigma_scale_sd = sigma_scale_sd,
    global_scale = global_scale,
    slab_scale = slab_scale,
    slab_df = slab_df,
    priors_only = as.numeric(prior_only)
  )

  pars <- c("sigma_proc", "sigma_obs", "Bmat", "x0")
  if(save_log_lik==TRUE) pars = c(pars,"log_lik")
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
    pars <- c(pars, "c2_hs", "lambda")
  }

  if (!is.null(pars_to_monitor)) {
    pars <- pars_to_monitor
    # if all then monitor everything
    if (pars_to_monitor == "all") {
      pars <- NULL
    }
  }


  if (map_estimation == FALSE) {
    out <- rstan::sampling(
      object = stanmodels$varlasso,
      data = data_list,
      pars = pars,
      warmup = warmup,
      iter = iter,
      thin = thin,
      chains = chains, ...
    )
  } else {
    sampling_args <- list(
      object = stanmodels$varlasso,
      data = data_list,
      #verbose = verbose,
      hessian = est_hessian,
      ...
    )
    out <- do.call(optimizing, sampling_args)
  }

  return(list(fit=out, data = data_list, pars = pars))
}
