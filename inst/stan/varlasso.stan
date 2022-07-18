functions {
  /* Efficient computation of the horseshoe prior
   * (generated with brms 2.16.3)
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slap regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
}
data {
  int<lower=0> n_time;           // years
  int<lower=0> n_spp;            // species
  int<lower=0> n_off;            // non-zero off-diagonals
  int<lower=1> n_q;              // unique proc SD
  int<lower=1> n_r;              // unique obs SD
  int id_q[n_spp];                 // IDs for proc SD
  int id_r[n_spp];                 // IDs for obs SD
  int<lower=0> rc_off[n_off,2];  // indices of non-zero off-diags
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> time_index[n_pos];
  int<lower=0> species_index[n_pos];
  real yy[n_pos];       // data
  //real b_mu[n_off]; // prior on mean of offdiagonal
  real b_sd[n_off]; // prior on sd offdiagonal
  real b_diag_min; // minimum value of estimated B diag
  real b_mu_diag[n_spp];// prior on mean of diagonal
  real b_sd_diag[n_spp];// prior on sd of diagonal
  real fixed_r[n_spp]; // fixed_r is optional, will be all 0s if estimated
  int off_diag_priors; // 0 if normal, 1 if Student-t, 2 if Laplace, 3 if horseshoe
  real<lower=0> sigma_proc_mu; // optional
  real<lower=0> sigma_obs_mu; // optional
  real<lower=0> sigma_proc_sd; // optional
  real<lower=0> sigma_obs_sd; // optional
  real<lower=0> nu_known; // optional
  real<lower=0> sigma_scale_df;
  real<lower=0> sigma_scale_sd;
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
  // options
  int priors_only; // whether to sample from priors only
  int est_trend; // whether to estimate a trend for each species
  int process_only; // whether to fit a VAR and estimate process error only
  int est_unique_reg;
}
transformed data {
  int est_nu; // whether to estimate student-t parameters
  int est_hs; // whether to estimate hs parameters
  int est_lambda; // whether to estimate laplace/student-t prior parameters
  int est_sigma_obs;
  int unique_reg;
  real dummy;
  int est_sigma_scale;
  matrix[n_spp,n_time] ymat;  // matrix of y for VAR model
  // initialize
  est_nu = 0;
  est_hs = 0;
  est_lambda = 0;
  est_sigma_obs = 1; // by default, estimate it
  dummy = 0; // for sampling from priors
  est_sigma_scale = 0;
  // indicators
  if(off_diag_priors==1) {
    est_nu = 1;
    if(nu_known > 0) est_nu = 0;
    est_lambda = 1;
    est_sigma_scale = 1;
  }
  if(off_diag_priors==2) est_lambda = 1;
  if(off_diag_priors==3) est_hs = 1;
  if(off_diag_priors==4) est_sigma_scale = 1;

  for(i in 1:n_spp) {
    if(fixed_r[i] != 0) est_sigma_obs = 0;
  }
  if(process_only==1) est_sigma_obs = 0;

  // convert vector to matrix for VAR
  for(i in 1:n_spp) {
    for(j in 1:n_time) {ymat[i,j]=0;}
  }
  for(i in 1:n_pos) {
    ymat[species_index[i],time_index[i]] = yy[i];
  }
  // regularization parameter fixed at 1, but more flexible aproach
  // is to estimate coefficient specific values
  unique_reg = 1;
  if(est_unique_reg==1) {
    unique_reg = n_off;
  }

}
parameters {
  real<lower=0> sigma_obs[est_sigma_obs * n_r];
  vector<lower=0>[est_sigma_scale] sigma_scale; // variance for shrinkage / hierarchical B off diags
  vector[n_off] B_z;  // off-diags of B, in normal (0,1) space
  vector<lower=b_diag_min,upper=1>[n_spp] Bdiag;   // diag of B
  vector[n_spp] x0;
  vector[est_trend*n_spp] U;
  vector<lower=2>[est_nu] nu; // student-t df parameters for Student-t priors
  vector<lower=0>[est_lambda*unique_reg] lambda2; // parameters for Student-t and laplace priors
  real<lower=0> sigma_proc[n_q];           // proc SD
  matrix[n_spp*(1-process_only),(n_time-1)*(1-process_only)] devs;       // states
  // paramters specific to hs
  vector<lower=0>[est_hs*n_off] hs_local; // parameters for horseshoe priors
  real<lower=0> hs_global[est_hs];  // global shrinkage parameters
  real<lower=0> hs_slab[est_hs];  // slab regularization parameter
}
transformed parameters {
  matrix[n_spp,n_spp] Bmat;
  matrix[n_spp,n_time] x;       // states
  vector<lower=0>[n_spp] sigma;
  vector<lower=0>[n_r * (1-process_only)] sigma_r;
  vector[n_off] Boffd;  // off-diags of B

  // B off-diagonals
  if(off_diag_priors == 0) {
     //Normal priors
     //sigma_scale ~ student_t(3,0,2); // not estimated
     for(i in 1:n_off) {
       Boffd[i] = B_z[i] * b_sd[i];
     }
  }
  if(off_diag_priors == 1) {
    if(est_unique_reg ==1) {
      Boffd = sigma_scale[1] * sqrt(lambda2) .* B_z;
    } else {
      Boffd = sigma_scale[1] * sqrt(lambda2[1]) * B_z;
    }
  }
  if(off_diag_priors == 2) {
    if(est_unique_reg == 1) {
      Boffd = sqrt(lambda2) .* B_z;
    } else {
      Boffd = sqrt(lambda2[1]) * B_z;
    }
  }
  if(off_diag_priors == 3) {
    Boffd = horseshoe(B_z, hs_local, hs_global[1], hs_scale_slab^2 * hs_slab[1]);
  }
  if(off_diag_priors == 4) {
     //Normal priors with sd estimated
     Boffd = B_z*sigma_scale[1];
  }

  // construct B matrix, starting with diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // this maps shared variances
  for(i in 1:n_spp) {
    sigma[i] = sigma_proc[id_q[i]];
    //if(off_diag_priors> 0) sigma[i] = sigma[i] * sigma_scale;
  }
  // this is autoregression equation for VARSS
  if(process_only == 0) {
    x[,1] = x0;
    for(t in 2:n_time) {
      x[,t] = Bmat * x[,t-1] + sigma .* devs[,t-1];
      if(est_trend == 1) {
        x[,t] = x[,t] + U;
      }
    }
    // map potentially shared observation errors
    for(i in 1:n_r) {
      if(est_sigma_obs==1) {
        sigma_r[i] = sigma_obs[i];
      } else {
        sigma_r[i] = fixed_r[i];
      }
    }
  } else {
    x[,1] = x0;
    for(t in 2:n_time) {
      x[,t] = Bmat * ymat[,t-1];
      if(est_trend == 1) {
        x[,t] = x[,t] + U;
      }
    }
  }

}
model {
  // PRIORS
  x0 ~ std_normal();  // initial state
  if(process_only==0) {
    for(i in 1:n_spp) {
      devs[i] ~ std_normal(); // process devs
    }
    sigma_obs ~ student_t(3, sigma_obs_mu, sigma_obs_sd);  // obs SD
  }
  if(est_trend==1) {
    U ~ std_normal(); // linear trend
  }

  sigma_proc ~ student_t(3, sigma_proc_mu, sigma_proc_sd);  // process SD's
  Bdiag ~ normal(b_mu_diag,b_sd_diag);  // B diagonal

  if(off_diag_priors == 0) {
    // this model treats B elements as fixed effects
    B_z ~ std_normal();// std normal prior on B_z,
  }
  if(off_diag_priors == 1) {
    // this model treats B elements as random with shrinkage lambda
    B_z ~ std_normal();
    //Student t priors
    sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    if(est_nu==1) {
      nu[1] ~ gamma(2, 0.1);
      lambda2 ~ inv_gamma(nu[1]/2, nu[1]/2);
    } else {
      lambda2 ~ inv_gamma(nu_known/2, nu_known/2);
    }
  }
  if(off_diag_priors == 2) {
    // parameterization from brms
    //lambda2 ~ chi_square(sigma_scale_df);
    //B_z ~ double_exponential(0, lambda2[1] * sigma_scale_sd);

    // this parameterization is in Stan manual, Ding and Blitzstein 2018
    // a ~ exp(1/(2*sigma2)) B|a ~ N(0, sqrt(a))
    B_z ~ std_normal();
    //sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    lambda2 ~ exponential(0.5 * (1.0 / pow(sigma_scale_sd,2)));

    //B_z ~ std_normal();
    //sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    //lambda2 ~ exponential(0.5);
    //     // vector
    //sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    //B_z ~ double_exponential(0, sigma_scale);
  }
  if(off_diag_priors == 3) {
    B_z ~ std_normal();
    hs_local ~ student_t(hs_df, 0, 1);
    // this is scaled by residual sd in brms when autoscale=T (default)
    if(process_only==0) {
      hs_global ~ student_t(hs_df_global, 0, hs_scale_global * (sigma_proc[1]));
    } else {
      hs_global ~ student_t(hs_df_global, 0, hs_scale_global);
    }
    hs_slab ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab);
  }
  if(off_diag_priors == 4) {
    // std normal prior on B_z,
    B_z ~ std_normal();
    sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
  }

  // LIKELIHOOD
  if(priors_only==0) {
    if(process_only==0) {
      // VARSS
      for(i in 1:n_pos) {
        yy[i] ~ normal(x[species_index[i],time_index[i]], sigma_r[id_r[species_index[i]]]);
      }
    } else {
      // VAR
      for(i in 1:n_pos) {
        yy[i] ~ normal(x[species_index[i],time_index[i]], sigma[species_index[i]]);
      }
    }
  } else {
    dummy ~ normal(0,1);
  }

}
generated quantities {
  vector[n_pos] log_lik;
  if(priors_only==0) {
    if(process_only==0) {
      for(i in 1:n_pos) {
        log_lik[i] = normal_lpdf(yy[i] | x[species_index[i],time_index[i]], sigma_r[id_r[species_index[i]]]);
      }
    } else {
      for(i in 1:n_pos) {
        log_lik[i] = normal_lpdf(yy[i] | x[species_index[i],time_index[i]], sigma[species_index[i]]);
      }
    }
  }
}
