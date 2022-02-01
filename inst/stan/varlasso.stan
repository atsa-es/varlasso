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
  real b_mu[n_off]; // prior on mean of offdiagonal
  real b_sd[n_off]; // prior on sd offdiagonal
  real b_mu_diag[n_spp];// prior on mean of diagonal
  real b_sd_diag[n_spp];// prior on sd of diagonal
  real fixed_r[n_spp]; // fixed_r is optional, will be all 0s if estimated
  int off_diag_priors; // 0 if normal, 1 if Student-t, 2 if Laplace, 3 if horseshoe
  real<lower=0> sigma_proc_mu; // optional
  real<lower=0> sigma_obs_mu; // optional
  real<lower=0> sigma_proc_sd; // optional
  real<lower=0> sigma_obs_sd; // optional
  real<lower=0> nu_known; // optional
  real<lower=0> global_scale;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  real<lower=0> sigma_scale_df;
  real<lower=0> sigma_scale_sd;
  int priors_only; // whether to sample from priors only
  int est_trend; // whether to estimate a trend for each species
}
transformed data {
  int est_nu; // whether to estimate student-t parameters
  int est_hs; // whether to estimate hs parameters
  int est_lambda; // whether to estimate laplace/student-t prior parameters
  int est_sigma_obs;
  real dummy;
  // initialize
  est_nu = 0;
  est_hs = 0;
  est_lambda = 0;
  est_sigma_obs = 1; // by default, estimate it
  dummy = 0; // for sampling from priors
  // indicators
  if(off_diag_priors==1) {
    est_nu = 1;
    if(nu_known > 0) est_nu = 0;
  }
  if(off_diag_priors==1) est_lambda = 1;
  if(off_diag_priors==2) est_lambda = 1;
  if(off_diag_priors==3) est_hs = 1;
  for(i in 1:n_spp) {
    if(fixed_r[i] != 0) est_sigma_obs = 0;
  }
}
parameters {
  real<lower=0> sigma_obs[est_sigma_obs * n_r];
  real<lower=0> sigma_scale; // variance for shrinkage / hierarchical B off diags
  vector[n_off] B_z;  // off-diags of B, in normal (0,1) space
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  vector[n_spp] x0;
  vector[est_trend*n_spp] U;
  vector<lower=2>[est_nu] nu; // student-t df parameters for Student-t priors
  vector<lower=0>[est_lambda*n_off] lambda2; // parameters for Student-t and laplace priors
  vector<lower=0>[est_hs] c2_hs; // c parameters for horseshoe priors
  real<lower=0> lambda[est_hs*n_off]; // parameters for horseshoe priors
  real<lower=0> sigma_proc[n_q];           // proc SD
  matrix[n_spp,n_time-1] devs;       // states
}
transformed parameters {
  matrix[n_spp,n_spp] Bmat;
  matrix[n_spp,n_time] x;       // states
  real<lower=0> lambda_tilde[est_hs*n_off]; // parameters for horseshoe priors
  vector<lower=0>[n_spp] sigma;
  vector<lower=0>[n_r] sigma_r;
  vector[n_off] Boffd;  // off-diags of B

  // B off-diagonals
  if(off_diag_priors == 0) {
     //Normal priors
     //sigma_scale ~ student_t(3,0,2); // not estimated
     for(i in 1:n_off) {
        Boffd[i] = B_z[i]*b_sd[i] + b_mu[i];//~ normal(b_mu[i], b_sd[i]);
     }
  }
  if(off_diag_priors == 1) {
     //Student t priors
     for(i in 1:n_off) {
        Boffd[i] = sigma_scale * sqrt(lambda2[i]) * B_z[i];
        //Boffd[i] ~ normal(0, sigma_scale*sqrt(1/inv_lambda2[i]));
     }
  }
  if(off_diag_priors == 2) {
    for(i in 1:n_off) {
      // this is implementation by Jeffrey Arnold
      Boffd[i] = sigma_scale * sqrt(lambda2[i]) * B_z[i];
      //Boffd[i] = B_z[i];
     }
  }
  if(off_diag_priors == 3) {
    for(i in 1:n_off) {
      lambda_tilde[i] = sqrt((c2_hs[1] * lambda[i] * lambda[i]) / (c2_hs[1] + square(sigma_scale * lambda[i])));
      //lambda_tilde[i] = sqrt(c2[1])*lambda[i] / sqrt(c2[1] + sigma_scale*lambda[i]*lambda[i]);
      Boffd[i] = B_z[i]*sigma_scale*lambda_tilde[i]; //normal(0, sigma_scale*lambda_tilde[i]);
    }
  }
  if(off_diag_priors == 4) {
     //Normal priors with sd estimated
     for(i in 1:n_off) {
        Boffd[i] = B_z[i]*sigma_scale + b_mu[i];//~ normal(b_mu[i], sigma_scale);
     }
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
    if(off_diag_priors> 0) sigma[i] = sigma[i] * sigma_scale;
  }
  // this is autoregression equation
  x[,1] = x0;
  for(t in 2:n_time) {
    x[,t] = Bmat * x[,t-1] + sigma .* devs[,t-1];
    if(est_trend == 1) {
      x[,t] = x[,t] + U;
    }
  }

  for(i in 1:n_r) {
    if(est_sigma_obs==1) {
      sigma_r[i] = sigma_obs[i];
    } else {
      sigma_r[i] = fixed_r[i];
    }
  }
}
model {
  // PRIORS
  // initial state
  x0 ~ std_normal();
  for(i in 1:n_spp) {
    devs[i] ~ std_normal();
  }
  if(est_trend==1) {
    U ~ std_normal();
  }
  // process SD's
  sigma_proc ~ student_t(3, sigma_proc_mu, sigma_proc_sd);
  // obs SD
  sigma_obs ~ student_t(3, sigma_obs_mu, sigma_obs_sd);
  // B diagonal
  Bdiag ~ normal(b_mu_diag,b_sd_diag);

  if(off_diag_priors == 0) {
    // std normal prior on B_z,
    B_z ~ std_normal();
  }
  if(off_diag_priors == 1) {
    B_z ~ std_normal();
    //Student t priors
    sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    if(est_nu==1) {
      nu[1] ~ gamma(2, 0.1);
      lambda2 ~ inv_gamma(nu[1]/2, nu[1]/2);
      //for(i in 1:n_off) {
      //  lambda2[i] ~ inv_gamma(nu[1]/2, nu[1]/2); // vector of local variances
        //B_z[i] ~ student_t(nu[1], 0, sigma_scale);
      //}
    } else {
      lambda2 ~ inv_gamma(nu_known/2, nu_known/2);
      //for(i in 1:n_off) {
        //B_z[i] ~ student_t(nu_known, 0, sigma_scale);
        //lambda2[i] ~ inv_gamma(nu_known/2, nu_known/2); // vector of local variances
      //}
    }
  }
  if(off_diag_priors == 2) {
    B_z ~ std_normal();
    sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
    //B_z ~ double_exponential(0, sigma_scale);
    lambda2 ~ exponential(0.5); // vector
  }
  if(off_diag_priors == 3) {
    B_z ~ std_normal();
    // regularized horseshore prior in rstanarm
    // see 2.8 in Piironen and Vehtari 2017
    // global_scale argument equal to the ratio of the expected number of non-zero
    //coefficients to the expected number of zero coefficients, divided by the square
    //root of the number of observations
    //hs(df = 1, global_df = 1, global_scale = 0.01, slab_df = 4, slab_scale = 2.5)
    lambda ~ student_t(1, 0.0, 1.0); // df_lambda can also be passed in
    sigma_scale ~ student_t(3, 0, global_scale);//cauchy(0, global_scale);
    c2_hs ~ inv_gamma(slab_df/2.0, slab_df*slab_scale*slab_scale/2.0);
    //c_hs ~ student_t(slab_df, 0, slab_scale);
  }
  if(off_diag_priors == 4) {
    // std normal prior on B_z,
    B_z ~ std_normal();
    sigma_scale ~ student_t(sigma_scale_df,0,sigma_scale_sd);
  }

  // LIKELIHOOD
  if(priors_only==0) {
    for(i in 1:n_pos) {
      yy[i] ~ normal(x[species_index[i],time_index[i]], sigma_r[id_r[species_index[i]]]);
    }
  } else {
    dummy ~ normal(0,1);
  }

}
generated quantities {
  vector[n_pos] log_lik;
  if(priors_only==0) {
    for(i in 1:n_pos) {
      log_lik[i] = normal_lpdf(yy[i] | x[species_index[i],time_index[i]], sigma_r[id_r[species_index[i]]]);
    }
  }
}
