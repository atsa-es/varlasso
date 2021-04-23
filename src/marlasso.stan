data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=0> n_off;            // # non-zero off-diagonals
  int<lower=1> n_q;              // # unique proc SD
  int<lower=1> n_obs;            // # observed ts
  // vectors
  int id_q[n_spp];                 // IDs for proc SD
  // matrices
  int<lower=0> rc_off[n_off,2];  // indices of non-zero off-diags
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos];
  int<lower=0> col_indx_pos[n_pos];
  real yy[n_pos];       // data
  real pro_mu; // mean of process variance prior
  real pro_cv; // cv of process variance prior
  real obs_mu; // mean of obs variance prior
  real obs_cv; // cv of obs variance prior
  real b_mu[6]; // prior on mean of offdiagonal
  real b_sd[6]; // prior on sd offdiagonal
  real b_mu_diag[4];// prior on mean of diagonal
  real b_sd_diag[4];// prior on sd offdiagonal
}
parameters {
  real<lower=0> SD_obs;
  vector<lower=-1,upper=1>[n_off] Boffd;  // off-diags of B
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  vector[n_spp] xx0;
  real<lower=0> SD_proc;           // proc SD
  matrix[n_spp,n_year-1] devs;       // states
}
transformed parameters {
  // obs SD
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  matrix[n_obs,n_year] yymiss;
  matrix[n_spp,n_year] xx;       // states
  // diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }

  xx[,1] = xx0;
  for(t in 2:n_year) {
    xx[,t] = Bmat * xx[,t-1] + SD_proc*devs[,t-1];
    // col(xx,t) ~ multi_normal(Bmat * col(xx,t-1), QQ);
  }

  // Deal with missing and non-missing values separately
  for(i in 1:n_pos) {
    yymiss[row_indx_pos[i], col_indx_pos[i]] = yy[i];
  }
}
model {
  // PRIORS
  // initial state
  xx0 ~ normal(0,1);
  for(i in 1:n_spp) {
    devs[i] ~ std_normal();
  }
  // process SD's
  SD_proc ~ normal(pro_mu,pro_mu*pro_cv);
  // obs SD
  SD_obs ~ normal(obs_mu,obs_mu*obs_cv);
  // B diagonal
  for(i in 1:4) Bdiag[i] ~ normal(b_mu_diag[i],b_sd_diag[i]);
  // B off-diagonals
  for(i in 1:6) Boffd[i] ~ normal(b_mu[i],b_sd[i]);
  // LIKELIHOOD
  //for(t in 2:n_year) {
  //  xx[,t] ~ normal(Bmat * xx[,t-1], SD_proc);
    // col(xx,t) ~ multi_normal(Bmat * col(xx,t-1), QQ);
  //}
  to_vector(yymiss) ~ normal(to_vector(xx), SD_obs);
}
