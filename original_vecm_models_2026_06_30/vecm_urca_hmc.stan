data {
  int<lower=1> t; //no. of observations
  int<lower=1> N; //no. of variables
  int<lower=1> p; //no. of lags
  int<lower=1, upper=N - 1> h; //no. of cointegrations
  int<lower=0> num_exo;
  matrix[t, N] y; //all time series data
  matrix[t, num_exo] D;
}
transformed data {
  int pd = p - 1;
  int K = t - p;
  matrix[K, N] dy;
  matrix[K, N] y_lagp;
  matrix[K, N * pd] dX; //lag difference matrix
  matrix[K, num_exo] D_trunc;
  for (i in 1 : K) {
    int t_idx = i + p;
    dy[i,  : ] = y[t_idx,  : ] - y[t_idx - 1,  : ];
    y_lagp[i,  : ] = y[t_idx - p,  : ];
    D_trunc[i,  : ] = D[i + p,  : ];
    for (lag in 1 : pd) {
      dX[i, ((lag - 1) * N + 1) : (lag * N)] = y[t_idx - lag,  : ]
                                               - y[t_idx - lag - 1,  : ];
    }
  }
}
parameters {
  vector[N] mu;
  // cholesky_factor_cov[N] L;
  cholesky_factor_corr[N] L_Omega; // Better behavior than raw covariance
  vector<lower=0>[N] L_sigma;
  array[pd] matrix[N, N] gamma;
  matrix[N, num_exo] phi;
  matrix[N - h, h] B;
  matrix[N, h] alpha;
}
transformed parameters {
  matrix[h, N] beta;
  beta[1 : h, 1 : h] = diag_matrix(rep_vector(1.0, h));
  beta[1 : h, h + 1 : N] = -1 * B';
  matrix[N, N * pd] topblock;
  for (lag in 1 : pd) {
    topblock[ : , ((lag - 1) * N + 1) : (lag * N)] = gamma[lag];
  }
  cholesky_factor_cov[N] L = diag_pre_multiply(L_sigma, L_Omega);
  matrix[N, N] pi = beta' * alpha';
  matrix[K, N] ect = y_lagp * pi; // Cointegration/Error Correction part
  matrix[K, N] short_run = dX * topblock';
}
model {
  to_vector(topblock) ~ normal(0, 5);
  to_vector(B) ~ normal(0, 10);
  L_Omega ~ lkj_corr_cholesky(2.0); // Prior favoring low-to-moderate correlation
  L_sigma ~ exponential(0.01);
  
  for (i in 1 : K) {
    vector[N] current_mu = mu + short_run[i,  : ]' + ect[i,  : ]'
                           + (phi * D_trunc[i,  : ]');
    dy[i,  : ]' ~ multi_normal_cholesky(current_mu, L);
  }
}
generated quantities {
  array[p] matrix[N, N] A;
  A[1] = diag_matrix(rep_vector(1.0, N)) + pi + gamma[1];
  for (i in 2 : (p - 1)) {
    A[i] = gamma[i] - gamma[i - 1];
  }
  A[p] = -gamma[pd];
  matrix[N, N * p] inLevelsTop;
  for (lag in 1 : p) {
    inLevelsTop[ : , ((lag - 1) * N + 1) : (lag * N)] = A[lag];
  }
  
  matrix[N * p, N * p] inLevelsCompanion = rep_matrix(0.0, N * p, N * p);
  inLevelsCompanion[(N + 1) : (N * p), 1 : (N * (p - 1))] = diag_matrix(
                                                                    rep_vector
                                                                    (1.0,
                                                                    N
                                                                    * (
                                                                    p - 1)));
  complex_vector[N * p] lambdas = eigenvalues(inLevelsCompanion);
  vector[N * p] lambda_moduli = abs(lambdas);
  real max_lambda_modulus = max(lambda_moduli);
}
