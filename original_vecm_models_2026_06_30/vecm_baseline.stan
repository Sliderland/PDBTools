data{
    int<lower=1> t; //no. of observations
    int<lower=1> N; //no. of variables
    int<lower=1> p; //no. of lags
    int<lower=1, upper=N> h; //no. of cointegrations
    matrix[t, N] y; //all time series data
}
transformed data {
    int pd = p - 1;
    int K = t - p;
    matrix[K, N] dy;
    matrix[K, N] y_lag1;
    matrix[K, N*pd] dX; //lag difference matrix
    for(i in 1:K){
        int t_idx = i + p;
        dy[i, ] = y[t_idx, ] - y[t_idx - 1, ];
        y_lag1[i, ] = y[t_idx - 1, ];

        for(lag in 1:pd){
            dX[i, ((lag-1)*N+1):(lag*N)] = y[t_idx - lag, ] - y[t_idx - lag - 1, ];
        }
    }

    // matrix[t-1, N] one_lag = y[2:t, ];
    // matrix[t-1, N] first_diff = y[2:t, ] - y[1:(t-1), ];
    // int<lower=0> pd = p-1; //no. of difference lags
    // matrix[t-pd-1, N*pd] diff_lags;
    // for(lag in 1:pd){
    //     diff_lags[, ((lag-1)*N+1):(lag*N)] = first_diff[(pd-lag+1):(t-lag-1), ];
    // }
}
parameters {
    vector[N] mu;
    cholesky_factor_cov[N] L;
    array[pd] matrix[N, N] xi;
    // matrix[N, N] rho;
    matrix[N-h, h] B;
    matrix[N, h] alpha;
}
transformed parameters{
    matrix[h, N] beta;
    beta[1:h, 1:h] = diag_matrix(rep_vector(1.0, h));
    beta[1:h, h+1:N] = -1*B';
    matrix[N, N*pd] topblock;
    for(lag in 1:pd){
        topblock[, ((lag-1)*N+1):(lag*N)] = xi[lag];
    }
    matrix[N, N] pi = beta' * alpha';
    matrix[K, N] ect = y_lag1 * pi; // Cointegration/Error Correction part
    matrix[K, N] short_run = dX * topblock';
}
model{
    // to_vector(alpha) ~ normal(0, 100);
    // to_vector(topblock) ~ normal(0, );
    to_vector(B) ~ normal(0, 3);
    for(i in 1:K){
        // dy[i, ]' ~ multi_normal_cholesky(mu' + ect[i, ] + (dX*topblock')[i,1:N], L);
        // y[i, ] ~ multi_normal_cholesky(y_hat[i, ], L);
        vector[N] current_mu = mu + short_run[i, ]' + ect[i, ]';
        dy[i, ]' ~ multi_normal_cholesky(current_mu, L);
        // dy[i, ] ~ multi_normal_cholesky(y_hat[i, ], L);
        // diff_lags[i, 1:N]' ~ multi_normal_cholesky(mu + (topblock*diff_lags[i-pd,]')[1:N] + b * a' * y[i-1, ]', L);
        // y[i, ]' ~ multi_normal_cholesky(mu + (topblock*diff_lags[i-pd,]')[1:N] + b * a' * y[i-1, ]', L);
    }
}
generated quantities{
  array[p] matrix[N, N] A;
  A[1] = diag_matrix(rep_vector(1.0, N)) + pi + xi[1];
  for (i in 2 : (p - 1)) {
    A[i] = xi[i] - xi[i - 1];
  }
  A[p] = -xi[pd];
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
//     // matrix[t-pd-1, N*pd] diff_lag_matrix = diff_lags;
//     // array[p] matrix[N, N] phi;
//     // phi[1] = rho + xi[1];
//     // // phi[p] = xi[pd];
//     // for(lag in 2:p){
//     //     phi[lag] = xi[lag];
//     // }
//     // matrix[N, N] xi0 = rho - diag_matrix(rep_vector(1.0, N));
}
