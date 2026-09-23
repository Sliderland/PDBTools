// Local baked-prior variant of semiconj.stan.
// The original file remains unchanged; prior hyperparameters are fixed locally.
data {
  int<lower=1> m; // Dimension of observation vector
  int<lower=1> p; // Order of VAR model
  int<lower=p> N; // Length of time series; permits an empty likelihood at N = p
  array[N] vector[m] y; // Time series
}
transformed data {
  int n_cond = N - p;
  int n_offdiag = m * (m - 1);
  matrix[n_cond, m * p] X;
  matrix[n_cond, m] Y;
  array[n_offdiag] int offdiag_idx;
  int k = 1;

  // Rows of X hold lagged observations in lag-major blocks. Thus if each
  // block of B is phi[lag]', X * B gives sum_lag phi[lag] * y[t - lag].
  if (n_cond > 0) {
    for (n in 1:n_cond) {
      Y[n] = y[p + n]';
      for (lag in 1:p) {
        X[n, ((lag - 1) * m + 1):(lag * m)] = y[p + n - lag]';
      }
    }
  }

  // Indices into to_vector(phi[s]) for all off-diagonal matrix entries.
  // Matrix vectorization is column-major, matching (column - 1) * m + row.
  if (m > 1) {
    for (j in 1:m) {
      for (i in 1:m) {
        if (i != j) {
          offdiag_idx[k] = (j - 1) * m + i;
          k += 1;
        }
      }
    }
  }
}
parameters {
  array[p] matrix[m, m] phi; // The phi_i
  cov_matrix[m] Sigma; // Error variance, Sigma
}
model {
  matrix[m * p, m] B = rep_matrix(0.0, m * p, m);
  matrix[m, m] L_Sigma = cholesky_decompose(Sigma);
  matrix[n_cond, m] residual;
  matrix[m, n_cond] whitened;

  // Assemble coefficient blocks once, then evaluate the conditional VAR
  // likelihood for every t > p with a single matrix product and solve.
  for (lag in 1:p) {
    B[((lag - 1) * m + 1):(lag * m)] = phi[lag]';
  }
  residual = Y - X * B;
  whitened = mdivide_left_tri_low(L_Sigma, residual');
  // Match Stan's proportional multi_normal target by omitting its
  // parameter-independent Gaussian normalizing constant.
  target += - n_cond * sum(log(diagonal(L_Sigma)))
            - 0.5 * sum(square(whitened));

  Sigma ~ inv_wishart(m + 4, diag_matrix(rep_vector(1.0, m)));
  for (lag in 1:p) {
    diagonal(phi[lag]) ~ normal(1.0, 10.0);
    if (m > 1) {
      to_vector(phi[lag])[offdiag_idx] ~ normal(0.0, 1.0);
    }
  }
}
generated quantities {
  matrix[m, m * p] topblock;
  for (lag in 1:p) {
    topblock[:, ((lag - 1) * m + 1):(lag * m)] = phi[lag];
  }
  matrix[m * p, m * p] companion = rep_matrix(0.0, m * p, m * p);
  companion[1:m, :] = topblock;
  if (p > 1) {
    companion[(m + 1):(m * p), 1:(m * (p - 1))] =
      diag_matrix(rep_vector(1.0, m * (p - 1)));
  }
  complex_vector[m * p] lambdas = eigenvalues(companion);
  vector[m * p] lambda_moduli = abs(lambdas);
  real max_lambda_modulus = max(lambda_moduli);
}
