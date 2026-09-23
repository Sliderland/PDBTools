// Local baked-prior variant of statrmlprior.stan.
// The original file remains unchanged; prior hyperparameters are fixed locally.
functions {
  /* Function to compute the matrix square root */
  matrix sqrtm(matrix A) {
    int m = rows(A);
    // Use the principal symmetric square root. The extra square root in
    // `root_root_evals` is required because `tcrossprod` squares the factors.
    matrix[m, m] A_sym = 0.5 * (A + A');
    vector[m] root_root_evals = sqrt(sqrt(eigenvalues_sym(A_sym)));
    matrix[m, m] evecs = eigenvectors_sym(A_sym);
    matrix[m, m] eprod = diag_post_multiply(evecs, root_root_evals);
    return tcrossprod(eprod);
  }
  /* Function to perform Algorithm [VQ] from Roy et al. (2019).
     Returned: a (2 x p) array of (m x m) matrices; the (1, i)-th component
               of the array is phi_i and the (2, i)-th component of the array
               is Gamma_{i-1} (assuming M = Sigma) */
  array[,] matrix rev_mapping(array[] matrix C, array[] matrix V, matrix M) {
    int p = size(V);
    int m = rows(M);
    array[p + 1] matrix[m, m] U;
    array[p + 1] matrix[m * (p + 1), m * (p + 1)] Uunder;
    array[p] matrix[m * p, m] eps;
    array[p] matrix[m * p, m] kappa;
    array[p + 1] matrix[m, m] D;
    matrix[m, p * m] A;
    array[2, p] matrix[m, m] phiGamma;
    U[1] = M;
    for (i in 1 : p) {
      U[1] = U[1] + V[i];
    } // U(0)
    Uunder[1][1 : m, 1 : m] = U[1];
    D[1] = U[1];
    U[2] = C[1]' * sqrtm(U[1]); // U(1), etc.
    eps[1][1 : m,  : ] = U[2]';
    kappa[1][1 : m,  : ] = U[2];
    D[2] = U[1]
           - kappa[1][1 : m,  : ]'
             * mdivide_left_spd(Uunder[1][1 : m, 1 : m],
                                kappa[1][1 : m,  : ]);
    Uunder[2][1 : m, 1 : m] = U[1];
    Uunder[2][1 : m, (m + 1) : (2 * m)] = eps[1][1 : m,  : ]';
    Uunder[2][(m + 1) : (2 * m), 1 : m] = eps[1][1 : m,  : ];
    Uunder[2][(m + 1) : (2 * m), (m + 1) : (2 * m)] = Uunder[1][1 : m, 1 : m];
    if (p > 1) {
      for (i in 2 : p) {
        int end = (i - 1) * m;
        U[i + 1] = eps[i - 1][1 : end,  : ]'
                   * mdivide_left_spd(Uunder[i - 1][1 : end, 1 : end],
                                      kappa[i - 1][1 : end,  : ])
                   + C[i]' * sqrtm(D[i]);
        eps[i][1 : end,  : ] = eps[i - 1][1 : end,  : ];
        eps[i][(end + 1) : (end + m),  : ] = U[i + 1]';
        kappa[i][1 : m,  : ] = U[i + 1];
        kappa[i][(m + 1) : (end + m),  : ] = kappa[i - 1][1 : end,  : ];
        D[i + 1] = U[1]
                   - kappa[i][1 : (end + m),  : ]'
                     * mdivide_left_spd(Uunder[i][1 : (end + m), 1 : (
                                        end + m)],
                                        kappa[i][1 : (end + m),  : ]);
        Uunder[i + 1][1 : m, 1 : m] = U[1];
        Uunder[i + 1][1 : m, (m + 1) : (end + 2 * m)] = eps[i][1 : (end + m),  : ]';
        Uunder[i + 1][(m + 1) : (end + 2 * m), 1 : m] = eps[i][1 : (end + m),  : ];
        Uunder[i + 1][(m + 1) : (end + 2 * m), (m + 1) : (end + 2 * m)] = Uunder[i][1 : (end
                                                                    + m), 1 : (end
                                                                    + m)];
      }
    }
    A = mdivide_right_spd(eps[p]', Uunder[p][1 : (p * m), 1 : (p * m)]);
    for (i in 1 : p) 
      phiGamma[1, i] = A[ : , ((i - 1) * m + 1) : (i * m)];
    for (i in 1 : p) 
      phiGamma[2, i] = Uunder[p][1 : m, ((i - 1) * m + 1) : (i * m)]';
    return phiGamma;
  }
}
data {
  int<lower=1> m; // Dimension of observation vector
  int<lower=1> p; // Order of VAR model
  int<lower=1> N; // Length of time series
  array[N] vector[m] y; // Time series
}
transformed data {
  matrix[m, m] identity; // Identity matrix
  real df = m + 4;
  vector[p * m] y1top;
  matrix[N - p, p * m] lag_design;
  matrix[N - p, m] y_rest;
  for (t in 1 : p)
    y1top[((t - 1) * m + 1) : (t * m)] = y[t];
  for (t in (p + 1) : N) {
    y_rest[t - p] = y[t]';
    for (i in 1 : p)
      lag_design[t - p, ((i - 1) * m + 1) : (i * m)] = y[t - i]';
  }
  identity = diag_matrix(rep_vector(1.0, m));
}
parameters {
  array[p] matrix[m, m] C; // The C_i
  cov_matrix[m] Sigma; // Error variance, Sigma
}
transformed parameters {
  array[p] matrix[m, m] phi; // The phi_i
  cov_matrix[p * m] Gamma; // (Stationary) variance of (y_1, ..., y_p)
  {
    array[p] matrix[m, m] V;
    array[2, p] matrix[m, m] phiGamma;
    for (i in 1 : p) {
      V[i] = C[i]' * C[i];
    }
    phiGamma = rev_mapping(C, V, Sigma);
    phi = phiGamma[1];
    for (i in 1 : p) {
      for (j in 1 : p) {
        if (i <= j) 
          Gamma[((i - 1) * m + 1) : (i * m), ((j - 1) * m + 1) : (j * m)] = phiGamma[2, j
                                                                    - i + 1];
        else 
          Gamma[((i - 1) * m + 1) : (i * m), ((j - 1) * m + 1) : (j * m)] = phiGamma[2, i
                                                                    - j + 1]';
      }
    }
  }
}
model {
  vector[p * m] mut_init = rep_vector(0.0, p * m); // Marginal mean of (y_1^T, ..., y_p^T)^T
  matrix[p * m, m] B = rep_matrix(0.0, p * m, m);
  matrix[m, m] L_Sigma = cholesky_decompose(Sigma);
  matrix[m, N - p] conditional_residuals;
  for (i in 1 : p)
    B[((i - 1) * m + 1) : (i * m),  : ] = phi[i]';
  conditional_residuals = (y_rest - lag_design * B)';
  // Likelihood:
  y1top ~ multi_normal_cholesky(mut_init, cholesky_decompose(Gamma));
  target += -0.5 * dot_self(to_vector(mdivide_left_tri_low(L_Sigma, conditional_residuals)))
            - (N - p) * sum(log(diagonal(L_Sigma)));
  // Prior:
  Sigma ~ inv_wishart(df, identity);
  for (i in 1 : p) 
    to_vector(C[i]) ~ std_normal();
}

generated quantities {
  matrix[m, m * p] topblock;
  for (lag in 1 : p) {
    topblock[ : , ((lag - 1) * m + 1) : (lag * m)] = phi[lag];
  }
  matrix[m * p, m * p] companion = rep_matrix(0.0, m * p, m * p);
  companion[1 : m,  : ] = topblock;
  if (p > 1) {
    companion[(m + 1) : (m * p), 1 : (m * (p - 1))] = diag_matrix(
                                                                  rep_vector
                                                                  (1.0,
                                                                   m
                                                                   * (
                                                                   p - 1)));
  }
  complex_vector[m * p] lambdas = eigenvalues(companion);
  vector[m * p] lambda_moduli = abs(lambdas);
  real max_lambda_modulus = max(lambda_moduli);
}


