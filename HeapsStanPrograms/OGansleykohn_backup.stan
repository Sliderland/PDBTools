functions {
  /* Function to transform A to P (inverse of part 2 of reparameterisation) */
  matrix AtoP(matrix A) {
    int m = rows(A);
    matrix[m, m] B = tcrossprod(A);
    for(i in 1:m) B[i, i] += 1.0;
    return mdivide_left_tri_low(cholesky_decompose(B), A);
  }
  /* Function to perform the reverse mapping from the Appendix. The details of
     how to perform Step 1 are in Section S1.3 of the Supplementary Materials.
     Returned: a (2 x p) array of (m x m) matrices; the (1, i)-th component
               of the array is phi_i and the (2, i)-th component of the array
               is Gamma_{i-1}*/
  matrix[,] rev_mapping(matrix[] P, matrix Sigma) {
    int p = size(P);
    int m = rows(Sigma);
    matrix[m, m] phi_for[p, p];   matrix[m, m] phi_rev[p, p];
    matrix[m, m] Sigma_for[p+1];  matrix[m, m] Sigma_rev[p+1];
    matrix[m, m] L_for;           matrix[m, m] L_rev;
    matrix[m, m] Ll_for[p+1];
    matrix[m, m] Gamma_trans[p+1];
    matrix[m, m] phiGamma[2, p];
    // Step 1:
    Sigma_for[p+1] = Sigma;
    Ll_for[p+1] = cholesky_decompose(Sigma);
    for(s in 1:p) {
      // In this block of code L_rev is B^{-1} and L_for is a working matrix
      L_for = - tcrossprod(P[p-s+1]);
      for(i in 1:m) L_for[i, i] += 1.0;
      L_rev = cholesky_decompose(L_for);
      Ll_for[p-s+1] = mdivide_right_tri_low(Ll_for[p-s+2], L_rev);
      Sigma_for[p-s+1] = tcrossprod(Ll_for[p-s+1]);
    }
    // Step 2:
    Sigma_rev[1] = Sigma_for[1];
    Gamma_trans[1] = Sigma_for[1];
    for(s in 0:(p-1)) {
      L_for = Ll_for[s+1];
      L_rev = cholesky_decompose(Sigma_rev[s+1]);
      phi_for[s+1, s+1] = mdivide_right_tri_low(L_for * P[s+1], L_rev);
      phi_rev[s+1, s+1] = mdivide_right_tri_low(L_rev * P[s+1]', L_for);
      Gamma_trans[s+2] = phi_for[s+1, s+1] * Sigma_rev[s+1];
      if(s>=1) {
        for(k in 1:s) {
          phi_for[s+1, k] = phi_for[s, k] - phi_for[s+1, s+1] * phi_rev[s, s-k+1];
          phi_rev[s+1, k] = phi_rev[s, k] - phi_rev[s+1, s+1] * phi_for[s, s-k+1];
        }
        for(k in 1:s) Gamma_trans[s+2] = Gamma_trans[s+2] + phi_for[s, k] * Gamma_trans[s+2-k];
      }
      Sigma_rev[s+2] = Sigma_rev[s+1] - quad_form_sym(Sigma_for[s+1], phi_rev[s+1, s+1]');
    }
    for(i in 1:p) phiGamma[1, i] = phi_for[p, i];
    for(i in 1:p) phiGamma[2, i] = Gamma_trans[i]';
    return phiGamma;
  }
}
data {
  int<lower=1> m; // Dimension of observation vector
  int<lower=1> p; // Order of VAR model
  int<lower=1> N; // Length of time series
  vector[m] y[N]; // Time series
}
transformed data {
  vector[p*m] y1top;                 // y_1, ..., y_p
  vector[m] mu = rep_vector(0.0, m); // (Zero)-mean of VAR process
  for(t in 1:p) y1top[((t-1)*m+1):(t*m)] = y[t];
}
parameters {
  matrix[m, m] A[p];   // The A_i
  cov_matrix[m] Sigma; // Error variance, Sigma
}
transformed parameters {
  matrix[m, m] phi[p];   // The phi_i
  cov_matrix[p*m] Gamma; // (Stationary) variance of (y_1, ..., y_p)
  {
    matrix[m, m] P[p];
    matrix[m, m] phiGamma[2, p];
    for(i in 1:p) P[i] = AtoP(A[i]);
    phiGamma = rev_mapping(P, Sigma);
    phi = phiGamma[1];
    for(i in 1:p) {
      for(j in 1:p) {
        if(i<=j) Gamma[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] = phiGamma[2, j-i+1];
        else Gamma[((i-1)*m+1):(i*m), ((j-1)*m+1):(j*m)] = phiGamma[2, i-j+1]';
      }
    }
  }
}
model {
  vector[p*m] mut_init;    // Marginal mean of (y_1^T, ..., y_p^T)^T
  vector[m] mut_rest[N-p]; // Conditional means of y_{p+1}, ..., y_{N}
  // Likelihood:
  for(t in 1:p) mut_init[((t-1)*m+1):(t*m)] = mu;
  for(t in (p+1):N) {
    mut_rest[t-p] = mu;
    for(i in 1:p) {
      mut_rest[t-p] += phi[i] * (y[t-i] - mu);
    }
  }
  y1top ~ multi_normal(mut_init, Gamma);
  y[(p+1):N] ~  multi_normal(mut_rest, Sigma);
}
