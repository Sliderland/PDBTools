data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=p> t;
  matrix[t, N] y;
}
generated quantities {
  // N: number of variables
  // t: number of rows/time points
  // p: number of lags
  
  matrix[N, N * p] topblock;
  matrix[N * p, N * p] companion;
  if (p == 1) {
    companion = topblock;
  } else {
    companion[1 : (N), 1 : (N * p)] = topblock;
    companion[(N + 1) : (N * p), 1 : (N * (p - 1))] = diag_matrix(
                                                                  rep_vector
                                                                  (1.0,
                                                                   N
                                                                   * (
                                                                   p - 1)));
  }
}
