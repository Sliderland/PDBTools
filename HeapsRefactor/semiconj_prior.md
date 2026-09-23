# Prompt: `semiconj_prior.stan`

Modernize `HeapsStanPrograms/semiconj_prior.stan` conservatively. Preserve the zero-mean, **unrestricted** VAR(p), likelihood conditional on the first `p` observations, independent normal priors for each `phi` entry, and `Sigma ~ IW(m+4, I)`. Its fixed prior constants are already in `transformed data` (lines 9–35). Keep `data` as `m,p,N,y`; do not enforce stationarity.

## Hot spots and changes

1. The repeated `multi_normal` likelihood (lines 40–49) factors `Sigma` and builds an array of conditional means. Consider a single `cholesky_factor_cov[m] L_Sigma`, `Sigma = L_Sigma * L_Sigma'`, `L_Sigma ~ inv_wishart_cholesky(m+4, I)`, and vectorized `multi_normal_cholesky` for all `t>p`. This is the same prior on `Sigma` only if the Cholesky distribution's Jacobian is retained. If the installed Stan version lacks this built-in, leave the covariance parameter and prior intact and factor it once inside the model.
2. Precompute the lagged data arrangement in `transformed data` if it reduces repeated indexing; matrix multiplication for all conditional means is acceptable only if it reproduces `sum_i phi[i] * y[t-i]` with the same orientation. Keep the likelihood conditional: do not add a density for `y[1:p]`.
3. The prior loop (lines 51–60) is simple but its diagonal/off-diagonal scales differ. Preserve diagonal mean 1, standard deviation 10, off-diagonal mean 0, standard deviation 1 **for every lag**. Simplify the identity scale matrix, and remove unused intermediate vectors only if clarity improves.
4. Keep the companion eigenvalue calculation (lines 62–79) in `generated quantities`; no stationarity filter belongs in the target.

## Required equivalence evidence

At the same `phi`, `Sigma`, and `y`, compare conditional means, the likelihood, each prior contribution, and generated eigenvalue summary. If the covariance factor is sampled instead, compare after mapping `L_Sigma * L_Sigma'` to the old coordinates and accounting for the Jacobian. Show `p=1` and `p=4`; keep output names available to downstream code.
