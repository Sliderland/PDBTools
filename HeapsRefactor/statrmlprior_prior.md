# Prompt: `statrmlprior_prior.stan`

Refactor `HeapsStanPrograms/statrmlprior_prior.stan` for the Roy et al. stationary VAR(p) parameterization. Preserve `C[i] ~ std_normal()`, `Sigma ~ IW(m+4,I)`, zero mean, the Roy reverse map, stationary density for the first `p` observations, and the generated stability measure. The inverse-Wishart hyperparameters are already baked in (lines 98–109).

## Highest value change: remove an avoidable SVD

`CtoVQ` (lines 19–28) computes an SVD for every `C[i]`. With `C=U diag(s) V'`, its outputs are `V_i=V diag(s^2)V'=C'C` and `Q_i=VU'`. In `rev_mapping`, `Q_i` occurs only in `sqrtm(V_i) * Q_i * sqrtm(D_i)` (lines 49, 66). The first two factors multiply to **`C'`**, so the term is exactly `C' * sqrtm(D_i)` for full-rank `C`, with a continuous extension at rank deficiency. Pass `C` directly to `rev_mapping`, compute `V_i = C_i' * C_i`, and replace that product with `C_i'`. This removes SVD singular-vector sign/degeneracy behavior and one symmetric root per lag without changing the induced `phi` or prior on `C`. Check the transpose orientation carefully.

## Other issues

- Lines 99–108 declare `y_1top` but fill and later use `y1top`. Consolidate to one variable, declare it before executable statements for compatibility, and preserve the first-`p` observation order. The unused `y_1top` is misleading and the declaration placement may fail on older parsers.
- `rev_mapping` (lines 33–90) allocates several maximum-size arrays and repeatedly solves growing SPD systems. Allocate only needed blocks if practical, cache reusable factorizations, and preserve every lag orientation. The remaining `sqrtm(D_i)` is the principal symmetric root; a Cholesky factor is not an algebraically identical substitution here.
- Build and factor stationary `Gamma` once; use `multi_normal_cholesky` for its initial density and the conditional densities. A Cholesky covariance parameter plus `inv_wishart_cholesky` may also help while preserving the `IW(m+4,I)` distribution.
- Do not change the standard-normal `C` prior to a prior on `V` or `Q`; this would change the stationary VAR prior. Keep eigenvalue diagnostics in generated quantities.

## Required equivalence evidence

For several full-rank `C` matrices, including `p=1` and `p=4`, compare old `sqrtm(V)*Q` with `C'`, then compare `phi`, all `Gamma` blocks, and full log density after the substitution. Include a near rank-deficient case to assess numerical behavior, report any difference, and verify all output names are preserved.
