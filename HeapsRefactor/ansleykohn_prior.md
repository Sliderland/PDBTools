# Prompt: `ansleykohn_prior.stan`

Refactor `HeapsStanPrograms/ansleykohn_prior.stan` for lower HMC cost and clearer numerical behavior. Keep the existing posterior **exactly**: zero-mean stationary VAR(p), Ansley–Kohn Cholesky partial-autocorrelation map, stationary joint density for the first `p` observations, conditional Gaussian densities afterward, and no explicit prior on `A` or `Sigma`. Preserve the improper flat prior implied by the existing model. Do not import the exchangeable normal/gamma prior from `statpriorPDB_prior.stan` or change the transformation to principal symmetric roots.

## Hot spots

1. `AtoP` and `rev_mapping` (lines 5–69): each leapfrog step repeatedly performs Cholesky factorizations and triangular solves. The forward/reverse recursion and covariance blocks become ill conditioned near the unit-root boundary. Preserve the order, orientation, and indexing of these recursions.
2. `transformed parameters` (lines 88–108) builds dense `Gamma` of dimension `p*m`; `model` (lines 110–124) factors it again in `multi_normal`. This is expensive at `m=20`, `p=4`.
3. Generated companion eigenvalues (lines 125–143) are output-only; retain them there.

## Implementation plan

- First simplify allocations and repeated work without changing the map: cache Cholesky factors already computed within one `rev_mapping` call, avoid constructing matrices used once where a triangular solve directly computes the needed result, and remove redundant zero-mean arithmetic. Prove the matrix orientation for every triangular replacement; `mdivide_right_tri_low(B,L)` means `B * L^{-1}`.
- Form `L_Gamma = cholesky_decompose(Gamma)` once and use `multi_normal_cholesky` for the initial observation vector. Likewise factor `Sigma` once for the repeated conditional likelihood. Keep the same observation order and full stationary initial density. If reusing the factor of `Sigma` from the recursion, preserve the precise Cholesky convention.
- If moving `Sigma` to a `cholesky_factor_cov` parameter is explored, derive the density in the new coordinate system, including its Jacobian. Since the current `Sigma` has no explicit prior, this is easy to get wrong. Prefer leaving `Sigma` as `cov_matrix` in the first pass.
- Add a clear data-domain guard for `N > p` if needed by the chosen Stan version; the existing likelihood assumes at least one conditional observation.

## Required equivalence evidence

Compare `P`, `phi`, all lag covariance blocks of `Gamma`, initial/conditional log likelihood, and `max_lambda_modulus` at identical valid `A`, `Sigma`, and `y` for `p=1` and `p=4`. Any HMC improvement claim should report divergences, treedepth, ESS, and time per effective draw on the same data and sampler settings. If the flat prior causes persistent geometry problems, report that separately; do not repair it by changing the prior.
