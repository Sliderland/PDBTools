# Prompt: `statinvertprior_prior.stan`

Refactor `HeapsStanPrograms/statinvertprior_prior.stan` cautiously. It implements a zero-mean stationary and invertible VARMA(p,q), a joint stationary Gaussian density for latent pre-sample `y` and innovations, then Gaussian observation densities. Keep the AR and MA partial-autocorrelation maps, the minus sign in `theta` (lines 204–207), the exact latent-state indexing, `Sigma ~ IW(m+4,I)`, and the separate hierarchical normal/gamma priors for `A` and `D`. Fixed hyperparameters already live in `transformed data` (lines 146–176); `Amu/Aomega/Dmu/Domega` are sampled hierarchy parameters and must remain so.

## Main computational and geometric problems

1. `initial_joint_var` (lines 79–137) propagates a `(p+q)m` state for **30 matrix-doubling steps** on every gradient evaluation. It sums a finite covariance series with `2^30` terms, then symmetrizes. Near a unit root, truncation error can be material; replacing it with a direct Lyapunov solution is a model change relative to this implemented target. Do not silently change the step count, use a parameter-dependent stopping rule, or drop the initial density.
2. Symmetric matrix roots occur throughout `AtoP` and `rev_mapping` (lines 5–75), for both AR and MA maps. The separate eigenvalue/eigenvector calls duplicate work and can behave poorly near repeated eigenvalues or the stationarity boundary.
3. `init` (lines 177–181) has a covariance `Omega` determined by `phi`, `theta`, and `Sigma`, producing a strongly coupled high-dimensional latent state. The conditional innovations are recursively recovered in `mut` (lines 218–257), so a wrong index or MA sign changes the likelihood.

## Implementation plan

- Make one symmetric eigendecomposition per `sqrtm` via `eigendecompose_sym` on a modern Stan version. Preserve the **principal symmetric** roots and all right/left SPD solves; do not use Cholesky factors as substitutes for these roots. Cache roots used more than once in the reverse mapping.
- Keep the fixed 30-step covariance approximation while optimizing its implementation. Reuse constant state-layout work where possible, avoid redundant matrix construction, and retain the final symmetrization. If evaluating an exact discrete Lyapunov solver as a separate research option, report its discrepancy from the present 30-step `Omega` across posterior-relevant states before proposing adoption. Stan has no direct built-in Lyapunov solver.
- Explore noncentering the latent state: sample `z_init ~ std_normal()` and set `init = mut_init + cholesky_decompose(Omega) * z_init`. This reproduces the same joint model when mapped back to `init`; remove the old `init ~ multi_normal(...)` factor and account for the coordinate transformation. If full noncentering worsens geometry, consider a partial parameterization, still with a proof of density equivalence.
- Use `multi_normal_cholesky` with a shared factor of `Sigma` for the observation density, and for the initial density if it remains centered. A Cholesky covariance parameter plus `inv_wishart_cholesky` can preserve the exact `Sigma` prior.
- Preserve gamma rate parameters, normal scales `1/sqrt(precision)`, and the two separate AR/MA hierarchies. Add explicit valid-data checks for `N>q` and `p,q>=1` if needed for the existing loops. Keep the eigenvalue report in generated quantities.

## Required equivalence evidence

At identical mapped `A,D,Sigma,init` and `y`, compare `P`, `R`, `phi`, `theta`, every block of `Omega`, recovered conditional innovations, conditional means, and the full target density. Include `p=q=1` and the paper's `p=4,q=2`; examine a near-boundary draw for sensitivity of the 30-step covariance. Report numerical error and HMC diagnostics separately. Do not describe a change in `Omega` or the initial-state law as a pure optimization.
