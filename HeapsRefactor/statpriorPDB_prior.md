# Prompt: `statpriorPDB_prior.stan`

Refactor `HeapsStanPrograms/statpriorPDB_prior.stan` while preserving the paper's exchangeable stationary VAR(p) model. Keep the zero mean, principal **symmetric** square roots, the map `A -> P -> phi`, the stationary density of the first `p` observations, `Sigma ~ IW(m+4,I)`, and hierarchical normal/gamma priors on diagonal and off-diagonal `A` entries. Fixed `es=0`, `fs=sqrt(0.455)`, `gs=1.365`, `hs=0.071175` for both groups already live in `transformed data` (lines 100–129). `Amu` and `Aomega` remain sampled parameters; do not bake them into constants.

## Hot spots

1. `sqrtm` (lines 5–18) calls `eigenvalues_sym` and `eigenvectors_sym` separately, then builds the principal root. It runs repeatedly inside `AtoP` and both passes of `rev_mapping` (lines 20–91). Repeated or close eigenvalues and covariances near singularity can make gradients and SPD solves difficult.
2. `Gamma` construction (lines 138–159) plus the stationary initial likelihood (lines 160–173) creates and factors a dense `p*m` covariance for each gradient evaluation.
3. The hierarchical precision parameters (lines 130–137, 176–189) can have strong posterior dependence with the `A` entries. Any reparameterization must preserve the gamma **rate** convention and normal standard deviation `1/sqrt(Aomega)`.

## Implementation plan

- Replace each separate symmetric eigenvalue/eigenvector pair by one `eigendecompose_sym` call if Stan >=2.33. Keep the principal square-root result and its symmetry. Do not substitute `cholesky_decompose` for `sqrtm`: the paper's exchangeable map depends on symmetric roots. Do not clamp eigenvalues without proving that it changes only roundoff and recording that choice.
- Cache repeated roots or solves in `rev_mapping`; use Stan's specialized `mdivide_left_spd`, `mdivide_right_spd`, and triangular solves only where the dividend really has the declared property. Build `Gamma` with the same block orientation and use `multi_normal_cholesky` after one Cholesky factorization for the initial vector and conditional series.
- Consider noncentering `A` as `A_ij = group_mean + group_scale * z_ij`, `z_ij ~ N(0,1)` to weaken the hierarchy when data are weak. This is a change of coordinates, so derive the induced target correctly and retain `Amu`/`Aomega` with their original priors; check whether it actually helps before keeping it. A Cholesky covariance parameter with `inv_wishart_cholesky` is another exact option.
- Require `N>p` for valid input if the loop bounds need an explicit guard. Keep `max_lambda_modulus` in generated quantities.

## Required equivalence evidence

Compare `sqrtm`, `P`, every `phi`, `Gamma`, initial and conditional log densities, hyperprior and conditional-prior contributions, and generated quantities at identical mapped values for `p=1` and `p=4`. Report any eigenvalue or SPD failures rather than silently shifting the prior or stationarity boundary.
