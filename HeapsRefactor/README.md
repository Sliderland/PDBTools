# Heaps Stan model refactor briefs

These are implementation prompts for the five **baked prior** variants in `../HeapsStanPrograms/`. No Stan source is changed by this review. Apply each brief to its named `_prior.stan` file, preferably on an isolated branch or copy, and keep the existing model definition.

| File | Prior constants already in `transformed data`? | Main issue | Priority |
| --- | --- | --- | --- |
| `ansleykohn_prior.stan` | There are no numeric prior hyperparameters. No explicit prior is placed on `A` or `Sigma`. | Expensive recursion and stationary initial covariance; an improper flat prior is part of the present target. | Medium |
| `semiconj_prior.stan` | Yes: lagwise normal means/scales and `IW(m+4, I)`. | Dense covariance likelihood; this is the simplest reference model. | Low |
| `statinvertprior_prior.stan` | Yes: normal/gamma hyperpriors and `IW(m+4, I)`. The hierarchy still contains sampled `Amu`, `Aomega`, `Dmu`, `Domega`. | A 30-step Lyapunov doubling recursion, repeated symmetric roots, and a coupled initial state. | Very high |
| `statpriorPDB_prior.stan` | Yes: normal/gamma hyperpriors and `IW(m+4, I)`. The hierarchy still contains sampled `Amu`, `Aomega`. | Repeated eigendecompositions in the partial-autocorrelation recursion. | Medium |
| `statrmlprior_prior.stan` | Yes: `IW(m+4, I)`; `C` has a standard-normal prior. | The SVD is avoidable algebraically; the initial-observation variable has inconsistent names. | High |

All five `data` blocks contain only dimensions (`m`, `p`, optionally `q`, `N`) and the observed series `y`. Here “baked in” refers to **fixed hyperparameters**, not the elimination of sampled hierarchical parameters. The Ansley–Kohn file deserves special attention: its absence of explicit priors is intentional in the existing code, but it is not a proper prior. Adding priors would change the posterior.

## Shared preservation contract for the implementing agents

- Preserve the zero mean; the exact conditional versus stationary initial likelihood; all lag and MA sign conventions; all prior families and numeric values; and the generated quantity names and meanings. Keep the models' data interfaces limited to observations and dimensions.
- Preserve the original mathematical transformation. In particular, a Cholesky factor is **not** interchangeable with a principal symmetric matrix square root in `statpriorPDB_prior.stan` or `statinvertprior_prior.stan`. It changes the map from unconstrained matrices to VAR coefficients.
- Prefer algebraic identities, shared decompositions, and built-in Cholesky densities. If changing the parameterization of `Sigma`, use `inv_wishart_cholesky` with the same degrees of freedom and identity scale factor; it includes the requisite covariance-to-factor Jacobian. Do not replace the inverse-Wishart prior with an LKJ or independent-scale prior.
- Treat optimization as a sequence of small changes. For each change, compare the old and new deterministic outputs (`phi`, `theta` where present, stationary covariance, and generated quantities) and log density at identical valid parameter/data values. For a changed coordinate system, compare after mapping into common coordinates and account for its Jacobian. Report compiler version, numerical error, and any HMC diagnostics from a short representative fit if fitting is in scope. Do not label approximate agreement as exact equivalence.
- Do not move eigenvalue calculations from `generated quantities` into the HMC target. They currently report stability and do not affect gradients.

## Sources and built-ins

- Local paper: `../HeapsStanPrograms/Heaps - 2023 - Enforcing Stationarity through the Prior in Vector Autoregressions.pdf`; accessible [arXiv HTML](https://arxiv.org/html/2004.09455), especially Sections 2, 3.2, 3.5, and 6.
- [Stan matrix operations](https://mc-stan.org/docs/functions-reference/matrix_operations.html): `eigendecompose_sym`, `cholesky_decompose`, and specialized triangular/SPD solves. `eigendecompose_sym` shares work otherwise repeated by separate eigenvalue and eigenvector calls; it requires Stan 2.33 or later.
- [Stan multivariate normal Cholesky density](https://mc-stan.org/docs/functions-reference/distributions_over_unbounded_vectors.html): `multi_normal_cholesky` is equivalent to `multi_normal` given the same covariance factor.
- [Stan inverse-Wishart Cholesky density](https://mc-stan.org/docs/functions-reference/covariance_matrix_distributions.html): `inv_wishart_cholesky` is available from Stan 2.30 and includes the factorization Jacobian.

There is no simple built-in Stan function that replaces the complete VAR partial-autocorrelation recursion or the discrete Lyapunov covariance calculation. The latter's current fixed 30-step calculation is an approximation; replacing it with an exact solver would change the implemented target slightly.
