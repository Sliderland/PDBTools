# Paper and supplement audit (2026-09-23)

Sources: `HeapsStanPrograms/Heaps - 2023 - Enforcing Stationarity through the Prior in Vector Autoregressions.pdf` and `Heaps - Supplementary Information.pdf`, relative to the repository root. The latter is the 27-page author supplement from arXiv 2004.09455v3. Compare against the pre-refactor commits in Git history before changing any statistical model.

## Verdict on completed refactors

| File | Verdict | Paper anchor |
| --- | --- | --- |
| `statrmlprior_prior.stan` | Keep the SVD removal. With `C = U D W'`, the old `V = W D^2 W'` and `Q = W U'` imply `sqrt(V) Q = C'`. The `C` and `Sigma` priors are unchanged. | Supplement S3.3, Eq. S34 |
| `ansleykohn_prior.stan` | Keep. The Cholesky reverse mapping and stationary first-`p` density remain; the lag design has the same coefficient orientation. No explicit prior was added. | Supplement S1.3, especially Eq. S31 |
| `statpriorPDB_prior.stan` | Keep the refactor. The principal symmetric-root mapping, prior, stationary first-`p` density, and lag orientation are preserved. | Supplement S1.3 and S2 |
| `semiconj_prior.stan` | Keep the refactor. The conditional likelihood and baked priors match its pre-refactor file, including the zero-length off-diagonal case at `m=1`. | Main paper Section 4 (conditional likelihood) |
| `statinvertprior_prior.stan` | The new `N >= q+1` data bound excluded valid `N=q`; corrected to `N >= q` in commit `4013203`. The noncentered initial state and companion-shift changes otherwise preserve the pre-refactor distribution. | Supplement S7 |

The manually whitened likelihoods omit only Gaussian constants independent of all parameters. This changes some absolute `log_prob` values, not the posterior. The VARMA noncentering changes the sampler coordinates and therefore should be compared after transforming `z_init` back to `init`, rather than by equating log densities in different coordinates.

## Existing differences from the paper application

These differences **predate** the refactors. Changing them in place would change the posterior currently fitted by the baked variants.

1. `semiconj_prior.stan` uses diagonal `Normal(1, 10)` and off-diagonal `Normal(0, 1)`. The paper's application uses zero means and `W_{k3} = 10 I_{m^2}`, so every coefficient has `Normal(0, sqrt(10))`; see supplement S6.2.
2. `statpriorPDB_prior.stan` fixes the exchangeable hyperparameters to the paper's `m=3` values for all `m`: `f=sqrt(0.455)`, `g=1.365`, `h≈0.071`. The paper's `m=10,20` runs use `f=sqrt(0.700)`, `g=2.100`, `h=0.333`; see supplement S6.2. The generic `statpriorPDB.stan` accepts these as data.
3. `statinvertprior_prior.stan` uses 30 fixed Lyapunov-doubling steps to approximate its stationary initial covariance. This is already present in the original/baked baseline. For a scalar transition `r`, the omitted covariance fraction is `r^(2*2^30)`; at `r=0.9999999999` it is about 80.6%. The paper's stationary initialization in S7 is exact in principle. Do not replace this with an adaptive stop rule without checking HMC smoothness and computational cost.
4. The paper uses the Ansley–Kohn Cholesky parameterization for a **frequentist maximum-likelihood** comparison (main paper Section 5). Sampling `ansleykohn_prior.stan` with no explicit parameter priors instead defines a flat-prior Bayesian target; that file remains the same as its baseline, but HMC draws from it are not the paper's MLE results.

## Follow-up implementation prompts, only if paper reproduction is requested

### Semi-conjugate application prior

Create a separate paper-application baked-prior variant of `semiconj_prior.stan`, preserving its conditional VAR likelihood and generated quantities. Set every coefficient prior to `Normal(0, sqrt(10))` and retain `Sigma ~ IW(m+4, I)`, per supplement S6.2. Do not silently overwrite the existing baked variant. Compare log density and gradients with `semiconj.stan` supplied the same hyperparameters for `p=1` and `p=4`, including `m=1`.

### Exchangeable application prior for larger dimensions

Create a separate paper-application baked-prior variant of `statpriorPDB_prior.stan`. Use S6.2's `m=3` hyperparameters for `m=3`, and its `m=10,20` hyperparameters for those two dimensions. Define what to do for other dimensions explicitly. Preserve the symmetric-root map, hierarchical prior form, stationary initial density, conditional likelihood, and generated quantities. Compare same-point log density and gradients with generic `statpriorPDB.stan` supplied identical hyperparameters.

### VARMA stationary covariance accuracy

Evaluate the 30-step Lyapunov-doubling approximation against a direct Kronecker solve for small state dimensions over a grid of stable companion spectral radii, including values close to one. Record relative covariance error, minimum eigenvalue, autodiff gradient error, and runtime. Propose a differentiable Stan implementation that maintains the S7 stationary initial distribution without a prohibitively large matrix solve. Keep the current file unchanged until an exact or numerically controlled replacement is demonstrated.
