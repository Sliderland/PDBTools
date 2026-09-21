# Stan sources for `VECM_info.RDS`

These four files are byte-for-byte copies of the Stan sources at Git commit
`5acaec9` (2026-06-29), the latest checked-in versions before
`../VECM_info.RDS` was saved on 2026-06-30. The RDS does not record a Git
commit or the Stan source text, so an uncommitted edit at fitting time cannot
be ruled out. The mapping below is verified against `../unconstrainedVECM.R`.

| Saved fit (`readRDS("../VECM_info.RDS")$models`) | Archived Stan source | Original source |
| --- | --- | --- |
| `vecm_baseline` | `vecm_baseline.stan` | `VECM.stan` |
| `vecm_urca_hmc` | `vecm_urca_hmc.stan` | `VECM_URCA_test.stan` |
| `vecm_priors` | `vecm_priors.stan` | `VECMwPriors.stan` |
| `vecm_long_run` | `vecm_long_run.stan` | `VECMLongRun.stan` |

All fits used `p = 4`, cointegration rank `h = 2`, and four chains. The baseline,
priors, and long-run fits used the same 1,000-row synthetic series generated
from the five-variable `urca::denmark` data by `fit_simulate_VECM()` in
`../helper.R`. The URCA-style fit used the original 55-row Denmark series plus
three centered quarterly seasonal dummies. The synthetic series itself and its
random seed were not saved separately, so rerunning the fitting script will
not exactly reproduce those input observations.

For goodness-of-fit and specification work on the **stored** posterior fits,
load `../VECM_info.RDS` and select `models$vecm_baseline`,
`models$vecm_urca_hmc`, `models$vecm_priors`, or `models$vecm_long_run`.
`samps` contains their saved posterior draws. Loading this RDS requires
substantial RAM (the compressed file is about 1.7 GiB).

The current `../VECM.stan` and `../VECMLongRun.stan` differ substantively from
the archived versions in *generated quantities* (the VAR-in-levels/companion
matrix calculations were revised after the RDS was saved). Their likelihood
and prior blocks were not changed in those later revisions. The current
`../VECM_URCA_test.stan` and `../VECMwPriors.stan` have the same contents as
their archived versions, apart from file mode.

No source files outside this folder were modified for this archive.
