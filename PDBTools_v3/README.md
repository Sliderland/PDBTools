# PDBTools V3: unconstrained parameter counts

This folder contains an isolated V3 copy of the PDBTools builder and its
runnable V2 workflows. V3 stores one integer unconstrained-parameter count per
included model parameter in posterior `dimensions`, following the exploratory
contract implemented in the current `Forked_posteriordb-r` branch. V2 and its
workflow files remain unchanged.

`PDBEntryBuilder_v3.R` infers those counts when posterior dimensions are
omitted, checks the sampled fit's counts against the registered posterior
before accepting reference draws, and continues to use dimension names to
select variables for diagnostics and writing. RStan inference uses a zero-chain
fit; CmdStanR uses a short fit.

The copied workflow scripts resolve the builder and repository-local
`HeapsStanPrograms` paths from their own location. Their existing database,
Stan source, and environment-variable configuration still applies. Some
examples depend on objects supplied by the calling session, as they did in the
original workflow files.

Run the focused V3 count checks from the repository root with:

```sh
Rscript PDBTools_v3/tests/test_unconstrained_parameter_counts.R
```
