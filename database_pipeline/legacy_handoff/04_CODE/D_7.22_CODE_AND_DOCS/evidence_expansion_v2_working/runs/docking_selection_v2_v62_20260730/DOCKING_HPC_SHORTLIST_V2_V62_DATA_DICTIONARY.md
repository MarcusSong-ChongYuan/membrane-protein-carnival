# Docking HPC shortlist V2 (MemPro V6.2)

The source release is immutable. Positive-negative conflicts are excluded from
default HPC lists and retained in a separate review table.

## Workload choices

| Set | Rule | Rows |
|---|---|---:|
| Pilot | Per target: R0 1, D1 5, D2 3, D3 2 | 6,019 |
| Standard | Per target: R0 3, D1 20, D2 10, D3 5 | 14,713 |
| Full non-conflict | All R0-D3 pairs without V6.2 positive-negative conflict | 186,625 |

R0 is a redocking/protocol-validation set. D1-D3 are prospective or mechanistic
sets. A candidate receptor PDB is provided for triage only; receptor state,
chain, biological assembly, microstates, and the docking box still require
preparation and validation.

## Ranking

Ranking is deterministic within target and tier. It favors independent sources,
BE1/BE2 support, stronger quantitative potency, approved/clinical/probe/
endogenous status, and tractable physicochemical properties. Ranking is not a
binding-affinity prediction.
