# Hypnos external validation — propofol PK→BIS on VitalDB

> **NOT FOR CLINICAL USE.** Reproducible, envelope-stratified external validation of the propofol depth-of-anaesthesia stack against **measured BIS** on the open VitalDB cohort (v0.4 VE2/VE3). Computed by Hypnos's own reference kernels; kept strictly separate from any publisher-reported metric.

- **cohort:** vitaldb · 30 subjects scored (of 30 requested) · seed 7
- **target:** measured BIS · **shared PD link:** `pd_effect.propofol.bis_sigmoid`
- **tracks:** BIS/BIS, Orchestra/PPF20_RATE · PPF20 = 20.0 mg/mL

## Leaderboard (ranked by overall MDAPE — inaccuracy; lower is better)

| Rank | Model | Tier | MDPE % (bias) | MDAPE % | wobble % | in-env MDAPE % (n) | out-env MDAPE % (n) |
| ---: | --- | :---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `marsh_1991` | B | -32.3 | 32.3 | 7.3 | 32.3 (30) | — |
| 2 | `schnider_1998` | B | -33.4 | 33.7 | 7.6 | 33.1 (22) | 34.9 (8) |
| 3 | `eleveld_2018` | A | -34.4 | 34.4 | 7.1 | 34.4 (30) | — |

## Caveats (domain-review items — these numbers are evidence, not yet verified facts)

- Observation = MEASURED BIS (independent); the pump's predicted Ce is never used.
- Propofol-only PK->BIS stack: remifentanil synergy is NOT modelled, so BIS is systematically OVER-predicted (depth under-predicted) in balanced anaesthesia — a real, interpretable bias. The cross-model RANKING is still apples-to-apples (shared cohort + PD link).
- All PK models share one PD-BIS link, so differences reflect the PK model + that shared link.
- Track names + PPF20 vial conc + the shared PD link are domain-review items; verify per cohort.
- Raw VitalDB records are NOT committed (spec §3); only this aggregate leaderboard + manifest are.

## Reproduce

```bash
pip install vitaldb pandas
python scripts/fetch_vitaldb.py --use-cache --seed 7   # byte-identical from the local cache
```

*See the companion `vitaldb_bis_leaderboard.json` for the full machine-readable record + manifest (the pinned cohort, tracks, seed, and Hypnos version).*
