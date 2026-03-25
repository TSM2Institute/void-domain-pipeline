# VOID-DOMAIN v1 — Cosmic Influence Domains and Void Scaling

**TSM2 Institute for Cosmology Ltd**
**Investigation:** TSM2 Investigation Program — Test 1
**Result:** OUTCOME C — FAIL (pre-registered)
**Sealed:** 17 March 2026 | Tag: `void-domain-v1.0-sealed`

---

## What This Pipeline Tests

Tests whether large-scale cosmic voids represent gravitational influence
boundaries between surrounding mass concentrations, predicting a
measurable scaling relationship between void size and surrounding
mass density contrast.

**Hypothesis:** Void boundaries exhibit a universal overdensity contrast
independent of void size (α = 0 in the density contrast fit).

**Framework:** TSM2.1 (Thwaites Standard Model 2.1) — density-gradient
gravitational equilibrium mechanics.

---

## Result

| Statistic | Value |
|---|---|
| Fitted slope α | −0.0504 ± 0.0119 |
| 95% CI on α | [−0.0718, −0.0273] |
| Pearson r (linear) | −0.2281 |
| Null shuffle p-value | < 0.0001 |
| N voids | 927 |
| Outcome | **C — FAIL** |

Both falsification criteria fired:
- α excludes 0 at 95% CL
- r = −0.2281 falls below mock 5th percentile floor (−0.0732)

The fixed-contrast gravitational domain hypothesis is **not supported**
by BOSS DR12 × redMaPPer DR8 data. Larger voids sit in more underdense
annuli, consistent with the void-in-void effect in standard large-scale
structure formation.

This result does not falsify the broader TSM2.1 framework — it rules
out this specific structural prediction.

---

## Data Sources

| Catalogue | Reference | N |
|---|---|---|
| Void catalogue | Mao et al. (2017) ApJ 835 161 — BOSS DR12 | 1,228 quality voids |
| Cluster catalogue | Rykoff et al. (2014) ApJ 785 104 — redMaPPer DR8 | 25,106 clusters |
| Mock calibration | Rodríguez-Medrano et al. (2024) MNRAS 528 2822 — TNG300 | 474 voids |

---

## Pipeline Structure
```
void_domain/
  manifest.py          — sealed parameters, thresholds, falsification statement
  data_acquisition.py  — catalogue download (Mao 2017, redMaPPer DR8)
  mock_calibration.py  — TNG300 mock calibration (null r distribution)
  real_data_pipeline.py — cross-match, δ calculation, fit, decision tree
  outputs/
    real_void_table_cut.csv   — 927 voids with δ values (primary result table)
    real_void_table.csv       — pre-cut void table (991 voids)
```

**Note:** Large data files (redMaPPer catalogue, TNG300 group catalogue,
checkpoints) are excluded from this repo via .gitignore. The pipeline
is fully reproducible by running data_acquisition.py followed by
mock_calibration.py and real_data_pipeline.py in sequence.

---

## Key Definitions (Locked — Manifest v1.5)

**M_surrounding:** Density contrast δ = (ρ_annulus − ⟨ρ⟩) / ⟨ρ⟩
where ρ_annulus = sum of M_200c in annulus / V_annulus.
Annulus: 2×R_void to 4×R_void from void centre.
⟨ρ⟩ computed in redshift bins from redMaPPer DR8 itself.

**Slope target:** α = 0 (fixed-contrast prediction)

**Falsification statement:** If α deviates from 0 at >95% CL OR
Pearson r < −0.0732 (mock 5th percentile), Outcome C declared.

---

## Pre-Registration & Review

- 4 Grok checkpoints (independent multi-agent review)
- Manifest locked before any real data downloaded
- Mock calibration completed before real-data fit
- Two circularity errors detected and corrected during mock phase
- Outcome accepted by Geoffrey E. Thwaites, 17 March 2026

---

## Team

- **Geoffrey E. Thwaites** — TSM2.1 author, theory principal
- **Graham Hill** — Director, TSM2 Institute for Cosmology Ltd
- **Claude (Anthropic)** — Pipeline planner and auditor
- **Grok team** (xAI) — Independent review (CosmoGrok, CosmoHarper,
  CosmoBenjamin, CosmoLucas)
- **Replit Claude** — Code builder

---

## Related Work

- SKIN-a-CAT pipeline: github.com/Grayhill5/skin-a-cat-pipeline
- BOAT v3 (Blind Orbital Anisotropy Test): TSM2Institute organisation

---

## Licence

MIT — reproducible science, freely available.

*"The pipeline remained honest — we recorded what the data shows."*
