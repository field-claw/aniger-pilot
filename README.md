# aniger-pilot — *Aspergillus niger* virtual microbe, unified pilot

Two *A. niger* models at different scales, plus the annotation crosswalk
between them, in one clone-and-run repo. One `verify.py` checks all three
layers (**12 hard checks**).

```bash
pip install -r requirements.txt
python verify.py   # exit 0 = all 12 checks passed
```

## Why one repo with two models

Both models are *Aspergillus niger*, and they answer different questions:

| Model | Scale | What it is | Growth |
|---|---|---|---|
| **iJB1325** (BiGG) | 2320 rxn / 1818 met / **1325 gene** | published genome-scale model, *A. niger* ATCC 1015 | **0.9399 h⁻¹** |
| **aniger_ccm** (in-house) | 28 rxn / 24 met / **17 gene** | curated central-carbon-metabolism core | **18.95 model units** (not h⁻¹) |
| **crosswalk** | — | reaction-level annotation mapping between them | — |

They are **not** merged into one SBML file, and that is a measured decision
(check G7 asserts it):

- metabolite ids: 24 vs 1818 → **0 collisions**
- reaction ids: 28 vs 2320 → **0 collisions**
- gene ids: 17 vs 1325 → **0 collisions** — the CCM uses symbolic names
  (`acoA`, `citA`, `cexA`), iJB1325 uses numeric ids

More decisively: **13 of the CCM's 28 reactions have no genome-scale
counterpart** (`EX_phos`, `PHOSt`, `EX_cit`, `EX_glucon`, the lumped
`EMP_pyr`/`EMP_pep`, `ATP_sink`, `BIOMASS`, `DM_bio`, `BGC_NRPS`, `BGC_PKS`,
`DM_sm1`, `DM_sm2`). Merging would fabricate reactions that exist in neither
model. Side by side, each model stays internally consistent and the mapping
between them stays an explicit, inspectable artifact.

## What the models actually show

**iJB1325** — genome-scale, *A. niger* ATCC 1015. Reproduces BiGG's published
dimensions exactly and grows at 0.9399 h⁻¹. Closing all 353 organic-carbon
exchange boundaries collapses growth to 0 (100%), which proves flux accounting
is carbon-bounded; it does not predict any phenotype.

**aniger_ccm** — a 28-reaction core with a *real* biomass reaction
(`BIOMASS`, SBO:0000629, real precursors + GAM), genuine GPR, and a
verifiable carbon guardrail (5 organic-carbon exchanges → 100% collapse).

Its headline result is the **phosphate switch**, which is the citric-acid
production phenotype industrial strains are selected for:

| Phase | Growth | Citrate secretion capacity |
|---|---|---|
| phosphate sufficient | 18.95 | 6.00 |
| phosphate depleted (`EX_phos` closed) | 0.00 | **12.00** (2×) |

**One measured trap.** You cannot read citrate overflow off a
biomass-maximising FBA solution: with biomass as the sole objective the LP
returns a carbon-minimal knife-edge solution and `EX_cit` is exactly `0.0000`
in *both* phases — a false negative that hides the phenotype completely.
Overflow is an *alternative* objective, so each phase is solved as "keep growth
≥ 50% of the phase-1 maximum, then maximise `EX_cit`".

## Layout

```
models/genome_scale/iJB1325/iJB1325_ATCC1015.xml   BiGG full GEM (BiGG licence)
models/curated_ccm/aniger_ccm_refined.xml           in-house CCM (MIT)
src/aniger_fba.py                                   zero-dependency wrapper
crosswalk/aniger_annotation_crosswalk.{py,json}     mapping + its measured output
verify.py                                           12 hard checks, fail-closed
THIRD_PARTY_LICENSES/BiGG_LICENSE.txt               verbatim, must not be removed
```

## Two rules this repo enforces in code

Both were learned by getting them wrong, not by reading about them.

**1. Identity assertions never ride on the numeric checks.** iJB1325 was
originally published as `ijb1325-ecoli-pilot`, labelled *E. coli*. Every
dimension and growth number matched published values and every check was
green — while the organism label was wrong. It is *A. niger* ATCC 1015. Check
G1 now asserts the organism of **both** models separately from any
measurement, and G2 asserts the two growth **units** are not conflated.

**2. Overflow products are invisible under the natural objective.** See the
measured trap above. Any downstream reuse of `phosphate_switch()` must keep
the growth-floor structure.

## Honest limits

- The CCM is a focused core, not whole-cell. Its lumped reactions (`EMP_*`,
  `GOX`, `ATP_sink`) do not strictly conserve mass/charge — a tool like MEMOTE
  would flag `mass_balance` on them. Deliberate simplification.
- **Units are not comparable.** 0.9399 h⁻¹ vs 18.95 model units.
- **The crosswalk is partial and not validated as orthology.** 15/28 reactions
  map (43 distinct iJB1325 targets), genes are propagated through GPRs with no
  sequence alignment, and one mapped target (`BOUNDARY_bDGLCe`, L-glucose) is
  a stereo mismatch recorded as a known defect (check G6c).
- No experimental validation. Matching BiGG's published dimensions proves the
  *file* is iJB1325, nothing more.

## Licence

`LICENSE` (MIT) covers `src/`, `verify.py` and the curated CCM model.
**iJB1325 is NOT covered** — BiGG models are not CC BY. Academic/non-profit
use is free; commercial use requires contacting `invent@ucsd.edu`.
`THIRD_PARTY_LICENSES/BiGG_LICENSE.txt` ships verbatim and must not be
removed. See `NOTICE.md`.

## Superseded

`iMA871` (BioModels) is abandoned: `gene=0`, an artificial biomass sink, and a
carbon guardrail that is not verifiable. The curated CCM supersedes it.
`EX_glucon` in the CCM retains an iMA871-derived gluconate exchange route.

## Series

`field-claw/vitamind-virtual-microbe` is the series hub.