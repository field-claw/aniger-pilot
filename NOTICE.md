# NOTICE — unified *Aspergillus niger* pilot

This repository bundles **two models at different scales** for the same
organism, plus the annotation crosswalk between them. They answer different
questions and are deliberately **not** merged into one SBML file.

## Contents and licences

| Path | What | Scale | Licence |
|---|---|---|---|
| `models/genome_scale/iJB1325/iJB1325_ATCC1015.xml` | iJB1325, *A. niger* ATCC 1015 genome-scale model (BiGG) | 2320 rxn / 1818 met / **1325 gene** | **BiGG non-profit academic licence** — `THIRD_PARTY_LICENSES/BiGG_LICENSE.txt` |
| `models/curated_ccm/aniger_ccm_refined.xml` | curated central-carbon-metabolism model, built in-house | 28 rxn / 24 met / **17 gene** | **MIT** (in-house, no third-party GEM licence) |
| `src/`, `verify.py`, `crosswalk/*.py` | tooling | — | **MIT** (`LICENSE`) |

**The `LICENSE` file covers `src/`, `verify.py` and the curated CCM model only.
It does NOT cover iJB1325.** BiGG models are **not** CC BY: academic and
non-profit use is free, commercial use requires contacting `invent@ucsd.edu`,
and redistribution must ship the complete copyright notice together with the
three paragraphs that follow it. The verbatim licence text is in
`THIRD_PARTY_LICENSES/BiGG_LICENSE.txt` and **must not be removed**.

## Why two models instead of one merged file

Measured, not assumed (asserted by check G7):

- metabolite id: 24 vs 1818 — **0 collisions**
- reaction id: 28 vs 2320 — **0 collisions**
- gene id: 17 vs 1325 — **0 collisions**, and the CCM uses symbolic names
  (`acoA`, `citA`, `cexA`) where iJB1325 uses numeric ids

There is additionally no reaction-level merge that would be honest here: 13 of
the CCM's 28 reactions have no genome-scale counterpart (`EX_phos`, `PHOSt`,
`EX_cit`, `EX_glucon`, the lumped `EMP_pyr`/`EMP_pep`, `ATP_sink`, `BIOMASS`,
`DM_bio`, `BGC_NRPS`, `BGC_PKS`, `DM_sm1`, `DM_sm2`). Merging would fabricate
reactions that exist in neither model.

## Known limitations (disclosed, not hidden)

1. **The curated CCM is a focused core, not whole-cell.** 28 reactions. No
   genome-wide coverage.
2. **Lumped approximations.** Glycolysis (`EMP_pyr`/`EMP_pep`) and
   `GOX`/`ATP_sink` are simplifying lumps; mass and charge are not strictly
   conserved there, so a tool such as MEMOTE would flag `mass_balance` on those
   reactions. Deliberate simplification, not a defect of a published GEM.
3. **Units are not comparable across the two models.** iJB1325 grows at
   0.9399 h⁻¹; the CCM's pFBA value is ~18.95 in **lumped model units**. Never
   cross-compare them (asserted by check G2).
4. **The crosswalk is partial and not stereo-aware.** 15 of 28 CCM reactions
   map to iJB1325 (43 distinct targets); 13 do not. Mapping uses stoichiometric
   fingerprints plus name keywords, and genes are propagated through GPRs with
   **no sequence alignment**. Fingerprints cannot distinguish D- from
   L-glucose: one mapped target, `BOUNDARY_bDGLCe` (L-glucose), is a
   chemically wrong annotation and is recorded as a known defect in
   `verify.py` (check G6c asserts the shipped set is exactly that one known
   entry, so the defect is disclosed rather than shipped silently, and any
   change in either direction is caught). **Do not treat the crosswalk as a
   validated orthology map.**
5. **No experimental validation.** iJB1325's dimensions match BiGG's published
   values, which proves the *file* is iJB1325. It does not prove the model is
   experimentally validated, and it does not predict any real phenotype.

## Honest boundary of `verify.py`

`verify.py` reproduces specific numbers on one machine with one solver
(GLPK). Passing it means: the shipped files load, their dimensions and growth
values reproduce, carbon is genuinely bounded in both models, the phosphate
switch is present, and the crosswalk references ids that really exist in the
shipped models. It does **not** mean the biology is correct.

## Superseded: iMA871

An earlier `iMA871-aniger-pilot` (BioModels) was abandoned: `gene=0`, an
artificial biomass sink, and a carbon guardrail that is not verifiable. The
curated CCM supersedes it. `EX_glucon` in the CCM retains an iMA871-derived
gluconate exchange route.