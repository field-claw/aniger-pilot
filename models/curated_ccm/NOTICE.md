# NOTICE — *A. niger* CCM (curated) pilot

## Model provenance
`data/aniger_ccm_refined.xml` is a **curated central-carbon-metabolism (CCM)
model built in-house** as part of the VitaMind virtual-microbe project
(Phase 1 refined). It is **not** a published full-genome GEM (unlike iJB1325
from BiGG or ecYeastGEM). Source builder:
`vitamind_core/assets/gem/aniger_ccm_refined_build.py`.

## License
The model and all `src/` tooling are released under **MIT** (see `LICENSE`).
No third-party GEM license applies.

## Known limitations (disclosed, not hidden)
1. **Curated core, not whole-cell.** 28 rxn / 24 met / 17 gene — a focused CCM
   core. No full genome coverage, no full metabolite pool.
2. **Lumped approximations.** Glycolysis (`EMP_pyr`/`EMP_pep`) and
   `GOX`/`ATP_sink` are simplifying lumps; mass/charge are not strictly
   conserved there. MEMOTE would flag `mass_balance` fails on those reactions.
   This is a deliberate Phase-1 simplification, not an error in a published GEM.
3. **Biomass is in model units, not h⁻¹.** The pFBA value ~18.95 reflects
   lumped biomass coefficients, NOT a calibrated per-hour growth rate. Do NOT
   compare to iJB1325 (0.9399 h⁻¹) or ecYeastGEM (0.087974 h⁻¹).
4. **What IS solid:** a *real* biomass reaction (`BIOMASS`, SBO:0000629) with
   real precursors + GAM; a *verifiable* carbon guardrail (closing 5
   organic-carbon exchanges collapses growth to 0, drop ~100%); genuine GPR
   rules (17 genes); and the native *A. niger* phosphate-switch phenotype,
   asserted by hard check L7 — closing `EX_phos` collapses growth 18.95 -> 0
   while citrate secretion capacity doubles 6.00 -> 12.00.
5. **Citrate overflow is invisible under a biomass objective.** Solved with
   biomass as the sole objective, `EX_cit` is exactly 0.0000 in both phosphate
   phases (carbon-minimal knife-edge solution). The phenotype only appears when
   citrate secretion is the objective subject to a growth floor. Any downstream
   reuse of `phosphate_switch()` must keep that structure.

## Why this replaces iMA871
`iMA871` (BioModels) ships `gene=0`, an artificial biomass sink, and a carbon
guardrail that is **NOT** verifiable — it cannot serve as a trustworthy
reference. This curated CCM model is the honest, working *A. niger* pilot for
the series. The `iMA871-aniger-pilot` repository was created empty by mistake
and is superseded by this one.
