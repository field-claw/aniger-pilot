"""Build an annotation crosswalk between the two in-house/published A. niger
models WITHOUT merging the SBML files (merging would entangle a BiGG-licensed
model with an MIT-licensed one; a crosswalk keeps both license-clean).

Direction: refined CCM (28 rxn, gene symbols) -> iJB1325 (2320 rxn, numeric
gene IDs). Mapping evidence, in priority order:
  E1 stoichiometric fingerprint : multiset of (formula, coefficient) equal
  E2 name keywords              : curated keyword hit in iJB1325 rxn name
Genes propagate: refined gene --(refined GPR)--> refined rxn --(map)-->
iJB1325 rxn --(iJB1325 GPR)--> candidate iJB1325 gene IDs.

Honest boundaries: lumped reactions (EMP_pyr/EMP_pep/BIOMASS/ATP_sink) map to
MULTIPLE iJB1325 reactions or none; BGC clusters have no iJB1325 counterpart.
Every entry records its evidence; unmapped entries stay unmapped.
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WS, "aniger_ccm_pilot", "src"))

import aniger_fba as af  # noqa: E402  (zero-dep wrapper from the pilot repo)

from cobra.io import read_sbml_model  # noqa: E402

REFINED_XML = os.path.join(WS, "aniger_ccm_pilot", "data",
                           "aniger_ccm_refined.xml")
IJB_XML = os.path.join(WS, "ecoli_pilot", "data", "iJB1325_ATCC1015.xml")
OUT_JSON = os.path.join(WS, "demo_community_out",
                        "aniger_annotation_crosswalk.json")

# curated keyword hints: refined rxn id -> keywords that must appear in the
# iJB1325 reaction name (case-insensitive). Empty = fingerprint only.
NAME_HINTS = {
    "CS": ["citrate synthase"],
    "ACO": ["aconitase", "aconitate hydratase"],
    "ICDH": ["isocitrate dehydrogenase"],
    "SDH": ["succinate dehydrogenase"],
    "PDH": ["pyruvate dehydrogenase"],
    "GOX": ["glucose oxidase"],
    "GLCt": ["hexokinase", "glucose transport"],
    "CITex": ["citrate"],
    "PPC": ["carboxylase"],
    "PHOSt": ["phosphate transport"],
}


def fingerprint(rxn):
    """Multiset of (formula, coefficient) ignoring compartment suffixes."""
    sig = Counter()
    for met, coef in rxn.metabolites.items():
        f = (met.formula or "?").strip()
        sig[(f, round(coef, 6))] += 1
    return frozenset(sig.items())


def load():
    refined = read_sbml_model(REFINED_XML)
    ijb = read_sbml_model(IJB_XML)
    return refined, ijb


def build_crosswalk():
    refined, ijb = load()

    # index iJB1325 by fingerprint (only reactions with full formulas)
    fp_index = {}
    no_formula = 0
    for r in ijb.reactions:
        mets = list(r.metabolites.items())
        if not mets or any((m.formula or "").strip() in ("", "?")
                           for m, _ in mets):
            no_formula += 1
            continue
        fp_index.setdefault(fingerprint(r), []).append(r)

    entries = []
    for rr in refined.reactions:
        cands = []
        # E1 fingerprint
        exact = fp_index.get(fingerprint(rr), [])
        for hit in exact:
            cands.append({"iJB1325": hit.id, "name": hit.name or "",
                          "evidence": "E1_fingerprint"})
        # E2 name keywords
        hints = NAME_HINTS.get(rr.id, [])
        for kw in hints:
            for r in ijb.reactions:
                if kw in (r.name or "").lower():
                    cands.append({"iJB1325": r.id, "name": r.name or "",
                                  "evidence": "E2_name:%s" % kw})
        # dedupe, prefer exact fingerprint hits
        seen, uniq = set(), []
        for c in sorted(cands, key=lambda c: 0 if c["evidence"].startswith("E1") else 1):
            if c["iJB1325"] not in seen:
                seen.add(c["iJB1325"])
                uniq.append(c)
        entries.append({
            "refined": rr.id,
            "refined_gpr": rr.gene_reaction_rule or "",
            "matches": uniq[:6],
            "n_matches": len(uniq),
        })
    return refined, ijb, entries


def propagate_genes(refined, ijb, entries):
    """refined gene -> refined rxn -> matched iJB1325 rxn -> iJB1325 genes."""
    gene_map = {}
    for e in entries:
        rgpr = e["refined_gpr"]
        if not rgpr:
            continue
        for g in rgpr.replace("(", " ").replace(")", " ").split():
            if g in ("and", "or"):
                continue
            slot = gene_map.setdefault(g, [])
            for m in e["matches"]:
                r = ijb.reactions.get_by_id(m["iJB1325"])
                for gg in r.genes:
                    slot.append({"via": m["iJB1325"],
                                 "evidence": m["evidence"],
                                 "gene": gg.id})
    return {g: v for g, v in gene_map.items()}


def main():
    refined, ijb, entries = build_crosswalk()
    gm = propagate_genes(refined, ijb, entries)

    n_exact = sum(1 for e in entries
                  if any(m["evidence"].startswith("E1") for m in e["matches"]))
    n_any = sum(1 for e in entries if e["matches"])
    n_unmapped = [e["refined"] for e in entries if not e["matches"]]

    out = {
        "direction": "refined_ccm(28rxn,17genes) -> iJB1325(2320rxn,1325genes)",
        "method": "E1 stoichiometric fingerprint; E2 curated name keywords; "
                  "genes propagated through GPRs (no sequence alignment)",
        "counts": {
            "refined_reactions": len(refined.reactions),
            "mapped_exact_fingerprint": n_exact,
            "mapped_any": n_any,
            "unmapped": n_unmapped,
            "genes_with_candidates": len(gm),
        },
        "reaction_crosswalk": entries,
        "gene_candidates": gm,
    }
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)

    print("=== crosswalk written -> %s" % OUT_JSON)
    print("exact-fingerprint mapped : %d/28" % n_exact)
    print("any-evidence mapped      : %d/28" % n_any)
    print("unmapped                 : %s" % n_unmapped)
    print("genes with candidates    : %d/17" % len(gm))
    print("\n--- key spot checks ---")
    for want in ("CS", "ICDH", "SDH", "PDH", "GOX", "CITex"):
        e = next(x for x in entries if x["refined"] == want)
        tops = ", ".join("%s(%s)" % (m["iJB1325"], m["evidence"][:2])
                         for m in e["matches"][:3]) or "-"
        print("  %-6s -> %s" % (want, tops))
    for g in ("citA", "cexA", "goxC", "sdhA"):
        c = gm.get(g, [])
        print("  gene %-6s -> %d candidates, e.g. %s"
              % (g, len(c), [x["gene"] for x in c[:4]]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
