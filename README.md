# Linear coordinate presentations of Lotka–Volterra tree systems

Author: Anonymous. Paper: `lotka_volterra_tree_classification.pdf` (source `.tex`, 22 pages).

## In plain English

A Lotka–Volterra system is a quadratic differential equation used for interacting
populations. Some of these equations have extra "invariant planes", given by linear
Darboux polynomials. For tree systems there is one extra plane for each edge of a tree.
If you change coordinates, the same equation can look like a different picture
(a hypergraph instead of a tree). Van der Kamp conjectured that two different trees
can never describe the same equation, for generic parameters. This paper proves that
conjecture, in a stronger form that needs genericity on one side only. It also shows by
an explicit example that the statement does not extend to special parameter values. Along the way it lists every
coordinate system in which a tree system keeps its form, and it counts how many
genuinely different pictures there are.

## What is proved, what is computed, what is open

- **Proved for every tree** (full proofs in the paper):
  - the complete list of linear Darboux polynomials, with a genericity condition that is necessary and sufficient (Thm 2.3);
  - the codimension-one exceptional strata (Prop 2.5);
  - the tangent-field description (Prop 3.2);
  - one-sided tree rigidity, i.e. the conjecture in readings (G) and (O) (Thm 3.6, Cor 3.7);
  - intersecting tree families: the all-parameter extension (A) of the conjecture fails (Prop 3.8); the conjecture itself, which concerns generic parameters, is unaffected;
  - all coordinate presentations (Thm 4.2);
  - the symmetry group (Thm 5.2) and the core-forest reduction (Lemma 6.2, Thm 6.3);
  - the path and caterpillar counts (Thms 7.1, 7.2);
  - the recognition test (Thm 8.1) and the Kahan-map transfer (Prop 8.3); constant growth rate b with the modifications of Remark 2.6.
- **Exact computation** (finite, exact arithmetic):
  - all trees n ≤ 9 via cone bases, and n ≤ 10 via core forests;
  - agreement with van der Kamp's Table 7 for n ≤ 8;
  - checks of every new statement for small n.
- **Exploratory**: van der Kamp's 13-parameter system (34). At 5 sample points it has 9 rays, but its matroid is not a tree-cone matroid.
- **Open**:
  - multi-factor exceptional strata;
  - whether the general class of each hypergraph equals the tree presentations, for all n;
  - formulas beyond caterpillars;
  - quasi-monomial equivalence;
  - completeness for Kahan maps.

## Replay

Requires Python 3 with `sympy`. Tested with Python 3.13 and sympy (venv `analysis-py313`).

    python3 verify_all.py          # ~25 s; reruns every verifier and compares bytes with data/
    python3 verify_all.py --write  # regenerate data/

Outputs are deterministic: fixed seeds, sorted keys, no timings or paths.

| Script | Output | Content |
|---|---|---|
| `code/verify.py` | `data/verification.json` | Darboux supports, cone bases, orbit counts, Table 7 (n ≤ 9) |
| `code/verify_hypergraphs.py` | `data/direct_hypergraph_verification.json` | direct n! hypergraph canonicalisation (n ≤ 6) |
| `code/lv_core_count.py` | `data/lv_core_counts.json` | decorated core-forest counts (n ≤ 10) |
| `code/verify_extensions.py` | `data/extensions_verification.json` | checks (a)–(k) of the new results (Section 9 of the paper) |

Integrity: `shasum -a 256 -c SHA256SUMS`.

Data and code DOI: [10.5281/zenodo.23564794](https://doi.org/10.5281/zenodo.23564794) (version 1.0.0-candidate); all versions: [10.5281/zenodo.23564793](https://doi.org/10.5281/zenodo.23564793). Repository: https://github.com/ipitchford/lotka-volterra-tree-presentations

## Release records

Version 1.0.0-candidate (unrefereed). [STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md), [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md), [PROVENANCE.md](PROVENANCE.md), [AI_INDEX.md](AI_INDEX.md), [CLAIMS.json](CLAIMS.json), [RESEARCH_GATES.json](RESEARCH_GATES.json), [LICENSES.md](LICENSES.md) (CC0-1.0 prose and data, MIT code), [CITATION.cff](CITATION.cff). Semantic negative controls: `python3 -B -O negative_controls.py` (standard library). Evidence Press page: https://evidencepress.org/releases/lotka-volterra-tree-presentations/
