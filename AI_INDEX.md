# Agent navigation

## Identity

- Release slug: `lotka-volterra-tree-presentations`; version `1.0.0-candidate`; DOI 10.5281/zenodo.23564794 (concept DOI 10.5281/zenodo.23564793)
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md))
- Manuscript: [PDF](lotka_volterra_tree_classification.pdf), [LaTeX source](lotka_volterra_tree_classification.tex)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json); research target: [target.yaml](target.yaml)
- Assurance limits: [ASSURANCE.md](ASSURANCE.md); internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance: [PROVENANCE.md](PROVENANCE.md); prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json)
- Citation: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json); licences: [LICENSES.md](LICENSES.md)
- Human-readable page: https://evidencepress.org/releases/lotka-volterra-tree-presentations/

## Scientific map

Object: homogeneous Lotka–Volterra systems x_i' = x_i (Ax)_i with A in the tree family L_T of a tree T on n vertices (3n−2 parameters: A_jj = a_j, A_ij = b_(v→j) with v the neighbour of j towards i). Genericity polynomial Delta_T: product of all d_(v→j) = b_(v→j) − a_j and of (d_(u→j) − d_(v→j))(d_(u→j) + d_(v→j)) over neighbour pairs. Source: van der Kamp, J. Dyn. Diff. Equat. 38 (2026), arXiv:2411.18264v4 numbering; Conjecture 12 under the source's generic-parameter convention.

- **Complete Darboux list** (Theorem 2.3, proof): exactly the n coordinates and one ray per edge (2n−1 rays) iff Delta_T ≠ 0. Codimension-one strata: Proposition 2.5.
- **Conjecture 12 proved, one-sided** (Theorem 3.6, Corollary 3.7, proof): U_T-system conjugate to any L_T'-system implies T ≅ T'. All-parameter extension fails (Proposition 3.8); the conjecture is unaffected.
- **Presentations** (Theorem 4.2): spanning trees of Cone(T), det(I + L_T) of them; supports are fundamental paths. **Symmetry** (Theorem 5.2): Aut M(Cone T) = (C_2)^leaves ⋊ Aut(T). **Core forests** (Theorem 6.3). **Counts**: paths (F_(2n−2) + F_(n−1))/2 (Theorem 7.1), caterpillars (Theorem 7.2).
- **Consequences**: tangent fields of dimension 2n−1 (Proposition 3.2), recognition test (Theorem 8.1), Kahan maps (Proposition 8.3), positivity (Remark 8.4), constant growth rate (Remark 2.6).
- **Exact computation**: all trees n ≤ 9 (cone bases), n ≤ 10 (core forests); van der Kamp's Table 7 matched for n ≤ 8. **Exploratory**: the 13-parameter system (34), five sample points.

Do not read: "the conjecture is false" (only the all-parameter extension fails); the counts as proved equal to the source's classes for all n (open); finite computations as proofs of the general theorems; results as applying to quasi-monomial or nonlinear changes of variables.

## Replay

Python 3 with SymPy:

```sh
shasum -a 256 -c SHA256SUMS
python3 -B -O negative_controls.py   # exact, standard library; corruptions rejected
python3 -B verify_all.py             # four verifiers, every output byte for byte, about 25 s
```

Programs: [code/](code/); outputs: [data/](data/).

## Useful next research question

Classify the linear Darboux sets where two or more factors of Delta_T vanish, and decide which special systems of different tree families are linearly conjugate.
