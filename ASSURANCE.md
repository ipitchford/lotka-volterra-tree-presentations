# Assurance

| Dimension | State | Note |
|---|---|---|
| Written proofs | present | Every theorem, all finite trees |
| Exact computation | passed | `verify_all.py` regenerates every output in `data/` byte for byte (about 25 s) |
| Semantic negative controls | passed | `negative_controls.py` (standard library): true instances accepted, plausible corruptions rejected |
| Clean replay | passed | A clean full replay of the final package in a separate copy reproduced every output; the released code and data were bound to it by SHA-256. Publication changed only the manuscript text (DOI statement), the PDF and README; those checks were rerun |
| Internal review | done | Internal AI-assisted review in several roles, an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol) and a verification re-review |
| External peer review | not assessed | |
| Independent reproduction | not assessed | The core-forest enumeration is algorithmically distinct but written in the same project |
| Formal verification | not assessed | |
| Priority | bounded | Dated bounded searches, see [RESEARCH_GATES.json](RESEARCH_GATES.json); not exhaustive |

'proof' means a written proof in the manuscript, valid for all finite trees, not formal verification; 'exact computation' means finite exact-arithmetic computation confirming instances; 'exploratory' is computation at sample points supporting a statement not proved in general; 'open' is not claimed. Internal AI-assisted review is not peer review or independent reproduction.
