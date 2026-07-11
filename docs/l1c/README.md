# L1C — Geometric Processing

Georeferencing, orthorectification and band co-registration: L1B radiance in sensor
coordinates onto a projected map grid. **Partly implemented** (`spp-l1c`) — georeferencing
and orthorectification work; band co-registration reaches 2–3 px, not the sub-pixel
target, and the absolute accuracy is unvalidated. Read
[limitations.md](limitations.md) before trusting a product.

| Document | Contents |
|---|---|
| [spec.md](spec.md) | The design: sensor model, terrain, resampling, refinement, error budget |
| [findings.md](findings.md) | **What measurement overturned** — six assumptions, and how each was caught |
| [architecture.md](architecture.md) | Component design and data flow |
| [validation.md](validation.md) | What was validated, against what — and what was not |
| [qa_report.md](qa_report.md) | `qa_report_l1c.json` schema: metrics, flags, thresholds |
| [performance.md](performance.md) | Measured runtime and memory |
| [limitations.md](limitations.md) | **The product's honest state** |

If you read only one, read [findings.md](findings.md). Five of the specification's
load-bearing assumptions were wrong, four of them were indistinguishable from an
irreducible platform error, and every one was caught by a check that could have failed.

Cross-cutting documents live one level up: [product hierarchy](../product_hierarchy.md),
[decision log](../decision_log.md), [generalisation](../generalisation.md),
[references](../references.md).
