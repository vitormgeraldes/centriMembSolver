# Claude Code — Implement Abstract RO/NF Membrane Classes and Commercial Module Selection in Choupo

**Implementation directive | 31 October 2026 | Project: Choupo**

## 0. Mission and scope

Implement a **two-layer membrane workflow**:

1. **Conceptual design and teaching:** a user selects an **abstract functional membrane class**, not a brand or a part number. Choupo simulates and optimises the RO/NF process against a clearly declared, hypothetical/reference performance specification.
2. **Equipment design / procurement support:** after the process is specified, Choupo compares the design requirements with **separately curated, traceable commercial membranes and modules**, simulates promising candidates under the **actual process conditions**, and recommends the closest documented commercial options with a transparent ranking and limitations.

**The classification and compatibility of chemical components and thermodynamic systems have ALREADY been implemented.** Do **not** reimplement them, create another chemistry ontology, or change the existing component classifier except for a narrowly justified integration fix. Reuse the existing aqueous eligibility checks and the actual thermophysical capabilities of the selected package.

The major engineering objective is to couple **water/solute transport + concentration polarisation + thermodynamics + scaling**. Membrane transport models should be as simple as the data justify. Do not block the useful workflow waiting for an exhaustive proprietary membrane catalogue.

**Act on this specification:** inspect the repository, implement, test, document, and report results. Only stop for a genuinely irreversible or owner-reserved architectural decision. Do not stop at a proposal or create a demonstration disconnected from the existing code.

## 1. Source of truth: inspect first, reuse what Choupo already ships

Read **in this order**:

1. `CLAUDE.md` and the relevant current architectural decisions.
2. Existing dictionaries, schema definitions and `docs/ai/schemas-reference.md`.
3. `docs/ai/unit-ops.md`, especially `spiralWoundModule`.
4. `data/standards/README.md` and existing `data/standards/assets/` membrane, membrane-module and feed-spacer records.
5. Existing implementation in `src/thermo/membrane/`, `src/unitOperations/membrane/`, the GUI and the sizing/design post-processor.
6. Existing regression/curation gates, especially the membrane module catalogue check and spiral-array gates.
7. **Only then**, if needed, consult tutorials as executable examples; do not start by copying a tutorial.

At the time of drafting, the public repository documents these **existing** capabilities (verify their exact current names before editing):

- `spiralWoundModule` (`membraneSW`), one stage with parallel `vessels` and series `elements`;
- `module <name>;` hardware references (`kind membraneModule`), `membrane <name>;` transport records, and `kind feedSpacer` records;
- generic module sizes (`Generic-2540`, `Generic-4040`, `Generic-8040`, `Generic-16040`);
- transport laws `solutionDiffusion`, `SDEM`, `DSPM-DE`, and `constantRejection`;
- `massTransfer` and `pressureDrop` models based on Schock–Miquel, plus ion-coupled polarisation;
- existing property/speciation infrastructure, a rating-test concept (`ratedTest {}`), explicit provenance, and the `data/standards` vs `data/local` separation.

**Do not create duplicate membrane registries, solvers, unit operations, or geometric calculators.** Extend existing record kinds, resolvers, dictionaries and schemas. Keep existing runnable cases compatible unless an existing behaviour is demonstrably incorrect.

## 2. Architectural contract: three distinct things

Do not conflate:

| Layer | Purpose | Owner of data | Permitted content |
|---|---|---|---|
| **Functional class** | Student's process objective / initial idealised membrane | generic model/class registry | abstract performance targets, permitted model choices, default assumptions and uncertainty |
| **Transport membrane** | Physical/selectivity model used by the unit | existing membrane transport model / record | `A_w`, optional fitted transport parameters, model provenance, validity, uncertainties |
| **Commercial element/module** | Hardware candidate offered by a manufacturer | separate `kind membraneModule` + manufacturer-rated tests | brand, model number, area, geometry, operating limits, published test results and sources |

The commercial catalogue must **not** control a conceptual simulation. A conceptual case must run **without** choosing any manufacturer. Commercial selection is an **explicit design-stage conversion** from abstract requirements to candidates, not hidden substitution in a student's original simulation.

Keep the *raw manufacturer specifications* independent of any transport model. Never store a fitted `B` as if it were a manufacturer specification.

## 3. Abstract functional classes (first implementation)

Add a small, **data-driven, editable** class catalogue, with at least these six conceptual service profiles:

| Class ID | Teaching intent | Example screening windows — NOT universal guarantees |
|---|---|---|
| `NF_OPEN` | Pass monovalent salts; partial selective retention | NaCl rejection about 0–30% |
| `NF_DIVALENT_SELECTIVE` | High sulfate/divalent-salt retention with NaCl passage | NaCl rejection about 10–50%; MgSO4 typically >90% |
| `NF_TIGHT` | Stronger salt retention, partial desalination | NaCl rejection about 50–90% |
| `BWRO_HIGH_FLUX` | Brackish RO with relatively high water permeability | NaCl rejection generally >95%; context-dependent pressure |
| `BWRO_HIGH_REJECTION` | Brackish RO for high product purity | NaCl rejection generally >99% |
| `SWRO` | Seawater reverse osmosis | NaCl rejection generally >99.5% at SWRO reference conditions |

These are **pedagogical starting categories**, deliberately not a manufacturer taxonomy. Do not classify only by a product's marketing name or blindly apply the numerical bands to unlike test conditions. A product can plausibly match more than one functional class; present a ranking rather than forcing exclusivity.

A class record must contain:

- stable `classId`, display name, concise teaching description;
- reference *observed* rejections by **solute/salt**, with explicit **reference conditions** and whether a target is point, range, lower bound or illustrative;
- reference water productivity/flux or allowable design ranges **only when defensible**;
- application intent (softening, sulfate control, partial desalination, SWRO, etc.);
- reference/scenario operating envelope and confidence/provenance;
- supported transport modes (prefer simple conceptual approximation first; expose advanced models as optional);
- uncertainty/scenario policy for unspecified solutes, including ions relevant to scaling.

Do not present a class as a physical membrane with fully known selectivity. The student may override the abstract target values. The GUI must make every such override visible.

### 3.1 Conceptual calculation policy

The current `constantRejection` transport law is a suitable **teaching starting point** but has **intrinsic/wall** rejection as its input (`c_p=(1-R_i)c_m`), whereas manufacturer/catalogue values are normally **observed/feed-to-permeate** rejections. These are **not interchangeable** if there is concentration polarisation.

Implement or reuse an explicit **nominal-rating inversion** that maps an observed class target to the intrinsic constant-rejection parameter at a stated reference condition, using the **same film/hydraulic solver**. If insufficient reference inputs exist, declare the value a *hypothetical intrinsic assumption* and do **not** pretend it reproduces a rated test.

For more accurate multi-ion simulations, prefer the existing `SDEM` or `DSPM-DE` solver **where the required independent data actually exist**. Preserve zero-current/electroneutrality and the unit's coupled ion-film model. Do not insert `R_MgSO4` as a universal rejection of each magnesium or sulfate ion; salt tests and multi-ion transport are distinct.

Missing calcium, sulfate, boron or silica selectivity must trigger a documented assumption / sensitivity range, not made-up manufacturer values or a default rejection of 100%.

## 4. Commercial catalogue: separate, source-traceable and multi-test

Leverage the existing `kind membraneModule`, `kind RO`/`kind NF`, and `kind feedSpacer` asset records. Introduce a separate **curated manufacturer-test dataset** only as needed; do not force many incompatible tests into a single `ratedTest` field. Maintain compatibility with the shipped single `ratedTest {}` records; use an optional repeated-test/reference mechanism through a single authoritative schema.

A commercial element record should provide, when known:

- stable internal ID; manufacturer; exact model/variant; technology (RO/NF); element geometry/format (flat sheet, 1812, 2540, 4040, 8040, 16040, etc.);
- active area `m2`, active channel height/spacer thickness `m`, spacer reference, effective leaf length `m`, vessel/element dimensions, operational limits;
- manufacturer spec/document URL and document revision/date;
- geometry provenance (`datasheet`, `other official source`, `measured`, `literature`, `estimated`, `missing`), confidence and status;
- a one-to-many collection of **rated tests** per membrane, each linked to its original document.

Each rated test should store the **original meaning** as well as SI conversions:

- `solute` / analytical salt (e.g. NaCl, MgSO4, CaCl2, Na2SO4, boron as B vs boric acid, silica as SiO2), including basis;
- observed rejection `R_obs` as a fraction, and `valueType`: typical / minimum / maximum / range / approximate;
- nominal feed composition or concentration, pressure and **its meaning** (feed gauge/absolute vs net driving pressure), temperature, pH, recovery, feed flow and permeate flow when reported;
- whether flux is measured in the **same test** or is an independent catalogue productivity figure;
- source URL, revision/date, transcription date, verification status and notes.

Preserve each source result separately. A NaCl test at 2 g/L is **not** the same operating point as an MgSO4 test at 2 g/L. A boron rejection at pH 8 must **not** be transferred to another pH without a transport/speciation calculation. Never average minimum, typical and range values as if they were identical measurements.

### 4.1 Available seed data — import carefully

The working files already assembled for the project are:

- `Choupo_RO_NF_catalog_extended_Synder_TriSep_2026-10-28.xlsx` — legacy/expanded manufacturer catalogue, including loose/selective NF;
- `Choupo_RO_NF_Simulator_Input_SI_extended.csv` — 119 test rows in the working snapshot, with SI conversions and source URLs;
- `Choupo_RO_NF_English_Simulator_Inputs.xlsx` — preceding SI reference work;
- `Choupo_catalogo_RO_NF_consolidado_SI.xlsx` — earlier working catalogue.

These are **staging inputs, not verified standards**. Counts/records were partly inherited from earlier lookups, may contain duplicates and incomplete or outdated manufacturer fields. Locate and inspect the attached files if available; if not, implement the parser and schema and explicitly report that import was not executed, rather than inventing rows.

Include and distinguish loose NF catalogues, notably **Synder** (NFG, NDF, NFDC, NFW, NFDE, NFS, NFD, NFX, NDX) and **TriSep** (XN45, TS40), alongside **DuPont, Hydranautics, LG, Toray** and other verified manufacturer data. A manufacturer page or a reseller is a source, **not an independent experimental validation**.

Implement a deterministic, idempotent import/QA tool that:

1. validates schema and SI dimensions;
2. normalises exact manufacturer/model variants and IDs without silently merging distinct revisions;
3. detects duplicate entries or contradictory rated data;
4. records provenance per field, not merely per row;
5. marks unverified values and conflicting specifications without fabricating resolutions;
6. emits an import/rejection report and counts **computed from real records**;
7. stages unverified data under the project's existing `data/local/` workflow, and promotes to `data/standards/` **only after human curation and according to the existing read-only standards doctrine**.

Do not scrape manufacturers during simulation, embed PDFs/images, reproduce entire proprietary data sheets or publish an unreviewed bulk copy of a third-party database. Store minimal technical facts, citations/URLs, and Choupo's independent calculations. Legal/licence review of a published bulk catalogue is required; do not claim legal clearance automatically.

## 5. Engineering units and the per-element/per-module SI summary

Keep **one compact user-facing SI summary** with only the parameters the Choupo model actually consumes or computes. The simulator must use the actual component feed, not a hardcoded 'seawater TDS'.

**Per membrane element, with clear input/output distinction:**

| Field | Unit | Status and interpretation |
|---|---|---|
| `A_active` | `m2` | source/geometry input |
| `L_leaf` | `m` | source or explicitly estimated |
| `h_feed_channel` | `m` | source or explicitly estimated |
| `epsilon_spacer` | `1` | source or explicitly estimated; **not universally 0.8** |
| `d_h` | `m` | calculated Schock–Miquel, with geometry provenance |
| `A_cross_section` | `m2` | calculated, with the two-membrane-face convention |
| `P_feed_in`, `P_feed_out`, `P_permeate` | `Pa` | inlet + calculated/declared outlet; identify absolute vs gauge |
| `T_feed` | `K` | inlet/input |
| `Q_feed_in`, `Q_feed_out`, `Q_permeate` | `m3/s` | specified or solved; not conflated with membrane-area flux |
| `u_superficial_in`, `u_superficial_mean` | `m/s` | calculated from channel gross cross-section |
| `u_interstitial_mean` | `m/s` | superficial velocity divided by porosity |
| `Jw_mean` | `m/s` (=`m3/(m2 s)`) | area-averaged permeate volume flux |
| `C_feed,i`, `C_bulk_mean,i`, `C_wall_max,i`, `C_permeate,i` | `kg/m3` | per component, never a universal scalar substitute |
| `R_obs,i` | `1` | predicted vs cited nominal, distinguished |
| `SI_bulk,m`, `SI_wall,m` | `1` | mineral-specific saturation indices with model provenance |
| `Re`, `Sc_i`, `Sh_i`, `k_i` | `1`, `1`, `1`, `m/s` | hydraulic/film diagnostics |

The GUI may show a compact subset by default with a details panel for per-solute and per-mineral results. The exported table should have stable machine-readable names and units.

**Important:** nominal test conditions are **inputs/observations**, whereas mean crossflow velocity, mean bulk concentration and wall concentration are usually **model outputs**. Do not store guessed mean module quantities as manufacturer facts. For a one-point provisional rating estimate, mark each derived mean as `estimatedFromNominalTest`.

### 5.1 Spiral-channel geometry: factor 2 is mandatory

A spiral-wound feed channel has membrane area on **both faces**. When `A_active` is the sum of active area on both sides, and `L_leaf` is the effective feed path length:

```text
W_total = A_active / (2 * L_leaf)
A_cross_section = W_total * h_feed_channel
                = A_active * h_feed_channel / (2 * L_leaf)
u_superficial(z) = Q_feed(z) / A_cross_section
u_interstitial(z) = u_superficial(z) / epsilon_spacer
```

Do not accidentally double-count or omit that factor of 2. For **flat-sheet single-face cells** the corresponding divisor is 1, as already documented in Choupo. Reuse the existing geometry code rather than reimplementing the formula in a new GUI helper.

For an approximate linear flow decrease along an element, a provisional average is:

```text
Q_feed_mean ≈ Q_feed_in - Q_permeate / 2
u_superficial_mean ≈ Q_feed_mean / A_cross_section
```

In production, compute averages from the existing axial solution, preferably length/area-integrated, **not** from the arithmetic formula if the axial profile is available. In a train of elements, use local `Q_feed` for each element and divide by the number of parallel vessels in the existing stage formulation.

### 5.2 Schock–Miquel geometry: do not make porosity universal

For a spacer of thickness `h`, porosity `epsilon` and cylindrical filament diameter `d_f`:

```text
d_h = 4 * epsilon / (2/h + (1-epsilon) * 4/d_f)
```

Prefer the actual feed-spacer record and existing Schock–Miquel implementation. If no geometric characterisation exists, allow a **labelled estimation policy** (e.g., `d_f/h = 0.5` and an explicit assumed epsilon), with `source`, `confidence`, `estimated`, and uncertainty/sensitivity range. The existing verified paper's RO spacers have porosities close to 0.88–0.91; **do not globally hardwire epsilon = 0.80**. Use manufacturer-specific or engineering-profile assumptions as appropriate.

If channel height itself is missing, no geometrically justified `d_h` can be derived. Leave it unresolved or use an explicitly selected *generic design module*, never pass an invented value as a commercial fact. Preserve correlation validity-window warnings (particularly Reynolds range).

### 5.3 Unit/basis checks

- LMH -> `m/s`: divide by `3.6e6`.
- `1 m3/(m2 s) = 1 m/s` for volumetric permeate flux.
- `bar -> Pa`: multiply by `1e5`.
- `25 degC -> 298.15 K`.
- `mg/L -> kg/m3`: multiply by `1e-3` **for concentrations explicitly stated per solution volume**; do **not** confuse this with ppm by mass (`mg/kg`) without density.
- Keep osmotic pressure, gauge/absolute pressure and transmembrane pressure difference distinct.
- An observed rejection fraction is **not** an intrinsic `B` or a pore parameter.

## 6. Thermodynamics and scaling are the distinctive capability

Use the chemical classifications **already implemented**; validate aqueous-model compatibility through the existing system, not via hardcoded water-component name tests.

For each local axial station/element, couple:

- solvent chemical potential / osmotic pressure (appropriate aqueous electrolyte activity model);
- water transport and ion/solute transport;
- concentration polarisation and ion coupling;
- aqueous speciation and mineral saturation **in the bulk and at the membrane surface**;
- concentration and pressure changes along each element, vessel and stage.

Report potential scaling for minerals supported by the selected thermo package, e.g. calcite/aragonite, gypsum/anhydrite, barite, celestite, silica and relevant phosphate phases where data exist. **Do not claim that SI > 0 implies immediate deposition**: report supersaturation/risk, and distinguish equilibrium from precipitation kinetics.

For high-recovery NF/RO networks, permit NaCl-permeable NF to be chosen deliberately (reduced NaCl osmotic driving penalty, but not desalinated NF permeate). The actual maximum recovery must be constrained by predicted scaling, hydraulic limits, concentrate flow and required product quality.

If activity parameters or mineral equilibria are missing/outside their validity range, propagate explicit uncertainty or refusal. Do not quietly switch to ideal thermodynamics merely to obtain a design result.

## 7. Workflow: Conceptual -> simulate/optimise -> commercial shortlist -> re-rate

**Stage A — Conceptual:**

1. User chooses `Functional membrane class` (six initial choices) and enters the feed analysis, T/P and process targets.
2. Choupo shows the assumed reference selectivities, flux/range and uncertainty, editable by the user.
3. Solve using existing `spiralWoundModule` for a generic module size/spacing if no commercial model is chosen. Keep abstraction separate from commercial hardware.
4. The student adjusts stage/vessel/element count, recovery, energy and scaling constraints. Display whether the conceptual target is feasible under the declared assumptions.

**Stage B — Candidate search (design/procurement support):**

1. Query **only the independently curated commercial catalogue** by fit to the conceptual duty, desired performance, format and applicable operating envelope.
2. Reject candidates violating hard limits (pressure, temperature, dimensions, declared chemical tolerance, rated max feed flow and hydraulic/pressure-drop limits where known).
3. Use rated tests to identify or estimate a compatible transport model **under the candidate's original test conditions**; then re-simulate the candidate at the student's process conditions. Do not compare raw nominal rejections at incompatible P/T/concentrations as though they were directly commensurate.
4. Rank feasible candidates by process performance and confidence: required product composition, flux/area, pressure/energy, minimum concentrate flow, bulk/wall scaling margin, necessary module/vessel count. Apply user-visible weights or constraints, never a magic hidden score.
5. Show **top three documented candidates** (or fewer when insufficient data), exact manufacturer/model IDs, URLs, test conditions, predicted operating point, deviations from conceptual target, important missing fields and confidence.
6. Keep supplier part numbers/trademarks in this final selection only. No fabricated current price, vendor availability, purchase link or certification.

**Stage C — Committed commercial rating/design:**

On selecting a candidate, bind the existing `module <commercial-id>` and `membrane <transport-id>` in a **new design/rating version** of the case. Run the existing solver again; show the impact relative to the abstract ideal: permeability, salt passage, scaling margin, recovery, pressure/energy, dimensions and number of elements. Do not silently rewrite the original conceptual case.

If the published data are too sparse to identify a model, give a *screening-only ranking* with a visible warning and a range of outcomes. Do not present a numerical fit as uniquely identified from a single NaCl rejection datum.

## 8. Ranking formula — transparent, not self-fulfilling

Build a dimensionless score from **normalised performance at a common operating point**, with hard eligibility filters first. An illustrative structure:

```text
score(candidate) =
    w_quality * normalised_product_quality_penalty
  + w_flux    * normalised_area_or_flux_penalty
  + w_energy  * normalised_pressure_or_energy_penalty
  + w_scaling * normalised_scaling_risk_penalty
  + w_fit     * normalised_target_selectivity_penalty
  + w_data    * documented_uncertainty_penalty
```

Publish the score breakdown and weights. Missing data must reduce confidence, not silently become zero penalty or an impossibly good prediction. Respect constraints before optimising the score. Never fit parameters on the same experimental points that are then called *independent validation*.

Support comparison of an optimised **abstract** NF/RO system with the best available **commercial** alternative, rather than merely returning the nearest name based on two raw rejection percentages.

## 9. GUI, API, F1 and agent-facing documentation

Expose this as a coherent extension of the existing sequential-modular flowsheet UX:

- **Membrane selection mode:** `Functional class` (default in teaching) or `Specific commercial module` (advanced/rating), not a confusing mixture of both.
- Compact functional-class selector with a one-line explanation and reference performance assumptions.
- A clearly separated **"Find commercial equivalents"** action in the **design** stage, not during basic composition entry.
- One concise SI operating-point table and a small number of performance plots (axial flux, pressure, major species, bulk/wall SI).
- Commercial comparison table: candidate, provenance, test condition, predicted process results, constraints, uncertainty and datasheet URL.
- Explain each value's state: `published`, `measured`, `estimated`, `assumed`, `fitted`, `simulated`, `unverified`. Show editable assumptions without advertising them as product facts.
- F1 examples: why NaCl passage can help selective NF, why high MgSO4 rejection does **not** guarantee CaCO3/CaSO4 control, why an observed manufacturer rejection differs from intrinsic rejection.
- Generate schemas/help from the same sources of truth; update agent-facing documentation, not a second handwritten option list.

## 10. Regression, physical checks and acceptance criteria

Create minimal focused tests and a small teaching case without duplicating existing tutorial suites. Cover at least:

1. **Abstract stage independent of vendor.** A complete conceptual NF selective-sulfate example runs with no manufacturer or part number.
2. **Old cases stay valid.** Shipped `spiralWoundModule` cases, known `NF270-4040`, `SW30HR-380` and generic modules retain valid behaviour unless an explicit physics correction changes them with independently justified baselines.
3. **Factor-of-two check.** With `A_active=35 m2`, `h=0.0007112 m`, `L_leaf=0.95 m`, the computed gross feed cross-section is approximately `0.013102 m2`; a flat-sheet, single-face interpretation must be twice that under otherwise identical inputs.
4. **Schock–Miquel check.** For `h=0.7112 mm`, `epsilon=0.90`, `d_f=h/2`, verify `d_h≈0.9144 mm` and use the already-implemented spacer test as authority.
5. **No epsilon universal default.** Two spacer records with different porosities produce distinct `d_h`, interstitial velocity, mass transfer and pressure drop.
6. **Observed vs intrinsic rejection.** Nominal `R_obs` is not directly injected as `constantRejection` intrinsic R when the film is active; inversion/assumptions are reported.
7. **Salt tests not interchangeable.** Separate NaCl and MgSO4 conditions remain separate, and boron requires the specified pH/basis; 97% MgSO4 must not imply 97% CaCl2.
8. **Aqueous compatibility.** Reuse the completed component classification system: an aqueous RO unit refuses a completely organic feed, whereas another supported operation may process it.
9. **SI and balances.** All table fields have dimensional tests; stage and element mass balances/electroneutrality pass; pressure/flow/area do not get swapped.
10. **Scaling at two locations.** Report `SI_bulk` and `SI_wall`, demonstrate a case where they differ, and never equate supersaturation with measured deposition.
11. **Commercial lookup.** Deterministic top candidates; rejection of out-of-envelope products; unavailable information clearly flagged; links trace back to exact manufacturer tests.
12. **Idempotent import / catalogue QA.** Repeated ingestion does not duplicate tests; conflicting source revisions are surfaced, not overwritten.
13. **Independent validation doctrine.** Preserve distinction between regression, model verification, rating consistency and independent experimental validation.
14. **Full suite:** run the repo's official test gates and accept only `FAIL 0`; regenerate inventory/schema artefacts by their prescribed generators, never maintain copied counters by hand.

One pedagogical end-to-end case should run:

```text
Aqueous feed -> generic NF_DIVALENT_SELECTIVE -> thermodynamic scaling and recovery
-> commercial shortlist (e.g. a curated Synder/TriSep candidate)
-> rerun with selected actual module -> compare designs.
```

Use actual source-vetted product facts in the final part; if the attached catalogue is not accessible, supply the workflow and tests with explicitly hypothetical fixture records rather than claiming a manufacturer endorsement.

## 11. Implementation order and outputs

Execute as a series of narrow, verifiable commits (not a broad rewrite):

**A. Inspection and mapping:** identify current record and GUI interfaces, existing generic classes/models, resolver and relevant gates; record what can be reused.

**B. Data model/import:** add the conceptual functional-class records and the multi-test commercial metadata adapter; import staged CSV/Excel when available; emit field-level QA report.

**C. Solver integration:** expose class targets through the existing transport law/rating workflow, enforce geometry and SI summary, integrate thermodynamic scaling diagnostics; do not destabilise the existing physics.

**D. Design-stage matchmaker:** candidate hard filters, re-rating at actual conditions, uncertainty-conscious ranking, recommendation output.

**E. GUI + help:** simple conceptual entry, separate commercial selection, clear warnings, SI summary, F1 and agent/schema help.

**F. End-to-end gates and documentation:** pass the complete official suite, distinguish verified from unverified catalogue records, and demonstrate the pedagogical case.

### Final report to the project owner

Deliver:

- paths and short purpose of each modified/new file;
- number of **actually imported** unique membranes/modules and independently sourced rated tests (computed, not copied from a conversation);
- number of vendor records that are verified, partially verified, and unverified;
- explicit treatment of porosity, spacer geometry and nominal vs intrinsic rejection;
- a working conceptual case and a working commercial-selection case;
- before/after comparison of the best available commercial alternative, clearly noting unvalidated estimates;
- exact tests run and results, with `FAIL 0`, or precise uncompleted blockers.

**Definition of done:** an engineering student can select an **abstract membrane behaviour**, simulate a chemically valid RO/NF process with thermodynamic scaling predictions, then — only during design — receive a **traceable, uncertainty-aware shortlist of commercial modules**, choose one and re-simulate the real equipment without modifying the original conceptual model.