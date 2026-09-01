# Session log

## 2026-09-02 (later) — Reorganization executed (Tiew syllabus arrived)

Tiew's ECON-GA 3002 syllabus (local, untracked) settled the open questions:
she covers entry, moment inequalities, dynamic games, firm learning — so those
stay pointers here (L15 Dickstein guest slot: decide per availability).
Executed with Chris's rulings:

- **Pass-through → Week 7**: new `passthrough.tex` (29 frames) ported from the
  Leuven deck with notation conversions and TWO source sign-error fixes in the
  Villas-Boas/IFT algebra; `Week 15b/pt_notes.tex` fully subsumed and retired.
- **Common ownership → half-lecture**: `conduct_new.tex` gained a 13-frame
  BCS application section (hypothesis, Big Three, κ derivation, reduced-form
  critiques, cereal data/κ variation, implementation callbacks, results,
  mechanism), integrating the old 3-frame stub. 46→57 pages.
- **`syllabus-2026.tex` drafted**: corrected Other Courses rotation +
  division-of-labor statement (Nguyen 1802 / Tiew 3002 / Waldinger EMD);
  Week 5 retitled with a Pass-Through reading block (Weyl–Fabinger, MRRS,
  Mrázová–Neary); 2023 Merger Guidelines; zeros subsection; CG 2020/2023;
  Small–Rosen; BCS AEJ:Micro. Flag: cereal BCS still cited "(2021, WP)".
- **Website (chrisconlon.github.io, local edits, not pushed)**: fixed the
  broken Common Ownership link (pointed at the DP notes); added the Lecture 7
  Pass-Through entry; repointed Other-Notes Pass Through to the Week 7 deck.
  Fall-2026 header/dates refresh still pending.
- Updated `other_courses/course-overlap-notes.md` with the Tiew section and
  the full sequence division-of-labor table.

## 2026-09-02 — PS0 expansion + integration notes update

- **PS0 → 5 parts**: Part 2 gained the 4-D curse-of-dimensionality exercise
  (1-D ridge reduction as truth; nquad timed at moderate tolerance as the
  "adaptive doesn't scale" lesson; tensor GH vs sparse grids — both routes:
  sparse-grids.de GQN download or hand-rolled Smolyak — vs matched-node MC;
  monomial validation required; new (f) on why tensor wins at d=4 and dies at
  d=10). New Part 3: Berry-inversion fixed point (contraction vs Newton;
  returns as BLP inner loop + Rust Bellman). New Part 4: numerical derivatives
  (forward/central vs analytic; truncation-vs-roundoff log-log plot).
- **Solution keys (local, gitignored)**: ps0_solution.py (Smolyak combination
  technique validated against published Heiss–Winschel GQN node counts — d=10
  K=5 gives exactly 8,761; runs ~3 min incl. the deliberate 184s nquad row)
  and ps1_solution.py + generate_schools.py from 2026-09-01.
- **Design findings baked in**: sparse grids initially LOST to MC because the
  original sigmas made the integrand a wide ridge (negative-weight
  cancellation, sum|w| up to 833); rescaled sigma=(1,.5,.75,.5) gives the
  clean sparse-beats-MC table while tensor's d=4 dominance becomes question (f).
- **Extra Notes/integration.tex updated**: Simpson frame had the TRAPEZOID
  figure (proper figure generated); date/typo/N(mu,Sigma) errata; QMC frame
  gains MLHS + scipy.stats.qmc + pyblp.Integration; three new frames — Sparse
  Grids (Smolyak), Sparse Grids: Caveats, Know Your Integrand (ridge
  reduction) — carrying the PS0 lessons with cross-references both ways.

## 2026-09-01 (later) — Assignments overhaul

All five active psets repaired, modernized, and extended; all build. Uncommitted
pending review.

- **Solutions privatized**: PS0's comment-block answers and PS4's
  `rust.py`/`generate_data.py` (the DGP!) moved to gitignored
  `Assignments/solutions/`; public tree now ships data only.
- **Retired**: MATLAB PS0 variant (stale, contradicted the convexity fix),
  Su–Judd Knitro/AMPL folders (PS4 now points at scipy/CasADi/Pyomo/JuMP),
  orphan `BLP_hw_2024.pdf`.
- **Repairs**: PS4's malformed Hotz–Miller displays fixed (log-sum-exp over
  actions, primed states, Euler-constant comment) + β=0.975 stated; PS2's
  `\calJ` compile error, WTP typo, Part-IV sign, λ→ρ relabel (file's λ was
  already Berry's parameter — exponents kept `1/(1−ρ)`), continuous numbering
  1–22, `Σ_η` notation; PS1's ξ₁=0 normalization + data-description fixes;
  PS0's URL/hermgauss/notation refresh.
- **New material**: SE questions (inverse Hessian/OPG; GMM sandwich with
  (1+1/S)) and a school-closure ΔCS question (Small–Rosen, miles-as-money,
  log(1−σ) comparison) in PS1; zeros reframe of PS2 Part I (it *is*
  Quan–Williams) + selection question + second-choice diversion question +
  CES coda; sampling-zeros part (7a–7c: multinomial N_t=500, drop vs Laplace,
  bias directions) and micro-moment question (9a, salvaged from ps2_2016's
  commented draft, CG-2023 standardized form) in BLP_hw; diversion-object
  (D_jk vs D_{j→k}) reporting requirements in both diversion tables.
- **PS1 identification flaw → RESOLVED as pedagogy** (Chris's ruling): Q4 now
  estimates the characteristics model (ξ=0), Q6 the FE model, with an explicit
  nesting/LR structure (2 df); new Q7 asks why joint estimation fails (FEs
  absorb school-level characteristics — exactly BLP's ξ) and runs the two-step
  ξ̂-on-x projection with a small-J standard-errors caveat. RC question builds
  on the characteristics model.
- Due dates deliberately untouched per Chris.

## 2026-09-01 — Weeks 1–5 streamlining (see `REDESIGN-WEEKS-1-5.md` for full detail)

Three phases, all complete and verified building; uncommitted pending review:

1. **Errata + retirement**: `demand2.tex` retired (superseded ancestor of Week 3;
   one salvage — the corrected heteroskedasticity display — ported to mc1) with
   its ghost copies; 2023 Merger Guidelines thresholds; bertrand reordered +
   merged; rep_consumer CES exponent fix + figure-frame titles; AIG initials;
   part1b title; Vives 1999; assorted brace/typo fixes.
2. **Within-deck merges**: mc1/mc2/mc3 consolidated (mc3: 482→203 lines,
   commented graveyard cleared); demand_new review section compressed, pseudocode
   /optimal-IV merges, cross-equation supply frame restored, argmax→argmin GMM
   fix; micro_data citation tables merged (identical paper sets!), Comparisons
   frame restored.
3. **New frames**: WDZ theorem + Small–Rosen CS/CV worked removal example
   (ΔCS = (1/α)log(1−σ_j) identity; virtual-price demystifier per Chris) +
   CES↔logit bridge in mc1; nested-CES bullet (mc2); preview bullet
   (rep_consumer); zeros-in-shares ×2 (Gandhi–Lu–Shi, DHJ, Quan–Williams),
   PCM + vertical-model asides (archive ports, CDF bug fixed), "Taking Stock"
   IIOC retrospective (Good/Bad News), search/consideration preview — all in
   demand_new (66→75 pp).

Separately: syllabus-overlap analysis of ECON-GA 1802 vintages (S2025 Waldinger,
S2026 Nguyen) in `other_courses/course-overlap-notes.md`; reorganization ON HOLD
pending third syllabus.

## 2026-08-31 — Notation unification + infrastructure

Work done with Claude Code; commits `ca89c47`, `c7b52f7`, `c13affa` (all pushed).

## Housekeeping (`ca89c47`, `c7b52f7`)

- Deleted 433 untracked LaTeX aux artifacts and `.DS_Store` files; removed the
  stray macOS partial-transfer file `.!22825!switching_costs2.tex`, scratch
  `test.tex`/`test.pdf` (Week 12), and the empty `Week 5- Instruments/` folder.
- Untracked `build_log.txt` (now gitignored); deduped `.gitignore`.

## Notation unification (`c13affa`) — see `NOTATION.md` for the conventions

Berry–Haile (2014) conventions applied across Weeks 2–5, 7, 8, Assignments 2–3,
and the Leuven demand/conduct decks (Week 4-Old archives and `BLP_hw_2019.tex`
deliberately untouched):

- `σ_ijt` / `σ_jt(·)` for individual and market-level choice probabilities;
  `s_jt` strictly observed shares (retired `𝔰`, `\calS`-as-shares, `\tilde s`,
  `P_ij`, `S_ij`, `π`-as-probability).
- `ρ` for nesting (Berry's σ); **fixed the nested-logit inversion sign**
  (`−ρ ln s_{j|gt}`, so ρ ∈ (0,1) is coherent; iioc deck was right, the two
  `demand_new` copies were wrong).
- Capital `Σ`/`Π` for RC and demographic parameters (`β_i = β + Σν_i`);
  `w_i` type weights; `π_f` profits only.
- Diversion: `D_jk` (small-price-change) vs `D_{j→k}` (second-choice) for
  **aggregate** objects; individual diversion is intervention-invariant
  `D_{jk,i}` (no arrow).
- Scale parameter `σ → a`; Ho admission/diagnosis probabilities `p → φ`
  (wtp deck); `P_j`-as-price → `p_j`; `Prob(` → `\Pr`; Manski–Lerman bullet
  harmonized to `A_j` (population) / `H_j` (sample rate).
- Semester dates → Fall 2026 (incl. the "Fall 2029" typo); course-intro
  placements updated (Nano Ochoa → Johns Hopkins, Helena Pedrotti → FTC).

## Infrastructure (`c13affa`)

- **`teaching_slides.sty` fork healed**: builds resolve the sty via
  `~/Library/texmf`, which held a *diverged copy* (eucal + R listings vs the
  repo's never-active mathrsfs experiment). Reconciled to the texmf behavior +
  the `\Pr`→ℙ fix and new `\prob{}` macro; **texmf entries are now symlinks to
  the repo files** (recipe in `resources/README.md`).
- **`\input{preamble}` pattern retired**: all 22 such decks (Leuven, Extra
  Notes, Week 15/15b/15c, old Week 7 `conduct.tex`, Week 4-Old
  `instruments`/`supply`) migrated to `\usepackage{teaching_slides}`; root and
  folder-local preambles deleted; `resources/README.md` rewritten. (Root
  preamble had required a `gradio-preamble.sty` that existed nowhere — these
  decks had been unbuildable.)
- Latent defects fixed en route: missing `\end{document}` (Leuven
  `micro_data.tex`), dutchcal/unicode-math XeTeX crash, undefined
  `\custombullet`, stray `A` for `\calH` in the Leuven FOC, `Ω → Δ`.
- **Full build: 61/61 decks pass** `build_all.sh`; all PDFs regenerated.

## Leuven / content findings

- The **"30 Years of BLP" talk** is `Leuven-Lectures/demand/demand_iioc.tex`
  (IIOC 2024) — nondescript filename is why it was lost.
- Incorporation status: micro_data and conduct fully incorporated (weekly decks
  are supersets); ML/low-rank section → `machine_learning.tex`; CMA
  diversion-as-MTE → `diversion.tex`. **Not incorporated:** pass-through (~22
  theory frames; only 6-frame `pt_notes.tex` in Wk 15b), common ownership (~30
  BCS cereal frames; only the 3-slide stub), and the IIOC "Good News / Bad News
  about BLP" retrospective frames (natural Week 4 additions).

## Syllabus comparison

See `other_courses/course-overlap-notes.md`. Headline: Nguyen's S2026 ECON-GA
1802 dropped entry / moment inequalities / dynamic games and added vertical
contracts & bargaining — colliding with my planned vertical expansion and
orphaning entry/MI. **On hold until the third syllabus (Tiew?) arrives.**

## Open items

1. **Dead copy-paste / commented-out block cleanup** across decks (explicitly
   queued; kept separate from the notation pass for diff reviewability).
2. `demand_iioc.tex` line ~832: the IPDL row mixes `σ_j^{-1}` with
   `−x_jt β − α p_jt … + ξ_jt` — looks like ξ-recovery mislabeled as the
   inversion; needs a substantive look.
3. Stale assignment due dates: `ps0-matlab.tex` (Sept 2019), `ps4.tex`
   (Nov 2020).
4. GMPS `\hat{Π}` GMM-penalty vs `Π`-demographic-matrix collision (micro_data
   decks + iioc) — decided to let slide; noted here for completeness.
5. Content phase: pass-through + common-ownership incorporation, IIOC
   retrospective frames into Week 4, Vertical/Bargaining (pending overlap
   decision), Production Functions (Week 6 has no .tex source).
6. `other_courses/` syllabi PDFs are **untracked** — decide whether colleagues'
   syllabi should live in this public repo before committing them.
7. On a new machine: recreate the texmf symlinks per `resources/README.md`.
