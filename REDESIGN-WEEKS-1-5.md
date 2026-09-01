# Redesign: Weeks 1–5 (demand sequence)

*Working design doc, started 2026-09-01. Companion to `SESSION-LOG.md`.*

**Budget:** Lectures 1–5 of 15 (per the course webpage): L1 = Week 1's three decks,
L2 = bertrand + representative_consumer, L3 = multinomial_choice 1–3 (one slot!),
L4–5 = demand_new + micro_data + machine_learning. The binding constraint is L3.

**Phases (agreed):** (1) errata + demand2 retirement → (2) within-deck merges →
(3) new frames. Commit at phase boundaries.

---

## Phase 1 — errata + retirement  ✅ (in progress)

- ✅ `demand2.tex` retired (+ PDF, ghost `resources/multinomial-choice.tex/.pdf`,
  stale `mpec.pdf`, duplicate `hardmax/softmax/nesting.png`). It was the ancestor
  of all of Week 3, superseded frame-for-frame; only salvage was the corrected
  heteroskedasticity display, ported to `multinomial_choice1.tex`.
- ✅ `bertrand.tex`: reordered (roadmap frames first, technical core last), merged
  duplicate frames (10 → 8), simultaneity bullets → one-line Working (1927) callback.
- ✅ `representative_consumer.tex`: CES price-index exponent fixed to $1/(1-\sigma)$;
  desiderata duplicate → callback; Hausman figure frames retitled by table
  (missing `beer4.png` = Light-segment table; Premium + Popular now read as the
  intended examples); CGJ frame titles disambiguated.
- ✅ Week 1: 2023 Merger Guidelines thresholds (HHI 1800 / ΔHHI 100; 30% share
  presumption, Philadelphia Nat'l Bank); brace fix `part1a:131`; Vives 1999 in
  both decks; part1b retitled ("Structure, Conduct, Performance (and Conjectural
  Variations)") + duplicate HHI regression dropped; AIG initials normalized +
  Takeaways/Nonlinear-IV/Simultaneity merges in the homogeneous-products deck.
  All Week 1 decks rebuild cleanly.

## Phase 2 — within-deck merges  ✅ (completed 2026-09-01; all decks build)

**multinomial_choice1** — merge Proportional Substitution + Diversion Ratios;
merge the two elasticity frames; de-duplicate scale material between the two
scale frames (keep both locations, cut the restatement); delete commented
duplicates (incl. the retired MLE/aggregation block — `demand_new` owns that
material live).

**multinomial_choice2** — delete the IIA-recap opener (mc1 now ends there); keep
the probit "can we do better?" frame as the true opener with a one-line recall;
retitle the mis-titled "Alternative Interpretation" tree frame; merge the two
probability-derivation frames if genuinely redundant.

**multinomial_choice3** — consolidate the triple statement of the RC
specification; delete commented blocks owned elsewhere (FKRB → demand_new,
convexity/FIML/derivatives → `Extra Notes/mpec.tex`); fold
quadrature/curse-of-dimensionality content into the live Recommendations frame;
keep the Chamberlain-scores bullet with an explicit "much more on this next
week" pointer.

**demand_new** — compress the opening review section (McFadden/MLE/aggregation,
~8 frames) to a ~2-frame bridge (L3 immediately precedes); merge
Pseudocode+Estimation and the two Optimal-IV frames; single Bertrand-FOC
statement with callback to the bertrand deck; de-duplicate the ρ→1 bullet and
inversion re-listing; **restore** the commented cross-equation "what supply
tells us about demand" frame (notation-checked per `NOTATION.md`); add a
"(defined next lecture)" pointer where the Conlon–Rao case study cites micro
moments.

**micro_data** — merge the two large citation tables; **restore** the commented
MDPLE-vs-microBLP "Comparisons" frame (notation-checked). *Leuven twins are
archived talks — do NOT sync content restructuring to them (notation only).*

## Phase 3 — new frames  ✅ (completed 2026-09-01; all decks build)

Implementation notes: mc1 gained WDZ + CS/CV pair (incl. the $\Delta CS =
\frac{1}{\alpha}\log(1-\sigma_j)$ removal identity) + CES↔logit bridge closer,
and the Inclusive Value frame's "globally concave" claim was corrected to
convex; mc2 gained the nested-CES bullet; representative_consumer gained the
preview bullet; demand_new gained 8 frames (zeros ×2, PCM + vertical asides —
the archived vertical frame's swapped CDF cutoffs were corrected in the port —
and the Taking Stock section: Good News, Bad News ×2, search preview).
**Flag for Chris:** the search-preview frame cites Honka–Hortaçsu–Vitorino
(2017) as banking/awareness (correct for that paper) — say the word if you
intended Honka (2014, insurance) instead.

### Original rulings (2026-09-01)

1. **CES↔logit bridge — YES**, plus nested/mixed CES if possible.
   Sketch: one frame closing mc1 ("CES and Logit: the same model?" — ADT (1992)
   representative-consumer foundation; both deliver constant markups; table
   comparing p = mc/ρ with η = −αp(1−σ)); one bullet/frame in mc2 pairing
   nested logit ↔ nested CES; forward pointer from representative_consumer's
   "better left for Trade and Macro" frame. Mixed CES = random-σ analogue (trade
   lit) as a remark.
2. **Zeros in market shares — YES**, into demand_new near the inversion.
   Sketch: 2 frames — the problem (log 0 undefined; sampling vs. structural
   zeros; selection on ξ; the drop-the-zeros and "add 1/2n" hacks and why they
   bias) and the solutions (Gandhi–Lu–Shi (QE 2023) bounds, Dubé–Hortaçsu–Joo
   (2021), Quan–Williams (RAND 2018) aggregation; practical advice).
   *(Correction: earlier drafts said "Gandhi–Ho–Nevo" — the zeros paper is
   Gandhi–Lu–Shi.)*
3. **Worked CS/CV computation — YES (approved 2026-09-01)**, in mc1 right after
   Logit Inclusive Value: (i) Small–Rosen $E[CS] = \frac{1}{\alpha}\ln\sum_k
   e^{V_k} + C$, one line off the preceding frame; (ii) worked 3-good **removal**
   example — chosen over a price change because it threads Hausman (1997)/Week 2,
   HW2, Week 8's $WTP_i(j)$, and outside-good diversion; (iii) caveat bullets:
   no income effects (constant α; random $\alpha_i$ breaks the clean form) and
   the new-goods overstatement critique with a Petrin (2002) pointer (replicated
   in Week 5). **Pedagogical note from Chris: students find the Hausman
   virtual-price construction confusing — frame the removal example explicitly
   as "what Hausman was computing, without the virtual price": drop $j$ from the
   choice set and the log-sum-exp handles it.** Applied echo in the demand_new
   case study: deferred (deck already dense).
4. **Pure characteristics + vertical model — one slide each**, resurrected in
   condensed form from `Week 4-Old/extensions.tex` (lines ~181, ~138/162) into
   demand_new's "Other Analytic Inverses" neighborhood. Framing: "these exist,
   rarely used today."
5. **IIOC "Good News / Bad News" retrospective — YES**, into the main chain:
   the Good News (BLP works; tuning answered Knittel–Metaxoglou), Bad News:
   Substitution (rank-one logit, low-dim projection), Bad News: Elasticities
   frames from `Leuven-Lectures/demand/demand_iioc.tex` become demand_new's
   closing section.
6. **Surplus function + Williams–Daly–Zachary theorem — YES (Chris, 2026-09-01):
   missing from the old McFadden/Train material.** One frame in mc1 between
   Logit Inclusive Value and the new CS/CV frames: define the social surplus
   function $G(\symbf{V}) = \E[\max_j (V_j + \varepsilon_j)]$; WDZ:
   $\sigma_j = \partial G/\partial V_j$ (convex $G$; the ARUM analogue of Roy's
   identity/Shephard's lemma); logit case: $\nabla \log\sum_k e^{V_k} =$
   softmax — note the deck's existing smooth-max/softmax frames are this theorem
   in disguise; remarks: chain rule through $V_j(p_j)$ links to the CS/CV frames
   that follow; forward pointers to Hotz–Miller (Week 11 inverts the map) and
   convex-duality treatments (Chiong–Galichon–Shum). Sequence becomes: inclusive
   value → WDZ → worked ΔCS.
7. **Search/consideration — mention + preview only.** One frame (likely
   alongside the retrospective closer): what dropping full-information/full-
   consideration does, pointers (Goeree 2008; Honka et al.), and "covered in
   depth in ECON-GA 1802 (search weeks)".

## Open questions

- Spring-course pointer in `course_intro.tex` ("take Waldinger's course") —
  pending the 1802 rotation / third syllabus.
- Webpage (chrisconlon.github.io): broken Common Ownership link (points at the
  DP notes); demand2 links n/a after retirement (it wasn't linked — nothing to
  do); update lecture list if deck set changes.

## Deliberately NOT doing

- No dynamic-games/entry promotion decision until the third syllabus arrives
  (see `other_courses/course-overlap-notes.md`).
- No consideration-set/search *content* beyond the preview frame (1802 owns it).
- Leuven decks stay frozen as historical talks.
