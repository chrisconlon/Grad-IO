# Overlap with ECON-GA 1802 (econ-dept IO II) — notes toward reorganization

*Drafted 2026-08-31. Status: **on hold** pending the third syllabus (likely Audrey
Tiew's — the S2026 description names the sequence as Jovanovic + Conlon + Tiew).*

## The two vintages

- **Spring 2025 — Daniel Waldinger.** Entry (2 wks, incl. multiple equilibria),
  dynamic games (2 wks), auctions (1 wk), market design/matching (3 wks),
  selection markets (2 wks), search/transportation in the tail.
- **Spring 2026 — Thi Mai Anh Nguyen.** Substantially redesigned around
  transactions and frictions: auctions (wks 1–2), consumer + spatial search
  (wks 3–5), contracts incl. selection markets and contract design (wks 6–8),
  long-term/relational contracts (wks 10–12, her research area), **vertical
  contracts & bargaining (wks 13–14)**.

Both vintages explicitly presume my course delivers BLP, Rust, and Hotz–Miller —
that role in the sequence is fixed.

## Coverage by vintage

| Topic | Waldinger S2025 | Nguyen S2026 |
| --- | --- | --- |
| Entry models (+ multiple equilibria) | 2 weeks | **dropped** |
| Moment inequalities (Ciliberto–Tamer, PPHI) | ~1 week (with entry) | **dropped** |
| Dynamic games | 2 weeks | **dropped** |
| Auctions | 1 week | 2 weeks |
| Market design / matching | 3 weeks | **dropped** |
| Selection markets | 2 weeks | folded into Contracts (wks 6–8) |
| Consumer + spatial search | (TBD tail) | 3 weeks (new) |
| Contract design, moral hazard | — | wks 6–8 (new) |
| Long-term/relational contracts | — | 3 weeks (new) |
| Vertical contracts & bargaining | — | 2 weeks (new) |

## Findings that matter for my reorganization

1. **The planned "flesh out Vertical Markets and Bargaining" collides head-on
   with Nguyen wks 13–14.** Her reading list is exactly what I'd teach:
   Collard-Wexler–Gowrisankaran–Lee Nash-in-Nash, CLWY, Ho–Lee, Grennan,
   Villas-Boas, Ho–Ho–Mortimer, Conlon–Mortimer (2021). Proposed carve at the
   demand/transaction boundary: I keep what leans on demand estimation and
   pricing — double marginalization, foreclosure with estimated demand (CLWY
   deck), rebates/AUD, and the unincorporated **pass-through** block (~22
   Leuven frames she doesn't touch, which connects to my merger material) —
   plus at most one overview lecture on bargaining econometrics pointing into
   1802.

2. **Entry, moment inequalities, and dynamic games are now orphaned.**
   Waldinger covered them; Nguyen dropped all three. Unless Tiew picks them up,
   students only see Berry (1992), Bresnahan–Reiss, Ciliberto–Tamer, BBL/POB in
   my Week 15/15c leftovers → argues for promoting entry + moment inequalities
   to proper slots (`entry.tex` and `moment_inequalities.tex` build again after
   the teaching_slides migration). Dynamic games would be a new build — **wait
   for Tiew's syllabus** (industry dynamics is a plausible topic for her).

3. **Smaller notes.**
   - Week 8 (WTP/Health) only lightly overlaps her selection-markets weeks: I
     come at hospitals from demand/WTP (CDS, Ho 2009), she from asymmetric
     info/contract design. Complementary; worth cross-references (Ho 2009 and
     Ho–Lee appear on both sides).
   - Auctions and search are cleanly hers — no reason for me to enter.
   - Production functions (my Week 6, currently PDF-only) appear in nobody's
     course — the planned flesh-out fills a real hole in the sequence.

## Organizational ideas worth stealing

- Nguyen's **mini-estimation project** (pick a model → derive estimator →
  simulate → Monte Carlo, four milestones, ≤5-page write-up) as an alternative
  to a fifth problem set.
- Waldinger's **referee report** assignment (matches infrastructure I already
  have in `referee-reports` guidelines).

## Next step

When the third syllabus arrives: settle the entry/dynamic-games question, then
draft a revised week-by-week Fall 2026 syllabus with the vertical trim and
entry/MI promotion penciled in.
