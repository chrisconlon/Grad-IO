# Notation conventions (Grad-IO)

Following Berry & Haile (2014). Adopted 2026-08 for the demand material; new decks should conform.

## Choice probabilities and shares

| Object | Symbol | Notes |
|---|---|---|
| Individual choice probability (model) | `\sigma_{ijt}` (or `\sigma_{ij}` in static settings) | Conditionals: `\sigma_{ij\|g}`, nest-level `\sigma_{ig}` |
| Market-level choice probability / share function | `\sigma_{jt}(\cdot)`, e.g. `\sigma_j(\delta_t,\theta)` | Inversion: `\sigma_j^{-1}(\cdot)` |
| Observed market share (data) | `s_{jt}` **only** | Vector: `\symbf{s}_t`. Empirical frequencies may wear hats: `\hat{s}_{ij}` |
| Estimation equation | `\sigma_j(\delta_t,\theta) = s_{jt}` | The model matches the data |

**Retired symbols** (do not use): `\mathfrak{s}_{jt}`, `\mathcal{S}_{jt}`/`\calS` for shares, `\tilde{s}_{jt}`, `P_{ij}`, `S_{ij}`, and `\pi` as a choice probability.

## Parameters and other objects

- `\rho` — nesting parameter (Berry 1994's $\sigma$; matches PyBLP). Where a deck deliberately contrasts parameterizations, keep McFadden's `\lambda` with $\rho = 1-\lambda$ and a remark that Berry (1994) writes $\sigma$.
- `\pi_f` — profits, and nothing else. Type/mixture weights are `w_i` / `w_{it}` (so `s_{jt} = \sum_i w_{it}\,\sigma_{ijt}`).
- `\Pi`, `\Sigma` — demographic-interaction and random-coefficient parameter matrices (PyBLP convention). Random coefficients use **capital** `\Sigma` whenever possible: write `\beta_i = \beta + \Sigma \nu_i` (with $\Sigma$ diagonal if need be) and refer to elements of $\Sigma$, rather than lowercase `\sigma_x \nu_i` / `\sigma \nu_i`. Never use bare `\sigma_0` for an RC value at the truth (it now reads as the outside-good choice probability) — write `\Sigma_0` etc.
- Prices are always lowercase `p_{jt}`; capital `P` is never a price.
- `\Delta_{(j,k)}(\symbf{p}) = -\partial \sigma_j / \partial p_k` — the demand-derivative matrix (not `\Omega`).
- Diversion ratios — two distinct **aggregate** objects, kept notationally separate:
  - `D_{jk}(\symbf{p})` — **small-price-change** (derivative-based) diversion: $D_{jk} = \frac{\partial q_k/\partial p_j}{|\partial q_j/\partial p_j|}$.
  - `D_{j \rightarrow k}` — **second-choice** (choice-set-removal) diversion.

  Classify aggregate uses by the underlying experiment (marginal price change vs. removal of $j$). **Individual-level** diversion ratios are always plain `D_{jk,i}` (no arrow): $D_{jk,i} = \frac{\sigma_{ik}}{1-\sigma_{ij}}$ does not depend on the intervention — only the aggregation weights do.
- `a` — the T1EV/error **scale parameter** (formerly $\sigma$): $\sigma_{ij} = e^{V_{ij}/a}/\sum_k e^{V_{ik}/a}$, only $\beta/a$ identified.
- Probability statements: use `\Pr(\cdot)` or `\prob{\cdot}` (both render blackboard $\mathbb{P}$ via `teaching_slides.sty`), not literal `\mathbb{P}`.
