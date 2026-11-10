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

## Production functions (Week 6) — decided 2026-11-09

Applies to `Week 6- Production Functions/production.tex` and the production problem set. Lowercase letters are logs of uppercase levels, **except prices**, which stay lowercase levels as everywhere else in the course (logs written $\ln p$, $\ln w^X$).

| Object | Symbol | Note |
|---|---|---|
| Producer, period | $j$, $t$ | not $i$ (investment is $i_{jt}$; consumers are $i$ elsewhere) |
| Output, inputs | $Y,K,L,M$ / $y,k,l,m$; vector $\mathbf x_{jt}$; generic variable input $X$; variable inputs $\mathcal V$ | $q$ stays a demand quantity |
| Production function | $y = f(\mathbf x) + \omega + \varepsilon$ | Hicks-neutral $\omega$; $\varepsilon$ ex post shock / measurement error |
| Markov law | $\omega_{jt} = g(\omega_{jt-1}) + \xi_{jt}$; AR(1) persistence $\rho$ | **$\xi$ here is the productivity innovation, not BLP's demand unobservable**; demand shifters are $\delta_{jt}$. Say so on the slide where $\xi$ first appears. |
| Control function | $\omega = h_t(\mathbf o_{jt})$; step-1 observables $\mathbf o_{jt}$; extra input-demand shifters $\mathbf d_{jt}$; first stage $\Phi_t$; instruments $\mathbf z_{jt}$ | DLW's $\mathbf z$ is renamed $\mathbf d$ so $\mathbf z$ means instruments, as in the conduct deck |
| Prices | output price $p_{jt}$, input prices $w^X_{jt}$ (levels) | |
| Elasticities and shares | $\theta^X = \partial f/\partial x$; revenue share $\alpha^X = w^X X/(pY)$; $\alpha^{\mathcal V} = VC/(pY)$ | not $s$ (observed market shares) |
| Scale, markup | $\nu = \sum_{X\in\mathcal V}\theta^X$; $mc$; $\mu = p/mc$ | $\mu$ matches the pass-through deck |
| Deck-local reuse (label on first use) | $\rho$ = AR(1) persistence; $\nu$ = short-run scale elasticity; $\sigma$ = elasticity of substitution; $\eta > 1$ = CES demand elasticity, $\mu = \eta/(\eta-1)$ | course-wide meanings: nesting parameter, RC draws, choice probability |
| Controlled Markov process; multi-product firms | $\omega_{jt} = g(\omega_{jt-1}, a_{jt-1}) + \xi_{jt}$ with $a$ = export $e$, R&D $\mathit{rd}$, quota $\mathit{qr}$; product $n$ of firm $j$, input share $\phi_{jnt}$ | added 2026-11-09 (round 3); $g$ and $\rho$ were taken, so not DGKP's product index $j$ or share $\rho$ |
