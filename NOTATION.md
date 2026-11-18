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

## Conduct deck (Week 7): labor-market applications — added 2026-11-13

Deck-local to the labor frames of `Week 7- Mergers and Conduct/conduct_new.tex` (the monopsony-vs-oligopsony primer, Scuderi–Roussille, Sharma). Label on first use.

| Object | Symbol | Note |
|---|---|---|
| Wage, employment | $w_j$, $n_j$ | the papers' notation; italic $w_j$ is not the upright cost shifters $\textrm{w}$ of this deck, nor the type weights $w_i$ |
| Markdown | $w_j = \text{MRPL}_j/\psi_j$, $\psi_j = 1 + 1/\varepsilon_j$ | Week 6's input-market wedge $\psi$; Sharma writes the markdown as $\mu = 1/\psi$, but $\mu^m$ is the conditional moment in this deck |
| Labor-supply elasticity | $\varepsilon_j = \partial\ln n_j/\partial\ln w_j$ | Week 6's $\varepsilon$ is the ex post production shock |
| Nest elasticities | $e$ (across employers within a market), $e_M < e$ (market as a whole) | the literature's $\eta$, $\theta$; $\eta$ is the markup in this deck |
| Shares | $s_j$ = wage-bill share of the local labor market | an observed data share, consistent with $s$ = observed shares |

## Single-agent dynamics (Weeks 10–11: Rust, MPEC, Hotz–Miller; PS4) — decided 2026-11-17

Applies to `rust_part1.tex`, `rust_part2.tex`, `mpec.tex` (Week 10), `hotz_miller.tex` (Week 11) and `Assignments/Assignment 4- Rust/ps4.tex`. Follows **Su and Judd (2012)** with two exceptions: transitions are $f$/$\mathbf F$, not Su–Judd's $p$ (lowercase $p$ is always a price, and prices are state variables in the dynamic-demand weeks), and the stationary distribution is $\mu$, not $\pi$ ($\pi_f$ is profits only). Rust's own symbols ($i_t$ for the decision, $\theta_1, \theta_2, \theta_3$, $\pi(x,i)$ for payoffs) are retired; a slide that reproduces one of Rust's tables says which of his symbols map to ours.

| Object | Symbol | Note |
|---|---|---|
| Agent (bus), period | $i$, $t$ | $i$ is never the decision |
| State, next state | $x_t \in \{0,\ldots,N-1\}$ (mileage bins), $x'$ | costs may use a rescaled $x/100$; say so |
| Decision | $d_t \in \{0,1\}$: $0$ = keep, $1$ = replace; general decks $d \in \mathcal D$ | Su–Judd; data column `d_it` |
| Payoff shocks | $\varepsilon_t(d)$, iid T1EV, location $0$, scale $1$ (mean $\gamma$) | |
| Flow payoff | $u(x,d;\theta)$: $u(x,0) = -c(x;\theta^c)$, $u(x,1) = -RC - c(0;\theta^c)$ | deterministic part; $\varepsilon_t(d)$ added separately |
| Parameters | $\theta = (\theta^c, RC, \theta^f)$; $c(x;\theta^c) = \theta^c_1 x + \theta^c_2 x^2$; $\theta^f_k = \Pr(\Delta = k)$ | Rust's $\theta_1$ = our $\theta^c$, his $\theta_3$ = our $\theta^f$; $\beta$ is fixed, not in $\theta$ |
| Timing (standard) | $d_t$ chosen after seeing $(x_t, \varepsilon_t)$; $x_{t+1} = (1-d_t)\,x_t + \Delta_{t+1}$ | so $f(x' \mid x,1) = f(x' \mid 0,0)$ |
| Transitions | density $f(x' \mid x,d;\theta^f)$; matrices $\mathbf F_d$ (rows = current state) | $\mathbf F_1$ has every row equal to row $0$ of $\mathbf F_0$ |
| Expected value function | $EV_\theta(x,d) = \E[V_\theta(x') \mid x,d]$; Bellman operator $EV = T_\theta(EV)$ | $EV(x,1) = EV(0,0)$, so one function $EV(x) \equiv EV(x,0)$ on the grid |
| Choice-specific value | $v_\theta(x,d) = u(x,d;\theta) + \beta\, EV_\theta(x,d)$ | replaces $\tilde V$ |
| Ex ante (integrated) value | $V(x) = \ln \sum_d \exp v(x,d)$ | **Euler's $\gamma$ is dropped** from $V$ and $EV$ (it adds $\gamma/(1-\beta)$ and leaves CCPs unchanged); keep it, and say so, only where $\E[\varepsilon \mid d,x] = \gamma - \ln P(d \mid x)$ itself appears (forward simulation) |
| CCPs | model $P(d \mid x;\theta)$; first stage $\hat P(d \mid x)$; vector $\mathbf P(x)$ | capital $P$ is a probability here; $\sigma$ stays the static-demand choice probability |
| Hotz–Miller, Arcidiacono–Miller | $\phi_{dd'}(\mathbf P(x)) = v(x,d) - v(x,d')$; $\psi_d(\mathbf P(x)) = V(x) - v(x,d)$ | logit: $\ln P(d\mid x) - \ln P(d'\mid x)$ and $-\ln P(d \mid x)$ (no $\gamma$, per the row above) |
| Normalizations (KSS) | "replace normalization" $u(x,1) = 0$; "keep normalization" $u(x,0) = 0$ | the first is a level shift in Rust because $u(x,1)$ is constant in $x$; the second is not |
| Cost of never replacing | $C = (\mathbf I - \beta \mathbf F_0)^{-1} c$ | keep normalization gives $\tilde u(x,1) = -RC + C(x) - C(0)$ |
| Long run | stationary distribution $\mu(x)$ of the controlled chain; engine demand $\bar d = \sum_x \mu(x) P(1 \mid x)$ | replacements per bus per period |
| Counterfactuals | engine-price change $\tau$: $u^\tau(x,1) = u(x,1) - \tau$; new transitions $\mathbf F'_d$ | $\Delta$ is the mileage increment, so not $\Delta$ for the price change |

### Euler-equation CCPs (Week 11 `eccp.tex`, replaces Paul Scott's deck) — decided 2026-11-17

Extends the single-agent block above; everything there still applies (scale 1, no $\gamma$, $\psi_d = -\ln P(d \mid x)$). Chris decided $\eta$ and $d$; the rest is the proposal he did not object to and may be revised once the deck is built. A slide that reproduces one of Scott's tables says which of his symbols map to ours ($j \to d$, $p_{jkt} \to P_t(d \mid k)$, $\nu \to \varepsilon$, $\varepsilon^{CV} \to \eta$, $\pi \to u$).

| Object | Symbol | Note |
|---|---|---|
| Decision (land use) | $d_{it} \in \{0,1\}$: $1$ = crops (the renewal action, like "replace"), $0$ = other | Scott's $j$ |
| Field state | $k_{it} \in \{0,\ldots,\bar k\}$, years since last in crops; $k_{i,t+1} = \kappa(d_{it}, k_{it})$ with $\kappa(1,k) = 0$, $\kappa(0,k) = \min\{k+1, \bar k\}$ | deterministic, so $\kappa$ rather than $f$ |
| Market state | $\omega_t$, the farmers' information set; state $x_{it} = (k_{it}, \omega_t)$ | **deck-local**: $\omega$ is productivity in Week 6 |
| Flow payoff | $u_t(k,d) = \alpha_0(k,d) + \alpha_R R_t(d) + \xi_t(k,d)$, plus $\varepsilon_{it}(d)$ | $R$ observed returns; $\xi$ unobserved aggregate shock; normalization $\alpha_0(k,0) = 0$ for all $k$ |
| CCPs, values | $P_t(d \mid k)$; $v_t(k,d)$, $V_t(k)$, all carrying $t$ because $\omega_t$ makes the problem non-stationary | |
| Expectational error | $\eta_{t+1}(k') = V_{t+1}(k') - \E_t[V_{t+1}(k')]$, $\E[\eta_{t+1} \mid \omega_t] = 0$ | **deck-local**: $\eta$ is a demand elasticity or markup elsewhere; $\varepsilon$ stays the logit shock |
