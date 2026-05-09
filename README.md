# IPI Pharmaceutical Industry — Linear Time Series Analysis

> **Linear Time Series Assignment** &nbsp;|&nbsp; *Mathieu Tassel & Ralph Nader*

![Language](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Data](https://img.shields.io/badge/Data-INSEE%20BDM-003189?style=flat-square)
![Model](https://img.shields.io/badge/Model-ARIMA(0%2C1%2C1)-2e7d32?style=flat-square)
![Period](https://img.shields.io/badge/Period-Jan%201990%20–%20Feb%202026-8e44ad?style=flat-square)

---

A rigorous time series study of the French **Industrial Production Index (IPI)** for the pharmaceutical sector (NAF rev. 2, division 21). Starting from raw monthly INSEE data spanning **434 observations** (Jan 1990 – Feb 2026), we walk through the full econometric pipeline: stationarisation, ARMA identification and validation, two-step-ahead forecasting with joint confidence regions, and Granger-causality tests against early-available indicators.

---

## Table of Contents

1. [Data](#1-data)
2. [Stationarisation](#2-stationarisation)
3. [ARMA Identification & Estimation](#3-arma-identification--estimation)
4. [Forecasting (March – April 2026)](#4-forecasting-march--april-2026)
5. [Granger Causality](#5-granger-causality)
6. [Repository Structure](#6-repository-structure)
7. [How to Run](#7-how-to-run)

---

## 1. Data

**Source:** INSEE BDM — IPI 2021, series filtered on NAF division 21, monthly frequency, CVS-CJO adjustment, base 100 in 2021.

Division 21 covers basic pharmaceutical substances and finished medicinal products. The series is already seasonally and working-day adjusted, making it ready for direct modelling.

| Statistic | Value |
|---|---|
| Period | Jan. 1990 – Feb. 2026 |
| Observations (T) | 434 |
| Mean | 69.36 |
| Median | 70.92 |
| Std. deviation | 25.09 |
| Min | 25.11 (Feb. 1990) |
| Max | 133.20 (Aug. 2024) |

The long-run upward trajectory is driven by population ageing and rising public health expenditure. The most visible structural break is the **Covid-19 shock (2020)**, followed by a vaccine-driven rebound. Since the series is strictly positive with growing variance, we apply the log transformation `Yₜ = log Xₜ` — this simultaneously stabilises variance and converts multiplicative growth into additive increments.

---

## 2. Stationarisation

Unit root tests on `Yₜ = log(IPIₜ)` **in level**:

| Test | Specification | Statistic | p-value |
|---|---|---|---|
| ADF | Const. + trend | −1.890 | 0.624 |
| PP | Const. + trend | −31.734 | < 0.01 |
| KPSS | Trend stat. (H₀) | 1.179 | < 0.01 |
| KPSS | Level stat. (H₀) | 6.749 | < 0.01 |

The PP rejection can be attributed to its tendency to over-reject in the presence of structural breaks (visible around 2020). ADF and KPSS agree: **Yₜ ∉ I(0)**.

After applying a first difference `ΔYₜ = log(Xₜ / Xₜ₋₁)`, all three tests confirm stationarity:

| Test | Specification | Statistic | p-value |
|---|---|---|---|
| ADF | Constant | −9.063 | < 0.01 |
| PP | Constant | −545.6 | < 0.01 |
| KPSS | Level stat. (H₀) | 0.109 | > 0.10 |

**Conclusion:** `Yₜ ~ I(1)`, so `d = 1` is sufficient. We model `ΔYₜ` directly.

---

## 3. ARMA Identification & Estimation

### 3.1 ACF / PACF Reading

The ACF cuts off sharply after lag 1 (`ρ̂(1) ≈ −0.68`) while the PACF decays geometrically over the first four lags — the canonical signature of a **MA(1)** process. We nonetheless estimate a full set of candidate models to confirm.

### 3.2 Model Selection

| Model | AIC | BIC | Elimination rationale |
|---|---|---|---|
| **MA(1)** | **−1503.820** | **−1491.608** | **All coefficients significant — retained** |
| ARMA(1,2) | −1502.518 | −1482.165 | All coef. sig. but BIC +9.4 vs MA(1) |
| ARMA(1,1) | −1501.924 | −1485.642 | φ̂₁ insig. (t = −0.32) |
| MA(2) | −1501.923 | −1485.640 | θ̂₂ insig. (t = 0.32) |
| ARMA(2,1) | −1499.933 | −1479.580 | φ̂₁, φ̂₂ insig. |
| ARMA(4,1) | −1499.562 | −1471.067 | φ̂₁, φ̂₂, φ̂₃ insig. |
| ARMA(3,1) | −1498.144 | −1473.719 | φ̂₁, φ̂₂, φ̂₃ insig. |
| AR(2) | −1480.477 | −1464.194 | Dominated on both criteria |
| AR(1) | −1446.609 | −1434.397 | Dominated on both criteria |

Both AIC and BIC unanimously select **MA(1)**. Every ARMA(p,1) with p ≥ 1 is dominated on both criteria, confirming the AR component is not needed.

### 3.3 Estimated Model

The retained **ARIMA(0,1,1)** (equivalently, EWMA) on the log-series:

```
(1 − B) Yₜ = 0.0033 + (1 + 0.6841 B) εₜ,    εₜ ~ⁱⁱᵈ WN(0, 0.001789)
```

| Parameter | Estimate | Std. Error | t-stat | p-value |
|---|---|---|---|---|
| θ̂₁ | −0.6841 | 0.0348 | −19.64 | < 0.001 |
| μ̂ | 0.0033 | 0.0006 | 5.15 | < 0.001 |
| σ̂² | 0.001789 | — | — | — |

The drift `μ̂ = 0.0033` implies an average monthly growth rate of **0.33%**.

### 3.4 Residual Validation

Ljung-Box test on MA(1) residuals:

| Lags H | Q′ | df | p-value |
|---|---|---|---|
| 5 | 4.526 | 4 | 0.339 |
| 10 | 17.520 | 9 | 0.041 |
| 15 | 32.145 | 14 | 0.004 |
| 20 | 35.092 | 19 | 0.014 |

H₀ is not rejected at lag 5 (p = 0.339). Rejections at longer horizons reflect **conditional heteroskedasticity** (ARCH effects) driven by the Covid-19 outlier, not residual autocorrelation in the mean — as confirmed by the ACF of squared residuals. Since ARIMA targets the conditional mean, this is outside the model's scope and the MA(1) is considered valid.

---

## 4. Forecasting (March – April 2026)

### 4.1 Point Forecasts & Marginal Confidence Intervals

From the ARIMA(0,1,1) with T = Feb. 2026:

| Horizon | Point forecast (IPI) | 95% Marginal CI |
|---|---|---|
| T+1 (March 2026) | **120.7** | [111.1 ; 131.1] |
| T+2 (April 2026) | **121.1** | [109.5 ; 133.9] |

The CI for T+2 is wider than for T+1, reflecting `Var(eₜ₊₂) = σ̂²(1 + θ̂₁²) = 0.002626 > σ̂² = 0.001789`. Both forecasts lie well below the August 2024 peak (133.2), consistent with the moderate downward correction observed since mid-2024. The MA(1) memory is exhausted after one period, so `X̂ₜ₊₂ ≈ X̂ₜ₊₁`.

### 4.2 Joint Confidence Region

The joint 95% confidence ellipse for `(Xₜ₊₁, Xₜ₊₂)` is derived from the bivariate Gaussian structure of the prediction errors:

```
Σ̂ = σ̂² [ 1       θ̂₁      ]  =  [ 0.001789   0.001224 ]
         [ θ̂₁   1 + θ̂₁²  ]     [ 0.001224   0.002626 ]
```

The joint ellipse is **strictly smaller** than the Cartesian product of the two marginal intervals because `Cov(eₜ₊₁, eₜ₊₂) = σ̂²θ̂₁ = 0.001224 > 0`: a positive surprise in March 2026 propagates into April through the MA(1) structure.

---

## 5. Granger Causality

We test two series published before the official IPI release, using a bivariate VAR framework on `(Δlog Xₜ, Yₜ)ᵀ`. Two conditions must hold jointly for `Yₜ₊₁` to improve the forecast of `Xₜ₊₁`:

1. **Granger causality** (F-test): does the *past* of Y help predict X beyond X's own past?
2. **Instantaneous causality** (χ²(1)): does `Yₜ₊₁` itself carry additional information?

| Candidate | Model | Test | Stat. | p-value | ρ̂ |
|---|---|---|---|---|---|
| IPI Chimie (NAF 20) | VAR(7) | Granger (F) | 1.033 | 0.406 | −0.001 |
| IPI Chimie (NAF 20) | VAR(7) | Instantaneous (χ²) | 0.000 | 0.985 | — |
| Balance of opinion | VAR(10) | Granger (F) | 1.793 | 0.058 | 0.103 |
| Balance of opinion | VAR(10) | Instantaneous (χ²) | 4.422 | **0.035** | — |

**IPI Chimie (NAF 20):** Neither condition holds (p = 0.406; p = 0.985, ρ̂ ≈ 0). Pharmaceutical output is **demand-driven** (prescriptions, health policy) rather than input-constrained by the chemical industry.

**Balance of opinion (INSEE ENQ-CONJ-ACT-IND):** Granger causality is marginal (p = 0.058), but instantaneous causality is significant (p = 0.035, ρ̂ = 0.103): condition 2 holds. The survey and the IPI capture the same monthly production cycle from different measurement angles, yielding contemporaneous co-movement. Knowing the survey outcome for T+1 yields a marginal variance reduction of `1 − 0.103² ≈ 98.9%` — a modest but statistically valid improvement.

---

## 6. Repository Structure

```
.
├── src/
│   └── mtassel_rnader.R        # Full analysis script (data retrieval to Granger tests)
├── data/
│   ├── ipi_pharma_cvs_cjo.csv  # IPI Pharma monthly series (Jan 1990 – Feb 2026)
│   ├── ipi_chimie.csv          # IPI Chimie NAF 20 (Granger candidate 1)
│   └── balance_of_opinion.csv  # INSEE business survey balance (Granger candidate 2)
├── img/
│   ├── figure_q3.pdf           # Raw and stationary series
│   ├── figure_acf_pacf.pdf     # ACF and PACF of Δlog(IPI)
│   ├── figure_residuals_MA1.pdf# MA(1) residual diagnostics
│   ├── figure_forecast_q8.pdf  # Forecast for March–April 2026
│   ├── figure_ellipse_q8.pdf   # Joint 95% confidence ellipse
│   ├── figure_granger_q9.pdf   # Granger test: IPI Chimie
│   └── figure_granger_enq_q9.pdf # Granger test: Balance of opinion
└── DM_times_series.pdf         # Full written report
```

---

## 7. How to Run

The entire pipeline — data retrieval, stationarisation, model selection, forecasting, and causality tests — is contained in a single self-contained R script.

**Prerequisites:** R ≥ 4.1, with internet access to query the INSEE API.

```r
# All required packages are installed automatically on first run:
# insee, dplyr, stringr, tseries, forecast, vars, ellipse

source("src/mtassel_rnader.R")
```

The script will:
1. Fetch the IPI Pharma series live from the INSEE BDM API
2. Run unit root tests and print descriptive statistics
3. Identify, estimate, and validate the MA(1) model
4. Produce and export all figures to `img/`
5. Compute forecasts and joint confidence regions
6. Run Granger causality tests for both candidate series

---

## Authors

| | |
|---|---|
| **Mathieu Tassel** | M2 Economics |
| **Ralph Nader** | M2 Economics |

*Linear Time Series — Graduate Econometrics Assignment*
