# DM Linear Time Series
# IPI Pharmaceutical Industry (NAF 21)
# done by Mathieu Tassel & Ralph Nader


# loading all required Packages 
pkgs <- c("insee", "dplyr", "stringr", "tseries", "forecast", "vars","ellipse")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(insee); library(dplyr); library(stringr)
library(tseries); library(forecast); library(vars)

filter <- dplyr::filter
select <- dplyr::select

# Question 1) Data retrieval and descriptive statistics

idbank_list <- get_idbank_list("IPI-2021") |>
  filter(FREQ == "M")        |>
  filter(NATURE == "INDICE")   |>
  filter(CORRECTION == "CVS-CJO")  |>
  filter(str_detect(NAF2, "^21$")) |>
  add_insee_title()

idbank_pharma <- idbank_list |> pull(idbank)

raw       <- get_insee_idbank(idbank_pharma)
ipi_pharma <- raw |>
  transmute(date = as.Date(paste0(TIME_PERIOD, "-01")),
            ipi  = as.numeric(OBS_VALUE)) |>
  arrange(date) |>
  filter(!is.na(ipi))

cat("\n Question 1) descriptive statistics \n")
cat("Period:", format(min(ipi_pharma$date), "%Y-%m"),
    "to", format(max(ipi_pharma$date), "%Y-%m"), "\n")
cat("T:", nrow(ipi_pharma), "\n")
cat("Mean:", round(mean(ipi_pharma$ipi), 2), "\n")
cat("Median:", round(median(ipi_pharma$ipi), 2), "\n")
cat("Std dev:", round(sd(ipi_pharma$ipi), 2), "\n")
cat("Min:", round(min(ipi_pharma$ipi), 2),
    "(", format(ipi_pharma$date[which.min(ipi_pharma$ipi)], "%Y-%m"), ")\n")
cat("Max:", round(max(ipi_pharma$ipi), 2),
    "(", format(ipi_pharma$date[which.max(ipi_pharma$ipi)], "%Y-%m"), ")\n")
cat("CV:", round(sd(ipi_pharma$ipi) / mean(ipi_pharma$ipi) * 100, 1), "%\n")

write.csv(ipi_pharma, "ipi_pharma_cvs_cjo.csv", row.names = FALSE)

# Question 2) Stationarisation

ipi_ts    <- ts(ipi_pharma$ipi, start = c(1990, 1), frequency = 12)
log_ipi   <- log(ipi_ts)
d_log_ipi <- diff(log_ipi)

cat("\n Question 2) Unit root tests on log(IPI) in level \n")
print(adf.test(log_ipi, alternative = "stationary"))
print(pp.test(log_ipi, alternative = "stationary"))
print(kpss.test(log_ipi, null = "Trend"))
print(kpss.test(log_ipi, null = "Level"))

cat("\n Question 2) Unit root tests on diff(log(IPI)) \n")
print(adf.test(d_log_ipi,  alternative = "stationary"))
print(pp.test(d_log_ipi,   alternative = "stationary"))
print(kpss.test(d_log_ipi, null = "Level"))


# Question 3) Graphical representation

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

plot(ipi_ts,
     main = "Raw series: IPI Pharmaceutical industry (NAF 21)",
     ylab = "Index (base 100 = 2021)", xlab = "",
     col = "#B22222", lwd = 1.5)
abline(h = 100, lty = 2, col = "grey50")
abline(v = 2020, lty = 3, col = "red")
text(2015, 30, "Covid-19 (2020)", col = "red", cex = 0.8)
grid()

plot(d_log_ipi,
     main = "Stationary series: delta log(IPI)",
     ylab = "delta log(IPI)", xlab = "",
     col = "#B22222", lwd = 1.2)
abline(h = 0,    lty = 2, col = "grey50")
abline(v = 2020, lty = 3, col = "red")
text(2015, min(d_log_ipi) * 0.85, "Covid-19 (2020)",
     col = "red", cex = 0.8)
grid()

dev.copy(pdf, "figure_q3.pdf", width = 8, height = 7); dev.off()

# Question 4) Identification

# ACF / PACF
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(as.numeric(d_log_ipi),  lag.max = 24, main = "ACF  -  delta log(IPI)")
pacf(as.numeric(d_log_ipi), lag.max = 24, main = "PACF -  delta log(IPI)")
dev.copy(pdf, "figure_acf_pacf.pdf", width = 8, height = 4); dev.off()

# Candidate models (from ACF/PACF reading)
res_01 <- arima(d_log_ipi, order = c(0, 0, 1))  # MA(1) i.e. the baseline
res_02 <- arima(d_log_ipi, order = c(0, 0, 2))  # MA(2)
res_11 <- arima(d_log_ipi, order = c(1, 0, 1))  # ARMA(1,1)
res_21 <- arima(d_log_ipi, order = c(2, 0, 1))  # ARMA(2,1)
res_31 <- arima(d_log_ipi, order = c(3, 0, 1))  # ARMA(3,1)
res_41 <- arima(d_log_ipi, order = c(4, 0, 1))  # ARMA(4,1)
res_12 <- arima(d_log_ipi, order = c(1, 0, 2))  # ARMA(1,2)
res_20 <- arima(d_log_ipi, order = c(2, 0, 0))  # AR(2)
res_10 <- arima(d_log_ipi, order = c(1, 0, 0))  # AR(1)

# AIC / BIC table
models <- list("MA(1)"     = res_01, "ARMA(1,2)" = res_12,
               "ARMA(1,1)" = res_11, "MA(2)"     = res_02,
               "ARMA(2,1)" = res_21, "ARMA(4,1)" = res_41,
               "ARMA(3,1)" = res_31, "AR(2)"     = res_20,
               "AR(1)"     = res_10)

ic_table <- data.frame(
  Model = names(models),
  AIC   = round(sapply(models, AIC), 3),
  BIC   = round(sapply(models, BIC), 3)
)
cat("\n Qurestion 4) AIC / BIC (sorted by AIC) \n")
print(ic_table[order(ic_table$AIC), ])

# Coefficient significance check
check_coef <- function(nom, mod) {
  coefs   <- mod$coef
  ses     <- sqrt(diag(mod$var.coef))
  t_stats <- coefs / ses
  p_vals  <- 2 * (1 - pnorm(abs(t_stats)))
  cat("\n====", nom, "====\n")
  print(round(data.frame(Estimate = coefs, SE = ses,
                         t_stat = t_stats, p_value = p_vals), 4))
  cat("AIC:", round(AIC(mod), 3), "| BIC:", round(BIC(mod), 3), "\n")
}

check_coef("MA(1)",     res_01)
check_coef("ARMA(1,2)", res_12)
check_coef("ARMA(1,1)", res_11)
check_coef("MA(2)",     res_02)
check_coef("ARMA(2,1)", res_21)
check_coef("ARMA(4,1)", res_41)
check_coef("ARMA(3,1)", res_31)

# Retained model
best_model <- res_01
cat("\n>>> RETAINED MODEL: MA(1) <<<\n")

# Ljung-Box on MA(1) residuals
cat("\n Question 4) Ljung-Box test on MA(1) residuals \n")
for (h in c(5, 10, 15, 20)) {
  lb  <- Box.test(residuals(best_model), lag = h,
                  type = "Ljung-Box", fitdf = 1)
  sig <- ifelse(lb$p.value < 0.05, " *", "")
  cat("lag =", formatC(h, width = 2), "| Q' =",
      round(lb$statistic, 3), "| df =", lb$parameter,
      "| p =", round(lb$p.value, 4), sig, "\n")
}

# Residual diagnostics figure (ACF + ACF squared)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(as.numeric(residuals(best_model)),
    lag.max = 24, main = "ACF - residuals")
acf(as.numeric(residuals(best_model))^2,
    lag.max = 24, main = "ACF - residuals squared")
dev.copy(pdf, "figure_residuals_MA1.pdf", width = 8, height = 3.5)
dev.off()

# Final parameter recap 
coefs   <- best_model$coef
ses     <- sqrt(diag(best_model$var.coef))
t_stats <- coefs / ses

cat("\n Question 4) MA(1) parameter recap\n")
cat("theta1    :", round(coefs["ma1"], 4),
    "| SE:", round(ses["ma1"], 4),
    "| t:", round(t_stats["ma1"], 4), "\n")
cat("mu        :", round(coefs["intercept"], 4),
    "| SE:", round(ses["intercept"], 4),
    "| t:", round(t_stats["intercept"], 4), "\n")
cat("sigma^2   :", round(best_model$sigma2, 6), "\n")
cat("AIC       :", round(AIC(best_model), 3), "\n")
cat("BIC       :", round(BIC(best_model), 3), "\n")

# Q8) Forecast and joint confidence region


theta1  <- as.numeric(-best_model$coef["ma1"])
mu      <- as.numeric(best_model$coef["intercept"])
sig2    <- as.numeric(best_model$sigma2)
eps_T   <- as.numeric(tail(residuals(best_model), 1))
Y_T     <- as.numeric(tail(log_ipi, 1))

Yhat_T1 <- Y_T + mu + theta1 * eps_T
Yhat_T2 <- Y_T + 2 * mu + theta1 * eps_T

cat("\n=== Q8  Point forecasts ===\n")
cat("Yhat_{T+1|T} :", round(Yhat_T1, 4), "-> IPI:", round(exp(Yhat_T1), 2), "\n")
cat("Yhat_{T+2|T} :", round(Yhat_T2, 4), "-> IPI:", round(exp(Yhat_T2), 2), "\n")

# Covariance matrix
Sigma <- sig2 * matrix(c(1, theta1, theta1, 1 + theta1^2), nrow = 2, byrow = TRUE)
cat("Sigma:\n"); print(round(Sigma, 6))

# Marginal 95% CIs (log scale then exponentiated)
ic1_lo <- Yhat_T1 - qnorm(0.975) * sqrt(sig2)
ic1_hi <- Yhat_T1 + qnorm(0.975) * sqrt(sig2)
ic2_lo <- Yhat_T2 - qnorm(0.975) * sqrt(sig2 * (1 + theta1^2))
ic2_hi <- Yhat_T2 + qnorm(0.975) * sqrt(sig2 * (1 + theta1^2))

cat("95% CI X_{T+1}: [", round(exp(ic1_lo), 2), ";", round(exp(ic1_hi), 2), "]\n")
cat("95% CI X_{T+2}: [", round(exp(ic2_lo), 2), ";", round(exp(ic2_hi), 2), "]\n")

Xhat_T1 <- exp(Yhat_T1);  Xhat_T2 <- exp(Yhat_T2)
Xic1_lo <- exp(ic1_lo);   Xic1_hi <- exp(ic1_hi)
Xic2_lo <- exp(ic2_lo);   Xic2_hi <- exp(ic2_hi)

# ---- Figure 1: time-series forecast ----
n_show    <- 36
ipi_show  <- tail(ipi_ts, n_show)
dates_num <- as.numeric(time(tail(ipi_ts, n_show)))
date_T1   <- max(time(ipi_ts)) + 1/12
date_T2   <- max(time(ipi_ts)) + 2/12

ylim_range <- range(c(ipi_show, Xic1_lo, Xic1_hi, Xic2_lo, Xic2_hi))

pdf("figure_forecast_q8.pdf", width = 9, height = 5)
par(mfrow = c(1, 1), mar = c(6, 5, 4, 2))

plot(dates_num, as.numeric(ipi_show),
     type = "o", pch = 19, cex = 0.6, lwd = 1.5, col = "black",
     xlim = c(min(dates_num) - 0.1, date_T2 + 0.12),
     ylim = c(ylim_range[1] - 15, ylim_range[2] + 10),
     xlab = "", ylab = "IPI (base 100 = 2021)",
     main = "IPI Pharmaceutical Industry: forecast March-April 2026",
     xaxt = "n")

axis(1, at = seq(2023, 2026.5, by = 0.5),
     labels = c("Jan 2023","Jul 2023","Jan 2024","Jul 2024",
                "Jan 2025","Jul 2025","Jan 2026","Jul 2026"),
     cex.axis = 0.75, las = 2)

abline(v = max(dates_num), lty = 2, col = "grey50")

arrows(date_T1, Xic1_lo, date_T1, Xic1_hi, angle = 90, code = 3,
       length = 0.07, col = "#FFAAAA", lwd = 2.5)
arrows(date_T2, Xic2_lo, date_T2, Xic2_hi, angle = 90, code = 3,
       length = 0.07, col = "#FFAAAA", lwd = 2.5)

lines(c(max(dates_num), date_T1, date_T2),
      c(as.numeric(tail(ipi_ts, 1)), Xhat_T1, Xhat_T2),
      col = "#B22222", lwd = 2, lty = 2)
points(c(date_T1, date_T2), c(Xhat_T1, Xhat_T2),
       pch = 17, col = "#B22222", cex = 1.8)

# T+1 label below, T+2 label above
text(date_T1, Xic1_lo - 5,
     paste0(round(Xhat_T1, 1), "\n[", round(Xic1_lo, 1), "; ", round(Xic1_hi, 1), "]"),
     col = "#FFAAAA", cex = 0.78, font = 2, adj = 0.5)
text(date_T2, Xic2_hi + 4,
     paste0(round(Xhat_T2, 1), "\n[", round(Xic2_lo, 1), "; ", round(Xic2_hi, 1), "]"),
     col = "#FFAAAA", cex = 0.78, font = 2, adj = 0.5)

legend("topleft",
       legend = c("Observed IPI", "Point forecast", "95% marginal CI"),
       col    = c("black", "#B22222", "#FFAAAA"),
       lty    = c(1, 2, NA), pch = c(19, 17, NA),
       lwd    = c(1.5, 2, 2.5), cex = 0.82, bg = "white")
grid()
dev.off()

ell     <- ellipse(Sigma, centre = c(Yhat_T1, Yhat_T2), level = 0.95)
ell_ipi <- exp(ell)

# Cartesian product of marginal CIs (rectangle)
rect_x <- c(Xic1_lo, Xic1_hi, Xic1_hi, Xic1_lo, Xic1_lo)
rect_y <- c(Xic2_lo, Xic2_lo, Xic2_hi, Xic2_hi, Xic2_lo)

pdf("figure_ellipse_q8.pdf", width = 7, height = 7)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

plot(ell_ipi[, 1], ell_ipi[, 2],
     type = "l", lwd = 2, col = "#B22222",
     xlim = range(c(ell_ipi[, 1], rect_x)) + c(-0.5, 0.5),
     ylim = range(c(ell_ipi[, 2], rect_y)) + c(-0.5, 0.5),
     xlab = expression(X[T+1] ~ "(IPI, base 100 = 2021)"),
     ylab = expression(X[T+2] ~ "(IPI, base 100 = 2021)"),
     main = "Joint 95% Confidence Region for\n(X_{T+1}, X_{T+2})")

# Marginal CI rectangle
lines(rect_x, rect_y, lty = 2, col = "grey40", lwd = 1.5)

# Point forecast
points(Xhat_T1, Xhat_T2, pch = 21, bg = "#B22222", col = "#B22222", cex = 1.8)
text(Xhat_T1 + 0.2, Xhat_T2 + 0.2,
     expression(paste("(", hat(X)[T+1], ", ", hat(X)[T+2], ")")),
     col = "#B22222", cex = 0.85, font = 2)

legend("topleft",
       legend = c("95% joint ellipse", "Cartesian product of marginal CIs", "Point forecast"),
       col    = c("steelblue", "grey40", "#B22222"),
       lty    = c(1, 2, NA), pch = c(NA, NA, 17),
       lwd    = c(2, 1.5, NA), cex = 0.82, bg = "white")
grid()
dev.off()



# Question 9) Granger causality tests

# Candidaten nb 1: IPI Chimie (NAF 20)
raw_ch <- get_insee_idbank("010767783")

ipi_chimie <- raw_ch |>
  transmute(date = as.Date(paste0(TIME_PERIOD, "-01")),
            ipi  = as.numeric(OBS_VALUE)) |>
  arrange(date) |>
  filter(!is.na(ipi))

merged <- inner_join(ipi_pharma |> rename(pharma = ipi),
                     ipi_chimie  |> rename(chimie = ipi),
                     by = "date")

d_log_pharma <- diff(log(merged$pharma))
d_log_chimie  <- diff(log(merged$chimie))
bivar <- cbind(d_log_pharma, d_log_chimie)
colnames(bivar) <- c("pharma", "chimie")

p1   <- as.integer(VARselect(bivar, lag.max = 12,
                             type = "const")$selection["AIC(n)"])
var1 <- VAR(bivar, p = p1, type = "const")
g1   <- causality(var1, cause = "chimie")
s1   <- cov(resid(var1))
rho1 <- s1[1,2] / sqrt(s1[1,1] * s1[2,2])

cat("\n Question 9) Candidate nb 1: IPI Chimie, VAR(", p1, ") \n")
print(g1$Granger); print(g1$Instant)
cat("Rho =", round(rho1, 4), "\n")

# Candidate nb 2: balance of opinion (ENQ-CONJ-ACT-IND) 
raw_enq <- get_insee_idbank("001586064")

enq_prod <- raw_enq |>
  transmute(date  = as.Date(paste0(TIME_PERIOD, "-01")),
            solde = as.numeric(OBS_VALUE)) |>
  arrange(date) |>
  filter(!is.na(solde))

merged2 <- inner_join(ipi_pharma |> rename(pharma = ipi),
                      enq_prod, by = "date")

d_log_pharma2 <- diff(log(merged2$pharma))
solde_aligned  <- merged2$solde[-1]
bivar2 <- cbind(d_log_pharma2, solde_aligned)
colnames(bivar2) <- c("pharma", "solde")

p2   <- as.integer(VARselect(bivar2, lag.max = 12,
                             type = "const")$selection["AIC(n)"])
var2 <- VAR(bivar2, p = p2, type = "const")
g2   <- causality(var2, cause = "solde")
s2   <- cov(resid(var2))
rho2 <- s2[1,2] / sqrt(s2[1,1] * s2[2,2])

cat("\n Question 9) Candidate nb 2: Balance of opinion, VAR(", p2, ") \n")
print(g2$Granger); print(g2$Instant)
cat("Rho =", round(rho2, 4), "\n")
cat("Variance reduction (1 - rho^2) =", round(1 - rho2^2, 4), "\n")

