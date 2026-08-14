# =============================================================================
# Extreme value analysis of storm surge from a NOAA CO-OPS tide gauge
#
# v2 - supersedes fit_surge_gev.R. Run this file top to bottom.
#      source("fit_surge_gev_v2.R")   or paste in the section order below.
#
# CHANGE FROM v1: the annual mean residual is now removed before extracting
# annual maxima. NOAA predictions are referenced to the 1983-2001 National Tidal
# Datum Epoch, centred near 1992, so observed water in recent years sits above
# the prediction by the sea level rise accumulated since then. Portland shows
# +0.12 m in a January 2020 test. Left in, that secular offset would (a) be
# double-counted when sea level rise is added separately in the flood model, and
# (b) make the stationarity test fire on datum drift rather than storm climate.
#
# The removed offset is kept, because its slope is itself a local sea level rise
# estimate and serves as an independent check on the whole pipeline.
#
# STATUS: the download and column-matching logic has been verified against the
# live API. The fitting and plotting sections have NOT been executed by their
# author. Run section by section and read the output at each stage.
# =============================================================================

# ---- configuration ----------------------------------------------------------

STATION      <- "8418150"      # Portland, Maine - verified against the live API
YEAR_FROM    <- 1980
YEAR_TO      <- 2025
DATUM        <- "MSL"
CACHE_DIR    <- "noaa_cache"
MIN_COVERAGE <- 0.80           # reject years with less than this fraction of hours

RETURN_PERIODS <- c(1.5, 2, 3, 5, 10, 15, 25, 50, 100, 200, 350, 500)

dir.create(CACHE_DIR, showWarnings = FALSE)
cat("Working directory:", getwd(), "\n")


# ---- section 1: download and build the annual record ------------------------
# CO-OPS caps hourly requests at one year, so this loops year by year and caches
# to disk. Re-running is cheap. Delete noaa_cache to force a refresh.

coops_url <- function(product, year, interval = NULL) {
  u <- paste0("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter",
              "?product=", product, "&application=research",
              "&begin_date=", year, "0101", "&end_date=", year, "1231",
              "&datum=", DATUM, "&station=", STATION,
              "&time_zone=GMT&units=metric&format=csv")
  if (!is.null(interval)) u <- paste0(u, "&interval=", interval)
  u
}

fetch_year <- function(product, year, interval = NULL) {
  f <- file.path(CACHE_DIR, sprintf("%s_%s_%d.csv", STATION, product, year))
  if (!file.exists(f)) {
    ok <- tryCatch({ download.file(coops_url(product, year, interval), f,
                                   quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (!ok) { if (file.exists(f)) unlink(f); return(NULL) }
  }
  d <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  # A failed API call returns a CSV containing an error string, not data.
  if (is.null(d) || nrow(d) < 100) { unlink(f); return(NULL) }
  d
}

pick_col <- function(df, pattern) {
  i <- grep(pattern, names(df), ignore.case = TRUE)
  if (length(i) == 0) stop("no column matching '", pattern, "' in: ",
                           paste(names(df), collapse = ", "))
  names(df)[i[1]]
}

cat("\nDownloading", YEAR_FROM, "to", YEAR_TO, "for station", STATION, "\n")
records <- list()

for (yr in YEAR_FROM:YEAR_TO) {
  obs  <- fetch_year("hourly_height", yr)
  pred <- fetch_year("predictions",   yr, interval = "h")
  if (is.null(obs) || is.null(pred)) {
    cat(sprintf("  %d  skipped (no data)\n", yr)); next
  }

  o <- data.frame(t   = trimws(obs[[pick_col(obs, "date")]]),
                  obs = as.numeric(obs[[pick_col(obs, "water")]]))
  p <- data.frame(t    = trimws(pred[[pick_col(pred, "date")]]),
                  pred = as.numeric(pred[[pick_col(pred, "pred")]]))
  m <- merge(o, p, by = "t")
  m <- m[is.finite(m$obs) & is.finite(m$pred), ]
  if (nrow(m) < 1000) { cat(sprintf("  %d  skipped (only %d hours)\n",
                                    yr, nrow(m))); next }

  m$residual <- m$obs - m$pred
  offset     <- mean(m$residual)        # datum drift since epoch + seasonal mean
  m$surge    <- m$residual - offset     # meteorological surge for this year

  records[[as.character(yr)]] <- list(
    year          = yr,
    coverage      = nrow(m) / 8760,
    max_surge     = max(m$surge),
    mean_residual = offset,
    max_residual  = max(m$residual),
    max_observed  = max(m$obs))

  cat(sprintf("  %d  n=%5d  coverage=%4.0f%%  offset=%+.3f m  max surge=%5.2f m\n",
              yr, nrow(m), 100 * nrow(m) / 8760, offset, max(m$surge)))
}

if (length(records) == 0) stop("No data retrieved. Check STATION and connectivity.")
cat(sprintf("\n%d years retrieved.\n", length(records)))


# ---- section 2: annual maxima and sea level rise check ----------------------

am <- do.call(rbind, lapply(records, function(r)
  data.frame(year          = r$year,
             coverage      = r$coverage,
             surge         = r$max_surge,
             mean_residual = r$mean_residual,
             raw_residual  = r$max_residual,
             observed      = r$max_observed)))
am <- am[order(am$year), ]

n_all <- nrow(am)
am <- am[am$coverage >= MIN_COVERAGE, ]
cat(sprintf("\n%d years, %d retained after the %.0f%% coverage filter\n",
            n_all, nrow(am), 100 * MIN_COVERAGE))
cat(sprintf("Annual maximum surge: mean %.2f m, sd %.2f m, range %.2f-%.2f m\n",
            mean(am$surge), sd(am$surge), min(am$surge), max(am$surge)))

if (nrow(am) < 20)
  warning("Fewer than 20 annual maxima. The GEV shape estimate will be unstable.")

# Independent validation: the slope of the removed offset is a local sea level
# rise estimate relative to the 1983-2001 tidal epoch. Portland's published rate
# is roughly 2 mm/yr. A result near that confirms the pipeline is sound.
slr <- lm(mean_residual ~ year, data = am)
slr_rate <- unname(coef(slr)["year"]) * 1000
slr_p    <- summary(slr)$coefficients["year", "Pr(>|t|)"]
cat(sprintf("\nSEA LEVEL CHECK: mean annual residual trend %+.2f mm/yr (p = %.4f)\n",
            slr_rate, slr_p))
cat("  Compare against the published rate for this gauge at\n")
cat("  tidesandcurrents.noaa.gov/sltrends/. A close match validates the pipeline.\n")


# ---- section 3: stationarity of the surge itself ----------------------------
# With the datum drift removed, any remaining trend is a change in storm climate
# rather than in mean sea level. This does not prove stationarity - it only tests
# for a detectable linear trend.

trend    <- lm(surge ~ year, data = am)
trend_mm <- unname(coef(trend)["year"]) * 1000
trend_p  <- summary(trend)$coefficients["year", "Pr(>|t|)"]

cat(sprintf("\nTrend in annual maximum surge: %+.2f mm/yr (p = %.3f)\n",
            trend_mm, trend_p))
if (trend_p < 0.05) {
  cat("  Detectable. Options, in rough order of preference:\n")
  cat("   (a) raise YEAR_FROM to a recent window and refit\n")
  cat("   (b) fit a non-stationary GEV with a time-varying location\n")
  cat("   (c) proceed stationary and state the limitation explicitly\n")
} else {
  cat("  No trend detected at the 5% level. A stationary GEV is defensible.\n")
  cat("  Sea level rise still matters for TOTAL water level and must be added\n")
  cat("  separately in the risk model - it is not in this surge distribution.\n")
}


# ---- section 4: GEV fit by maximum likelihood -------------------------------
# Coles parameterisation: xi > 0 heavy (Frechet), xi = 0 Gumbel,
# xi < 0 bounded (Weibull).
# WARNING: scipy.stats.genextreme uses c = -xi. Do not mix conventions.

gev_nll <- function(p, x) {
  mu <- p[1]; sigma <- exp(p[2]); xi <- p[3]
  z <- (x - mu) / sigma
  if (abs(xi) < 1e-8) return(sum(log(sigma) + z + exp(-z)))
  t <- 1 + xi * z
  if (any(t <= 0)) return(1e10)
  sum(log(sigma) + (1 + 1/xi) * log(t) + t^(-1/xi))
}

fit_gev <- function(x) {
  s0  <- sd(x) * sqrt(6) / pi                 # method-of-moments Gumbel start
  mu0 <- mean(x) - 0.5772 * s0
  o <- optim(c(mu0, log(s0), 0.05), gev_nll, x = x, method = "Nelder-Mead",
             control = list(maxit = 5000, reltol = 1e-10))
  list(mu = o$par[1], sigma = exp(o$par[2]), xi = o$par[3],
       nll = o$value, convergence = o$convergence)
}

return_level <- function(T, mu, sigma, xi) {
  y <- -log(1 - 1/T)
  if (abs(xi) < 1e-8) return(mu - sigma * log(y))
  mu + (sigma / xi) * (y^(-xi) - 1)
}

fit <- fit_gev(am$surge)
if (fit$convergence != 0) warning("optim did not converge cleanly")

cat(sprintf("\nGEV fit (n = %d years)\n  location mu    %+.4f m\n  scale    sigma %+.4f m\n  shape    xi    %+.4f\n",
            nrow(am), fit$mu, fit$sigma, fit$xi))
cat(if (fit$xi >  0.05) "  Heavy tail: rare surges exceed what a Gumbel predicts.\n"
    else if (fit$xi < -0.05) "  Bounded tail: the fit implies a physical upper limit.\n"
    else "  Shape near zero; Gumbel is a reasonable simplification here.\n")

if (requireNamespace("extRemes", quietly = TRUE)) {
  ev <- extRemes::fevd(am$surge, type = "GEV")
  cat("  extRemes cross-check:", sprintf("%+.4f ", ev$results$par), "\n")
} else {
  cat("  (install.packages('extRemes') for an independent cross-check)\n")
}


# ---- section 5: return levels with bootstrap confidence intervals -----------
# Parametric bootstrap captures PARAMETER uncertainty only. It assumes the GEV
# family is correct, so the intervals are narrower than the true uncertainty.

rgev <- function(n, mu, sigma, xi) {
  u <- runif(n)
  if (abs(xi) < 1e-8) return(mu - sigma * log(-log(u)))
  mu + (sigma / xi) * ((-log(u))^(-xi) - 1)
}

set.seed(42)
NBOOT <- 1000
boot <- matrix(NA_real_, NBOOT, length(RETURN_PERIODS))
for (b in seq_len(NBOOT)) {
  fb <- tryCatch(fit_gev(rgev(nrow(am), fit$mu, fit$sigma, fit$xi)),
                 error = function(e) NULL)
  if (is.null(fb) || fb$convergence != 0) next
  boot[b, ] <- sapply(RETURN_PERIODS, return_level, fb$mu, fb$sigma, fb$xi)
}

rl <- data.frame(
  return_period_yr = RETURN_PERIODS,
  surge_m  = sapply(RETURN_PERIODS, return_level, fit$mu, fit$sigma, fit$xi),
  lower_95 = apply(boot, 2, quantile, 0.025, na.rm = TRUE),
  upper_95 = apply(boot, 2, quantile, 0.975, na.rm = TRUE),
  gumbel_m = sapply(RETURN_PERIODS, return_level, fit$mu, fit$sigma, 0))

cat("\nReturn levels (meteorological surge, m)\n")
print(rl, row.names = FALSE, digits = 3)


# ---- section 6: diagnostics -------------------------------------------------
# A fit that has not been checked against the data is not a result.

png("gev_diagnostics.png", width = 1400, height = 1100, res = 130)
par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))

plot(am$year, am$surge, type = "b", pch = 19, cex = 0.7,
     xlab = "Year", ylab = "Annual max surge (m)",
     main = sprintf("Annual maxima, station %s", STATION))
abline(trend, col = "firebrick", lwd = 2, lty = 2)

n <- nrow(am); xs <- sort(am$surge); emp <- (seq_len(n) - 0.5) / n
theo <- sapply(emp, function(q) return_level(1/(1-q), fit$mu, fit$sigma, fit$xi))
plot(theo, xs, pch = 19, cex = 0.7, xlab = "Model quantile (m)",
     ylab = "Empirical quantile (m)", main = "QQ plot")
abline(0, 1, col = "firebrick", lwd = 2)

plot(am$year, am$mean_residual, type = "b", pch = 19, cex = 0.7,
     xlab = "Year", ylab = "Mean annual residual (m)",
     main = sprintf("Removed offset: %+.2f mm/yr", slr_rate))
abline(slr, col = "steelblue", lwd = 2)

plot(rl$return_period_yr, rl$surge_m, log = "x", type = "b", pch = 19,
     ylim = range(rl$lower_95, rl$upper_95, am$surge),
     xlab = "Return period (years)", ylab = "Surge (m)",
     main = "Return level plot")
polygon(c(rl$return_period_yr, rev(rl$return_period_yr)),
        c(rl$lower_95, rev(rl$upper_95)),
        col = adjustcolor("steelblue", 0.20), border = NA)
lines(rl$return_period_yr, rl$gumbel_m, col = "darkgreen", lty = 3, lwd = 2)
points(1/(1-emp), xs, pch = 4, col = "grey30", cex = 0.8)
legend("topleft", bty = "n", cex = 0.8,
       legend = c("GEV fit", "95% bootstrap", "Gumbel (xi=0)", "observed"),
       col = c("black", "steelblue", "darkgreen", "grey30"),
       lty = c(1, 1, 3, NA), pch = c(19, NA, NA, 4), lwd = c(2, 8, 2, NA))

dev.off()
cat("\nDiagnostics written to gev_diagnostics.png\n")


# ---- section 7: export for the Python flood model ---------------------------

write.csv(rl, "surge_return_levels.csv", row.names = FALSE)
write.csv(am, "annual_maxima.csv", row.names = FALSE)

cat('\n--- paste into damage.py, replacing the illustrative values ---\n')
cat(sprintf('GEV_LOC   = %.4f   # NOAA station %s, %d-%d (n = %d years)\n',
            fit$mu, STATION, min(am$year), max(am$year), nrow(am)))
cat(sprintf('GEV_SCALE = %.4f\n', fit$sigma))
cat(sprintf('GEV_SHAPE = %.4f   # Coles convention; scipy genextreme uses c = -xi\n',
            fit$xi))
cat('GEV_IS_FITTED = True\n')
cat('---------------------------------------------------------------\n')

writeLines(c(
  sprintf("station: %s", STATION),
  sprintf("datum: %s", DATUM),
  sprintf("years: %d-%d (n = %d after coverage filter)",
          min(am$year), max(am$year), nrow(am)),
  sprintf("GEV location: %.4f m", fit$mu),
  sprintf("GEV scale: %.4f m", fit$sigma),
  sprintf("GEV shape: %.4f", fit$xi),
  sprintf("surge trend: %+.2f mm/yr (p = %.4f)", trend_mm, trend_p),
  sprintf("sea level trend from removed offset: %+.2f mm/yr (p = %.4f)",
          slr_rate, slr_p),
  "NOTE: this is METEOROLOGICAL SURGE only. Mean sea level rise has been",
  "removed and must be added separately in the risk model."
), "gev_fit_summary.txt")

cat("\nWritten: surge_return_levels.csv, annual_maxima.csv, gev_fit_summary.txt\n")
