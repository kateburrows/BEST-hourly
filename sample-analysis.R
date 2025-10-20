# ============================================================================
# Hourly Temperature and Psychiatric Emergencies: Sample Code
# ============================================================================

# Required packages
library(dplyr)
library(lubridate)
library(dlnm)
library(survival)
library(splines)
library(ggplot2)

# Load data (replace with real data)
# df_sample should contain:
# - case: binary indicator (1 = case, 0 = control)
# - temp_C.0 through temp_C.23: hourly temperatures for 24 hours before encounter
# - holid: federal holiday indicator
# - strata: matching strata for case-crossover design
# - rh_max: daily maximum relative humidity (for sensitivity analysis)

df_sample <- readRDS("df_sample.rds")

# Create temperature matrix (24 hourly temperatures)
temp_matrix <- as.matrix(df_sample[, paste0("temp_C.", 0:23)])

# ============================================================================
# FIGURE 1: Cumulative (0-24 hour) Exposure-Response Curve
# ============================================================================

# Create crossbasis with distributed lag non-linear model
# Natural spline (df=3) for temperature, quadratic polynomial for lag structure
cb_24lag <- crossbasis(temp_matrix, 
                       lag = c(0, 23),
                       argvar = list(fun = "ns", df = 3),
                       arglag = list(fun = "poly", degree = 2))

# Fit conditional logistic regression
model_24lag <- clogit(case ~ cb_24lag + holid + strata(strata), 
                      data = df_sample)

# Calculate minimum morbidity temperature (MMT)
perct <- seq(0.01, 0.99, by = 0.01)
temp_percent <- quantile(df_sample$temp_C.0, perct, na.rm = TRUE)

# Initial prediction to find MMT
pred_initial <- crosspred(cb_24lag, model_24lag, 
                         at = temp_percent, 
                         cen = median(df_sample$temp_C.0, na.rm = TRUE))

pmin <- perct[which.min(pred_initial$allfit)]
mmt <- quantile(df_sample$temp_C.0, pmin, na.rm = TRUE)

# Final prediction centered at MMT
pred_fig1 <- crosspred(cb_24lag, model_24lag,
                       at = temp_percent,
                       cen = mmt)

# Plot Figure 1
plot(pred_fig1, "overall",
     xlab = "Temperature (°C)",
     ylab = "Relative Risk",
     main = "Cumulative 24-Hour Temperature Effect",
     ylim = c(0.9, 1.25))

# Extract key estimates for Table S1
extreme_temp <- quantile(df_sample$temp_C.0, 0.99, na.rm = TRUE)
temps_of_interest <- c(10, 15, 20, 25, 30, extreme_temp)
pred_extreme <- crosspred(cb_24lag, model_24lag,
                         at = extreme_temp,
                         cen = mmt)

cat("\nExtreme heat (99th percentile) vs MMT:\n")
cat("OR:", round(pred_extreme$allRRfit[1], 2), 
    "95% CI:", round(pred_extreme$allRRlow[1], 2), "-", 
    round(pred_extreme$allRRhigh[1], 2), "\n")


# ============================================================================
# FIGURE 2: Lag-Specific Effects (Slices at 0, 6, 12, 18, 24 hours)
# ============================================================================

# Plot lag slices
plot(pred_fig1, "slices", 
     lag = c(0, 6, 12, 18, 24),
     xlab = "Temperature (°C)",
     ylab = "Relative Risk",
     main = "Lag-Specific Temperature Effects")

# For supplement: Extract all 24 lag-specific estimates
# (corresponds to Table S3 and Figure S1)


# ============================================================================
# FIGURE 3: Stratified by Daily Temperature Variability
# ============================================================================

# Calculate daily temperature range
temp_range <- apply(temp_matrix, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
df_sample$temp_range <- temp_range

# Create tertiles of variability
df_sample$range_tertile <- cut(temp_range,
                               breaks = quantile(temp_range, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                               labels = c("Low", "Mid", "High"))

# Model for each tertile
# LOW VARIABILITY
df_low <- subset(df_sample, range_tertile == "Low")
temp_matrix_low <- as.matrix(df_low[, paste0("temp_C.", 0:23)])
cb_low <- crossbasis(temp_matrix_low, 
                    lag = c(0, 23),
                    argvar = list(fun = "ns", df = 3),
                    arglag = list(fun = "poly", degree = 2))
model_low <- clogit(case ~ cb_low + holid + strata(strata), data = df_low)
pred_low <- crosspred(cb_low, model_low, at = temp_percent, cen = mmt)

# MEDIUM VARIABILITY
df_mid <- subset(df_sample, range_tertile == "Mid")
temp_matrix_mid <- as.matrix(df_mid[, paste0("temp_C.", 0:23)])
cb_mid <- crossbasis(temp_matrix_mid, 
                    lag = c(0, 23),
                    argvar = list(fun = "ns", df = 3),
                    arglag = list(fun = "poly", degree = 2))
model_mid <- clogit(case ~ cb_mid + holid + strata(strata), data = df_mid)
pred_mid <- crosspred(cb_mid, model_mid, at = temp_percent, cen = mmt)

# HIGH VARIABILITY
df_high <- subset(df_sample, range_tertile == "High")
temp_matrix_high <- as.matrix(df_high[, paste0("temp_C.", 0:23)])
cb_high <- crossbasis(temp_matrix_high, 
                     lag = c(0, 23),
                     argvar = list(fun = "ns", df = 3),
                     arglag = list(fun = "poly", degree = 2))
model_high <- clogit(case ~ cb_high + holid + strata(strata), data = df_high)
pred_high <- crosspred(cb_high, model_high, at = temp_percent, cen = mmt)

# Plot all three
par(mfrow = c(1, 3))
plot(pred_low, "overall", xlab = "Temperature (°C)", ylab = "Relative Risk",
     main = "Low Variability (<7.72°C)", col = "blue")
plot(pred_mid, "overall", xlab = "Temperature (°C)", ylab = "Relative Risk",
     main = "Mid Variability (7.72-10.59°C)", col = "darkgreen")
plot(pred_high, "overall", xlab = "Temperature (°C)", ylab = "Relative Risk",
     main = "High Variability (>10.59°C)", col = "red")
par(mfrow = c(1, 1))


# ============================================================================
# FIGURE 4: Consecutive Hours of Extreme Heat
# ============================================================================

# Calculate consecutive hours of extreme heat (≥99th percentile)
extreme_temp <- quantile(df_sample$temp_C.0, 0.99, na.rm = TRUE)

consecutive_heat <- apply(temp_matrix >= extreme_temp, 1, function(x) {
  if(all(!x)) return(0)  # No extreme heat
  max(rle(x)$lengths[rle(x)$values])  # Maximum consecutive TRUE values
})

# Create crossbasis for consecutive hours
cb_consecutive <- crossbasis(consecutive_heat, 
                            lag = c(0, 0), 
                            argvar = list(fun = "poly", degree = 2))

# Fit model
model_consecutive <- clogit(case ~ cb_consecutive + holid + strata(strata), 
                           data = df_sample)

# Prediction centered at 0 consecutive hours
pred_consecutive <- crosspred(cb_consecutive, model_consecutive, cen = 0)

# Plot
plot(pred_consecutive, "overall",
     xlab = "Consecutive hours of extreme heat",
     ylab = "Relative Risk",
     main = "Effect of Sustained Extreme Heat")

# Extract specific estimates
cat("\nConsecutive hours of extreme heat:\n")
cat("1 hour OR:", round(pred_consecutive$allRRfit[2], 2), "\n")
cat("10 hours OR:", round(pred_consecutive$allRRfit[11], 2), "\n")
cat("14 hours OR:", round(pred_consecutive$allRRfit[15], 2), "\n")


# ============================================================================
# SENSITIVITY ANALYSIS: Adjusted for Relative Humidity (Table S2)
# ============================================================================

# Refit main model with relative humidity adjustment
model_rh <- clogit(case ~ cb_24lag + holid + rh_max + strata(strata), 
                   data = df_sample)

pred_rh <- crosspred(cb_24lag, model_rh,
                    at = temp_percent,
                    cen = mmt)

pred_extreme_rh <- crosspred(cb_24lag, model_rh,
                            at = extreme_temp,
                            cen = mmt)

cat("\nWith RH adjustment - Extreme heat vs MMT:\n")
cat("OR:", round(pred_extreme_rh$allRRfit[1], 2), 
    "95% CI:", round(pred_extreme_rh$allRRlow[1], 2), "-", 
    round(pred_extreme_rh$allRRhigh[1], 2), "\n")

# ============================================================================
# Session Info for Reproducibility
# ============================================================================
sessionInfo()
