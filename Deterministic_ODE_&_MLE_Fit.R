
library(deSolve)
library(ggplot2)
library(ellipse)
library(tidyr)

is_positive_scalar <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0
}

setwd("C:/Users/Ennie Matlhanya/MMEDGit/MMF")
data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)


data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time),
                            format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))

# Location == 1 → IS (symptomatic), Location == 0 → IA (asymptomatic)
data$Status <- ifelse(data$Location == 1, "IS", "IA")

# counts
n_IA <- sum(data$Status == "IA")
n_IS <- sum(data$Status == "IS")
cat("Total cases    :", nrow(data), "\n")
cat("Asymptomatic (IA):", n_IA, "\n")
cat("Symptomatic  (IS):", n_IS, "\n")

#daily incidence
daily_total <- aggregate(Infection_number ~ Day, data = data, FUN = length)
names(daily_total)[2] <- "new_cases"

daily_IA <- aggregate(Infection_number ~ Day,
                      data = subset(data, Status == "IA"), FUN = length)
daily_IS <- aggregate(Infection_number ~ Day,
                      data = subset(data, Status == "IS"), FUN = length)
names(daily_IA)[2] <- "new_IA"
names(daily_IS)[2] <- "new_IS"

all_days <- data.frame(Day = 0:max(data$Day))
obsDat   <- merge(all_days, daily_total, by = "Day", all.x = TRUE)
obsDat   <- merge(obsDat,   daily_IA,    by = "Day", all.x = TRUE)
obsDat   <- merge(obsDat,   daily_IS,    by = "Day", all.x = TRUE)
obsDat[is.na(obsDat)] <- 0

cat("\nObserved daily incidence:\n")
print(obsDat)

disease_params <- function(
    beta   = 0.6,   
    pA     = 0.4,   
    gammaA = 1,     
    gammaS = 1,     
    N      = 34     
) {
  stopifnot("beta must be a positive scalar" = is_positive_scalar(beta))
  return(as.list(environment()))
}

index_case <- data[data$Infected.by == "Started", ][1, ]
I0_IS <- as.integer(index_case$Location == 1)
I0_IA <- as.integer(index_case$Location == 0)
N0    <- nrow(data)

pop.siar <- c(
  S  = N0 - I0_IA - I0_IS,
  IA = I0_IA,
  IS = I0_IS,
  R  = 0
)

cat("\nInitial conditions:\n")
print(pop.siar)

# assumed beta
values <- c(
  beta   = 0.6,
  pA     = 0.4,
  gammaA = 1,
  gammaS = 1,
  N      = N0
)

cat("\nParameter values:\n")
print(values)

# R0 = beta when gammaA = gammaS = 1
R0_value <- with(as.list(values), {
  beta * (pA / gammaA + (1 - pA) / gammaS)
})
cat("\nBasic Reproduction Number (R0):", round(R0_value, 3), "\n")

siar <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    
    lambda         <- beta * (IA + IS) / N
    new_infections <- lambda * S
    
    dSdt  <- -new_infections
    dIAdt <-  pA      * new_infections - gammaA * IA
    dISdt <- (1 - pA) * new_infections - gammaS * IS
    dRdt  <-  gammaA * IA + gammaS * IS
    
    return(list(c(dSdt, dIAdt, dISdt, dRdt)))
  })
}

siar_CI <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    
    lambda         <- beta * (IA + IS) / N
    new_infections <- lambda * S
    
    dSdt  <- -new_infections
    dIAdt <-  pA      * new_infections - gammaA * IA
    dISdt <- (1 - pA) * new_infections - gammaS * IS
    dRdt  <-  gammaA * IA + gammaS * IS
    dCIdt <-  new_infections             # cumulative incidence tracker
    
    return(list(c(dSdt, dIAdt, dISdt, dRdt, dCIdt)))
  })
}

# DETERMINISTIC MODEL — ASSUMED BETA = 0.6

time.out <- seq(0, 10, by = 0.1)

ts.siar <- data.frame(lsoda(
  y     = pop.siar,
  times = time.out,
  func  = siar,
  parms = values
))

cat("\nModel output at day 4:\n");  print(subset(ts.siar, time == 4))
cat("\nModel output at day 10:\n"); print(subset(ts.siar, time == 10))


daily_obs <- aggregate(Infection_number ~ Day + Status, data = data, FUN = length)
names(daily_obs)[3] <- "Count"
cat("\nObserved daily incidence by type:\n")
print(daily_obs)


# Plotting all four compartments over time
ts_long <- pivot_longer(ts.siar,
                        cols      = c("S", "IA", "IS", "R"),
                        names_to  = "Compartment",
                        values_to = "Count")

ggplot(ts_long, aes(x = time, y = Count, colour = Compartment)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(
    values = c(S = "steelblue", IA = "orange", IS = "red", R = "forestgreen"),
    labels = c(S  = "Susceptible",
               IA = "Infectious (Asymptomatic)",
               IS = "Infectious (Symptomatic)",
               R  = "Recovered")
  ) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("SIAR Deterministic Model (assumed beta = 0.6)",
          subtitle = paste0("R0 = ", round(R0_value, 2),
                            "  |  pA = ", round(values["pA"], 2))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Plotting model infectious (IA + IS) vs observed daily cases
cumulative_obs <- aggregate(Infection_number ~ Day, data = data, FUN = length)

ggplot() +
  geom_line(data = ts.siar,
            aes(x = time, y = IA + IS, colour = "Model (IA+IS)"),
            linewidth = 1.1) +
  geom_point(data = cumulative_obs,
             aes(x = Day, y = Infection_number, colour = "Observed daily cases"),
             size = 3) +
  scale_colour_manual(values = c("Model (IA+IS)"        = "red",
                                 "Observed daily cases" = "black")) +
  xlab("Days since index case") +
  ylab("Number of active infectious individuals") +
  ggtitle("SIAR Model (assumed beta) vs Observed Outbreak Incidence") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

# NEGATIVE LOG-LIKELIHOOD (POISSON)

initialize_pop_CI <- function(N = N0, IA0 = I0_IA, IS0 = I0_IS) {
  c(S = N - IA0 - IS0, IA = IA0, IS = IS0, R = 0, CI = 0)
}

simEpidemic <- function(parms = disease_params(),
                        tseq  = seq(0, 4, by = 0.01)) {
  init <- initialize_pop_CI(N = parms$N)
  as.data.frame(lsoda(init, tseq, siar_CI, parms = parms))
}

get_daily_incidence <- function(simDat, days = 0:4) {
  ci_vals   <- sapply(days, function(d) simDat$CI[which.min(abs(simDat$time - d))])
  daily_inc <- c(ci_vals[1], diff(ci_vals))
  data.frame(Day = days, predicted = daily_inc)
}

nllikelihood <- function(parms = disease_params(), obsDat = obsDat) {
  sim    <- simEpidemic(parms = parms, tseq = seq(0, max(obsDat$Day), by = 0.01))
  pred   <- get_daily_incidence(sim, days = obsDat$Day)
  lambda <- pmax(pred$predicted, 1e-10)   # avoid log(0)
  nlls   <- -dpois(round(obsDat$new_cases), lambda = lambda, log = TRUE)
  return(sum(nlls))
}

# Test: NLL at assumed beta = 0.6
cat("\nNLL at assumed beta = 0.6:", round(nllikelihood(disease_params(), obsDat), 4), "\n")

# Fitting beta on log scale 
subsParms <- function(fit.params, fixed.params = disease_params()) {
  within(fixed.params, {
    loggedParms   <- names(fit.params)[grepl("log_", names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl("log_", names(fit.params))]
    for (nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for (nm in loggedParms)   assign(gsub("log_", "", nm),
                                     exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}

objFXN <- function(fit.params, fixed.params = disease_params(), obsDat = obsDat) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat)
}

# Sanity check at a starting guess
cat("NLL at beta = 2 (guess):",
    round(objFXN(c(log_beta = log(2)), disease_params(), obsDat), 4), "\n")

#  Fitting beta via MLE

# (stochastic)
init.pars  <- c(log_beta = log(1))

optim.vals <- optim(
  par          = init.pars,
  fn           = objFXN,
  fixed.params = disease_params(),
  obsDat       = obsDat,
  control      = list(trace = 3, maxit = 150),
  method       = "SANN"
)
cat("\nSANN result — beta =", round(exp(optim.vals$par), 4), "\n")

# deterministic
optim.vals <- optim(
  par          = optim.vals$par,
  fn           = objFXN,
  fixed.params = disease_params(),
  obsDat       = obsDat,
  control      = list(trace = 3, maxit = 800, reltol = 1e-7),
  method       = "Nelder-Mead",
  hessian      = TRUE
)

MLEfits  <- optim.vals$par
beta.MLE <- exp(MLEfits["log_beta"])

cat("\n--- MLE Results ---\n")
cat("Fitted beta         :", round(beta.MLE, 4), "\n")
cat("R0 = beta           :", round(beta.MLE, 4), "\n")
cat("NLL at MLE          :", round(optim.vals$value, 4), "\n")
cat("Convergence (0=good):", optim.vals$convergence, "\n")

# Confidence Intervals

fisherInfMatrix <- solve(optim.vals$hessian)
se_log_beta     <- sqrt(diag(fisherInfMatrix))

ci_log_beta <- MLEfits["log_beta"] + c(-1, 1) * 1.96 * se_log_beta
ci_beta     <- exp(ci_log_beta)

cat("\n95% CI for beta / R0: [", round(ci_beta[1], 4), ",", round(ci_beta[2], 4), "]\n")

# Simulate with MLE beta and with assumed beta = 0.6 for comparison
fitParms  <- subsParms(MLEfits, disease_params())
fitDat    <- simEpidemic(parms = fitParms,      tseq = seq(0, 4, by = 0.01))
assumeDat <- simEpidemic(parms = disease_params(), tseq = seq(0, 4, by = 0.01))

fit_daily    <- get_daily_incidence(fitDat,    days = 0:4)
assumed_daily <- get_daily_incidence(assumeDat, days = 0:4)

ggplot() +
  geom_col(data = obsDat,
           aes(x = Day, y = new_cases, fill = "Observed"),
           alpha = 0.5, width = 0.4) +
  geom_line(data = assumed_daily,
            aes(x = Day, y = predicted, colour = "Assumed beta = 0.6"),
            linewidth = 1, linetype = "dashed") +
  geom_line(data = fit_daily,
            aes(x = Day, y = predicted, colour = "MLE fit"),
            linewidth = 1.2) +
  geom_point(data = fit_daily,
             aes(x = Day, y = predicted, colour = "MLE fit"),
             size = 3) +
  scale_colour_manual(values = c("MLE fit"            = "steelblue",
                                 "Assumed beta = 0.6" = "grey40")) +
  scale_fill_manual(values = c("Observed" = "tomato")) +
  labs(
    x        = "Days since index case (15 June 2026)",
    y        = "New cases per day",
    title    = "SIAR Model — Assumed Beta vs MLE Fit",
    subtitle = paste0("Assumed beta = 0.6  |  MLE beta = ", round(beta.MLE, 3),
                      "  |  R0 = ",         round(beta.MLE, 3),
                      "  |  95% CI: [",     round(ci_beta[1], 3),
                      ", ",                 round(ci_beta[2], 3), "]"),
    colour   = NULL, fill = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# Likelihood profile over beta

beta.seq    <- exp(seq(log(0.1), log(10), length.out = 100))
nll.seq     <- sapply(beta.seq, function(b) {
  objFXN(c(log_beta = log(b)), disease_params(), obsDat = obsDat)
})
conf.cutoff <- optim.vals$value + qchisq(0.95, df = 1) / 2
profile_df  <- data.frame(beta = beta.seq, nll = nll.seq)

ggplot(profile_df, aes(x = beta, y = nll)) +
  geom_line(linewidth = 1.1, colour = "steelblue") +
  geom_hline(aes(yintercept = conf.cutoff, linetype = "95% CI cutoff"),
             colour = "red") +
  geom_vline(aes(xintercept = beta.MLE, linetype = "MLE beta"),
             colour = "black") +
  geom_vline(aes(xintercept = 0.6, linetype = "Assumed beta = 0.6"),
             colour = "grey40") +
  scale_linetype_manual(values = c("95% CI cutoff"    = "dashed",
                                   "MLE beta"         = "solid",
                                   "Assumed beta = 0.6" = "dotted")) +
  scale_x_log10() +
  labs(
    x        = expression(beta ~ "(log scale)"),
    y        = "Negative Log-Likelihood",
    title    = "Likelihood Profile for beta",
    subtitle = paste0("MLE = ", round(beta.MLE, 3),
                      "  |  95% CI: [", round(ci_beta[1], 3),
                      ", ",             round(ci_beta[2], 3), "]"),
    linetype = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

