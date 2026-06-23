# TODO: 
library(deSolve)
library(ggplot2)
library(ellipse)
library(tidyr)
library(chron)

is_positive_scalar <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0
}

data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)


data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time),
                            format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))

# Location == 1 → IE (Empire), Location == 0 → IA (AIMS)
data$Status <- ifelse(data$Location == 1, "IE", "IA")

# counts
n_IA <- sum(data$Status == "IA")
n_IE <- sum(data$Status == "IE")
cat("Total cases    :", nrow(data), "\n")
cat("AIMS (IA):", n_IA, "\n")
cat("Empire (IE):", n_IE, "\n")

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


betaMixprior <- 1.8
betaNightEprior <- 0
betaNightAprior <-1.5
#disease_params <- function(
    # Inside your SIR derivative function:
 #   gammaA = 1,     
  #  gammaS = 1,     
#    N      = 46     
#) {
 # stopifnot("beta must be a positive scalar" = is_positive_scalar(beta))
  #return(as.list(environment()))
#}

index_cases <- data[data$Infected.by == "Started", ]
I0_IE <- sum(data$Location == 1 & data$Infected.by == "Started") #sum(as.integer(index_case$Location == 1))
I0_IA <- sum(data$Location == 0 & data$Infected.by == "Started")
#I0_IA <- sum(as.integer(index_case$Location == 0))
S0_E <- sum(data$Location) - I0_IE
N0    <- 42
S0_A <- N0 - S0_E - I0_IA
N_A <- 25 
N_E <- 17

pop.ssiirr <- c(
  SE  = S0_E,
  SA = S0_A,
  IA = I0_IA,
  IE = I0_IE,
  RA  = 0,
  RE = 0
)

cat("\nInitial conditions:\n")
print(pop.ssiirr)

# assumed beta
values <- c(
  gammaA = 1,
  gammaE = 1,
  N      = N0
)

cat("\nParameter values:\n")
print(values)

R0_prior <- betaMixprior
cat("\nBasic Reproduction Number (R0):", round(R0_value, 3), "\n")

ssiirr <- function(t, y, parms) {
  with(c(as.list(y), parms), {
    beta_M  <- ifelse((t >= as.integer(t) + 0.3328) & (t <= as.integer(t) + 0.77071), betaMixprior, 0)
    beta_NA <- ifelse((t >= as.integer(t) + 0.3328) & (t <= as.integer(t) + 0.77071), 0, betaNightAprior)
    beta_NE <- ifelse((t >= as.integer(t) + 0.3328) & (t <= as.integer(t) + 0.77071), 0, betaNightEprior)
    
    lambdaA        <- (beta_M*(IA + IE)/N0 + beta_NA*IA/N_A)*SA
    lambdaE       <- (beta_M*(IA + IE)/N0 + beta_NE*IE/N_E)*SE
    new_infectionsA <- lambdaA
    new_infectionsE <- lambdaE
    
    dSAdt  <- -lambdaA
    dSEdt <- -lambdaE
    dIAdt <-  new_infectionsA - gammaA * IA
    dIEdt <-  new_infectionsE - gammaE * IE
    dRAdt  <-  gammaA * IA 
    dREdt  <-  gammaE * IE 
    
    return(list(c(dSAdt, dSEdt, dIAdt, dIEdt, dRAdt, dREdt)))
  })
}

# TODO: did we even plot confidence intervals? Do we need them?
#ssiirr_CI <- function(t, y, parms) {
 # with(c(as.list(y), parms), {
    
  #  lambda         <- beta * (IA + IS) / N
   # new_infections <- lambda * S
    
    #dSdt  <- -new_infections
  #  dIAdt <-  pA      * new_infections - gammaA * IA
   # dISdt <- (1 - pA) * new_infections - gammaS * IS
  #  dRdt  <-  gammaA * IA + gammaS * IS
   # dCIdt <-  new_infections             # cumulative incidence tracker
    
 #   return(list(c(dSdt, dIAdt, dISdt, dRdt, dCIdt)))
#  })
#}

# DETERMINISTIC MODEL — ASSUMED BETA = 1.8

time.out <- seq(0, 5, by = 0.02083)

ts.ssiirr <- data.frame(lsoda(
  y     = pop.ssiirr,
  times = time.out,
  func  = ssiirr,
  parms = values
))

cat("\nModel output at day 2:\n");  print(subset(ts.ssiirr, time == 2))
cat("\nModel output at day 5:\n"); print(subset(ts.ssiirr, time == 5))

#TODO: split data by location for this thing
daily_obs <- aggregate(Infection_number ~ Day + Status, data = data, FUN = length)
names(daily_obs)[3] <- "Count"
cat("\nObserved daily incidence by type:\n")
print(daily_obs)

# TODO: Redo compartments
# Plotting all four compartments over time
ts_long <- pivot_longer(ts.ssiirr,
                        cols      = c("SA","SE", "IA", "IE", "RA", "RE"),
                        names_to  = "Compartment",
                        values_to = "Count")

ggplot(ts_long, aes(x = time, y = Count, colour = Compartment)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(
    values = c(SA = "steelblue", SE = "magenta", IA = "orange", IE = "red", RA = "forestgreen", RE = "black"),
    labels = c(SA  = "Susceptible living at AIMS",
               SE = "Susceptible living at Empire",
               IA = "Infectious (living at AIMS)",
               IE = "Infectious (living at Empire)",
               RA  = "Recovered at AIMS", 
               RE = "Recovered at Empire")
  ) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("ssiirr Deterministic Model (assumed beta = 1.8)",
          subtitle = paste0("R0 = ", round(R0_value, 2),
                            "  |  pA = ", round(values["pA"], 2))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

# Plotting model infectious (IA + IS) vs observed daily cases
cumulative_obs <- aggregate(Infection_number ~ Day, data = data, FUN = length)

ggplot() +
  geom_line(data = ts.ssiirr,
            aes(x = time, y = IA + IS, colour = "Model (IA+IE)"),
            linewidth = 1.1) +
  geom_point(data = cumulative_obs,
             aes(x = Day, y = Infection_number, colour = "Observed daily cases"),
             size = 3) +
  scale_colour_manual(values = c("Model (IA+IE)"        = "red",
                                 "Observed daily cases" = "black")) +
  xlab("Days since index case") +
  ylab("Number of active infectious individuals") +
  ggtitle("ssiirr Model (assumed beta) vs Observed Outbreak Incidence") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

# NEGATIVE LOG-LIKELIHOOD (POISSON)
# TODO: change this to include all original populations
initialize_pop_CI <- function(N = N0, IA0 = I0_IA, IS0 = I0_IS) {
  c(S = N - IA0 - IS0, IA = IA0, IS = IS0, R = 0, CI = 0)
}

simEpidemic <- function(parms = disease_params(),
                        tseq  = seq(0, 5, by = 0.02083)) {
  init <- initialize_pop_CI(N = parms$N)
  as.data.frame(lsoda(init, tseq, ssiirr_CI, parms = parms))
}
# TODO: need to figure out how to change this to simulate per half hour
get_daily_incidence <- function(simDat, days = 0:4) {
  ci_vals   <- sapply(days, function(d) simDat$CI[which.min(abs(simDat$time - d))])
  daily_inc <- c(ci_vals[1], diff(ci_vals))
  data.frame(Day = days, predicted = daily_inc)
}

nllikelihood <- function(parms = disease_params(), obsDat = obsDat) {
  sim    <- simEpidemic(parms = parms, tseq = seq(0, max(obsDat$Day)+1, by = 0.02083))
  pred   <- get_daily_incidence(sim, days = obsDat$Day)
  lambda <- pmax(pred$predicted, 1e-10)   # avoid log(0)
  nlls   <- -dpois(round(obsDat$new_cases), lambda = lambda, log = TRUE)
  return(sum(nlls))
}

# Test: NLL at assumed beta = 1.8
cat("\nNLL at assumed beta = 1.8:", round(nllikelihood(disease_params(), obsDat), 5), "\n")

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
cat("\nSANN result — beta =", round(exp(optim.vals$par), 5), "\n")

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
cat("Fitted beta         :", round(beta.MLE, 5), "\n")
cat("R0 = beta           :", round(beta.MLE, 5), "\n")
cat("NLL at MLE          :", round(optim.vals$value, 5), "\n")
cat("Convergence (0=good):", optim.vals$convergence, "\n")

# Confidence Intervals

fisherInfMatrix <- solve(optim.vals$hessian)
se_log_beta     <- sqrt(diag(fisherInfMatrix))

ci_log_beta <- MLEfits["log_beta"] + c(-1, 1) * 1.96 * se_log_beta
ci_beta     <- exp(ci_log_beta)

cat("\n95% CI for beta / R0: [", round(ci_beta[1], 5), ",", round(ci_beta[2], 5), "]\n")

# Simulate with MLE beta and with assumed beta = 1.8 for comparison
fitParms  <- subsParms(MLEfits, disease_params())
fitDat    <- simEpidemic(parms = fitParms,      tseq = seq(0, 5, by = 0.02083))
assumeDat <- simEpidemic(parms = disease_params(), tseq = seq(0, 5, by = 0.02083))

# TODO: double check how days works
fit_daily    <- get_daily_incidence(fitDat,    days = 0:5)
assumed_daily <- get_daily_incidence(assumeDat, days = 0:5)

ggplot() +
  geom_col(data = obsDat,
           aes(x = Day, y = new_cases, fill = "Observed"),
           alpha = 0.5, width = 0.4) +
  geom_line(data = assumed_daily,
            aes(x = Day, y = predicted, colour = "Assumed beta = 1.8"),
            linewidth = 1, linetype = "dashed") +
  geom_line(data = fit_daily,
            aes(x = Day, y = predicted, colour = "MLE fit"),
            linewidth = 1.2) +
  geom_point(data = fit_daily,
             aes(x = Day, y = predicted, colour = "MLE fit"),
             size = 3) +
  scale_colour_manual(values = c("MLE fit"            = "steelblue",
                                 "Assumed beta = 1.8" = "grey40")) +
  scale_fill_manual(values = c("Observed" = "tomato")) +
  labs(
    x        = "Days since index case (15 June 2026)",
    y        = "New cases per day",
    title    = "ssiirr Model — Assumed Beta vs MLE Fit",
    subtitle = paste0("Assumed beta = 1.8  |  MLE beta = ", round(beta.MLE, 3),
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
  geom_vline(aes(xintercept = 1.8, linetype = "Assumed beta = 1.8"),
             colour = "grey40") +
  scale_linetype_manual(values = c("95% CI cutoff"    = "dashed",
                                   "MLE beta"         = "solid",
                                   "Assumed beta = 1.8" = "dotted")) +
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

