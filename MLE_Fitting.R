library(deSolve)
library(ggplot2)
library(tidyr)

data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)

data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time), format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))
data$Status   <- ifelse(data$Location == 1, "IE", "IA")


## Split observed infections into DAY and NIGHT periods


day_start <- 0.3328
day_end   <- 0.77071

origin <- as.POSIXct(paste(min(data$Date), "00:00:00"))

data$t_inf <- as.numeric(difftime(data$DateTime, origin, units = "days"))
data$frac  <- data$t_inf %% 1

data$Period <- ifelse(
  data$frac >= day_start & data$frac <= day_end,
  "Day",
  "Night"
)

# IMPORTANT:
# Started cases are initial conditions, not new infections to fit.
fit_data <- subset(data, Infected.by != "Started")

# Assign each infection to a day/night fitting interval
fit_data$PeriodDay <- ifelse(
  fit_data$Period == "Day",
  floor(fit_data$t_inf),
  ifelse(
    fit_data$frac > day_end,
    floor(fit_data$t_inf),      # night after evening
    floor(fit_data$t_inf) - 1   # night before morning
  )
)

fit_data <- subset(fit_data, PeriodDay >= 0)

# Create all possible day/night periods, including zero-count periods
max_t <- max(fit_data$t_inf)
max_period_day <- floor(max_t)

obs_grid <- do.call(rbind, lapply(0:max_period_day, function(d) {
  data.frame(
    PeriodDay = d,
    Period    = c("Day", "Night"),
    bin_start = c(d + day_start, d + day_end),
    bin_end   = c(d + day_end, d + 1 + day_start)
  )
}))

# Keep only periods that start before the last observed infection
obs_grid <- subset(obs_grid, bin_start <= max_t)

obs_counts <- aggregate(
  Infection_number ~ PeriodDay + Period,
  data = fit_data,
  FUN = length
)

names(obs_counts)[names(obs_counts) == "Infection_number"] <- "new_cases"

obsPeriod <- merge(
  obs_grid,
  obs_counts,
  by = c("PeriodDay", "Period"),
  all.x = TRUE
)

obsPeriod$new_cases[is.na(obsPeriod$new_cases)] <- 0
obsPeriod <- obsPeriod[order(obsPeriod$bin_start), ]

############################################################
## Observed day/night counts split by AIMS vs Empire
############################################################

fit_data$Population <- ifelse(fit_data$Location == 0, "AIMS", "Empire")

obs_grid_loc <- merge(
  obs_grid,
  data.frame(Population = c("AIMS", "Empire"))
)

obs_counts_loc <- aggregate(
  Infection_number ~ PeriodDay + Period + Population,
  data = fit_data,
  FUN = length
)

names(obs_counts_loc)[names(obs_counts_loc) == "Infection_number"] <- "new_cases"

obsPeriodLoc <- merge(
  obs_grid_loc,
  obs_counts_loc,
  by = c("PeriodDay", "Period", "Population"),
  all.x = TRUE
)

obsPeriodLoc$new_cases[is.na(obsPeriodLoc$new_cases)] <- 0
obsPeriodLoc <- obsPeriodLoc[order(obsPeriodLoc$Population, obsPeriodLoc$bin_start), ]

print(obsPeriodLoc)

PRIOR_betaM      <- 4.5
PRIOR_betaNightA <- 1.5   
PRIOR_betaNightE <- 0     

I0_IE <- sum(data$Location == 1 & data$Infected.by == "Started")  
I0_IA <- sum(data$Location == 0 & data$Infected.by == "Started")  
N0    <- 42
N_A   <- 25
N_E   <- 17
S0_E  <- N_E - I0_IE   
S0_A  <- N_A - I0_IA   

stopifnot("Initial conditions don't sum to N0" = (S0_A + S0_E + I0_IA + I0_IE) == N0)

initialize_pop_CI <- function() {
  c(SE = S0_E, SA = S0_A,
    IA = I0_IA, IE = I0_IE,
    RA = 0, RE = 0,
    CI = 0,
    CIA = 0,   # cumulative infections in AIMS
    CIE = 0)   # cumulative infections in Empire
}

ssiirr_CI <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    beta_M  <- ifelse(is_daytime, betaM,      0)
    beta_NA <- ifelse(is_daytime, 0,          betaNightA)
    beta_NE <- ifelse(is_daytime, 0,          PRIOR_betaNightE)
    
    inf_mix <- beta_M  * (IA + IE) / N0
    inf_NA  <- beta_NA * IA / N_A
    inf_NE  <- beta_NE * IE / N_E
    
    # New infections by population
    new_A <- (inf_mix + inf_NA) * SA
    new_E <- (inf_mix + inf_NE) * SE
    
    dSEdt <- -new_E
    dSAdt <- -new_A
    
    dIAdt <-  new_A - gammaA * IA
    dIEdt <-  new_E - gammaE * IE
    
    dRAdt <- gammaA * IA
    dREdt <- gammaE * IE
    
    # Cumulative incidence
    dCIdt  <- new_A + new_E
    dCIAdt <- new_A
    dCIEdt <- new_E
    
    return(list(c(dSEdt, dSAdt,
                  dIAdt, dIEdt,
                  dRAdt, dREdt,
                  dCIdt, dCIAdt, dCIEdt)))
  })
}
simEpidemic <- function(betaM, betaNightA,
                        tseq = seq(0, max(obsPeriod$bin_end), by = 0.02083)) {
  init  <- initialize_pop_CI()
  parms <- c(
    betaM      = betaM,
    betaNightA = betaNightA,
    gammaA     = 1,
    gammaE     = 1,
    N0         = N0,
    N_A        = N_A,
    N_E        = N_E
  )
  as.data.frame(lsoda(init, tseq, ssiirr_CI, parms = parms))
}

get_period_incidence <- function(simDat, obsPeriod) {
  
  get_CI <- function(tt) {
    approx(simDat$time, simDat$CI, xout = tt, rule = 2)$y
  }
  
  CI_start <- sapply(obsPeriod$bin_start, get_CI)
  CI_end   <- sapply(obsPeriod$bin_end, get_CI)
  
  predicted <- CI_end - CI_start
  
  data.frame(
    obsPeriod,
    predicted = pmax(predicted, 1e-10)
  )
}

nllikelihood <- function(betaM, betaNightA) {
  sim  <- simEpidemic(betaM, betaNightA)
  pred <- get_period_incidence(sim, obsPeriod)
  
  nlls <- -dpois(
    x      = round(pred$new_cases),
    lambda = pred$predicted,
    log    = TRUE
  )
  
  return(sum(nlls))
}

cat("\nNLL at priors (beta_M=1.8, betaNightA=1.5):",
    round(nllikelihood(PRIOR_betaM, PRIOR_betaNightA), 5), "\n")

objFXN <- function(log_pars) {
  betaM      <- exp(log_pars[1])
  betaNightA <- exp(log_pars[2])
  nllikelihood(betaM, betaNightA)
}

set.seed(123)  

init.pars <- c(log(PRIOR_betaM), log(PRIOR_betaNightA))

optim.vals <- optim(
  par     = init.pars,
  fn      = objFXN,
  control = list(trace = 0, maxit = 200), 
  method  = "SANN"
)
cat("\nSANN result — beta_M =", round(exp(optim.vals$par[1]), 5),
    "| betaNightA =", round(exp(optim.vals$par[2]), 5), "\n")

optim.vals <- optim(
  par     = optim.vals$par,
  fn      = objFXN,
  control = list(trace = 0, maxit = 1000, reltol = 1e-7), 
  method  = "Nelder-Mead",
  hessian = TRUE
)

betaM.MLE      <- exp(optim.vals$par[1])
betaNightA.MLE <- exp(optim.vals$par[2])

cat("\n--- MLE Results ---\n")
cat("Prior beta_M          :", PRIOR_betaM,      "(unchanged)\n")
cat("Prior betaNightA      :", PRIOR_betaNightA, "(unchanged)\n")
cat("Fitted beta_M (MLE)   :", round(betaM.MLE, 5), "\n")
cat("Fitted betaNightA(MLE):", round(betaNightA.MLE, 5), "\n")
cat("NLL at MLE            :", round(optim.vals$value, 5), "\n")
cat("Convergence (0=good)  :", optim.vals$convergence, "\n")

fisherInfMatrix <- solve(optim.vals$hessian)
se_log_pars     <- sqrt(diag(fisherInfMatrix))

ci_betaM      <- exp(optim.vals$par[1] + c(-1, 1) * 1.96 * se_log_pars[1])
ci_betaNightA <- exp(optim.vals$par[2] + c(-1, 1) * 1.96 * se_log_pars[2])

cat("\n95% CI for beta_M     : [", round(ci_betaM[1], 5), ",", round(ci_betaM[2], 5), "]\n")
cat("95% CI for betaNightA : [", round(ci_betaNightA[1], 5), ",", round(ci_betaNightA[2], 5), "]\n")

fit_sim   <- simEpidemic(betaM.MLE,   betaNightA.MLE)
prior_sim <- simEpidemic(PRIOR_betaM, PRIOR_betaNightA)

## Visualization 

library(ggplot2)
library(tidyr)

# After fitting:
fit_sim <- simEpidemic(betaM.MLE, betaNightA.MLE)

# Get predicted incidence for each day/night bin
plot_dat <- get_period_incidence(fit_sim, obsPeriod)

# Make nicer labels
plot_dat$Period <- factor(plot_dat$Period, levels = c("Day", "Night"))

plot_dat$interval_label <- paste0(
  "Day ", plot_dat$PeriodDay, "\n", plot_dat$Period
)

############################################################
## Alternating Day/Night plot with prior + fitted lines
############################################################

library(ggplot2)

# Simulations
fit_sim   <- simEpidemic(betaM.MLE, betaNightA.MLE)
prior_sim <- simEpidemic(PRIOR_betaM, PRIOR_betaNightA)

# Predicted incidence by fitting interval
fit_pred   <- get_period_incidence(fit_sim, obsPeriod)
prior_pred <- get_period_incidence(prior_sim, obsPeriod)

# Build plotting data
plot_dat <- obsPeriod
plot_dat$pred_fit   <- fit_pred$predicted
plot_dat$pred_prior <- prior_pred$predicted

# Order intervals by actual time
plot_dat <- plot_dat[order(plot_dat$bin_start), ]
plot_dat$x_id <- seq_len(nrow(plot_dat))

# Nice alternating labels
plot_dat$interval_label <- paste0("Day ", plot_dat$PeriodDay, "\n", plot_dat$Period)

# Plot
############################################################
## Alternating Day/Night plot with prior + fitted lines
############################################################

library(ggplot2)

# Simulations
fit_sim   <- simEpidemic(betaM.MLE, betaNightA.MLE)
prior_sim <- simEpidemic(PRIOR_betaM, PRIOR_betaNightA)

# Predicted incidence by fitting interval
fit_pred   <- get_period_incidence(fit_sim, obsPeriod)
prior_pred <- get_period_incidence(prior_sim, obsPeriod)

# Build plotting data
plot_dat <- obsPeriod
plot_dat$pred_fit   <- fit_pred$predicted
plot_dat$pred_prior <- prior_pred$predicted

# Order intervals by actual time
plot_dat <- plot_dat[order(plot_dat$bin_start), ]
plot_dat$x_id <- seq_len(nrow(plot_dat))

# Nice alternating labels
plot_dat$interval_label <- paste0("Day ", plot_dat$PeriodDay, "\n", plot_dat$Period)

# Plot
############################################################
## Alternating Day/Night plot: observed + prior + fitted
############################################################

library(ggplot2)

# 1. Make sure obsPeriod is in true chronological order
obsPeriod <- obsPeriod[order(obsPeriod$bin_start), ]

# 2. Simulate prior and fitted models
fit_sim   <- simEpidemic(betaM.MLE, betaNightA.MLE)
prior_sim <- simEpidemic(PRIOR_betaM, PRIOR_betaNightA)

# 3. Get predicted incidence in the same day/night bins
fit_pred   <- get_period_incidence(fit_sim, obsPeriod)
prior_pred <- get_period_incidence(prior_sim, obsPeriod)

# 4. Build plotting dataframe
plot_dat <- obsPeriod

plot_dat$pred_fit   <- fit_pred$predicted
plot_dat$pred_prior <- prior_pred$predicted

# 5. Create x-axis ordering variable
plot_dat$x_id <- seq_len(nrow(plot_dat))

# 6. Cleaner labels
plot_dat$interval_label <- ifelse(
  plot_dat$Period == "Day",
  paste0("D", plot_dat$PeriodDay),
  paste0("N", plot_dat$PeriodDay)
)

# Quick check
print(plot_dat[, c("x_id", "interval_label", "Period", "new_cases", 
                   "pred_prior", "pred_fit")])

# 7. Plot
ggplot(plot_dat, aes(x = x_id)) +
  geom_col(
    aes(y = new_cases, fill = Period),
    alpha = 0.7,
    width = 0.8
  ) +
  
  geom_line(
    aes(y = pred_prior, colour = "Prior", group = 1),
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_point(
    aes(y = pred_prior, colour = "Prior"),
    size = 2.5
  ) +
  
  geom_line(
    aes(y = pred_fit, colour = "Fitted", group = 1),
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = pred_fit, colour = "Fitted"),
    size = 3
  ) +
  
  scale_x_continuous(
    breaks = plot_dat$x_id,
    labels = plot_dat$interval_label
  ) +
  
  labs(
    title = "Observed vs predicted infections by day/night interval",
    subtitle = "Bars = observed infections; dashed line = prior; solid line = fitted",
    x = "Fitting interval",
    y = "Number of new infections",
    fill = "Observed period",
    colour = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

get_period_location_incidence <- function(simDat, obsPeriodLoc) {
  
  get_val <- function(var, tt) {
    approx(simDat$time, simDat[[var]], xout = tt, rule = 2)$y
  }
  
  predicted <- numeric(nrow(obsPeriodLoc))
  
  for (i in seq_len(nrow(obsPeriodLoc))) {
    
    pop <- obsPeriodLoc$Population[i]
    
    if (pop == "AIMS") {
      var <- "CIA"
    } else {
      var <- "CIE"
    }
    
    start_val <- get_val(var, obsPeriodLoc$bin_start[i])
    end_val   <- get_val(var, obsPeriodLoc$bin_end[i])
    
    predicted[i] <- end_val - start_val
  }
  
  data.frame(
    obsPeriodLoc,
    predicted = pmax(predicted, 1e-10)
  )
}

############################################################
## Plot: observed vs predicted by day/night and population
############################################################

library(ggplot2)

# Simulate prior and fitted models
fit_sim   <- simEpidemic(betaM.MLE, betaNightA.MLE)
prior_sim <- simEpidemic(PRIOR_betaM, PRIOR_betaNightA)

# Get location-specific predictions
fit_pred_loc   <- get_period_location_incidence(fit_sim, obsPeriodLoc)
prior_pred_loc <- get_period_location_incidence(prior_sim, obsPeriodLoc)

# Build plotting data
plot_loc <- obsPeriodLoc

plot_loc$pred_fit   <- fit_pred_loc$predicted
plot_loc$pred_prior <- prior_pred_loc$predicted

# Order properly within each population
plot_loc <- plot_loc[order(plot_loc$Population, plot_loc$bin_start), ]

# Create x-axis id separately for each population
plot_loc$x_id <- ave(
  plot_loc$bin_start,
  plot_loc$Population,
  FUN = function(x) seq_along(x)
)

# Nice labels
plot_loc$interval_label <- ifelse(
  plot_loc$Period == "Day",
  paste0("D", plot_loc$PeriodDay),
  paste0("N", plot_loc$PeriodDay)
)

# Plot
ggplot(plot_loc, aes(x = x_id)) +
  geom_col(
    aes(y = new_cases, fill = Period),
    alpha = 0.7,
    width = 0.8
  ) +
  
  geom_line(
    aes(y = pred_prior, colour = "Prior", group = 1),
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_point(
    aes(y = pred_prior, colour = "Prior"),
    size = 2.5
  ) +
  
  geom_line(
    aes(y = pred_fit, colour = "Fitted", group = 1),
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = pred_fit, colour = "Fitted"),
    size = 3
  ) +
  
  facet_wrap(~ Population, ncol = 1) +
  
  scale_x_continuous(
    breaks = unique(plot_loc$x_id),
    labels = unique(plot_loc$interval_label)
  ) +
  
  labs(
    title = "Observed vs predicted infections by population",
    subtitle = "Bars = observed infections; dashed line = prior; solid line = fitted",
    x = "Day/night interval",
    y = "Number of new infections",
    fill = "Observed period",
    colour = "Model"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )