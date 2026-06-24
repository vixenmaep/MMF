library(deSolve)
library(ggplot2)
library(tidyr)

data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)

data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time), format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))
data$Status   <- ifelse(data$Location == 1, "IE", "IA")

daily_total           <- aggregate(Infection_number ~ Day, data = data, FUN = length)
names(daily_total)[2] <- "new_cases"

all_days <- data.frame(Day = 0:max(data$Day))
obsDat   <- merge(all_days, daily_total, by = "Day", all.x = TRUE)
obsDat[is.na(obsDat)] <- 0


PRIOR_betaM      <- 1.8   
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
    CI = 0)
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
    
    dSEdt <- -(inf_mix + inf_NE) * SE
    dSAdt <- -(inf_mix + inf_NA) * SA
    dIAdt <-  (inf_mix + inf_NA) * SA - gammaA * IA
    dIEdt <-  (inf_mix + inf_NE) * SE - gammaE * IE
    dRAdt <-   gammaA * IA
    dREdt <-   gammaE * IE
    dCIdt <-  (inf_mix + inf_NA) * SA + (inf_mix + inf_NE) * SE
    
    return(list(c(dSEdt, dSAdt, dIAdt, dIEdt, dRAdt, dREdt, dCIdt)))
  })
}

simEpidemic <- function(betaM, betaNightA,
                        tseq = seq(0, max(obsDat$Day), by = 0.02083)) {
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

get_daily_incidence <- function(simDat, days = obsDat$Day) {
  ci_vals   <- sapply(days, function(d) simDat$CI[which.min(abs(simDat$time - d))])
  daily_inc <- c(ci_vals[1], diff(ci_vals))
  data.frame(Day = days, predicted = daily_inc)
}

nllikelihood <- function(betaM, betaNightA) {
  sim    <- simEpidemic(betaM, betaNightA)
  pred   <- get_daily_incidence(sim)
  lambda <- pmax(pred$predicted, 1e-10)
  nlls   <- -dpois(round(obsDat$new_cases), lambda = lambda, log = TRUE)
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

fit_daily   <- get_daily_incidence(fit_sim)
prior_daily <- get_daily_incidence(prior_sim)



