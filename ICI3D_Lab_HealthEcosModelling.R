
#####################################################################

## Lab: Introduction to health economics in dynamical modelling 

#####################################################################
## Clinic on the Meaningful Modeling of Epidemiological Data
## International Clinics on Infectious Disease Dynamics and Data (ICI3D) Program
## https://www.ici3d.org
##
## Attribution: Mmamapudi Kubjane (2026)
##             
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)
##
#####################################################################

##		The goal of this to lab is to complete a simple health 
##		economics modelling exercise in R.
##
##		In this lab, you will: 
##    - Define and simulate alternative vaccination intervention scenarios
##    - Estimate the epidemiological impact of vaccination on disease incidence
##    - Attach costs to key health outcomes and interventions
##    - Calculate and interpret the Incremental Cost-Effectiveness Ratio (ICER)
##
##    In this lab, we will build on 'Lab: ODE models in R'
##      (https://github.com/ICI3D/RTutorials/blob/master/ICI3D_Lab_ODEmodels.R). 
##    We will extend the SIR model to include a vaccination intervention.
##		

#####################################################################
# Install packages if haven't got them
# install.packages("tidyverse")
# install.packages("deSolve")
# install.packages("reshape2")
 
# Load libraries
library(readr)
library(dplyr)
library(deSolve)

# -----------------------------
# Step 1: Read and prepare data


#Data prep

#the  CSV (MMF-Final+Duplicates.csv).

#Groups infections by date to get obs_cases (daily counts).

#That’s our real data
# -----------------------------
mmf <- read_csv("MMF-Final+Duplicates.csv")

# Convert Date to proper format
mmf <- mmf %>% mutate(Date = as.Date(Date, format="%d/%m/%Y"))

# Count daily infections
daily_cases <- mmf %>%
  group_by(Date) %>%
  summarise(cases = n()) %>%
  arrange(Date)

obs_cases <- daily_cases$cases
times <- seq(0, length(obs_cases), by=1)   # include day 0

# -----------------------------
# Step 2: SIIR model equations

#Model simulation

#Runs your SIIR equations with deSolve::ode.

#Produces cumcases (cumulative infections).

#Takes diff(cumcases) to get daily incidence — the model’s prediction.
# -----------------------------
sir_equations <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I_A + I_S + R
    beta <- P_D * C
    
    dS   <- -beta * S * (I_A + I_S)/N
    dI_S <-  beta * S * (I_A + I_S)/N * p - gamma * I_S
    dI_A <-  beta * S * (I_A + I_S)/N * (1 - p) - gamma * I_A
    dR   <-  gamma * I_A + gamma * I_S
    
    cumcases <- beta * S * (I_A + I_S)/N
    
    return(list(c(dS, dI_S, dI_A, dR, cumcases)))
  })
}

# -----------------------------
# Step 3: Likelihood function


#Compares obs_cases (real daily counts) to pred (model daily incidence).

#Uses a Poisson likelihood to measure how well the model matches the data.

#optim() adjusts parameters ( P_𝐷,𝐶,𝑝,gammA) to maximize the fit.
# -----------------------------
loglik <- function(par, data, times, init) {
  names(par) <- c("P_D", "C", "p", "gamma")
  
  out <- ode(y = init, times = times, func = sir_equations, parms = par)
  pred <- diff(out[ , "cumcases"])  # daily incidence
  
  # Ensure lengths match observed data
  pred <- pred[1:length(data)]
  
  sum(dpois(data, lambda = pmax(pred, 1e-6), log = TRUE))
}

# -----------------------------
# Step 4: Fit parameters
# -----------------------------
init <- c(S = 1000, I_A = 1, I_S = 1, R = 0, cumcases = 0)

fit <- optim(par = c(P_D = 0.05, C = 10, p = 0.3, gamma = 0.1),
             fn = function(par) -loglik(par, obs_cases, times, init),
             method = "Nelder-Mead")

print(fit$par)

# -----------------------------
# Step 5: Visual check
# -----------------------------
out <- ode(y = init, times = times, func = sir_equations, parms = fit$par)
pred <- diff(out[ , "cumcases"])
pred <- pred[1:length(obs_cases)]   # match lengths

plot(seq_along(obs_cases), obs_cases, type="b", col="red", pch=19,
     ylab="Daily cases", xlab="Day", main="MMF Outbreak: Observed vs Fitted")
lines(seq_along(pred), pred, col="blue", lwd=2)
legend("topright", legend=c("Observed","Fitted"),
       col=c("red","blue"), lty=c(1,1), pch=c(19,NA))


# ------------------------------------------------------------------------------
# MCMC
# -----------------------------------------------------------------------------
# Load libraries
library(readr)
library(dplyr)
library(deSolve)
library(ggplot2)

# -----------------------------
# Step 1: Read and prepare data
# -----------------------------
mmf <- read_csv("MMF-Final+Duplicates.csv")
mmf <- mmf %>% mutate(Date = as.Date(Date, format="%d/%m/%Y"))

daily_cases <- mmf %>%
  group_by(Date) %>%
  summarise(cases = n()) %>%
  arrange(Date)

obs_cases <- daily_cases$cases
times <- seq(0, length(obs_cases), by=1)

# -----------------------------
# Step 2: SIIR model equations
# -----------------------------
sir_equations <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I_A + I_S + R
    beta <- P_D * C
    
    dS   <- -beta * S * (I_A + I_S)/N
    dI_S <-  beta * S * (I_A + I_S)/N * p - gamma * I_S
    dI_A <-  beta * S * (I_A + I_S)/N * (1 - p) - gamma * I_A
    dR   <-  gamma * I_A + gamma * I_S
    
    cumcases <- beta * S * (I_A + I_S)/N
    
    return(list(c(dS, dI_S, dI_A, dR, cumcases)))
  })
}

# -----------------------------
# Step 3: Likelihood function
# -----------------------------
likelihood <- function(par, data, times, init) {
  names(par) <- c("P_D","C","p","gamma")
  out <- ode(y = init, times = times, func = sir_equations, parms = par)
  pred <- diff(out[ , "cumcases"])
  pred <- pred[1:length(data)]
  sum(dpois(data, lambda = pmax(pred, 1e-6), log = TRUE))
}

# -----------------------------
# Step 4: Priors
# -----------------------------
prior <- function(par) {
  P_D <- par[1]; C <- par[2]; p <- par[3]; gamma <- par[4]
  dunif(P_D, 0, 1, log=TRUE) +
    dunif(C, 0, 50, log=TRUE) +
    dbeta(p, 2, 2, log=TRUE) +
    dunif(gamma, 0, 1, log=TRUE)
}

posterior <- function(par, data, times, init) {
  prior(par) + likelihood(par, data, times, init)
}

# -----------------------------
# Step 5: MCMC loop
# -----------------------------

# Load libraries
library(readr)
library(dplyr)
library(deSolve)
library(ggplot2)

# -----------------------------
# Step 1: Read and prepare data
# -----------------------------
mmf <- read_csv("MMF-Final+Duplicates.csv")
mmf <- mmf %>% mutate(Date = as.Date(Date, format="%d/%m/%Y"))

daily_cases <- mmf %>%
  group_by(Date) %>%
  summarise(cases = n()) %>%
  arrange(Date)

obs_cases <- daily_cases$cases
times <- seq(0, length(obs_cases), by=1)

# -----------------------------
# Step 2: SIIR model equations
# -----------------------------
sir_equations <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + I_A + I_S + R
    beta <- P_D * C
    
    dS   <- -beta * S * (I_A + I_S)/N
    dI_S <-  beta * S * (I_A + I_S)/N * p - gamma * I_S
    dI_A <-  beta * S * (I_A + I_S)/N * (1 - p) - gamma * I_A
    dR   <-  gamma * I_A + gamma * I_S
    
    cumcases <- beta * S * (I_A + I_S)/N
    
    return(list(c(dS, dI_S, dI_A, dR, cumcases)))
  })
}

# -----------------------------
# Step 3: Likelihood function
# -----------------------------
likelihood <- function(par, data, times, init) {
  names(par) <- c("P_D","C","p","gamma")
  out <- ode(y = init, times = times, func = sir_equations, parms = par)
  pred <- diff(out[ , "cumcases"])
  pred <- pred[1:length(data)]
  sum(dpois(data, lambda = pmax(pred, 1e-6), log = TRUE))
}

# -----------------------------
# Step 4: Priors
# -----------------------------
prior <- function(par) {
  P_D <- par[1]; C <- par[2]; p <- par[3]; gamma <- par[4]
  dunif(P_D, 0, 1, log=TRUE) +
    dunif(C, 0, 50, log=TRUE) +
    dbeta(p, 2, 2, log=TRUE) +
    dunif(gamma, 0, 1, log=TRUE)
}

posterior <- function(par, data, times, init) {
  prior(par) + likelihood(par, data, times, init)
}

# -----------------------------
# Step 5: MCMC loop
# -----------------------------
iterations <- 5000
chain <- matrix(NA, nrow=iterations, ncol=4)
ll_chain <- numeric(iterations)

init <- c(S = 1000, I_A = 1, I_S = 1, R = 0, cumcases = 0)
par <- c(P_D=0.05, C=10, p=0.3, gamma=0.1)
ll <- posterior(par, obs_cases, times, init)

for (i in 1:iterations) {
  proposal <- par + rnorm(4, 0, c(0.01, 1, 0.01, 0.01)) # step sizes
  ll_prop <- posterior(proposal, obs_cases, times, init)
  ratio <- exp(ll_prop - ll)
  
  if (runif(1) < ratio) {
    par <- proposal
    ll <- ll_prop
  }
  
  chain[i, ] <- par
  ll_chain[i] <- ll
}

# -----------------------------
# Step 6: Diagnostics
# -----------------------------
colnames(chain) <- c("P_D","C","p","gamma")

# Trace plots
matplot(chain, type="l", lty=1, main="MCMC Trace Plots")

# Posterior histograms
par(mfrow=c(2,2))
hist(chain[,1], main="Posterior of P_D", xlab="")
hist(chain[,2], main="Posterior of C", xlab="")
hist(chain[,3], main="Posterior of p", xlab="")
hist(chain[,4], main="Posterior of gamma", xlab="")



















# -----------------------------
# Step 6: Diagnostics
# -----------------------------
colnames(chain) <- c("P_D","C","p","gamma")

# Trace plots
matplot(chain, type="l", lty=1, main="MCMC Trace Plots")

# Posterior histograms
par(mfrow=c(2,2))
hist(chain[,1], main="Posterior of P_D", xlab="")
hist(chain[,2], main="Posterior of C", xlab="")
hist(chain[,3], main="Posterior of p", xlab="")
hist(chain[,4], main="Posterior of gamma", xlab="")

ScenarioList <- c(
  "Baseline (no vaccine)",
  "Intervention (vaccine)"
)

  
 
# ------------------------------------------------------------------------------
# Specify cost parameters 
# ------------------------------------------------------------------------------
cost_vaccine   <- 10    # USD, vaccine course per person
cost_treatment <- 200   # USD, treatment course per case
#   Note: For simplicity test/diagnosis cost is bundled into cost_treatment;
#   only individuals who are diagnosed and treated incur this cost.
#   Treatment uptake: proportion of incident cases that reach treatment.
#   Accounts for lost to follow-up, non-diagnosis, refusal, incomplete treatment.
#   
treatment_uptake <- 0.85   # 85% treated, 15% lost individuals 

# ------------------------------------------------------------------------------
# Specify intervention (vaccine) parameters 
# ------------------------------------------------------------------------------
vaccine_coverage <- 0.40   # Proportion of susceptible offered vaccine
vaccine_efficacy <- 1.00   # 100% protection in those vaccinated
# Individuals not vaccinated remain in S, fully susceptible.

# ------------------------------------------------------------------------------
# SIR Differential equations as in Lab 1 (re-review if have forgotten)
# ------------------------------------------------------------------------------

# View all the calculations performed
Analysis
# what information is this data frame showing?

# ==============================================================================
# Discussion Questions
# ==============================================================================
#
# Work through these questions by modifying the parameter values above and
# re-running the relevant sections of the script.
#
# ------------------------------------------------------------------------------
# Q: What if vaccine coverage is lower?
# ------------------------------------------------------------------------------
#   Coverage reflects the proportion of susceptibles who are offered and accept
#   the vaccine. Low coverage may result from supply constraints, hesitancy, or
#   access barriers.
#
#   a) Change vaccine_coverage to 20% and 50% in turn, and re-run the model.
#      How does this affect total cases and the ICER?
#
#   b) Compare the effect of halving coverage vs. halving efficacy.
#      Which has a bigger impact on cases averted and on the ICER? Why might
#      that be?
#
# ------------------------------------------------------------------------------
