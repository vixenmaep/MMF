# Date : 22 July 2026
# Author : Mandita Rakei
# Purpose : Create Gillespie Model for SIR_baseline 

# Model Description

# We have a closed population of individuals at AIMS who are susceptible to MMF.
# Once infected some show symptoms (with probability Ps).
# Other infected individuals do are asymptotic (with probability (1-Ps)).
# We assume that they have the same contact rate (beta).
# Both sets of individuals recover at a rate of (gamma).

## Compartments:
## (S,(I_A, I_S),R) = (susceptible, infectious (symptoms & none), recovered)

## Transitions:
## Event                           Change        	  	          Rate
## Infection (I_A)                 (S,I_A,R)->(S-1,I_A+1,R)     beta*(I_A + I_S)*(1-Ps)S/N
## Infection(I_S)                  (S, I_S,R)->(S-1,I_S+1,R)    beta*(I_A + I_S)*(Ps)S/N
## Recovery/Removal (I_A)          (S,I_A,R)->(S,I_A-1,R+1)     gamma*(I_A)
## Recovery/Removal (I_A)          (S,I_S,R)->(S,I_S-1,R+1)     gamma*(I_S)

## 
# TODO LIST:
# make it hourly, so day*24hrs
# Make population 46
# have to use the MLE beta from Deterministic ODE MLE Fit
# Step 0: Load essential packages 

library(tidyverse) # has ggplot and other fancy stuff has more libraries rolled in one.
library(lubridate)
library(hms)

# Step 1: Load in Data for comparison to Gillespie 

df <- read_csv("MMF-Final+Locations+Doctors.csv")


head(df)

# Step 2: Function to step forward in time to next event and update states:
event_sir_base <- function(time, S, I_A, I_S, R, params, t_end) {
  
  with(as.list(params), {
    N <- S + I_A + I_S + R
    
    # Define 4 separate events for clean Gillespie state transitions
    rates <- c(
      infect_A  = beta * (I_A + I_S) * (1 - Ps) * S / N,
      infect_S  = beta * (I_A + I_S) * Ps * S / N,
      recover_A = gamma * I_A,
      recover_S = gamma * I_S
    )
    
    total_rate <- sum(rates)
    
    if (total_rate == 0) {
      count.inf <- 0
      event_time <- t_end
    } else {
      event_time <- time + rexp(1, total_rate)
      event_type <- sample(names(rates), 1, prob = rates / total_rate)
      
      switch(event_type,
             "infect_A" = {
               S <- S - 1
               I_A <- I_A + 1
               count.inf <- 1
             },
             "infect_S" = {
               S <- S - 1
               I_S <- I_S + 1
               count.inf <- 1
             },
             "recover_A" = {
               I_A <- I_A - 1
               R <- R + 1
               count.inf <- 0
             },
             "recover_S" = {
               I_S <- I_S - 1
               R <- R + 1
               count.inf <- 0
             })
    }
    
    return(data.frame(time = event_time, S = S, I_A = I_A, I_S = I_S, R = R, count.inf = count.inf))
  })
}

# Step 3: Function to simulate states from time 0 to t_end:
simulate_sir_base <- function(t_end, y, params) {
  with(as.list(y), {
    
    cum.inf <- 0
    # Initialize data frame with the correct columns matching your equations
    ts <- data.frame(time = 0, S = S, I_A = I_A, I_S = I_S, R = R, count.inf = 0) 
    next_event <- ts
    
    while (next_event$time < t_end) {
      next_event <- event_sir_base(
        time = next_event$time, 
        S = next_event$S, 
        I_A = next_event$I_A, 
        I_S = next_event$I_S, 
        R = next_event$R, 
        params = params, 
        t_end = t_end
      )
      
      cum.inf <- cum.inf + next_event$count.inf
      next_event$count.inf <- cum.inf
      ts <- rbind(ts, next_event)
    }
    
    return(ts)
  })
}

# Step 4: Run the model execution block
pop <- 50                                     # Population size (Should this be total MMED participants)
params <- c(beta = 0.3, gamma = 0.1, Ps = 0.6) # Parameters
final_time <- 400                              # Simulation timeline
y0 <- c(S = pop - 1, I_A = 1, I_S = 0, R = 0)  # Clean split of initial conditions

# Call the correct base simulation function
ts1 <- simulate_sir_base(final_time, y0, params)

# Step 5: Pivot and Plot
ts1_long <- (ts1 
             |> pivot_longer(cols = c(S, I_A, I_S, R), names_to = "Compartment", values_to = "count")
             |> mutate(Compartment = factor(Compartment, levels = c('S', 'I_A', 'I_S', 'R'))) 
)

ggplot(ts1_long, aes(x = time, y = count, color = Compartment)) +
  geom_step(linewidth = 1.2) +
  labs(title = "Gillespie SIR Dynamics with Asymptomatic Split", 
       y = "Count", 
       x = "Time") +
  theme_minimal(base_size = 14)


