######################################################################
## SIAR Model - MMF Outbreak (34 cases, closed population)
## State variables: S (Susceptible), IA (Asymptomatic Infectious),
##                  IS (Symptomatic Infectious), R (Recovered)
## Symptom status proxy: Location == 1 → IS (symptomatic)
##                       Location == 0 → IA (asymptomatic)
######################################################################

library(deSolve)
library(ggplot2)
library(tidyr)     

setwd("C:/Users/Ennie Matlhanya/MMEDGit/MMF")
data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)

MMF_DoctorVisits

data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time),
                             format = "%Y-%m-%d %H:%M:%S")

data$Status <- ifelse(data$Location == 1, "IS", "IA")
# summary of counts
n_IA <- sum(data$Status == "IA")  
n_IS <- sum(data$Status == "IS")  
cat("Total cases:", nrow(data), "\n")
cat("Asymptomatic (IA):", n_IA, "\n")
cat("Symptomatic  (IS):", n_IS, "\n")

#Model
siar <- function(t, y, parms) {
  with(c(as.list(y), parms), {

    lambda <- beta * (IA + IS) / N

    dSdt  <-  -lambda * S                     
    dIAdt <-  pA * lambda * S - gammaA * IA   
    dISdt <- (1 - pA) * lambda * S - gammaS * IS  
    dRdt  <-  gammaA * IA + gammaS * IS       

    return(list(c(dSdt, dIAdt, dISdt, dRdt)))
  })
}


N0 <- nrow(data)  #total population = 34 observed cases


I0_IA <- 0   
I0_IS <- 1  
S0    <- N0 - I0_IA - I0_IS
R0_init <- 0

pop.siar <- c(
  S  = S0,
  IA = I0_IA,
  IS = I0_IS,
  R  = R0_init
)

cat("\nInitial conditions:\n")
print(pop.siar)

######################################################################
## Section 4: Parameter values
######################################################################

# Epidemic spans 15/06/2026 to 19/06/2026 = ~4 days
# Infectious period: 24 hours for both IA and IS → gammaA = gammaS = 1
# After 1 day an infectious person can no longer infect others and moves to R
# beta is a placeholder (0.6) (MODEL FITTING)
# pA = proportion of infections that are asymptomatic = 1 - P(IS) = 1 - 0.6 = 0.4

values <- c(
  beta   = 0.6,           
  pA     = 0.4,          
  gammaA = 1,             
  gammaS = 1,             
  N      = N0            
)

cat("\nParameter values:\n")
print(values)

R0_value <- with(as.list(values), {
  beta * (pA / gammaA + (1 - pA) / gammaS)
})
cat("\nBasic Reproduction Number (R0):", round(R0_value, 3), "\n")


# Outbreak lasted ~4 days; model over 10 days to see full trajectory
time.out <- seq(0, 10, by = 0.1)

ts.siar <- data.frame(lsoda(
  y     = pop.siar,
  times = time.out,
  func  = siar,
  parms = values
))

head(ts.siar)
subset(ts.siar, time == 4)   
subset(ts.siar, time == 10)  

data$Day <- as.numeric(data$Date - min(data$Date))  # Days since index case

daily_obs <- aggregate(Infection_number ~ Day + Status, data = data, FUN = length)
names(daily_obs)[3] <- "Count"

cat("\nObserved daily incidence by type:\n")
print(daily_obs)

#Plotting
ts_long <- pivot_longer(ts.siar,
                        cols = c("S", "IA", "IS", "R"),
                        names_to = "Compartment",
                        values_to = "Count")

ggplot(ts_long, aes(x = time, y = Count, colour = Compartment)) +
  geom_line(size = 1.1) +
  scale_colour_manual(
    values = c(S = "steelblue", IA = "orange", IS = "red", R = "forestgreen"),
    labels = c(S = "Susceptible", IA = "Infectious (Asymptomatic)",
               IS = "Infectious (Symptomatic)", R = "Recovered")
  ) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("SIAR Model - MMF Outbreak (N = 34)",
          subtitle = paste0("R0 = ", round(R0_value, 2),
                            ";  pA (asymptomatic fraction) = ",
                            round(values["pA"], 2))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Plot observed incidence overlay


# cases per day (cumulative)
cumulative_obs <- aggregate(Infection_number ~ Day, data = data, FUN = length)
cumulative_obs$Cumulative <- cumsum(cumulative_obs$Infection_number)

ggplot() +
  geom_line(data = ts.siar, aes(x = time, y = IA + IS, colour = "Model (IA+IS)"),
            size = 1.1) +
  geom_point(data = cumulative_obs,
             aes(x = Day, y = Infection_number, colour = "Observed daily cases"),
             size = 3) +
  scale_colour_manual(values = c("Model (IA+IS)" = "red",
                                  "Observed daily cases" = "black")) +
  xlab("Days since index case") +
  ylab("Number of active infectious individuals") +
  ggtitle("SIAR Model vs Observed Outbreak Incidence") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

