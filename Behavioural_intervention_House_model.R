# Date    : 24 July 2026
# Authors : Vix Pretorius, edited from Lebogang Mathlanya's code
# Purpose : ODE (Deterministic) Interactive SIR model — two sub-populations (AIMS & Empire)
#           with time-varying (daytime / nighttime) betas

# ---------------------------------------------------------------------------
# Model Description
# ---------------------------------------------------------------------------
# Two sub-populations sharing AIMS during the day:
#   Population A — residents of AIMS     (N_A = 25)
#   Population E — residents of Empire   (N_E = 17)
#
# Compartments: SA, SE, IA, IE, RA, RE
#
# Force of infection mirrors the ODE exactly:
#   DAYTIME   (day_fraction in [0.3328, 0.77071]):
#     inf_mix  = betaMix  * (IA + IE) / N0   <- shared AIMS mixing
#     inf_NA   = 0                             <- no AIMS-only night contact
#     inf_NE   = 0                             <- no Empire-only night contact
#
#   NIGHTTIME (all other hours):
#     inf_mix  = 0
#     inf_NA   = betaNightA * IA / N_A        <- AIMS residents only
#     inf_NE   = betaNightE * IE / N_E        <- Empire residents only
#
# Transitions and rates (exactly matching ODE derivatives):
#   Event               Change                  Rate
#   Infect SA (day)     SA -> IA                (inf_mix) * SA
#   Infect SA (night)   SA -> IA                (inf_NA)  * SA
#   Infect SE (day)     SE -> IE                (inf_mix) * SE
#   Infect SE (night)   SE -> IE                (inf_NE)  * SE
#   Recover IA          IA -> RA                gammaA * IA
#   Recover IE          IE -> RE                gammaE * IE

library(deSolve)
library(ggplot2)
library(ellipse)
library(tidyr)
library(chron)

data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)


data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time),
                            format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))

# Location == 1 → IE (Empire), Location == 0 → IA (AIMS)
data$Status <- ifelse(data$Location == 1, "IE", "IA")

betaMixpost <- 3.1
betaNightEpost <- 0
betaNightApost <-0.9

index_cases <- data[data$Infected.by == "Started", ]
I0_IE <- sum(data$Location == 1 & data$Infected.by == "Started") 
I0_IA <- sum(data$Location == 0 & data$Infected.by == "Started")
#I0_IA <- sum(as.integer(index_case$Location == 0))
N0  <- 42
N_A <- 25 
N_E <- 17
S0_E <- N_E - I0_IE
S0_A <- N_A - I0_IA

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

# Update your values vector to include all constants needed by the function
values <- c(
  gammaA          = 1,
  gammaE          = 1,
  N0              = 42,
  N_A             = 25,
  N_E             = 17,
  betaMixpost    = betaMixpost,
  betaNightApost = betaNightApost,
  betaNightEpost = betaNightEpost
)

ssiirr_new <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    # Time-of-day logic
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    p_day   <- 0.85
    p_night <- 0.30
    
    # Calculate transmission reduction factors based on the joint probabilities
    factor_day   <- 0.03*p_day^2   + 0.21*p_day*(1-p_day)     + (1-p_day)^2
    factor_night <- 0.03*p_night^2 + 0.21*p_night*(1-p_night) + (1-p_night)^2
   
    # Apply factors to your transmission coefficients
    beta_M  <- ifelse(is_daytime, betaMixpost * factor_day, 0)
    beta_NA <- ifelse(is_daytime, 0, betaNightApost * factor_night)
    beta_NE <- ifelse(is_daytime, 0, betaNightEpost * factor_night)
    # Force of infection (per house)
    inf_mix <- beta_M * (IA + IE) / N0
    inf_NA   <- beta_NA * IA / N_A
    inf_NE   <- beta_NE * IE / N_E
    
    # Derivatives
    dSAdt <- -(inf_mix + inf_NA) * SA
    dSEdt <- -(inf_mix + inf_NE) * SE
    
    dIAdt <- (inf_mix + inf_NA) * SA - gammaA * IA
    dIEdt <- (inf_mix + inf_NE) * SE - gammaE * IE
    
    dRAdt <- gammaA * IA 
    dREdt <- gammaE * IE 
    
    # CRITICAL FIX: The order here now perfectly matches pop.ssiirr
    # SE, SA, IA, IE, RA, RE
    return(list(c(dSEdt, dSAdt, dIAdt, dIEdt, dRAdt, dREdt)))
  })
}

# 3. Run the model
time.out <- seq(0, 10, by = 0.02083)

ts.ssiirr_new <- data.frame(lsoda(
  y     = pop.ssiirr,
  times = time.out,
  func  = ssiirr_new,
  parms = values
))

ssiirr_base <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    # Time-of-day logic
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    # Dynamic betas
    beta_M  <- ifelse(is_daytime, betaMixpost, 0)
    beta_NA <- ifelse(is_daytime, 0, betaNightApost)
    beta_NE <- ifelse(is_daytime, 0, betaNightEpost)
    
    # Force of infection (per house)
    inf_mix <- beta_M * (IA + IE) / N0
    inf_NA   <- beta_NA * IA / N_A
    inf_NE   <- beta_NE * IE / N_E
    
    # Derivatives
    dSAdt <- -(inf_mix + inf_NA) * SA
    dSEdt <- -(inf_mix + inf_NE) * SE
    
    dIAdt <- (inf_mix + inf_NA) * SA - gammaA * IA
    dIEdt <- (inf_mix + inf_NE) * SE - gammaE * IE
    
    dRAdt <- gammaA * IA 
    dREdt <- gammaE * IE 
    
    return(list(c(dSEdt, dSAdt, dIAdt, dIEdt, dRAdt, dREdt)))
  })
}

# 3. Run the model
time.out <- seq(0, 10, by = 0.02083)

ts.ssiirr_base <- data.frame(lsoda(
  y     = pop.ssiirr,
  times = time.out,
  func  = ssiirr_base,
  parms = values
))

pop_colours <- c(
  SA = "lightblue4", IA = "lightblue4", RA = "lightblue4",  # AIMS
  SE = "orange", IE = "orange", RE = "orange" # Empire
)

pop_labels <- c(
  SA = "Susceptible living at AIMS",
  SE = "Susceptible living at Empire",
  IA = "Infectious (living at AIMS)",
  IE = "Infectious (living at Empire)",
  RA = "Recovered at AIMS", 
  RE = "Recovered at Empire"
)

state_linetypes <- c(
  SA = "dotted",    SE = "dotted",
  IA = "solid",     IE = "solid",
  RA = "dotdash",   RE = "dotdash"
)


# ==============================================================================
# Side-by-Side Model Comparison: IA and IE trajectories (Color-Coded by Model)
# ==============================================================================

# Add scenario identifiers to each model output
ts.ssiirr_base$Scenario <- "Baseline Model"
ts.ssiirr_new$Scenario  <- "Behavioural Change"

# Combine the datasets
ts_combined <- rbind(ts.ssiirr_base, ts.ssiirr_new)

# Pivot long selecting only IA and IE compartments
ts_infectious_long <- pivot_longer(
  ts_combined,
  cols      = c("IA", "IE"),
  names_to  = "Compartment",
  values_to = "Count"
)

# Create clean display labels for the side-by-side facet panels
facet_labels <- c(
  IA = "AIMS",
  IE = "Empire"
)

# Color mapping: Base model is magenta, New model is black
model_colours <- c(
  "Baseline Model"     = "magenta", 
  "Behavioural Change" = "violetred4"
)

plotInter <- ggplot(ts_infectious_long, aes(x = time, y = Count, colour = Scenario)) +
  geom_line(linewidth = 1.2, linetype = "solid") + # Forces all lines to be solid
  # Facet side-by-side by compartment (IA next to IE)
  facet_wrap(~ Compartment, labeller = as_labeller(facet_labels)) +
  scale_colour_manual(values = model_colours) +
  xlab("Days since 15 June 2026") +
  ylab("Number of Infectious Individuals") +
  ggtitle("AIMS vs Empire Epidemic Curves",
          subtitle = "Magenta = Baseline, Violet Red = Behavioural Change") +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(face = "bold"),
    axis.title = element_text(size = 15), 
    axis.text = element_text(size = 13))
ggsave("Intervention Plot - Behavioural Change.png", plotInter, width = 8, height = 6, dpi = 300)



