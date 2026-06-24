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


betaMixprior <- 4.82234
betaNightEprior <- 0
betaNightAprior <-0.12757

I0_IE <- sum(data$Location == 1 & data$Infected.by == "Started") #sum(as.integer(index_case$Location == 1))
I0_IA <- sum(data$Location == 0 & data$Infected.by == "Started")
N0    <- 42
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

values <- c(
  gammaA          = 1,
  gammaE          = 1,
  N0              = 42,
  N_A             = 25,
  N_E             = 17,
  betaMixprior    = betaMixprior,
  betaNightAprior = betaNightAprior,
  betaNightEprior = betaNightEprior
)

ssiirr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    # Time-of-day logic
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    # Dynamic betas
    beta_M  <- ifelse(is_daytime, betaMixprior, 0)
    beta_NA <- ifelse(is_daytime, 0, betaNightAprior)
    beta_NE <- ifelse(is_daytime, 0, betaNightEprior)
    
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
time.out <- seq(0, 5, by = 0.02083)

ts.ssiirr <- data.frame(lsoda(
  y     = pop.ssiirr,
  times = time.out,
  func  = ssiirr,
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
# 1. Combined Plot (All Compartments)
# ==============================================================================
ts_long <- pivot_longer(ts.ssiirr,
                        cols      = c("SA","SE", "IA", "IE", "RA", "RE"),
                        names_to  = "Compartment",
                        values_to = "Count")

plotTest <- ggplot(ts_long, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since 15 June 2026") +
  ylab("Number of individuals") +
  ggtitle("Interacting SIR Deterministic Model (assumed beta = 1.8)",) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13))
ggsave("test_axis_size.png", plotTest, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 2. AIMS Plot (Population A)
# ==============================================================================
ts_long_A <- pivot_longer(ts.ssiirr,
                          cols      = c("SA", "IA", "RA"),
                          names_to  = "Compartment",
                          values_to = "Count")

ggplot(ts_long_A, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("AIMS (Population A) Deterministic Model",
          subtitle = paste0("Daytime Beta = ", betaMixprior, " | Nighttime Beta = ", betaNightAprior)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())


# ==============================================================================
# 3. Empire Plot (Population E)
# ==============================================================================
ts_long_E <- pivot_longer(ts.ssiirr,
                          cols      = c("SE", "IE", "RE"),
                          names_to  = "Compartment",
                          values_to = "Count")

ggplot(ts_long_E, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("Empire (Population E) Deterministic Model",
          subtitle = paste0("Daytime Beta = ", betaMixprior, " | Nighttime Beta = ", betaNightEprior)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())



