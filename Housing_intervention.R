

# --- Libraries ---
library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# --- Initial conditions ---
pop.ssiirr <- c(SE = 17, SA = 25, IA = 1, IE = 1, RA = 0, RE = 0)

values <- c(
  gammaA = 1 , gammaE = 1,
  N0 = 42, N_A = 25, N_E = 17,
  betaMixprior = 1.8, betaNightAprior = 1.5, betaNightEprior = 0
)

# --- ODE function (ssiirr) ---
ssiirr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    beta_M  <- ifelse(is_daytime, betaMixprior, 0)
    beta_NA <- ifelse(is_daytime, 0, betaNightAprior)
    beta_NE <- ifelse(is_daytime, 0, betaNightEprior)
    
    inf_mix <- beta_M * (IA + IE) / N0
    inf_NA  <- beta_NA * IA / N_A
    inf_NE  <- beta_NE * IE / N_E
    
    dSAdt <- -(inf_mix + inf_NA) * SA
    dSEdt <- -(inf_mix + inf_NE) * SE
    dIAdt <- (inf_mix + inf_NA) * SA - gammaA * IA
    dIEdt <- (inf_mix + inf_NE) * SE - gammaE * IE
    dRAdt <- gammaA * IA 
    dREdt <- gammaE * IE 
    
    return(list(c(dSEdt, dSAdt, dIAdt, dIEdt, dRAdt, dREdt)))
  })
}

# --- Helper function to run scenarios ---
run_scenario <- function(betaMix, betaNightA, betaNightE, label) {
  values["betaMixprior"]    <- betaMix
  values["betaNightAprior"] <- betaNightA
  values["betaNightEprior"] <- betaNightE
  time.out <- seq(0, 10, by = 0.02083)
  out <- data.frame(lsoda(y = pop.ssiirr, times = time.out, func = ssiirr, parms = values))
  out$Scenario <- label
  return(out)
}

# --- Run scenarios ---
out_mix    <- run_scenario(1.8, 1.8, 1.8, "Continuous Mixing")
out_iso    <- run_scenario(1.8, 1.5, 1.2, "Nighttime Isolation")
out_strict <- run_scenario(1.8, 0, 0, "Strict Isolation")

# --- Combine results ---
out_all <- rbind(out_mix, out_iso, out_strict)

# --- Summarise final recovered counts ---
final_counts <- out_all %>%
  group_by(Scenario) %>%
  summarise(
    Final_RA = max(RA),
    Final_RE = max(RE),
    Total_R  = Final_RA + Final_RE
  )

print(final_counts)

# --- Plot cumulative infections (RA + RE) ---
ggplot(final_counts, aes(x = Scenario, y = Total_R, fill = Scenario)) +
  geom_col(width = 0.6) +
  labs(
    x = "Scenario",
    y = "Total cumulative infections (RA + RE)",
    title = "Total Outbreak Size Across Housing Policies"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# --- Plot total infectious (IA + IE) across scenarios ---
ggplot(out_all, aes(x = time, y = IA + IE, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("Continuous Mixing"   = "purple",
                                 "Nighttime Isolation" = "magenta",
                                 "Strict Isolation"    = "darkorchid4")) +
  labs(
    x = "Days since index case",
    y = "Total infectious individuals",
    title = "Impact of Housing Interventions on Epidemic Spread",
    subtitle = "Continuous mixing vs nighttime isolation vs strict isolation"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# --- Separate plots for AIMS vs Empire ---
out_long <- pivot_longer(out_all, cols = c("IA","IE"), names_to = "Group", values_to = "Infectious")

ggplot(out_long, aes(x = time, y = Infectious,
                     colour = Scenario, linetype = Group)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Group, scales = "free_y") +
  scale_colour_manual(values = c("Continuous Mixing"   = "purple",
                                 "Nighttime Isolation" = "magenta",
                                 "Strict Isolation"    = "darkorchid4")) +
  scale_linetype_manual(values = c("IA" = "solid", "IE" = "solid")) +
  labs(
    x = "Days since index case",
    y = "Infectious individuals",
    title = "AIMS vs Empire Epidemic Curves",
    subtitle = "Purple = Continuous Mixing, Magenta = Nighttime Isolation, Dark Purple = Strict Isolation"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

