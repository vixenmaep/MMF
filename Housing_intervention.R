library(deSolve)
library(ggplot2)
library(tidyr)

# --- Helper function to run your ODEs under different assumptions ---
run_scenario <- function(betaMix, betaNightA, betaNightE, label) {
  values <- c(
    gammaA          = 1,
    gammaE          = 1,
    N0              = 42,
    N_A             = 25,
    N_E             = 17,
    betaMixprior    = betaMix,
    betaNightAprior = betaNightA,
    betaNightEprior = betaNightE
  )
  
  time.out <- seq(0, 10, by = 0.02083)
  out <- data.frame(lsoda(y = pop.ssiirr, times = time.out, func = ssiirr, parms = values))
  out$Scenario <- label
  return(out)
}


# Extract final recovered counts
final_counts <- out_all %>%
  group_by(Scenario) %>%
  summarise(
    Final_RA = max(RA),
    Final_RE = max(RE),
    Total_R  = Final_RA + Final_RE
  )

print(final_counts)

# Plot cumulative infections
ggplot(final_counts, aes(x = Scenario, y = Total_R, fill = Scenario)) +
  geom_col(width = 0.6) +
  labs(
    x = "Scenario",
    y = "Total cumulative infections (RA + RE)",
    title = "Total Outbreak Size Across Housing Policies"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")


# --- Scenario 1: Continuous mixing (day + night) ---
out_mix <- run_scenario(betaMix = 1.8, betaNightA = 1.8, betaNightE = 1.8, label = "Continuous Mixing")

# --- Scenario 2: Nighttime isolation (within-group spread only) ---
out_iso <- run_scenario(betaMix = 1.8, betaNightA = 1.5, betaNightE = 1.2, label = "Nighttime Isolation")

# --- Scenario 3: Strict isolation (no spread at night) ---
out_strict <- run_scenario(betaMix = 1.8, betaNightA = 0, betaNightE = 0, label = "Strict Isolation")

# --- Combine results ---
out_all <- rbind(out_mix, out_iso, out_strict)

# --- Plot total infectious (IA + IE) across scenarios ---
ggplot(out_all, aes(x = time, y = IA + IE, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Days since index case",
    y = "Total infectious individuals",
    title = "Impact of Housing Interventions on Epidemic Spread",
    subtitle = "Continuous mixing vs nighttime isolation vs strict isolation"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# --- Optional: Separate plots for AIMS vs Empire ---
out_long <- pivot_longer(out_all, cols = c("IA","IE"), names_to = "Group", values_to = "Infectious")

ggplot(out_long, aes(x = time, y = Infectious, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Group, scales = "free_y") +
  labs(
    x = "Days since index case",
    y = "Infectious individuals",
    title = "AIMS vs Empire Epidemic Curves",
    subtitle = "Comparing housing policies"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

library(ggplot2)

# Reshape data for plotting IA (AIMS) and IE (Empire)
out_long <- pivot_longer(out_all,
                         cols = c("IA","IE"),
                         names_to = "Group",
                         values_to = "Infectious")

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





