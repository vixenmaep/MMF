library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# --- Initial conditions ---
pop.ssiirr <- c(SE = 16, SA = 23, IA = 2, IE = 1, RA = 0, RE = 0)

values <- c(
  gammaA = 1 , gammaE = 1,
  N0 = 42, N_A = 25, N_E = 17,
  betaMixprior = 3.1, betaNightAprior = 0.9, betaNightEprior = 0
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
out_mix    <- run_scenario(3.1, 3.1, 3.1, "Continuous Mixing")
out_iso    <- run_scenario(3.1, 0.9, 0, "Baseline")
out_strict <- run_scenario(3.1, 0, 0, "Strict Isolation")

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
                                 "Baseline" = "magenta",
                                 "Strict Isolation"    = "darkorchid4")) +
  labs(
    x = "Days since index case",
    y = "Total infectious individuals",
    title = "Impact of Housing Interventions on Epidemic Spread",
    subtitle = "Continuous mixing vs Baseline vs strict isolation"
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
                                 "Baseline" = "magenta",
                                 "Strict Isolation"    = "darkorchid4")) +
  scale_linetype_manual(values = c("IA" = "solid", "IE" = "solid")) +
  labs(
    x = "Days since index case",
    y = "Infectious individuals",
    title = "AIMS vs Empire Epidemic Curves",
    subtitle = "Purple = Continuous Mixing, Magenta = Baseline, Dark Purple = Strict Isolation"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())


ssiirr_swap <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    day_fraction <- t %% 1
    is_daytime   <- (day_fraction >= 0.3328) & (day_fraction <= 0.77071)
    
    if (is_daytime) {
      # Daytime mixing
      beta_M  <- betaMixprior
      inf_mix <- beta_M * (IA + IE) / N0
      
      dSAdt <- -inf_mix * SA
      dSEdt <- -inf_mix * SE
      dIAdt <- inf_mix * SA - gammaA * IA
      dIEdt <- inf_mix * SE - gammaE * IE
    } else {
      # Nighttime swapping
      beta_swap <- betaNightAprior  # use same parameter for swap intensity
      
      inf_A_swap <- beta_swap * IE / N_E
      inf_E_swap <- beta_swap * IA / N_A
      
      dSAdt <- -inf_A_swap * SA
      dSEdt <- -inf_E_swap * SE
      dIAdt <- inf_A_swap * SA - gammaA * IA
      dIEdt <- inf_E_swap * SE - gammaE * IE
    }
    
    dRAdt <- gammaA * IA
    dREdt <- gammaE * IE
    
    return(list(c(dSEdt, dSAdt, dIAdt, dIEdt, dRAdt, dREdt)))
  })
}

# --- Run swap scenario ---
run_swap <- function(betaMix, betaSwap, label) {
  values["betaMixprior"]    <- betaMix
  values["betaNightAprior"] <- betaSwap
  time.out <- seq(0, 10, by = 0.02083)
  out <- data.frame(lsoda(y = pop.ssiirr, times = time.out, func = ssiirr_swap, parms = values))
  out$Scenario <- label
  return(out)
}

out_swap <- run_swap(3.1, 0.9 , "Nighttime Swapping")
out_all <- rbind(out_mix, out_iso, out_strict, out_swap)

ggplot(out_all, aes(x = time, y = IA + IE, colour = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(values = c("Continuous Mixing"   = "purple",
                                 "Baseline"           = "magenta",
                                 "Strict Isolation"   = "darkorchid4",
                                 "Nighttime Swapping" = "plum3")) +
  labs(
    x = "Days since index case",
    y = "Total infectious individuals",
    title = "Impact of Housing Interventions on Epidemic Spread",
    subtitle = "Continuous mixing vs baseline vs strict isolation vs swapping"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

out_long <- pivot_longer(out_all, cols = c("IA","IE"), names_to = "Group", values_to = "Infectious")
facet_labels <- c(
  IA = "AIMS",
  IE = "Empire"
)
plotInter <- ggplot(out_long, aes(x = time, y = Infectious,
                     colour = Scenario, linetype = Group)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Group, scales = "free_y", ~ Compartment, labeller = as_labeller(facet_labels)) +
   +
  scale_colour_manual(values = c("Continuous Mixing"   = "purple",
                                 "Baseline" = "magenta",
                                 "Strict Isolation"    = "darkorchid4",
                                 "Nighttime Swapping" = "plum3")) +
  scale_linetype_manual(values = c("IA" = "solid", "IE" = "solid")) +
  labs(
    x = "Days since index case",
    y = "Infectious individuals",
    title = "AIMS vs Empire Epidemic Curves",
    subtitle = "Housing/Quarantine Intervention"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 13))
ggsave("Intervention Plot - Housing Change.png", plotInter, width = 8, height = 6, dpi = 300)

 