# Date    : 22 July 2026
# Authors : Mandita Rakei (original), updated to match SSIIRR deterministic ODE
# Purpose : Gillespie (stochastic) SIR model — two sub-populations (AIMS & Empire)
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
#     inf_mix  = betaMix  * (IA + IE) / N0  <- shared AIMS mixing
#     inf_NA   = 0                             <- no AIMS-only night contact
#     inf_NE   = 0                             <- no Empire-only night contact
#
#   NIGHTTIME (all other hours):
#     inf_mix  = 0
#     inf_NA   = betaNightA * IA / N_A         <- AIMS residents only
#     inf_NE   = betaNightE * IE / N_E         <- Empire residents only
#
# Transitions and rates (exactly matching ODE derivatives):
#   Event                Change                 Rate
#   Infect SA (day)      SA -> IA                (inf_mix) * SA
#   Infect SA (night)    SA -> IA                (inf_NA)  * SA
#   Infect SE (day)      SE -> IE                (inf_mix) * SE
#   Infect SE (night)    SE -> IE                (inf_NE)  * SE
#   Recover IA           IA -> RA                gammaA * IA
#   Recover IE           IE -> RE                gammaE * IE

# ---------------------------------------------------------------------------
# Step 0: Load packages
# ---------------------------------------------------------------------------

library(tidyverse)
library(lubridate)

# ---------------------------------------------------------------------------
# Step 1: Load data and set initial conditions (mirrors ODE script exactly)
# ---------------------------------------------------------------------------

data <- read.csv("MMF-Final+Locations+Doctors.csv", stringsAsFactors = FALSE)

data$Date     <- as.Date(data$Date, format = "%d/%m/%Y")
data$DateTime <- as.POSIXct(paste(data$Date, data$Time),
                            format = "%Y-%m-%d %H:%M:%S")
data$Day      <- as.numeric(data$Date - min(data$Date))

# Location == 1 -> IE (Empire), Location == 0 -> IA (AIMS)
data$Status <- ifelse(data$Location == 1, "IE", "IA")

# Index cases — same logic as ODE script
I0_IE <- sum(data$Location == 1 & data$Infected.by == "Started")
I0_IA <- sum(data$Location == 0 & data$Infected.by == "Started")

# Population sizes — identical to ODE script
N0  <- 42
N_A <- 25
N_E <- 17
S0_A <- N_A - I0_IA
S0_E <- N_E - I0_IE

# Initial state vector — compartment names and order match pop.ssiirr in ODE
y0 <- c(
  SE = S0_E,
  SA = S0_A,
  IA = I0_IA,
  IE = I0_IE,
  RA = 0,
  RE = 0
)

# ---------------------------------------------------------------------------
# Step 2: Parameters — identical to the ODE script
# ---------------------------------------------------------------------------

betaMixprior    <- 4.82234
betaNightAprior <- 0.12757
betaNightEprior <- 0       # Empire has no night-only transmission in the ODE

params <- c(
  betaMix    = betaMixprior,
  betaNightA = betaNightAprior,
  betaNightE = betaNightEprior,
  gammaA     = 1,
  gammaE     = 1,
  N0         = N0,
  N_A        = N_A,
  N_E        = N_E
)

# Daytime window (fraction of day) — matches ODE exactly
DAYTIME_START <- 0.3328
DAYTIME_END   <- 0.77071

# ---------------------------------------------------------------------------
# Step 3: Single Gillespie step — six possible events, time-dependent rates
# ---------------------------------------------------------------------------

event_ssiirr <- function(time, SA, SE, IA, IE, RA, RE, params, t_end) {
  with(as.list(params), {
    is_daytime   <- (time >= floor(time) + 0.3328) & (time <= floor(time) + 0.77071)
    
    # Dynamic betas
    beta_M  <- ifelse(is_daytime, betaMixprior, 0)
    beta_NA <- ifelse(is_daytime, 0, betaNightAprior)
    beta_NE <- ifelse(is_daytime, 0, betaNightEprior)
    
    inf_mix <- beta_M  * (IA + IE) / N0
    inf_NA  <- beta_NA * IA / N_A
    inf_NE  <- beta_NE * IE / N_E
    
    # Six event rates — exactly the terms in the ODE derivatives
    rates <- c(
      infect_A_mix   = inf_mix * SA,          # daytime AIMS infection
      infect_A_night = inf_NA  * SA,           # nighttime AIMS-only infection
      infect_E_mix   = inf_mix * SE,           # daytime Empire infection
      infect_E_night = inf_NE  * SE,           # nighttime Empire-only infection
      recover_A      = gammaA  * IA,
      recover_E      = gammaE  * IE
    )
    
    total_rate <- sum(rates)
    
    new_infection <- 0
    
    if (total_rate == 0) {
      # Epidemic is over — jump to end of simulation
      return(data.frame(
        time = t_end,
        SA = SA, SE = SE, IA = IA, IE = IE, RA = RA, RE = RE,
        new_infection = 0
      ))
    }
    
    # Time to next event (exponential inter-arrival)
    event_time <- time + rexp(1, total_rate)
    
    # Which event fires?
    event_type <- sample(names(rates), 1, prob = rates / total_rate)
    
    switch(event_type,
           "infect_A_mix"   = { SA <- SA - 1; IA <- IA + 1; new_infection <- 1 },
           "infect_A_night" = { SA <- SA - 1; IA <- IA + 1; new_infection <- 1 },
           "infect_E_mix"   = { SE <- SE - 1; IE <- IE + 1; new_infection <- 1 },
           "infect_E_night" = { SE <- SE - 1; IE <- IE + 1; new_infection <- 1 },
           "recover_A"      = { IA <- IA - 1; RA <- RA + 1 },
           "recover_E"      = { IE <- IE - 1; RE <- RE + 1 }
    )
    
    return(data.frame(
      time = event_time,
      SA = SA, SE = SE, IA = IA, IE = IE, RA = RA, RE = RE,
      new_infection = new_infection
    ))
  })
}

# ---------------------------------------------------------------------------
# Step 4: Full Gillespie simulation from t = 0 to t_end
#         Includes a CI (cumulative incidence) column to mirror ssiirr_CI
# ---------------------------------------------------------------------------

simulate_ssiirr <- function(t_end, y, params) {
  with(as.list(y), {
    
    CI <- 0
    
    ts <- data.frame(
      time = 0,
      SA = SA, SE = SE, IA = IA, IE = IE, RA = RA, RE = RE,
      CI = CI
    )
    current <- ts[1, ]
    
    while (current$time < t_end) {
      
      nxt <- event_ssiirr(
        time = current$time,
        SA   = current$SA,
        SE   = current$SE,
        IA   = current$IA,
        IE   = current$IE,
        RA   = current$RA,
        RE   = current$RE,
        params = params,
        t_end  = t_end
      )
      
      # Stop if next event falls beyond t_end
      if (nxt$time >= t_end) break
      
      CI <- CI + nxt$new_infection
      
      current <- data.frame(
        time = nxt$time,
        SA   = nxt$SA, SE = nxt$SE,
        IA   = nxt$IA, IE = nxt$IE,
        RA   = nxt$RA, RE = nxt$RE,
        CI   = CI
      )
      ts <- rbind(ts, current)
    }
    
    # Final row at exactly t_end for clean downstream indexing
    ts <- rbind(ts, data.frame(
      time = t_end,
      SA   = current$SA, SE = current$SE,
      IA   = current$IA, IE = current$IE,
      RA   = current$RA, RE = current$RE,
      CI   = CI
    ))
    
    return(ts)
  })
}

# ---------------------------------------------------------------------------
# Step 5: Run one realisation
# ---------------------------------------------------------------------------

final_time <- 10     # days — same horizon as the ODE script
ts1 <- simulate_ssiirr(final_time, y0, params)

cat("\nGillespie trajectory (last 6 rows):\n"); print(tail(ts1))

# ---------------------------------------------------------------------------
# Global Aesthetic Mappings (Styles & Color Definitions)
# ---------------------------------------------------------------------------

pop_colours <- c(
  SA = "lightblue4", IA = "lightblue4", RA = "lightblue4",
  SE = "orange",     IE = "orange",     RE = "orange"
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
  SA = "dotted",  SE = "dotted",
  IA = "solid",   IE = "solid",
  RA = "dotdash", RE = "dotdash"
)

# ---------------------------------------------------------------------------
# Step 6a: Plot all six compartments (mirrors ODE compartment plot)
# ---------------------------------------------------------------------------

ts1_long <- ts1 %>%
  pivot_longer(cols      = c(SA, SE, IA, IE, RA, RE),
               names_to  = "Compartment",
               values_to = "Count") %>%
  mutate(Compartment = factor(Compartment,
                              levels = c("SA", "SE", "IA", "IE", "RA", "RE")))

Combination <- ggplot(ts1_long, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_step(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since index case (15 June 2026)") +
  ylab("Number of individuals") +
  ggtitle("Interacting SIR Gillespie Model — (single run)",
          subtitle = paste0("Transmission Day = ", round(params["betaMix"], 1),
                            "  |  Trans. Night AIMS = ", round(params["betaNightA"], 1),
                            "  |  Trans. Night Empire = ", round(params["betaNightE"], 1))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13))

ggsave("Combination_Gillespie.png", Combination, width = 8, height = 6, dpi = 300)



# ---------------------------------------------------------------------------
# Step 6b: AIMS compartments only (mirrors ODE Empire plot)
# ---------------------------------------------------------------------------


ts1_long_A <- ts1 %>%
  pivot_longer(cols = c(SA, IA, RA),
               names_to  = "Compartment",
               values_to = "Count")

A <- ggplot(ts1_long_A, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_step(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since 15 June 2026") +
  ylab("Number of individuals") +
  ggtitle("AIMS (Population A) — Interacting SIR Gillespie Model (single run)",
          subtitle = paste0("Trans. Day = ", round(params["betaMix"], 1),
                            " | Trans. Night = ", round(params["betaNightA"], 1))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13))

ggsave("AIMS_Gillespie.png", A, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Step 6c: Empire compartments only (mirrors ODE Empire plot)
# ---------------------------------------------------------------------------

ts1_long_E <- ts1 %>%
  pivot_longer(cols = c(SE, IE, RE),
               names_to  = "Compartment",
               values_to = "Count")

E<-ggplot(ts1_long_E, aes(x = time, y = Count, colour = Compartment, linetype = Compartment)) +
  geom_step(linewidth = 1.1) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_linetype_manual(values = state_linetypes, labels = pop_labels) +
  xlab("Days since 15 June 2026") +
  ylab("Number of individuals") +
  ggtitle("Empire (Population E) — Interacting SIR Gillespie Model (single run)",
          subtitle = paste0("Trans. Day = ", round(params["betaMix"], 1),
                            " | Trans. Night = ", round(params["betaNightE"], 1))) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13))

ggsave("Empire_Gillespie.png", E, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# Step 7: Multiple realisations — median + 95 % envelope (Original Simulation Method)
# ---------------------------------------------------------------------------

n_sims   <- 1000
time_out <- seq(0, final_time, by = 0.02083)   # ~30-min grid, matches ODE step

# Helper: look up compartment value at each grid point (last-value carry-forward)
get_compartment_at <- function(traj, tgrid, col) {
  sapply(tgrid, function(t) {
    idx <- max(which(traj$time <= t))
    traj[[col]][idx]
  })
}

sim_list <- lapply(seq_len(n_sims), function(i) {
  traj <- simulate_ssiirr(final_time, y0, params)
  data.frame(
    sim  = i,
    time = time_out,
    SA   = get_compartment_at(traj, time_out, "SA"),
    SE   = get_compartment_at(traj, time_out, "SE"),
    IA   = get_compartment_at(traj, time_out, "IA"),
    IE   = get_compartment_at(traj, time_out, "IE"),
    RA   = get_compartment_at(traj, time_out, "RA"),
    RE   = get_compartment_at(traj, time_out, "RE"),
    CI   = get_compartment_at(traj, time_out, "CI")
  )
})

sim_all <- bind_rows(sim_list)

# Summarise per time point
ensemble <- sim_all %>%
  group_by(time) %>%
  summarise(
    IA_med = median(IA), IA_lo = quantile(IA, 0.025), IA_hi = quantile(IA, 0.975),
    IE_med = median(IE), IE_lo = quantile(IE, 0.025), IE_hi = quantile(IE, 0.975),
    SA_med = median(SA), SA_lo = quantile(SA, 0.025), SA_hi = quantile(SA, 0.975),
    SE_med = median(SE), SE_lo = quantile(SE, 0.025), SE_hi = quantile(SE, 0.975),
    RA_med = median(RA), RA_lo = quantile(RA, 0.025), RA_hi = quantile(RA, 0.975),
    RE_med = median(RE), RE_lo = quantile(RE, 0.025), RE_hi = quantile(RE, 0.975),
    .groups = "drop"
  )

# --- Plot: total infectious (IA + IE) ensemble (Saved using plotTest dimensions) ---
plotEnsembleTotal <- ggplot(ensemble, aes(x = time)) +
  geom_ribbon(aes(ymin = IA_lo + IE_lo, ymax = IA_hi + IE_hi),
              fill = "tomato", alpha = 0.25) +
  geom_line(aes(y = IA_med + IE_med, colour = "Gillespie median (IA + IE)"),
            linewidth = 1.1) +
  scale_colour_manual(values = c("Gillespie median (IA + IE)" = "red")) +
  xlab("Days since index case") +
  ylab("Active infectious individuals") +
  ggtitle(
    paste0("Interacting SIR Gillespie — ", n_sims, " runs (total infectious)"),
    subtitle = paste0("Shaded: 95 % simulation envelope  |  ",
                      "Trans. Day = ", round(params["betaMix"], 1),
                      "  |  Trans. Night AIMS = ", round(params["betaNightA"], 1))
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13))

ggsave("ensemble_total_infectious.png", plotEnsembleTotal, width = 8, height = 6, dpi = 300)

# --- Plot: AIMS vs Empire infectious separately (Saved using plotTest dimensions) ---
ensemble_long <- ensemble %>%
  select(time, IA_med, IE_med, IA_lo, IE_lo, IA_hi, IE_hi) %>%
  pivot_longer(
    cols      = -time,
    names_to  = c("Compartment", ".value"),
    names_sep = "_"
  )

plotEnsembleSplit <- ggplot(ensemble_long, aes(x = time, colour = Compartment, fill = Compartment)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(aes(y = med), linewidth = 1.1) +
  scale_colour_manual(values = c(IA = "orange", IE = "lightblue4"),
                      labels = c(IA = "AIMS infectious",
                                 IE = "Empire infectious")) +
  scale_fill_manual(values   = c(IA = "orange", IE = "lightblue4"),
                    labels   = c(IA = "AIMS infectious",
                                 IE = "Empire infectious")) +
  xlab("Days since index case") +
  ylab("Infectious individuals") +
  ggtitle(
    paste0("Interacting SIR Gillespie — ", n_sims, " runs (by sub-population)"),
    subtitle = "Shaded: 95 % simulation envelope"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        axis.title = element_text(size = 15),
        axis.text  = element_text(size = 13))

ggsave("ensemble_split_infectious.png", plotEnsembleSplit, width = 8, height = 6, dpi = 300)



