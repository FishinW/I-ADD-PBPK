# =============================================================================
# File: 07_simulate_rat_kidney_iadd_mc.R
# Purpose: Simulate the rat kidney I-ADD scenario and Monte Carlo prediction intervals.
# Inputs: Kidney I-ADD model and 7-day device data.
# Outputs: Observed-versus-predicted summaries and Monte Carlo plots.
# Run from the repository root after installing the packages listed in README.md.
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FME)
library(minpack.lm)
library(gridExtra)
library(patchwork)
library(data.table)
library(purrr)

###############################################################################
### Kidney I-ADD Model
###############################################################################

source(file.path('R', 'models', '05_model_rat_kidney_iadd.R'))

mod_60Days <- mcode("DexPBPK", DexPBPK) %>%
  param(
    PR          = 0.13,
    PABRC_VAS   = 0.11,
    PABRC_INT   = 0.375,
    fu          = 0.175,
    fucel_brain = 0.24,
    fuint_brain = 0.29,
    BP          = 0.725,
    CLint       = 0.04
  )
out_60Days <- mod_60Days %>%
  mrgsim(end = 1000, delta = 0.5) %>%
  as.data.frame()
predicted <- out_60Days %>%
  select(ID, CPLASMA_out, CLIVER_out, CKID_out, time) %>%
  rename(
    Time     = time,
    Study_ID = ID,
    Plasma   = CPLASMA_out,
    Liver    = CLIVER_out,
    Kidney   = CKID_out
  ) %>%
  pivot_longer(cols = c(Plasma, Liver, Kidney),
               names_to = "Compartments", values_to = "Predicted")

observed <- read.csv("data/animal/Kidney-7Days.csv",
                     check.names = FALSE, strip.white = TRUE) %>%
  pivot_longer(cols = c(Rat1, Rat2, Rat3),
               names_to = "Subject", values_to = "Observed", values_drop_na = TRUE) %>%
  mutate(Observed = Observed / 1000,
         Time = 24 * 7) %>%
  mutate(Compartments = gsub("Left Kidney",  "Kidney", Compartments),
         Compartments = gsub("Right Kidney", "Kidney", Compartments))

obs_sum_kid <- observed %>%
  group_by(Compartments, Time) %>%
  summarise(
    Mean = mean(Observed, na.rm = TRUE),
    SD   = sd(Observed,   na.rm = TRUE),
    N    = sum(!is.na(Observed)),
    .groups = "drop"
  )

obs_time_kid <- unique(obs_sum_kid$Time)

p_study_kid <- ggplot() +
  geom_line(data = predicted, aes(x = Time, y = Predicted, color = "Predicted"),
            linewidth = 1.3, alpha = 0.95) +
  geom_vline(xintercept = obs_time_kid, linetype = "dashed",
             linewidth = 0.7, color = "#1f78b4", alpha = 0.65) +
  geom_errorbar(data = obs_sum_kid,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-5), ymax = Mean + SD,
                    color = "Observed (mean \u00b1 SD)"),
                width = 12, linewidth = 1.0) +
  geom_point(data = obs_sum_kid,
             aes(x = Time, y = Mean, color = "Observed (mean \u00b1 SD)"),
             size = 3.8, alpha = 0.95) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(limits = c(0.0001, 10),
                breaks = c(10, 1, 0.1, 0.01, 0.001, 0.0001),
                labels = c("10", "1", "0.1", "0.01", "0.001", "0.0001")) +
  scale_color_manual(values = c("Predicted" = "#1f78b4",
                                "Observed (mean \u00b1 SD)" = "#ff7f00")) +
  facet_wrap(~ Compartments, scales = "fixed", ncol = 3) +
  labs(title = "PBPK-simulated DEX exposure following I-ADD (Kidney)",
       x = "Time (h)", y = "DEX concentration (\u00b5g/mL)", color = NULL) +
  theme_bw(base_size = 13) +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray85", linewidth = 0.35),
    panel.grid.minor = element_line(color = "gray92", linewidth = 0.2),
    axis.text        = element_text(size = 12, colour = "black", face = "bold"),
    axis.title       = element_text(size = 15, colour = "black", face = "bold"),
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.background = element_rect(colour = "black", linewidth = 1.1, fill = "gray90"),
    strip.text       = element_text(face = "bold", size = 13),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.key.width = unit(1.2, "cm"),
    panel.spacing    = unit(1.2, "lines"),
    plot.margin      = margin(8, 8, 8, 8)
  )

print(p_study_kid)

pred_at_obs <- predicted %>%
  semi_join(obs_sum_kid %>% select(Compartments, Time), by = c("Compartments", "Time")) %>%
  group_by(Compartments, Time) %>%
  summarise(Predicted = mean(Predicted, na.rm = TRUE), .groups = "drop")

compare_tbl <- obs_sum_kid %>%
  select(Compartments, Time, Observed = Mean) %>%
  left_join(pred_at_obs, by = c("Compartments", "Time")) %>%
  mutate(`Fold error` = Predicted / Observed) %>%
  mutate(Compartments = recode(Compartments,
                               "Kidney" = "I-ADD target tissue (Kidney)",
                               "Liver"  = "Liver",
                               "Plasma" = "Plasma"))

summary_display <- compare_tbl %>%
  select(Compartments, Observed, Predicted, `Fold error`) %>%
  pivot_longer(cols = c(Observed, Predicted, `Fold error`),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed", "Predicted", "Fold error")),
    Compartments = factor(Compartments,
                          levels = c("I-ADD target tissue (Kidney)", "Liver", "Plasma"))
  ) %>%
  pivot_wider(names_from = Compartments, values_from = Value) %>%
  arrange(Metric)

summary_display

###############################################################################
### Monte Carlo Simulation - Kidney I-ADD Model
###############################################################################
set.seed(123)

# Parameter table (lognormal only, all CV = 0.2)
mc_params_kid <- data.frame(
  name = c("PK",    "PR",    "PABRC_VAS", "PABRC_INT", "fu",    "BP",    "CLint", "KREL"),
  mean = c(1.519,   0.134,   0.11,       0.375,       0.175,   0.725,   0.040,   0.0126),
  CV   = c(0.2,     0.2,     0.2,         0.2,         0.2,     0.2,     0.2,     0.2)
)

# Generate MC parameters (lognormal, truncated at 2.5th-97.5th percentiles)
generate_mc_pars_kid <- function(N = 1000, params_table) {
  
  mc <- data.frame(MC_ID = 1:N)
  
  for (i in seq_len(nrow(params_table))) {
    pname   <- params_table$name[i]
    mu_nat  <- params_table$mean[i]
    cv      <- params_table$CV[i]
    sdlog   <- sqrt(log(1 + cv^2))
    meanlog <- log(mu_nat) - 0.5 * sdlog^2
    lower_bound <- qlnorm(0.025, meanlog = meanlog, sdlog = sdlog)
    upper_bound <- qlnorm(0.975, meanlog = meanlog, sdlog = sdlog)
    samples <- rlnorm(N * 10, meanlog = meanlog, sdlog = sdlog)
    samples <- samples[samples >= lower_bound & samples <= upper_bound]
    while (length(samples) < N) {
      extra   <- rlnorm(N * 10, meanlog = meanlog, sdlog = sdlog)
      samples <- c(samples, extra[extra >= lower_bound & extra <= upper_bound])
    }
    mc[[pname]] <- samples[1:N]
  }
  
  mc
}

mc_pars_kid <- generate_mc_pars_kid(N = 1000, params_table = mc_params_kid)

cat("MC parameter summary (mean | CV):\n")
print(mc_pars_kid %>% select(-MC_ID) %>%
        summarise(across(everything(), list(mean = mean, cv = ~sd(.) / mean(.)))))

# One-subject simulation function
run_one_mc_kid <- function(par_row, end = 1000, delta = 0.5) {
  pars_nat <- as.list(par_row[, -1])   # drop MC_ID
  out <- mod_60Days %>%
    param(pars_nat) %>%
    mrgsim(end = end, delta = delta) %>%
    as.data.frame()
  out$MC_ID <- par_row$MC_ID
  out
}

# Summarise helper
summarise_mc_kid <- function(sim_df, compartments_map) {
  sim_df %>%
    select(time, MC_ID, all_of(names(compartments_map))) %>%
    pivot_longer(cols = all_of(names(compartments_map)),
                 names_to = "Compartments", values_to = "value") %>%
    mutate(Compartments = recode(Compartments, !!!compartments_map)) %>%
    group_by(Compartments, time) %>%
    summarise(
      P50 = median(value, na.rm = TRUE),
      P05 = quantile(value, 0.05, na.rm = TRUE),
      P95 = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Time = time)
}

# Compartment mapping and colors
compartments_map_kid <- c(
  CPLASMA_out = "Plasma",
  CLIVER_out  = "Liver",
  CKID_out    = "Kidney"
)

organ_colors_kid <- c(
  "Plasma" = "black",
  "Liver"  = "forestgreen",
  "Kidney" = "firebrick"
)

# Run MC
cat("\nRunning Monte Carlo for Kidney I-ADD model (N=1000)...\n")

sim_mc_kid <- map_dfr(
  split(mc_pars_kid, mc_pars_kid$MC_ID),
  run_one_mc_kid,
  end   = 1000,
  delta = 0.5
)

mc_summary_kid <- summarise_mc_kid(sim_mc_kid, compartments_map_kid)

# Plot - Kidney MC
mc_summary_kid <- mc_summary_kid %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_kid)))
obs_sum_kid_plot <- obs_sum_kid %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_kid)))

obs_time_mc_kid <- unique(obs_sum_kid_plot$Time)

p_mc_kid <- ggplot() +
  geom_ribbon(data = mc_summary_kid,
              aes(x = Time, ymin = pmax(P05, 1e-4), ymax = P95, fill = Compartments),
              alpha = 0.3) +
  geom_line(data = mc_summary_kid,
            aes(x = Time, y = P50, color = Compartments),
            linewidth = 1) +
  geom_vline(xintercept = obs_time_mc_kid,
             linetype = "dashed", color = "#ff7f00",
             linewidth = 0.7, alpha = 0.6) +
  geom_errorbar(data = obs_sum_kid_plot,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-4), ymax = Mean + SD),
                width = 6, linewidth = 1, color = "#ff7f00") +
  geom_point(data = obs_sum_kid_plot,
             aes(x = Time, y = Mean),
             size = 4, color = "#ff7f00") +
  facet_wrap(~ Compartments, scales = "fixed", ncol = 3) +
  scale_color_manual(values = organ_colors_kid) +
  scale_fill_manual(values = organ_colors_kid) +
  scale_x_continuous(limits = c(0, 250),
                     breaks = seq(0, 250, 50),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(limits = c(1e-4, 10),
                breaks = c(10, 1, 0.1, 0.01, 0.001, 0.0001),
                labels = c("10", "1", "0.1", "0.01", "0.001", "0.0001")) +
  theme_bw(base_size = 13) +
  labs(title = "Monte Carlo Simulation (5th\u201395th Percentile) vs Observed (I-ADD Kidney)",
       x = "Time (h)", y = "DEX concentration (\u00b5g/mL)") +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "gray90"),
    strip.text       = element_text(face = "bold", size = 13),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = 12, colour = "black", face = "bold"),
    axis.title       = element_text(size = 14, colour = "black", face = "bold"),
    panel.spacing    = unit(0.8, "lines")
  )

print(p_mc_kid)


# =============================================================================
# Save kidney I-ADD simulation outputs
# =============================================================================
write.csv(mc_pars_kid, file.path("results", "tables", "kidney_iadd_monte_carlo_parameters.csv"), row.names = FALSE)
if (exists("mc_summary_kid")) write.csv(mc_summary_kid, file.path("results", "tables", "kidney_iadd_mc_summary.csv"), row.names = FALSE)

ggsave(file.path("results", "figures", "kidney_iadd_deterministic.png"),
       p_study_kid, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
ggsave(file.path("results", "figures", "kidney_iadd_monte_carlo.png"),
       p_mc_kid, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
