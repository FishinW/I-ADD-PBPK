# =============================================================================
# File: 06_simulate_rat_brain_iadd_mc.R
# Purpose: Simulate rat brain I-ADD scenarios and Monte Carlo prediction intervals.
# Inputs: Brain I-ADD models and 7-day device data.
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
### 30-Day Brain I-ADD Model
###############################################################################

source(file.path('R', 'models', '03_model_rat_brain_iadd_30day.R'))

mod_30Days <- mcode("DexPBPK", DexPBPK) %>%
  param(
    DOSE_IADD   = 0.37,
    PR          = 0.13,
    PABRC_VAS   = 0.11,
    PABRC_INT   = 0.375,
    fu          = 0.175,
    fucel_brain = 0.24,
    fuint_brain = 0.29,
    BP          = 0.725,
    CLint       = 0.040
  )
out_30Days <- mod_30Days %>%
  mrgsim(end = 1000, delta = 0.5) %>%
  as.data.frame()
predicted <- out_30Days %>%
  select(ID, CPLASMA_out, CLIVER_out, CBR_out, time) %>%
  rename(
    Time     = time,
    Study_ID = ID,
    Plasma   = CPLASMA_out,
    Liver    = CLIVER_out,
    Brain    = CBR_out
  ) %>%
  pivot_longer(cols = c(Plasma, Liver, Brain),
               names_to = "Compartments", values_to = "Predicted")

observed <- read.csv("data/animal/Brain-7Days.csv",
                     check.names = FALSE, strip.white = TRUE) %>%
  pivot_longer(cols = c(Rat1, Rat2, Rat3, Rat4, Rat5),
               names_to = "Subject", values_to = "Observed", values_drop_na = TRUE) %>%
  mutate(Compartments = gsub("Serum", "Plasma", Compartments),
         Observed = Observed / 1000,
         Time = 24 * 7)

obs_sum_30 <- observed %>%
  group_by(Compartments, Time) %>%
  summarise(
    Mean = mean(Observed, na.rm = TRUE),
    SD   = sd(Observed,   na.rm = TRUE),
    N    = sum(!is.na(Observed)),
    .groups = "drop"
  )

obs_time <- unique(obs_sum_30$Time)

p_study_30 <- ggplot() +
  geom_line(data = predicted, aes(x = Time, y = Predicted, color = "Predicted"),
            linewidth = 1.3, alpha = 0.95) +
  geom_vline(xintercept = obs_time, linetype = "dashed",
             linewidth = 0.7, color = "#1f78b4", alpha = 0.65) +
  geom_errorbar(data = obs_sum_30,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-5), ymax = Mean + SD,
                    color = "Observed (mean \u00b1 SD)"),
                width = 12, linewidth = 1.0) +
  geom_point(data = obs_sum_30,
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
  labs(title = "PBPK-simulated DEX exposure following I-ADD (30-Day)",
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

print(p_study_30)

pred_at_obs <- predicted %>%
  semi_join(obs_sum_30 %>% select(Compartments, Time), by = c("Compartments", "Time")) %>%
  group_by(Compartments, Time) %>%
  summarise(Predicted = mean(Predicted, na.rm = TRUE), .groups = "drop")

compare_tbl <- obs_sum_30 %>%
  select(Compartments, Time, Observed = Mean) %>%
  left_join(pred_at_obs, by = c("Compartments", "Time")) %>%
  mutate(`Fold error` = Predicted / Observed) %>%
  mutate(Compartments = recode(Compartments,
                               "Brain"  = "I-ADD target tissue (Brain)",
                               "Liver"  = "Liver",
                               "Plasma" = "Plasma"))

summary_display <- compare_tbl %>%
  select(Compartments, Observed, Predicted, `Fold error`) %>%
  pivot_longer(cols = c(Observed, Predicted, `Fold error`),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed", "Predicted", "Fold error")),
    Compartments = factor(Compartments,
                          levels = c("I-ADD target tissue (Brain)", "Liver", "Plasma"))
  ) %>%
  pivot_wider(names_from = Compartments, values_from = Value) %>%
  arrange(Metric)

summary_display


###############################################################################
### Monte Carlo Simulation - 30-Day Brain I-ADD Model
###############################################################################
set.seed(123)

# Parameter table
mc_params_brain <- list(
  normal = data.frame(
    name = c("BW",  "QCC", "QBRC"),
    mean = c(0.25,  14,    0.02),
    CV   = c(0.2,   0.2,   0.2)
  ),
  lognormal = data.frame(
    name = c("PR",  "PBR", "PABRC_VAS", "PABRC_INT", "PABRC_OUT",
             "PABRC_IN", "fucel_brain", "fu",   "BP",   "CLint"),
    mean = c(0.134, 1.29,  0.11,       0.375,       1.0,
             0.24,       0.24,         0.175,  0.725,  0.041),
    CV   = c(0.2,   0.2,   0.2,         0.2,         0.2,
             0.2,        0.2,          0.2,    0.2,    0.2)
  )
)

# Generate MC parameters (truncated at 2.5th-97.5th percentiles)
generate_mc_pars_brain <- function(N = 500, params = mc_params_brain) {
  
  mc <- data.frame(MC_ID = 1:N)
  
  for (i in seq_len(nrow(params$normal))) {
    pname <- params$normal$name[i]
    mu    <- params$normal$mean[i]
    sigma <- mu * params$normal$CV[i]
    lower_bound <- qnorm(0.025, mean = mu, sd = sigma)
    upper_bound <- qnorm(0.975, mean = mu, sd = sigma)
    samples <- rnorm(N * 10, mean = mu, sd = sigma)
    samples <- samples[samples >= lower_bound & samples <= upper_bound]
    while (length(samples) < N) {
      extra <- rnorm(N * 10, mean = mu, sd = sigma)
      samples <- c(samples, extra[extra >= lower_bound & extra <= upper_bound])
    }
    mc[[pname]] <- samples[1:N]
  }
  
  for (i in seq_len(nrow(params$lognormal))) {
    pname   <- params$lognormal$name[i]
    mu_nat  <- params$lognormal$mean[i]
    sdlog   <- sqrt(log(1 + params$lognormal$CV[i]^2))
    meanlog <- log(mu_nat) - 0.5 * sdlog^2
    lower_bound <- qlnorm(0.025, meanlog = meanlog, sdlog = sdlog)
    upper_bound <- qlnorm(0.975, meanlog = meanlog, sdlog = sdlog)
    samples <- rlnorm(N * 10, meanlog = meanlog, sdlog = sdlog)
    samples <- samples[samples >= lower_bound & samples <= upper_bound]
    while (length(samples) < N) {
      extra <- rlnorm(N * 10, meanlog = meanlog, sdlog = sdlog)
      samples <- c(samples, extra[extra >= lower_bound & extra <= upper_bound])
    }
    mc[[pname]] <- samples[1:N]
  }
  
  mc
}

mc_pars_brain <- generate_mc_pars_brain(N = 1000, params = mc_params_brain)

cat("MC parameter summary (mean | CV):\n")
print(mc_pars_brain %>% select(-MC_ID) %>%
        summarise(across(everything(), list(mean = mean, cv = ~sd(.) / mean(.)))))

# One-subject simulation function
run_one_mc_brain <- function(par_row, model, end = 1000, delta = 0.5) {
  pars_nat <- as.list(par_row[, -1])
  out <- model %>%
    param(pars_nat) %>%
    mrgsim(end = end, delta = delta) %>%
    as.data.frame()
  out$MC_ID <- par_row$MC_ID
  out
}

# Summarise helper
summarise_mc_brain <- function(sim_df, compartments_map) {
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

# Compartment mappings
compartments_map_30 <- c(CPLASMA_out = "Plasma", CLIVER_out = "Liver", CBR_out = "Brain")
compartments_map_60 <- c(CPLASMA_out = "Plasma", CLIVER_out = "Liver",
                         CKID_out = "Kidney", CBR_out = "Brain")

organ_colors_brain <- c("Plasma" = "black", "Liver" = "forestgreen", "Brain" = "steelblue")
organ_colors_brain60 <- c("Plasma" = "black", "Liver" = "forestgreen",
                          "Kidney" = "firebrick", "Brain" = "steelblue")

# Run MC - 30-Day
cat("\nRunning Monte Carlo for 30-Day Brain I-ADD model (N=1000)...\n")

sim_mc_brain30 <- map_dfr(
  split(mc_pars_brain, mc_pars_brain$MC_ID),
  run_one_mc_brain,
  model = mod_30Days,
  end   = 1000,
  delta = 0.5
)

mc_summary_brain30 <- summarise_mc_brain(sim_mc_brain30, compartments_map_30)

# Plot - 30-Day
mc_summary_brain30 <- mc_summary_brain30 %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_brain)))
obs_sum_30_plot <- obs_sum_30 %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_brain)))

obs_time_30 <- unique(obs_sum_30_plot$Time)

p_mc_brain30 <- ggplot() +
  geom_ribbon(data = mc_summary_brain30,
              aes(x = Time, ymin = pmax(P05, 1e-4), ymax = P95, fill = Compartments),
              alpha = 0.3) +
  geom_line(data = mc_summary_brain30,
            aes(x = Time, y = P50, color = Compartments),
            linewidth = 1) +
  geom_vline(xintercept = obs_time_30,
             linetype = "dashed", color = "#ff7f00",
             linewidth = 0.7, alpha = 0.6) +
  geom_errorbar(data = obs_sum_30_plot,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-4), ymax = Mean + SD),
                width = 6, linewidth = 1, color = "#ff7f00") +
  geom_point(data = obs_sum_30_plot,
             aes(x = Time, y = Mean),
             size = 4, color = "#ff7f00") +
  facet_wrap(~ Compartments, scales = "fixed", ncol = 3) +
  scale_color_manual(values = organ_colors_brain) +
  scale_fill_manual(values = organ_colors_brain) +
  scale_x_continuous(limits = c(0, 250),
                     breaks = seq(0, 250, 50),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(limits = c(1e-4, 10),
                breaks = c(10, 1, 0.1, 0.01, 0.001, 0.0001),
                labels = c("10", "1", "0.1", "0.01", "0.001", "0.0001")) +
  theme_bw(base_size = 13) +
  labs(title = "Monte Carlo Simulation (5th\u201395th Percentile) vs Observed (Brain I-ADD 30-Day)",
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

print(p_mc_brain30)

###############################################################################
### 60-Day Brain I-ADD Model
###############################################################################

source(file.path('R', 'models', '04_model_rat_brain_iadd_60day.R'))

mod_60Days <- mcode("DexPBPK", DexPBPK) %>%
  param(
    PR          = 0.13,
    PABRC_VAS   = 0.15,
    PABRC_INT   = 0.30,
    fu          = 0.175,
    fucel_brain = 0.24,
    fuint_brain = 0.29,
    BP          = 0.725,
    CLint       = 0.040,
    FMAX_REL    = 0.928,
    KREL        = 0.0133
  )
out_60Days <- mod_60Days %>%
  mrgsim(end = 1000, delta = 0.5) %>%
  as.data.frame()
predicted <- out_60Days %>%
  select(ID, CPLASMA_out, CLIVER_out, CKID_out, CBR_out, time) %>%
  rename(
    Time     = time,
    Study_ID = ID,
    Plasma   = CPLASMA_out,
    Liver    = CLIVER_out,
    Kidney   = CKID_out,
    Brain    = CBR_out
  ) %>%
  pivot_longer(cols = c(Plasma, Liver, Brain, Kidney),
               names_to = "Compartments", values_to = "Predicted")

observed <- read.csv("data/animal/Brain-7Days-60Days.csv",
                     check.names = FALSE, strip.white = TRUE) %>%
  pivot_longer(cols = c(Rat1, Rat2, Rat3, Rat4),
               names_to = "Subject", values_to = "Observed", values_drop_na = TRUE) %>%
  mutate(Compartments = gsub("Serum", "Plasma", Compartments),
         Observed = Observed / 1000,
         Time = 24 * 7)

obs_sum_60 <- observed %>%
  group_by(Compartments, Time) %>%
  summarise(
    Mean = mean(Observed, na.rm = TRUE),
    SD   = sd(Observed,   na.rm = TRUE),
    N    = sum(!is.na(Observed)),
    .groups = "drop"
  )

obs_time <- unique(obs_sum_60$Time)

p_study_60 <- ggplot() +
  geom_line(data = predicted, aes(x = Time, y = Predicted, color = "Predicted"),
            linewidth = 1.3, alpha = 0.95) +
  geom_vline(xintercept = obs_time, linetype = "dashed",
             linewidth = 0.7, color = "#1f78b4", alpha = 0.65) +
  geom_errorbar(data = obs_sum_60,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-5), ymax = Mean + SD,
                    color = "Observed (mean \u00b1 SD)"),
                width = 12, linewidth = 1.0) +
  geom_point(data = obs_sum_60,
             aes(x = Time, y = Mean, color = "Observed (mean \u00b1 SD)"),
             size = 3.8, alpha = 0.95) +
  scale_x_continuous(limits = c(0, 500), breaks = seq(0, 500, 100),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(limits = c(0.0001, 10),
                breaks = c(10, 1, 0.1, 0.01, 0.001, 0.0001),
                labels = c("10", "1", "0.1", "0.01", "0.001", "0.0001")) +
  scale_color_manual(values = c("Predicted" = "#1f78b4",
                                "Observed (mean \u00b1 SD)" = "#ff7f00")) +
  facet_wrap(~ Compartments, scales = "fixed", ncol = 4) +
  labs(title = "PBPK-simulated DEX exposure following I-ADD (60-Day)",
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

print(p_study_60)

pred_at_obs <- predicted %>%
  semi_join(obs_sum_60 %>% select(Compartments, Time), by = c("Compartments", "Time")) %>%
  group_by(Compartments, Time) %>%
  summarise(Predicted = mean(Predicted, na.rm = TRUE), .groups = "drop")

compare_tbl <- obs_sum_60 %>%
  select(Compartments, Time, Observed = Mean) %>%
  left_join(pred_at_obs, by = c("Compartments", "Time")) %>%
  mutate(`Fold error` = Predicted / Observed) %>%
  mutate(Compartments = recode(Compartments,
                               "Brain"  = "I-ADD target tissue (Brain)",
                               "Kidney" = "Kidney",
                               "Liver"  = "Liver",
                               "Plasma" = "Plasma"))

summary_display <- compare_tbl %>%
  select(Compartments, Observed, Predicted, `Fold error`) %>%
  pivot_longer(cols = c(Observed, Predicted, `Fold error`),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed", "Predicted", "Fold error")),
    Compartments = factor(Compartments,
                          levels = c("I-ADD target tissue (Brain)", "Kidney", "Liver", "Plasma"))
  ) %>%
  pivot_wider(names_from = Compartments, values_from = Value) %>%
  arrange(Metric)

summary_display


###############################################################################
### Monte Carlo Simulation - 60-Day Brain I-ADD Model
###############################################################################

# Run MC - 60-Day (reuse same mc_pars_brain population)
cat("\nRunning Monte Carlo for 60-Day Brain I-ADD model (N=1000)...\n")

sim_mc_brain60 <- map_dfr(
  split(mc_pars_brain, mc_pars_brain$MC_ID),
  run_one_mc_brain,
  model = mod_60Days,
  end   = 1000,
  delta = 0.5
)

mc_summary_brain60 <- summarise_mc_brain(sim_mc_brain60, compartments_map_60)

# Plot - 60-Day
mc_summary_brain60 <- mc_summary_brain60 %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_brain60)))
obs_sum_60_plot <- obs_sum_60 %>%
  mutate(Compartments = factor(Compartments, levels = names(organ_colors_brain60)))

obs_time_60 <- unique(obs_sum_60_plot$Time)

p_mc_brain60 <- ggplot() +
  geom_ribbon(data = mc_summary_brain60,
              aes(x = Time, ymin = pmax(P05, 1e-4), ymax = P95, fill = Compartments),
              alpha = 0.3) +
  geom_line(data = mc_summary_brain60,
            aes(x = Time, y = P50, color = Compartments),
            linewidth = 1) +
  geom_vline(xintercept = obs_time_60,
             linetype = "dashed", color = "#ff7f00",
             linewidth = 0.7, alpha = 0.6) +
  geom_errorbar(data = obs_sum_60_plot,
                aes(x = Time, ymin = pmax(Mean - SD, 1e-4), ymax = Mean + SD),
                width = 6, linewidth = 1, color = "#ff7f00") +
  geom_point(data = obs_sum_60_plot,
             aes(x = Time, y = Mean),
             size = 4, color = "#ff7f00") +
  facet_wrap(~ Compartments, scales = "fixed", ncol = 4) +
  scale_color_manual(values = organ_colors_brain60) +
  scale_fill_manual(values = organ_colors_brain60) +
  scale_x_continuous(limits = c(0, 250),
                     breaks = seq(0, 250, 50),
                     expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_log10(limits = c(1e-4, 10),
                breaks = c(10, 1, 0.1, 0.01, 0.001, 0.0001),
                labels = c("10", "1", "0.1", "0.01", "0.001", "0.0001")) +
  theme_bw(base_size = 13) +
  labs(title = "Monte Carlo Simulation (5th\u201395th Percentile) vs Observed (Brain I-ADD 60-Day)",
       x = "Time (h)", y = "DEX concentration (\u00b5g/mL)") +
  theme(
    legend.position  = "bottom",
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

print(p_mc_brain60)


# =============================================================================
# Save brain I-ADD simulation outputs
# =============================================================================
write.csv(mc_pars_brain, file.path("results", "tables", "brain_iadd_monte_carlo_parameters.csv"), row.names = FALSE)
if (exists("mc_summary_30")) write.csv(mc_summary_30, file.path("results", "tables", "brain_iadd_30day_mc_summary.csv"), row.names = FALSE)
if (exists("mc_summary_60")) write.csv(mc_summary_60, file.path("results", "tables", "brain_iadd_60day_mc_summary.csv"), row.names = FALSE)

ggsave(file.path("results", "figures", "brain_iadd_30day_deterministic.png"),
       p_study_30, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
ggsave(file.path("results", "figures", "brain_iadd_30day_monte_carlo.png"),
       p_mc_brain30, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
ggsave(file.path("results", "figures", "brain_iadd_60day_deterministic.png"),
       p_study_60, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
ggsave(file.path("results", "figures", "brain_iadd_60day_monte_carlo.png"),
       p_mc_brain60, width = 10.5, height = 6.5, units = "in", dpi = 600, bg = "white")
