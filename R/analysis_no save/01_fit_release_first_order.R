# =============================================================================
# File: 01_fit_release_first_order.R
# Purpose: Fit the in vitro I-ADD release profile.
# Input: data/release/release_profile.csv
# Run from the repository root.
# =============================================================================


library(dplyr)
library(tidyr)
library(ggplot2)

setwd("..")

release_data <- read.csv(file.path("data", "release", "release_profile.csv"))

release_30_day_brain <- release_data %>% 
  filter(Device == "30-day in vitro release",
         Organ  == "Brain") 

release_60_day_brain <- release_data %>% 
  filter(Device == "60-day in vitro release",
         Organ  == "Brain") 


release_60_day_kidney <- release_data %>% 
  filter(Device == "60-day in vitro release",
         Organ  == "Kidney")


fit_exp_30day_brain <- nls(
  Propotion ~ Fmax * (1 - exp(-k * Time)),
  data  = release_30_day_brain,
  start = list(Fmax = 0.99, k = 0.01)
)

summary(fit_exp_30day_brain)
coef_30_brain <- coef(fit_exp_30day_brain)

fit_exp_60day_kidney <- nls(
  Propotion ~ Fmax * (1 - exp(-k * Time)),
  data  = release_60_day_kidney,
  start = list(Fmax = 0.99, k = 0.01)
)

summary(fit_exp_60day_kidney)
coef_60_kidney <- coef(fit_exp_60day_kidney)


fit_exp_60day_brain <- nls(
  Propotion ~ Fmax * (1 - exp(-k * Time)),
  data  = release_60_day_brain,
  start = list(Fmax = 0.99, k = 0.01)
)

summary(fit_exp_60day_brain)
coef_60_brain <- coef(fit_exp_60day_brain)

# ---- Plot for Kidney----
Fmax_hat_kidney   <- coef_60_kidney["Fmax"]
k_hat_kidney      <- coef_60_kidney["k"]

# MODIFIED: Simplified formula notation
eq_label_60_kidney <- sprintf(
  "F(t) = %.3f × (1 - exp(-t * %.4f))",
  Fmax_hat_kidney, k_hat_kidney
)


pred_grid_60_kidney <- tibble(
  Time = seq(0, max(release_60_day_kidney$Time), length.out = 300)
)

pred_grid_60_kidney$Pred <- predict(fit_exp_60day_kidney, newdata = pred_grid_60_kidney)



ggplot(release_60_day_kidney, aes(Time, Propotion)) +
  geom_line(data = pred_grid_60_kidney, aes(Time, Pred),
            linewidth = 1.2, color = "blue") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)) +
  geom_point(size = 2) +
  # MODIFIED: Position formula inside the plot area
  annotate(
    "text",
    x = max(release_60_day_kidney$Time) * 0.7, 
    y = 0.1,  
    label = eq_label_60_kidney,
    hjust = 0.5, vjust = 0.5,
    size = 4.5,
    fontface = "italic"
  ) +
  labs(
    x = "Time (h)",
    y = "Fraction released",
    title = "First order release formula: 60-day in vitro release data for kidney I-ADD"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25),
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(colour = "black", linewidth = 1.5, fill = "gray90"),
    strip.text = element_text(face = "bold", size = 15),
    legend.position = "bottom"
  )


# ---- Plot for Brain (60 days)----
Fmax_hat_brain   <- coef_60_brain["Fmax"]
k_hat_brain      <- coef_60_brain["k"]

# MODIFIED: Simplified formula notation
eq_label_60_brain <- sprintf(
  "F(t) = %.3f × (1 - exp(-t * %.4f))",
  Fmax_hat_brain, k_hat_brain
)


pred_grid_60_brain <- tibble(
  Time = seq(0, max(release_60_day_brain$Time), length.out = 300)
)

pred_grid_60_brain$Pred <- predict(fit_exp_60day_brain, newdata = pred_grid_60_brain)


ggplot(release_60_day_brain, aes(Time, Propotion)) +
  geom_line(
    data = pred_grid_60_brain,
    aes(Time, Pred),
    linewidth = 1.2,
    color = "blue"
  ) +
  geom_point(size = 2) +
  # MODIFIED: Position formula inside the plot area
  annotate(
    "text",
    x = max(release_60_day_brain$Time) * 0.7,  
    y = 0.1,  
    label = eq_label_60_brain,
    hjust = 0.5, vjust = 0.5,
    size = 4.5,
    fontface = "italic"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    x = "Time (h)",
    y = "Fraction released",
    title = "First order release formula: 60-day in vitro release data for brain I-ADD"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, colour = "black", face = "bold"),
    axis.title = element_text(size = 18, colour = "black", face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

