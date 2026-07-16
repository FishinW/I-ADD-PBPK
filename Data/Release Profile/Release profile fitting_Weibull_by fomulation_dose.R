

rm(list = ls())
curr.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)

library(dplyr)
library(tidyr)
library(ggplot2)
library(minpack.lm)


release_data <- read.csv("Data for fitting release profile.csv")

weibull_fun <- function(Time, Fmax, lambda, k) {
  Fmax * (1 - exp(-(Time / lambda)^k))
}

# --- 30-days release---

release_30_day <- release_data %>% 
  filter(Device == "30-day in vitro release")

fit_weibull_30 <- nlsLM(
  Propotion ~ weibull_fun(Time, Fmax, lambda, k),
  data = release_30_day,
  start = list(
    Fmax   = max(release_30_day$Propotion),  # ~0.99
    lambda = 100,                              # hours
    k      = 0.5                               # burst-like
  ),
  lower = c(Fmax = 0.8, lambda = 1,   k = 0.01),
  upper = c(Fmax = 1.0, lambda = 1000, k = 5),
  control = nls.lm.control(maxiter = 500)
)

summary(fit_weibull_30)
coef(fit_weibull_30)

pred_grid <- tibble(
  Time = seq(0, max(release_30_day$Time), length.out = 300)
)

pred_grid$Pred <- predict(fit_weibull_30, newdata = pred_grid)

coef_30 <- coef(fit_weibull_30)

Fmax_hat   <- coef_30["Fmax"]
lambda_hat <- coef_30["lambda"]
k_hat      <- coef_30["k"]

eq_label_30 <- sprintf(
  "F(t) = %.3f × (1 - exp(-(t / %.1f)^%.2f))",
  Fmax_hat, lambda_hat, k_hat
)

ggplot(release_30_day, aes(Time, Propotion)) +
  geom_line(data = pred_grid, aes(Time, Pred),
            linewidth = 1.2, color = "blue") +
  geom_point(size = 2) +
  annotate(
    "text",
    x = max(release_30_day$Time) * 0.7,  
    y = 0.1,
    label = eq_label_30,
    hjust = 0.5, vjust = 0.5,
    size = 4.5,
    fontface = "italic"
  ) + 
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)) +
  labs(
    x = "Time (h)",
    y = "Fraction released",
    title = "Weibull release formula: 30-day in vitro release data"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25),
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(colour = "black", linewidth = 1.5, fill = "gray90"),
    strip.text = element_text(face = "bold", size = 15),
    legend.position = "bottom"
  )

