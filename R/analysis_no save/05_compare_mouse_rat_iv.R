# =============================================================================
# File: 05_compare_mouse_rat_iv.R
# Purpose: Run the combined mouse and rat IV model-evaluation workflow.
# Inputs: Mouse and rat PK datasets and PBPK models.
# Outputs: Combined MRD summaries and diagnostic plots.
# Run from the repository root after installing the packages listed in README.md.
# =============================================================================

# ----------------------------
# Clean history and Set data road
# ----------------------------
# ----------------------------
# Load required libraries
# ----------------------------

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FME)
library(minpack.lm)
library(gridExtra)
library(patchwork)
library(data.table)


# ----------------------------
# Helper to force numeric columns
# ----------------------------

# ----------------------
# Mouse Model
# ----------------------


source(file.path('R', 'models', '01_model_mouse_iv.R'))  # Use the new version of I-ADD model

# Compile model
mod <- mcode("DexPBPK", DexPBPK)
MODEL_PARAM_NAMES <- names(param(mod))

# ----------------------------
# Load and clean observed data
# ----------------------------
data1 <- read.csv("data/pk/Dex_PK_Dataset_Mouse_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
  filter(Route == "IV") %>% 
  filter(Compartments %in% c("Plasma", "Liver", "Kidney", "Brain")) 

data2 <- read.csv("data/pk/Dex_PK_Dataset2_Mouse_IV_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
  filter(Route == "IV")

data_comb <- rbind(data1, data2) %>% 
  select(-1) %>%                              # remove the first column
  mutate(Number = dense_rank( 
    interaction(F_Author, Study_ID, drop = TRUE))) %>%   # same author -> same ID, starts at 1
  select(Number, everything())  %>% 
  filter (Number %in% c(2, 6, 9)) %>% 
  filter (Concentration_Mean != 0)

data_list <- 
  rbind(data1, data2) %>% 
  select(-1) %>%                              # remove the first column
  mutate(Number = dense_rank( 
    interaction(F_Author, Study_ID, drop = TRUE))) %>%   # same author -> same ID, starts at 1
  select(Number, everything()) %>% 
  select(Number, DOI, F_Author, Year, Compound, Species, Route, Condition, Dosing, Compartments) %>%
  distinct()

names(data_comb) <- make.unique(names(data_comb))

data_brain <- data_comb %>%
  filter(Compartments == "Brain")

data_liver_plasma_kidney <- data_comb %>%
  filter(Compartments %in% c("Liver", "Plasma", "Kidney"))

data <- data_comb


names(data) <- make.unique(names(data))


# Seperate the sheet

study_list <- data %>%
  group_by(Number, Study_ID) %>%
  group_split() %>%
  set_names(
    data %>%
      distinct(Number, Study_ID) %>%
      arrange(Number, Study_ID) %>%
      transmute(name = paste0("[", Number, "] Study ", Study_ID)) %>%
      pull(name)
  ) %>%
  map(~ .x %>%
        select(Compartments, Time, Concentration_Mean, Dosing) %>%
        filter(Compartments %in% c("Plasma", "Liver", "Kidney", "Brain")) %>%
        mutate(
          Compartments = recode(
            Compartments,
            Plasma = "CPLASMA_out",
            Liver  = "CLIVER_out",
            Kidney = "CKID_out",
            Brain  = "CBR_out"
          )
        )
  )


# ----------------------------
# Helper to force numeric columns
# ----------------------------

force_numeric_df <- function(df) {
  df <- as.data.frame(df)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) df[[nm]] <- unlist(df[[nm]])
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}


# ----------------------------
# Prediction function
# ----------------------------

pred <- function(pars = NULL,
                 BW = 0.02,
                 dose_per_kg,
                 tinterval_h = 24,
                 ndoses = 1,
                 t_end_h = 24,
                 dt_h = 0.01,
                 atol = 1e-10,
                 rtol = 5e-10) {
  
  # log -> natural; keep only model params by name
  pars_nat <- if (is.null(pars)) numeric(0) else {
    stopifnot(!is.null(names(pars)))  # must be named
    exp(pars)
  }
  pars_nat <- pars_nat[names(pars_nat) %in% names(param(mod))]
  
  # dosing & time grid
  ex    <- ev(amt = dose_per_kg * BW, ii = tinterval_h, addl = ndoses - 1, cmt = "AV")
  tsamp <- tgrid(0, t_end_h, dt_h)
  
  # apply params to a copy, and set BW
  simmod <- mod
  if (length(pars_nat)) simmod <- param(simmod, pars_nat)
  simmod <- param(simmod, BW = BW)
  
  out <- simmod %>%
    update(atol = atol, rtol = rtol) %>%
    mrgsim_d(data = ex, tgrid = tsamp) %>%
    as.data.frame() %>%
    dplyr::rename(Time = time) %>%
    dplyr::distinct(Time, .keep_all = TRUE)
  
  out
}

# ----------------------------
# Prediction for each study
# ----------------------------
theta.int <- log(c(PL = 2.32, 
                   PK = 0.498,
                   PR = 0.133,
                   PBR = 2.22,
                   PABRC_VAS =  0.11,
                   PABRC_INT = 0.375,
                   fu = 0.175,
                   BP = 0.725,
                   CLint = 0.020
))


sim_list <- imap(study_list, function(df, nm) {
  dose_val <- unique(df$Dosing)
  dose_val <- dose_val[!is.na(dose_val)]
  
  if (length(dose_val) == 0) return(NULL)
  
  num_val   <- if ("Number"   %in% names(df)) unique(df$Number)[1] else NA
  study_val <- if ("Study_ID" %in% names(df)) unique(df$Study_ID)[1] else NA
  
  pred(pars = theta.int, dose_per_kg = dose_val[1]) %>%
    mutate(
      Study_Name = nm,
      Dose_mgkg = dose_val[1]
    )
})



pars_log_fit <- theta.int

# -----------------------------------------------------
# 7.2 Match predictions to observations (fitted model)
# -----------------------------------------------------

pred_list <- imap(study_list, function(df, nm) {
  dose_val <- unique(df$Dosing)
  dose_val <- dose_val[!is.na(dose_val)]
  
  if (length(dose_val) == 0) return(NULL)
  
  num_val   <- if ("Number"   %in% names(df)) unique(df$Number)[1] else NA
  study_val <- if ("Study_ID" %in% names(df)) unique(df$Study_ID)[1] else NA
  
  pred(pars = pars_log_fit, dose_per_kg = dose_val[1]) %>%
    mutate(
      Study_Name = nm,
      Dose_mgkg = dose_val[1]
    )
})



combine_obs_pred <- function(study_list, pred_list,
                             model_vars = c("CPLASMA_out","CLIVER_out","CKID_out","CBR_out")) {
  
  common_studies <- intersect(names(study_list), names(pred_list))
  if (length(common_studies) == 0) {
    stop("No matching study names between study_list and pred_list")
  }
  
  purrr::map_dfr(common_studies, function(nm) {
    obs_df  <- study_list[[nm]]
    pred_df <- pred_list[[nm]]
    
    if (is.null(obs_df) || is.null(pred_df) || nrow(obs_df) == 0 || nrow(pred_df) == 0) {
      return(NULL)
    }
    
    # observed
    obs_use <- obs_df %>%
      dplyr::select(Compartments, Time, Concentration_Mean, dplyr::any_of("Dosing")) %>%
      dplyr::rename(Observed = Concentration_Mean) %>%
      dplyr::mutate(
        Study_Name = nm,
        Compartments = as.character(Compartments),
        Time = as.numeric(Time),
        Observed = as.numeric(Observed)
      ) %>%
      dplyr::filter(!is.na(Compartments), !is.na(Time), !is.na(Observed))
    
    if (nrow(obs_use) == 0) return(NULL)
    
    # prediction wide -> long (your new pred_list is wide, this is expected)
    need_cols <- c("Time", model_vars)
    miss_cols <- setdiff(need_cols, names(pred_df))
    if (length(miss_cols) > 0) {
      warning("Skipping ", nm, " missing pred columns: ", paste(miss_cols, collapse = ", "))
      return(NULL)
    }
    
    pred_long <- pred_df %>%
      dplyr::select(dplyr::all_of(need_cols)) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(model_vars),
        names_to = "Compartments",
        values_to = "Predicted"
      ) %>%
      dplyr::mutate(
        Compartments = as.character(Compartments),
        Predicted = as.numeric(Predicted)
      ) %>%
      dplyr::filter(!is.na(Time), !is.na(Compartments), !is.na(Predicted))
    
    # nearest-time match within each compartment
    comp_names <- intersect(unique(obs_use$Compartments), unique(pred_long$Compartments))
    if (length(comp_names) == 0) return(NULL)
    
    out_list <- lapply(comp_names, function(comp) {
      obs_sub  <- obs_use  %>% dplyr::filter(Compartments == comp)
      pred_sub <- pred_long %>% dplyr::filter(Compartments == comp)
      
      if (nrow(obs_sub) == 0 || nrow(pred_sub) == 0) return(NULL)
      
      idx <- vapply(obs_sub$Time, function(tt) which.min(abs(pred_sub$Time - tt)), integer(1))
      
      obs_sub %>%
        dplyr::mutate(
          Pred_Time = pred_sub$Time[idx],
          Predicted = pred_sub$Predicted[idx]
        )
    })
    
    joined <- dplyr::bind_rows(out_list)
    if (nrow(joined) == 0) return(NULL)
    
    joined %>%
      dplyr::mutate(
        Residual = Observed - Predicted,
        Log.OBS  = ifelse(Observed  > 0, log10(Observed),  NA_real_),
        Log.PRE  = ifelse(Predicted > 0, log10(Predicted), NA_real_),
        OPR      = Predicted / Observed,
        logOPR   = ifelse(OPR > 0, log10(OPR), NA_real_)
      )
  })
}


combined_mouse <- combine_obs_pred(study_list, pred_list) %>% 
  mutate(Species = "Mouse")

cat("Total matched observations:", nrow(combined_mouse), "\n")
# ------------------
# Rat Model
# ------------------
source(file.path('R', 'models', '02_model_rat_iv.R'))  # Use the new version of I-ADD model

# Compile model
mod <- mcode("DexPBPK", DexPBPK)
MODEL_PARAM_NAMES <- names(param(mod))


data3 <- read.csv("data/pk/Dex_PK_Dataset_Rat_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
  filter(Route == "IV") %>% 
  filter(Compartments %in% c("Plasma", "Liver", "Kidney", "Brain")) 

data4 <- read.csv("data/pk/Dex_PK_Dataset2_Rat_IV_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
  filter(Route == "IV") %>% 
  rename("Repeat.Dose" = "Repeat Dose")



data <- rbind(data3, data4) %>% 
  select(-1) %>%                              # remove the first column
  mutate(Number = dense_rank( 
    interaction(F_Author, Study_ID, drop = TRUE))) %>%   # same author -> same ID, starts at 1
  select(Number, everything()) %>% 
  filter(Concentration_Mean != 0) %>% 
  filter(!Number %in% c(2, 3, 5, 6 )) %>% 
  filter(Time != 0)

names(data) <- make.unique(names(data))



study_list <- data %>%
  group_by(Number, Study_ID) %>%
  group_split() %>%
  set_names(
    data %>%
      distinct(Number, Study_ID) %>%
      arrange(Number, Study_ID) %>%
      transmute(name = paste0("[", Number, "] Study ", Study_ID)) %>%
      pull(name)
  ) %>%
  map(~ .x %>%
        select(Compartments, Time, Concentration_Mean, Dosing) %>%
        filter(Compartments %in% c("Plasma", "Liver", "Kidney", "Brain")) %>%
        mutate(
          Compartments = recode(
            Compartments,
            Plasma = "CPLASMA_out",
            Liver  = "CLIVER_out",
            Kidney = "CKID_out",
            Brain  = "CBR_out"
          )
        )
  )


# ----------------------------
# Helper to force numeric columns
# ----------------------------

force_numeric_df <- function(df) {
  df <- as.data.frame(df)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) df[[nm]] <- unlist(df[[nm]])
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
  }
  df
}


# ----------------------------
# Prediction function
# ----------------------------

pred <- function(pars = NULL,
                 BW = 0.02,
                 dose_per_kg,
                 tinterval_h = 24,
                 ndoses = 1,
                 t_end_h = 24,
                 dt_h = 0.01,
                 atol = 1e-10,
                 rtol = 5e-10) {
  
  # log -> natural; keep only model params by name
  pars_nat <- if (is.null(pars)) numeric(0) else {
    stopifnot(!is.null(names(pars)))  # must be named
    exp(pars)
  }
  pars_nat <- pars_nat[names(pars_nat) %in% names(param(mod))]
  
  # dosing & time grid
  ex    <- ev(amt = dose_per_kg * BW, ii = tinterval_h, addl = ndoses - 1, cmt = "AV")
  tsamp <- tgrid(0, t_end_h, dt_h)
  
  # apply params to a copy, and set BW
  simmod <- mod
  if (length(pars_nat)) simmod <- param(simmod, pars_nat)
  simmod <- param(simmod, BW = BW)
  
  out <- simmod %>%
    update(atol = atol, rtol = rtol) %>%
    mrgsim_d(data = ex, tgrid = tsamp) %>%
    as.data.frame() %>%
    dplyr::rename(Time = time) %>%
    dplyr::distinct(Time, .keep_all = TRUE)
  
  out
}

# ----------------------------
# Prediction for each study
# ----------------------------
theta.int <- log(c( PL = 6.76, 
                    PK = 1.51,
                    PR = 0.13,
                    PABRC_VAS = 0.11,
                    PABRC_INT = 0.375,
                    fu = 0.175,
                    BP = 0.725,
                    CLint = 0.041
))


sim_list <- imap(study_list, function(df, nm) {
  dose_val <- unique(df$Dosing)
  dose_val <- dose_val[!is.na(dose_val)]
  
  if (length(dose_val) == 0) return(NULL)
  
  num_val   <- if ("Number"   %in% names(df)) unique(df$Number)[1] else NA
  study_val <- if ("Study_ID" %in% names(df)) unique(df$Study_ID)[1] else NA
  
  pred(pars = theta.int, dose_per_kg = dose_val[1]) %>%
    mutate(
      Study_Name = nm,
      Dose_mgkg = dose_val[1]
    )
})



pars_log_fit <- theta.int

# -----------------------------------------------------
# 7.2 Match predictions to observations (fitted model)
# -----------------------------------------------------

pred_list <- imap(study_list, function(df, nm) {
  dose_val <- unique(df$Dosing)
  dose_val <- dose_val[!is.na(dose_val)]
  
  if (length(dose_val) == 0) return(NULL)
  
  num_val   <- if ("Number"   %in% names(df)) unique(df$Number)[1] else NA
  study_val <- if ("Study_ID" %in% names(df)) unique(df$Study_ID)[1] else NA
  
  pred(pars = pars_log_fit, dose_per_kg = dose_val[1]) %>%
    mutate(
      Study_Name = nm,
      Dose_mgkg = dose_val[1]
    )
})



combine_obs_pred <- function(study_list, pred_list,
                             model_vars = c("CPLASMA_out","CLIVER_out","CKID_out","CBR_out")) {
  
  common_studies <- intersect(names(study_list), names(pred_list))
  if (length(common_studies) == 0) {
    stop("No matching study names between study_list and pred_list")
  }
  
  purrr::map_dfr(common_studies, function(nm) {
    obs_df  <- study_list[[nm]]
    pred_df <- pred_list[[nm]]
    
    if (is.null(obs_df) || is.null(pred_df) || nrow(obs_df) == 0 || nrow(pred_df) == 0) {
      return(NULL)
    }
    
    # observed
    obs_use <- obs_df %>%
      dplyr::select(Compartments, Time, Concentration_Mean, dplyr::any_of("Dosing")) %>%
      dplyr::rename(Observed = Concentration_Mean) %>%
      dplyr::mutate(
        Study_Name = nm,
        Compartments = as.character(Compartments),
        Time = as.numeric(Time),
        Observed = as.numeric(Observed)
      ) %>%
      dplyr::filter(!is.na(Compartments), !is.na(Time), !is.na(Observed))
    
    if (nrow(obs_use) == 0) return(NULL)
    
    # prediction wide -> long (your new pred_list is wide, this is expected)
    need_cols <- c("Time", model_vars)
    miss_cols <- setdiff(need_cols, names(pred_df))
    if (length(miss_cols) > 0) {
      warning("Skipping ", nm, " missing pred columns: ", paste(miss_cols, collapse = ", "))
      return(NULL)
    }
    
    pred_long <- pred_df %>%
      dplyr::select(dplyr::all_of(need_cols)) %>%
      dplyr::mutate(Time = as.numeric(Time)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(model_vars),
        names_to = "Compartments",
        values_to = "Predicted"
      ) %>%
      dplyr::mutate(
        Compartments = as.character(Compartments),
        Predicted = as.numeric(Predicted)
      ) %>%
      dplyr::filter(!is.na(Time), !is.na(Compartments), !is.na(Predicted))
    
    # nearest-time match within each compartment
    comp_names <- intersect(unique(obs_use$Compartments), unique(pred_long$Compartments))
    if (length(comp_names) == 0) return(NULL)
    
    out_list <- lapply(comp_names, function(comp) {
      obs_sub  <- obs_use  %>% dplyr::filter(Compartments == comp)
      pred_sub <- pred_long %>% dplyr::filter(Compartments == comp)
      
      if (nrow(obs_sub) == 0 || nrow(pred_sub) == 0) return(NULL)
      
      idx <- vapply(obs_sub$Time, function(tt) which.min(abs(pred_sub$Time - tt)), integer(1))
      
      obs_sub %>%
        dplyr::mutate(
          Pred_Time = pred_sub$Time[idx],
          Predicted = pred_sub$Predicted[idx]
        )
    })
    
    joined <- dplyr::bind_rows(out_list)
    if (nrow(joined) == 0) return(NULL)
    
    joined %>%
      dplyr::mutate(
        Residual = Observed - Predicted,
        Log.OBS  = ifelse(Observed  > 0, log10(Observed),  NA_real_),
        Log.PRE  = ifelse(Predicted > 0, log10(Predicted), NA_real_),
        OPR      = Predicted / Observed,
        logOPR   = ifelse(OPR > 0, log10(OPR), NA_real_)
      )
  })
}


combined_rat <- combine_obs_pred(study_list, pred_list) %>% 
  mutate(Species = "Rat")

cat("Total matched observations:", nrow(combined_rat), "\n")
combined <- rbind(combined_mouse, combined_rat) %>% 
  filter(Observed != 0) %>% 
  mutate(
    LogOBS = log10(Observed),
    LogPRE = log10(Predicted),
    Compartments = recode(
      Compartments,
      "CPLASMA_out" = "Plasma",
      "CLIVER_out"  = "Liver",
      "CBR_out" = "Brain",
      "CKID_out" = "Kidney"
    )
  )


# ---------------------------------------------------
# PLOT 1: Overall Goodness of Fit (fitted model)
# ---------------------------------------------------

# Calculate overall R²
# Natural-scale R²
R2_overall_nat <- cor(combined$Observed, combined$Predicted, use = "complete.obs")^2

# Log-scale R² (drop non-positive values first)
combined_log <- combined %>%
  filter(Observed > 0, Predicted > 0) %>%
  mutate(
    logObs  = log10(Observed),
    logPred = log10(Predicted)
  )

R2_overall_log <- cor(combined_log$logObs, combined_log$logPred, use = "complete.obs")^2

# By organ, both scales
R2_by_organ <- combined_log %>%
  group_by(Compartments) %>%
  summarise(
    R2_natural = cor(Observed, Predicted, use = "complete.obs")^2,
    R2_log     = cor(logObs, logPred, use = "complete.obs")^2,
    N    = n(),
    RMSE_natural = sqrt(mean((Observed - Predicted)^2, na.rm = TRUE)),
    RMSE_log     = sqrt(mean((logObs - logPred)^2, na.rm = TRUE)),
    .groups = "drop"
  )

print(R2_overall_nat)
print(R2_overall_log)
print(R2_by_organ)


# Calculate overall R²
R2_overall <- R2_overall_log

# R² by organ
R2_by_organ <- combined %>%
  group_by(Compartments) %>%
  summarise(
    R2 = cor(Observed, Predicted, use = "complete.obs")^2,
    N = n(),
    RMSE = sqrt(mean(Residual^2)),
    .groups = "drop"
  )

cat("\n=== Goodness of Fit Metrics ===\n")
cat("Overall R² =", round(R2_overall, 3), "\n\n")
cat("R² by organ:\n")


print(R2_by_organ)

# Plot
windowsFonts("Times" = windowsFont("Times New Roman"))


p_gof <- ggplot(combined, aes(x = Log.OBS, y = Log.PRE, color = Compartments, shape = Species)) +
  geom_point(size = 3.5, stroke = 1.1, alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_abline(intercept = log(2, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.5, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(3, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.33, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  scale_shape_manual(
    values = c("Mouse" = 16, "Rat"   = 17)) +
  scale_color_manual(
    values = c(
      "Brain"  = "#1f78b4",
      "Kidney" = "#33a02c",
      "Liver"  = "#e31a1c",
      "Plasma" = "grey40"
    )
  ) +
  scale_x_continuous(
    limits = c(-1, 2),
    labels = scales::math_format(10^.x)
  ) +
  scale_y_continuous(
    limits = c(-1, 2),
    labels = scales::math_format(10^.x)
  ) +
  labs(
    title = sprintf("Goodness of Fit Across All Studies (R² = %.3f)", R2_overall),
    x = "Observed ug/mL",
    y = "Predicted ug/mL",
    color = "Organ"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # text = element_text(family = "Times"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    plot.background         = element_rect (fill ="White"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.box = "vertical"
  )



print(p_gof)



# ------------------------------------------------
# Two fold plots (fitting model)
# ------------------------------------------------
# Calculate overall %2e and %3e
E2_overall <- mean(combined$OPR >= 0.5 & combined$OPR <= 2.0, na.rm = TRUE) * 100

# E3: proportion of rows where 1/3 ≤ OPB ≤ 3.0
E3_overall <- mean(combined$OPR >= 1/3 & combined$OPR <= 3.0, na.rm = TRUE) * 100

# Print neatly
cat(sprintf("E2_overall = %.2f%%\n", E2_overall))
cat(sprintf("E3_overall = %.2f%%\n", E3_overall))



# R² by organ
R2_by_organ <- combined %>%
  group_by(Compartments) %>%
  summarise(
    R2 = cor(Observed, Predicted, use = "complete.obs")^2,
    N = n(),
    RMSE = sqrt(mean(Residual^2)),
    .groups = "drop"
  )

p <- ggplot(combined, aes(Log.OBS, Log.PRE, color = Species, shape = Compartments)) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(3, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.33, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(2, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.5, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.33, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_point(
    position = position_jitter(width = 0.05, height = 0.05),
    size = 3.2, alpha = 0.9, stroke = 1.0
  ) +
  scale_color_manual(
    values = c(
      "Mouse"  = "blue",
      "Rat"  = "red"
    )
  ) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2, 2), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-2, 2), labels = scales::math_format(10^.x))

## Set up your theme and font
windowsFonts("Times" = windowsFont("Times New Roman"))

p1 <- p + 
  theme (
    plot.background         = element_rect (fill="White"),
    # text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 20, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 20, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "Observed values (μg/mL)",  y = "Predicted values (μg/mL)") +
  annotate("text", x = Inf, y = -1.8, label = paste0("R² = ", round(R2_overall, 3)), size = 5, fontface = "bold", hjust = 1.1, vjust = 1.1)

p1

p2 <-
  ggplot(combined, aes(Log.PRE, logOPR, color = Species, shape = Compartments)) +
  geom_hline(yintercept = log10(2),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(3),linetype = 3, color   = "grey", size =1) +
  geom_hline(yintercept = log10(0.33),linetype = 3, color   = "grey", size =1) +
  geom_point(
    position = position_jitter(width = 0.05, height = 0.05),
    size = 3.2, alpha = 0.9, stroke = 1.0
  ) +
   scale_color_manual(
    values = c(
      "Mouse"  = "blue",
      "Rat"  = "red"
    )
  ) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x))

p2 <- p2 +
  theme (
    plot.background         = element_rect (fill ="White"),
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill ="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 20, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 20, colour = "black", face = "bold"),   # label of axes
    legend.position         = 'bottom',
    legend.box              = "vertical",
    legend.text             = element_text(size = 15),   # label text
    legend.title            = element_text(size = 18, face = "bold"),   # if title exists
    legend.key.size  = unit(1.2, "lines")) +                        # box size) +
  labs (x = "Predicted values (μg/mL)" ,y = "Predicted/Observed") +
  annotate("text", x = Inf, y = -1.8, label = paste0("%2e: ", round(E2_overall, 2)), size = 5, fontface = "bold", hjust = 1.1, vjust = 1.1) +
  annotate("text", x = Inf, y = -1.4, label = paste0("%3e: ", round(E3_overall, 2)), size = 5, fontface = "bold", hjust = 1.1, vjust = 1.1)


p2
combined_plot <- p1 + p2 
combined_plot

# ----------------------------
# MRD Calculation
# MRD = 10^x, where x = sqrt(1/N * sum((log10(Cpred) - log10(Cobs))^2))
# ----------------------------
calc_MRD <- function(combined_df, label = "") {
  
  df_valid <- combined_df %>%
    filter(Observed > 0, Predicted > 0,
           !is.na(Log.OBS), !is.na(Log.PRE))
  
  # Overall MRD
  N_overall <- nrow(df_valid)
  x_overall <- sqrt(sum((df_valid$Log.PRE - df_valid$Log.OBS)^2) / N_overall)
  MRD_overall <- 10^x_overall
  
  cat("\n============================================================\n")
  cat("MRD Results:", label, "\n")
  cat("============================================================\n")
  cat(sprintf("Overall  | N = %d | x = %.4f | MRD = %.4f\n",
              N_overall, x_overall, MRD_overall))
  
  # MRD by organ/compartment
  MRD_by_organ <- df_valid %>%
    group_by(Compartments) %>%
    summarise(
      N   = n(),
      x   = sqrt(sum((Log.PRE - Log.OBS)^2) / n()),
      MRD = 10^x,
      .groups = "drop"
    ) %>%
    arrange(Compartments)
  
  cat("\nBy compartment:\n")
  print(MRD_by_organ, n = Inf)
  
  # MRD by study
  MRD_by_study <- df_valid %>%
    group_by(Study_Name, Compartments) %>%
    summarise(
      N   = n(),
      x   = sqrt(sum((Log.PRE - Log.OBS)^2) / n()),
      MRD = 10^x,
      .groups = "drop"
    ) %>%
    arrange(Study_Name, Compartments)
  
  cat("\nBy study and compartment:\n")
  print(MRD_by_study, n = Inf)
  
  list(
    overall      = data.frame(Label = label, N = N_overall, x = x_overall, MRD = MRD_overall),
    by_organ     = MRD_by_organ %>% mutate(Label = label),
    by_study     = MRD_by_study %>% mutate(Label = label)
  )
}

# Run for mouse and rat
MRD_mouse <- calc_MRD(combined_mouse, label = "Mouse")
MRD_rat   <- calc_MRD(combined_rat,   label = "Rat")
MRD_all   <- calc_MRD(combined_log)


# Combined summary table
MRD_summary <- bind_rows(
  MRD_mouse$overall,
  MRD_rat$overall,
  MRD_all$overall
)

print(MRD_summary)

# Save results
MRD_organ_all <- bind_rows(MRD_mouse$by_organ, MRD_rat$by_organ)
MRD_study_all <- bind_rows(MRD_mouse$by_study, MRD_rat$by_study)

print(MRD_organ_all)
print(MRD_study_all)


write.csv(MRD_organ_all, file.path("results", "tables", "MRD_by_organ.csv"),  row.names = FALSE)
write.csv(MRD_study_all, file.path("results", "tables", "MRD_by_study.csv"),  row.names = FALSE)
write.csv(MRD_summary, file.path("results", "tables", "MRD_overall.csv"),   row.names = FALSE)
