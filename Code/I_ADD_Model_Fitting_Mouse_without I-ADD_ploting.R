############################################################
# Clean and robust I-ADD PBPK fitting workflow
# Xiaowen Wang, 2026 April mouse
############################################################

# ----------------------------
# Clean history and Set data road
# ----------------------------
rm(list = ls())
curr.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)


# ----------------------------
# Load required libraries
# ----------------------------

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(FME)
library(minpack.lm)
library(gridExtra)
library(patchwork)
library(data.table)



source (file = 'I_ADD_PBPK_Model_2026May22_Brain Three Layer Kideny One Layer_mouse_without I-ADD_with fuorgan.R')  # Use the new version of I-ADD model

# Compile model
mod <- mcode("DexPBPK", DexPBPK)
MODEL_PARAM_NAMES <- names(param(mod))

# ----------------------------
# Load and clean observed data
# ----------------------------
data1 <- read.csv("../../Database/PK_Data/Dex/Dex_PK_Dataset_Mouse_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
  filter(Route == "IV") %>% 
  filter(Compartments %in% c("Plasma", "Liver", "Kidney", "Brain")) 

data2 <- read.csv("../../Database/PK_Data/03_Data searching for DEX/Mouse from Pia/Dex_PK_Dataset2_Mouse_IV_Unit_Consist_Round4.csv", check.names = FALSE, strip.white = TRUE) %>% 
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

# make.unique can rename duplicate columns (e.g. F_Author -> F_Author.1);
# restore F_Author to its original name so downstream selects find it
names(data_comb) <- sub("^F_Author\\..*$", "F_Author", names(data_comb))

# ----------------------------
# Build human-readable study labels: "First Author et al. (Year), X mg/kg"
# Used instead of "[N] Study ID" in all facet strip labels
# ----------------------------
study_label_map <- rbind(data1, data2) %>%
  select(-1) %>%
  mutate(Number = dense_rank(interaction(F_Author, Study_ID, drop = TRUE))) %>%
  filter(Number %in% c(2, 6, 9)) %>%
  group_by(Number, Study_ID, F_Author, Year) %>%
  summarise(Dosing = first(na.omit(Dosing)), .groups = "drop") %>%
  mutate(
    # Format author: keep only last name + "et al." if there is more than one word
    Author_short = sub("^(\\S+).*", "\\1", F_Author),
    Study_Label  = paste0(Author_short, " et al. (", Year, "), ", Dosing, " mg/kg"),
    Name_key     = paste0("[", Number, "] Study ", Study_ID)   # old key to match on
  )

# Named vector: old key -> new label
label_lookup <- setNames(study_label_map$Study_Label, study_label_map$Name_key)


data_brain <- data_comb %>%
  filter(Compartments == "Brain")

data_liver_plasma_kidney <- data_comb %>%
  filter(Compartments %in% c("Liver", "Plasma", "Kidney"))

data <- data_liver_plasma_kidney

# ---------------------------------
# Study_list
# ---------------------------------
study_list <- data %>%
  group_by(Number, Study_ID) %>%
  group_split() %>%
  set_names(
    data %>%
      distinct(Number, Study_ID) %>%
      arrange(Number, Study_ID) %>%
      transmute(name = label_lookup[paste0("[", Number, "] Study ", Study_ID)]) %>%
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


theta.int <- log(c(PL          = 6.76,
                   PK          = 1.51,
                   PR          = 1.0,
                   PBR         = 1.2,
                   PABRC_VAS   = 0.05,# initial 0.6
                   PABRC_INT   = 0.8, # initial 0.6
                   fu          = 0.175,
                   fucel_brain = 0.24,
                   fuint_brain = 0.29,   
                   BP          = 0.725,
                   CLint       = 0.4     # initial 0.4
))


# ---------------------------------
# Function: sim_list
# ---------------------------------

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

# ----------------------------
# Residual function for sensFun / modFit (list-based, automatic doses)
# ----------------------------

total_residuals <- function(pars_log) {
  model_vars <- c("CPLASMA_out","CLIVER_out","CKID_out","CBR_out")
  
  # IMPORTANT:
  # If modFit passes only a subset (theta_sel), merge into full theta.int baseline
  pars_full <- theta.int
  pars_full[names(pars_log)] <- pars_log
  
  res_all <- list()
  k <- 1
  
  for (nm in names(study_list)) {
    obs_df <- study_list[[nm]]
    if (is.null(obs_df) || nrow(obs_df) == 0) next
    
    # get dose automatically from study data
    dose_vals <- unique(obs_df$Dosing)
    dose_vals <- dose_vals[!is.na(dose_vals)]
    if (length(dose_vals) == 0) next
    dose_val <- dose_vals[1]
    
    # simulate to observed max time (or use 24 h fallback)
    t_end <- suppressWarnings(max(obs_df$Time, na.rm = TRUE))
    if (!is.finite(t_end)) t_end <- 24
    
    pred_df <- pred(pars = pars_full, dose_per_kg = dose_val, t_end_h = t_end)
    
    # observation table for FME (must be base data.frame)
    obs_use <- obs_df %>%
      dplyr::filter(!is.na(Time), !is.na(Compartments), !is.na(Concentration_Mean)) %>%
      dplyr::select(Compartments, Time, Concentration_Mean) %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    if (nrow(obs_use) == 0) next
    
    # prediction table for FME
    need_cols <- c("Time", model_vars)
    miss_cols <- setdiff(need_cols, names(pred_df))
    if (length(miss_cols) > 0) next
    
    pred_use <- pred_df[, need_cols, drop = FALSE] %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    # type safety (important for FME::modCost)
    obs_use$Compartments <- as.character(obs_use$Compartments)
    obs_use$Time <- as.numeric(obs_use$Time)
    obs_use$Concentration_Mean <- as.numeric(obs_use$Concentration_Mean)
    
    pred_use$Time <- as.numeric(pred_use$Time)
    for (v in model_vars) pred_use[[v]] <- as.numeric(pred_use[[v]])
    
    # modCost for this study
    cst <- FME::modCost(
      model = pred_use,
      obs   = obs_use,
      x     = "Time",
      y     = "Concentration_Mean",
      var   = "Compartments"
    )
    
    # extract NUMERIC residuals (this is what sensFun/modFit need)
    rtab <- cst$residuals
    rvec <- rtab$res
    
    # optional: informative names
    names(rvec) <- paste0(nm, "::", rtab$var, "::t=", rtab$x)
    
    res_all[[k]] <- rvec
    k <- k + 1
  }
  
  if (length(res_all) == 0) stop("No residuals generated from study_list.")
  
  unlist(res_all, use.names = TRUE)
}

# ----------------------------
# Sensitivity analysis
# ----------------------------

Sns <- FME::sensFun(func = total_residuals, parms = theta.int)
Sen <- summary(Sns)

# Rank robustly
Sen_rank <- Sen
names(Sen_rank) <- tolower(names(Sen_rank))
if (!"mean" %in% names(Sen_rank)) stop("summary(Sns) missing 'Mean' column")

Sen_rank$meanabs <- abs(Sen_rank$mean)
Sen_rank <- Sen_rank[order(-Sen_rank$meanabs), , drop = FALSE]

par_names <- rownames(Sen_rank)
top_params <- head(par_names, 12)
top_params

# Subset starting vector
theta_sel <- theta.int[c("PL", "PK", "PR", "CLint")]

cat("\nSelected parameters for fitting:\n"); print(theta_sel)

# ----------------------------
# Bounded fitting (±5× around starts)
# ----------------------------
theta_nat <- exp(theta_sel)
lower_nat <- theta_nat / 100
upper_nat <- theta_nat * 10

lower <- log(lower_nat)
upper <- log(upper_nat)

Fit <- modFit(
  f       = total_residuals,
  p       = theta_sel,       # log-scale starts
  lower   = lower,           # log-scale bounds
  upper   = upper,           # log-scale bounds
  method  = "Port",          # honors bounds
  control = list(eval.max = 300, iter.max = 300, trace = 1)
)

print(summary(Fit))
cat("\nFitted parameters (natural scale):\n"); print(exp(Fit$par))

pars_log_fit <- theta.int
pars_log_fit[names(Fit$par)] <- Fit$par 

# -------------------------------
# Step wise model fitting: Step 2
# -------------------------------

data <- data_brain
theta.int <- pars_log_fit
theta_sel <- theta.int[c("PBR", "PABRC_INT", "PABRC_VAS")]

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

# ----------------------------
# Bounded fitting (±5× around starts)
# ----------------------------
theta_nat <- exp(theta_sel)
lower_nat <- theta_nat / 10
upper_nat <- theta_nat * 10

lower <- log(lower_nat)
upper <- log(upper_nat)

Fit <- modFit(
  f       = total_residuals,
  p       = theta_sel,       # log-scale starts
  lower   = lower,           # log-scale bounds
  upper   = upper,           # log-scale bounds
  method  = "Port",          # honors bounds
  control = list(eval.max = 300, iter.max = 300, trace = 1)
)

print(summary(Fit))
cat("\nFitted parameters (natural scale):\n"); print(exp(Fit$par))

pars_log_fit <- theta.int
pars_log_fit[names(Fit$par)] <- Fit$par 


# -----------------------------------------------------
# 7.2 Match predictions to observations (fitted model)
# -----------------------------------------------------
data <- data_comb

study_list <- data %>%
  group_by(Number, Study_ID) %>%
  group_split() %>%
  set_names(
    data %>%
      distinct(Number, Study_ID) %>%
      arrange(Number, Study_ID) %>%
      transmute(name = label_lookup[paste0("[", Number, "] Study ", Study_ID)]) %>%
      pull(name)
  ) %>%
  map(~ .x %>%
        select(Compartments, Time, Concentration_Mean, Dosing, Year, any_of("F_Author")) %>%
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
      dplyr::select(Compartments, Time, Concentration_Mean, dplyr::any_of(c("Dosing", "F_Author"))) %>%
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


combined <- combine_obs_pred(study_list, pred_list)

cat("Total matched observations:", nrow(combined), "\n")
head(combined)


# ---------------------------------------------------
# 7.3 PLOT 1: Overall Goodness of Fit (fitted model)
# ---------------------------------------------------

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
p_gof <- ggplot(combined, aes(x = Observed, y = Predicted, color = Compartments)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", 
              linewidth = 1) +
  scale_x_log10(
    limits = c(0.01, 100),
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")) +
  scale_y_log10(
    limits = c(0.01, 100),
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(
    title = sprintf("Goodness of Fit Across All Studies (R² = %.3f)", R2_overall),
    x = "Observed ug/mL",
    y = "Predicted ug/mL",
    color = "Organ"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_gof)

# -------------------------------------------------------------------
# 7.4 PLOT 2: Obs vs Pred by Study (faceted by organ) (initial model)
# -------------------------------------------------------------------
# Unique study IDs present in BOTH tables

combine_obs_pred_list <- function(study_list, pred_list,
                                  model_vars = c("CPLASMA_out","CLIVER_out","CKID_out","CBR_out")) {
  common_studies <- intersect(names(study_list), names(pred_list))
  
  out <- setNames(vector("list", length(common_studies)), common_studies)
  
  for (nm in common_studies) {
    obs_df  <- study_list[[nm]]
    pred_df <- pred_list[[nm]]
    
    if (is.null(obs_df) || is.null(pred_df) || nrow(obs_df) == 0 || nrow(pred_df) == 0) {
      out[[nm]] <- NULL
      next
    }
    
    # observed
    obs_use <- obs_df %>%
      dplyr::select(Compartments, Time, Concentration_Mean, dplyr::any_of(c("Dosing", "F_Author", "Year"))) %>%
      dplyr::rename(Observed = Concentration_Mean) %>%
      dplyr::mutate(
        Study_Name = nm,
        Compartments = as.character(Compartments),
        Time = as.numeric(Time),
        Observed = as.numeric(Observed)
      ) %>%
      filter(!is.na(Compartments), !is.na(Time), !is.na(Observed))
    
    # prediction (wide -> long)
    need_cols <- c("Time", model_vars)
    miss_cols <- setdiff(need_cols, names(pred_df))
    if (length(miss_cols) > 0) {
      out[[nm]] <- NULL
      next
    }
    
    pred_long <- pred_df %>%
      select(all_of(need_cols)) %>%
      mutate(Time = as.numeric(Time)) %>%
      pivot_longer(
        cols = all_of(model_vars),
        names_to = "Compartments",
        values_to = "Predicted"
      ) %>%
      mutate(
        Compartments = as.character(Compartments),
        Predicted = as.numeric(Predicted),
        Study_Name = nm
      )
    
    # nearest-time matched observed + predicted points (for R² / residuals)
    comp_names <- intersect(unique(obs_use$Compartments), unique(pred_long$Compartments))
    matched_list <- lapply(comp_names, function(comp) {
      obs_sub  <- obs_use  %>% filter(Compartments == comp)
      pred_sub <- pred_long %>% filter(Compartments == comp)
      if (nrow(obs_sub) == 0 || nrow(pred_sub) == 0) return(NULL)
      
      idx <- vapply(obs_sub$Time, function(tt) which.min(abs(pred_sub$Time - tt)), integer(1))
      obs_sub %>%
        mutate(
          Pred_Time = pred_sub$Time[idx],
          Predicted = pred_sub$Predicted[idx]
        )
    })
    
    matched_df <- bind_rows(matched_list)
    
    out[[nm]] <- list(
      obs = obs_use,
      pred = pred_long,
      matched = matched_df
    )
  }
  
  out
}

plot_list <- combine_obs_pred_list(study_list, pred_list)

study_names <- names(plot_list)


# Optional prettier organ labels
organ_labs <- c(
  CPLASMA_out  = "Plasma",
  CLIVER_out   = "Liver",
  CKID_out = "Kidney",
  CBR_out  = "Brain"
)

pred_all_list <- list()
obs_all_list  <- list()

for (study_nm in names(plot_list)) {
  
  item <- plot_list[[study_nm]]
  
  if (is.null(item)) {
    message("Skipping ", study_nm, " - no data available")
    next
  }
  
  sim_fit  <- item$pred
  obs_plot <- item$obs
  
  if (nrow(sim_fit) == 0 || nrow(obs_plot) == 0) {
    message("Skipping ", study_nm, " - empty obs/pred")
    next
  }
  
  f_author_val  <- if ("F_Author" %in% names(obs_plot)) unique(obs_plot$F_Author)[1] else study_nm
  year_val      <- if ("Year"     %in% names(obs_plot)) unique(na.omit(obs_plot$Year))[1]  else NA
  dose_val_plot <- if ("Dosing"   %in% names(obs_plot)) unique(na.omit(obs_plot$Dosing))[1] else NA
  row_label <- study_nm   # study_nm already = "Author et al. (Year), X mg/kg" from label_lookup
  
  sim_fit <- sim_fit %>%
    mutate(
      Study      = study_nm,
      F_Author   = f_author_val,
      Year       = year_val,
      Row_Label  = row_label,
      Organ      = recode(Compartments, !!!organ_labs, .default = Compartments)
    )
  
  obs_plot <- obs_plot %>%
    mutate(
      Study      = study_nm,
      F_Author   = f_author_val,
      Year       = year_val,
      Row_Label  = row_label,
      Organ      = recode(Compartments, !!!organ_labs, .default = Compartments)
    )
  
  pred_all_list[[study_nm]] <- sim_fit
  obs_all_list[[study_nm]]  <- obs_plot
}

# Combine all studies
pred_all <- bind_rows(pred_all_list)
obs_all  <- bind_rows(obs_all_list)

# Keep only Row_Label + Organ panels that have BOTH prediction and observation
valid_panels <- inner_join(
  pred_all %>% distinct(Row_Label, Organ),
  obs_all  %>% distinct(Row_Label, Organ),
  by = c("Row_Label", "Organ")
)

pred_all <- pred_all %>%
  semi_join(valid_panels, by = c("Row_Label", "Organ"))

obs_all <- obs_all %>%
  semi_join(valid_panels, by = c("Row_Label", "Organ"))

# Optional: control organ order
organ_order <- c( "Liver", "Kidney", "Brain","Plasma")
year_order <- c(2011, 2010, 2007)


row_order <- bind_rows(
  obs_all  %>% distinct(Row_Label, Year),
  pred_all %>% distinct(Row_Label, Year)
) %>%
  distinct(Row_Label, Year) %>%
  mutate(
    Year = as.numeric(as.character(Year)),
    Year_order = match(Year, year_order)
  ) %>%
  arrange(Year_order) %>%
  pull(Row_Label)

# 3. Apply factor order to both datasets
obs_all <- obs_all %>%
  mutate(
    Row_Label = factor(Row_Label, levels = row_order),
    Organ     = factor(Organ, levels = organ_order)
  )

pred_all <- pred_all %>%
  mutate(
    Row_Label = factor(Row_Label, levels = row_order),
    Organ     = factor(Organ, levels = organ_order)
  )


# Organ colour palette — print-safe, colour-blind friendly
organ_colors <- c(
  Plasma = "#1a1a1a",
  Liver  = "#2e7d32",
  Kidney = "#c62828",
  Brain  = "#1565c0"
)

# Drop panels where there is no observed data (e.g. Study 12 Plasma/Kidney/Brain)
valid_panels2 <- obs_all %>% distinct(Row_Label, Organ)
pred_all <- pred_all %>% semi_join(valid_panels2, by = c("Row_Label", "Organ"))
obs_all  <- obs_all  %>% semi_join(valid_panels2, by = c("Row_Label", "Organ"))

# install.packages("ggh4x")  # run once if not installed
library(ggplot2)
library(ggh4x)

p_all <- ggplot() +
  geom_line(
    data = pred_all,
    aes(x = Time, y = Predicted, colour = Organ),
    linewidth = 0.7,
    alpha = 0.9
  ) +
  geom_point(
    data = obs_all,
    aes(x = Time, y = Observed, colour = Organ),
    shape = 16,
    size  = 1.8,
    alpha = 0.9
  ) +
  ggh4x::facet_grid2(
    rows = vars(Row_Label),
    cols = vars(Organ),
    scales = "free_y",
    independent = "y",
    space = "fixed",
    drop = TRUE,
    axes = "x",              # show x axes for all panels
    remove_labels = "none"   # keep x-axis tick labels for all panels
  )+
  scale_colour_manual(values = organ_colors, name = "Organ") +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = c(0, 2, 4, 6, 8, 10),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = NULL,
    x = "Time (hours)",
    y = expression("Concentration (" * mu * "g/mL)")
  ) +
  theme_classic(base_size = 9, base_family = "sans") +
  theme(
    plot.background  = element_blank(),
    panel.background = element_blank(),
    
    # remove grid/background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # keep panel border
    panel.border = element_rect(
      colour = "grey40",
      linewidth = 0.5,
      fill = NA
    ),
    
    panel.spacing = unit(0.5, "lines"),
    
    # remove background of facet labels, including right-side row labels
    strip.background   = element_blank(),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    
    strip.text.x = element_text(
      size = 8,
      face = "bold",
      margin = margin(2, 2, 2, 2)
    ),
    strip.text.y = element_text(
      size = 7,
      face = "bold",
      angle = 0,
      hjust = 0,
      margin = margin(2, 4, 2, 4)
    ),
    
    axis.text = element_text(size = 7, colour = "grey20"),
    axis.title = element_text(size = 8.5, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
    axis.ticks.length = unit(0.15, "cm"),
    
    legend.position = "bottom",
    plot.margin = margin(6, 8, 4, 6)
  )

print(p_all)
# PNG preview
ggsave(
  filename = "Figure_DEX_PBPK_Mouse_preview.png",
  plot     = p_all,
  width    = 178,
  height   = 140,
  units    = "mm",
  dpi      = 300
)

