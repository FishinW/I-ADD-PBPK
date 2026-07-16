############################################################
# Full PBPK + AUC + Normalized Sensitivity Analysis (NSC)
# Two separate mrgsolve models:
#   1) brain model
#   2) kidney model
############################################################

rm(list = ls())
curr.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)


# ============================================================
# Libraries
# ============================================================
library(mrgsolve)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

# ============================================================
# USER SETTINGS
# ============================================================

# simulation settings
sim_end   <- 168      # h
sim_delta <- 0.1      # h
BW_use    <- 0.25

# NSC settings
perturb_frac <- 0.01  # 1% perturbation
top_n_plot   <- 10

# parameters excluded from SA
# adjust this list to match your final model parameter names
exclude_pars <- c()   #"BW", "DOSE_IADD", "k", "FMAX", "KREL", "FMAX_REL"

# model files
brain_model_file  <- "I_ADD_PBPK_Model_2026Jul10_Brain Three Layer Kidney One Layer_rat_with I-ADD in brain_with fuorgan_AUC.R"
kidney_model_file <- "I_ADD_PBPK_Model_2026Jul10_Brain Three Layer Kidney One Layer_rat_with I-ADD in kidney_with fuorgan_AUC.R"

# compiled model names
brain_model_name  <- "DexPBPK_brain"
kidney_model_name <- "DexPBPK_kidney"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

compile_model_from_file <- function(model_file, model_name = "DexPBPK") {
  source(file = model_file, local = TRUE)
  
  if (!exists("DexPBPK")) {
    stop(paste0("Object 'DexPBPK' was not found after sourcing: ", model_file))
  }
  
  mod <- mcode(model_name, DexPBPK)
  return(mod)
}

run_sim <- function(mod, pars = NULL, end = sim_end, delta = sim_delta, BW = BW_use) {
  simmod <- mod
  
  if (!is.null(pars)) {
    simmod <- param(simmod, pars)
  }
  
  # simmod <- param(simmod, BW = BW)
  
  out <- simmod %>%
    mrgsim(end = end, delta = delta) %>%
    as.data.frame()
  
  out
}

get_final_auc <- function(sim_df, target = c("brain", "kidney", "plasma")) {
  target <- match.arg(target)
  last_row <- sim_df %>% slice_tail(n = 1)
  
  if (target == "brain") {
    if (!"AUC_BRAIN_out" %in% names(last_row)) {
      stop("AUC_BRAIN_out not found in simulation output.")
    }
    return(last_row$AUC_BRAIN_out)
  }
  
  if (target == "kidney") {
    if (!"AUC_KIDNEY_out" %in% names(last_row)) {
      stop("AUC_KIDNEY_out not found in simulation output.")
    }
    return(last_row$AUC_KIDNEY_out)
  }
  
  if (target == "plasma") {
    if (!"AUC_PLASMA_out" %in% names(last_row)) {
      stop("AUC_PLASMA_out not found in simulation output.")
    }
    return(last_row$AUC_PLASMA_out)
  }
}

calc_nsc <- function(mod,
                     base_pars,
                     target = c("brain", "kidney", "plasma"),
                     frac = 0.01,
                     end = sim_end,
                     delta = sim_delta,
                     BW = BW_use,
                     eps_y = 1e-12) {
  
  target <- match.arg(target)
  
  base_sim <- run_sim(mod, pars = base_pars, end = end, delta = delta, BW = BW)
  y0 <- get_final_auc(base_sim, target = target)
  
  if (is.na(y0) || abs(y0) < eps_y) {
    stop(paste0("Baseline AUC for ", target, " is too small or NA."))
  }
  
  res <- lapply(names(base_pars), function(pn) {
    p0 <- base_pars[[pn]]
    
    pars_up <- base_pars
    pars_dn <- base_pars
    
    pars_up[[pn]] <- p0 * (1 + frac)
    pars_dn[[pn]] <- p0 * (1 - frac)
    
    sim_up <- run_sim(mod, pars = pars_up, end = end, delta = delta, BW = BW)
    sim_dn <- run_sim(mod, pars = pars_dn, end = end, delta = delta, BW = BW)
    
    y_up <- get_final_auc(sim_up, target = target)
    y_dn <- get_final_auc(sim_dn, target = target)
    
    nsc <- (y_up - y_dn) / (2 * frac * y0)
    
    data.frame(
      Parameter = pn,
      BaselineParameter = p0,
      BaselineAUC = y0,
      AUC_up = y_up,
      AUC_down = y_dn,
      NSC = nsc,
      absNSC = abs(nsc),
      Direction = case_when(
        nsc > 0  ~ "Positive",
        nsc < 0  ~ "Negative",
        TRUE     ~ "Zero"
      ),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(res) %>%
    arrange(desc(absNSC))
}

plot_nsc <- function(nsc_df, title_text, top_n = 10, use_relative = FALSE) {
  
  plot_df <- nsc_df %>%
    slice_max(order_by = absNSC, n = top_n) %>%
    arrange(NSC) %>%
    mutate(
      PlotValue = if (use_relative) NSC / max(abs(NSC)) * 100 else NSC,
      Parameter = factor(Parameter, levels = Parameter),
      Direction = factor(Direction, levels = c("Negative", "Zero", "Positive"))
    )
  
  xlab_text <- if (use_relative) {
    "Relative normalized sensitivity (% of max signed NSC)"
  } else {
    "Normalized sensitivity coefficient (NSC)"
  }
  
  max_abs <- max(abs(plot_df$PlotValue), na.rm = TRUE)
  
  ggplot(plot_df, aes(x = Parameter, y = PlotValue, fill = Direction)) +
    geom_col(width = 0.75) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    scale_y_continuous(limits = c(-1.1 * max_abs, 1.1 * max_abs)) +
    scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue", "Zero" = "grey70")) +
    labs(
      title = title_text,
      x = NULL,
      y = xlab_text,
      fill = "Effect"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.text = element_text(size = 13, colour = "black", face = "bold"),
      axis.title = element_text(size = 15, colour = "black", face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 12, face = "bold")
    )
}

prepare_base_pars <- function(mod, exclude_pars = NULL) {
  all_pars <- param(mod)
  sa_par_names <- setdiff(names(all_pars), exclude_pars)
  sa_par_vals  <- unlist(all_pars[sa_par_names])
  
  if (any(sa_par_vals <= 0)) {
    bad_pars <- names(sa_par_vals)[sa_par_vals <= 0]
    stop(
      paste0(
        "Some selected parameters are <= 0 and cannot be perturbed by percent change: ",
        paste(bad_pars, collapse = ", ")
      )
    )
  }
  
  return(sa_par_vals)
}

run_model_nsc_workflow <- function(model_file,
                                   model_name,
                                   target_main = c("brain", "kidney"),
                                   include_plasma = TRUE,
                                   exclude_pars = exclude_pars,
                                   sim_end = sim_end,
                                   sim_delta = sim_delta,
                                   BW_use = BW_use,
                                   perturb_frac = perturb_frac,
                                   top_n_plot = top_n_plot) {
  
  target_main <- match.arg(target_main)
  
  mod <- compile_model_from_file(model_file = model_file, model_name = model_name)
  base_pars <- prepare_base_pars(mod, exclude_pars = exclude_pars)
  
  baseline_out <- run_sim(mod, pars = base_pars, end = sim_end, delta = sim_delta, BW = BW_use)
  
  cat("\n============================================================\n")
  cat("Model file:", model_file, "\n")
  cat("Primary target:", target_main, "\n")
  cat("============================================================\n")
  
  print(
    tail(
      baseline_out %>%
        dplyr::select(any_of(c(
          "time",
          "CPLASMA_out", "CBRAIN_out", "CKIDNEY_out",
          "AUC_PLASMA_out", "AUC_BRAIN_out", "AUC_KIDNEY_out"
        ))),
      5
    )
  )
  
  auc_cols <- intersect(
    c("AUC_PLASMA_out", "AUC_BRAIN_out", "AUC_KIDNEY_out"),
    names(baseline_out)
  )
  
  print(
    baseline_out %>%
      slice_tail(n = 1) %>%
      select(all_of(auc_cols))
  )
  
  main_nsc <- calc_nsc(
    mod = mod,
    base_pars = base_pars,
    target = target_main,
    frac = perturb_frac,
    end = sim_end,
    delta = sim_delta,
    BW = BW_use
  )
  
  main_plot <- plot_nsc(
    nsc_df = main_nsc,
    title_text = paste0(
      "Normalized sensitivity of ", target_main,
      " AUC (0-", sim_end, " h)"
    ),
    top_n = top_n_plot,
    use_relative = FALSE
  )
  
  plasma_nsc <- NULL
  plasma_plot <- NULL
  
  if (include_plasma && "AUC_PLASMA_out" %in% names(baseline_out)) {
    plasma_nsc <- calc_nsc(
      mod = mod,
      base_pars = base_pars,
      target = "plasma",
      frac = perturb_frac,
      end = sim_end,
      delta = sim_delta,
      BW = BW_use
    )
    
    plasma_plot <- plot_nsc(
      nsc_df = plasma_nsc,
      title_text = paste0(
        "Normalized sensitivity of plasma AUC (0-", sim_end, " h)"
      ),
      top_n = top_n_plot,
      use_relative = FALSE
    )
  }
  
  list(
    mod = mod,
    base_pars = base_pars,
    baseline_out = baseline_out,
    main_nsc = main_nsc,
    main_plot = main_plot,
    plasma_nsc = plasma_nsc,
    plasma_plot = plasma_plot
  )
}

# ============================================================
# RUN BRAIN MODEL
# ============================================================

brain_res <- run_model_nsc_workflow(
  model_file = brain_model_file,
  model_name = brain_model_name,
  target_main = "brain",
  include_plasma = TRUE,
  exclude_pars = exclude_pars,
  sim_end = sim_end,
  sim_delta = sim_delta,
  BW_use = BW_use,
  perturb_frac = perturb_frac,
  top_n_plot = top_n_plot
)

brain_nsc <- brain_res$main_nsc
p_brain   <- brain_res$main_plot

brain_plasma_nsc <- brain_res$plasma_nsc
p_brain_plasma   <- brain_res$plasma_plot

print(head(brain_nsc, 15))
print(p_brain)

if (!is.null(brain_plasma_nsc)) {
  print(head(brain_plasma_nsc, 15))
  print(p_brain_plasma)
}

# ============================================================
# RUN KIDNEY MODEL
# ============================================================

kidney_res <- run_model_nsc_workflow(
  model_file = kidney_model_file,
  model_name = kidney_model_name,
  target_main = "kidney",
  include_plasma = TRUE,
  exclude_pars = exclude_pars,
  sim_end = sim_end,
  sim_delta = sim_delta,
  BW_use = BW_use,
  perturb_frac = perturb_frac,
  top_n_plot = top_n_plot
)

kidney_nsc <- kidney_res$main_nsc
p_kidney   <- kidney_res$main_plot

kidney_plasma_nsc <- kidney_res$plasma_nsc
p_kidney_plasma   <- kidney_res$plasma_plot

print(head(kidney_nsc, 15))
print(p_kidney)

if (!is.null(kidney_plasma_nsc)) {
  print(head(kidney_plasma_nsc, 15))
  print(p_kidney_plasma)
}

# ============================================================
# OPTIONAL COMBINED DISPLAY
# ============================================================

plot_list <- list(p_brain, p_kidney)

if (!is.null(p_brain_plasma))  plot_list <- append(plot_list, list(p_brain_plasma))
if (!is.null(p_kidney_plasma)) plot_list <- append(plot_list, list(p_kidney_plasma))

wrap_plots(plotlist = plot_list, ncol = 2)


# ============================================================
# Circular NSC plot
# ============================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)
library(grid)

plot_circular_nsc <- function(
    nsc_list,
    cutoff = 0.30,
    top_n_per_group = NULL,
    empty_bar = 3,
    ring_breaks = c(0.3, 0.6, 0.9),
    label_size = 3.5,
    group_label_size = 5,
    title = "Normalized sensitivity coefficients"
) {
  
  # ----------------------------------------------------------
  # 1. Check input
  # ----------------------------------------------------------
  if (is.null(names(nsc_list)) || any(names(nsc_list) == "")) {
    stop("nsc_list must be a named list, for example list(Brain = brain_nsc).")
  }
  
  required_cols <- c("Parameter", "NSC")
  
  bad_input <- purrr::keep(
    nsc_list,
    ~ !all(required_cols %in% names(.x))
  )
  
  if (length(bad_input) > 0) {
    stop("Each NSC table must contain columns named 'Parameter' and 'NSC'.")
  }
  
  # Remove NULL objects
  nsc_list <- purrr::compact(nsc_list)
  
  if (length(nsc_list) == 0) {
    stop("No valid NSC tables were supplied.")
  }
  
  group_order <- names(nsc_list)
  
  # ----------------------------------------------------------
  # 2. Combine NSC results
  # ----------------------------------------------------------
  plot_data <- purrr::imap_dfr(
    nsc_list,
    function(df, group_name) {
      
      out <- df %>%
        transmute(
          Group = group_name,
          Parameter = as.character(Parameter),
          NSC = as.numeric(NSC),
          absNSC = abs(NSC),
          Direction = case_when(
            NSC > 0 ~ "Positive",
            NSC < 0 ~ "Negative",
            TRUE    ~ "Zero"
          )
        ) %>%
        filter(
          is.finite(NSC),
          absNSC >= cutoff
        ) %>%
        arrange(desc(absNSC))
      
      if (!is.null(top_n_per_group)) {
        out <- out %>%
          slice_head(n = top_n_per_group)
      }
      
      out
    }
  )
  
  if (nrow(plot_data) == 0) {
    stop(
      paste0(
        "No parameters met the cutoff of |NSC| >= ",
        cutoff,
        ". Reduce the cutoff."
      )
    )
  }
  
  plot_data <- plot_data %>%
    mutate(
      Group = factor(Group, levels = group_order),
      is_empty = FALSE
    )
  
  # ----------------------------------------------------------
  # 3. Add empty bars between sectors
  # ----------------------------------------------------------
  empty_data <- tibble(
    Group = factor(
      rep(group_order, each = empty_bar),
      levels = group_order
    ),
    Parameter = NA_character_,
    NSC = NA_real_,
    absNSC = 0,
    Direction = NA_character_,
    is_empty = TRUE
  )
  
  circle_data <- bind_rows(plot_data, empty_data) %>%
    arrange(Group, is_empty, desc(absNSC)) %>%
    mutate(id = row_number())
  
  number_of_bars <- nrow(circle_data)
  
  # ----------------------------------------------------------
  # 4. Parameter-label positions and angles
  # ----------------------------------------------------------
  label_offset <- max(circle_data$absNSC, na.rm = TRUE) * 0.08
  
  label_data <- circle_data %>%
    filter(!is_empty) %>%
    mutate(
      angle = 90 - 360 * (id - 0.5) / number_of_bars,
      hjust = ifelse(angle < -90, 1, 0),
      angle = ifelse(angle < -90, angle + 180, angle),
      label_y = absNSC + label_offset
    )
  
  # ----------------------------------------------------------
  # 5. Group-sector positions
  # ----------------------------------------------------------
  base_data <- circle_data %>%
    group_by(Group) %>%
    summarise(
      start = min(id),
      end = max(id[!is_empty]),
      title_position = mean(c(start, end)),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------
  # 6. Reference-ring segments
  # ----------------------------------------------------------
  ring_data <- tidyr::crossing(
    base_data %>% select(Group, start, end),
    Ring = ring_breaks
  )
  
  # Dynamic radial limits
  maximum_value <- max(
    circle_data$absNSC,
    ring_breaks,
    na.rm = TRUE
  )
  
  outer_limit <- maximum_value * 1.35
  inner_limit <- -maximum_value * 0.45
  
  # Position used to label rings
  ring_label_x <- max(circle_data$id)
  
  # ----------------------------------------------------------
  # 7. Circular bar plot
  # ----------------------------------------------------------
  p <- ggplot(
    circle_data,
    aes(
      x = id,
      y = absNSC,
      fill = Direction
    )
  ) +
    
    # Reference rings, drawn before bars
    geom_segment(
      data = ring_data,
      aes(
        x = start - 0.45,
        xend = end + 0.45,
        y = Ring,
        yend = Ring
      ),
      inherit.aes = FALSE,
      linewidth = 0.35,
      color = "grey70"
    ) +
    
    # Bars
    geom_col(
      width = 0.86,
      alpha = 0.85,
      na.rm = TRUE
    ) +
    
    # Parameter labels
    geom_text(
      data = label_data,
      aes(
        x = id,
        y = label_y,
        label = Parameter,
        angle = angle,
        hjust = hjust
      ),
      inherit.aes = FALSE,
      size = label_size,
      color = "black",
      fontface = "bold"
    ) +
    
    # Black arc under each group
    geom_segment(
      data = base_data,
      aes(
        x = start - 0.35,
        xend = end + 0.35,
        y = -maximum_value * 0.035,
        yend = -maximum_value * 0.035
      ),
      inherit.aes = FALSE,
      linewidth = 0.9,
      color = "black"
    ) +
    
    # Group labels
    geom_text(
      data = base_data,
      aes(
        x = title_position,
        y = -maximum_value * 0.26,
        label = Group
      ),
      inherit.aes = FALSE,
      size = group_label_size,
      fontface = "bold",
      color = "black"
    ) +
    
    # Ring-value labels
    annotate(
      "text",
      x = ring_label_x,
      y = ring_breaks,
      label = paste0(ring_breaks*100, "%"),
      size = 5,
      fontface = "bold",
      colour = "red",
      hjust = 1
    ) +
    
    scale_fill_manual(
      values = c(
        "Positive" = "red",
        "Negative" = "blue",
        "Zero" = "grey70"
      ),
      breaks = c("Negative", "Positive"),
      drop = FALSE
    ) +
    
    scale_x_continuous(
      limits = c(0.5, number_of_bars + 0.5),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(inner_limit, outer_limit),
      expand = c(0, 0)
    ) +
    
    coord_polar(start = 0) +
    
    labs(
      title = title,
      fill = "Effect"
    ) +
    
    theme_void(base_size = 13) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16,
        margin = margin(b = 10)
      ),
      legend.position = "bottom",
      legend.margin = margin(t = -10),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}



p_circle <- plot_circular_nsc(
  nsc_list = list(
    Brain  = brain_nsc,
    Kidney = kidney_nsc,
    Plasma = brain_plasma_nsc
  ),
  cutoff = 0.30,
  top_n_per_group = NULL,
  empty_bar = 3,
  ring_breaks = c(0.3, 0.6, 0.9),
  title = "Normalized sensitivity of AUC (0–168 h)"
)

print(p_circle)

p_circle_BK <- plot_circular_nsc(
  nsc_list = list(
    Brain  = brain_nsc,
    Kidney = kidney_nsc
  ),
  cutoff = 0.30,
  empty_bar = 4,
  ring_breaks = c(0.3, 0.6, 0.9),
  title = "Normalized sensitivity of brain and kidney AUC (0–168 h)"
)

print(p_circle_BK)


p_circle_top10 <- plot_circular_nsc(
  nsc_list = list(
    Brain  = brain_nsc,
    Kidney = kidney_nsc,
    Plasma = brain_plasma_nsc
  ),
  cutoff = 0,
  top_n_per_group = 10,
  empty_bar = 3,
  ring_breaks = c(0.3, 0.6, 0.9),
  title = "Top 10 normalized sensitivity coefficients"
)

print(p_circle_top10)
# ============================================================
# OPTIONAL SAVE RESULTS
# ============================================================

# write.csv(brain_nsc, "BrainModel_BrainAUC_NSC.csv", row.names = FALSE)
# write.csv(kidney_nsc, "KidneyModel_KidneyAUC_NSC.csv", row.names = FALSE)

# if (!is.null(brain_plasma_nsc)) {
#   write.csv(brain_plasma_nsc, "BrainModel_PlasmaAUC_NSC.csv", row.names = FALSE)
# }
# if (!is.null(kidney_plasma_nsc)) {
#   write.csv(kidney_plasma_nsc, "KidneyModel_PlasmaAUC_NSC.csv", row.names = FALSE)
# }

# ggsave("BrainModel_BrainAUC_NSC.png", p_brain, width = 8, height = 6, dpi = 300)
# ggsave("KidneyModel_KidneyAUC_NSC.png", p_kidney, width = 8, height = 6, dpi = 300)