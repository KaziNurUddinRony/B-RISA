# ================================================================
# B-RISA Stationary Simulation Data Generation
#
# Generates the 576-scenario clinically calibrated factorial design used
# in the manuscript. The script saves one generated-data CSV per sample
# size (n = 100, 500, 1000), together with scenario/calibration summaries.
#
# Run from the repository root, from the R/ directory, or set BRISA_ROOT
# to the repository root before running.
# ================================================================

rm(list = ls())
gc()

# ------------------------------------------------
# 0. Project directory
# ------------------------------------------------

get_brisa_root <- function() {
  env_root <- Sys.getenv("BRISA_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }

  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (basename(wd) == "R") dirname(wd) else wd
}

PROJECT_ROOT <- get_brisa_root()

PROJECT_DIR <- PROJECT_ROOT
setwd(PROJECT_DIR)

# ------------------------------------------------
# 1. Output folders
# ------------------------------------------------

RESULT_ROOT <- "results/clinical_calibrated_simulation_by_n"
FIGURE_ROOT <- "figures/clinical_calibrated_simulation_by_n"

DATA_DIR <- file.path(RESULT_ROOT, "data_by_sample_size")
SUMMARY_DIR <- file.path(RESULT_ROOT, "summaries")
SCENARIO_DIR <- file.path(RESULT_ROOT, "scenario_grid")
PARAMETER_DIR <- file.path(RESULT_ROOT, "parameter_tables")
EXAMPLE_DIR <- file.path(RESULT_ROOT, "examples")

FIG_EXAMPLE_DIR <- file.path(FIGURE_ROOT, "example_pairplots")
FIG_PREVALENCE_DIR <- file.path(FIGURE_ROOT, "prevalence_checks")
FIG_MISSINGNESS_DIR <- file.path(FIGURE_ROOT, "missingness_checks")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SCENARIO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PARAMETER_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EXAMPLE_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(FIG_EXAMPLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_PREVALENCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_MISSINGNESS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# 2. Packages
# ------------------------------------------------

required_packages <- c(
  "MASS",
  "Matrix",
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "GGally",
  "tibble",
  "stringr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Restore the repository R environment before running this script."
  )
}

library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)
library(tibble)
library(stringr)

# ------------------------------------------------
# 3. User settings
# ------------------------------------------------

BASE_SEED <- 20260616

# Keep this 1 first.
# If you set this to 100, files will become very large.
N_REPS_TO_SAVE <- 100

MAX_P <- 20

# Final generated-data files
FILE_N100 <- file.path(DATA_DIR, "generated_data_n100.csv")
FILE_N500 <- file.path(DATA_DIR, "generated_data_n500.csv")
FILE_N1000 <- file.path(DATA_DIR, "generated_data_n1000.csv")

# ------------------------------------------------
# 4. Derived calibration-input paths
# ------------------------------------------------

PARAM_DIR <- file.path(PROJECT_ROOT, "data", "derived")

adult_mu_file <- file.path(PARAM_DIR, "adult_corrected_mean_vectors.csv")
adult_sigma0_file <- file.path(PARAM_DIR, "adult_corrected_Sigma_Y0_no_appendicitis.csv")
adult_sigma1_file <- file.path(PARAM_DIR, "adult_corrected_Sigma_Y1_appendicitis.csv")
adult_sigmap_file <- file.path(PARAM_DIR, "adult_corrected_Sigma_pooled.csv")

rpa_mu_file <- file.path(PARAM_DIR, "rpa_standardized_corrected_mean_vectors.csv")
rpa_sigma0_file <- file.path(PARAM_DIR, "rpa_standardized_corrected_Sigma_Y0_no_appendicitis.csv")
rpa_sigma1_file <- file.path(PARAM_DIR, "rpa_standardized_corrected_Sigma_Y1_appendicitis.csv")
rpa_sigmap_file <- file.path(PARAM_DIR, "rpa_standardized_corrected_Sigma_pooled.csv")

needed_files <- c(
  adult_mu_file,
  adult_sigma0_file,
  adult_sigma1_file,
  adult_sigmap_file,
  rpa_mu_file,
  rpa_sigma0_file,
  rpa_sigma1_file,
  rpa_sigmap_file
)

missing_needed <- needed_files[!file.exists(needed_files)]

if (length(missing_needed) > 0) {
  stop(
    "The following corrected EDA files are missing:\n",
    paste(missing_needed, collapse = "\n"),
    "\nCreate or restore the derived calibration files before running this script."
  )
}

# ------------------------------------------------
# 5. Helper functions for reading mu and Sigma
# ------------------------------------------------

read_mean_vectors <- function(path) {
  df <- data.table::fread(
    path,
    data.table = FALSE
  )
  
  required_cols <- c(
    "Feature",
    "mu_Y0_no_appendicitis",
    "mu_Y1_appendicitis"
  )
  
  if (!all(required_cols %in% names(df))) {
    stop("Mean vector file does not contain required columns: ", path)
  }
  
  list(
    features = df$Feature,
    mu0 = as.numeric(df$mu_Y0_no_appendicitis),
    mu1 = as.numeric(df$mu_Y1_appendicitis)
  )
}

read_matrix_csv <- function(path) {
  df <- data.table::fread(
    path,
    data.table = FALSE
  )
  
  feature_names <- df[[1]]
  
  mat <- as.matrix(
    df[, -1, drop = FALSE]
  )
  
  storage.mode(mat) <- "numeric"
  
  rownames(mat) <- feature_names
  colnames(mat) <- names(df)[-1]
  
  mat
}

nearest_pd <- function(Sigma, ridge = 1e-6) {
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  
  eig <- tryCatch(
    eigen(
      Sigma,
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    error = function(e) NA_real_
  )
  
  if (all(is.finite(eig)) && min(eig) > 1e-8) {
    return(Sigma)
  }
  
  out <- as.matrix(
    Matrix::nearPD(
      Sigma,
      corr = FALSE
    )$mat
  )
  
  out <- out + diag(ridge, ncol(out))
  
  out
}

load_empirical_parameters <- function(
    dataset_name,
    mu_file,
    sigma0_file,
    sigma1_file,
    sigmap_file
) {
  mu <- read_mean_vectors(mu_file)
  
  Sigma0 <- read_matrix_csv(sigma0_file)
  Sigma1 <- read_matrix_csv(sigma1_file)
  Sigmap <- read_matrix_csv(sigmap_file)
  
  Sigma0 <- nearest_pd(Sigma0)
  Sigma1 <- nearest_pd(Sigma1)
  Sigmap <- nearest_pd(Sigmap)
  
  list(
    dataset_name = dataset_name,
    features = mu$features,
    p_base = length(mu$features),
    mu0 = mu$mu0,
    mu1 = mu$mu1,
    Sigma0 = Sigma0,
    Sigma1 = Sigma1,
    Sigmap = Sigmap
  )
}

adult_params <- load_empirical_parameters(
  dataset_name = "Adult",
  mu_file = adult_mu_file,
  sigma0_file = adult_sigma0_file,
  sigma1_file = adult_sigma1_file,
  sigmap_file = adult_sigmap_file
)

rpa_params <- load_empirical_parameters(
  dataset_name = "RPA_standardized",
  mu_file = rpa_mu_file,
  sigma0_file = rpa_sigma0_file,
  sigma1_file = rpa_sigma1_file,
  sigmap_file = rpa_sigmap_file
)

cat("\nAdult empirical parameters loaded:\n")
print(adult_params$features)

cat("\nRPA standardized empirical parameters loaded:\n")
print(rpa_params$features)

# ------------------------------------------------
# 6. Expand empirical adult parameters to p = 10 and p = 20
# ------------------------------------------------

expand_empirical_params <- function(
    base_params,
    p_target,
    cross_block_attenuation = 0.30,
    diagonal_jitter = 1e-6
) {
  m <- base_params$p_base
  
  if (p_target == m) {
    return(base_params)
  }
  
  index_map <- ((seq_len(p_target) - 1) %% m) + 1
  block_id <- ((seq_len(p_target) - 1) %/% m) + 1
  
  expand_vector <- function(v) {
    v[index_map]
  }
  
  expand_matrix <- function(Sigma_base) {
    Sigma_p <- matrix(
      0,
      nrow = p_target,
      ncol = p_target
    )
    
    for (i in seq_len(p_target)) {
      for (j in seq_len(p_target)) {
        original_i <- index_map[i]
        original_j <- index_map[j]
        
        block_distance <- abs(block_id[i] - block_id[j])
        
        attenuation <- cross_block_attenuation ^ block_distance
        
        Sigma_p[i, j] <- Sigma_base[original_i, original_j] * attenuation
      }
    }
    
    Sigma_p <- nearest_pd(
      Sigma_p + diag(diagonal_jitter, p_target)
    )
    
    Sigma_p
  }
  
  expanded_features <- paste0(
    base_params$features[index_map],
    "_copy",
    block_id
  )
  
  list(
    dataset_name = paste0(base_params$dataset_name, "_expanded_p", p_target),
    features = expanded_features,
    p_base = p_target,
    mu0 = expand_vector(base_params$mu0),
    mu1 = expand_vector(base_params$mu1),
    Sigma0 = expand_matrix(base_params$Sigma0),
    Sigma1 = expand_matrix(base_params$Sigma1),
    Sigmap = expand_matrix(base_params$Sigmap)
  )
}

get_params_for_p <- function(p) {
  if (p == 3) {
    return(rpa_params)
  }
  
  if (p == 5) {
    return(adult_params)
  }
  
  if (p %in% c(10, 20)) {
    return(
      expand_empirical_params(
        base_params = adult_params,
        p_target = p,
        cross_block_attenuation = 0.30
      )
    )
  }
  
  stop("Unsupported p: ", p)
}

# ------------------------------------------------
# 7. Scenario grid: 576 scenarios
# ------------------------------------------------

scenario_grid <- tidyr::crossing(
  n_train = c(100, 500, 1000),
  p = c(3, 5, 10, 20),
  prevalence = c(0.20, 0.40, 0.60),
  boundary = c("linear", "nonlinear"),
  covariance = c("pooled", "class_specific"),
  label_noise = c(0, 0.10),
  missing_rate = c(0, 0.15)
) |>
  dplyr::mutate(
    scenario_id = dplyr::row_number(),
    calibration_source = dplyr::case_when(
      p == 3 ~ "RPA_standardized",
      p == 5 ~ "Adult",
      p %in% c(10, 20) ~ "Adult_expanded",
      TRUE ~ "Unknown"
    ),
    scenario_name = paste0(
      "S", scenario_id,
      "_n", n_train,
      "_p", p,
      "_prev", prevalence,
      "_", boundary,
      "_", covariance,
      "_ln", label_noise,
      "_miss", missing_rate,
      "_", calibration_source
    )
  ) |>
  dplyr::select(
    scenario_id,
    scenario_name,
    calibration_source,
    dplyr::everything()
  )

stopifnot(nrow(scenario_grid) == 576)

data.table::fwrite(
  scenario_grid,
  file.path(SCENARIO_DIR, "scenario_grid_576_clinical_calibrated.csv")
)

cat("\nScenario grid created with", nrow(scenario_grid), "scenarios.\n")
print(head(scenario_grid, 12))

# ------------------------------------------------
# 8. Core data-generation functions
# ------------------------------------------------

logit <- function(p) {
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

make_exact_binary_labels <- function(n, prevalence, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n1 <- round(n * prevalence)
  n0 <- n - n1
  
  y <- c(
    rep(0L, n0),
    rep(1L, n1)
  )
  
  sample(
    y,
    size = n,
    replace = FALSE
  )
}

generate_gaussian_by_class <- function(
    n,
    prevalence,
    mu0,
    mu1,
    Sigma0,
    Sigma1,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  p <- length(mu0)
  
  y <- make_exact_binary_labels(
    n = n,
    prevalence = prevalence,
    seed = seed + 1
  )
  
  X <- matrix(
    NA_real_,
    nrow = n,
    ncol = p
  )
  
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  
  if (length(idx0) > 0) {
    X[idx0, ] <- MASS::mvrnorm(
      n = length(idx0),
      mu = mu0,
      Sigma = Sigma0
    )
  }
  
  if (length(idx1) > 0) {
    X[idx1, ] <- MASS::mvrnorm(
      n = length(idx1),
      mu = mu1,
      Sigma = Sigma1
    )
  }
  
  list(
    X = X,
    y = as.integer(y)
  )
}

apply_risa_style_nonlinearity <- function(X, y) {
  X_new <- X
  
  p <- ncol(X_new)
  
  if (p >= 1) {
    idx1 <- which(y == 1)
    X_new[idx1, 1] <- X_new[idx1, 1]^2
  }
  
  if (p >= 2) {
    idx0 <- which(y == 0)
    X_new[idx0, 2] <- sin(X_new[idx0, 2])
  }
  
  if (p >= 3) {
    X_new[, 3] <- X_new[, 3] + 0.25 * X_new[, 1] * X_new[, 2]
  }
  
  X_new
}

apply_label_noise <- function(y, label_noise, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (label_noise <= 0) {
    return(as.integer(y))
  }
  
  flip <- rbinom(
    n = length(y),
    size = 1,
    prob = label_noise
  )
  
  y_noisy <- ifelse(
    flip == 1,
    1 - y,
    y
  )
  
  as.integer(y_noisy)
}

apply_missingness <- function(X, y, missing_rate, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (missing_rate <= 0) {
    return(X)
  }
  
  X_miss <- X
  
  p <- ncol(X_miss)
  affected_features <- seq_len(min(3, p))
  
  z1 <- as.numeric(scale(X_miss[, 1]))
  z1[is.na(z1)] <- 0
  
  for (j in affected_features) {
    prob_missing <- inv_logit(
      logit(missing_rate) +
        0.35 * z1 -
        0.30 * y
    )
    
    r_missing <- rbinom(
      n = nrow(X_miss),
      size = 1,
      prob = prob_missing
    )
    
    X_miss[r_missing == 1, j] <- NA_real_
  }
  
  X_miss
}

generate_one_dataset <- function(
    n,
    p,
    prevalence,
    boundary,
    covariance,
    label_noise,
    missing_rate,
    seed = 1
) {
  params <- get_params_for_p(p)
  
  if (covariance == "pooled") {
    Sigma0 <- params$Sigmap
    Sigma1 <- params$Sigmap
  } else if (covariance == "class_specific") {
    Sigma0 <- params$Sigma0
    Sigma1 <- params$Sigma1
  } else {
    stop("Unknown covariance setting: ", covariance)
  }
  
  generated <- generate_gaussian_by_class(
    n = n,
    prevalence = prevalence,
    mu0 = params$mu0,
    mu1 = params$mu1,
    Sigma0 = Sigma0,
    Sigma1 = Sigma1,
    seed = seed + 1
  )
  
  X <- generated$X
  y_true <- generated$y
  
  if (boundary == "nonlinear") {
    X <- apply_risa_style_nonlinearity(
      X = X,
      y = y_true
    )
  }
  
  y_observed <- apply_label_noise(
    y = y_true,
    label_noise = label_noise,
    seed = seed + 2
  )
  
  X <- apply_missingness(
    X = X,
    y = y_observed,
    missing_rate = missing_rate,
    seed = seed + 3
  )
  
  colnames(X) <- paste0("X", seq_len(p))
  
  df <- as.data.frame(X)
  
  df$Y_true <- y_true
  df$Y_observed <- y_observed
  df$Outcome <- factor(
    y_observed,
    levels = c(0, 1),
    labels = c("No_appendicitis", "Appendicitis")
  )
  
  df
}

# ------------------------------------------------
# 9. Helper: pad feature columns to X1,...,X20
# ------------------------------------------------

pad_feature_columns <- function(df, max_p = MAX_P) {
  all_feature_cols <- paste0("X", seq_len(max_p))
  
  for (feature_name in all_feature_cols) {
    if (!feature_name %in% names(df)) {
      df[[feature_name]] <- NA_real_
    }
  }
  
  df |>
    dplyr::select(
      dplyr::all_of(all_feature_cols),
      dplyr::everything()
    )
}

add_scenario_metadata <- function(df, scenario_row, rep_id) {
  df <- pad_feature_columns(
    df = df,
    max_p = MAX_P
  )
  
  df <- df |>
    dplyr::mutate(
      scenario_id = scenario_row$scenario_id,
      scenario_name = scenario_row$scenario_name,
      rep_id = rep_id,
      n_train = scenario_row$n_train,
      p = scenario_row$p,
      prevalence_target = scenario_row$prevalence,
      boundary = scenario_row$boundary,
      covariance = scenario_row$covariance,
      label_noise = scenario_row$label_noise,
      missing_rate_target = scenario_row$missing_rate,
      calibration_source = scenario_row$calibration_source,
      row_id_within_scenario = dplyr::row_number(),
      .before = 1
    )
  
  df
}

summarize_dataset <- function(df, scenario_row, rep_id) {
  feature_cols <- paste0("X", seq_len(scenario_row$p))
  affected_cols <- paste0("X", seq_len(min(3, scenario_row$p)))
  
  feature_matrix <- as.matrix(
    df[, feature_cols, drop = FALSE]
  )
  
  affected_matrix <- as.matrix(
    df[, affected_cols, drop = FALSE]
  )
  
  tibble::tibble(
    scenario_id = scenario_row$scenario_id,
    scenario_name = scenario_row$scenario_name,
    rep_id = rep_id,
    n = nrow(df),
    p = scenario_row$p,
    prevalence_target_true = scenario_row$prevalence,
    prevalence_true_observed = mean(df$Y_true == 1),
    prevalence_after_label_noise = mean(df$Y_observed == 1),
    expected_noisy_prevalence = scenario_row$prevalence * (1 - scenario_row$label_noise) +
      (1 - scenario_row$prevalence) * scenario_row$label_noise,
    missing_rate_target = scenario_row$missing_rate,
    missing_rate_observed_all_features = mean(is.na(feature_matrix)),
    missing_rate_observed_affected_features = mean(is.na(affected_matrix)),
    boundary = scenario_row$boundary,
    covariance = scenario_row$covariance,
    label_noise = scenario_row$label_noise,
    calibration_source = scenario_row$calibration_source
  )
}

# ------------------------------------------------
# 10. Generate all scenarios and save THREE CSV files by sample size
# ------------------------------------------------

all_generated_data <- list()
all_summaries <- list()

counter <- 1

cat("\n================================================\n")
cat("Generating all scenarios and storing in memory\n")
cat("================================================\n")

for (i in seq_len(nrow(scenario_grid))) {
  scenario_i <- scenario_grid[i, ]
  
  if (i %% 25 == 0) {
    cat("Generating scenario", i, "of", nrow(scenario_grid), "\n")
  }
  
  for (rep_id in seq_len(N_REPS_TO_SAVE)) {
    seed_i <- BASE_SEED + scenario_i$scenario_id * 10000 + rep_id * 100
    
    df_i <- generate_one_dataset(
      n = scenario_i$n_train,
      p = scenario_i$p,
      prevalence = scenario_i$prevalence,
      boundary = scenario_i$boundary,
      covariance = scenario_i$covariance,
      label_noise = scenario_i$label_noise,
      missing_rate = scenario_i$missing_rate,
      seed = seed_i
    )
    
    all_summaries[[counter]] <- summarize_dataset(
      df = df_i,
      scenario_row = scenario_i,
      rep_id = rep_id
    )
    
    df_i_meta <- add_scenario_metadata(
      df = df_i,
      scenario_row = scenario_i,
      rep_id = rep_id
    )
    
    all_generated_data[[counter]] <- data.table::as.data.table(df_i_meta)
    
    counter <- counter + 1
  }
}

summary_df <- dplyr::bind_rows(all_summaries)

combined_data <- data.table::rbindlist(
  all_generated_data,
  use.names = TRUE,
  fill = TRUE
)

# Reorder columns
metadata_cols <- c(
  "scenario_id",
  "scenario_name",
  "rep_id",
  "row_id_within_scenario",
  "n_train",
  "p",
  "prevalence_target",
  "boundary",
  "covariance",
  "label_noise",
  "missing_rate_target",
  "calibration_source"
)

feature_cols <- paste0("X", seq_len(MAX_P))

outcome_cols <- c(
  "Y_true",
  "Y_observed",
  "Outcome"
)

combined_data <- combined_data |>
  dplyr::select(
    dplyr::all_of(metadata_cols),
    dplyr::all_of(feature_cols),
    dplyr::all_of(outcome_cols)
  )

cat("\nCombined data dimensions before splitting:\n")
print(dim(combined_data))

# Split by n_train
data_n100 <- combined_data |>
  dplyr::filter(n_train == 100)

data_n500 <- combined_data |>
  dplyr::filter(n_train == 500)

data_n1000 <- combined_data |>
  dplyr::filter(n_train == 1000)

# Save THREE CSV files
data.table::fwrite(
  data_n100,
  FILE_N100
)

data.table::fwrite(
  data_n500,
  FILE_N500
)

data.table::fwrite(
  data_n1000,
  FILE_N1000
)

# Save summary CSV
data.table::fwrite(
  summary_df,
  file.path(SUMMARY_DIR, "generated_data_summary.csv")
)

cat("\nSaved generated datasets by sample size:\n")
cat("  ", FILE_N100, "\n")
cat("  ", FILE_N500, "\n")
cat("  ", FILE_N1000, "\n")

cat("\nGenerated data summary head:\n")
print(head(summary_df, 20))

# ------------------------------------------------
# 11. Save parameter tables used for p = 3, 5, 10, 20
# ------------------------------------------------

save_parameter_table_for_p <- function(p) {
  params <- get_params_for_p(p)
  
  feature_names <- paste0("X", seq_len(p))
  
  mu_table <- tibble::tibble(
    Feature = feature_names,
    Source_Feature = params$features,
    mu_Y0_no_appendicitis = params$mu0,
    mu_Y1_appendicitis = params$mu1,
    delta_appendicitis_minus_no_appendicitis = params$mu1 - params$mu0
  )
  
  data.table::fwrite(
    mu_table,
    file.path(
      PARAMETER_DIR,
      paste0("p_", p, "_mean_vectors.csv")
    )
  )
  
  save_mat <- function(mat, name) {
    rownames(mat) <- feature_names
    colnames(mat) <- feature_names
    
    out <- as.data.frame(mat) |>
      tibble::rownames_to_column("Feature")
    
    data.table::fwrite(
      out,
      file.path(
        PARAMETER_DIR,
        paste0("p_", p, "_", name, ".csv")
      )
    )
  }
  
  save_mat(params$Sigma0, "Sigma_Y0_no_appendicitis")
  save_mat(params$Sigma1, "Sigma_Y1_appendicitis")
  save_mat(params$Sigmap, "Sigma_pooled")
  
  invisible(mu_table)
}

param_tables <- lapply(
  c(3, 5, 10, 20),
  save_parameter_table_for_p
)

# ------------------------------------------------
# 12. p = 3 example shown in output cell
# ------------------------------------------------

example_p3 <- scenario_grid |>
  dplyr::filter(
    n_train == 100,
    p == 3,
    prevalence == 0.20,
    boundary == "linear",
    covariance == "class_specific",
    label_noise == 0,
    missing_rate == 0
  ) |>
  dplyr::slice(1)

example_p3_data <- generate_one_dataset(
  n = example_p3$n_train,
  p = example_p3$p,
  prevalence = example_p3$prevalence,
  boundary = example_p3$boundary,
  covariance = example_p3$covariance,
  label_noise = example_p3$label_noise,
  missing_rate = example_p3$missing_rate,
  seed = 12345
)

cat("\n================================================\n")
cat("p = 3 example scenario\n")
cat("================================================\n")
print(example_p3)

cat("\nFirst 10 rows of generated p = 3 data:\n")
print(head(example_p3_data, 10))

cat("\nObserved class counts for p = 3 example:\n")
print(table(example_p3_data$Outcome))

cat("\nEmpirical mean by generated class for p = 3 example:\n")
print(
  example_p3_data |>
    dplyr::group_by(Outcome) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::starts_with("X"),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
)

cat("\nEmpirical covariance for generated no-appendicitis group, p = 3 example:\n")
print(
  round(
    cov(
      example_p3_data |>
        dplyr::filter(Y_observed == 0) |>
        dplyr::select(X1, X2, X3)
    ),
    4
  )
)

cat("\nEmpirical covariance for generated appendicitis group, p = 3 example:\n")
print(
  round(
    cov(
      example_p3_data |>
        dplyr::filter(Y_observed == 1) |>
        dplyr::select(X1, X2, X3)
    ),
    4
  )
)

data.table::fwrite(
  example_p3_data,
  file.path(EXAMPLE_DIR, "example_p3_generated_dataset.csv")
)

# ------------------------------------------------
# 13. p = 3 pairplot shown and saved
# ------------------------------------------------

plot_p3 <- example_p3_data |>
  dplyr::select(X1, X2, X3, Outcome) |>
  GGally::ggpairs(
    columns = 1:3,
    mapping = ggplot2::aes(linetype = Outcome),
    upper = list(
      continuous = GGally::wrap("cor", size = 3)
    ),
    lower = list(
      continuous = GGally::wrap("points", alpha = 0.6, size = 1)
    ),
    diag = list(
      continuous = GGally::wrap("densityDiag", alpha = 0.4)
    ),
    title = "Example p = 3 Clinical-Calibrated Simulation"
  ) +
  ggplot2::theme_bw(base_size = 10)

ggplot2::ggsave(
  file.path(FIG_EXAMPLE_DIR, "example_p3_pairplot.pdf"),
  plot_p3,
  width = 8,
  height = 8,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIG_EXAMPLE_DIR, "example_p3_pairplot.png"),
  plot_p3,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(plot_p3)

# ------------------------------------------------
# 14. Prevalence validation plot
# ------------------------------------------------

plot_prev <- summary_df |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = prevalence_target_true,
      y = prevalence_after_label_noise,
      shape = calibration_source
    )
  ) +
  ggplot2::geom_point(alpha = 0.45) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::facet_wrap(~ n) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Observed prevalence after label-noise mechanism",
    subtitle = "Dashed line represents true target prevalence before label noise",
    x = "Target true prevalence before label noise",
    y = "Observed prevalence after label noise"
  )

ggplot2::ggsave(
  file.path(FIG_PREVALENCE_DIR, "observed_prevalence_summary.pdf"),
  plot_prev,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIG_PREVALENCE_DIR, "observed_prevalence_summary.png"),
  plot_prev,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(plot_prev)

# ------------------------------------------------
# 15. Missingness validation plots
# ------------------------------------------------

plot_missing_all <- summary_df |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(missing_rate_target),
      y = missing_rate_observed_all_features,
      shape = calibration_source
    )
  ) +
  ggplot2::geom_boxplot(outlier.alpha = 0.3) +
  ggplot2::facet_wrap(~ n) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Observed missingness across all feature columns",
    x = "Target missingness",
    y = "Observed missingness across all features"
  )

ggplot2::ggsave(
  file.path(FIG_MISSINGNESS_DIR, "observed_missingness_all_features_summary.pdf"),
  plot_missing_all,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIG_MISSINGNESS_DIR, "observed_missingness_all_features_summary.png"),
  plot_missing_all,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(plot_missing_all)

plot_missing_affected <- summary_df |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(missing_rate_target),
      y = missing_rate_observed_affected_features,
      shape = calibration_source
    )
  ) +
  ggplot2::geom_boxplot(outlier.alpha = 0.3) +
  ggplot2::facet_wrap(~ n) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title = "Observed missingness among affected features only",
    subtitle = "Missingness was applied only to X1, ..., Xmin(3,p)",
    x = "Target missingness",
    y = "Observed missingness among affected features"
  )

ggplot2::ggsave(
  file.path(FIG_MISSINGNESS_DIR, "observed_missingness_affected_features_summary.pdf"),
  plot_missing_affected,
  width = 10,
  height = 7,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIG_MISSINGNESS_DIR, "observed_missingness_affected_features_summary.png"),
  plot_missing_affected,
  width = 10,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

print(plot_missing_affected)

# ------------------------------------------------
# 16. Final checks
# ------------------------------------------------

cat("\n================================================\n")
cat("FINAL CHECKS\n")
cat("================================================\n")

cat("\nNumber of scenarios:\n")
print(nrow(scenario_grid))

cat("\nRows in n=100 file:\n")
print(nrow(data_n100))

cat("\nRows in n=500 file:\n")
print(nrow(data_n500))

cat("\nRows in n=1000 file:\n")
print(nrow(data_n1000))

expected_n100 <- 192 * 100 * N_REPS_TO_SAVE
expected_n500 <- 192 * 500 * N_REPS_TO_SAVE
expected_n1000 <- 192 * 1000 * N_REPS_TO_SAVE

cat("\nExpected rows:\n")
cat("n=100:  ", expected_n100, "\n")
cat("n=500:  ", expected_n500, "\n")
cat("n=1000: ", expected_n1000, "\n")

if (nrow(data_n100) != expected_n100) {
  warning("n=100 file row count does not match expected.")
} else {
  cat("n=100 row count is correct.\n")
}

if (nrow(data_n500) != expected_n500) {
  warning("n=500 file row count does not match expected.")
} else {
  cat("n=500 row count is correct.\n")
}

if (nrow(data_n1000) != expected_n1000) {
  warning("n=1000 file row count does not match expected.")
} else {
  cat("n=1000 row count is correct.\n")
}

cat("\nSummary of generated rows by n_train:\n")
print(
  combined_data |>
    dplyr::count(n_train)
)

cat("\nSummary of generated rows by p:\n")
print(
  combined_data |>
    dplyr::count(p)
)

cat("\nSummary of generated rows by target prevalence:\n")
print(
  combined_data |>
    dplyr::count(prevalence_target)
)

# ------------------------------------------------
# 17. Final message
# ------------------------------------------------

cat("\n================================================\n")
cat("CLINICAL-CALIBRATED DATA GENERATION COMPLETE\n")
cat("================================================\n")
cat("Generated data CSV files by sample size:\n")
cat("  ", FILE_N100, "\n")
cat("  ", FILE_N500, "\n")
cat("  ", FILE_N1000, "\n\n")

cat("Scenario grid:\n")
cat("  ", file.path(SCENARIO_DIR, "scenario_grid_576_clinical_calibrated.csv"), "\n\n")

cat("Summary:\n")
cat("  ", file.path(SUMMARY_DIR, "generated_data_summary.csv"), "\n\n")

cat("Parameter tables:\n")
cat("  ", PARAMETER_DIR, "\n\n")

cat("Example p=3 dataset:\n")
cat("  ", file.path(EXAMPLE_DIR, "example_p3_generated_dataset.csv"), "\n\n")

cat("Figures:\n")
cat("  ", FIG_EXAMPLE_DIR, "\n")
cat("  ", FIG_PREVALENCE_DIR, "\n")
cat("  ", FIG_MISSINGNESS_DIR, "\n")
cat("================================================\n")

