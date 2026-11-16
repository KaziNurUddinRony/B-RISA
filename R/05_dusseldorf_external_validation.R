# ================================================================
# B-RISA Locked Düsseldorf External Validation
#
# Fits each prespecified Bayesian model once to the full eligible RPA
# development cohort, applies RPA-derived preprocessing unchanged, and
# evaluates locked predictions in the 301-patient Düsseldorf cohort.
# The 81 complete cases are retained as a sensitivity analysis.
#
# Run from the repository root, from the R/ directory, or set BRISA_ROOT
# to the repository root before running.
# ================================================================

rm(list = ls())
graphics.off()
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 180,
  dplyr.summarise.inform = FALSE
)

# ------------------------------------------------
# 0. Project paths
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

RPA_STD_CLEAN_FILE <- file.path(
  PROJECT_ROOT, "data", "derived",
  "rpa_standardized_selected_cleaned_data.csv"
)

RPA_SCALING_PARAMETER_FILE <- file.path(
  PROJECT_ROOT, "data", "derived", "rpa_standardization_parameters.csv"
)

DUSSELDORF_RAW_FILE <- file.path(
  PROJECT_ROOT, "data", "raw", "app_data_ext.csv"
)

RESULT_ROOT <- file.path(PROJECT_ROOT, "analysis_outputs", "dusseldorf_external_validation")

TABLE_ROOT <- file.path(
  RESULT_ROOT,
  "tables"
)

R_OBJECT_ROOT <- file.path(
  RESULT_ROOT,
  "r_objects"
)

FIGURE_STEM_ROOT <- file.path(
  RESULT_ROOT,
  "figures"
)

FIGURE_PDF_ROOT <- file.path(
  FIGURE_STEM_ROOT,
  "pdf"
)

FIGURE_PNG_ROOT <- file.path(
  FIGURE_STEM_ROOT,
  "png"
)

DATASET_ID <- "dusseldorf_external"
DATASET_LABEL <- "Düsseldorf External Validation"

RPA_STD_FEATURES <- c(
  "Body_Temperature",
  "WBC_Count",
  "Neutrophil_Percentage"
)

# Exact upstream RPA standardization parameters used to create the
# standardized RPA dataset. These were estimated from the 677 complete
# RPA patients and are used only as a validated fallback when the saved
# parameter CSV is unavailable.
RPA_SCALING_FALLBACK <- data.frame(
  Feature = RPA_STD_FEATURES,
  mean_overall = c(
    37.39734,
    12.76388,
    71.83604
  ),
  sd_overall = c(
    0.916253,
    5.344710,
    14.459023
  )
)

# ------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------

CLEAN_OLD_OUTPUT <- FALSE
INSTALL_MISSING_PACKAGES <- FALSE
PRINT_LARGE_SUMMARY_TABLES <- TRUE

SEED <- 20260618L
EXTERNAL_FIT_SEED <- 20260730L
BOOTSTRAP_SEED <- 20260619L
BOOTSTRAP_REPS <- 2000L
CONFIDENCE_LEVEL <- 0.95
ALPHA <- 0.05
P_ADJUST_METHOD <- "holm"

MODEL_SET <- c("BLR", "BLDA", "BART")
MODEL_ORDER <- MODEL_SET

N_POST_DRAWS_BLR <- 300L
N_POST_DRAWS_BLDA <- 300L
N_POST_DRAWS_BART <- 300L
N_BURN_BART <- 150L
BART_NTREE <- 50L

BLR_INTERCEPT_PRIOR_SD <- 2.5
BLR_COEFFICIENT_PRIOR_SD <- 2.5
BLR_OPTIM_MAXIT <- 2000L
BLR_OPTIM_RELTOL <- 1e-10
BLR_MAX_ABS_GRADIENT <- 1e-5
BLR_HESSIAN_EIGEN_FLOOR <- 1e-8

# Genuine Bayesian LDA prior settings.
# All predictors have already undergone locked RPA-derived preprocessing,
# so an identity-centered covariance prior is interpretable on this scale.
BLDA_MEAN_PRIOR_KAPPA <- 0.01
BLDA_IW_DF_OFFSET <- 2L
BLDA_PRIOR_COVARIANCE_SCALE <- 1.0
BLDA_CLASS_PRIOR_ALPHA <- 0.5
BLDA_CLASS_PRIOR_BETA <- 0.5
BLDA_PRECISION_EIGEN_FLOOR <- 1e-10

CLASSIFICATION_THRESHOLD <- 0.50
SELECTIVE_COVERAGE <- 0.80
DCA_THRESHOLDS <- seq(0.01, 0.99, by = 0.01)
CALIBRATION_BINS <- 10L
EPS <- 1e-6

DEGENERATE_PREDICTION_TOL <- 1e-12
STOP_ON_MODEL_FAILURE <- TRUE
STOP_ON_DEGENERATE_PREDICTIONS <- TRUE

PRIMARY_ANALYSIS_ID <- "primary_imputed_all"
PRIMARY_ANALYSIS_LABEL <- "Primary: all patients with RPA-trained imputation"
COMPLETE_CASE_ANALYSIS_ID <- "complete_case"
COMPLETE_CASE_ANALYSIS_LABEL <- "Sensitivity: complete cases"

INFERENTIAL_METRICS <- c(
  "accuracy",
  "sensitivity",
  "specificity",
  "ppv",
  "npv",
  "f1",
  "auc",
  "auprc",
  "brier",
  "logloss",
  "ece",
  "selective_risk"
)

DESCRIPTIVE_METRICS <- c(
  "calibration_intercept",
  "calibration_slope",
  "mean_uncertainty",
  "mean_interval_width",
  "selective_coverage"
)

BOOTSTRAP_MODEL_METRICS <- unique(c(
  INFERENTIAL_METRICS,
  DESCRIPTIVE_METRICS
))

BOOTSTRAP_PAIRWISE_METRICS <- INFERENTIAL_METRICS

HIGHER_IS_BETTER <- c(
  "accuracy",
  "sensitivity",
  "specificity",
  "ppv",
  "npv",
  "f1",
  "auc",
  "auprc"
)

LOWER_IS_BETTER <- c(
  "brier",
  "logloss",
  "ece",
  "selective_risk"
)

model_colors <- c(
  "BLR" = "#1B9AAA",
  "BLDA" = "#F46036",
  "BART" = "#6A4C93",
  "Treat all" = "#444444",
  "Treat none" = "#999999"
)

outcome_colors <- c(
  "No_appendicitis" = "#2E86AB",
  "Appendicitis" = "#F18F01"
)

cohort_colors <- c(
  "RPA development" = "#4C78A8",
  "Düsseldorf external" = "#E45756"
)

# ------------------------------------------------
# 2. Packages
# ------------------------------------------------

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble",
  "scales",
  "pROC",
  "PRROC",
  "MASS",
  "janitor",
  "readr"
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
  if (!INSTALL_MISSING_PACKAGES) {
    stop(
      "Missing packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }

  install.packages(
    missing_packages,
    dependencies = TRUE
  )
}

if (!requireNamespace("BART", quietly = TRUE)) {
  if (INSTALL_MISSING_PACKAGES) {
    install.packages("BART", dependencies = TRUE)
  }
}

if (!requireNamespace("BART", quietly = TRUE)) {
  stop("The BART package is required.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(pROC)
  library(PRROC)
  library(MASS)
  library(janitor)
  library(readr)
})

# ------------------------------------------------
# 3. Output folders
# ------------------------------------------------

if (CLEAN_OLD_OUTPUT) {
  unlink(RESULT_ROOT, recursive = TRUE, force = TRUE)
  unlink(FIGURE_STEM_ROOT, recursive = TRUE, force = TRUE)
}

dirs <- list(
  result_base = RESULT_ROOT,
  figure_base = FIGURE_STEM_ROOT,
  dataset = file.path(TABLE_ROOT, "dataset"),
  preprocessing = file.path(TABLE_ROOT, "preprocessing"),
  cohort_shift = file.path(TABLE_ROOT, "cohort_shift"),
  metrics = file.path(TABLE_ROOT, "metrics"),
  predictions = file.path(TABLE_ROOT, "predictions"),
  status = file.path(TABLE_ROOT, "status"),
  calibration = file.path(TABLE_ROOT, "calibration"),
  selective = file.path(TABLE_ROOT, "selective_prediction"),
  dca = file.path(TABLE_ROOT, "decision_curve"),
  bootstrap = file.path(TABLE_ROOT, "bootstrap"),
  manifests = file.path(TABLE_ROOT, "manifests"),
  r_objects = R_OBJECT_ROOT,
  performance_figs = file.path(FIGURE_STEM_ROOT, "performance"),
  calibration_figs = file.path(FIGURE_STEM_ROOT, "calibration"),
  uncertainty_figs = file.path(FIGURE_STEM_ROOT, "uncertainty"),
  selective_figs = file.path(FIGURE_STEM_ROOT, "selective_prediction"),
  dca_figs = file.path(FIGURE_STEM_ROOT, "decision_curve"),
  cohort_shift_figs = file.path(FIGURE_STEM_ROOT, "cohort_shift"),
  bootstrap_figs = file.path(FIGURE_STEM_ROOT, "bootstrap")
)

table_directories <- c(
  dirs$dataset,
  dirs$preprocessing,
  dirs$cohort_shift,
  dirs$metrics,
  dirs$predictions,
  dirs$status,
  dirs$calibration,
  dirs$selective,
  dirs$dca,
  dirs$bootstrap,
  dirs$manifests,
  dirs$r_objects
)

for (directory in unique(c(
  RESULT_ROOT,
  TABLE_ROOT,
  table_directories,
  FIGURE_PDF_ROOT,
  FIGURE_PNG_ROOT
))) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------

# 4. General helper functions
# ------------------------------------------------

theme_b_risa <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 3
      ),
      plot.subtitle = element_text(size = base_size),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

save_plot_pdf_png <- function(
  plot,
  path_without_extension,
  width = 10,
  height = 7,
  dpi = 600
) {
  expected_prefix <- paste0(
    FIGURE_STEM_ROOT,
    .Platform$file.sep
  )

  if (!startsWith(path_without_extension, expected_prefix)) {
    stop(
      "Figure path must begin with FIGURE_STEM_ROOT: ",
      path_without_extension
    )
  }

  relative_path <- substring(
    path_without_extension,
    nchar(expected_prefix) + 1L
  )

  pdf_path <- file.path(
    FIGURE_PDF_ROOT,
    paste0(relative_path, ".pdf")
  )

  png_path <- file.path(
    FIGURE_PNG_ROOT,
    paste0(relative_path, ".png")
  )

  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    pdf_path,
    plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )

  ggplot2::ggsave(
    png_path,
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  print(plot)
}

safe_fwrite <- function(x, path) {
  data.table::fwrite(
    as.data.frame(x),
    path,
    na = "NA"
  )
}

# This function deliberately preserves dimensions and dimnames.
# A direct as.numeric(p) conversion would flatten posterior
# prediction matrices and caused one common probability to be recycled
# across all test patients.
clip_prob <- function(p, eps = EPS) {
  original_dim <- dim(p)
  original_dimnames <- dimnames(p)

  out <- pmin(
    pmax(p, eps),
    1 - eps
  )

  storage.mode(out) <- "double"

  if (!is.null(original_dim)) {
    dim(out) <- original_dim
    dimnames(out) <- original_dimnames
  }

  out
}

validate_prediction_draws <- function(
  draws,
  expected_draws,
  expected_test_size,
  model_name
) {
  draws <- as.matrix(draws)

  observed_dimension <- as.integer(dim(draws))
  expected_dimension <- c(
    as.integer(expected_draws),
    as.integer(expected_test_size)
  )

  if (
    length(observed_dimension) != 2L ||
    !all(observed_dimension == expected_dimension)
  ) {
    stop(
      model_name,
      " returned posterior prediction dimensions ",
      paste(observed_dimension, collapse = " x "),
      "; expected ",
      paste(expected_dimension, collapse = " x "),
      "."
    )
  }

  if (any(!is.finite(draws))) {
    stop(model_name, " returned non-finite posterior predictions.")
  }

  if (any(draws < 0 | draws > 1)) {
    stop(
      model_name,
      " returned values outside the probability interval [0, 1]."
    )
  }

  draws
}

summarise_draws <- function(draws) {
  draws <- as.matrix(draws)

  tibble::tibble(
    p_mean = colMeans(draws),
    p_sd = apply(draws, 2L, stats::sd),
    p_lwr = apply(
      draws,
      2L,
      stats::quantile,
      probs = 0.025,
      names = FALSE,
      type = 8
    ),
    p_upr = apply(
      draws,
      2L,
      stats::quantile,
      probs = 0.975,
      names = FALSE,
      type = 8
    )
  )
}

safe_auc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)

  valid <- is.finite(y) & is.finite(p)
  y <- y[valid]
  p <- p[valid]

  if (length(unique(y)) < 2L) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = y,
          predictor = p,
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        )
      )
    ),
    error = function(e) NA_real_
  )
}

safe_auprc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)

  valid <- is.finite(y) & is.finite(p)
  y <- y[valid]
  p <- p[valid]

  if (length(unique(y)) < 2L) {
    return(NA_real_)
  }

  positive_scores <- p[y == 1L]
  negative_scores <- p[y == 0L]

  if (
    length(positive_scores) == 0L ||
    length(negative_scores) == 0L
  ) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(
      PRROC::pr.curve(
        scores.class0 = positive_scores,
        scores.class1 = negative_scores,
        curve = FALSE
      )$auc.integral
    ),
    error = function(e) NA_real_
  )
}

calc_ece <- function(y, p, n_bins = CALIBRATION_BINS) {
  y <- as.integer(y)
  p <- clip_prob(as.numeric(p))

  bins <- cut(
    p,
    breaks = seq(0, 1, length.out = n_bins + 1L),
    include.lowest = TRUE,
    labels = FALSE
  )

  total_n <- length(y)
  ece <- 0

  for (bin_id in seq_len(n_bins)) {
    index <- which(bins == bin_id)

    if (length(index) == 0L) {
      next
    }

    ece <- ece +
      length(index) / total_n *
      abs(
        mean(y[index]) -
          mean(p[index])
      )
  }

  ece
}

calc_calibration_intercept_slope <- function(y, p) {
  y <- as.integer(y)
  p <- clip_prob(as.numeric(p))

  if (length(unique(y)) < 2L) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  eta <- qlogis(p)

  fit <- tryCatch(
    suppressWarnings(
      stats::glm(
        y ~ eta,
        family = stats::binomial()
      )
    ),
    error = function(e) NULL
  )

  if (
    is.null(fit) ||
    any(!is.finite(stats::coef(fit)))
  ) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  c(
    intercept = unname(stats::coef(fit)[1L]),
    slope = unname(stats::coef(fit)[2L])
  )
}

select_least_uncertain <- function(
  uncertainty,
  target_coverage
) {
  uncertainty <- as.numeric(uncertainty)

  if (
    length(uncertainty) == 0L ||
    any(!is.finite(uncertainty))
  ) {
    return(integer(0))
  }

  number_to_retain <- max(
    1L,
    min(
      length(uncertainty),
      as.integer(ceiling(
        target_coverage * length(uncertainty)
      ))
    )
  )

  order(uncertainty)[seq_len(number_to_retain)]
}

calc_metrics <- function(
  y,
  p_mean,
  p_sd,
  p_lwr,
  p_upr,
  threshold = CLASSIFICATION_THRESHOLD,
  coverage = SELECTIVE_COVERAGE
) {
  y <- as.integer(y)
  p_mean <- clip_prob(as.numeric(p_mean))
  p_sd <- as.numeric(p_sd)
  p_lwr <- as.numeric(p_lwr)
  p_upr <- as.numeric(p_upr)

  valid <- (
    is.finite(y) &
      is.finite(p_mean) &
      is.finite(p_sd) &
      is.finite(p_lwr) &
      is.finite(p_upr)
  )

  y <- y[valid]
  p_mean <- p_mean[valid]
  p_sd <- p_sd[valid]
  p_lwr <- p_lwr[valid]
  p_upr <- p_upr[valid]

  if (length(y) == 0L) {
    stop("No finite observations were available for metric calculation.")
  }

  yhat <- as.integer(p_mean >= threshold)

  tp <- sum(yhat == 1L & y == 1L)
  tn <- sum(yhat == 0L & y == 0L)
  fp <- sum(yhat == 1L & y == 0L)
  fn <- sum(yhat == 0L & y == 1L)

  sensitivity <- if (tp + fn > 0L) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  specificity <- if (tn + fp > 0L) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  ppv <- if (tp + fp > 0L) {
    tp / (tp + fp)
  } else {
    NA_real_
  }

  npv <- if (tn + fn > 0L) {
    tn / (tn + fn)
  } else {
    NA_real_
  }

  f1 <- if (
    is.finite(ppv) &&
    is.finite(sensitivity) &&
    ppv + sensitivity > 0
  ) {
    2 * ppv * sensitivity / (ppv + sensitivity)
  } else {
    NA_real_
  }

  retained_index <- select_least_uncertain(
    p_sd,
    coverage
  )

  selective_risk <- if (length(retained_index) > 0L) {
    mean(
      yhat[retained_index] != y[retained_index]
    )
  } else {
    NA_real_
  }

  calibration <- calc_calibration_intercept_slope(
    y,
    p_mean
  )

  tibble::tibble(
    n_test = length(y),
    prevalence_test = mean(y),
    threshold = threshold,
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn,
    accuracy = mean(yhat == y),
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    f1 = f1,
    auc = safe_auc(y, p_mean),
    auprc = safe_auprc(y, p_mean),
    brier = mean((p_mean - y)^2),
    logloss = -mean(
      y * log(p_mean) +
        (1 - y) * log(1 - p_mean)
    ),
    calibration_intercept = unname(
      calibration["intercept"]
    ),
    calibration_slope = unname(
      calibration["slope"]
    ),
    ece = calc_ece(y, p_mean),
    mean_uncertainty = mean(p_sd),
    mean_interval_width = mean(p_upr - p_lwr),
    selective_coverage = (
      length(retained_index) / length(y)
    ),
    selective_risk = selective_risk
  )
}

# ------------------------------------------------

# ------------------------------------------------
# 4A. External-validation helper functions
# ------------------------------------------------

print_table <- function(x, title = NULL) {
  if (!is.null(title)) {
    cat(
      "\n",
      paste(rep("=", 110), collapse = ""),
      "\n",
      title,
      "\n",
      paste(rep("=", 110), collapse = ""),
      "\n",
      sep = ""
    )
  }

  if (inherits(x, "tbl_df") || inherits(x, "tbl")) {
    print(x, n = Inf, width = Inf)
  } else {
    print(x)
  }
}

save_and_print_table <- function(
  x,
  directory,
  filename,
  title = NULL,
  print_table_object = TRUE
) {
  safe_fwrite(
    x,
    file.path(directory, paste0(filename, ".csv"))
  )

  if (print_table_object) {
    print_table(x, title)
  }

  invisible(x)
}

metric_direction <- function(metric_name) {
  if (metric_name %in% HIGHER_IS_BETTER) {
    return("Higher is better")
  }

  if (metric_name %in% LOWER_IS_BETTER) {
    return("Lower is better")
  }

  "Descriptive; magnitude is model-specific"
}

favored_model <- function(
  metric_name,
  model_1,
  model_2,
  difference
) {
  if (!is.finite(difference)) {
    return(NA_character_)
  }

  if (difference == 0) {
    return("Tie")
  }

  if (metric_name %in% HIGHER_IS_BETTER) {
    return(if (difference > 0) model_1 else model_2)
  }

  if (metric_name %in% LOWER_IS_BETTER) {
    return(if (difference < 0) model_1 else model_2)
  }

  "Interpret descriptively"
}

metric_label <- function(metric_name) {
  labels <- c(
    accuracy = "Accuracy",
    sensitivity = "Sensitivity",
    specificity = "Specificity",
    ppv = "PPV",
    npv = "NPV",
    f1 = "F1 score",
    auc = "AUROC",
    auprc = "AUPRC",
    brier = "Brier score",
    logloss = "Log loss",
    calibration_intercept = "Calibration intercept",
    calibration_slope = "Calibration slope",
    ece = "ECE",
    mean_uncertainty = "Mean predictive SD",
    mean_interval_width = "Mean 95% interval width",
    selective_coverage = "Selective coverage",
    selective_risk = "Selective risk"
  )

  unname(labels[metric_name])
}

read_rpa_scaling_parameters <- function(path) {
  fallback <- RPA_SCALING_FALLBACK

  if (!file.exists(path)) {
    warning(
      "RPA scaling parameter file was not found. ",
      "Using the validated embedded parameters."
    )
    return(fallback)
  }

  candidate <- data.table::fread(
    path,
    data.table = FALSE
  )

  required <- c(
    "Feature",
    "mean_overall",
    "sd_overall"
  )

  missing_columns <- setdiff(required, names(candidate))

  if (length(missing_columns) > 0L) {
    stop(
      "RPA scaling parameter file is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  candidate <- candidate |>
    dplyr::filter(Feature %in% RPA_STD_FEATURES) |>
    dplyr::arrange(match(Feature, RPA_STD_FEATURES)) |>
    dplyr::select(
      Feature,
      mean_overall,
      sd_overall
    )

  if (nrow(candidate) != length(RPA_STD_FEATURES)) {
    stop("The RPA scaling parameter file does not contain all three predictors.")
  }

  difference_mean <- max(abs(
    candidate$mean_overall -
      fallback$mean_overall
  ))

  difference_sd <- max(abs(
    candidate$sd_overall -
      fallback$sd_overall
  ))

  if (
    !is.finite(difference_mean) ||
    !is.finite(difference_sd) ||
    difference_mean > 1e-4 ||
    difference_sd > 1e-4
  ) {
    stop(
      "Saved RPA scaling parameters differ from the validated parameters. ",
      "Do not continue until the standardized RPA source is reconciled."
    )
  }

  candidate
}

load_dusseldorf_dataset <- function(path) {
  if (!file.exists(path)) {
    stop("Düsseldorf dataset file does not exist:\n", path)
  }

  raw <- readr::read_csv(
    path,
    show_col_types = FALSE,
    trim_ws = TRUE,
    na = c(
      "",
      "NA",
      "N/A",
      "na",
      "n/a",
      "NaN",
      "nan",
      "NAN",
      "NULL",
      "null",
      ".",
      "-999",
      "?"
    )
  ) |>
    janitor::clean_names()

  required_columns <- c(
    "diagnosis_by_criteria",
    "body_temp",
    "wbc_count",
    "neutrophil_perc"
  )

  missing_columns <- setdiff(
    required_columns,
    names(raw)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Düsseldorf dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  cleaned <- raw |>
    dplyr::mutate(
      diagnosis_original = trimws(
        tolower(
          as.character(diagnosis_by_criteria)
        )
      ),
      Y = dplyr::case_when(
        diagnosis_original %in% c(
          "appendicitis",
          "acute appendicitis",
          "positive",
          "yes"
        ) ~ 1L,
        diagnosis_original %in% c(
          "no appendicitis",
          "no_appendicitis",
          "non-appendicitis",
          "non appendicitis",
          "nonappendicitis",
          "noappendicitis",
          "negative",
          "no"
        ) ~ 0L,
        TRUE ~ NA_integer_
      ),
      Body_Temperature_Raw = suppressWarnings(
        as.numeric(body_temp)
      ),
      WBC_Count_Original = suppressWarnings(
        as.numeric(wbc_count)
      ),
      WBC_Count_Raw = dplyr::case_when(
        is.na(WBC_Count_Original) ~ NA_real_,
        WBC_Count_Original > 1000 ~
          WBC_Count_Original / 1000,
        TRUE ~ WBC_Count_Original
      ),
      Neutrophil_Percentage_Raw = suppressWarnings(
        as.numeric(neutrophil_perc)
      )
    )

  if (any(is.na(cleaned$Y))) {
    unmapped <- cleaned |>
      dplyr::filter(is.na(Y)) |>
      dplyr::count(diagnosis_original, name = "n")

    print_table(
      unmapped,
      "Unmapped Düsseldorf outcome values"
    )

    stop("At least one Düsseldorf diagnosis could not be mapped.")
  }

  if (dplyr::n_distinct(cleaned$Y) != 2L) {
    stop("The Düsseldorf outcome must contain both diagnostic classes.")
  }

  cleaned <- cleaned |>
    dplyr::mutate(
      Row_ID = dplyr::row_number(),
      Outcome = factor(
        Y,
        levels = c(0, 1),
        labels = c(
          "No_appendicitis",
          "Appendicitis"
        )
      ),
      Complete_Three_Predictors = stats::complete.cases(
        Body_Temperature_Raw,
        WBC_Count_Raw,
        Neutrophil_Percentage_Raw
      ),
      Missing_Predictor_Count = rowSums(
        cbind(
          is.na(Body_Temperature_Raw),
          is.na(WBC_Count_Raw),
          is.na(Neutrophil_Percentage_Raw)
        )
      )
    )

  cleaned
}

convert_dusseldorf_to_rpa_upstream_scale <- function(
  dusseldorf_data,
  scaling_parameters
) {
  output <- tibble::tibble(
    Y = dusseldorf_data$Y,
    Body_Temperature = dusseldorf_data$Body_Temperature_Raw,
    WBC_Count = dusseldorf_data$WBC_Count_Raw,
    Neutrophil_Percentage = dusseldorf_data$Neutrophil_Percentage_Raw,
    Row_ID = dusseldorf_data$Row_ID,
    Outcome = dusseldorf_data$Outcome,
    Complete_Three_Predictors =
      dusseldorf_data$Complete_Three_Predictors,
    Missing_Predictor_Count =
      dusseldorf_data$Missing_Predictor_Count
  )

  for (feature_name in RPA_STD_FEATURES) {
    feature_mean <- scaling_parameters$mean_overall[
      scaling_parameters$Feature == feature_name
    ]

    feature_sd <- scaling_parameters$sd_overall[
      scaling_parameters$Feature == feature_name
    ]

    output[[feature_name]] <- (
      output[[feature_name]] -
        feature_mean
    ) / feature_sd
  }

  output
}

reconstruct_rpa_raw_scale <- function(
  rpa_standardized_data,
  scaling_parameters
) {
  output <- rpa_standardized_data

  for (feature_name in RPA_STD_FEATURES) {
    feature_mean <- scaling_parameters$mean_overall[
      scaling_parameters$Feature == feature_name
    ]

    feature_sd <- scaling_parameters$sd_overall[
      scaling_parameters$Feature == feature_name
    ]

    output[[feature_name]] <-
      output[[feature_name]] * feature_sd +
      feature_mean
  }

  output
}

make_external_fit_seed <- function(model_name) {
  model_code <- match(model_name, MODEL_SET)

  as.integer(
    EXTERNAL_FIT_SEED +
      1000L * model_code
  )
}

calculate_metric_vector <- function(
  y,
  p_mean,
  p_sd,
  p_lwr,
  p_upr
) {
  metric_row <- calc_metrics(
    y = y,
    p_mean = p_mean,
    p_sd = p_sd,
    p_lwr = p_lwr,
    p_upr = p_upr
  )

  values <- unlist(
    metric_row[
      c(
        INFERENTIAL_METRICS,
        DESCRIPTIVE_METRICS
      )
    ],
    use.names = TRUE
  )

  as.numeric_values <- as.numeric(values)
  names(as.numeric_values) <- names(values)
  as.numeric_values
}

summarise_predictor <- function(x) {
  x_finite <- x[is.finite(x)]

  tibble::tibble(
    n_available = length(x_finite),
    n_missing = sum(is.na(x)),
    missing_percent = 100 * mean(is.na(x)),
    mean = if (length(x_finite) > 0L) mean(x_finite) else NA_real_,
    sd = if (length(x_finite) > 1L) stats::sd(x_finite) else NA_real_,
    median = if (length(x_finite) > 0L) stats::median(x_finite) else NA_real_,
    q1 = if (length(x_finite) > 0L) as.numeric(stats::quantile(x_finite, 0.25, names = FALSE)) else NA_real_,
    q3 = if (length(x_finite) > 0L) as.numeric(stats::quantile(x_finite, 0.75, names = FALSE)) else NA_real_,
    minimum = if (length(x_finite) > 0L) min(x_finite) else NA_real_,
    maximum = if (length(x_finite) > 0L) max(x_finite) else NA_real_
  )
}

calculate_cohort_shift <- function(
  rpa_raw,
  dusseldorf_raw
) {
  feature_map <- tibble::tibble(
    Feature = RPA_STD_FEATURES,
    Dusseldorf_column = c(
      "Body_Temperature_Raw",
      "WBC_Count_Raw",
      "Neutrophil_Percentage_Raw"
    )
  )

  summary_rows <- list()
  shift_rows <- list()

  for (row_id in seq_len(nrow(feature_map))) {
    feature_name <- feature_map$Feature[row_id]
    dusseldorf_column <- feature_map$Dusseldorf_column[row_id]

    rpa_values <- rpa_raw[[feature_name]]
    dusseldorf_values <- dusseldorf_raw[[dusseldorf_column]]

    summary_rows[[length(summary_rows) + 1L]] <-
      summarise_predictor(rpa_values) |>
      dplyr::mutate(
        Cohort = "RPA development",
        Feature = feature_name,
        .before = 1
      )

    summary_rows[[length(summary_rows) + 1L]] <-
      summarise_predictor(dusseldorf_values) |>
      dplyr::mutate(
        Cohort = "Düsseldorf external",
        Feature = feature_name,
        .before = 1
      )

    rpa_finite <- rpa_values[is.finite(rpa_values)]
    dusseldorf_finite <- dusseldorf_values[is.finite(dusseldorf_values)]

    pooled_sd <- sqrt(
      (
        (length(rpa_finite) - 1L) * stats::var(rpa_finite) +
          (length(dusseldorf_finite) - 1L) * stats::var(dusseldorf_finite)
      ) /
        (
          length(rpa_finite) +
            length(dusseldorf_finite) -
            2L
        )
    )

    smd <- if (
      is.finite(pooled_sd) &&
      pooled_sd > 0
    ) {
      (
        mean(dusseldorf_finite) -
          mean(rpa_finite)
      ) / pooled_sd
    } else {
      NA_real_
    }

    shift_rows[[length(shift_rows) + 1L]] <-
      tibble::tibble(
        Feature = feature_name,
        RPA_mean = mean(rpa_finite),
        Dusseldorf_mean = mean(dusseldorf_finite),
        Mean_difference_Dusseldorf_minus_RPA =
          mean(dusseldorf_finite) -
          mean(rpa_finite),
        Pooled_SD = pooled_sd,
        Standardized_mean_difference = smd,
        Absolute_SMD = abs(smd),
        RPA_missing_percent = 100 * mean(is.na(rpa_values)),
        Dusseldorf_missing_percent = 100 * mean(is.na(dusseldorf_values)),
        Missingness_difference_percentage_points =
          100 * mean(is.na(dusseldorf_values)) -
          100 * mean(is.na(rpa_values))
      )
  }

  list(
    summary = dplyr::bind_rows(summary_rows),
    shift = dplyr::bind_rows(shift_rows)
  )
}

make_confusion_table <- function(prediction_data) {
  prediction_data |>
    dplyr::mutate(
      Predicted = factor(
        as.integer(p_mean >= CLASSIFICATION_THRESHOLD),
        levels = c(0, 1),
        labels = c(
          "No_appendicitis",
          "Appendicitis"
        )
      )
    ) |>
    dplyr::count(
      Analysis_Set,
      Analysis_Label,
      Model,
      Outcome_true,
      Predicted,
      name = "n",
      .drop = FALSE
    )
}

make_calibration_bins <- function(prediction_data) {
  prediction_data |>
    dplyr::mutate(
      Bin = cut(
        p_mean,
        breaks = seq(
          0,
          1,
          length.out = CALIBRATION_BINS + 1L
        ),
        include.lowest = TRUE,
        labels = FALSE
      )
    ) |>
    dplyr::group_by(
      Analysis_Set,
      Analysis_Label,
      Model,
      Bin
    ) |>
    dplyr::summarise(
      Mean_Predicted = mean(p_mean),
      Observed_Rate = mean(Y_true),
      n = dplyr::n(),
      .groups = "drop"
    )
}

make_risk_coverage_data <- function(prediction_data) {
  coverage_grid <- seq(0.50, 1.00, by = 0.05)
  rows <- list()

  for (analysis_set in unique(prediction_data$Analysis_Set)) {
    analysis_data <- prediction_data |>
      dplyr::filter(Analysis_Set == analysis_set)

    analysis_label <- unique(analysis_data$Analysis_Label)

    for (model_name in MODEL_ORDER) {
      model_data <- analysis_data |>
        dplyr::filter(Model == model_name) |>
        dplyr::arrange(Row_ID)

      predicted_class <- as.integer(
        model_data$p_mean >= CLASSIFICATION_THRESHOLD
      )

      for (target_coverage in coverage_grid) {
        retained_index <- select_least_uncertain(
          model_data$p_sd,
          target_coverage
        )

        rows[[length(rows) + 1L]] <- tibble::tibble(
          Analysis_Set = analysis_set,
          Analysis_Label = analysis_label,
          Model = model_name,
          Target_Coverage = target_coverage,
          Coverage = length(retained_index) / nrow(model_data),
          Retained_n = length(retained_index),
          Selective_Risk = mean(
            predicted_class[retained_index] !=
              model_data$Y_true[retained_index]
          )
        )
      }
    }
  }

  dplyr::bind_rows(rows)
}

make_decision_curve_data <- function(prediction_data) {
  rows <- list()

  for (analysis_set in unique(prediction_data$Analysis_Set)) {
    analysis_data <- prediction_data |>
      dplyr::filter(Analysis_Set == analysis_set)

    analysis_label <- unique(analysis_data$Analysis_Label)

    dca <- decision_curve_table(
      analysis_data,
      DCA_THRESHOLDS
    ) |>
      dplyr::mutate(
        Analysis_Set = analysis_set,
        Analysis_Label = analysis_label,
        .before = 1
      )

    rows[[length(rows) + 1L]] <- dca
  }

  dplyr::bind_rows(rows) |>
    dplyr::mutate(
      Standardized_Net_Benefit =
        Net_Benefit / Prevalence
    )
}

bootstrap_external_predictions <- function(
  prediction_data,
  bootstrap_reps = BOOTSTRAP_REPS,
  bootstrap_seed = BOOTSTRAP_SEED
) {
  analysis_sets <- unique(prediction_data$Analysis_Set)
  model_estimate_rows <- list()
  pairwise_rows <- list()
  draw_rows <- list()

  set.seed(bootstrap_seed)

  for (analysis_set in analysis_sets) {
    analysis_data <- prediction_data |>
      dplyr::filter(Analysis_Set == analysis_set)

    analysis_label <- unique(analysis_data$Analysis_Label)

    prediction_by_model <- lapply(
      MODEL_ORDER,
      function(model_name) {
        analysis_data |>
          dplyr::filter(Model == model_name) |>
          dplyr::arrange(Row_ID)
      }
    )

    names(prediction_by_model) <- MODEL_ORDER

    reference_ids <- prediction_by_model[[1L]]$Row_ID
    reference_y <- prediction_by_model[[1L]]$Y_true

    for (model_name in MODEL_ORDER) {
      model_data <- prediction_by_model[[model_name]]

      if (!identical(model_data$Row_ID, reference_ids)) {
        stop(
          "Patient ordering differs across models in analysis set ",
          analysis_set,
          "."
        )
      }

      if (!identical(
        as.integer(model_data$Y_true),
        as.integer(reference_y)
      )) {
        stop(
          "Outcome ordering differs across models in analysis set ",
          analysis_set,
          "."
        )
      }
    }

    point_rows <- list()

    for (model_name in MODEL_ORDER) {
      model_data <- prediction_by_model[[model_name]]

      values <- calculate_metric_vector(
        y = model_data$Y_true,
        p_mean = model_data$p_mean,
        p_sd = model_data$p_sd,
        p_lwr = model_data$p_lwr,
        p_upr = model_data$p_upr
      )

      point_rows[[length(point_rows) + 1L]] <-
        tibble::tibble(
          Analysis_Set = analysis_set,
          Analysis_Label = analysis_label,
          Model = model_name,
          Metric = names(values),
          Point_estimate = as.numeric(values)
        )
    }

    point_estimates <- dplyr::bind_rows(point_rows) |>
      dplyr::filter(Metric %in% BOOTSTRAP_MODEL_METRICS)

    analysis_bootstrap_rows <- vector(
      mode = "list",
      length = bootstrap_reps * length(MODEL_ORDER)
    )

    row_index <- 0L
    patient_count <- length(reference_ids)

    for (bootstrap_id in seq_len(bootstrap_reps)) {
      sampled_index <- sample.int(
        patient_count,
        size = patient_count,
        replace = TRUE
      )

      for (model_name in MODEL_ORDER) {
        row_index <- row_index + 1L
        model_data <- prediction_by_model[[model_name]]

        values <- calculate_metric_vector(
          y = model_data$Y_true[sampled_index],
          p_mean = model_data$p_mean[sampled_index],
          p_sd = model_data$p_sd[sampled_index],
          p_lwr = model_data$p_lwr[sampled_index],
          p_upr = model_data$p_upr[sampled_index]
        )

        analysis_bootstrap_rows[[row_index]] <-
          tibble::as_tibble_row(
            c(
              list(
                Analysis_Set = analysis_set,
                Analysis_Label = analysis_label,
                Bootstrap = bootstrap_id,
                Model = model_name
              ),
              as.list(values)
            )
          )
      }

      if (bootstrap_id %% 100L == 0L) {
        cat(
          sprintf(
            "\r%s patient bootstrap %d/%d",
            analysis_set,
            bootstrap_id,
            bootstrap_reps
          )
        )
      }
    }

    cat("\n")

    analysis_draws <- dplyr::bind_rows(
      analysis_bootstrap_rows
    )

    analysis_draws_long <- analysis_draws |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(
          BOOTSTRAP_MODEL_METRICS
        ),
        names_to = "Metric",
        values_to = "Value"
      )

    draw_rows[[length(draw_rows) + 1L]] <-
      analysis_draws_long

    model_estimates <- analysis_draws_long |>
      dplyr::filter(is.finite(Value)) |>
      dplyr::group_by(
        Analysis_Set,
        Analysis_Label,
        Model,
        Metric
      ) |>
      dplyr::summarise(
        Bootstrap_finite_replicates = dplyr::n(),
        Bootstrap_mean = mean(Value),
        Bootstrap_standard_error = stats::sd(Value),
        CI_lower = as.numeric(
          stats::quantile(
            Value,
            probs = (1 - CONFIDENCE_LEVEL) / 2,
            names = FALSE
          )
        ),
        CI_upper = as.numeric(
          stats::quantile(
            Value,
            probs = 1 - (1 - CONFIDENCE_LEVEL) / 2,
            names = FALSE
          )
        ),
        .groups = "drop"
      ) |>
      dplyr::left_join(
        point_estimates,
        by = c(
          "Analysis_Set",
          "Analysis_Label",
          "Model",
          "Metric"
        )
      ) |>
      dplyr::mutate(
        Metric_Label = vapply(
          Metric,
          metric_label,
          character(1)
        ),
        Direction = vapply(
          Metric,
          metric_direction,
          character(1)
        )
      )

    model_estimate_rows[[length(model_estimate_rows) + 1L]] <-
      model_estimates

    model_pairs <- combn(
      MODEL_ORDER,
      2L,
      simplify = FALSE
    )

    for (metric_name in BOOTSTRAP_PAIRWISE_METRICS) {
      metric_draws <- analysis_draws_long |>
        dplyr::filter(Metric == metric_name) |>
        dplyr::select(
          Bootstrap,
          Model,
          Value
        ) |>
        tidyr::pivot_wider(
          names_from = Model,
          values_from = Value
        )

      metric_points <- point_estimates |>
        dplyr::filter(Metric == metric_name)

      metric_pair_rows <- list()

      for (model_pair in model_pairs) {
        model_1 <- model_pair[1L]
        model_2 <- model_pair[2L]

        difference <-
          metric_draws[[model_1]] -
          metric_draws[[model_2]]

        difference <- difference[
          is.finite(difference)
        ]

        point_1 <- metric_points |>
          dplyr::filter(Model == model_1) |>
          dplyr::pull(Point_estimate)

        point_2 <- metric_points |>
          dplyr::filter(Model == model_2) |>
          dplyr::pull(Point_estimate)

        point_difference <- point_1 - point_2

        if (length(difference) == 0L) {
          ci_lower <- NA_real_
          ci_upper <- NA_real_
          p_value_raw <- NA_real_
        } else {
          ci_lower <- as.numeric(
            stats::quantile(
              difference,
              probs = (1 - CONFIDENCE_LEVEL) / 2,
              names = FALSE
            )
          )

          ci_upper <- as.numeric(
            stats::quantile(
              difference,
              probs = 1 - (1 - CONFIDENCE_LEVEL) / 2,
              names = FALSE
            )
          )

          p_value_raw <- min(
            1,
            2 * min(
              mean(difference <= 0),
              mean(difference >= 0)
            )
          )
        }

        metric_pair_rows[[length(metric_pair_rows) + 1L]] <-
          tibble::tibble(
            Analysis_Set = analysis_set,
            Analysis_Label = analysis_label,
            Metric = metric_name,
            Metric_Label = metric_label(metric_name),
            Comparison = paste0(model_1, " - ", model_2),
            Model_1 = model_1,
            Model_2 = model_2,
            Model_1_point_estimate = point_1,
            Model_2_point_estimate = point_2,
            Point_difference = point_difference,
            Bootstrap_finite_replicates = length(difference),
            CI_lower = ci_lower,
            CI_upper = ci_upper,
            p_value_raw = p_value_raw,
            Direction = metric_direction(metric_name),
            Favored_model = favored_model(
              metric_name,
              model_1,
              model_2,
              point_difference
            )
          )
      }

      metric_pair_data <- dplyr::bind_rows(
        metric_pair_rows
      ) |>
        dplyr::mutate(
          p_value_adjusted = stats::p.adjust(
            p_value_raw,
            method = P_ADJUST_METHOD
          ),
          CI_excludes_zero =
            is.finite(CI_lower) &
            is.finite(CI_upper) &
            (
              CI_lower > 0 |
              CI_upper < 0
            ),
          Significant_adjusted =
            is.finite(p_value_adjusted) &
            p_value_adjusted < ALPHA,
          Interpretation = paste0(
            "Paired patient bootstrap with fixed locked ",
            "RPA-trained predictions."
          )
        )

      pairwise_rows[[length(pairwise_rows) + 1L]] <-
        metric_pair_data
    }
  }

  list(
    model_estimates = dplyr::bind_rows(model_estimate_rows),
    pairwise_results = dplyr::bind_rows(pairwise_rows),
    bootstrap_draws_long = dplyr::bind_rows(draw_rows)
  )
}


# 5. Data loading and preprocessing
# ------------------------------------------------

load_rpa_standardized_dataset <- function(path) {
  if (!file.exists(path)) {
    stop("Standardized RPA dataset file does not exist:\n", path)
  }

  data <- data.table::fread(
    path,
    data.table = FALSE
  )

  required_columns <- c("Y", RPA_STD_FEATURES)
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "The standardized RPA dataset is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  data <- data[, required_columns, drop = FALSE]
  data$Y <- as.integer(data$Y)

  data <- data[
    !is.na(data$Y) &
      data$Y %in% c(0L, 1L),
    ,
    drop = FALSE
  ]

  if (length(unique(data$Y)) != 2L) {
    stop("The outcome Y must contain both 0 and 1.")
  }

  for (feature_name in RPA_STD_FEATURES) {
    data[[feature_name]] <- suppressWarnings(
      as.numeric(data[[feature_name]])
    )
  }

  data$Row_ID <- seq_len(nrow(data))
  data$Outcome <- factor(
    data$Y,
    levels = c(0, 1),
    labels = c(
      "No_appendicitis",
      "Appendicitis"
    )
  )

  data
}

fit_preprocess <- function(train_data) {
  medians <- vapply(
    train_data[, RPA_STD_FEATURES, drop = FALSE],
    function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }

      stats::median(x, na.rm = TRUE)
    },
    numeric(1)
  )

  if (any(!is.finite(medians))) {
    stop(
      "At least one predictor is entirely missing in a training fold."
    )
  }

  imputed <- train_data[
    ,
    RPA_STD_FEATURES,
    drop = FALSE
  ]

  for (feature_name in RPA_STD_FEATURES) {
    missing_index <- is.na(imputed[[feature_name]])

    if (any(missing_index)) {
      imputed[[feature_name]][missing_index] <-
        medians[[feature_name]]
    }
  }

  centers <- vapply(
    imputed,
    mean,
    numeric(1)
  )

  scales <- vapply(
    imputed,
    stats::sd,
    numeric(1)
  )

  scales[
    !is.finite(scales) |
      scales <= 0
  ] <- 1

  list(
    medians = medians,
    centers = centers,
    scales = scales
  )
}

apply_preprocess <- function(data, preprocessing) {
  x <- data[
    ,
    RPA_STD_FEATURES,
    drop = FALSE
  ]

  for (feature_name in RPA_STD_FEATURES) {
    missing_index <- is.na(x[[feature_name]])

    if (any(missing_index)) {
      x[[feature_name]][missing_index] <-
        preprocessing$medians[[feature_name]]
    }

    x[[feature_name]] <- (
      x[[feature_name]] -
        preprocessing$centers[[feature_name]]
    ) / preprocessing$scales[[feature_name]]
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  colnames(x) <- RPA_STD_FEATURES

  if (any(!is.finite(x))) {
    stop("Non-finite values remain after preprocessing.")
  }

  x
}

# ------------------------------------------------
# 7. BLR
# ------------------------------------------------

softplus_stable <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

nearest_positive_definite <- function(
  covariance_matrix,
  eigen_floor = BLR_HESSIAN_EIGEN_FLOOR
) {
  covariance_matrix <- as.matrix(covariance_matrix)
  covariance_matrix <- (
    covariance_matrix +
      t(covariance_matrix)
  ) / 2

  decomposition <- eigen(
    covariance_matrix,
    symmetric = TRUE
  )

  decomposition$values[
    decomposition$values < eigen_floor
  ] <- eigen_floor

  repaired <- decomposition$vectors %*%
    diag(
      decomposition$values,
      nrow = length(decomposition$values)
    ) %*%
    t(decomposition$vectors)

  (
    repaired +
      t(repaired)
  ) / 2
}

predict_blr <- function(
  x_train,
  y_train,
  x_test,
  n_draws
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  if (length(unique(y_train)) < 2L) {
    stop("BLR training data contain only one outcome class.")
  }

  X_train <- cbind(
    "(Intercept)" = 1,
    x_train
  )

  X_test <- cbind(
    "(Intercept)" = 1,
    x_test
  )

  parameter_count <- ncol(X_train)

  smoothed_prevalence <- (
    sum(y_train) + 0.5
  ) / (
    length(y_train) + 1
  )

  prior_mean <- c(
    qlogis(smoothed_prevalence),
    rep(0, parameter_count - 1L)
  )

  prior_sd <- c(
    BLR_INTERCEPT_PRIOR_SD,
    rep(
      BLR_COEFFICIENT_PRIOR_SD,
      parameter_count - 1L
    )
  )

  prior_precision <- 1 / prior_sd^2

  negative_log_posterior <- function(beta) {
    eta <- as.vector(X_train %*% beta)

    negative_log_likelihood <- sum(
      softplus_stable(eta) -
        y_train * eta
    )

    negative_log_prior <- 0.5 * sum(
      prior_precision *
        (beta - prior_mean)^2
    )

    negative_log_likelihood +
      negative_log_prior
  }

  gradient_negative_log_posterior <- function(beta) {
    eta <- as.vector(X_train %*% beta)
    probability <- plogis(eta)

    as.vector(
      crossprod(
        X_train,
        probability - y_train
      )
    ) +
      prior_precision *
      (beta - prior_mean)
  }

  optimization <- stats::optim(
    par = prior_mean,
    fn = negative_log_posterior,
    gr = gradient_negative_log_posterior,
    method = "BFGS",
    control = list(
      maxit = BLR_OPTIM_MAXIT,
      reltol = BLR_OPTIM_RELTOL
    ),
    hessian = FALSE
  )

  beta_map <- optimization$par
  gradient_at_map <- gradient_negative_log_posterior(
    beta_map
  )

  maximum_absolute_gradient <- max(
    abs(gradient_at_map)
  )

  converged <- (
    optimization$convergence == 0L ||
      maximum_absolute_gradient <= BLR_MAX_ABS_GRADIENT
  )

  if (!converged) {
    stop(
      "BLR optimization did not converge. Code = ",
      optimization$convergence,
      "; max absolute gradient = ",
      signif(maximum_absolute_gradient, 6),
      "."
    )
  }

  fitted_eta <- as.vector(
    X_train %*% beta_map
  )

  fitted_probability <- plogis(fitted_eta)

  logistic_weights <- pmax(
    fitted_probability *
      (1 - fitted_probability),
    1e-8
  )

  hessian <- crossprod(
    X_train,
    X_train * logistic_weights
  ) +
    diag(
      prior_precision,
      nrow = parameter_count
    )

  hessian <- (
    hessian +
      t(hessian)
  ) / 2

  hessian_eigenvalues <- eigen(
    hessian,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (any(!is.finite(hessian_eigenvalues))) {
    stop("BLR Hessian has non-finite eigenvalues.")
  }

  posterior_covariance <- tryCatch(
    solve(hessian),
    error = function(e) NULL
  )

  if (
    is.null(posterior_covariance) ||
    any(!is.finite(posterior_covariance))
  ) {
    stop("BLR Hessian inversion failed.")
  }

  posterior_covariance <- nearest_positive_definite(
    posterior_covariance
  )

  beta_draws <- MASS::mvrnorm(
    n = n_draws,
    mu = beta_map,
    Sigma = posterior_covariance
  )

  beta_draws <- matrix(
    beta_draws,
    nrow = n_draws,
    ncol = parameter_count
  )

  linear_predictor_draws <- beta_draws %*%
    t(X_test)

  probability_draws <- plogis(
    linear_predictor_draws
  )

  probability_draws <- clip_prob(
    probability_draws
  )

  diagnostics <- list(
    Optimizer_convergence_code =
      as.integer(optimization$convergence),
    Optimizer_objective =
      as.numeric(optimization$value),
    Max_absolute_gradient =
      as.numeric(maximum_absolute_gradient),
    Coefficient_L2_norm =
      sqrt(sum(beta_map[-1L]^2)),
    Intercept_MAP =
      as.numeric(beta_map[1L]),
    Minimum_Hessian_eigenvalue =
      min(hessian_eigenvalues),
    Maximum_Hessian_eigenvalue =
      max(hessian_eigenvalues),
    Hessian_condition_number =
      max(hessian_eigenvalues) /
      min(hessian_eigenvalues),
    Training_prevalence =
      mean(y_train),
    Prior_center_prevalence =
      smoothed_prevalence
  )

  list(
    draws = probability_draws,
    diagnostics = diagnostics
  )
}

# ------------------------------------------------
# 8. Genuine Bayesian LDA
# ------------------------------------------------

predict_blda <- function(
  x_train,
  y_train,
  x_test,
  n_draws = N_POST_DRAWS_BLDA
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.integer(y_train)

  index_0 <- which(y_train == 0L)
  index_1 <- which(y_train == 1L)

  n_0 <- length(index_0)
  n_1 <- length(index_1)
  n_total <- n_0 + n_1
  predictor_count <- ncol(x_train)
  test_count <- nrow(x_test)

  if (n_0 < 2L || n_1 < 2L) {
    stop(
      "BLDA requires at least two training observations per observed class."
    )
  }

  x_0 <- x_train[index_0, , drop = FALSE]
  x_1 <- x_train[index_1, , drop = FALSE]

  xbar_0 <- colMeans(x_0)
  xbar_1 <- colMeans(x_1)

  centered_0 <- sweep(
    x_0,
    2L,
    xbar_0,
    FUN = "-"
  )

  centered_1 <- sweep(
    x_1,
    2L,
    xbar_1,
    FUN = "-"
  )

  scatter_0 <- crossprod(centered_0)
  scatter_1 <- crossprod(centered_1)

  prior_mean_0 <- rep(0, predictor_count)
  prior_mean_1 <- rep(0, predictor_count)

  kappa_prior_0 <- BLDA_MEAN_PRIOR_KAPPA
  kappa_prior_1 <- BLDA_MEAN_PRIOR_KAPPA

  kappa_post_0 <- kappa_prior_0 + n_0
  kappa_post_1 <- kappa_prior_1 + n_1

  posterior_mean_0 <- (
    kappa_prior_0 * prior_mean_0 +
      n_0 * xbar_0
  ) / kappa_post_0

  posterior_mean_1 <- (
    kappa_prior_1 * prior_mean_1 +
      n_1 * xbar_1
  ) / kappa_post_1

  nu_prior <- predictor_count + BLDA_IW_DF_OFFSET

  if (nu_prior <= predictor_count + 1L) {
    stop(
      "BLDA inverse-Wishart prior requires nu0 > p + 1."
    )
  }

  prior_scale <- diag(
    BLDA_PRIOR_COVARIANCE_SCALE,
    nrow = predictor_count,
    ncol = predictor_count
  )

  mean_shift_0 <- xbar_0 - prior_mean_0
  mean_shift_1 <- xbar_1 - prior_mean_1

  posterior_scale <- (
    prior_scale +
      scatter_0 +
      scatter_1 +
      (
        kappa_prior_0 * n_0 /
          kappa_post_0
      ) * tcrossprod(mean_shift_0) +
      (
        kappa_prior_1 * n_1 /
          kappa_post_1
      ) * tcrossprod(mean_shift_1)
  )

  posterior_scale <- nearest_positive_definite(
    posterior_scale,
    eigen_floor = BLDA_PRECISION_EIGEN_FLOOR
  )

  nu_post <- nu_prior + n_total

  posterior_precision_scale <- tryCatch(
    solve(posterior_scale),
    error = function(e) NULL
  )

  if (
    is.null(posterior_precision_scale) ||
      any(!is.finite(posterior_precision_scale))
  ) {
    stop(
      "BLDA posterior inverse-Wishart scale inversion failed."
    )
  }

  posterior_precision_scale <- nearest_positive_definite(
    posterior_precision_scale,
    eigen_floor = BLDA_PRECISION_EIGEN_FLOOR
  )

  # If Sigma | data ~ Inv-Wishart(nu_post, posterior_scale),
  # then Omega = Sigma^{-1} | data ~ Wishart(
  #   nu_post, posterior_scale^{-1}
  # ).
  precision_draw_array <- stats::rWishart(
    n = n_draws,
    df = nu_post,
    Sigma = posterior_precision_scale
  )

  class_1_prior_draws <- stats::rbeta(
    n_draws,
    shape1 = BLDA_CLASS_PRIOR_ALPHA + n_1,
    shape2 = BLDA_CLASS_PRIOR_BETA + n_0
  )

  probability_draws <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = test_count
  )

  sampled_precision_min_eigen <- numeric(n_draws)
  sampled_precision_condition <- numeric(n_draws)

  for (draw_id in seq_len(n_draws)) {
    precision_draw <- precision_draw_array[
      ,
      ,
      draw_id,
      drop = FALSE
    ][, , 1L]

    precision_draw <- (
      precision_draw +
        t(precision_draw)
    ) / 2

    precision_eigenvalues <- eigen(
      precision_draw,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    if (
      any(!is.finite(precision_eigenvalues)) ||
        min(precision_eigenvalues) <= 0
    ) {
      stop(
        "BLDA sampled a non-positive-definite precision matrix at draw ",
        draw_id,
        "."
      )
    }

    sampled_precision_min_eigen[draw_id] <-
      min(precision_eigenvalues)

    sampled_precision_condition[draw_id] <-
      max(precision_eigenvalues) /
      min(precision_eigenvalues)

    # R's chol() returns upper-triangular U such that
    # t(U) %*% U = precision. Therefore solve(U, z)
    # has covariance precision^{-1}.
    precision_chol <- chol(precision_draw)

    mean_0_draw <- posterior_mean_0 +
      as.vector(
        backsolve(
          precision_chol,
          stats::rnorm(predictor_count)
        )
      ) /
      sqrt(kappa_post_0)

    mean_1_draw <- posterior_mean_1 +
      as.vector(
        backsolve(
          precision_chol,
          stats::rnorm(predictor_count)
        )
      ) /
      sqrt(kappa_post_1)

    class_1_prior <- class_1_prior_draws[draw_id]
    class_0_prior <- 1 - class_1_prior

    precision_mean_0 <- as.vector(
      precision_draw %*% mean_0_draw
    )

    precision_mean_1 <- as.vector(
      precision_draw %*% mean_1_draw
    )

    discriminant_0 <- (
      as.vector(
        x_test %*%
          precision_mean_0
      ) -
        0.5 *
        sum(
          mean_0_draw *
            precision_mean_0
        ) +
        log(class_0_prior)
    )

    discriminant_1 <- (
      as.vector(
        x_test %*%
          precision_mean_1
      ) -
        0.5 *
        sum(
          mean_1_draw *
            precision_mean_1
        ) +
        log(class_1_prior)
    )

    probability_draws[draw_id, ] <- plogis(
      discriminant_1 -
        discriminant_0
    )
  }

  probability_draws <- clip_prob(
    probability_draws
  )

  probability_draws <- validate_prediction_draws(
    draws = probability_draws,
    expected_draws = n_draws,
    expected_test_size = nrow(x_test),
    model_name = "BLDA"
  )

  diagnostics <- list(
    BLDA_method =
      "Normal-Inverse-Wishart shared-covariance Bayesian LDA",
    BLDA_prior_kappa =
      BLDA_MEAN_PRIOR_KAPPA,
    BLDA_prior_df =
      nu_prior,
    BLDA_posterior_df =
      nu_post,
    BLDA_prior_covariance_scale =
      BLDA_PRIOR_COVARIANCE_SCALE,
    BLDA_beta_prior_alpha =
      BLDA_CLASS_PRIOR_ALPHA,
    BLDA_beta_prior_beta =
      BLDA_CLASS_PRIOR_BETA,
    BLDA_mean_sampled_precision_min_eigen =
      mean(sampled_precision_min_eigen),
    BLDA_max_sampled_precision_condition =
      max(sampled_precision_condition),
    Training_prevalence =
      mean(y_train),
    Posterior_mean_class_prevalence =
      mean(class_1_prior_draws)
  )

  list(
    draws = probability_draws,
    diagnostics = diagnostics
  )
}

# ------------------------------------------------
# 9. BART
# ------------------------------------------------

predict_bart <- function(
  x_train,
  y_train,
  x_test,
  n_draws
) {
  fit <- BART::pbart(
    x.train = as.matrix(x_train),
    y.train = as.numeric(y_train),
    x.test = as.matrix(x_test),
    ntree = BART_NTREE,
    ndpost = n_draws,
    nskip = N_BURN_BART,
    printevery = 100000L
  )

  if (is.null(fit$yhat.test)) {
    stop("BART did not return yhat.test.")
  }

  latent_draws <- as.matrix(
    fit$yhat.test
  )

  # BART::pbart returns f(x) + binaryOffset on the latent probit scale.
  # Posterior event-probability draws are Phi(f(x) + binaryOffset).
  probability_draws <- stats::pnorm(
    latent_draws
  )

  expected_test_size <- nrow(x_test)

  if (
    nrow(probability_draws) == n_draws &&
    ncol(probability_draws) == expected_test_size
  ) {
    # Documented BART orientation: posterior draws in rows, patients in columns.
    probability_draws <- probability_draws
  } else if (
    nrow(probability_draws) == expected_test_size &&
    ncol(probability_draws) == n_draws
  ) {
    probability_draws <- t(probability_draws)
  } else {
    stop(
      "Unexpected BART prediction dimensions: ",
      paste(dim(probability_draws), collapse = " x "),
      ". Expected ",
      n_draws,
      " x ",
      expected_test_size,
      "."
    )
  }

  probability_draws <- clip_prob(
    probability_draws
  )

  list(
    draws = probability_draws,
    diagnostics = list(
      BART_binary_offset = if (
        length(fit$binaryOffset) == 0L
      ) {
        NA_real_
      } else {
        as.numeric(fit$binaryOffset[1L])
      },
      BART_number_of_trees =
        as.integer(BART_NTREE),
      Training_prevalence =
        mean(y_train)
    )
  )
}

fit_model <- function(
  model_name,
  x_train,
  y_train,
  x_test
) {
  if (model_name == "BLR") {
    return(
      predict_blr(
        x_train,
        y_train,
        x_test,
        N_POST_DRAWS_BLR
      )
    )
  }

  if (model_name == "BLDA") {
    return(
      predict_blda(
        x_train,
        y_train,
        x_test,
        N_POST_DRAWS_BLDA
      )
    )
  }

  if (model_name == "BART") {
    return(
      predict_bart(
        x_train,
        y_train,
        x_test,
        N_POST_DRAWS_BART
      )
    )
  }

  stop("Unknown model: ", model_name)
}

expected_draw_count <- function(model_name) {
  switch(
    model_name,
    BLR = N_POST_DRAWS_BLR,
    BLDA = N_POST_DRAWS_BLDA,
    BART = N_POST_DRAWS_BART,
    stop("Unknown model: ", model_name)
  )
}

# ------------------------------------------------

# ------------------------------------------------
# 9. Decision-curve helper
# ------------------------------------------------

decision_curve_table <- function(
  prediction_data,
  thresholds
) {
  model_rows <- lapply(
    unique(as.character(prediction_data$Model)),
    function(model_name) {
      model_data <- prediction_data[
        as.character(prediction_data$Model) == model_name,
        ,
        drop = FALSE
      ]

      y <- as.integer(model_data$Y_true)
      p <- clip_prob(model_data$p_mean)
      n <- length(y)
      prevalence <- mean(y)

      dplyr::bind_rows(
        lapply(
          thresholds,
          function(threshold) {
            predicted_positive <- p >= threshold

            true_positive <- sum(
              predicted_positive & y == 1L
            )

            false_positive <- sum(
              predicted_positive & y == 0L
            )

            tibble::tibble(
              Model = model_name,
              Threshold = threshold,
              Net_Benefit =
                true_positive / n -
                false_positive / n *
                threshold / (1 - threshold),
              Prevalence = prevalence
            )
          }
        )
      )
    }
  )

  unique_patient_data <- prediction_data |>
    dplyr::distinct(Row_ID, Y_true)

  prevalence <- mean(unique_patient_data$Y_true)

  reference_rows <- dplyr::bind_rows(
    lapply(
      thresholds,
      function(threshold) {
        tibble::tibble(
          Model = c("Treat all", "Treat none"),
          Threshold = threshold,
          Net_Benefit = c(
            prevalence -
              (1 - prevalence) *
              threshold / (1 - threshold),
            0
          ),
          Prevalence = prevalence
        )
      }
    )
  )

  dplyr::bind_rows(
    model_rows,
    list(reference_rows)
  )
}


# ------------------------------------------------
# 10. External-validation plotting functions
# ------------------------------------------------

plot_cohort_predictor_distributions <- function(
  rpa_raw,
  dusseldorf_raw
) {
  rpa_long <- rpa_raw |>
    dplyr::select(dplyr::all_of(RPA_STD_FEATURES)) |>
    dplyr::mutate(Cohort = "RPA development") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(RPA_STD_FEATURES),
      names_to = "Feature",
      values_to = "Value"
    )

  dusseldorf_long <- dusseldorf_raw |>
    dplyr::transmute(
      Body_Temperature = Body_Temperature_Raw,
      WBC_Count = WBC_Count_Raw,
      Neutrophil_Percentage = Neutrophil_Percentage_Raw,
      Cohort = "Düsseldorf external"
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(RPA_STD_FEATURES),
      names_to = "Feature",
      values_to = "Value"
    )

  plot_data <- dplyr::bind_rows(
    rpa_long,
    dusseldorf_long
  ) |>
    dplyr::filter(is.finite(Value))

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Value,
      color = Cohort,
      fill = Cohort
    )
  ) +
    ggplot2::geom_density(
      alpha = 0.22,
      linewidth = 1.0
    ) +
    ggplot2::facet_wrap(
      ~ Feature,
      scales = "free",
      ncol = 1
    ) +
    ggplot2::scale_color_manual(values = cohort_colors) +
    ggplot2::scale_fill_manual(values = cohort_colors) +
    ggplot2::labs(
      title = "RPA Development and Düsseldorf External Predictor Distributions",
      subtitle = "Observed raw clinical measurements before RPA-trained imputation",
      x = "Clinical measurement",
      y = "Density",
      color = "Cohort",
      fill = "Cohort"
    ) +
    theme_b_risa() +
    ggplot2::theme(legend.position = "bottom")

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$cohort_shift_figs,
      "rpa_dusseldorf_predictor_density"
    ),
    width = 10,
    height = 11
  )
}

plot_standardized_mean_differences <- function(shift_table) {
  plot_object <- ggplot2::ggplot(
    shift_table,
    ggplot2::aes(
      x = Standardized_mean_difference,
      y = reorder(Feature, Standardized_mean_difference)
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dashed",
      color = "gray50"
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "Predictor Distribution Shift from RPA to Düsseldorf",
      subtitle = "Positive values indicate a larger mean in Düsseldorf",
      x = "Standardized mean difference",
      y = NULL
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$cohort_shift_figs,
      "rpa_dusseldorf_standardized_mean_differences"
    ),
    width = 9,
    height = 5.5
  )
}

plot_roc_curves_external <- function(
  prediction_data,
  analysis_id,
  analysis_label
) {
  curve_rows <- list()

  for (model_name in MODEL_ORDER) {
    model_data <- prediction_data |>
      dplyr::filter(Model == model_name)

    roc_object <- pROC::roc(
      response = model_data$Y_true,
      predictor = model_data$p_mean,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )

    curve_rows[[length(curve_rows) + 1L]] <-
      tibble::tibble(
        Model = model_name,
        False_Positive_Rate = 1 - roc_object$specificities,
        True_Positive_Rate = roc_object$sensitivities
      )
  }

  curve_data <- dplyr::bind_rows(curve_rows)

  save_and_print_table(
    curve_data,
    dirs$metrics,
    paste0(analysis_id, "_roc_curve_coordinates"),
    paste0(analysis_label, ": ROC curve coordinates"),
    print_table_object = PRINT_LARGE_SUMMARY_TABLES
  )

  plot_object <- ggplot2::ggplot(
    curve_data,
    ggplot2::aes(
      x = False_Positive_Rate,
      y = True_Positive_Rate,
      color = Model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray50"
    ) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0("Düsseldorf External Validation: ROC Curves"),
      subtitle = analysis_label,
      x = "False positive rate",
      y = "True positive rate",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(analysis_id, "_roc_curves")
    ),
    width = 8,
    height = 7
  )
}

plot_precision_recall_curves_external <- function(
  prediction_data,
  analysis_id,
  analysis_label
) {
  curve_rows <- list()

  for (model_name in MODEL_ORDER) {
    model_data <- prediction_data |>
      dplyr::filter(Model == model_name)

    thresholds <- sort(
      unique(model_data$p_mean),
      decreasing = TRUE
    )

    model_curve <- dplyr::bind_rows(
      lapply(
        thresholds,
        function(threshold) {
          predicted_positive <- model_data$p_mean >= threshold

          tp <- sum(predicted_positive & model_data$Y_true == 1L)
          fp <- sum(predicted_positive & model_data$Y_true == 0L)
          fn <- sum(!predicted_positive & model_data$Y_true == 1L)

          tibble::tibble(
            Model = model_name,
            Threshold = threshold,
            Precision = if (tp + fp > 0L) tp / (tp + fp) else NA_real_,
            Recall = if (tp + fn > 0L) tp / (tp + fn) else NA_real_
          )
        }
      )
    ) |>
      tidyr::drop_na(Precision, Recall)

    curve_rows[[length(curve_rows) + 1L]] <- model_curve
  }

  curve_data <- dplyr::bind_rows(curve_rows)

  save_and_print_table(
    curve_data,
    dirs$metrics,
    paste0(analysis_id, "_precision_recall_curve_coordinates"),
    paste0(analysis_label, ": precision-recall curve coordinates"),
    print_table_object = PRINT_LARGE_SUMMARY_TABLES
  )

  plot_object <- ggplot2::ggplot(
    curve_data,
    ggplot2::aes(
      x = Recall,
      y = Precision,
      color = Model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Precision-Recall Curves",
      subtitle = analysis_label,
      x = "Recall / sensitivity",
      y = "Precision / PPV",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(analysis_id, "_precision_recall_curves")
    ),
    width = 8,
    height = 7
  )
}

plot_calibration_external <- function(
  calibration_data,
  analysis_id,
  analysis_label
) {
  plot_data <- calibration_data |>
    dplyr::filter(Analysis_Set == analysis_id)

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Mean_Predicted,
      y = Observed_Rate,
      color = Model
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray45"
    ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(
      ggplot2::aes(size = n),
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::scale_size_continuous(range = c(2, 6)) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Calibration Curves",
      subtitle = analysis_label,
      x = "Mean predicted probability",
      y = "Observed appendicitis rate",
      color = "Model",
      size = "Bin n"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$calibration_figs,
      paste0(analysis_id, "_calibration_curves")
    ),
    width = 8,
    height = 7
  )
}

plot_uncertainty_external <- function(
  prediction_data,
  analysis_id,
  analysis_label
) {
  plot_data <- prediction_data |>
    dplyr::filter(Analysis_Set == analysis_id) |>
    dplyr::mutate(
      Interval_Width = p_upr - p_lwr
    )

  uncertainty_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Model,
      y = p_sd,
      fill = Model
    )
  ) +
    ggplot2::geom_violin(alpha = 0.75, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.25) +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Predictive Uncertainty",
      subtitle = analysis_label,
      x = "Model",
      y = "Posterior predictive SD",
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    uncertainty_plot,
    file.path(
      dirs$uncertainty_figs,
      paste0(analysis_id, "_posterior_predictive_uncertainty")
    ),
    width = 8,
    height = 7
  )

  interval_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Model,
      y = Interval_Width,
      fill = Model
    )
  ) +
    ggplot2::geom_violin(alpha = 0.75, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.25) +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: 95% Predictive Interval Width",
      subtitle = analysis_label,
      x = "Model",
      y = "95% posterior interval width",
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    interval_plot,
    file.path(
      dirs$uncertainty_figs,
      paste0(analysis_id, "_posterior_interval_width")
    ),
    width = 8,
    height = 7
  )
}

plot_probability_density_external <- function(
  prediction_data,
  analysis_id,
  analysis_label
) {
  plot_data <- prediction_data |>
    dplyr::filter(Analysis_Set == analysis_id)

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = p_mean,
      fill = Outcome_true,
      color = Outcome_true
    )
  ) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 1.1) +
    ggplot2::facet_wrap(~ Model, ncol = 1) +
    ggplot2::scale_color_manual(values = outcome_colors) +
    ggplot2::scale_fill_manual(values = outcome_colors) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Predicted Probability Distributions",
      subtitle = analysis_label,
      x = "Predicted probability of appendicitis",
      y = "Density",
      fill = "True outcome",
      color = "True outcome"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(analysis_id, "_predicted_probability_density")
    ),
    width = 8,
    height = 9
  )
}

plot_risk_coverage_external <- function(
  risk_coverage_data,
  analysis_id,
  analysis_label
) {
  plot_data <- risk_coverage_data |>
    dplyr::filter(Analysis_Set == analysis_id)

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Coverage,
      y = Selective_Risk,
      color = Model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Risk-Coverage Curves",
      subtitle = paste0(analysis_label, "; lower risk is better"),
      x = "Coverage retained",
      y = "Selective risk",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$selective_figs,
      paste0(analysis_id, "_selective_risk_coverage_curve")
    ),
    width = 8,
    height = 7
  )
}

plot_decision_curves_external <- function(
  dca_data,
  analysis_id,
  analysis_label
) {
  plot_data <- dca_data |>
    dplyr::filter(Analysis_Set == analysis_id)

  full_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Threshold,
      y = Net_Benefit,
      color = Model,
      linetype = Model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Decision Curve Analysis",
      subtitle = analysis_label,
      x = "Threshold probability",
      y = "Net benefit",
      color = "Strategy",
      linetype = "Strategy"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    full_plot,
    file.path(
      dirs$dca_figs,
      paste0(analysis_id, "_decision_curve_net_benefit")
    ),
    width = 9,
    height = 7
  )

  zoom_data <- plot_data |>
    dplyr::filter(
      Threshold >= 0.05,
      Threshold <= 0.80
    )

  zoom_plot <- ggplot2::ggplot(
    zoom_data,
    ggplot2::aes(
      x = Threshold,
      y = Net_Benefit,
      color = Model,
      linetype = Model
    )
  ) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Decision Curve, Clinical Range",
      subtitle = paste0(analysis_label, "; thresholds 0.05 to 0.80"),
      x = "Threshold probability",
      y = "Net benefit",
      color = "Strategy",
      linetype = "Strategy"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    zoom_plot,
    file.path(
      dirs$dca_figs,
      paste0(analysis_id, "_decision_curve_net_benefit_zoom")
    ),
    width = 9,
    height = 7
  )
}

plot_confusion_matrices_external <- function(confusion_table) {
  plot_object <- ggplot2::ggplot(
    confusion_table,
    ggplot2::aes(
      x = Predicted,
      y = Outcome_true,
      fill = n
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = n),
      size = 4.5
    ) +
    ggplot2::facet_grid(
      Analysis_Label ~ Model,
      scales = "free"
    ) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Confusion Matrices",
      x = "Predicted outcome",
      y = "Observed outcome",
      fill = "Count"
    ) +
    theme_b_risa(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
    )

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      "dusseldorf_confusion_matrices"
    ),
    width = 13,
    height = 7
  )
}

plot_bootstrap_model_estimates <- function(bootstrap_model_estimates) {
  plot_data <- bootstrap_model_estimates |>
    dplyr::filter(Metric %in% INFERENTIAL_METRICS) |>
    dplyr::mutate(
      Model = factor(Model, levels = MODEL_ORDER),
      Metric_Label = factor(
        Metric_Label,
        levels = unique(Metric_Label)
      )
    )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Model,
      y = Point_estimate,
      color = Model
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      width = 0.2
    ) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::facet_grid(
      Analysis_Label ~ Metric_Label,
      scales = "free_y"
    ) +
    ggplot2::scale_color_manual(values = model_colors) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Bootstrap Model Estimates",
      subtitle = paste0(
        round(CONFIDENCE_LEVEL * 100),
        "% paired patient-bootstrap percentile intervals"
      ),
      x = "Model",
      y = "Metric estimate",
      color = "Model"
    ) +
    theme_b_risa(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$bootstrap_figs,
      "dusseldorf_bootstrap_model_estimates"
    ),
    width = 22,
    height = 8
  )
}

plot_bootstrap_pairwise_forest <- function(bootstrap_pairwise_results) {
  plot_data <- bootstrap_pairwise_results |>
    dplyr::mutate(
      Metric_Label = factor(
        Metric_Label,
        levels = unique(Metric_Label)
      )
    )

  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Point_difference,
      y = Comparison
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray45"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::facet_grid(
      Analysis_Label ~ Metric_Label,
      scales = "free_x"
    ) +
    ggplot2::labs(
      title = "Düsseldorf External Validation: Paired Bootstrap Differences",
      subtitle = "Difference = first model minus second model",
      x = "Metric difference",
      y = "Model comparison"
    ) +
    theme_b_risa(base_size = 10)

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$bootstrap_figs,
      "dusseldorf_bootstrap_pairwise_forest"
    ),
    width = 22,
    height = 9
  )
}



# ------------------------------------------------
# 11. Load and validate the two cohorts
# ------------------------------------------------

rpa_standardized_data <- load_rpa_standardized_dataset(
  RPA_STD_CLEAN_FILE
)

rpa_scaling_parameters <- read_rpa_scaling_parameters(
  RPA_SCALING_PARAMETER_FILE
)

dusseldorf_raw_data <- load_dusseldorf_dataset(
  DUSSELDORF_RAW_FILE
)

dusseldorf_standardized_data <-
  convert_dusseldorf_to_rpa_upstream_scale(
    dusseldorf_raw_data,
    rpa_scaling_parameters
  )

rpa_raw_reconstructed <- reconstruct_rpa_raw_scale(
  rpa_standardized_data,
  rpa_scaling_parameters
)

# ------------------------------------------------
# 12. Dataset, unit, and preprocessing audits
# ------------------------------------------------

outcome_mapping_audit <- dusseldorf_raw_data |>
  dplyr::count(
    diagnosis_original,
    Y,
    Outcome,
    name = "n"
  ) |>
  dplyr::arrange(Y)

wbc_harmonization_audit <- dusseldorf_raw_data |>
  dplyr::mutate(
    Conversion = dplyr::case_when(
      is.na(WBC_Count_Original) ~ "Missing",
      WBC_Count_Original > 1000 ~ "Divided by 1000",
      TRUE ~ "No conversion"
    )
  ) |>
  dplyr::count(Conversion, name = "n") |>
  dplyr::mutate(
    percentage = 100 * n / sum(n)
  )

external_missingness_audit <- tibble::tibble(
  Feature = RPA_STD_FEATURES,
  Missing_n = c(
    sum(is.na(dusseldorf_raw_data$Body_Temperature_Raw)),
    sum(is.na(dusseldorf_raw_data$WBC_Count_Raw)),
    sum(is.na(dusseldorf_raw_data$Neutrophil_Percentage_Raw))
  )
) |>
  dplyr::mutate(
    Total_n = nrow(dusseldorf_raw_data),
    Missing_percent = 100 * Missing_n / Total_n
  )

cohort_summary <- tibble::tibble(
  Cohort = c(
    "RPA development",
    "Düsseldorf external"
  ),
  n = c(
    nrow(rpa_standardized_data),
    nrow(dusseldorf_raw_data)
  ),
  n_no_appendicitis = c(
    sum(rpa_standardized_data$Y == 0L),
    sum(dusseldorf_raw_data$Y == 0L)
  ),
  n_appendicitis = c(
    sum(rpa_standardized_data$Y == 1L),
    sum(dusseldorf_raw_data$Y == 1L)
  ),
  prevalence_appendicitis = c(
    mean(rpa_standardized_data$Y),
    mean(dusseldorf_raw_data$Y)
  )
)

analysis_sample_summary <- tibble::tibble(
  Analysis_Set = c(
    PRIMARY_ANALYSIS_ID,
    COMPLETE_CASE_ANALYSIS_ID
  ),
  Analysis_Label = c(
    PRIMARY_ANALYSIS_LABEL,
    COMPLETE_CASE_ANALYSIS_LABEL
  ),
  n = c(
    nrow(dusseldorf_raw_data),
    sum(dusseldorf_raw_data$Complete_Three_Predictors)
  ),
  n_no_appendicitis = c(
    sum(dusseldorf_raw_data$Y == 0L),
    sum(
      dusseldorf_raw_data$Y == 0L &
        dusseldorf_raw_data$Complete_Three_Predictors
    )
  ),
  n_appendicitis = c(
    sum(dusseldorf_raw_data$Y == 1L),
    sum(
      dusseldorf_raw_data$Y == 1L &
        dusseldorf_raw_data$Complete_Three_Predictors
    )
  )
) |>
  dplyr::mutate(
    prevalence_appendicitis =
      n_appendicitis / n
  )

cohort_shift_objects <- calculate_cohort_shift(
  rpa_raw_reconstructed,
  dusseldorf_raw_data
)

cohort_predictor_summary <-
  cohort_shift_objects$summary

cohort_shift_table <-
  cohort_shift_objects$shift

# Fit the final preprocessing object using the complete RPA development cohort.
# This repeats the exact median-imputation and training-only centering/scaling
# principle used inside the repeated-CV analysis, but now the complete
# RPA cohort is the locked development sample.
final_preprocessing <- fit_preprocess(
  rpa_standardized_data
)

x_rpa_final <- apply_preprocess(
  rpa_standardized_data,
  final_preprocessing
)

x_dusseldorf_final <- apply_preprocess(
  dusseldorf_standardized_data,
  final_preprocessing
)

rpa_preprocessing_audit <- tibble::tibble(
  Feature = RPA_STD_FEATURES,
  Upstream_Raw_Mean = rpa_scaling_parameters$mean_overall,
  Upstream_Raw_SD = rpa_scaling_parameters$sd_overall,
  RPA_Imputation_Median_Upstream_Standardized =
    as.numeric(final_preprocessing$medians[RPA_STD_FEATURES]),
  RPA_Imputation_Median_Raw =
    as.numeric(final_preprocessing$medians[RPA_STD_FEATURES]) *
      rpa_scaling_parameters$sd_overall +
      rpa_scaling_parameters$mean_overall,
  Final_RPA_Center_After_Imputation =
    as.numeric(final_preprocessing$centers[RPA_STD_FEATURES]),
  Final_RPA_Scale_After_Imputation =
    as.numeric(final_preprocessing$scales[RPA_STD_FEATURES])
)

external_imputation_audit <- external_missingness_audit |>
  dplyr::left_join(
    rpa_preprocessing_audit |>
      dplyr::select(
        Feature,
        RPA_Imputation_Median_Upstream_Standardized,
        RPA_Imputation_Median_Raw
      ),
    by = "Feature"
  ) |>
  dplyr::mutate(
    Imputation_source = "RPA development cohort only"
  )

# Verify that the external final matrix is finite and that no Düsseldorf
# outcome was used in preprocessing.
if (any(!is.finite(x_rpa_final))) {
  stop("Non-finite values remain in the final RPA training matrix.")
}

if (any(!is.finite(x_dusseldorf_final))) {
  stop("Non-finite values remain in the final Düsseldorf matrix.")
}

if (nrow(x_dusseldorf_final) != nrow(dusseldorf_raw_data)) {
  stop("Düsseldorf preprocessing changed the patient count.")
}

# Save cleaned and transformed analytic datasets. These patient-level files
# are not printed because they are data outputs rather than summary tables.
safe_fwrite(
  dusseldorf_raw_data,
  file.path(
    dirs$dataset,
    "dusseldorf_cleaned_raw_external_dataset.csv"
  )
)

safe_fwrite(
  dusseldorf_standardized_data,
  file.path(
    dirs$preprocessing,
    "dusseldorf_rpa_upstream_standardized_dataset.csv"
  )
)

safe_fwrite(
  tibble::as_tibble(x_dusseldorf_final) |>
    dplyr::mutate(
      Row_ID = dusseldorf_raw_data$Row_ID,
      Y = dusseldorf_raw_data$Y,
      .before = 1
    ),
  file.path(
    dirs$preprocessing,
    "dusseldorf_final_locked_model_matrix.csv"
  )
)

safe_fwrite(
  tibble::as_tibble(x_rpa_final) |>
    dplyr::mutate(
      Row_ID = rpa_standardized_data$Row_ID,
      Y = rpa_standardized_data$Y,
      .before = 1
    ),
  file.path(
    dirs$preprocessing,
    "rpa_final_locked_training_matrix.csv"
  )
)

# Save and print cohort/preprocessing tables before fitting.
save_and_print_table(
  cohort_summary,
  dirs$dataset,
  "cohort_summary",
  "Development and external cohort summary"
)

save_and_print_table(
  analysis_sample_summary,
  dirs$dataset,
  "external_analysis_sample_summary",
  "Düsseldorf external-validation analysis samples"
)

save_and_print_table(
  outcome_mapping_audit,
  dirs$dataset,
  "dusseldorf_outcome_mapping_audit",
  "Düsseldorf outcome mapping audit"
)

save_and_print_table(
  wbc_harmonization_audit,
  dirs$preprocessing,
  "dusseldorf_wbc_harmonization_audit",
  "Düsseldorf WBC unit harmonization audit"
)

save_and_print_table(
  rpa_scaling_parameters,
  dirs$preprocessing,
  "rpa_upstream_standardization_parameters",
  "RPA upstream raw-to-standardized parameters"
)

save_and_print_table(
  rpa_preprocessing_audit,
  dirs$preprocessing,
  "final_rpa_preprocessing_parameters",
  "Final locked RPA preprocessing parameters"
)

save_and_print_table(
  external_imputation_audit,
  dirs$preprocessing,
  "dusseldorf_external_imputation_audit",
  "Düsseldorf missingness and RPA-trained imputation audit"
)

save_and_print_table(
  cohort_predictor_summary,
  dirs$cohort_shift,
  "rpa_dusseldorf_predictor_summary",
  "RPA and Düsseldorf predictor summaries"
)

save_and_print_table(
  cohort_shift_table,
  dirs$cohort_shift,
  "rpa_dusseldorf_distribution_shift",
  "RPA-to-Düsseldorf predictor distribution shift"
)

# ------------------------------------------------
# 13. Fit locked RPA models and predict Düsseldorf
# ------------------------------------------------

y_rpa <- rpa_standardized_data$Y
y_dusseldorf <- dusseldorf_raw_data$Y

prediction_list <- vector(
  mode = "list",
  length = length(MODEL_SET)
)

status_list <- vector(
  mode = "list",
  length = length(MODEL_SET)
)

posterior_draws_by_model <- vector(
  mode = "list",
  length = length(MODEL_SET)
)

names(posterior_draws_by_model) <- MODEL_SET

for (model_index in seq_along(MODEL_SET)) {
  model_name <- MODEL_SET[model_index]
  fit_seed <- make_external_fit_seed(model_name)

  cat(
    sprintf(
      "\nFitting locked %s model on all RPA patients and predicting Düsseldorf...\n",
      model_name
    )
  )

  set.seed(fit_seed)
  start_time <- Sys.time()

  fit_result <- tryCatch(
    fit_model(
      model_name,
      x_rpa_final,
      y_rpa,
      x_dusseldorf_final
    ),
    error = function(e) {
      structure(
        list(message = e$message),
        class = "model_fit_error"
      )
    }
  )

  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )

  if (inherits(fit_result, "model_fit_error")) {
    failure_row <- tibble::tibble(
      Dataset = DATASET_LABEL,
      Model = model_name,
      Fit_seed = fit_seed,
      Status = "FAILED",
      Failure_message = fit_result$message,
      Runtime_seconds = runtime_seconds
    )

    status_list[[model_index]] <- failure_row

    safe_fwrite(
      dplyr::bind_rows(status_list),
      file.path(
        dirs$status,
        "dusseldorf_external_model_status_partial.csv"
      )
    )

    if (STOP_ON_MODEL_FAILURE) {
      stop(
        "Locked external model failure for ",
        model_name,
        ": ",
        fit_result$message
      )
    }

    next
  }

  posterior_draws <- validate_prediction_draws(
    draws = fit_result$draws,
    expected_draws = expected_draw_count(model_name),
    expected_test_size = nrow(x_dusseldorf_final),
    model_name = model_name
  )

  posterior_draws_by_model[[model_name]] <- posterior_draws

  posterior_summary <- summarise_draws(
    posterior_draws
  )

  prediction_sd_across_patients <- stats::sd(
    posterior_summary$p_mean
  )

  number_unique_predictions <- dplyr::n_distinct(
    round(posterior_summary$p_mean, digits = 12)
  )

  degenerate_prediction <-
    !is.finite(prediction_sd_across_patients) ||
    prediction_sd_across_patients <= DEGENERATE_PREDICTION_TOL ||
    number_unique_predictions <= 1L

  prediction_data <- tibble::tibble(
    Dataset = DATASET_LABEL,
    Row_ID = dusseldorf_raw_data$Row_ID,
    Model = model_name,
    Y_true = y_dusseldorf,
    Outcome_true = dusseldorf_raw_data$Outcome,
    Complete_Three_Predictors =
      dusseldorf_raw_data$Complete_Three_Predictors,
    Missing_Predictor_Count =
      dusseldorf_raw_data$Missing_Predictor_Count,
    p_mean = posterior_summary$p_mean,
    p_sd = posterior_summary$p_sd,
    p_lwr = posterior_summary$p_lwr,
    p_upr = posterior_summary$p_upr
  )

  status_base <- tibble::tibble(
    Dataset = DATASET_LABEL,
    Model = model_name,
    Fit_seed = fit_seed,
    Status = if (degenerate_prediction) "DEGENERATE" else "OK",
    Failure_message = NA_character_,
    Runtime_seconds = runtime_seconds,
    RPA_training_n = nrow(x_rpa_final),
    Dusseldorf_test_n = nrow(x_dusseldorf_final),
    Posterior_draw_rows = nrow(posterior_draws),
    Posterior_draw_columns = ncol(posterior_draws),
    Probability_min = min(prediction_data$p_mean),
    Probability_max = max(prediction_data$p_mean),
    Probability_mean = mean(prediction_data$p_mean),
    Probability_sd_across_patients =
      prediction_sd_across_patients,
    Number_unique_predictions =
      number_unique_predictions,
    Number_predicted_positive = sum(
      prediction_data$p_mean >= CLASSIFICATION_THRESHOLD
    ),
    Degenerate_prediction = degenerate_prediction
  )

  diagnostic_row <- tibble::as_tibble(
    fit_result$diagnostics
  )

  status_list[[model_index]] <- dplyr::bind_cols(
    status_base,
    diagnostic_row
  )

  prediction_list[[model_index]] <- prediction_data

  if (
    degenerate_prediction &&
    STOP_ON_DEGENERATE_PREDICTIONS
  ) {
    safe_fwrite(
      dplyr::bind_rows(status_list),
      file.path(
        dirs$status,
        "dusseldorf_external_model_status_partial.csv"
      )
    )

    stop(
      "Degenerate Düsseldorf predictions detected for ",
      model_name,
      "."
    )
  }

  rm(
    fit_result,
    posterior_draws,
    posterior_summary,
    prediction_data,
    status_base,
    diagnostic_row
  )

  gc(FALSE)
}

locked_predictions <- dplyr::bind_rows(
  prediction_list
)

model_status <- dplyr::bind_rows(
  status_list
)

if (nrow(locked_predictions) != nrow(dusseldorf_raw_data) * length(MODEL_SET)) {
  stop("The locked external prediction row count is incorrect.")
}

prediction_completeness <- locked_predictions |>
  dplyr::count(Row_ID, name = "Model_count")

if (any(prediction_completeness$Model_count != length(MODEL_SET))) {
  stop("Each Düsseldorf patient must have one prediction from every model.")
}

# ------------------------------------------------
# 14. Primary and complete-case analyses
# ------------------------------------------------

analysis_predictions <- dplyr::bind_rows(
  locked_predictions |>
    dplyr::mutate(
      Analysis_Set = PRIMARY_ANALYSIS_ID,
      Analysis_Label = PRIMARY_ANALYSIS_LABEL,
      .before = 1
    ),
  locked_predictions |>
    dplyr::filter(Complete_Three_Predictors) |>
    dplyr::mutate(
      Analysis_Set = COMPLETE_CASE_ANALYSIS_ID,
      Analysis_Label = COMPLETE_CASE_ANALYSIS_LABEL,
      .before = 1
    )
) |>
  dplyr::mutate(
    Analysis_Set = factor(
      Analysis_Set,
      levels = c(
        PRIMARY_ANALYSIS_ID,
        COMPLETE_CASE_ANALYSIS_ID
      )
    ),
    Model = factor(Model, levels = MODEL_ORDER)
  )

overall_metrics <- analysis_predictions |>
  dplyr::group_by(
    Analysis_Set,
    Analysis_Label,
    Dataset,
    Model
  ) |>
  dplyr::group_modify(
    ~ calc_metrics(
      y = .x$Y_true,
      p_mean = .x$p_mean,
      p_sd = .x$p_sd,
      p_lwr = .x$p_lwr,
      p_upr = .x$p_upr
    )
  ) |>
  dplyr::ungroup()

confusion_table <- make_confusion_table(
  analysis_predictions
)

calibration_bins <- make_calibration_bins(
  analysis_predictions
)

risk_coverage_data <- make_risk_coverage_data(
  analysis_predictions
)

decision_curve_data <- make_decision_curve_data(
  analysis_predictions
)

# ------------------------------------------------
# 15. Paired patient bootstrap
# ------------------------------------------------

bootstrap_results <- bootstrap_external_predictions(
  prediction_data = analysis_predictions,
  bootstrap_reps = BOOTSTRAP_REPS,
  bootstrap_seed = BOOTSTRAP_SEED
)

bootstrap_model_estimates <-
  bootstrap_results$model_estimates

bootstrap_pairwise_results <-
  bootstrap_results$pairwise_results

bootstrap_draws_long <-
  bootstrap_results$bootstrap_draws_long

# ------------------------------------------------
# 16. Save and print external-validation outputs
# ------------------------------------------------

save_and_print_table(
  model_status,
  dirs$status,
  "dusseldorf_external_model_status",
  "Locked external model-status audit"
)

save_and_print_table(
  prediction_completeness,
  dirs$status,
  "dusseldorf_external_prediction_completeness",
  "Düsseldorf prediction completeness audit",
  print_table_object = PRINT_LARGE_SUMMARY_TABLES
)

save_and_print_table(
  overall_metrics,
  dirs$metrics,
  "dusseldorf_external_overall_metrics",
  "Düsseldorf external-validation metrics"
)

save_and_print_table(
  confusion_table,
  dirs$metrics,
  "dusseldorf_external_confusion_matrices",
  "Düsseldorf confusion matrices"
)

save_and_print_table(
  calibration_bins,
  dirs$calibration,
  "dusseldorf_external_calibration_bins",
  "Düsseldorf calibration-bin table",
  print_table_object = PRINT_LARGE_SUMMARY_TABLES
)

save_and_print_table(
  risk_coverage_data,
  dirs$selective,
  "dusseldorf_external_risk_coverage_curve",
  "Düsseldorf risk-coverage table",
  print_table_object = PRINT_LARGE_SUMMARY_TABLES
)

save_and_print_table(
  decision_curve_data,
  dirs$dca,
  "dusseldorf_external_decision_curve_net_benefit",
  "Düsseldorf decision-curve table",
  print_table_object = PRINT_LARGE_SUMMARY_TABLES
)

save_and_print_table(
  bootstrap_model_estimates,
  dirs$bootstrap,
  "dusseldorf_external_patient_bootstrap_model_estimates",
  "Düsseldorf patient-bootstrap model estimates"
)

save_and_print_table(
  bootstrap_pairwise_results,
  dirs$bootstrap,
  "dusseldorf_external_patient_bootstrap_pairwise_comparisons",
  "Düsseldorf paired patient-bootstrap model comparisons"
)

# The per-patient prediction file is saved but not printed because it is a
# patient-level output rather than a summary table.
safe_fwrite(
  analysis_predictions,
  file.path(
    dirs$predictions,
    "dusseldorf_external_locked_predictions.csv"
  )
)

# Bootstrap draws are saved for audit but not printed because they contain
# many thousands of rows.
safe_fwrite(
  bootstrap_draws_long,
  file.path(
    dirs$bootstrap,
    "dusseldorf_external_patient_bootstrap_draws_long.csv"
  )
)

# ------------------------------------------------
# 17. Figures
# ------------------------------------------------

plot_cohort_predictor_distributions(
  rpa_raw_reconstructed,
  dusseldorf_raw_data
)

plot_standardized_mean_differences(
  cohort_shift_table
)

for (analysis_id in c(
  PRIMARY_ANALYSIS_ID,
  COMPLETE_CASE_ANALYSIS_ID
)) {
  analysis_label <- if (
    analysis_id == PRIMARY_ANALYSIS_ID
  ) {
    PRIMARY_ANALYSIS_LABEL
  } else {
    COMPLETE_CASE_ANALYSIS_LABEL
  }

  analysis_data <- analysis_predictions |>
    dplyr::filter(
      as.character(Analysis_Set) == analysis_id
    )

  plot_roc_curves_external(
    analysis_data,
    analysis_id,
    analysis_label
  )

  plot_precision_recall_curves_external(
    analysis_data,
    analysis_id,
    analysis_label
  )

  plot_calibration_external(
    calibration_bins,
    analysis_id,
    analysis_label
  )

  plot_uncertainty_external(
    analysis_predictions,
    analysis_id,
    analysis_label
  )

  plot_probability_density_external(
    analysis_predictions,
    analysis_id,
    analysis_label
  )

  plot_risk_coverage_external(
    risk_coverage_data,
    analysis_id,
    analysis_label
  )

  plot_decision_curves_external(
    decision_curve_data,
    analysis_id,
    analysis_label
  )
}

plot_confusion_matrices_external(
  confusion_table
)

plot_bootstrap_model_estimates(
  bootstrap_model_estimates
)

plot_bootstrap_pairwise_forest(
  bootstrap_pairwise_results
)

# ------------------------------------------------
# 18. Manifest, R objects, and session information
# ------------------------------------------------

method_manifest <- tibble::tibble(
  Item = c(
    "Development dataset",
    "External dataset",
    "Development predictors",
    "RPA upstream scaling source",
    "External WBC harmonization",
    "External primary analysis",
    "External sensitivity analysis",
    "External preprocessing",
    "External-outcome use in preprocessing",
    "Models",
    "BLR method",
    "BLR priors",
    "BLDA method",
    "BLDA posterior draws",
    "BLDA mean prior kappa",
    "BLDA inverse-Wishart prior df",
    "BLDA prior covariance scale",
    "BLDA class-prior Beta hyperparameters",
    "BART method",
    "BART posterior draws",
    "BART burn-in",
    "BART number of trees",
    "Classification threshold",
    "Selective coverage",
    "Decision-curve thresholds",
    "Bootstrap replicates",
    "Bootstrap seed",
    "External fit seed",
    "Primary inferential note"
  ),
  Value = c(
    RPA_STD_CLEAN_FILE,
    DUSSELDORF_RAW_FILE,
    paste(RPA_STD_FEATURES, collapse = ", "),
    if (file.exists(RPA_SCALING_PARAMETER_FILE)) {
      RPA_SCALING_PARAMETER_FILE
    } else {
      "Validated embedded fallback parameters"
    },
    "Values greater than 100 divided by 1000",
    paste0(
      nrow(dusseldorf_raw_data),
      " patients after RPA-trained median imputation"
    ),
    paste0(
      sum(dusseldorf_raw_data$Complete_Three_Predictors),
      " complete cases"
    ),
    paste0(
      "Raw Düsseldorf values transformed with RPA upstream means/SDs; ",
      "missing values imputed with full-RPA medians; final centering and ",
      "scaling estimated from full RPA and applied unchanged"
    ),
    "None",
    paste(MODEL_SET, collapse = ", "),
    "MAP plus Laplace posterior approximation",
    paste0(
      "Intercept SD = ",
      BLR_INTERCEPT_PRIOR_SD,
      "; coefficient SD = ",
      BLR_COEFFICIENT_PRIOR_SD
    ),
    "Normal-Inverse-Wishart shared-covariance Bayesian LDA",
    N_POST_DRAWS_BLDA,
    BLDA_MEAN_PRIOR_KAPPA,
    paste0("p + ", BLDA_IW_DF_OFFSET),
    BLDA_PRIOR_COVARIANCE_SCALE,
    paste0(
      "Beta(",
      BLDA_CLASS_PRIOR_ALPHA,
      ", ",
      BLDA_CLASS_PRIOR_BETA,
      ")"
    ),
    "Probit Bayesian additive regression trees",
    N_POST_DRAWS_BART,
    N_BURN_BART,
    BART_NTREE,
    CLASSIFICATION_THRESHOLD,
    SELECTIVE_COVERAGE,
    paste0(
      min(DCA_THRESHOLDS),
      " to ",
      max(DCA_THRESHOLDS),
      " by ",
      DCA_THRESHOLDS[2] - DCA_THRESHOLDS[1]
    ),
    BOOTSTRAP_REPS,
    BOOTSTRAP_SEED,
    EXTERNAL_FIT_SEED,
    paste0(
      "Paired patient bootstrap conditions on locked predictions and ",
      "does not refit the RPA models"
    )
  )
)

save_and_print_table(
  method_manifest,
  dirs$manifests,
  "dusseldorf_external_method_manifest",
  "Düsseldorf external-validation method manifest"
)

saveRDS(
  list(
    rpa_scaling_parameters = rpa_scaling_parameters,
    final_preprocessing = final_preprocessing,
    cohort_summary = cohort_summary,
    analysis_sample_summary = analysis_sample_summary,
    cohort_predictor_summary = cohort_predictor_summary,
    cohort_shift_table = cohort_shift_table,
    model_status = model_status,
    locked_predictions = locked_predictions,
    analysis_predictions = analysis_predictions,
    overall_metrics = overall_metrics,
    confusion_table = confusion_table,
    calibration_bins = calibration_bins,
    risk_coverage_data = risk_coverage_data,
    decision_curve_data = decision_curve_data,
    bootstrap_model_estimates = bootstrap_model_estimates,
    bootstrap_pairwise_results = bootstrap_pairwise_results,
    posterior_draws_by_model = posterior_draws_by_model,
    method_manifest = method_manifest
  ),
  file.path(
    dirs$r_objects,
    "dusseldorf_external_validation_complete_analysis.rds"
  ),
  compress = "xz"
)

saveRDS(
  posterior_draws_by_model,
  file.path(
    dirs$r_objects,
    "dusseldorf_external_posterior_prediction_draws.rds"
  ),
  compress = "xz"
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    RESULT_ROOT,
    "sessionInfo.txt"
  )
)

# ------------------------------------------------
# 19. Final validation and summary
# ------------------------------------------------

if (any(model_status$Status != "OK")) {
  stop("At least one locked external model has a non-OK status.")
}

if (any(!is.finite(overall_metrics$auc))) {
  stop("At least one external AUC is non-finite.")
}

cat("\n============================================================\n")
cat("B-RISA DÜSSELDORF LOCKED EXTERNAL VALIDATION COMPLETE\n")
cat("============================================================\n")

cat("\nExternal analysis samples:\n")
print(analysis_sample_summary)

cat("\nOverall external-validation metrics:\n")
print(overall_metrics, n = Inf, width = Inf)

cat("\nModel-status audit:\n")
print(model_status, n = Inf, width = Inf)

cat("\nTables saved in:\n")
cat(TABLE_ROOT, "\n")

cat("\nPDF figures saved in:\n")
cat(FIGURE_PDF_ROOT, "\n")

cat("\nPNG figures saved in:\n")
cat(FIGURE_PNG_ROOT, "\n")

cat("\nR objects saved in:\n")
cat(R_OBJECT_ROOT, "\n")

cat("============================================================\n")

