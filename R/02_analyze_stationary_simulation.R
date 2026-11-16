# ================================================================
# B-RISA Stationary Simulation Analysis
#
# Fits BLR, BLDA, and BART to the generated stationary factorial datasets.
# Training uses Y_observed; evaluation uses Y_true. The script produces
# the simulation metrics, audits, figures, checkpoints, and summaries used
# in the manuscript.
#
# Run from the repository root, from the R/ directory, or set BRISA_ROOT
# to the repository root before running.
# ================================================================

rm(list = ls())
gc()

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

SOURCE_PROJECT_DIR <- PROJECT_ROOT
INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "results",
  "clinical_calibrated_simulation_by_n",
  "data_by_sample_size"
)

OUTPUT_ROOT <- file.path(PROJECT_ROOT, "analysis_outputs", "stationary_simulation")
RESULT_ROOT <- file.path(OUTPUT_ROOT, "results")
FIGURE_ROOT <- file.path(OUTPUT_ROOT, "figures")
CHECKPOINT_ROOT <- file.path(RESULT_ROOT, "checkpoints")
REPRO_ROOT <- file.path(RESULT_ROOT, "reproducibility")

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------

ANALYSIS_VERSION <- "B-RISA_stationary_simulation_final"

INSTALL_MISSING_PACKAGES <- FALSE

# Set TRUE only when you intentionally want to erase current analysis
# checkpoints and restart the entire simulation analysis.
CLEAN_OUTPUT_ROOT <- FALSE

# TRUE allows the script to resume from current analysis RDS checkpoints.
RESUME_ANALYSIS <- TRUE

BASE_SEED <- 20260618L

N_FILES <- c(100L, 500L, 1000L)
TRAIN_PROP <- 0.70

MODEL_SET <- c("BLR", "BLDA", "BART")

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

# Genuine Bayesian LDA priors.
#
# Because every training partition is standardized using training-only
# moments, a zero-centered prior for both class means and an identity-
# centered prior for the shared covariance are natural weak defaults.
#
# nu0 = p + 2 and Lambda0 = I imply E[Sigma] = I under the
# inverse-Wishart parameterization used here.
BLDA_MEAN_PRIOR_KAPPA <- 0.01
BLDA_IW_DF_OFFSET <- 2L
BLDA_PRIOR_COVARIANCE_SCALE <- 1.0
BLDA_CLASS_PRIOR_ALPHA <- 0.5
BLDA_CLASS_PRIOR_BETA <- 0.5
BLDA_PRECISION_EIGEN_FLOOR <- 1e-10

CLASSIFICATION_THRESHOLD <- 0.50
SELECTIVE_COVERAGE <- 0.80
SELECTIVE_COVERAGE_GRID <- seq(0.50, 1.00, by = 0.05)
DCA_THRESHOLDS <- seq(0.01, 0.99, by = 0.01)
CALIBRATION_BINS <- 10L
EPS <- 1e-6

DEGENERATE_PREDICTION_TOL <- 1e-12

STOP_ON_MODEL_FAILURE <- FALSE
STOP_ON_DEGENERATE_PREDICTIONS <- FALSE

SAVE_PREDICTION_SAMPLE <- TRUE
MAX_PREDICTIONS_PER_UNIT <- 40L
CHECKPOINT_EVERY_UNITS <- 100L
SAVE_FULL_POSTERIOR_DRAWS <- FALSE

# For a short test, set an integer such as 10L.
# For the full analysis, retain Inf.
MAX_UNITS_PER_N <- Inf

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
  "viridis"
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
    install.packages(
      "BART",
      dependencies = TRUE
    )
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
  library(viridis)
})

# ------------------------------------------------
# 3. Output folders
# ------------------------------------------------

if (CLEAN_OUTPUT_ROOT && dir.exists(OUTPUT_ROOT)) {
  unlink(
    OUTPUT_ROOT,
    recursive = TRUE,
    force = TRUE
  )
}

COMBINED_FIGURE_ROOT <- file.path(
  FIGURE_ROOT,
  "combined"
)

make_sample_figure_dirs <- function(n_value) {
  sample_root <- file.path(
    FIGURE_ROOT,
    paste0("n_", n_value)
  )

  list(
    root = sample_root,
    performance = file.path(sample_root, "performance"),
    calibration = file.path(sample_root, "calibration"),
    uncertainty = file.path(sample_root, "uncertainty"),
    selective_prediction = file.path(sample_root, "selective_prediction"),
    decision_curve = file.path(sample_root, "decision_curve"),
    robustness = file.path(sample_root, "robustness")
  )
}

sample_figure_dirs <- setNames(
  lapply(N_FILES, make_sample_figure_dirs),
  as.character(N_FILES)
)

dirs <- list(
  result_base = RESULT_ROOT,
  checkpoints = CHECKPOINT_ROOT,
  manifests = REPRO_ROOT,
  performance_figs = file.path(COMBINED_FIGURE_ROOT, "performance"),
  calibration_figs = file.path(COMBINED_FIGURE_ROOT, "calibration"),
  uncertainty_figs = file.path(COMBINED_FIGURE_ROOT, "uncertainty"),
  selective_figs = file.path(COMBINED_FIGURE_ROOT, "selective_prediction"),
  dca_figs = file.path(COMBINED_FIGURE_ROOT, "decision_curve"),
  robustness_figs = file.path(COMBINED_FIGURE_ROOT, "robustness")
)

for (directory in unlist(dirs, use.names = FALSE)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

for (sample_dirs in sample_figure_dirs) {
  for (directory in unlist(sample_dirs, use.names = FALSE)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

cat(
  "Input directory:\n",
  INPUT_DIR,
  "\n\nOutput root:\n",
  OUTPUT_ROOT,
  "\n"
)


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
  ggplot2::ggsave(
    paste0(path_without_extension, ".pdf"),
    plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )

  ggplot2::ggsave(
    paste0(path_without_extension, ".png"),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

compute_dca_zoom_limits <- function(
  dca_data,
  threshold_min = 0.05,
  threshold_max = 0.80,
  excluded_strategy = "Treat all",
  padding_fraction = 0.08,
  minimum_span = 0.05
) {
  reference_values <- dca_data |>
    dplyr::filter(
      Threshold >= threshold_min,
      Threshold <= threshold_max,
      Strategy != excluded_strategy,
      is.finite(Mean_Net_Benefit)
    ) |>
    dplyr::pull(
      Mean_Net_Benefit
    )

  if (length(reference_values) == 0L) {
    return(
      c(
        -0.05,
        0.50
      )
    )
  }

  lower_limit <- min(
    c(
      0,
      reference_values
    ),
    na.rm = TRUE
  )

  upper_limit <- max(
    c(
      0,
      reference_values
    ),
    na.rm = TRUE
  )

  observed_span <- max(
    upper_limit - lower_limit,
    minimum_span
  )

  c(
    lower_limit -
      padding_fraction *
      observed_span,
    upper_limit +
      padding_fraction *
      observed_span
  )
}

safe_fwrite <- function(x, path) {
  data.table::fwrite(
    as.data.frame(x),
    path,
    na = "NA"
  )
}

append_deduplicated_rds <- function(
  new_data,
  path,
  key_columns
) {
  if (is.null(new_data) || nrow(new_data) == 0L) {
    return(invisible(NULL))
  }

  combined <- as.data.frame(new_data)

  if (file.exists(path)) {
    old_data <- readRDS(path)
    combined <- dplyr::bind_rows(
      old_data,
      combined
    )
  }

  if (length(key_columns) > 0L) {
    duplicate_index <- duplicated(
      combined[, key_columns, drop = FALSE],
      fromLast = TRUE
    )
    combined <- combined[!duplicate_index, , drop = FALSE]
  }

  saveRDS(
    combined,
    path,
    compress = "xz"
  )

  invisible(NULL)
}

read_rds_rows <- function(paths) {
  existing_paths <- paths[file.exists(paths)]
  if (length(existing_paths) == 0L) {
    return(tibble::tibble())
  }

  dplyr::bind_rows(
    lapply(
      existing_paths,
      readRDS
    )
  )
}

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
    stop(
      model_name,
      " returned non-finite posterior predictions."
    )
  }

  if (any(draws < 0 | draws > 1)) {
    stop(
      model_name,
      " returned values outside [0, 1]."
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

# ------------------------------------------------
# 5. Clinical metric functions, unchanged in method
# ------------------------------------------------

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

calc_ece <- function(
  y,
  p,
  n_bins = CALIBRATION_BINS
) {
  y <- as.integer(y)
  p <- clip_prob(as.numeric(p))

  bins <- cut(
    p,
    breaks = seq(
      0,
      1,
      length.out = n_bins + 1L
    ),
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
      as.integer(
        ceiling(
          target_coverage *
            length(uncertainty)
        )
      )
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
    stop(
      "No finite observations were available for metric calculation."
    )
  }

  yhat <- as.integer(
    p_mean >= threshold
  )

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
    2 * ppv * sensitivity /
      (ppv + sensitivity)
  } else {
    NA_real_
  }

  retained_index <- select_least_uncertain(
    p_sd,
    coverage
  )

  selective_risk <- if (
    length(retained_index) > 0L
  ) {
    mean(
      yhat[retained_index] !=
        y[retained_index]
    )
  } else {
    NA_real_
  }

  calibration <-
    calc_calibration_intercept_slope(
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
        (1 - y) *
        log(1 - p_mean)
    ),
    calibration_intercept =
      unname(calibration["intercept"]),
    calibration_slope =
      unname(calibration["slope"]),
    ece = calc_ece(y, p_mean),
    mean_uncertainty = mean(p_sd),
    mean_interval_width =
      mean(p_upr - p_lwr),
    selective_coverage =
      length(retained_index) /
      length(y),
    selective_risk =
      selective_risk
  )
}

risk_coverage_table <- function(
  y,
  p_mean,
  p_sd,
  coverage_grid =
    SELECTIVE_COVERAGE_GRID
) {
  y <- as.integer(y)
  p_mean <- clip_prob(as.numeric(p_mean))
  p_sd <- as.numeric(p_sd)

  predicted_class <- as.integer(
    p_mean >= CLASSIFICATION_THRESHOLD
  )

  dplyr::bind_rows(
    lapply(
      coverage_grid,
      function(target_coverage) {
        retained_index <-
          select_least_uncertain(
            p_sd,
            target_coverage
          )

        tibble::tibble(
          Target_Coverage =
            target_coverage,
          Coverage =
            length(retained_index) /
            length(y),
          Selective_Risk =
            mean(
              predicted_class[
                retained_index
              ] !=
                y[retained_index]
            )
        )
      }
    )
  )
}

decision_curve_table_multi <- function(
  prediction_list,
  y,
  thresholds = DCA_THRESHOLDS
) {
  y <- as.integer(y)
  n <- length(y)
  prevalence <- mean(y)

  model_rows <- dplyr::bind_rows(
    lapply(
      names(prediction_list),
      function(model_name) {
        p <- clip_prob(
          prediction_list[[model_name]]
        )

        dplyr::bind_rows(
          lapply(
            thresholds,
            function(threshold) {
              predicted_positive <-
                p >= threshold

              true_positive <- sum(
                predicted_positive &
                  y == 1L
              )

              false_positive <- sum(
                predicted_positive &
                  y == 0L
              )

              tibble::tibble(
                Strategy = model_name,
                Threshold = threshold,
                Net_Benefit =
                  true_positive / n -
                  false_positive / n *
                  threshold /
                  (1 - threshold),
                Prevalence = prevalence
              )
            }
          )
        )
      }
    )
  )

  reference_rows <- dplyr::bind_rows(
    lapply(
      thresholds,
      function(threshold) {
        tibble::tibble(
          Strategy = c(
            "Treat all",
            "Treat none"
          ),
          Threshold = threshold,
          Net_Benefit = c(
            prevalence -
              (1 - prevalence) *
              threshold /
              (1 - threshold),
            0
          ),
          Prevalence = prevalence
        )
      }
    )
  )

  dplyr::bind_rows(
    model_rows,
    reference_rows
  ) |>
    dplyr::mutate(
      Standardized_Net_Benefit =
        Net_Benefit / Prevalence
    )
}

# ------------------------------------------------
# 6. Input validation and preprocessing
# ------------------------------------------------

input_files <- tibble::tibble(
  n_generated = N_FILES,
  file_path = file.path(
    INPUT_DIR,
    paste0(
      "generated_data_n",
      N_FILES,
      ".csv"
    )
  )
)

missing_input <- input_files$file_path[
  !file.exists(input_files$file_path)
]

if (length(missing_input) > 0L) {
  stop(
    "Missing input files:\n",
    paste(
      missing_input,
      collapse = "\n"
    )
  )
}

required_metadata_columns <- c(
  "scenario_id",
  "scenario_name",
  "rep_id",
  "n_train",
  "p",
  "prevalence_target",
  "boundary",
  "covariance",
  "label_noise",
  "missing_rate_target",
  "calibration_source",
  "Y_true",
  "Y_observed"
)

validate_input_schema <- function(data) {
  missing_columns <- setdiff(
    required_metadata_columns,
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Simulation input is missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

fit_preprocess <- function(
  train_data,
  feature_columns
) {
  medians <- vapply(
    train_data[
      ,
      feature_columns,
      drop = FALSE
    ],
    function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }

      stats::median(
        x,
        na.rm = TRUE
      )
    },
    numeric(1)
  )

  if (any(!is.finite(medians))) {
    stop(
      "At least one predictor is entirely missing in the training data."
    )
  }

  imputed <- train_data[
    ,
    feature_columns,
    drop = FALSE
  ]

  for (feature_name in feature_columns) {
    missing_index <-
      is.na(imputed[[feature_name]])

    if (any(missing_index)) {
      imputed[[feature_name]][
        missing_index
      ] <- medians[[feature_name]]
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

apply_preprocess <- function(
  data,
  feature_columns,
  preprocessing
) {
  x <- data[
    ,
    feature_columns,
    drop = FALSE
  ]

  for (feature_name in feature_columns) {
    missing_index <-
      is.na(x[[feature_name]])

    if (any(missing_index)) {
      x[[feature_name]][
        missing_index
      ] <-
        preprocessing$medians[[feature_name]]
    }

    x[[feature_name]] <- (
      x[[feature_name]] -
        preprocessing$centers[[feature_name]]
    ) /
      preprocessing$scales[[feature_name]]
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  colnames(x) <- feature_columns

  if (any(!is.finite(x))) {
    stop(
      "Non-finite values remain after preprocessing."
    )
  }

  x
}

stratified_split <- function(
  data,
  stratify_column = "Y_true",
  train_prop = TRAIN_PROP,
  seed = BASE_SEED
) {
  y <- as.integer(
    data[[stratify_column]]
  )

  if (
    length(unique(y)) != 2L ||
    !all(y %in% c(0L, 1L))
  ) {
    stop(
      "The stratification outcome must contain both 0 and 1."
    )
  }

  set.seed(seed)

  index_0 <- which(y == 0L)
  index_1 <- which(y == 1L)

  if (
    length(index_0) < 3L ||
    length(index_1) < 3L
  ) {
    stop(
      "Each true outcome class requires at least three observations."
    )
  }

  training_count <- function(class_count) {
    max(
      2L,
      min(
        class_count - 1L,
        as.integer(
          floor(
            class_count *
              train_prop
          )
        )
      )
    )
  }

  n_train_0 <- training_count(
    length(index_0)
  )

  n_train_1 <- training_count(
    length(index_1)
  )

  train_index <- c(
    sample(
      index_0,
      size = n_train_0,
      replace = FALSE
    ),
    sample(
      index_1,
      size = n_train_1,
      replace = FALSE
    )
  )

  train_index <- sample(train_index)
  test_index <- setdiff(
    seq_len(nrow(data)),
    train_index
  )

  if (
    length(unique(y[train_index])) < 2L ||
    length(unique(y[test_index])) < 2L
  ) {
    stop(
      "The stratified split did not retain both true classes."
    )
  }

  list(
    train_index = train_index,
    test_index = test_index
  )
}

# ------------------------------------------------
# 7. Exact corrected clinical BLR
# ------------------------------------------------

softplus_stable <- function(x) {
  pmax(x, 0) +
    log1p(exp(-abs(x)))
}

nearest_positive_definite <- function(
  covariance_matrix,
  eigen_floor =
    BLR_HESSIAN_EIGEN_FLOOR
) {
  covariance_matrix <-
    as.matrix(covariance_matrix)

  covariance_matrix <- (
    covariance_matrix +
      t(covariance_matrix)
  ) / 2

  decomposition <- eigen(
    covariance_matrix,
    symmetric = TRUE
  )

  decomposition$values[
    decomposition$values <
      eigen_floor
  ] <- eigen_floor

  repaired <-
    decomposition$vectors %*%
    diag(
      decomposition$values,
      nrow =
        length(
          decomposition$values
        )
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
  n_draws =
    N_POST_DRAWS_BLR
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  if (
    length(unique(y_train)) <
      2L
  ) {
    stop(
      "BLR training labels contain only one class."
    )
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
    rep(
      0,
      parameter_count - 1L
    )
  )

  prior_sd <- c(
    BLR_INTERCEPT_PRIOR_SD,
    rep(
      BLR_COEFFICIENT_PRIOR_SD,
      parameter_count - 1L
    )
  )

  prior_precision <-
    1 / prior_sd^2

  negative_log_posterior <-
    function(beta) {
      eta <- as.vector(
        X_train %*% beta
      )

      negative_log_likelihood <-
        sum(
          softplus_stable(eta) -
            y_train * eta
        )

      negative_log_prior <-
        0.5 *
        sum(
          prior_precision *
            (
              beta -
                prior_mean
            )^2
        )

      negative_log_likelihood +
        negative_log_prior
    }

  gradient_negative_log_posterior <-
    function(beta) {
      eta <- as.vector(
        X_train %*% beta
      )

      probability <- plogis(eta)

      as.vector(
        crossprod(
          X_train,
          probability -
            y_train
        )
      ) +
        prior_precision *
        (
          beta -
            prior_mean
        )
    }

  optimization <- stats::optim(
    par = prior_mean,
    fn = negative_log_posterior,
    gr =
      gradient_negative_log_posterior,
    method = "BFGS",
    control = list(
      maxit = BLR_OPTIM_MAXIT,
      reltol =
        BLR_OPTIM_RELTOL
    ),
    hessian = FALSE
  )

  beta_map <- optimization$par

  gradient_at_map <-
    gradient_negative_log_posterior(
      beta_map
    )

  maximum_absolute_gradient <-
    max(
      abs(gradient_at_map)
    )

  converged <- (
    optimization$convergence == 0L ||
      maximum_absolute_gradient <=
        BLR_MAX_ABS_GRADIENT
  )

  if (!converged) {
    stop(
      "BLR optimization did not converge. Code = ",
      optimization$convergence,
      "; max absolute gradient = ",
      signif(
        maximum_absolute_gradient,
        6
      ),
      "."
    )
  }

  fitted_eta <- as.vector(
    X_train %*% beta_map
  )

  fitted_probability <-
    plogis(fitted_eta)

  logistic_weights <- pmax(
    fitted_probability *
      (
        1 -
          fitted_probability
      ),
    1e-8
  )

  hessian <- crossprod(
    X_train,
    X_train *
      logistic_weights
  ) +
    diag(
      prior_precision,
      nrow =
        parameter_count
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

  if (
    any(
      !is.finite(
        hessian_eigenvalues
      )
    )
  ) {
    stop(
      "BLR Hessian has non-finite eigenvalues."
    )
  }

  posterior_covariance <-
    tryCatch(
      solve(hessian),
      error = function(e) NULL
    )

  if (
    is.null(
      posterior_covariance
    ) ||
    any(
      !is.finite(
        posterior_covariance
      )
    )
  ) {
    stop(
      "BLR Hessian inversion failed."
    )
  }

  posterior_covariance <-
    nearest_positive_definite(
      posterior_covariance
    )

  beta_draws <- MASS::mvrnorm(
    n = n_draws,
    mu = beta_map,
    Sigma =
      posterior_covariance
  )

  beta_draws <- matrix(
    beta_draws,
    nrow = n_draws,
    ncol =
      parameter_count
  )

  linear_predictor_draws <-
    beta_draws %*%
    t(X_test)

  probability_draws <-
    plogis(
      linear_predictor_draws
    )

  probability_draws <-
    clip_prob(
      probability_draws
    )

  diagnostics <- list(
    Optimizer_convergence_code =
      as.integer(
        optimization$convergence
      ),
    Optimizer_objective =
      as.numeric(
        optimization$value
      ),
    Max_absolute_gradient =
      as.numeric(
        maximum_absolute_gradient
      ),
    Coefficient_L2_norm =
      sqrt(
        sum(
          beta_map[-1L]^2
        )
      ),
    Intercept_MAP =
      as.numeric(
        beta_map[1L]
      ),
    Minimum_Hessian_eigenvalue =
      min(
        hessian_eigenvalues
      ),
    Maximum_Hessian_eigenvalue =
      max(
        hessian_eigenvalues
      ),
    Hessian_condition_number =
      max(
        hessian_eigenvalues
      ) /
      min(
        hessian_eigenvalues
      ),
    Training_prevalence =
      mean(y_train),
    Prior_center_prevalence =
      smoothed_prevalence
  )

  list(
    draws =
      probability_draws,
    diagnostics =
      diagnostics
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
    # t(U) %*% U = precision.  Therefore
    # solve(U, z) has covariance precision^{-1}.
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
# 9. Exact corrected clinical BART
# ------------------------------------------------

predict_bart <- function(
  x_train,
  y_train,
  x_test,
  n_draws =
    N_POST_DRAWS_BART
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  if (
    length(unique(y_train)) <
      2L
  ) {
    stop(
      "BART training labels contain only one class."
    )
  }

  fit <- BART::pbart(
    x.train =
      x_train,
    y.train =
      y_train,
    x.test =
      x_test,
    ntree =
      BART_NTREE,
    ndpost =
      n_draws,
    nskip =
      N_BURN_BART,
    nkeeptest =
      n_draws,
    printevery =
      100000L
  )

  if (
    is.null(
      fit$yhat.test
    )
  ) {
    stop(
      "BART did not return yhat.test."
    )
  }

  latent_draws <- as.matrix(
    fit$yhat.test
  )

  observed_dimension <-
    as.integer(
      dim(
        latent_draws
      )
    )

  expected_dimension <- c(
    as.integer(n_draws),
    as.integer(nrow(x_test))
  )

  reversed_dimension <- rev(
    expected_dimension
  )

  orientation_action <-
    "accepted_documented_draw_by_patient_orientation"

  # BART::pbart returns yhat.test in posterior-draw by test-patient
  # orientation. The documented orientation must be checked first.
  # This is essential when n_draws = nrow(x_test) = 300 because both
  # orientations have dimensions 300 x 300. In that square case, the
  # matrix is accepted as returned and is never transposed.
  if (
    identical(
      observed_dimension,
      expected_dimension
    )
  ) {
    latent_draws <- latent_draws
  } else if (
    !identical(
      expected_dimension,
      reversed_dimension
    ) &&
      identical(
        observed_dimension,
        reversed_dimension
      )
  ) {
    latent_draws <- t(
      latent_draws
    )

    orientation_action <-
      "transposed_unambiguously_reversed_non_square_matrix"
  } else {
    stop(
      "BART returned yhat.test dimensions ",
      paste(
        observed_dimension,
        collapse = " x "
      ),
      "; expected ",
      paste(
        expected_dimension,
        collapse = " x "
      ),
      " in posterior-draw by test-patient orientation."
    )
  }

  probability_draws <-
    stats::pnorm(
      latent_draws
    )

  probability_draws <-
    clip_prob(
      probability_draws
    )

  probability_draws <-
    validate_prediction_draws(
      draws =
        probability_draws,
      expected_draws =
        n_draws,
      expected_test_size =
        nrow(x_test),
      model_name =
        "BART"
    )

  list(
    draws =
      probability_draws,
    diagnostics = list(
      BART_binary_offset =
        if (
          length(
            fit$binaryOffset
          ) == 0L
        ) {
          NA_real_
        } else {
          as.numeric(
            fit$binaryOffset[[1L]]
          )
        },
      BART_number_of_trees =
        as.integer(
          BART_NTREE
        ),
      BART_draw_rows =
        as.integer(
          nrow(
            probability_draws
          )
        ),
      BART_test_columns =
        as.integer(
          ncol(
            probability_draws
          )
        ),
      BART_orientation_action =
        orientation_action,
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
  switch(
    model_name,
    BLR = predict_blr(
      x_train,
      y_train,
      x_test
    ),
    BLDA = predict_blda(
      x_train,
      y_train,
      x_test
    ),
    BART = predict_bart(
      x_train,
      y_train,
      x_test
    ),
    stop(
      "Unknown model: ",
      model_name
    )
  )
}

expected_draw_count <- function(model_name) {
  switch(
    model_name,
    BLR =
      N_POST_DRAWS_BLR,
    BLDA =
      N_POST_DRAWS_BLDA,
    BART =
      N_POST_DRAWS_BART,
    stop(
      "Unknown model: ",
      model_name
    )
  )
}

normalize_seed <- function(raw_seed) {
  # R's integer seed range is limited to 2^31 - 1. Calculating the
  # original seed formula with integer literals caused overflow for
  # n = 500 and n = 1000. We calculate with doubles and reduce only
  # values outside the safe range. This preserves the exact seeds
  # already used for the completed n = 100 analysis.
  modulus <- 2147483646

  raw_seed <- as.double(raw_seed)

  if (
    length(raw_seed) != 1L ||
    is.na(raw_seed) ||
    !is.finite(raw_seed)
  ) {
    stop(
      "The raw seed is not finite."
    )
  }

  normalized_seed <- (
    (raw_seed - 1) %% modulus
  ) + 1

  seed_value <- as.integer(
    floor(normalized_seed)
  )

  if (
    is.na(seed_value) ||
    !is.finite(seed_value) ||
    seed_value <= 0L
  ) {
    stop(
      "Unable to construct a valid deterministic random seed."
    )
  }

  seed_value
}

make_fit_seed <- function(
  n_generated,
  scenario_id,
  rep_id,
  model_name
) {
  model_code <- match(
    model_name,
    MODEL_SET
  )

  raw_seed <- (
    as.double(BASE_SEED) +
      as.double(n_generated) * 10000000 +
      as.double(scenario_id) * 10000 +
      as.double(rep_id) * 10 +
      as.double(model_code)
  )

  normalize_seed(raw_seed)
}

make_split_seed <- function(
  n_generated,
  scenario_id,
  rep_id
) {
  raw_seed <- (
    as.double(BASE_SEED) +
      as.double(n_generated) * 100000 +
      as.double(scenario_id) * 100 +
      as.double(rep_id)
  )

  normalize_seed(raw_seed)
}

# ------------------------------------------------
# 10. Online curve-summary accumulators
# ------------------------------------------------

update_online_accumulator <- function(
  accumulator,
  new_data,
  key_columns,
  value_column
) {
  if (
    is.null(new_data) ||
    nrow(new_data) == 0L
  ) {
    return(accumulator)
  }

  chunk <- new_data |>
    dplyr::filter(
      is.finite(
        .data[[value_column]]
      )
    ) |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          key_columns
        )
      )
    ) |>
    dplyr::summarise(
      Count = dplyr::n(),
      Sum = sum(
        .data[[value_column]]
      ),
      SumSq = sum(
        .data[[value_column]]^2
      ),
      .groups = "drop"
    )

  if (
    is.null(accumulator) ||
    nrow(accumulator) == 0L
  ) {
    return(chunk)
  }

  dplyr::bind_rows(
    accumulator,
    chunk
  ) |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          key_columns
        )
      )
    ) |>
    dplyr::summarise(
      Count = sum(Count),
      Sum = sum(Sum),
      SumSq = sum(SumSq),
      .groups = "drop"
    )
}

finalize_online_accumulator <- function(
  accumulator,
  mean_name,
  sd_name
) {
  if (
    is.null(accumulator) ||
    nrow(accumulator) == 0L
  ) {
    return(tibble::tibble())
  }

  output <- accumulator |>
    dplyr::mutate(
      Mean_Value =
        Sum / Count,
      Variance_Value =
        dplyr::if_else(
          Count > 1L,
          pmax(
            (
              SumSq -
                Sum^2 /
                Count
            ) /
            (
              Count - 1L
            ),
            0
          ),
          NA_real_
        ),
      SD_Value =
        sqrt(
          Variance_Value
        )
    ) |>
    dplyr::select(
      -Sum,
      -SumSq,
      -Variance_Value
    )

  names(output)[
    names(output) ==
      "Mean_Value"
  ] <- mean_name

  names(output)[
    names(output) ==
      "SD_Value"
  ] <- sd_name

  output
}

design_key_columns <- c(
  "n_generated",
  "p",
  "prevalence_target",
  "boundary",
  "covariance",
  "label_noise",
  "missing_rate_target",
  "calibration_source"
)

dca_key_columns <- c(
  design_key_columns,
  "Strategy",
  "Threshold"
)

risk_key_columns <- c(
  design_key_columns,
  "Model",
  "Target_Coverage"
)

# ------------------------------------------------
# 11. Analyze one scenario-replication
# ------------------------------------------------

analyze_one_unit <- function(
  unit_data,
  unit_metadata
) {
  n_generated <-
    as.integer(
      unit_metadata$n_generated[[1L]]
    )

  scenario_id <-
    as.integer(
      unit_metadata$scenario_id[[1L]]
    )

  rep_id <-
    as.integer(
      unit_metadata$rep_id[[1L]]
    )

  predictor_count <-
    as.integer(
      unit_metadata$p[[1L]]
    )

  feature_columns <- paste0(
    "X",
    seq_len(
      predictor_count
    )
  )

  missing_features <- setdiff(
    feature_columns,
    names(unit_data)
  )

  if (
    length(missing_features) >
      0L
  ) {
    stop(
      "Missing feature columns: ",
      paste(
        missing_features,
        collapse = ", "
      )
    )
  }

  if (
    !("Simulation_Row_ID" %in%
      names(unit_data))
  ) {
    unit_data$Simulation_Row_ID <-
      seq_len(
        nrow(unit_data)
      )
  }

  split_seed <- make_split_seed(
    n_generated,
    scenario_id,
    rep_id
  )

  if (
    is.na(split_seed) ||
    !is.finite(split_seed) ||
    split_seed <= 0L
  ) {
    stop(
      "Invalid split seed for n=",
      n_generated,
      ", scenario=",
      scenario_id,
      ", replication=",
      rep_id,
      ": ",
      split_seed,
      "."
    )
  }

  split_object <- stratified_split(
    data = unit_data,
    stratify_column = "Y_true",
    train_prop = TRAIN_PROP,
    seed = split_seed
  )

  training_data <- unit_data[
    split_object$train_index,
    ,
    drop = FALSE
  ]

  test_data <- unit_data[
    split_object$test_index,
    ,
    drop = FALSE
  ]

  preprocessing <- fit_preprocess(
    training_data,
    feature_columns
  )

  x_train <- apply_preprocess(
    training_data,
    feature_columns,
    preprocessing
  )

  x_test <- apply_preprocess(
    test_data,
    feature_columns,
    preprocessing
  )

  y_train_observed <- as.integer(
    training_data$Y_observed
  )

  y_test_true <- as.integer(
    test_data$Y_true
  )

  if (
    length(
      unique(
        y_train_observed
      )
    ) < 2L
  ) {
    stop(
      "Observed training labels contain only one class."
    )
  }

  unit_id <- paste(
    n_generated,
    scenario_id,
    rep_id,
    sep = "_"
  )

  metadata_row <- tibble::tibble(
    Unit_ID = unit_id,
    n_generated =
      n_generated,
    scenario_id =
      scenario_id,
    scenario_name =
      as.character(
        unit_metadata$scenario_name[[1L]]
      ),
    rep_id =
      rep_id,
    p =
      predictor_count,
    prevalence_target =
      as.numeric(
        unit_metadata$prevalence_target[[1L]]
      ),
    boundary =
      as.character(
        unit_metadata$boundary[[1L]]
      ),
    covariance =
      as.character(
        unit_metadata$covariance[[1L]]
      ),
    label_noise =
      as.numeric(
        unit_metadata$label_noise[[1L]]
      ),
    missing_rate_target =
      as.numeric(
        unit_metadata$missing_rate_target[[1L]]
      ),
    calibration_source =
      as.character(
        unit_metadata$calibration_source[[1L]]
      ),
    n_total_unit =
      nrow(unit_data),
    n_train_internal =
      nrow(training_data),
    n_test_internal =
      nrow(test_data),
    training_observed_prevalence =
      mean(
        y_train_observed
      ),
    training_true_prevalence =
      mean(
        as.integer(
          training_data$Y_true
        )
      ),
    test_true_prevalence =
      mean(
        y_test_true
      ),
    split_seed =
      split_seed
  )

  sample_size <- if (
    SAVE_PREDICTION_SAMPLE
  ) {
    min(
      MAX_PREDICTIONS_PER_UNIT,
      nrow(test_data)
    )
  } else {
    0L
  }

  if (sample_size > 0L) {
    set.seed(
      split_seed + 999L
    )

    sampled_test_positions <- sample(
      seq_len(
        nrow(test_data)
      ),
      size = sample_size,
      replace = FALSE
    )
  } else {
    sampled_test_positions <-
      integer(0)
  }

  metric_rows <- list()
  status_rows <- list()
  prediction_sample_rows <- list()
  risk_rows <- list()
  successful_prediction_means <- list()

  unit_problem <- FALSE

  for (
    model_name in MODEL_SET
  ) {
    fit_seed <- make_fit_seed(
      n_generated,
      scenario_id,
      rep_id,
      model_name
    )

    if (
      is.na(fit_seed) ||
      !is.finite(fit_seed) ||
      fit_seed <= 0L
    ) {
      stop(
        "Invalid model seed for unit ",
        unit_id,
        ", model ",
        model_name,
        ": ",
        fit_seed,
        "."
      )
    }

    set.seed(fit_seed)

    start_time <- Sys.time()

    fit_result <- tryCatch(
      fit_model(
        model_name,
        x_train,
        y_train_observed,
        x_test
      ),
      error = function(e) {
        structure(
          list(
            message = e$message
          ),
          class =
            "model_fit_error"
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

    if (
      inherits(
        fit_result,
        "model_fit_error"
      )
    ) {
      unit_problem <- TRUE

      status_rows[[
        length(status_rows) + 1L
      ]] <- dplyr::bind_cols(
        metadata_row,
        tibble::tibble(
          Model = model_name,
          Fit_seed = fit_seed,
          Status = "FAILED",
          Failure_message =
            fit_result$message,
          Runtime_seconds =
            runtime_seconds
        )
      )

      next
    }

    posterior_draws <-
      tryCatch(
        validate_prediction_draws(
          draws =
            fit_result$draws,
          expected_draws =
            expected_draw_count(
              model_name
            ),
          expected_test_size =
            nrow(x_test),
          model_name =
            model_name
        ),
        error = function(e) {
          structure(
            list(
              message = e$message
            ),
            class =
              "model_fit_error"
          )
        }
      )

    if (
      inherits(
        posterior_draws,
        "model_fit_error"
      )
    ) {
      unit_problem <- TRUE

      status_rows[[
        length(status_rows) + 1L
      ]] <- dplyr::bind_cols(
        metadata_row,
        tibble::tibble(
          Model = model_name,
          Fit_seed = fit_seed,
          Status = "FAILED",
          Failure_message =
            posterior_draws$message,
          Runtime_seconds =
            runtime_seconds
        )
      )

      next
    }

    posterior_summary <-
      summarise_draws(
        posterior_draws
      )

    prediction_sd_across_patients <-
      if (
        nrow(
          posterior_summary
        ) >= 2L
      ) {
        stats::sd(
          posterior_summary$p_mean
        )
      } else {
        NA_real_
      }

    number_unique_predictions <-
      dplyr::n_distinct(
        round(
          posterior_summary$p_mean,
          digits = 12
        )
      )

    degenerate_prediction <- (
      nrow(
        posterior_summary
      ) > 1L &&
        (
          !is.finite(
            prediction_sd_across_patients
          ) ||
            prediction_sd_across_patients <=
              DEGENERATE_PREDICTION_TOL ||
            number_unique_predictions <=
              1L
        )
    )

    if (degenerate_prediction) {
      unit_problem <- TRUE
    }

    model_status <- if (
      degenerate_prediction
    ) {
      "DEGENERATE"
    } else {
      "OK"
    }

    status_base <- dplyr::bind_cols(
      metadata_row,
      tibble::tibble(
        Model = model_name,
        Fit_seed = fit_seed,
        Status = model_status,
        Failure_message =
          NA_character_,
        Runtime_seconds =
          runtime_seconds,
        Posterior_draw_rows =
          nrow(
            posterior_draws
          ),
        Posterior_draw_columns =
          ncol(
            posterior_draws
          ),
        Probability_min =
          min(
            posterior_summary$p_mean
          ),
        Probability_max =
          max(
            posterior_summary$p_mean
          ),
        Probability_mean =
          mean(
            posterior_summary$p_mean
          ),
        Probability_sd_across_patients =
          prediction_sd_across_patients,
        Number_unique_predictions =
          number_unique_predictions,
        Number_predicted_positive =
          sum(
            posterior_summary$p_mean >=
              CLASSIFICATION_THRESHOLD
          ),
        Degenerate_prediction =
          degenerate_prediction
      )
    )

    diagnostic_row <- tibble::as_tibble(
      fit_result$diagnostics
    )

    status_rows[[
      length(status_rows) + 1L
    ]] <- dplyr::bind_cols(
      status_base,
      diagnostic_row
    )

    if (degenerate_prediction) {
      next
    }

    metric_row <- calc_metrics(
      y =
        y_test_true,
      p_mean =
        posterior_summary$p_mean,
      p_sd =
        posterior_summary$p_sd,
      p_lwr =
        posterior_summary$p_lwr,
      p_upr =
        posterior_summary$p_upr
    )

    metric_rows[[
      length(metric_rows) + 1L
    ]] <- dplyr::bind_cols(
      metadata_row,
      tibble::tibble(
        Model =
          model_name
      ),
      metric_row
    )

    successful_prediction_means[[model_name]] <- posterior_summary$p_mean

    risk_row <- risk_coverage_table(
      y =
        y_test_true,
      p_mean =
        posterior_summary$p_mean,
      p_sd =
        posterior_summary$p_sd
    ) |>
      dplyr::mutate(
        Model =
          model_name,
        .before = 1
      )

    risk_rows[[
      length(risk_rows) + 1L
    ]] <- dplyr::bind_cols(
      metadata_row[
        rep(
          1L,
          nrow(risk_row)
        ),
        ,
        drop = FALSE
      ],
      risk_row
    )

    if (
      sample_size > 0L
    ) {
      prediction_sample_rows[[
        length(
          prediction_sample_rows
        ) + 1L
      ]] <- dplyr::bind_cols(
        metadata_row[
          rep(
            1L,
            sample_size
          ),
          ,
          drop = FALSE
        ],
        tibble::tibble(
          Model = model_name,
          Simulation_Row_ID =
            test_data$
            Simulation_Row_ID[
              sampled_test_positions
            ],
          Y_true =
            y_test_true[
              sampled_test_positions
            ],
          Y_observed =
            as.integer(
              test_data$
              Y_observed[
                sampled_test_positions
              ]
            ),
          Outcome_true = factor(
            y_test_true[
              sampled_test_positions
            ],
            levels = c(0, 1),
            labels = c(
              "No_appendicitis",
              "Appendicitis"
            )
          ),
          p_mean =
            posterior_summary$p_mean[
              sampled_test_positions
            ],
          p_sd =
            posterior_summary$p_sd[
              sampled_test_positions
            ],
          p_lwr =
            posterior_summary$p_lwr[
              sampled_test_positions
            ],
          p_upr =
            posterior_summary$p_upr[
              sampled_test_positions
            ]
        )
      )
    }

    rm(
      fit_result,
      posterior_draws,
      posterior_summary,
      status_base,
      diagnostic_row
    )

    gc(FALSE)
  }

  if (
    length(
      successful_prediction_means
    ) > 0L
  ) {
    dca_rows <-
      decision_curve_table_multi(
        prediction_list =
          successful_prediction_means,
        y =
          y_test_true
      )

    dca_rows <- dplyr::bind_cols(
      metadata_row[
        rep(
          1L,
          nrow(dca_rows)
        ),
        ,
        drop = FALSE
      ],
      dca_rows
    )
  } else {
    dca_rows <- tibble::tibble()
  }

  status_data <- dplyr::bind_rows(
    status_rows
  )

  all_requested_models_ok <- (
    nrow(status_data) ==
      length(MODEL_SET) &&
      all(
        status_data$Status ==
          "OK"
      )
  )

  # Curve summaries are intended to compare all three methods on the
  # same scenario-replication units.  Therefore an incomplete unit
  # contributes no DCA or risk-coverage rows.  If it is successfully
  # rerun later, the completed rerun contributes normally.
  risk_data <- dplyr::bind_rows(
    risk_rows
  )

  if (!all_requested_models_ok) {
    dca_rows <- tibble::tibble()
    risk_data <- tibble::tibble()
  }

  unit_status <- dplyr::bind_cols(
    metadata_row,
    tibble::tibble(
      Unit_Status =
        if (
          all_requested_models_ok
        ) {
          "COMPLETED"
        } else {
          "COMPLETED_WITH_PROBLEMS"
        },
      Requested_models =
        length(MODEL_SET),
      Successful_models =
        sum(
          status_data$Status ==
            "OK"
        ),
      Failed_models =
        sum(
          status_data$Status ==
            "FAILED"
        ),
      Degenerate_models =
        sum(
          status_data$Status ==
            "DEGENERATE"
        )
    )
  )

  list(
    metrics =
      dplyr::bind_rows(
        metric_rows
      ),
    status =
      status_data,
    unit_status =
      unit_status,
    prediction_sample =
      dplyr::bind_rows(
        prediction_sample_rows
      ),
    risk_curve =
      risk_data,
    dca_curve =
      dca_rows,
    strict_problem =
      unit_problem
  )
}


# ------------------------------------------------
# 12. Analyze one generated sample-size file
# ------------------------------------------------

checkpoint_chunk_files <- function(n_generated) {
  sample_checkpoint_dir <- file.path(
    dirs$checkpoints,
    paste0("n_", n_generated)
  )

  dir.create(
    sample_checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  sort(
    list.files(
      sample_checkpoint_dir,
      pattern = "^chunk_[0-9]{6}\\.rds$",
      full.names = TRUE
    )
  )
}

atomic_save_rds <- function(
  object,
  path,
  compress = "xz"
) {
  temp_path <- paste0(
    path,
    ".tmp"
  )

  saveRDS(
    object,
    temp_path,
    compress = compress
  )

  if (file.exists(path)) {
    unlink(
      path,
      force = TRUE
    )
  }

  moved <- file.rename(
    temp_path,
    path
  )

  if (!moved) {
    unlink(
      temp_path,
      force = TRUE
    )
    stop(
      "Unable to atomically save RDS checkpoint: ",
      path
    )
  }

  invisible(path)
}

validate_checkpoint_version <- function(
  checkpoint_object,
  path
) {
  object_version <- checkpoint_object$analysis_version

  if (
    is.null(object_version) ||
      !identical(
        as.character(object_version),
        ANALYSIS_VERSION
      )
  ) {
    stop(
      "Checkpoint version mismatch at:\n",
      path,
      "\nExpected: ",
      ANALYSIS_VERSION,
      "\nFound: ",
      ifelse(
        is.null(object_version),
        "<missing>",
        as.character(object_version)
      ),
      "\nDelete the current analysis output folder or set CLEAN_OUTPUT_ROOT=TRUE."
    )
  }

  invisible(TRUE)
}

deduplicate_keep_latest <- function(
  data,
  key_columns
) {
  if (
    is.null(data) ||
      nrow(data) == 0L
  ) {
    return(
      tibble::tibble()
    )
  }

  data <- as.data.frame(data)

  duplicate_index <- duplicated(
    data[
      ,
      key_columns,
      drop = FALSE
    ],
    fromLast = TRUE
  )

  data[
    !duplicate_index,
    ,
    drop = FALSE
  ]
}

analyze_sample_size_file <- function(
  n_generated,
  input_file
) {
  cat("\n================================================\n")
  cat(
    "Analyzing generated n =",
    n_generated,
    "\n"
  )
  cat(
    "Input:",
    input_file,
    "\n"
  )
  cat("================================================\n")

  n_generated_value <- as.integer(
    n_generated
  )

  sample_checkpoint_dir <- file.path(
    dirs$checkpoints,
    paste0(
      "n_",
      n_generated_value
    )
  )

  dir.create(
    sample_checkpoint_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  input_audit_file <- file.path(
    sample_checkpoint_dir,
    "input_unit_size_audit.rds"
  )

  completion_audit_file <- file.path(
    sample_checkpoint_dir,
    "sample_size_completion_audit.rds"
  )

  input_header <- names(
    data.table::fread(
      input_file,
      nrows = 0L,
      showProgress = FALSE
    )
  )

  feature_columns_in_file <- grep(
    "^X[0-9]+$",
    input_header,
    value = TRUE
  )

  missing_required_header_columns <- setdiff(
    required_metadata_columns,
    input_header
  )

  if (
    length(
      missing_required_header_columns
    ) > 0L
  ) {
    stop(
      "Input file is missing required columns: ",
      paste(
        missing_required_header_columns,
        collapse = ", "
      )
    )
  }

  selected_columns <- unique(
    c(
      required_metadata_columns,
      intersect(
        "row_id_within_scenario",
        input_header
      ),
      feature_columns_in_file
    )
  )

  data <- data.table::fread(
    input_file,
    select = selected_columns,
    data.table = TRUE,
    showProgress = TRUE
  )

  validate_input_schema(data)

  if (
    "row_id_within_scenario" %in%
      names(data)
  ) {
    data[
      ,
      Simulation_Row_ID :=
        row_id_within_scenario
    ]
  } else {
    data[
      ,
      Simulation_Row_ID :=
        seq_len(.N),
      by = .(
        scenario_id,
        rep_id
      )
    ]
  }

  unit_columns <- c(
    "scenario_id",
    "scenario_name",
    "rep_id",
    "n_train",
    "p",
    "prevalence_target",
    "boundary",
    "covariance",
    "label_noise",
    "missing_rate_target",
    "calibration_source"
  )

  units <- unique(
    data[
      ,
      ..unit_columns
    ]
  )

  data.table::setorder(
    units,
    scenario_id,
    rep_id
  )

  units[
    ,
    n_generated :=
      n_generated_value
  ]

  units[
    ,
    Unit_ID :=
      paste(
        n_generated_value,
        scenario_id,
        rep_id,
        sep = "_"
      )
  ]

  data.table::setcolorder(
    units,
    c(
      "Unit_ID",
      "n_generated",
      setdiff(
        names(units),
        c(
          "Unit_ID",
          "n_generated"
        )
      )
    )
  )

  unit_size_audit <- data[
    ,
    .(
      Rows_in_unit = .N
    ),
    by = .(
      scenario_id,
      rep_id
    )
  ]

  unit_size_audit[
    ,
    `:=`(
      n_generated =
        n_generated_value,
      Expected_rows =
        n_generated_value,
      Unit_size_correct =
        Rows_in_unit ==
        n_generated_value
    )
  ]

  atomic_save_rds(
    list(
      analysis_version =
        ANALYSIS_VERSION,
      audit =
        as.data.frame(
          unit_size_audit
        )
    ),
    input_audit_file
  )

  if (
    any(
      !unit_size_audit$
      Unit_size_correct
    )
  ) {
    stop(
      "At least one unit in generated_data_n",
      n_generated_value,
      ".csv does not contain exactly ",
      n_generated_value,
      " rows."
    )
  }

  data.table::setkeyv(
    data,
    c(
      "scenario_id",
      "rep_id"
    )
  )

  if (
    is.finite(
      MAX_UNITS_PER_N
    )
  ) {
    units <- units |>
      dplyr::slice_head(
        n = min(
          as.integer(
            MAX_UNITS_PER_N
          ),
          nrow(units)
        )
      )
  }

  existing_chunk_files <-
    checkpoint_chunk_files(
      n_generated_value
    )

  completed_unit_ids <-
    character(0)

  if (
    RESUME_ANALYSIS &&
      length(
        existing_chunk_files
      ) > 0L
  ) {
    existing_chunk_objects <-
      lapply(
        existing_chunk_files,
        readRDS
      )

    invisible(
      mapply(
        validate_checkpoint_version,
        existing_chunk_objects,
        existing_chunk_files,
        SIMPLIFY = FALSE
      )
    )

    prior_unit_status <-
      dplyr::bind_rows(
        lapply(
          existing_chunk_objects,
          function(x) {
            x$unit_status
          }
        )
      ) |>
      deduplicate_keep_latest(
        "Unit_ID"
      )

    prior_metrics <-
      dplyr::bind_rows(
        lapply(
          existing_chunk_objects,
          function(x) {
            x$metrics
          }
        )
      ) |>
      deduplicate_keep_latest(
        c(
          "Unit_ID",
          "Model"
        )
      )

    prior_model_status <-
      dplyr::bind_rows(
        lapply(
          existing_chunk_objects,
          function(x) {
            x$status
          }
        )
      ) |>
      deduplicate_keep_latest(
        c(
          "Unit_ID",
          "Model"
        )
      )

    status_complete_ids <-
      prior_unit_status |>
      dplyr::filter(
        Unit_Status ==
          "COMPLETED"
      ) |>
      dplyr::pull(
        Unit_ID
      )

    metric_complete_ids <-
      prior_metrics |>
      dplyr::distinct(
        Unit_ID,
        Model
      ) |>
      dplyr::count(
        Unit_ID,
        name =
          "Model_count"
      ) |>
      dplyr::filter(
        Model_count ==
          length(
            MODEL_SET
          )
      ) |>
      dplyr::pull(
        Unit_ID
      )

    model_status_complete_ids <-
      prior_model_status |>
      dplyr::filter(
        Status ==
          "OK"
      ) |>
      dplyr::distinct(
        Unit_ID,
        Model
      ) |>
      dplyr::count(
        Unit_ID,
        name =
          "OK_model_count"
      ) |>
      dplyr::filter(
        OK_model_count ==
          length(
            MODEL_SET
          )
      ) |>
      dplyr::pull(
        Unit_ID
      )

    completed_unit_ids <- Reduce(
      intersect,
      list(
        status_complete_ids,
        metric_complete_ids,
        model_status_complete_ids
      )
    )

    cat(
      "Verified completed current analysis units:",
      length(
        completed_unit_ids
      ),
      "of",
      nrow(units),
      "\n"
    )
  }

  if (
    length(
      existing_chunk_files
    ) > 0L
  ) {
    existing_chunk_numbers <-
      as.integer(
        sub(
          "^chunk_([0-9]{6})\\.rds$",
          "\\1",
          basename(
            existing_chunk_files
          )
        )
      )

    next_chunk_number <-
      max(
        existing_chunk_numbers,
        na.rm = TRUE
      ) + 1L
  } else {
    next_chunk_number <- 1L
  }

  metric_buffer <- list()
  prediction_buffer <- list()
  status_buffer <- list()
  unit_status_buffer <- list()
  dca_buffer <- list()
  risk_buffer <- list()

  buffer_index <- 0L
  processed_since_checkpoint <- 0L

  flush_checkpoint <- function() {
    if (
      length(metric_buffer) == 0L &&
        length(unit_status_buffer) ==
        0L
    ) {
      return(
        invisible(NULL)
      )
    }

    metric_data <-
      dplyr::bind_rows(
        metric_buffer
      )

    prediction_data <-
      dplyr::bind_rows(
        prediction_buffer
      )

    status_data <-
      dplyr::bind_rows(
        status_buffer
      )

    unit_status_data <-
      dplyr::bind_rows(
        unit_status_buffer
      )

    dca_data <-
      dplyr::bind_rows(
        dca_buffer
      )

    risk_data <-
      dplyr::bind_rows(
        risk_buffer
      )

    dca_accumulator <-
      update_online_accumulator(
        accumulator =
          tibble::tibble(),
        new_data =
          dca_data,
        key_columns =
          dca_key_columns,
        value_column =
          "Net_Benefit"
      )

    risk_accumulator <-
      update_online_accumulator(
        accumulator =
          tibble::tibble(),
        new_data =
          risk_data,
        key_columns =
          risk_key_columns,
        value_column =
          "Selective_Risk"
      )

    chunk_path <- file.path(
      sample_checkpoint_dir,
      sprintf(
        "chunk_%06d.rds",
        next_chunk_number
      )
    )

    chunk_object <- list(
      analysis_version =
        ANALYSIS_VERSION,
      n_generated =
        n_generated_value,
      saved_at =
        as.character(
          Sys.time()
        ),
      metrics =
        metric_data,
      prediction_sample =
        prediction_data,
      status =
        status_data,
      unit_status =
        unit_status_data,
      dca_accumulator =
        dca_accumulator,
      risk_accumulator =
        risk_accumulator
    )

    atomic_save_rds(
      chunk_object,
      chunk_path
    )

    next_chunk_number <<-
      next_chunk_number + 1L

    metric_buffer <<- list()
    prediction_buffer <<- list()
    status_buffer <<- list()
    unit_status_buffer <<- list()
    dca_buffer <<- list()
    risk_buffer <<- list()

    gc(FALSE)

    invisible(NULL)
  }

  for (
    unit_position in
    seq_len(
      nrow(units)
    )
  ) {
    unit_metadata <-
      units[
        unit_position,
        ,
        drop = FALSE
      ]

    unit_id <-
      unit_metadata$
      Unit_ID[[1L]]

    if (
      unit_id %in%
        completed_unit_ids
    ) {
      next
    }

    scenario_value <-
      unit_metadata$
      scenario_id[[1L]]

    rep_value <-
      unit_metadata$
      rep_id[[1L]]

    unit_data <- data[
      list(
        scenario_value,
        rep_value
      ),
      nomatch = 0L
    ]

    data.table::setorder(
      unit_data,
      Simulation_Row_ID
    )

    unit_data <- as.data.frame(
      unit_data
    )

    if (
      nrow(unit_data) !=
        n_generated_value
    ) {
      stop(
        "Unit ",
        unit_id,
        " contains ",
        nrow(unit_data),
        " rows; expected ",
        n_generated_value,
        "."
      )
    }

    buffer_index <-
      buffer_index + 1L

    processed_since_checkpoint <-
      processed_since_checkpoint + 1L

    result <- tryCatch(
      analyze_one_unit(
        unit_data =
          unit_data,
        unit_metadata =
          unit_metadata
      ),
      error = function(e) {
        metadata_row <-
          tibble::tibble(
            Unit_ID =
              unit_id,
            n_generated =
              n_generated_value,
            scenario_id =
              as.integer(
                scenario_value
              ),
            scenario_name =
              as.character(
                unit_metadata$
                scenario_name[[1L]]
              ),
            rep_id =
              as.integer(
                rep_value
              ),
            p =
              as.integer(
                unit_metadata$
                p[[1L]]
              ),
            prevalence_target =
              as.numeric(
                unit_metadata$
                prevalence_target[[1L]]
              ),
            boundary =
              as.character(
                unit_metadata$
                boundary[[1L]]
              ),
            covariance =
              as.character(
                unit_metadata$
                covariance[[1L]]
              ),
            label_noise =
              as.numeric(
                unit_metadata$
                label_noise[[1L]]
              ),
            missing_rate_target =
              as.numeric(
                unit_metadata$
                missing_rate_target[[1L]]
              ),
            calibration_source =
              as.character(
                unit_metadata$
                calibration_source[[1L]]
              ),
            Unit_Status =
              "UNIT_ERROR",
            Requested_models =
              length(
                MODEL_SET
              ),
            Successful_models =
              0L,
            Failed_models =
              length(
                MODEL_SET
              ),
            Degenerate_models =
              0L,
            Failure_message =
              e$message
          )

        list(
          metrics =
            tibble::tibble(),
          status =
            tibble::tibble(),
          unit_status =
            metadata_row,
          prediction_sample =
            tibble::tibble(),
          risk_curve =
            tibble::tibble(),
          dca_curve =
            tibble::tibble(),
          strict_problem =
            TRUE
        )
      }
    )

    metric_buffer[[buffer_index]] <-
      result$metrics

    prediction_buffer[[buffer_index]] <-
      result$prediction_sample

    status_buffer[[buffer_index]] <-
      result$status

    unit_status_buffer[[buffer_index]] <-
      result$unit_status

    dca_buffer[[buffer_index]] <-
      result$dca_curve

    risk_buffer[[buffer_index]] <-
      result$risk_curve

    if (
      unit_position %% 25L ==
        0L
    ) {
      cat(
        "n =",
        n_generated_value,
        "| unit",
        unit_position,
        "of",
        nrow(units),
        "\n"
      )
    }

    must_stop <- (
      result$strict_problem &&
        (
          STOP_ON_MODEL_FAILURE ||
            STOP_ON_DEGENERATE_PREDICTIONS
        )
    )

    if (
      processed_since_checkpoint >=
        CHECKPOINT_EVERY_UNITS ||
        must_stop
    ) {
      flush_checkpoint()
      processed_since_checkpoint <- 0L
      buffer_index <- 0L
    }

    if (must_stop) {
      stop(
        "Strict quality control stopped after unit ",
        unit_id,
        ". Inspect RDS checkpoints under ",
        sample_checkpoint_dir,
        "."
      )
    }

    rm(
      unit_data,
      result
    )

    gc(FALSE)
  }

  if (
    length(metric_buffer) > 0L ||
      length(unit_status_buffer) >
      0L
  ) {
    flush_checkpoint()
  }

  all_chunk_files <-
    checkpoint_chunk_files(
      n_generated_value
    )

  all_chunk_objects <-
    lapply(
      all_chunk_files,
      readRDS
    )

  invisible(
    mapply(
      validate_checkpoint_version,
      all_chunk_objects,
      all_chunk_files,
      SIMPLIFY = FALSE
    )
  )

  latest_unit_status <-
    dplyr::bind_rows(
      lapply(
        all_chunk_objects,
        function(x) {
          x$unit_status
        }
      )
    ) |>
    deduplicate_keep_latest(
      "Unit_ID"
    )

  completed_count <-
    sum(
      latest_unit_status$
      Unit_Status ==
        "COMPLETED",
      na.rm = TRUE
    )

  problem_count <-
    sum(
      latest_unit_status$
      Unit_Status !=
        "COMPLETED",
      na.rm = TRUE
    )

  completion_audit <-
    tibble::tibble(
      n_generated =
        n_generated_value,
      Expected_units =
        nrow(units),
      Completed_units =
        completed_count,
      Problem_units =
        problem_count,
      Pending_units =
        max(
          nrow(units) -
            nrow(
              latest_unit_status
            ),
          0L
        ),
      Complete =
        completed_count ==
        nrow(units)
    )

  atomic_save_rds(
    list(
      analysis_version =
        ANALYSIS_VERSION,
      audit =
        completion_audit
    ),
    completion_audit_file
  )

  cat(
    "Finished generated n =",
    n_generated_value,
    "| completed units:",
    completed_count,
    "of",
    nrow(units),
    "| problem units:",
    problem_count,
    "\n"
  )

  rm(data)
  gc(FALSE)

  invisible(
    completion_audit
  )
}

# ------------------------------------------------
# 13. Run all generated sample sizes
# ------------------------------------------------

print(input_files)

sample_size_run_audits <- list()

for (
  file_index in
  seq_len(
    nrow(input_files)
  )
) {
  n_value <-
    input_files$
    n_generated[
      file_index
    ]

  sample_size_run_audits[[
    as.character(
      n_value
    )
  ]] <-
    analyze_sample_size_file(
      n_generated =
        n_value,
      input_file =
        input_files$
        file_path[
          file_index
        ]
    )
}

sample_size_run_audit <-
  dplyr::bind_rows(
    sample_size_run_audits
  )

# ------------------------------------------------
# 14. Combine completed checkpoint results
# ------------------------------------------------

all_chunk_files <- sort(
  list.files(
    dirs$checkpoints,
    pattern =
      "^chunk_[0-9]{6}\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
)

if (
  length(
    all_chunk_files
  ) == 0L
) {
  stop(
    "No current analysis checkpoint chunks were produced."
  )
}

all_chunk_objects <- lapply(
  all_chunk_files,
  readRDS
)

invisible(
  mapply(
    validate_checkpoint_version,
    all_chunk_objects,
    all_chunk_files,
    SIMPLIFY = FALSE
  )
)

all_metrics <- dplyr::bind_rows(
  lapply(
    all_chunk_objects,
    function(x) {
      x$metrics
    }
  )
) |>
  deduplicate_keep_latest(
    c(
      "Unit_ID",
      "Model"
    )
  )

all_status <- dplyr::bind_rows(
  lapply(
    all_chunk_objects,
    function(x) {
      x$status
    }
  )
) |>
  deduplicate_keep_latest(
    c(
      "Unit_ID",
      "Model"
    )
  )

all_unit_status <- dplyr::bind_rows(
  lapply(
    all_chunk_objects,
    function(x) {
      x$unit_status
    }
  )
) |>
  deduplicate_keep_latest(
    "Unit_ID"
  )

all_prediction_samples <-
  dplyr::bind_rows(
    lapply(
      all_chunk_objects,
      function(x) {
        x$prediction_sample
      }
    )
  )

if (
  nrow(
    all_prediction_samples
  ) > 0L
) {
  all_prediction_samples <-
    all_prediction_samples |>
    deduplicate_keep_latest(
      c(
        "Unit_ID",
        "Model",
        "Simulation_Row_ID"
      )
    )
}

complete_metric_unit_ids <-
  all_metrics |>
  dplyr::distinct(
    Unit_ID,
    Model
  ) |>
  dplyr::count(
    Unit_ID,
    name =
      "Model_count"
  ) |>
  dplyr::filter(
    Model_count ==
      length(
        MODEL_SET
      )
  ) |>
  dplyr::pull(
    Unit_ID
  )

complete_status_unit_ids <-
  all_status |>
  dplyr::filter(
    Status ==
      "OK"
  ) |>
  dplyr::distinct(
    Unit_ID,
    Model
  ) |>
  dplyr::count(
    Unit_ID,
    name =
      "OK_model_count"
  ) |>
  dplyr::filter(
    OK_model_count ==
      length(
        MODEL_SET
      )
  ) |>
  dplyr::pull(
    Unit_ID
  )

complete_unit_status_ids <-
  all_unit_status |>
  dplyr::filter(
    Unit_Status ==
      "COMPLETED"
  ) |>
  dplyr::pull(
    Unit_ID
  )

complete_unit_ids <- Reduce(
  intersect,
  list(
    complete_metric_unit_ids,
    complete_status_unit_ids,
    complete_unit_status_ids
  )
)

all_metrics <- all_metrics |>
  dplyr::filter(
    Unit_ID %in%
      complete_unit_ids
  )

if (
  nrow(
    all_prediction_samples
  ) > 0L
) {
  all_prediction_samples <-
    all_prediction_samples |>
    dplyr::filter(
      Unit_ID %in%
        complete_unit_ids
    )
}

if (
  nrow(all_metrics) == 0L
) {
  stop(
    "No complete three-model simulation units were available."
  )
}

# Input-unit audit objects are stored as RDS and included in the
# final audit CSV later.
input_audit_files <- sort(
  list.files(
    dirs$checkpoints,
    pattern =
      "^input_unit_size_audit\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
)

all_input_unit_audits <-
  dplyr::bind_rows(
    lapply(
      input_audit_files,
      function(path) {
        object <- readRDS(path)
        validate_checkpoint_version(
          object,
          path
        )
        object$audit
      }
    )
  )

# ------------------------------------------------
# 15. Finalize DCA and risk summaries
# ------------------------------------------------

combined_dca_accumulator <-
  tibble::tibble()

combined_risk_accumulator <-
  tibble::tibble()

for (
  chunk_object in
  all_chunk_objects
) {
  combined_dca_accumulator <-
    dplyr::bind_rows(
      combined_dca_accumulator,
      chunk_object$
      dca_accumulator
    ) |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          dca_key_columns
        )
      )
    ) |>
    dplyr::summarise(
      Count =
        sum(Count),
      Sum =
        sum(Sum),
      SumSq =
        sum(SumSq),
      .groups =
        "drop"
    )

  combined_risk_accumulator <-
    dplyr::bind_rows(
      combined_risk_accumulator,
      chunk_object$
      risk_accumulator
    ) |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          risk_key_columns
        )
      )
    ) |>
    dplyr::summarise(
      Count =
        sum(Count),
      Sum =
        sum(Sum),
      SumSq =
        sum(SumSq),
      .groups =
        "drop"
    )
}

dca_summary_by_design <-
  finalize_online_accumulator(
    combined_dca_accumulator,
    mean_name =
      "Mean_Net_Benefit",
    sd_name =
      "SD_Net_Benefit"
  )

risk_summary_by_design <-
  finalize_online_accumulator(
    combined_risk_accumulator,
    mean_name =
      "Mean_Selective_Risk",
    sd_name =
      "SD_Selective_Risk"
  )

dca_summary_by_n <-
  combined_dca_accumulator |>
  dplyr::group_by(
    n_generated,
    Strategy,
    Threshold
  ) |>
  dplyr::summarise(
    Count =
      sum(Count),
    Sum =
      sum(Sum),
    SumSq =
      sum(SumSq),
    .groups =
      "drop"
  ) |>
  finalize_online_accumulator(
    mean_name =
      "Mean_Net_Benefit",
    sd_name =
      "SD_Net_Benefit"
  )

risk_summary_by_n <-
  combined_risk_accumulator |>
  dplyr::group_by(
    n_generated,
    Model,
    Target_Coverage
  ) |>
  dplyr::summarise(
    Count =
      sum(Count),
    Sum =
      sum(Sum),
    SumSq =
      sum(SumSq),
    .groups =
      "drop"
  ) |>
  finalize_online_accumulator(
    mean_name =
      "Mean_Selective_Risk",
    sd_name =
      "SD_Selective_Risk"
  )

# ------------------------------------------------
# 16. Summary tables
# ------------------------------------------------

metric_columns <- c(
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
  "calibration_intercept",
  "calibration_slope",
  "ece",
  "mean_uncertainty",
  "mean_interval_width",
  "selective_coverage",
  "selective_risk"
)

summarise_metrics <- function(
  data,
  grouping_columns
) {
  data |>
    dplyr::group_by(
      dplyr::across(
        dplyr::all_of(
          grouping_columns
        )
      )
    ) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(
          metric_columns
        ),
        list(
          mean =
            ~ mean(
              .x,
              na.rm = TRUE
            ),
          sd =
            ~ stats::sd(
              .x,
              na.rm = TRUE
            ),
          median =
            ~ stats::median(
              .x,
              na.rm = TRUE
            ),
          q025 =
            ~ as.numeric(
              stats::quantile(
                .x,
                0.025,
                na.rm = TRUE,
                names = FALSE
              )
            ),
          q975 =
            ~ as.numeric(
              stats::quantile(
                .x,
                0.975,
                na.rm = TRUE,
                names = FALSE
              )
            )
        ),
        .names =
          "{.col}_{.fn}"
      ),
      n_results =
        dplyr::n(),
      .groups =
        "drop"
    )
}

overall_model_summary <-
  summarise_metrics(
    all_metrics,
    "Model"
  )

summary_by_sample_size <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "Model"
    )
  )

summary_by_dimension <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "p",
      "Model"
    )
  )

summary_by_boundary <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "boundary",
      "Model"
    )
  )

summary_by_design <-
  summarise_metrics(
    all_metrics,
    c(
      design_key_columns,
      "Model"
    )
  )

summary_by_boundary_covariance <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "boundary",
      "covariance",
      "Model"
    )
  )

summary_by_noise_missingness <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "label_noise",
      "missing_rate_target",
      "Model"
    )
  )

summary_by_calibration_source <-
  summarise_metrics(
    all_metrics,
    c(
      "n_generated",
      "calibration_source",
      "Model"
    )
  )

# Strict unique-best AUROC counts. Exact equality is treated as a tie,
# matching the manuscript definition.
auc_wide <- all_metrics |>
  dplyr::select(
    Unit_ID,
    n_generated,
    Model,
    auc
  ) |>
  tidyr::pivot_wider(
    names_from =
      Model,
    values_from =
      auc
  ) |>
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(
        MODEL_SET
      ),
      is.finite
    )
  )

classify_auc_winner <- function(
  blr,
  blda,
  bart
) {
  values <- c(
    BLR = blr,
    BLDA = blda,
    BART = bart
  )

  maximum_value <- max(
    values
  )

  winner_names <- names(
    values
  )[
    values ==
      maximum_value
  ]

  if (
    length(
      winner_names
    ) == 1L
  ) {
    winner_names
  } else {
    "Tied"
  }
}

if (
  nrow(
    auc_wide
  ) > 0L
) {
  auc_wide$Winner <- mapply(
    classify_auc_winner,
    auc_wide$BLR,
    auc_wide$BLDA,
    auc_wide$BART,
    USE.NAMES = FALSE
  )

  unique_best_auc_counts <-
    auc_wide |>
    dplyr::count(
      n_generated,
      Winner,
      name =
        "Count"
    ) |>
    tidyr::complete(
      n_generated =
        N_FILES,
      Winner =
        c(
          MODEL_SET,
          "Tied"
        ),
      fill =
        list(
          Count = 0L
        )
    ) |>
    dplyr::group_by(
      n_generated
    ) |>
    dplyr::mutate(
      Total =
        sum(Count)
    ) |>
    dplyr::ungroup()
} else {
  unique_best_auc_counts <-
    tibble::tibble()
}

# ------------------------------------------------
# 17. Status and completeness audit
# ------------------------------------------------

if (
  !(
    "Probability_sd_across_patients" %in%
      names(all_status)
  )
) {
  all_status$
  Probability_sd_across_patients <-
    NA_real_
}

if (
  !(
    "Degenerate_prediction" %in%
      names(all_status)
  )
) {
  all_status$
  Degenerate_prediction <-
    FALSE
}

model_status_audit <-
  all_status |>
  dplyr::group_by(
    n_generated,
    Model
  ) |>
  dplyr::summarise(
    Total_fits =
      dplyr::n(),
    OK_fits =
      sum(
        Status ==
          "OK",
        na.rm = TRUE
      ),
    Failed_fits =
      sum(
        Status ==
          "FAILED",
        na.rm = TRUE
      ),
    Degenerate_fits =
      sum(
        Status ==
          "DEGENERATE",
        na.rm = TRUE
      ),
    Minimum_prediction_SD =
      if (
        any(
          is.finite(
            Probability_sd_across_patients
          )
        )
      ) {
        min(
          Probability_sd_across_patients,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    Mean_prediction_SD =
      if (
        any(
          is.finite(
            Probability_sd_across_patients
          )
        )
      ) {
        mean(
          Probability_sd_across_patients,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    .groups =
      "drop"
  )

metric_completeness <-
  all_metrics |>
  dplyr::count(
    Unit_ID,
    name =
      "Observed_model_results"
  ) |>
  dplyr::mutate(
    Expected_model_results =
      length(
        MODEL_SET
      ),
    Complete =
      Observed_model_results ==
      Expected_model_results
  )

# ------------------------------------------------
# 18. Paired simulation model comparisons
# ------------------------------------------------

inferential_metrics <- c(
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

model_pairs <- combn(
  MODEL_SET,
  2L,
  simplify = FALSE
)

paired_comparison_rows <- list()
omnibus_rows <- list()

comparison_group_values <- c(
  "ALL",
  as.character(
    N_FILES
  )
)

for (
  group_value in
  comparison_group_values
) {
  if (
    group_value ==
      "ALL"
  ) {
    comparison_data <-
      all_metrics
  } else {
    comparison_data <-
      all_metrics |>
      dplyr::filter(
        n_generated ==
          as.integer(
            group_value
          )
      )
  }

  for (
    metric_name in
    inferential_metrics
  ) {
    wide_data <-
      comparison_data |>
      dplyr::select(
        Unit_ID,
        Model,
        dplyr::all_of(
          metric_name
        )
      ) |>
      tidyr::pivot_wider(
        names_from =
          Model,
        values_from =
          dplyr::all_of(
            metric_name
          )
      ) |>
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(
            MODEL_SET
          ),
          is.finite
        )
      )

    if (
      nrow(
        wide_data
      ) < 2L
    ) {
      next
    }

    long_data <-
      wide_data |>
      tidyr::pivot_longer(
        cols =
          dplyr::all_of(
            MODEL_SET
          ),
        names_to =
          "Model",
        values_to =
          "Value"
      )

    friedman_result <-
      tryCatch(
        stats::friedman.test(
          Value ~
            Model |
            Unit_ID,
          data =
            long_data
        ),
        error =
          function(e) NULL
      )

    if (
      !is.null(
        friedman_result
      )
    ) {
      omnibus_rows[[
        length(
          omnibus_rows
        ) + 1L
      ]] <-
        tibble::tibble(
          Sample_Size_Group =
            group_value,
          Metric =
            metric_name,
          Test =
            "Friedman test",
          Complete_units =
            nrow(
              wide_data
            ),
          Statistic =
            unname(
              friedman_result$
              statistic
            ),
          Degrees_of_freedom =
            unname(
              friedman_result$
              parameter
            ),
          p_value =
            unname(
              friedman_result$
              p.value
            )
        )
    }

    metric_pair_rows <- list()

    for (
      model_pair in
      model_pairs
    ) {
      model_1 <-
        model_pair[1L]

      model_2 <-
        model_pair[2L]

      difference <-
        wide_data[[
          model_1
        ]] -
        wide_data[[
          model_2
        ]]

      test_result <-
        tryCatch(
          suppressWarnings(
            stats::wilcox.test(
              wide_data[[
                model_1
              ]],
              wide_data[[
                model_2
              ]],
              paired = TRUE,
              exact = FALSE,
              correct = TRUE
            )
          ),
          error =
            function(e) NULL
        )

      metric_pair_rows[[
        length(
          metric_pair_rows
        ) + 1L
      ]] <-
        tibble::tibble(
          Sample_Size_Group =
            group_value,
          Metric =
            metric_name,
          Comparison =
            paste0(
              model_1,
              " - ",
              model_2
            ),
          Model_1 =
            model_1,
          Model_2 =
            model_2,
          Complete_pairs =
            length(
              difference
            ),
          Model_1_mean =
            mean(
              wide_data[[
                model_1
              ]]
            ),
          Model_2_mean =
            mean(
              wide_data[[
                model_2
              ]]
            ),
          Mean_difference =
            mean(
              difference
            ),
          Median_difference =
            stats::median(
              difference
            ),
          Difference_q025 =
            as.numeric(
              stats::quantile(
                difference,
                0.025,
                names = FALSE
              )
            ),
          Difference_q975 =
            as.numeric(
              stats::quantile(
                difference,
                0.975,
                names = FALSE
              )
            ),
          Statistic =
            if (
              is.null(
                test_result
              )
            ) {
              NA_real_
            } else {
              unname(
                test_result$
                statistic
              )
            },
          p_value_raw =
            if (
              is.null(
                test_result
              )
            ) {
              NA_real_
            } else {
              unname(
                test_result$
                p.value
              )
            }
        )
    }

    metric_pair_data <-
      dplyr::bind_rows(
        metric_pair_rows
      ) |>
      dplyr::mutate(
        p_value_adjusted =
          stats::p.adjust(
            p_value_raw,
            method =
              "holm"
          ),
        Significant_adjusted =
          p_value_adjusted <
          0.05
      )

    paired_comparison_rows[[
      length(
        paired_comparison_rows
      ) + 1L
    ]] <-
      metric_pair_data
  }
}

simulation_omnibus_comparisons <-
  dplyr::bind_rows(
    omnibus_rows
  )

simulation_pairwise_comparisons <-
  dplyr::bind_rows(
    paired_comparison_rows
  )

# ------------------------------------------------
# 19. Calibration bins from saved prediction samples
# ------------------------------------------------

if (
  nrow(
    all_prediction_samples
  ) > 0L
) {
  calibration_bins <-
    all_prediction_samples |>
    dplyr::mutate(
      Bin =
        cut(
          p_mean,
          breaks =
            seq(
              0,
              1,
              length.out =
                CALIBRATION_BINS +
                1L
            ),
          include.lowest =
            TRUE,
          labels =
            FALSE
        )
    ) |>
    dplyr::group_by(
      n_generated,
      Model,
      Bin
    ) |>
    dplyr::summarise(
      Mean_Predicted =
        mean(
          p_mean
        ),
      Observed_Rate =
        mean(
          Y_true
        ),
      n =
        dplyr::n(),
      .groups =
        "drop"
    )
} else {
  calibration_bins <-
    tibble::tibble()
}


# ------------------------------------------------
# 20. Figures
# ------------------------------------------------

plot_metric_box <- function(
  metric_name,
  y_label,
  title,
  file_name
) {
  plot_object <- ggplot(
    all_metrics,
    aes(
      x = factor(
        n_generated
      ),
      y =
        .data[[metric_name]],
      fill = Model
    )
  ) +
    geom_boxplot(
      alpha = 0.82,
      outlier.alpha = 0.12,
      position =
        position_dodge(
          width = 0.8
        )
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      aes(group = Model),
      position =
        position_dodge(
          width = 0.8
        ),
      shape = 23,
      size = 2.6,
      color = "black"
    ) +
    scale_fill_manual(
      values = model_colors
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "Clinical-method model implementations; ",
        "independent scenario-replication test sets"
      ),
      x = "Generated sample size",
      y = y_label,
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      file_name
    ),
    width = 11,
    height = 7
  )

  print(plot_object)
}

plot_metric_box(
  "auc",
  "AUC",
  "B-RISA Simulation Discrimination by Sample Size",
  "auc_by_sample_size"
)

plot_metric_box(
  "auprc",
  "AUPRC",
  "B-RISA Simulation Precision-Recall Performance",
  "auprc_by_sample_size"
)

plot_metric_box(
  "brier",
  "Brier score",
  "B-RISA Simulation Probability Error",
  "brier_by_sample_size"
)

plot_metric_box(
  "logloss",
  "Log loss",
  "B-RISA Simulation Logarithmic Loss",
  "logloss_by_sample_size"
)

plot_metric_box(
  "selective_risk",
  "Selective risk",
  "B-RISA Simulation Selective Prediction",
  "selective_risk_by_sample_size"
)

plot_metric_box(
  "mean_uncertainty",
  "Mean posterior predictive SD",
  "B-RISA Simulation Posterior Predictive Uncertainty",
  "uncertainty_by_sample_size"
)

plot_metric_box(
  "mean_interval_width",
  "Mean 95% interval width",
  "B-RISA Simulation Posterior Predictive Interval Width",
  "interval_width_by_sample_size"
)

classification_plot_data <- all_metrics |>
  dplyr::select(
    n_generated,
    Model,
    sensitivity,
    specificity
  ) |>
  tidyr::pivot_longer(
    cols = c(
      sensitivity,
      specificity
    ),
    names_to = "Metric",
    values_to = "Value"
  )

classification_plot <- ggplot(
  classification_plot_data,
  aes(
    x = factor(
      n_generated
    ),
    y = Value,
    fill = Model
  )
) +
  geom_boxplot(
    alpha = 0.82,
    outlier.alpha = 0.12
  ) +
  facet_wrap(
    ~ Metric,
    ncol = 2
  ) +
  scale_fill_manual(
    values = model_colors
  ) +
  labs(
    title =
      "Sensitivity and Specificity by Sample Size",
    x =
      "Generated sample size",
    y =
      "Metric value",
    fill =
      "Model"
  ) +
  theme_b_risa()

save_plot_pdf_png(
  classification_plot,
  file.path(
    dirs$performance_figs,
    "sensitivity_specificity_by_sample_size"
  ),
  width = 12,
  height = 7
)

print(classification_plot)

overall_summary_plot_data <-
  all_metrics |>
  dplyr::group_by(
    Model
  ) |>
  dplyr::summarise(
    AUC = mean(
      auc,
      na.rm = TRUE
    ),
    AUPRC = mean(
      auprc,
      na.rm = TRUE
    ),
    `Brier score` = mean(
      brier,
      na.rm = TRUE
    ),
    `Log loss` = mean(
      logloss,
      na.rm = TRUE
    ),
    `Mean uncertainty` =
      mean(
        mean_uncertainty,
        na.rm = TRUE
      ),
    `Selective risk` =
      mean(
        selective_risk,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  )

overall_summary_plot <- ggplot(
  overall_summary_plot_data,
  aes(
    x = Model,
    y = Value,
    fill = Model
  )
) +
  geom_col(
    width = 0.66,
    alpha = 0.9
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = model_colors
  ) +
  labs(
    title =
      "B-RISA Clinical-Calibrated Simulation: Model Summary",
    subtitle =
      "Mean results across all scenarios and replications",
    x = "Model",
    y = "Metric value",
    fill = "Model"
  ) +
  theme_b_risa()

save_plot_pdf_png(
  overall_summary_plot,
  file.path(
    dirs$performance_figs,
    "overall_model_summary"
  ),
  width = 12,
  height = 8
)

print(overall_summary_plot)

auc_heatmap_data <- all_metrics |>
  dplyr::filter(
    is.finite(auc)
  ) |>
  dplyr::mutate(
    Model =
      factor(
        Model,
        levels = MODEL_SET
      ),
    n_generated =
      factor(
        n_generated,
        levels = N_FILES
      )
  ) |>
  dplyr::group_by(
    n_generated,
    p,
    prevalence_target,
    Model
  ) |>
  dplyr::summarise(
    Mean_AUC =
      mean(
        auc,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Text_Color =
      dplyr::if_else(
        Mean_AUC >= 0.72,
        "black",
        "white"
      )
  )

auc_heatmap <- ggplot(
  auc_heatmap_data,
  aes(
    x = factor(p),
    y =
      factor(
        prevalence_target
      ),
    fill =
      Mean_AUC
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label =
        sprintf(
          "%.2f",
          Mean_AUC
        ),
      color =
        Text_Color
    ),
    fontface = "bold",
    size = 3
  ) +
  scale_color_identity() +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, 1),
    name = "Mean AUC"
  ) +
  facet_grid(
    Model ~
      n_generated,
    drop = FALSE
  ) +
  labs(
    title =
      "AUC Robustness Across Dimension and Prevalence",
    subtitle =
      paste0(
        "Means across boundary, covariance, label-noise, ",
        "missingness, and replication settings"
      ),
    x =
      "Number of predictors",
    y =
      "Target prevalence"
  ) +
  theme_b_risa()

save_plot_pdf_png(
  auc_heatmap,
  file.path(
    dirs$robustness_figs,
    "auc_heatmap_dimension_prevalence"
  ),
  width = 13,
  height = 8
)

print(auc_heatmap)

auc_heatmap_boundary_data <- all_metrics |>
  dplyr::filter(
    is.finite(auc)
  ) |>
  dplyr::mutate(
    Model =
      factor(
        Model,
        levels = MODEL_SET
      ),
    boundary =
      factor(
        boundary,
        levels = c(
          "linear",
          "nonlinear"
        )
      ),
    n_generated =
      factor(
        n_generated,
        levels = N_FILES
      )
  ) |>
  dplyr::group_by(
    n_generated,
    boundary,
    p,
    prevalence_target,
    Model
  ) |>
  dplyr::summarise(
    Mean_AUC =
      mean(
        auc,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Text_Color =
      dplyr::if_else(
        Mean_AUC >= 0.72,
        "black",
        "white"
      )
  )

auc_heatmap_boundary <- ggplot(
  auc_heatmap_boundary_data,
  aes(
    x =
      factor(p),
    y =
      factor(
        prevalence_target
      ),
    fill =
      Mean_AUC
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label =
        sprintf(
          "%.2f",
          Mean_AUC
        ),
      color =
        Text_Color
    ),
    fontface = "bold",
    size = 2.7
  ) +
  scale_color_identity() +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, 1),
    name = "Mean AUC"
  ) +
  facet_grid(
    Model ~
      n_generated +
      boundary,
    drop = FALSE
  ) +
  labs(
    title =
      paste0(
        "AUC Robustness Across Dimension, Prevalence, ",
        "and Boundary Type"
      ),
    subtitle =
      paste0(
        "Linear and nonlinear scenarios are displayed separately; ",
        "values are means across covariance, label-noise, ",
        "missingness, and replication settings"
      ),
    x =
      "Number of predictors",
    y =
      "Target prevalence"
  ) +
  theme_b_risa(
    base_size = 11
  ) +
  theme(
    strip.text =
      element_text(
        face = "bold",
        size = 8
      ),
    axis.text.x =
      element_text(
        size = 8
      )
  )

save_plot_pdf_png(
  auc_heatmap_boundary,
  file.path(
    dirs$robustness_figs,
    "auc_heatmap_dimension_prevalence_boundary"
  ),
  width = 20,
  height = 8
)

print(auc_heatmap_boundary)

boundary_covariance_plot_data <-
  all_metrics |>
  dplyr::group_by(
    n_generated,
    Model,
    boundary,
    covariance
  ) |>
  dplyr::summarise(
    Mean_AUC =
      mean(
        auc,
        na.rm = TRUE
      ),
    SD_AUC =
      stats::sd(
        auc,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

boundary_covariance_plot <- ggplot(
  boundary_covariance_plot_data,
  aes(
    x = boundary,
    y = Mean_AUC,
    fill = covariance
  )
) +
  geom_col(
    position =
      position_dodge(
        width = 0.8
      ),
    alpha = 0.9
  ) +
  geom_errorbar(
    aes(
      ymin =
        Mean_AUC -
        SD_AUC,
      ymax =
        Mean_AUC +
        SD_AUC
    ),
    position =
      position_dodge(
        width = 0.8
      ),
    width = 0.18
  ) +
  facet_grid(
    Model ~
      n_generated
  ) +
  labs(
    title =
      "AUC by Boundary and Covariance Structure",
    x =
      "Boundary type",
    y =
      "Mean AUC",
    fill =
      "Covariance"
  ) +
  theme_b_risa()

save_plot_pdf_png(
  boundary_covariance_plot,
  file.path(
    dirs$robustness_figs,
    "auc_by_boundary_covariance"
  ),
  width = 13,
  height = 8
)

print(boundary_covariance_plot)

noise_missingness_data <- all_metrics |>
  dplyr::group_by(
    n_generated,
    Model,
    label_noise,
    missing_rate_target
  ) |>
  dplyr::summarise(
    Mean_AUC =
      mean(
        auc,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

noise_missingness_plot <- ggplot(
  noise_missingness_data,
  aes(
    x =
      factor(
        missing_rate_target
      ),
    y =
      Mean_AUC,
    color =
      factor(
        label_noise
      ),
    group =
      factor(
        label_noise
      )
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_grid(
    Model ~
      n_generated
  ) +
  scale_color_viridis_d(
    option = "magma",
    name =
      "Label noise"
  ) +
  labs(
    title =
      "AUC Robustness to Missingness and Label Noise",
    x =
      "Target missingness rate",
    y =
      "Mean AUC"
  ) +
  theme_b_risa()

save_plot_pdf_png(
  noise_missingness_plot,
  file.path(
    dirs$robustness_figs,
    "auc_by_missingness_label_noise"
  ),
  width = 13,
  height = 8
)

print(noise_missingness_plot)

if (
  nrow(calibration_bins) >
    0L
) {
  calibration_plot <- ggplot(
    calibration_bins,
    aes(
      x =
        Mean_Predicted,
      y =
        Observed_Rate,
      color =
        Model
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray45"
    ) +
    geom_line(
      linewidth = 1
    ) +
    geom_point(
      aes(size = n),
      alpha = 0.85
    ) +
    facet_wrap(
      ~ n_generated,
      nrow = 1
    ) +
    scale_color_manual(
      values = model_colors
    ) +
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    labs(
      title =
        "Calibration Curves from Saved Prediction Samples",
      x =
        "Mean predicted probability",
      y =
        "Observed true event rate",
      color =
        "Model",
      size =
        "Bin n"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    calibration_plot,
    file.path(
      dirs$calibration_figs,
      "calibration_curves_by_sample_size"
    ),
    width = 12,
    height = 6
  )

  print(calibration_plot)

  density_plot_data <-
    all_prediction_samples |>
    dplyr::group_by(
      n_generated,
      Model
    ) |>
    dplyr::group_modify(
      ~ dplyr::slice_sample(
        .x,
        n = min(
          25000L,
          nrow(.x)
        )
      )
    ) |>
    dplyr::ungroup()

  density_plot <- ggplot(
    density_plot_data,
    aes(
      x = p_mean,
      fill = Outcome_true,
      color = Outcome_true
    )
  ) +
    geom_density(
      alpha = 0.25,
      linewidth = 0.9
    ) +
    facet_grid(
      Model ~
        n_generated
    ) +
    scale_fill_manual(
      values = outcome_colors
    ) +
    scale_color_manual(
      values = outcome_colors
    ) +
    coord_cartesian(
      xlim = c(0, 1)
    ) +
    labs(
      title =
        "Predicted Probability Distributions",
      x =
        "Predicted probability of true appendicitis",
      y =
        "Density",
      fill =
        "True outcome",
      color =
        "True outcome"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    density_plot,
    file.path(
      dirs$calibration_figs,
      "predicted_probability_density_by_sample_size"
    ),
    width = 13,
    height = 9
  )

  print(density_plot)
}

if (
  nrow(
    risk_summary_by_n
  ) > 0L
) {
  risk_plot <- ggplot(
    risk_summary_by_n,
    aes(
      x =
        Target_Coverage,
      y =
        Mean_Selective_Risk,
      color =
        Model
    )
  ) +
    geom_line(
      linewidth = 1.15
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~ n_generated,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = model_colors
    ) +
    scale_x_continuous(
      labels =
        scales::percent_format(
          accuracy = 1
        )
    ) +
    scale_y_continuous(
      labels =
        scales::percent_format(
          accuracy = 1
        )
    ) +
    labs(
      title =
        "Selective Prediction Risk-Coverage Curves",
      subtitle =
        "Mean risk across scenario-replication datasets; lower is better",
      x =
        "Coverage retained",
      y =
        "Mean selective risk",
      color =
        "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    risk_plot,
    file.path(
      dirs$selective_figs,
      "selective_risk_coverage_by_sample_size"
    ),
    width = 12,
    height = 6
  )

  print(risk_plot)
}

if (
  nrow(
    dca_summary_by_n
  ) > 0L
) {
  dca_plot <- ggplot(
    dca_summary_by_n,
    aes(
      x =
        Threshold,
      y =
        Mean_Net_Benefit,
      color =
        Strategy,
      linetype =
        Strategy
    )
  ) +
    geom_line(
      linewidth = 1
    ) +
    facet_wrap(
      ~ n_generated,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = model_colors
    ) +
    labs(
      title =
        "Simulation Decision-Curve Analysis",
      subtitle =
        "Mean net benefit across scenario-replication datasets",
      x =
        "Threshold probability",
      y =
        "Mean net benefit",
      color =
        "Strategy",
      linetype =
        "Strategy"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    dca_plot,
    file.path(
      dirs$dca_figs,
      "decision_curve_by_sample_size"
    ),
    width = 13,
    height = 6
  )

  print(dca_plot)

  dca_zoom_data <- dca_summary_by_n |>
    dplyr::filter(
      Threshold >= 0.05,
      Threshold <= 0.80
    )

  dca_zoom_limits <- compute_dca_zoom_limits(
    dca_data =
      dca_summary_by_n,
    threshold_min =
      0.05,
    threshold_max =
      0.80
  )

  dca_zoom_plot <- ggplot(
    dca_zoom_data,
    aes(
      x =
        Threshold,
      y =
        Mean_Net_Benefit,
      color =
        Strategy,
      linetype =
        Strategy
    )
  ) +
    geom_line(
      linewidth = 1
    ) +
    facet_wrap(
      ~ n_generated,
      nrow = 1
    ) +
    scale_color_manual(
      values = model_colors
    ) +
    coord_cartesian(
      xlim = c(
        0.05,
        0.80
      ),
      ylim =
        dca_zoom_limits,
      expand = FALSE
    ) +
    labs(
      title =
        "Simulation Decision Curves: Clinical Threshold Range",
      subtitle =
        paste0(
          "Thresholds 0.05 to 0.80; the y-axis is determined from ",
          "the model curves and Treat-none, while Treat-all is ",
          "retained and clipped outside the displayed range"
        ),
      x =
        "Threshold probability",
      y =
        "Mean net benefit",
      color =
        "Strategy",
      linetype =
        "Strategy"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    dca_zoom_plot,
    file.path(
      dirs$dca_figs,
      "decision_curve_by_sample_size_zoom"
    ),
    width = 13,
    height = 6
  )

  print(dca_zoom_plot)
}

# ------------------------------------------------
# 20.1 Separate figure suites for n = 100, 500, and 1000
# ------------------------------------------------

save_sample_plot <- function(
  plot_object,
  sample_dirs,
  category,
  file_name,
  width = 10,
  height = 7
) {
  output_directory <- sample_dirs[[category]]

  save_plot_pdf_png(
    plot_object,
    file.path(
      output_directory,
      file_name
    ),
    width = width,
    height = height
  )

  print(plot_object)
  invisible(plot_object)
}

generate_sample_size_figures <- function(
  n_value,
  metrics_data,
  prediction_data,
  calibration_data,
  risk_data,
  dca_data,
  sample_dirs
) {
  if (nrow(metrics_data) == 0L) {
    warning(
      "No metric rows are available for n = ",
      n_value,
      ". No separate figures were generated."
    )
    return(invisible(NULL))
  }

  plot_metric_for_sample <- function(
    metric_name,
    y_label,
    title_text,
    category,
    file_name
  ) {
    plot_object <- ggplot(
      metrics_data,
      aes(
        x = Model,
        y = .data[[metric_name]],
        fill = Model
      )
    ) +
      geom_boxplot(
        alpha = 0.82,
        outlier.alpha = 0.12,
        width = 0.68
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 23,
        size = 3,
        color = "black"
      ) +
      scale_fill_manual(
        values = model_colors[MODEL_SET]
      ) +
      labs(
        title = paste0(title_text, ": n = ", n_value),
        subtitle = "Independent scenario-replication test sets",
        x = "Model",
        y = y_label,
        fill = "Model"
      ) +
      theme_b_risa()

    save_sample_plot(
      plot_object,
      sample_dirs,
      category,
      paste0(file_name, "_n_", n_value),
      width = 9,
      height = 7
    )
  }

  plot_metric_for_sample(
    "auc",
    "AUC",
    "B-RISA Simulation Discrimination",
    "performance",
    "auc"
  )

  plot_metric_for_sample(
    "auprc",
    "AUPRC",
    "B-RISA Simulation Precision-Recall Performance",
    "performance",
    "auprc"
  )

  plot_metric_for_sample(
    "brier",
    "Brier score",
    "B-RISA Simulation Probability Error",
    "performance",
    "brier"
  )

  plot_metric_for_sample(
    "logloss",
    "Log loss",
    "B-RISA Simulation Logarithmic Loss",
    "performance",
    "logloss"
  )

  plot_metric_for_sample(
    "selective_risk",
    "Selective risk",
    "B-RISA Simulation Selective Prediction",
    "selective_prediction",
    "selective_risk"
  )

  plot_metric_for_sample(
    "mean_uncertainty",
    "Mean posterior predictive SD",
    "Posterior Predictive Uncertainty",
    "uncertainty",
    "uncertainty"
  )

  plot_metric_for_sample(
    "mean_interval_width",
    "Mean 95% interval width",
    "Posterior Predictive Interval Width",
    "uncertainty",
    "interval_width"
  )

  classification_data <- metrics_data |>
    dplyr::select(
      Model,
      sensitivity,
      specificity
    ) |>
    tidyr::pivot_longer(
      cols = c(sensitivity, specificity),
      names_to = "Metric",
      values_to = "Value"
    )

  classification_plot <- ggplot(
    classification_data,
    aes(
      x = Model,
      y = Value,
      fill = Model
    )
  ) +
    geom_boxplot(
      alpha = 0.82,
      outlier.alpha = 0.12,
      width = 0.68
    ) +
    facet_wrap(~ Metric, ncol = 2) +
    scale_fill_manual(values = model_colors[MODEL_SET]) +
    labs(
      title = paste0(
        "Sensitivity and Specificity: n = ",
        n_value
      ),
      x = "Model",
      y = "Metric value",
      fill = "Model"
    ) +
    theme_b_risa()

  save_sample_plot(
    classification_plot,
    sample_dirs,
    "performance",
    paste0(
      "sensitivity_specificity_n_",
      n_value
    ),
    width = 11,
    height = 7
  )

  overall_plot_data <- metrics_data |>
    dplyr::group_by(Model) |>
    dplyr::summarise(
      AUC = mean(auc, na.rm = TRUE),
      AUPRC = mean(auprc, na.rm = TRUE),
      `Brier score` = mean(brier, na.rm = TRUE),
      `Log loss` = mean(logloss, na.rm = TRUE),
      `Mean uncertainty` = mean(mean_uncertainty, na.rm = TRUE),
      `Selective risk` = mean(selective_risk, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      cols = -Model,
      names_to = "Metric",
      values_to = "Value"
    )

  overall_plot <- ggplot(
    overall_plot_data,
    aes(
      x = Model,
      y = Value,
      fill = Model
    )
  ) +
    geom_col(width = 0.66, alpha = 0.9) +
    facet_wrap(
      ~ Metric,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(values = model_colors[MODEL_SET]) +
    labs(
      title = paste0(
        "B-RISA Simulation Model Summary: n = ",
        n_value
      ),
      subtitle = "Mean results across scenarios and replications",
      x = "Model",
      y = "Metric value",
      fill = "Model"
    ) +
    theme_b_risa()

  save_sample_plot(
    overall_plot,
    sample_dirs,
    "performance",
    paste0(
      "overall_model_summary_n_",
      n_value
    ),
    width = 12,
    height = 8
  )

  heatmap_data <- metrics_data |>
    dplyr::filter(
      is.finite(auc)
    ) |>
    dplyr::mutate(
      Model =
        factor(
          Model,
          levels = MODEL_SET
        )
    ) |>
    dplyr::group_by(
      p,
      prevalence_target,
      Model
    ) |>
    dplyr::summarise(
      Mean_AUC =
        mean(
          auc,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Text_Color =
        dplyr::if_else(
          Mean_AUC >= 0.72,
          "black",
          "white"
        )
    )

  heatmap_plot <- ggplot(
    heatmap_data,
    aes(
      x =
        factor(p),
      y =
        factor(
          prevalence_target
        ),
      fill =
        Mean_AUC
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        label =
          sprintf(
            "%.2f",
            Mean_AUC
          ),
        color =
          Text_Color
      ),
      fontface = "bold",
      size = 3.2
    ) +
    scale_color_identity() +
    scale_fill_viridis_c(
      option = "plasma",
      limits = c(0, 1),
      name = "Mean AUC"
    ) +
    facet_wrap(
      ~ Model,
      ncol = 3,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "AUC Robustness Across Dimension and Prevalence: n = ",
        n_value
      ),
      subtitle =
        paste0(
          "Means across boundary, covariance, label-noise, ",
          "missingness, and replication settings"
        ),
      x =
        "Number of predictors",
      y =
        "Target prevalence"
    ) +
    theme_b_risa()

  save_sample_plot(
    heatmap_plot,
    sample_dirs,
    "robustness",
    paste0(
      "auc_heatmap_dimension_prevalence_n_",
      n_value
    ),
    width = 12,
    height = 6
  )

  heatmap_boundary_data <- metrics_data |>
    dplyr::filter(
      is.finite(auc)
    ) |>
    dplyr::mutate(
      Model =
        factor(
          Model,
          levels = MODEL_SET
        ),
      boundary =
        factor(
          boundary,
          levels = c(
            "linear",
            "nonlinear"
          )
        )
    ) |>
    dplyr::group_by(
      boundary,
      p,
      prevalence_target,
      Model
    ) |>
    dplyr::summarise(
      Mean_AUC =
        mean(
          auc,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Text_Color =
        dplyr::if_else(
          Mean_AUC >= 0.72,
          "black",
          "white"
        )
    )

  heatmap_boundary_plot <- ggplot(
    heatmap_boundary_data,
    aes(
      x =
        factor(p),
      y =
        factor(
          prevalence_target
        ),
      fill =
        Mean_AUC
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.5
    ) +
    geom_text(
      aes(
        label =
          sprintf(
            "%.2f",
            Mean_AUC
          ),
        color =
          Text_Color
      ),
      fontface = "bold",
      size = 3
    ) +
    scale_color_identity() +
    scale_fill_viridis_c(
      option = "plasma",
      limits = c(0, 1),
      name = "Mean AUC"
    ) +
    facet_grid(
      boundary ~
        Model,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "AUC Robustness by Boundary Type: n = ",
        n_value
      ),
      subtitle =
        paste0(
          "Linear and nonlinear scenarios are displayed separately; ",
          "values are means across covariance, label-noise, ",
          "missingness, and replication settings"
        ),
      x =
        "Number of predictors",
      y =
        "Target prevalence"
    ) +
    theme_b_risa()

  save_sample_plot(
    heatmap_boundary_plot,
    sample_dirs,
    "robustness",
    paste0(
      "auc_heatmap_dimension_prevalence_boundary_n_",
      n_value
    ),
    width = 12,
    height = 8
  )

  boundary_data <- metrics_data |>
    dplyr::group_by(
      Model,
      boundary,
      covariance
    ) |>
    dplyr::summarise(
      Mean_AUC = mean(auc, na.rm = TRUE),
      SD_AUC = stats::sd(auc, na.rm = TRUE),
      .groups = "drop"
    )

  boundary_plot <- ggplot(
    boundary_data,
    aes(
      x = boundary,
      y = Mean_AUC,
      fill = covariance
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      alpha = 0.9
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(Mean_AUC - SD_AUC, 0),
        ymax = pmin(Mean_AUC + SD_AUC, 1)
      ),
      position = position_dodge(width = 0.8),
      width = 0.18
    ) +
    facet_wrap(~ Model, ncol = 3) +
    labs(
      title = paste0(
        "AUC by Boundary and Covariance Structure: n = ",
        n_value
      ),
      x = "Boundary type",
      y = "Mean AUC",
      fill = "Covariance"
    ) +
    theme_b_risa()

  save_sample_plot(
    boundary_plot,
    sample_dirs,
    "robustness",
    paste0(
      "auc_by_boundary_covariance_n_",
      n_value
    ),
    width = 12,
    height = 7
  )

  noise_data <- metrics_data |>
    dplyr::group_by(
      Model,
      label_noise,
      missing_rate_target
    ) |>
    dplyr::summarise(
      Mean_AUC = mean(auc, na.rm = TRUE),
      .groups = "drop"
    )

  noise_plot <- ggplot(
    noise_data,
    aes(
      x = factor(missing_rate_target),
      y = Mean_AUC,
      color = factor(label_noise),
      group = factor(label_noise)
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.2) +
    facet_wrap(~ Model, ncol = 3) +
    scale_color_viridis_d(
      option = "magma",
      name = "Label noise"
    ) +
    labs(
      title = paste0(
        "AUC Robustness to Missingness and Label Noise: n = ",
        n_value
      ),
      x = "Target missingness rate",
      y = "Mean AUC"
    ) +
    theme_b_risa()

  save_sample_plot(
    noise_plot,
    sample_dirs,
    "robustness",
    paste0(
      "auc_by_missingness_label_noise_n_",
      n_value
    ),
    width = 12,
    height = 7
  )

  if (nrow(calibration_data) > 0L) {
    calibration_plot <- ggplot(
      calibration_data,
      aes(
        x = Mean_Predicted,
        y = Observed_Rate,
        color = Model,
        group = Model
      )
    ) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "gray45"
      ) +
      geom_line(linewidth = 1) +
      geom_point(
        aes(size = n),
        alpha = 0.85
      ) +
      scale_color_manual(values = model_colors[MODEL_SET]) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(
        title = paste0(
          "Calibration Curves: n = ",
          n_value
        ),
        x = "Mean predicted probability",
        y = "Observed true event rate",
        color = "Model",
        size = "Bin n"
      ) +
      theme_b_risa()

    save_sample_plot(
      calibration_plot,
      sample_dirs,
      "calibration",
      paste0(
        "calibration_curves_n_",
        n_value
      ),
      width = 8,
      height = 7
    )
  }

  if (nrow(prediction_data) > 0L) {
    set.seed(BASE_SEED + n_value)

    density_data <- prediction_data |>
      dplyr::group_by(Model) |>
      dplyr::group_modify(
        ~ dplyr::slice_sample(
          .x,
          n = min(25000L, nrow(.x))
        )
      ) |>
      dplyr::ungroup()

    density_plot <- ggplot(
      density_data,
      aes(
        x = p_mean,
        fill = Outcome_true,
        color = Outcome_true
      )
    ) +
      geom_density(alpha = 0.25, linewidth = 0.9) +
      facet_wrap(~ Model, ncol = 1) +
      scale_fill_manual(values = outcome_colors) +
      scale_color_manual(values = outcome_colors) +
      coord_cartesian(xlim = c(0, 1)) +
      labs(
        title = paste0(
          "Predicted Probability Distributions: n = ",
          n_value
        ),
        x = "Predicted probability of true appendicitis",
        y = "Density",
        fill = "True outcome",
        color = "True outcome"
      ) +
      theme_b_risa()

    save_sample_plot(
      density_plot,
      sample_dirs,
      "calibration",
      paste0(
        "predicted_probability_density_n_",
        n_value
      ),
      width = 8,
      height = 10
    )
  }

  if (nrow(risk_data) > 0L) {
    risk_plot <- ggplot(
      risk_data,
      aes(
        x = Target_Coverage,
        y = Mean_Selective_Risk,
        color = Model
      )
    ) +
      geom_line(linewidth = 1.15) +
      geom_point(size = 2) +
      scale_color_manual(values = model_colors[MODEL_SET]) +
      scale_x_continuous(
        labels = scales::percent_format(accuracy = 1)
      ) +
      scale_y_continuous(
        labels = scales::percent_format(accuracy = 1)
      ) +
      labs(
        title = paste0(
          "Selective Prediction Risk-Coverage: n = ",
          n_value
        ),
        subtitle = "Lower risk is better",
        x = "Coverage retained",
        y = "Mean selective risk",
        color = "Model"
      ) +
      theme_b_risa()

    save_sample_plot(
      risk_plot,
      sample_dirs,
      "selective_prediction",
      paste0(
        "selective_risk_coverage_n_",
        n_value
      ),
      width = 8,
      height = 6
    )
  }

  if (nrow(dca_data) > 0L) {
    dca_plot <- ggplot(
      dca_data,
      aes(
        x = Threshold,
        y = Mean_Net_Benefit,
        color = Strategy,
        linetype = Strategy
      )
    ) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = model_colors) +
      labs(
        title = paste0(
          "Simulation Decision-Curve Analysis: n = ",
          n_value
        ),
        x = "Threshold probability",
        y = "Mean net benefit",
        color = "Strategy",
        linetype = "Strategy"
      ) +
      theme_b_risa()

    save_sample_plot(
      dca_plot,
      sample_dirs,
      "decision_curve",
      paste0(
        "decision_curve_n_",
        n_value
      ),
      width = 9,
      height = 6
    )

    dca_zoom_data <- dca_data |>
      dplyr::filter(
        Threshold >= 0.05,
        Threshold <= 0.80
      )

    dca_zoom_limits <- compute_dca_zoom_limits(
      dca_data =
        dca_data,
      threshold_min =
        0.05,
      threshold_max =
        0.80
    )

    dca_zoom_plot <- ggplot(
      dca_zoom_data,
      aes(
        x =
          Threshold,
        y =
          Mean_Net_Benefit,
        color =
          Strategy,
        linetype =
          Strategy
      )
    ) +
      geom_line(
        linewidth = 1
      ) +
      scale_color_manual(
        values = model_colors
      ) +
      coord_cartesian(
        xlim = c(
          0.05,
          0.80
        ),
        ylim =
          dca_zoom_limits,
        expand = FALSE
      ) +
      labs(
        title = paste0(
          "Decision Curves, Clinical Threshold Range: n = ",
          n_value
        ),
        subtitle =
          paste0(
            "Thresholds 0.05 to 0.80; the y-axis is determined from ",
            "the model curves and Treat-none, while Treat-all is ",
            "retained and clipped outside the displayed range"
          ),
        x =
          "Threshold probability",
        y =
          "Mean net benefit",
        color =
          "Strategy",
        linetype =
          "Strategy"
      ) +
      theme_b_risa()

    save_sample_plot(
      dca_zoom_plot,
      sample_dirs,
      "decision_curve",
      paste0(
        "decision_curve_zoom_n_",
        n_value
      ),
      width = 9,
      height = 6
    )
  }

  invisible(TRUE)
}

for (n_value in N_FILES) {
  n_key <- as.character(n_value)

  generate_sample_size_figures(
    n_value = n_value,
    metrics_data = all_metrics |>
      dplyr::filter(n_generated == n_value),
    prediction_data = all_prediction_samples |>
      dplyr::filter(n_generated == n_value),
    calibration_data = calibration_bins |>
      dplyr::filter(n_generated == n_value),
    risk_data = risk_summary_by_n |>
      dplyr::filter(n_generated == n_value),
    dca_data = dca_summary_by_n |>
      dplyr::filter(n_generated == n_value),
    sample_dirs = sample_figure_dirs[[n_key]]
  )
}


# ------------------------------------------------
# 21. Methodology manifest, three CSV exports, and reproducibility
# ------------------------------------------------

method_manifest <- tibble::tibble(
  Item = c(
    "Analysis version",
    "Input directory",
    "Output root",
    "Generated sample sizes",
    "Validation design",
    "Training label",
    "Evaluation label",
    "Stratification label",
    "Training proportion",
    "Models",
    "BLR method",
    "BLR intercept prior SD",
    "BLR coefficient prior SD",
    "BLR optimizer relative tolerance",
    "BLR posterior draws",
    "BLDA method",
    "BLDA mean prior kappa",
    "BLDA inverse-Wishart df rule",
    "BLDA inverse-Wishart scale",
    "BLDA class-prior Beta alpha",
    "BLDA class-prior Beta beta",
    "BLDA posterior draws",
    "BART method",
    "BART probability transformation",
    "BART trees",
    "BART burn-in",
    "BART posterior draws",
    "Classification threshold",
    "Calibration bins",
    "Selective coverage",
    "Decision thresholds",
    "Probability clipping",
    "Prediction-draw validation",
    "Degenerate prediction tolerance",
    "Failure policy",
    "Random seed",
    "Seed construction",
    "Checkpoint format",
    "User-facing CSV count"
  ),
  Value = c(
    ANALYSIS_VERSION,
    INPUT_DIR,
    OUTPUT_ROOT,
    paste(
      N_FILES,
      collapse = ", "
    ),
    paste0(
      "One true-label-stratified ",
      round(
        TRAIN_PROP * 100
      ),
      "/",
      round(
        (1 - TRAIN_PROP) * 100
      ),
      " train/test split per independent scenario-replication"
    ),
    "Y_observed",
    "Y_true",
    "Y_true",
    as.character(
      TRAIN_PROP
    ),
    paste(
      MODEL_SET,
      collapse = ", "
    ),
    paste0(
      "Regularized Bayesian logistic regression; MAP plus Laplace posterior; ",
      "training-only median imputation and standardization"
    ),
    as.character(
      BLR_INTERCEPT_PRIOR_SD
    ),
    as.character(
      BLR_COEFFICIENT_PRIOR_SD
    ),
    as.character(
      BLR_OPTIM_RELTOL
    ),
    as.character(
      N_POST_DRAWS_BLR
    ),
    paste0(
      "Genuine Bayesian LDA with class-specific means, one shared covariance, ",
      "Normal-Inverse-Wishart parameter posterior and Beta class-prevalence posterior"
    ),
    as.character(
      BLDA_MEAN_PRIOR_KAPPA
    ),
    paste0(
      "nu0 = p + ",
      BLDA_IW_DF_OFFSET
    ),
    paste0(
      "Lambda0 = ",
      BLDA_PRIOR_COVARIANCE_SCALE,
      " * I_p"
    ),
    as.character(
      BLDA_CLASS_PRIOR_ALPHA
    ),
    as.character(
      BLDA_CLASS_PRIOR_BETA
    ),
    as.character(
      N_POST_DRAWS_BLDA
    ),
    "Probit Bayesian additive regression trees",
    "pnorm applied to latent yhat.test draws",
    as.character(
      BART_NTREE
    ),
    as.character(
      N_BURN_BART
    ),
    as.character(
      N_POST_DRAWS_BART
    ),
    as.character(
      CLASSIFICATION_THRESHOLD
    ),
    as.character(
      CALIBRATION_BINS
    ),
    as.character(
      SELECTIVE_COVERAGE
    ),
    paste0(
      min(
        DCA_THRESHOLDS
      ),
      " to ",
      max(
        DCA_THRESHOLDS
      ),
      " by ",
      DCA_THRESHOLDS[2L] -
        DCA_THRESHOLDS[1L]
    ),
    paste0(
      "Dimension-preserving truncation to [",
      EPS,
      ", ",
      1 - EPS,
      "]"
    ),
    paste0(
      "Expected posterior draws x test patients; ",
      "all finite and within [0,1]"
    ),
    as.character(
      DEGENERATE_PREDICTION_TOL
    ),
    paste0(
      "STOP_ON_MODEL_FAILURE=",
      STOP_ON_MODEL_FAILURE,
      "; STOP_ON_DEGENERATE_PREDICTIONS=",
      STOP_ON_DEGENERATE_PREDICTIONS
    ),
    as.character(
      BASE_SEED
    ),
    paste0(
      "Deterministic formula evaluated in double precision and reduced ",
      "modulo 2147483646 only when necessary"
    ),
    paste0(
      "RDS chunks every ",
      CHECKPOINT_EVERY_UNITS,
      " newly processed units"
    ),
    "3"
  )
)

tag_table <- function(
  data,
  table_id
) {
  if (
    is.null(data) ||
      nrow(data) == 0L
  ) {
    return(
      tibble::tibble(
        Table_ID =
          character(0)
      )
    )
  }

  output <- as.data.frame(
    data
  )

  output[] <- lapply(
    output,
    function(x) {
      if (is.factor(x)) {
        as.character(x)
      } else {
        x
      }
    }
  )

  output |>
    dplyr::mutate(
      Table_ID =
        table_id,
      .before = 1
    )
}

# CSV 1: one row per complete scenario-replication/model result.
csv_1_unit_metrics <- all_metrics |>
  dplyr::arrange(
    n_generated,
    scenario_id,
    rep_id,
    factor(
      Model,
      levels = MODEL_SET
    )
  )

# CSV 2: all user-facing summary tables are stacked in one file and
# identified by Table_ID. Filtering Table_ID reconstructs each table.
csv_2_summary_tables <- dplyr::bind_rows(
  tag_table(
    overall_model_summary,
    "overall_model_summary"
  ),
  tag_table(
    summary_by_sample_size,
    "summary_by_sample_size"
  ),
  tag_table(
    summary_by_dimension,
    "summary_by_dimension"
  ),
  tag_table(
    summary_by_boundary,
    "summary_by_boundary"
  ),
  tag_table(
    summary_by_design,
    "summary_by_design_condition"
  ),
  tag_table(
    summary_by_boundary_covariance,
    "summary_by_boundary_covariance"
  ),
  tag_table(
    summary_by_noise_missingness,
    "summary_by_noise_missingness"
  ),
  tag_table(
    summary_by_calibration_source,
    "summary_by_calibration_source"
  ),
  tag_table(
    dca_summary_by_design,
    "decision_curve_summary_by_design"
  ),
  tag_table(
    dca_summary_by_n,
    "decision_curve_summary_by_sample_size"
  ),
  tag_table(
    risk_summary_by_design,
    "risk_coverage_summary_by_design"
  ),
  tag_table(
    risk_summary_by_n,
    "risk_coverage_summary_by_sample_size"
  ),
  tag_table(
    calibration_bins,
    "calibration_bins_from_prediction_sample"
  ),
  tag_table(
    unique_best_auc_counts,
    "unique_best_auc_counts"
  )
)

# The manifest uses dedicated column names so its character values
# cannot collide with numerical result columns when tables are stacked.
manifest_for_export <- method_manifest |>
  dplyr::transmute(
    Manifest_Item =
      Item,
    Manifest_Value =
      Value
  )

# CSV 3: inferential tests, status tables, completeness tables and
# reproducibility/audit tables in one file.
csv_3_inference_and_audits <- dplyr::bind_rows(
  tag_table(
    simulation_omnibus_comparisons,
    "friedman_omnibus_comparisons"
  ),
  tag_table(
    simulation_pairwise_comparisons,
    "paired_wilcoxon_comparisons"
  ),
  tag_table(
    model_status_audit,
    "model_status_audit"
  ),
  tag_table(
    metric_completeness,
    "metric_completeness_by_unit"
  ),
  tag_table(
    sample_size_run_audit,
    "sample_size_completion_audit"
  ),
  tag_table(
    all_input_unit_audits,
    "input_unit_size_audit"
  ),
  tag_table(
    all_status,
    "model_fit_status_and_diagnostics"
  ),
  tag_table(
    all_unit_status,
    "unit_status"
  ),
  tag_table(
    manifest_for_export,
    "methodology_manifest"
  )
)

CSV_UNIT_METRICS <- file.path(
  RESULT_ROOT,
  "01_simulation_unit_model_metrics.csv"
)

CSV_SUMMARY_TABLES <- file.path(
  RESULT_ROOT,
  "02_simulation_summary_tables.csv"
)

CSV_INFERENCE_AUDITS <- file.path(
  RESULT_ROOT,
  "03_simulation_inference_and_audits.csv"
)

# These are the ONLY three CSV writes in this script.
safe_fwrite(
  csv_1_unit_metrics,
  CSV_UNIT_METRICS
)

safe_fwrite(
  csv_2_summary_tables,
  CSV_SUMMARY_TABLES
)

safe_fwrite(
  csv_3_inference_and_audits,
  CSV_INFERENCE_AUDITS
)

# Prediction samples are retained for calibration/density
# reproducibility as RDS rather than creating a fourth CSV.
atomic_save_rds(
  list(
    analysis_version =
      ANALYSIS_VERSION,
    prediction_samples =
      all_prediction_samples
  ),
  file.path(
    RESULT_ROOT,
    "simulation_prediction_samples_for_figures.rds"
  )
)

# Compact final analysis object.
atomic_save_rds(
  list(
    analysis_version =
      ANALYSIS_VERSION,
    method_manifest =
      method_manifest,
    overall_model_summary =
      overall_model_summary,
    summary_by_sample_size =
      summary_by_sample_size,
    summary_by_dimension =
      summary_by_dimension,
    summary_by_boundary =
      summary_by_boundary,
    summary_by_design =
      summary_by_design,
    summary_by_boundary_covariance =
      summary_by_boundary_covariance,
    summary_by_noise_missingness =
      summary_by_noise_missingness,
    summary_by_calibration_source =
      summary_by_calibration_source,
    unique_best_auc_counts =
      unique_best_auc_counts,
    model_status_audit =
      model_status_audit,
    sample_size_run_audit =
      sample_size_run_audit,
    dca_summary_by_design =
      dca_summary_by_design,
    dca_summary_by_n =
      dca_summary_by_n,
    risk_summary_by_design =
      risk_summary_by_design,
    risk_summary_by_n =
      risk_summary_by_n,
    calibration_bins =
      calibration_bins,
    simulation_omnibus_comparisons =
      simulation_omnibus_comparisons,
    simulation_pairwise_comparisons =
      simulation_pairwise_comparisons
  ),
  file.path(
    RESULT_ROOT,
    "b_risa_simulation_bayesian_lda_analysis.rds"
  )
)

session_information <- capture.output(
  sessionInfo()
)

writeLines(
  session_information,
  file.path(
    dirs$manifests,
    "session_info.txt"
  )
)

writeLines(
  c(
    paste0(
      "Analysis version: ",
      ANALYSIS_VERSION
    ),
    paste0(
      "Input directory: ",
      INPUT_DIR
    ),
    paste0(
      "Output root: ",
      OUTPUT_ROOT
    ),
    "",
    "User-facing CSV files:",
    basename(
      CSV_UNIT_METRICS
    ),
    basename(
      CSV_SUMMARY_TABLES
    ),
    basename(
      CSV_INFERENCE_AUDITS
    ),
    "",
    paste0(
      "Completed three-model units: ",
      length(
        complete_unit_ids
      )
    ),
    paste0(
      "Expected full-study units: ",
      576L * 100L
    )
  ),
  file.path(
    dirs$manifests,
    "README_outputs.txt"
  )
)

# ------------------------------------------------
# 22. Final checks and completion message
# ------------------------------------------------

if (
  nrow(
    model_status_audit
  ) > 0L &&
    any(
      model_status_audit$
      Failed_fits >
        0L |
        model_status_audit$
        Degenerate_fits >
        0L
    )
) {
  warning(
    "At least one failed or degenerate model fit is present. ",
    "Inspect Table_ID='model_fit_status_and_diagnostics' and ",
    "Table_ID='model_status_audit' in ",
    basename(
      CSV_INFERENCE_AUDITS
    ),
    "."
  )
}

expected_complete_units <- 576L * 100L

if (
  is.infinite(
    MAX_UNITS_PER_N
  ) &&
    length(
      complete_unit_ids
    ) !=
      expected_complete_units
) {
  warning(
    "The full run contains ",
    length(
      complete_unit_ids
    ),
    " complete units; expected ",
    expected_complete_units,
    ". Review the completion and status audits."
  )
}

cat("\n================================================\n")
cat("B-RISA SIMULATION ANALYSIS COMPLETE\n")
cat("GENUINE BAYESIAN LDA FINAL ANALYSIS\n")
cat("================================================\n")

cat("\nInput directory:\n")
cat(INPUT_DIR, "\n")

cat("\nOutput root:\n")
cat(OUTPUT_ROOT, "\n")

cat("\nExactly three user-facing CSV files were written:\n")
cat("  ", CSV_UNIT_METRICS, "\n")
cat("  ", CSV_SUMMARY_TABLES, "\n")
cat("  ", CSV_INFERENCE_AUDITS, "\n")

cat("\nInternal RDS checkpoints:\n")
cat("  ", CHECKPOINT_ROOT, "\n")

cat("\nFigure folders:\n")
cat(
  "  ",
  file.path(
    FIGURE_ROOT,
    "n_100"
  ),
  "\n"
)
cat(
  "  ",
  file.path(
    FIGURE_ROOT,
    "n_500"
  ),
  "\n"
)
cat(
  "  ",
  file.path(
    FIGURE_ROOT,
    "n_1000"
  ),
  "\n"
)
cat(
  "  ",
  file.path(
    FIGURE_ROOT,
    "combined"
  ),
  "\n"
)

cat("\nBLDA implementation:\n")
cat(
  "  Normal-Inverse-Wishart shared-covariance Bayesian LDA\n"
)
cat(
  "  + Beta posterior for class prevalence\n"
)
cat(
  "  + ",
  N_POST_DRAWS_BLDA,
  " posterior probability draws per fit\n",
  sep = ""
)

cat("\nThe generated simulation data were not regenerated.\n")
cat("================================================\n")

