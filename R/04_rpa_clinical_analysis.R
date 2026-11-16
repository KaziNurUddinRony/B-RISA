# ================================================================
# B-RISA RPA Clinical Analysis
#
# Reproduces the manuscript's repeated stratified 10-fold cross-validation
# analysis for standardized RPA data with body temperature, WBC count, and neutrophil percentage. BLR, BLDA, and BART use the final model
# definitions reported in the manuscript. The statistical-analysis stage is
# included in this same script and reads the saved clinical outputs.
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
RPA_STD_CLEAN_FILE <- file.path(
  PROJECT_ROOT,
  "data", "derived",
  "rpa_standardized_selected_cleaned_data.csv"
)

OUTPUT_ROOT <- file.path(PROJECT_ROOT, "analysis_outputs", "rpa_standardized_clinical")
RESULT_ROOT <- file.path(OUTPUT_ROOT, "results")
FIGURE_ROOT <- file.path(OUTPUT_ROOT, "figures")

dir.create(
  OUTPUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

DATASET_ID <- "rpa_standardized"
DATASET_LABEL <- "RPA Standardized"

# Use only the intended standardized RPA predictors.
# The input file is already standardized. The script nevertheless
# re-centers and re-scales each predictor using the training fold only.
# Algebraically, this removes dependence on the upstream global scaling
# constants and preserves fold-specific preprocessing.
RPA_STD_FEATURES <- c(
  "Body_Temperature",
  "WBC_Count",
  "Neutrophil_Percentage"
)

# ------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------

CLEAN_OLD_OUTPUT <- FALSE
INSTALL_MISSING_PACKAGES <- FALSE

SEED <- 20260618L
N_FOLDS <- 10L
N_REPEATS <- 50L

MODEL_SET <- c("BLR", "BLDA", "BART")

N_POST_DRAWS_BLR <- 300L
N_POST_DRAWS_BLDA <- 300L
N_POST_DRAWS_BART <- 300L
N_BURN_BART <- 150L
BART_NTREE <- 50L

# BLR prior and optimization settings.
BLR_INTERCEPT_PRIOR_SD <- 2.5
BLR_COEFFICIENT_PRIOR_SD <- 2.5
BLR_OPTIM_MAXIT <- 2000L
BLR_OPTIM_RELTOL <- 1e-10
BLR_MAX_ABS_GRADIENT <- 1e-5
BLR_HESSIAN_EIGEN_FLOOR <- 1e-8

# Genuine Bayesian LDA priors. These match the updated simulation analysis.
# Training folds are standardized using training-only moments, so a zero-centered
# prior on class means and identity-centered prior on shared covariance are used.
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

# A valid fitted model should not assign exactly the same posterior
# mean probability to every test patient when the test fold contains
# more than one patient.
DEGENERATE_PREDICTION_TOL <- 1e-12
STOP_ON_MODEL_FAILURE <- TRUE
STOP_ON_DEGENERATE_PREDICTIONS <- TRUE

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
  "MASS"
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

  install.packages(missing_packages, dependencies = TRUE)
}

if ("BART" %in% MODEL_SET &&
    !requireNamespace("BART", quietly = TRUE)) {
  if (INSTALL_MISSING_PACKAGES) {
    install.packages("BART", dependencies = TRUE)
  }
}

if ("BART" %in% MODEL_SET &&
    !requireNamespace("BART", quietly = TRUE)) {
  stop(
    "The BART package is required because BART is in MODEL_SET."
  )
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
})

# ------------------------------------------------
# 3. Output folders
# ------------------------------------------------

RESULT_DIR <- file.path(RESULT_ROOT, DATASET_ID)
FIGURE_DIR <- file.path(FIGURE_ROOT, DATASET_ID)

dirs <- list(
  result_base = RESULT_DIR,
  figure_base = FIGURE_DIR,
  metrics = file.path(RESULT_DIR, "metrics"),
  predictions = file.path(RESULT_DIR, "predictions"),
  summaries = file.path(RESULT_DIR, "summaries"),
  status = file.path(RESULT_DIR, "status"),
  dca = file.path(RESULT_DIR, "decision_curve"),
  performance_figs = file.path(FIGURE_DIR, "performance"),
  calibration_figs = file.path(FIGURE_DIR, "calibration"),
  uncertainty_figs = file.path(FIGURE_DIR, "uncertainty"),
  selective_figs = file.path(FIGURE_DIR, "selective_prediction"),
  dca_figs = file.path(FIGURE_DIR, "decision_curve")
)

if (CLEAN_OLD_OUTPUT) {
  unlink(RESULT_DIR, recursive = TRUE, force = TRUE)
  unlink(FIGURE_DIR, recursive = TRUE, force = TRUE)
}

for (directory in dirs) {
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
# 6. Repeated stratified folds
# ------------------------------------------------

make_repeated_stratified_folds <- function(
  y,
  k,
  repeats,
  seed
) {
  y <- as.integer(y)

  if (min(table(y)) < k) {
    stop(
      "Each outcome class must contain at least k observations."
    )
  }

  set.seed(seed)

  index_0 <- which(y == 0L)
  index_1 <- which(y == 1L)
  all_folds <- vector(
    mode = "list",
    length = repeats * k
  )

  list_index <- 0L

  for (repeat_id in seq_len(repeats)) {
    shuffled_0 <- sample(index_0)
    shuffled_1 <- sample(index_1)

    fold_0 <- split(
      shuffled_0,
      rep(seq_len(k), length.out = length(shuffled_0))
    )

    fold_1 <- split(
      shuffled_1,
      rep(seq_len(k), length.out = length(shuffled_1))
    )

    for (fold_id in seq_len(k)) {
      list_index <- list_index + 1L

      test_index <- sample(
        c(
          fold_0[[fold_id]],
          fold_1[[fold_id]]
        )
      )

      train_index <- setdiff(
        seq_along(y),
        test_index
      )

      if (
        length(unique(y[test_index])) < 2L ||
        length(unique(y[train_index])) < 2L
      ) {
        stop(
          "A stratified fold does not contain both classes."
        )
      }

      all_folds[[list_index]] <- list(
        repeat_id = repeat_id,
        fold_id = fold_id,
        train_index = train_index,
        test_index = test_index
      )
    }
  }

  all_folds
}

make_fit_seed <- function(
  repeat_id,
  fold_id,
  model_name
) {
  model_code <- match(model_name, MODEL_SET)

  as.integer(
    SEED +
      repeat_id * 10000L +
      fold_id * 100L +
      model_code
  )
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

  if (
    nrow(probability_draws) == nrow(x_test) &&
    ncol(probability_draws) == n_draws
  ) {
    probability_draws <- t(
      probability_draws
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
# 10. Prediction aggregation
# ------------------------------------------------

aggregate_repeated_predictions <- function(
  all_predictions
) {
  all_predictions |>
    dplyr::group_by(
      Row_ID,
      Model,
      Y_true,
      Outcome_true
    ) |>
    dplyr::summarise(
      count = dplyr::n(),
      mean_probability = mean(p_mean),
      mean_second_moment = mean(
        p_sd^2 +
          p_mean^2
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      predictive_variance = pmax(
        mean_second_moment -
          mean_probability^2,
        0
      ),
      p_mean = mean_probability,
      p_sd = sqrt(predictive_variance),
      p_lwr = pmax(
        0,
        p_mean -
          1.96 * p_sd
      ),
      p_upr = pmin(
        1,
        p_mean +
          1.96 * p_sd
      )
    ) |>
    dplyr::select(
      Row_ID,
      Model,
      Y_true,
      Outcome_true,
      count,
      p_mean,
      p_sd,
      p_lwr,
      p_upr
    )
}

decision_curve_table <- function(
  prediction_data,
  thresholds
) {
  model_rows <- lapply(
    unique(prediction_data$Model),
    function(model_name) {
      model_data <- prediction_data[
        prediction_data$Model == model_name,
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
              predicted_positive &
                y == 1L
            )

            false_positive <- sum(
              predicted_positive &
                y == 0L
            )

            tibble::tibble(
              Model = model_name,
              Threshold = threshold,
              Net_Benefit = (
                true_positive / n
              ) - (
                false_positive / n
              ) * (
                threshold /
                  (1 - threshold)
              ),
              Prevalence = prevalence
            )
          }
        )
      )
    }
  )

  unique_patient_data <- prediction_data |>
    dplyr::distinct(
      Row_ID,
      Y_true
    )

  prevalence <- mean(
    unique_patient_data$Y_true
  )

  reference_rows <- dplyr::bind_rows(
    lapply(
      thresholds,
      function(threshold) {
        tibble::tibble(
          Model = c(
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
    list(reference_rows)
  )
}

# ------------------------------------------------
# 11. Plotting
# ------------------------------------------------

plot_metric_summary <- function(overall_metrics) {
  plot_data <- overall_metrics |>
    dplyr::select(
      Model,
      auc,
      auprc,
      brier,
      logloss,
      mean_uncertainty,
      selective_risk
    ) |>
    tidyr::pivot_longer(
      cols = -Model,
      names_to = "Metric",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        auc = "AUC",
        auprc = "AUPRC",
        brier = "Brier score",
        logloss = "Log loss",
        mean_uncertainty = "Mean uncertainty",
        selective_risk = "Selective risk"
      )
    )

  plot_object <- ggplot(
    plot_data,
    aes(
      x = Model,
      y = Value,
      fill = Model
    )
  ) +
    geom_col(
      width = 0.65,
      alpha = 0.9
    ) +
    facet_wrap(
      ~ Metric,
      scales = "free_y",
      ncol = 3
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": B-RISA Model Summary"
      ),
      subtitle = paste0(
        N_REPEATS,
        " repeated stratified ",
        N_FOLDS,
        "-fold cross-validation"
      ),
      x = "Model",
      y = "Metric value",
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(
        DATASET_ID,
        "_overall_metric_summary"
      )
    ),
    width = 12,
    height = 8
  )

  print(plot_object)
}

plot_roc_curves <- function(prediction_data) {
  curve_rows <- list()

  for (model_name in unique(prediction_data$Model)) {
    model_data <- prediction_data[
      prediction_data$Model == model_name,
      ,
      drop = FALSE
    ]

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
        False_Positive_Rate =
          1 - roc_object$specificities,
        True_Positive_Rate =
          roc_object$sensitivities
      )
  }

  curve_data <- dplyr::bind_rows(
    curve_rows
  )

  plot_object <- ggplot(
    curve_data,
    aes(
      x = False_Positive_Rate,
      y = True_Positive_Rate,
      color = Model
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray50"
    ) +
    scale_color_manual(values = model_colors) +
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": ROC Curves"
      ),
      x = "False positive rate",
      y = "True positive rate",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(
        DATASET_ID,
        "_roc_curves"
      )
    ),
    width = 8,
    height = 7
  )

  print(plot_object)
}

plot_precision_recall_curves <- function(
  prediction_data
) {
  curve_rows <- list()

  for (model_name in unique(prediction_data$Model)) {
    model_data <- prediction_data[
      prediction_data$Model == model_name,
      ,
      drop = FALSE
    ]

    thresholds <- sort(
      unique(model_data$p_mean),
      decreasing = TRUE
    )

    model_curve <- dplyr::bind_rows(
      lapply(
        thresholds,
        function(threshold) {
          predicted_positive <-
            model_data$p_mean >= threshold

          true_positive <- sum(
            predicted_positive &
              model_data$Y_true == 1L
          )

          false_positive <- sum(
            predicted_positive &
              model_data$Y_true == 0L
          )

          false_negative <- sum(
            !predicted_positive &
              model_data$Y_true == 1L
          )

          tibble::tibble(
            Model = model_name,
            Threshold = threshold,
            Precision = if (
              true_positive +
                false_positive >
                0L
            ) {
              true_positive /
                (
                  true_positive +
                    false_positive
                )
            } else {
              NA_real_
            },
            Recall = if (
              true_positive +
                false_negative >
                0L
            ) {
              true_positive /
                (
                  true_positive +
                    false_negative
                )
            } else {
              NA_real_
            }
          )
        }
      )
    ) |>
      tidyr::drop_na(
        Precision,
        Recall
      )

    curve_rows[[length(curve_rows) + 1L]] <-
      model_curve
  }

  curve_data <- dplyr::bind_rows(
    curve_rows
  )

  plot_object <- ggplot(
    curve_data,
    aes(
      x = Recall,
      y = Precision,
      color = Model
    )
  ) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = model_colors) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Precision-Recall Curves"
      ),
      x = "Recall / Sensitivity",
      y = "Precision / PPV",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$performance_figs,
      paste0(
        DATASET_ID,
        "_precision_recall_curves"
      )
    ),
    width = 8,
    height = 7
  )

  print(plot_object)
}

plot_calibration <- function(prediction_data) {
  calibration_data <- prediction_data |>
    dplyr::mutate(
      Bin = cut(
        p_mean,
        breaks = seq(
          0,
          1,
          length.out =
            CALIBRATION_BINS + 1L
        ),
        include.lowest = TRUE,
        labels = FALSE
      )
    ) |>
    dplyr::group_by(
      Model,
      Bin
    ) |>
    dplyr::summarise(
      Mean_Predicted =
        mean(p_mean),
      Observed_Rate =
        mean(Y_true),
      n = dplyr::n(),
      .groups = "drop"
    )

  safe_fwrite(
    calibration_data,
    file.path(
      dirs$summaries,
      paste0(
        DATASET_ID,
        "_calibration_bins.csv"
      )
    )
  )

  plot_object <- ggplot(
    calibration_data,
    aes(
      x = Mean_Predicted,
      y = Observed_Rate,
      color = Model
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray45"
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(
      aes(size = n),
      alpha = 0.9
    ) +
    scale_color_manual(values = model_colors) +
    scale_size_continuous(
      range = c(2, 6)
    ) +
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Calibration Curves"
      ),
      x = "Mean predicted probability",
      y = "Observed event rate",
      color = "Model",
      size = "Bin n"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$calibration_figs,
      paste0(
        DATASET_ID,
        "_calibration_curves"
      )
    ),
    width = 8,
    height = 7
  )

  print(plot_object)
}

plot_uncertainty <- function(prediction_data) {
  uncertainty_plot <- ggplot(
    prediction_data,
    aes(
      x = Model,
      y = p_sd,
      fill = Model
    )
  ) +
    geom_violin(
      alpha = 0.75,
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.alpha = 0.25
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Posterior Predictive Uncertainty"
      ),
      x = "Model",
      y = "Posterior predictive SD",
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    uncertainty_plot,
    file.path(
      dirs$uncertainty_figs,
      paste0(
        DATASET_ID,
        "_posterior_predictive_uncertainty"
      )
    ),
    width = 8,
    height = 7
  )

  print(uncertainty_plot)

  interval_data <- prediction_data |>
    dplyr::mutate(
      interval_width =
        p_upr - p_lwr
    )

  interval_plot <- ggplot(
    interval_data,
    aes(
      x = Model,
      y = interval_width,
      fill = Model
    )
  ) +
    geom_violin(
      alpha = 0.75,
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.alpha = 0.25
    ) +
    scale_fill_manual(values = model_colors) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Posterior Predictive Interval Width"
      ),
      x = "Model",
      y = "95% posterior interval width",
      fill = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    interval_plot,
    file.path(
      dirs$uncertainty_figs,
      paste0(
        DATASET_ID,
        "_posterior_interval_width"
      )
    ),
    width = 8,
    height = 7
  )

  print(interval_plot)
}

plot_prediction_density <- function(
  prediction_data
) {
  plot_object <- ggplot(
    prediction_data,
    aes(
      x = p_mean,
      fill = Outcome_true,
      color = Outcome_true
    )
  ) +
    geom_density(
      alpha = 0.25,
      linewidth = 1.1
    ) +
    facet_wrap(
      ~ Model,
      ncol = 1
    ) +
    scale_color_manual(values = outcome_colors) +
    scale_fill_manual(values = outcome_colors) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Predicted Probability Distributions"
      ),
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
      paste0(
        DATASET_ID,
        "_predicted_probability_density"
      )
    ),
    width = 8,
    height = 9
  )

  print(plot_object)
}

plot_selective_risk <- function(
  prediction_data
) {
  coverage_grid <- seq(
    0.50,
    1.00,
    by = 0.05
  )

  curve_rows <- list()

  for (model_name in unique(prediction_data$Model)) {
    model_data <- prediction_data[
      prediction_data$Model == model_name,
      ,
      drop = FALSE
    ]

    predicted_class <- as.integer(
      model_data$p_mean >=
        CLASSIFICATION_THRESHOLD
    )

    for (target_coverage in coverage_grid) {
      retained_index <- select_least_uncertain(
        model_data$p_sd,
        target_coverage
      )

      curve_rows[[length(curve_rows) + 1L]] <-
        tibble::tibble(
          Model = model_name,
          Target_Coverage =
            target_coverage,
          Coverage =
            length(retained_index) /
            nrow(model_data),
          Selective_Risk = mean(
            predicted_class[retained_index] !=
              model_data$Y_true[retained_index]
          )
        )
    }
  }

  curve_data <- dplyr::bind_rows(
    curve_rows
  )

  safe_fwrite(
    curve_data,
    file.path(
      dirs$summaries,
      paste0(
        DATASET_ID,
        "_risk_coverage_curve.csv"
      )
    )
  )

  plot_object <- ggplot(
    curve_data,
    aes(
      x = Coverage,
      y = Selective_Risk,
      color = Model
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.2) +
    scale_color_manual(values = model_colors) +
    scale_x_continuous(
      labels =
        scales::percent_format(accuracy = 1)
    ) +
    scale_y_continuous(
      labels =
        scales::percent_format(accuracy = 1)
    ) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Selective Prediction Risk-Coverage Curve"
      ),
      subtitle = "Lower risk is better",
      x = "Coverage retained",
      y = "Selective risk",
      color = "Model"
    ) +
    theme_b_risa()

  save_plot_pdf_png(
    plot_object,
    file.path(
      dirs$selective_figs,
      paste0(
        DATASET_ID,
        "_selective_risk_coverage_curve"
      )
    ),
    width = 8,
    height = 7
  )

  print(plot_object)
}

plot_decision_curves <- function(dca_data) {
  full_plot <- ggplot(
    dca_data,
    aes(
      x = Threshold,
      y = Net_Benefit,
      color = Model,
      linetype = Model
    )
  ) +
    geom_line(linewidth = 1.15) +
    scale_color_manual(values = model_colors) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Decision Curve Analysis"
      ),
      subtitle = "Net benefit across threshold probabilities",
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
      paste0(
        DATASET_ID,
        "_decision_curve_net_benefit"
      )
    ),
    width = 9,
    height = 7
  )

  print(full_plot)

  zoom_data <- dca_data |>
    dplyr::filter(
      Threshold >= 0.05,
      Threshold <= 0.80
    )

  zoom_plot <- ggplot(
    zoom_data,
    aes(
      x = Threshold,
      y = Net_Benefit,
      color = Model,
      linetype = Model
    )
  ) +
    geom_line(linewidth = 1.15) +
    scale_color_manual(values = model_colors) +
    labs(
      title = paste0(
        DATASET_LABEL,
        ": Decision Curve Analysis, Clinical Threshold Range"
      ),
      subtitle = "Threshold probabilities 0.05 to 0.80",
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
      paste0(
        DATASET_ID,
        "_decision_curve_net_benefit_zoom"
      )
    ),
    width = 9,
    height = 7
  )

  print(zoom_plot)
}

# ------------------------------------------------
# 12. Run the repeated cross-validation
# ------------------------------------------------

rpa_standardized_data <- load_rpa_standardized_dataset(
  RPA_STD_CLEAN_FILE
)

dataset_summary <- tibble::tibble(
  Dataset = DATASET_LABEL,
  n = nrow(rpa_standardized_data),
  p = length(RPA_STD_FEATURES),
  Features = paste(
    RPA_STD_FEATURES,
    collapse = ", "
  ),
  n_no_appendicitis =
    sum(rpa_standardized_data$Y == 0L),
  n_appendicitis =
    sum(rpa_standardized_data$Y == 1L),
  prevalence_appendicitis =
    mean(rpa_standardized_data$Y),
  n_folds = N_FOLDS,
  n_repeats = N_REPEATS,
  model_set = paste(
    MODEL_SET,
    collapse = ", "
  )
)

safe_fwrite(
  dataset_summary,
  file.path(
    dirs$summaries,
    paste0(
      DATASET_ID,
      "_dataset_summary.csv"
    )
  )
)

cat("\n================================================\n")
cat("B-RISA STANDARDIZED RPA ANALYSIS\n")
cat("================================================\n")
print(dataset_summary)

folds <- make_repeated_stratified_folds(
  y = rpa_standardized_data$Y,
  k = N_FOLDS,
  repeats = N_REPEATS,
  seed = SEED
)

prediction_list <- vector(
  mode = "list",
  length = length(folds) *
    length(MODEL_SET)
)

fold_metric_list <- vector(
  mode = "list",
  length = length(prediction_list)
)

status_list <- vector(
  mode = "list",
  length = length(prediction_list)
)

list_index <- 0L
total_fits <- length(prediction_list)

for (fold_object in folds) {
  repeat_id <- fold_object$repeat_id
  fold_id <- fold_object$fold_id

  training_data <- rpa_standardized_data[
    fold_object$train_index,
    ,
    drop = FALSE
  ]

  test_data <- rpa_standardized_data[
    fold_object$test_index,
    ,
    drop = FALSE
  ]

  preprocessing <- fit_preprocess(
    training_data
  )

  x_train <- apply_preprocess(
    training_data,
    preprocessing
  )

  x_test <- apply_preprocess(
    test_data,
    preprocessing
  )

  y_train <- training_data$Y
  y_test <- test_data$Y

  for (model_name in MODEL_SET) {
    list_index <- list_index + 1L

    cat(
      sprintf(
        "\rRepeat %d/%d | Fold %d/%d | Model %s | Fit %d/%d",
        repeat_id,
        N_REPEATS,
        fold_id,
        N_FOLDS,
        model_name,
        list_index,
        total_fits
      )
    )

    fit_seed <- make_fit_seed(
      repeat_id,
      fold_id,
      model_name
    )

    set.seed(fit_seed)
    start_time <- Sys.time()

    fit_result <- tryCatch(
      fit_model(
        model_name,
        x_train,
        y_train,
        x_test
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
        Repeat = repeat_id,
        Fold = fold_id,
        Model = model_name,
        Fit_seed = fit_seed,
        Status = "FAILED",
        Failure_message =
          fit_result$message,
        Runtime_seconds =
          runtime_seconds
      )

      status_list[[list_index]] <- failure_row

      if (STOP_ON_MODEL_FAILURE) {
        safe_fwrite(
          dplyr::bind_rows(
            status_list[
              !vapply(
                status_list,
                is.null,
                logical(1)
              )
            ]
          ),
          file.path(
            dirs$status,
            paste0(
              DATASET_ID,
              "_model_status_partial.csv"
            )
          )
        )

        stop(
          "\nModel failure: ",
          model_name,
          ", repeat ",
          repeat_id,
          ", fold ",
          fold_id,
          ". ",
          fit_result$message
        )
      }

      next
    }

    posterior_draws <- validate_prediction_draws(
      draws = fit_result$draws,
      expected_draws =
        expected_draw_count(model_name),
      expected_test_size =
        nrow(x_test),
      model_name = model_name
    )

    posterior_summary <- summarise_draws(
      posterior_draws
    )

    prediction_sd_across_patients <- if (
      nrow(posterior_summary) >= 2L
    ) {
      stats::sd(
        posterior_summary$p_mean
      )
    } else {
      NA_real_
    }

    number_unique_predictions <- dplyr::n_distinct(
      round(
        posterior_summary$p_mean,
        digits = 12
      )
    )

    degenerate_prediction <- (
      nrow(posterior_summary) > 1L &&
        (
          !is.finite(
            prediction_sd_across_patients
          ) ||
            prediction_sd_across_patients <=
            DEGENERATE_PREDICTION_TOL ||
            number_unique_predictions <= 1L
        )
    )

    prediction_data <- tibble::tibble(
      Dataset = DATASET_LABEL,
      Repeat = repeat_id,
      Fold = fold_id,
      Row_ID = test_data$Row_ID,
      Model = model_name,
      Y_true = y_test,
      Outcome_true = factor(
        y_test,
        levels = c(0, 1),
        labels = c(
          "No_appendicitis",
          "Appendicitis"
        )
      ),
      p_mean = posterior_summary$p_mean,
      p_sd = posterior_summary$p_sd,
      p_lwr = posterior_summary$p_lwr,
      p_upr = posterior_summary$p_upr
    )

    fold_metrics <- calc_metrics(
      y = prediction_data$Y_true,
      p_mean = prediction_data$p_mean,
      p_sd = prediction_data$p_sd,
      p_lwr = prediction_data$p_lwr,
      p_upr = prediction_data$p_upr
    ) |>
      dplyr::mutate(
        Dataset = DATASET_LABEL,
        Repeat = repeat_id,
        Fold = fold_id,
        Model = model_name,
        .before = 1
      )

    status_base <- tibble::tibble(
      Dataset = DATASET_LABEL,
      Repeat = repeat_id,
      Fold = fold_id,
      Model = model_name,
      Fit_seed = fit_seed,
      Status = if (
        degenerate_prediction
      ) {
        "DEGENERATE"
      } else {
        "OK"
      },
      Failure_message = NA_character_,
      Runtime_seconds =
        runtime_seconds,
      Posterior_draw_rows =
        nrow(posterior_draws),
      Posterior_draw_columns =
        ncol(posterior_draws),
      Probability_min =
        min(prediction_data$p_mean),
      Probability_max =
        max(prediction_data$p_mean),
      Probability_mean =
        mean(prediction_data$p_mean),
      Probability_sd_across_patients =
        prediction_sd_across_patients,
      Number_unique_predictions =
        number_unique_predictions,
      Number_predicted_positive =
        sum(
          prediction_data$p_mean >=
            CLASSIFICATION_THRESHOLD
        ),
      Degenerate_prediction =
        degenerate_prediction
    )

    diagnostic_row <- tibble::as_tibble(
      fit_result$diagnostics
    )

    status_list[[list_index]] <- dplyr::bind_cols(
      status_base,
      diagnostic_row
    )

    prediction_list[[list_index]] <- prediction_data
    fold_metric_list[[list_index]] <- fold_metrics

    if (
      degenerate_prediction &&
      STOP_ON_DEGENERATE_PREDICTIONS
    ) {
      safe_fwrite(
        dplyr::bind_rows(
          status_list[
            !vapply(
              status_list,
              is.null,
              logical(1)
            )
          ]
        ),
        file.path(
          dirs$status,
          paste0(
            DATASET_ID,
            "_model_status_partial.csv"
          )
        )
      )

      stop(
        "\nDegenerate patient-level predictions detected for ",
        model_name,
        ", repeat ",
        repeat_id,
        ", fold ",
        fold_id,
        ". Prediction SD = ",
        signif(
          prediction_sd_across_patients,
          6
        ),
        "."
      )
    }

    rm(
      fit_result,
      posterior_draws,
      posterior_summary,
      prediction_data,
      fold_metrics,
      status_base,
      diagnostic_row
    )

    gc(FALSE)
  }

  rm(
    training_data,
    test_data,
    preprocessing,
    x_train,
    x_test,
    y_train,
    y_test
  )

  gc(FALSE)
}

cat("\n")

all_predictions <- dplyr::bind_rows(
  prediction_list
)

fold_level_metrics <- dplyr::bind_rows(
  fold_metric_list
)

model_status <- dplyr::bind_rows(
  status_list
)

# Verify exactly one out-of-fold prediction per patient, model, and repeat.
prediction_completeness <- all_predictions |>
  dplyr::count(
    Repeat,
    Model,
    Row_ID,
    name = "Prediction_count"
  )

if (any(prediction_completeness$Prediction_count != 1L)) {
  stop(
    "Each patient must have exactly one out-of-fold prediction ",
    "per model and repeat."
  )
}

expected_prediction_rows <- (
  nrow(rpa_standardized_data) *
    N_REPEATS *
    length(MODEL_SET)
)

if (nrow(all_predictions) != expected_prediction_rows) {
  stop(
    "Expected ",
    expected_prediction_rows,
    " out-of-fold prediction rows but found ",
    nrow(all_predictions),
    "."
  )
}

# Correct repeat-level metrics:
# pool all 10 out-of-fold test folds within a repeat before calculating
# AUC, AUPRC, threshold metrics, calibration, and proper scores.
repeat_level_metrics <- all_predictions |>
  dplyr::group_by(
    Dataset,
    Repeat,
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

aggregated_predictions <- aggregate_repeated_predictions(
  all_predictions
) |>
  dplyr::mutate(
    Dataset = DATASET_LABEL,
    .before = 1
  )

overall_metrics <- aggregated_predictions |>
  dplyr::group_by(
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

numeric_metric_columns <- setdiff(
  names(repeat_level_metrics),
  c(
    "Dataset",
    "Repeat",
    "Model"
  )
)

repeat_metric_summary <- repeat_level_metrics |>
  dplyr::group_by(
    Dataset,
    Model
  ) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(
        numeric_metric_columns
      ),
      list(
        mean = ~ {
          values <- .x[is.finite(.x)]
          if (
            length(values) == 0L
          ) {
            NA_real_
          } else {
            mean(values)
          }
        },
        sd = ~ {
          values <- .x[is.finite(.x)]
          if (
            length(values) < 2L
          ) {
            NA_real_
          } else {
            stats::sd(values)
          }
        }
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

decision_curve_data <- decision_curve_table(
  aggregated_predictions,
  DCA_THRESHOLDS
) |>
  dplyr::mutate(
    Dataset = DATASET_LABEL,
    .before = 1
  )

standardized_decision_curve_data <-
  decision_curve_data |>
  dplyr::mutate(
    Standardized_Net_Benefit =
      Net_Benefit / Prevalence
  )

# ------------------------------------------------
# 13. Save outputs
# ------------------------------------------------

safe_fwrite(
  fold_level_metrics,
  file.path(
    dirs$metrics,
    paste0(
      DATASET_ID,
      "_fold_level_metrics.csv"
    )
  )
)

safe_fwrite(
  repeat_level_metrics,
  file.path(
    dirs$metrics,
    paste0(
      DATASET_ID,
      "_repeat_level_metrics.csv"
    )
  )
)

safe_fwrite(
  overall_metrics,
  file.path(
    dirs$metrics,
    paste0(
      DATASET_ID,
      "_overall_metrics_from_cv_predictions.csv"
    )
  )
)

safe_fwrite(
  repeat_metric_summary,
  file.path(
    dirs$summaries,
    paste0(
      DATASET_ID,
      "_repeat_level_metric_summary.csv"
    )
  )
)

safe_fwrite(
  all_predictions,
  file.path(
    dirs$predictions,
    paste0(
      DATASET_ID,
      "_all_cv_predictions.csv"
    )
  )
)

safe_fwrite(
  aggregated_predictions,
  file.path(
    dirs$predictions,
    paste0(
      DATASET_ID,
      "_aggregated_cv_predictions.csv"
    )
  )
)

safe_fwrite(
  prediction_completeness,
  file.path(
    dirs$status,
    paste0(
      DATASET_ID,
      "_prediction_completeness.csv"
    )
  )
)

safe_fwrite(
  model_status,
  file.path(
    dirs$status,
    paste0(
      DATASET_ID,
      "_model_status.csv"
    )
  )
)

blr_diagnostics <- model_status |>
  dplyr::filter(Model == "BLR")

safe_fwrite(
  blr_diagnostics,
  file.path(
    dirs$status,
    paste0(
      DATASET_ID,
      "_blr_diagnostics.csv"
    )
  )
)

safe_fwrite(
  decision_curve_data,
  file.path(
    dirs$dca,
    paste0(
      DATASET_ID,
      "_decision_curve_net_benefit.csv"
    )
  )
)

safe_fwrite(
  standardized_decision_curve_data,
  file.path(
    dirs$dca,
    paste0(
      DATASET_ID,
      "_decision_curve_standardized_net_benefit.csv"
    )
  )
)

method_manifest <- tibble::tibble(
  Item = c(
    "Dataset",
    "Predictors",
    "Input scale",
    "Cross-validation",
    "Models",
    "BLR method",
    "BLR intercept prior SD",
    "BLR coefficient prior SD",
    "BLR optimizer tolerance",
    "BLDA method",
    "BLDA prior hyperparameters",
    "BART method",
    "BART probability transformation",
    "BART posterior draws",
    "BART burn-in",
    "BART number of trees",
    "Posterior draw validation",
    "Degenerate prediction rule",
    "Repeat-level metric calculation",
    "Repeated prediction aggregation",
    "Random seed"
  ),
  Value = c(
    DATASET_LABEL,
    paste(RPA_STD_FEATURES, collapse = ", "),
    paste0(
      "Pre-standardized input; median imputation and ",
      "training-fold re-centering/re-scaling applied inside CV"
    ),
    paste0(
      N_REPEATS,
      " repeated ",
      N_FOLDS,
      "-fold stratified CV"
    ),
    paste(MODEL_SET, collapse = ", "),
    paste0(
      "MAP plus Laplace approximation; ",
      "training-fold standardization"
    ),
    BLR_INTERCEPT_PRIOR_SD,
    BLR_COEFFICIENT_PRIOR_SD,
    BLR_OPTIM_RELTOL,
    "Normal-Inverse-Wishart shared-covariance Bayesian LDA",
    paste0(
      N_POST_DRAWS_BLDA, " posterior draws; kappa0=",
      BLDA_MEAN_PRIOR_KAPPA,
      "; nu0=p+", BLDA_IW_DF_OFFSET,
      "; prior covariance scale=", BLDA_PRIOR_COVARIANCE_SCALE,
      "; class prevalence prior=Beta(",
      BLDA_CLASS_PRIOR_ALPHA, ",", BLDA_CLASS_PRIOR_BETA, ")"
    ),
    "Probit Bayesian additive regression trees",
    "pnorm applied to latent yhat.test draws",
    N_POST_DRAWS_BART,
    N_BURN_BART,
    BART_NTREE,
    paste0(
      "Every matrix must equal posterior draws x test patients; ",
      "all entries finite and in [0,1]"
    ),
    paste0(
      "Stop if patient-level posterior mean SD <= ",
      DEGENERATE_PREDICTION_TOL,
      " within a test fold"
    ),
    paste0(
      "All 10 out-of-fold test predictions are pooled within ",
      "each repeat before metrics are calculated"
    ),
    "First and second moments combined by the law of total variance",
    SEED
  )
)

safe_fwrite(
  method_manifest,
  file.path(
    RESULT_ROOT,
    "rpa_standardized_corrected_methodology_manifest.csv"
  )
)

saveRDS(
  list(
    dataset_summary =
      dataset_summary,
    fold_level_metrics =
      fold_level_metrics,
    repeat_level_metrics =
      repeat_level_metrics,
    repeat_metric_summary =
      repeat_metric_summary,
    all_predictions =
      all_predictions,
    aggregated_predictions =
      aggregated_predictions,
    overall_metrics =
      overall_metrics,
    model_status =
      model_status,
    decision_curve_data =
      decision_curve_data,
    standardized_decision_curve_data =
      standardized_decision_curve_data,
    method_manifest =
      method_manifest
  ),
  file.path(
    dirs$result_base,
    paste0(
      DATASET_ID,
      "_updated_b_risa_analysis.rds"
    )
  ),
  compress = "xz"
)

# ------------------------------------------------
# 14. Figures
# ------------------------------------------------

plot_metric_summary(
  overall_metrics
)

plot_roc_curves(
  aggregated_predictions
)

plot_precision_recall_curves(
  aggregated_predictions
)

plot_calibration(
  aggregated_predictions
)

plot_uncertainty(
  aggregated_predictions
)

plot_prediction_density(
  aggregated_predictions
)

plot_selective_risk(
  aggregated_predictions
)

plot_decision_curves(
  decision_curve_data
)

# ------------------------------------------------
# 15. Final validation and summary
# ------------------------------------------------

if (any(model_status$Status != "OK")) {
  stop(
    "The analysis completed with at least one non-OK model status."
  )
}

if (
  "BLR" %in% model_status$Model &&
  any(
    model_status$Model == "BLR" &
      model_status$Degenerate_prediction
  )
) {
  stop("At least one BLR fold has degenerate predictions.")
}

cat("\n================================================\n")
cat("B-RISA STANDARDIZED RPA ANALYSIS COMPLETE\n")
cat("================================================\n")

cat("\nOverall metrics:\n")
print(overall_metrics)

cat("\nAll model fits completed successfully:\n")
print(
  model_status |>
    dplyr::count(
      Model,
      Status
    )
)

cat("\nResults directory:\n")
cat(dirs$result_base, "\n")

cat("\nFigures directory:\n")
cat(dirs$figure_base, "\n")

cat("\nClinical stage complete. Statistical stage starts next.\n")
cat("================================================\n")

# =====================================================================
# END OF CLINICAL MODEL-FITTING STAGE
# The statistical-analysis stage below intentionally starts from a clean
# R session and reads the clinical outputs generated above. This prevents
# accidental dependence on in-memory objects and makes the inferential
# stage independently auditable.
# =====================================================================

# ================================================================
# B-RISA Standardized RPA Statistical Analysis
# FINAL BAYESIAN MODEL IMPLEMENTATION
#
# Quality-control note:
#   Repeat-level descriptive summaries and model-wise Shapiro-Wilk
#   tests now extract the complete vector of repeat values. A repeat-value count audit now stops the analysis
#   unless every model contributes exactly EXPECTED_REPEATS values.
#
# Cohort:
#   Standardized Regensburg Pediatric Appendicitis dataset
#   Predictors: Body_Temperature, WBC_Count, and
#   Neutrophil_Percentage
#   Y = 1 appendicitis; Y = 0 no appendicitis
#
# Required inputs from:
#   rpa_standardized_updated_clinical_analysis.R
#
# Primary inputs:
#   1. rpa_standardized_repeat_level_metrics.csv
#   2. rpa_standardized_aggregated_cv_predictions.csv
#   3. rpa_standardized_fold_level_metrics.csv
#   4. rpa_standardized_model_status.csv
#
# Statistical strategy:
#
# A. Quality-control audit
#    - verifies 50 repeats and 10 folds;
#    - checks model failures and degenerate predictions;
#    - checks finite metric availability.
#
# B. Repeat-level stability analysis
#    - each repeat is calculated from all pooled out-of-fold
#      predictions, rather than the mean of fold-level AUCs;
#    - Shapiro-Wilk tests are applied to paired model differences;
#    - randomized-block ANOVA plus paired t-tests are used only when
#      all paired differences are compatible with normality;
#    - otherwise Friedman plus paired Wilcoxon tests are used;
#    - Holm correction is applied within each metric.
#
# C. Patient-level paired bootstrap
#    - uses the final aggregated out-of-fold predictions;
#    - resamples patients, preserving model pairing;
#    - reports percentile confidence intervals and pairwise
#      differences for clinically important metrics.
#
# Interpretation warning:
#   Cross-validation repeats reuse the same patients and are not
#   independent clinical datasets. Repeat-level p-values are therefore
#   reported as stability/sensitivity analyses. The paired patient
#   bootstrap is the primary inferential summary, but it conditions on
#   the already-fitted cross-validated predictions and does not fully
#   reproduce model-training uncertainty.
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
OUTPUT_ROOT <- file.path(PROJECT_ROOT, "analysis_outputs", "rpa_standardized_clinical")
RESULT_BASE <- file.path(OUTPUT_ROOT, "results", "rpa_standardized")
FIGURE_BASE <- file.path(OUTPUT_ROOT, "figures", "rpa_standardized")

REPEAT_METRIC_FILE <- file.path(
  RESULT_BASE, "metrics", "rpa_standardized_repeat_level_metrics.csv"
)
FOLD_METRIC_FILE <- file.path(
  RESULT_BASE, "metrics", "rpa_standardized_fold_level_metrics.csv"
)
AGGREGATED_PREDICTION_FILE <- file.path(
  RESULT_BASE, "predictions", "rpa_standardized_aggregated_cv_predictions.csv"
)
MODEL_STATUS_FILE <- file.path(
  RESULT_BASE, "status", "rpa_standardized_model_status.csv"
)

OUTPUT_DIR <- file.path(
  RESULT_BASE, "statistical_analysis", "final_analysis"
)
FIGURE_DIR <- file.path(
  FIGURE_BASE, "statistical_analysis", "final_analysis"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  REPEAT_METRIC_FILE,
  FOLD_METRIC_FILE,
  AGGREGATED_PREDICTION_FILE,
  MODEL_STATUS_FILE
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required clinical-analysis output files are missing:\n",
    paste(missing_files, collapse = "\n"),
    "\nThe clinical-analysis stage above must complete successfully."
  )
}

# ------------------------------------------------
# 1. Analysis settings
# ------------------------------------------------

INSTALL_MISSING_PACKAGES <- FALSE

ALPHA <- 0.05
CONFIDENCE_LEVEL <- 0.95
P_ADJUST_METHOD <- "holm"

EXPECTED_REPEATS <- 50L
EXPECTED_FOLDS <- 10L
MIN_COMPLETE_REPEATS <- 40L

MODEL_ORDER <- c(
  "BLR",
  "BLDA",
  "BART"
)

# "auto", "parametric", or "nonparametric".
REPEAT_TEST_FAMILY <- "auto"
RUN_REPEAT_LEVEL_TESTS <- TRUE

RUN_PATIENT_BOOTSTRAP <- TRUE
BOOTSTRAP_REPS <- 2000L
BOOTSTRAP_SEED <- 20260619L

CLASSIFICATION_THRESHOLD <- 0.50
SELECTIVE_COVERAGE <- 0.80
CALIBRATION_BINS <- 10L
EPS <- 1e-6

STOP_IF_MODEL_STATUS_PROBLEM <- TRUE
DEGENERATE_PREDICTION_TOL <- 1e-12

# Formal performance metrics.
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

# These are summarized, but narrower uncertainty is not automatically
# better and calibration coefficients are unstable in tiny folds.
DESCRIPTIVE_ONLY_METRICS <- c(
  "calibration_intercept",
  "calibration_slope",
  "mean_uncertainty",
  "mean_interval_width",
  "selective_coverage"
)

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
  "BART" = "#6A4C93"
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
  "PRROC"
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

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(pROC)
  library(PRROC)
})

# ------------------------------------------------
# 3. Helper functions
# ------------------------------------------------

safe_fwrite <- function(x, path) {
  data.table::fwrite(
    as.data.frame(x),
    path,
    na = "NA"
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

clip_prob <- function(p, eps = EPS) {
  pmin(
    pmax(as.numeric(p), eps),
    1 - eps
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
          y,
          p,
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
  p <- clip_prob(p)

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

  ece <- 0

  for (bin_id in seq_len(n_bins)) {
    index <- which(bins == bin_id)

    if (length(index) == 0L) {
      next
    }

    ece <- ece +
      length(index) / length(y) *
      abs(
        mean(y[index]) -
          mean(p[index])
      )
  }

  ece
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

  keep_n <- max(
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

  order(uncertainty)[seq_len(keep_n)]
}

calculate_metric_vector <- function(
  y,
  p_mean,
  p_sd,
  p_lwr,
  p_upr
) {
  y <- as.integer(y)
  p_mean <- clip_prob(p_mean)
  p_sd <- as.numeric(p_sd)
  p_lwr <- as.numeric(p_lwr)
  p_upr <- as.numeric(p_upr)

  predicted_class <- as.integer(
    p_mean >= CLASSIFICATION_THRESHOLD
  )

  tp <- sum(
    predicted_class == 1L &
      y == 1L
  )

  tn <- sum(
    predicted_class == 0L &
      y == 0L
  )

  fp <- sum(
    predicted_class == 1L &
      y == 0L
  )

  fn <- sum(
    predicted_class == 0L &
      y == 1L
  )

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
    SELECTIVE_COVERAGE
  )

  c(
    accuracy =
      mean(predicted_class == y),
    sensitivity =
      sensitivity,
    specificity =
      specificity,
    ppv =
      ppv,
    npv =
      npv,
    f1 =
      f1,
    auc =
      safe_auc(y, p_mean),
    auprc =
      safe_auprc(y, p_mean),
    brier =
      mean((p_mean - y)^2),
    logloss =
      -mean(
        y * log(p_mean) +
          (1 - y) *
          log(1 - p_mean)
      ),
    ece =
      calc_ece(y, p_mean),
    mean_uncertainty =
      mean(p_sd),
    mean_interval_width =
      mean(p_upr - p_lwr),
    selective_coverage =
      length(retained_index) /
      length(y),
    selective_risk =
      mean(
        predicted_class[retained_index] !=
          y[retained_index]
      )
  )
}

safe_descriptive_summary <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(
      tibble::tibble(
        n = 0L,
        mean = NA_real_,
        sd = NA_real_,
        standard_error = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        median = NA_real_,
        q1 = NA_real_,
        q3 = NA_real_,
        minimum = NA_real_,
        maximum = NA_real_
      )
    )
  }

  mean_value <- mean(x)
  sd_value <- if (n >= 2L) {
    stats::sd(x)
  } else {
    NA_real_
  }

  standard_error <- if (n >= 2L) {
    sd_value / sqrt(n)
  } else {
    NA_real_
  }

  if (n >= 2L) {
    critical_value <- stats::qt(
      1 - (1 - CONFIDENCE_LEVEL) / 2,
      df = n - 1L
    )

    ci_lower <- mean_value -
      critical_value *
      standard_error

    ci_upper <- mean_value +
      critical_value *
      standard_error
  } else {
    ci_lower <- NA_real_
    ci_upper <- NA_real_
  }

  tibble::tibble(
    n = n,
    mean = mean_value,
    sd = sd_value,
    standard_error = standard_error,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    median = stats::median(x),
    q1 = as.numeric(
      stats::quantile(
        x,
        0.25,
        names = FALSE
      )
    ),
    q3 = as.numeric(
      stats::quantile(
        x,
        0.75,
        names = FALSE
      )
    ),
    minimum = min(x),
    maximum = max(x)
  )
}


extract_metric_values <- function(
  data,
  metric_name,
  model_name,
  expected_repeats = EXPECTED_REPEATS
) {
  if (!(metric_name %in% names(data))) {
    stop(
      "Metric column is missing: ",
      metric_name
    )
  }

  if (!("Model" %in% names(data))) {
    stop("The repeat-level data do not contain a Model column.")
  }

  model_index <- as.character(data$Model) == model_name

  values <- data[[metric_name]][model_index]
  values <- suppressWarnings(as.numeric(values))

  if (length(values) != expected_repeats) {
    stop(
      "Expected ",
      expected_repeats,
      " repeat-level values for model ",
      model_name,
      " and metric ",
      metric_name,
      ", but found ",
      length(values),
      "."
    )
  }

  values
}

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n < 3L) {
    return(
      tibble::tibble(
        n = n,
        W = NA_real_,
        p_value = NA_real_,
        Decision = "Indeterminate",
        Reason = "Fewer than 3 finite observations"
      )
    )
  }

  if (
    !is.finite(stats::sd(x)) ||
    stats::sd(x) == 0
  ) {
    return(
      tibble::tibble(
        n = n,
        W = NA_real_,
        p_value = NA_real_,
        Decision = "Indeterminate",
        Reason = "All observations are identical"
      )
    )
  }

  test <- tryCatch(
    stats::shapiro.test(x),
    error = function(e) NULL
  )

  if (is.null(test)) {
    return(
      tibble::tibble(
        n = n,
        W = NA_real_,
        p_value = NA_real_,
        Decision = "Indeterminate",
        Reason = "Shapiro-Wilk test failed"
      )
    )
  }

  p_value <- unname(test$p.value)

  tibble::tibble(
    n = n,
    W = unname(test$statistic),
    p_value = p_value,
    Decision = if (
      p_value > ALPHA
    ) {
      "Fail to reject normality"
    } else {
      "Reject normality"
    },
    Reason = if (
      p_value > ALPHA
    ) {
      paste0("p > ", ALPHA)
    } else {
      paste0("p <= ", ALPHA)
    }
  )
}

cohens_dz <- function(difference) {
  difference <- difference[
    is.finite(difference)
  ]

  if (length(difference) < 2L) {
    return(NA_real_)
  }

  difference_sd <- stats::sd(difference)

  if (
    !is.finite(difference_sd) ||
    difference_sd == 0
  ) {
    return(NA_real_)
  }

  mean(difference) /
    difference_sd
}

rank_biserial_paired <- function(difference) {
  difference <- difference[
    is.finite(difference) &
      difference != 0
  ]

  if (length(difference) == 0L) {
    return(0)
  }

  absolute_ranks <- rank(
    abs(difference),
    ties.method = "average"
  )

  positive_sum <- sum(
    absolute_ranks[
      difference > 0
    ]
  )

  negative_sum <- sum(
    absolute_ranks[
      difference < 0
    ]
  )

  denominator <- positive_sum +
    negative_sum

  if (denominator == 0) {
    return(NA_real_)
  }

  (
    positive_sum -
      negative_sum
  ) / denominator
}

metric_direction <- function(metric_name) {
  if (metric_name %in% HIGHER_IS_BETTER) {
    return("Higher is better")
  }

  if (metric_name %in% LOWER_IS_BETTER) {
    return("Lower is better")
  }

  "Descriptive; smaller is not automatically better"
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
    if (difference > 0) {
      model_1
    } else {
      model_2
    }
  } else if (metric_name %in% LOWER_IS_BETTER) {
    if (difference < 0) {
      model_1
    } else {
      model_2
    }
  } else {
    "Interpret descriptively"
  }
}

# ------------------------------------------------
# 4. Read data
# ------------------------------------------------

repeat_metrics <- data.table::fread(
  REPEAT_METRIC_FILE,
  data.table = FALSE
)

fold_metrics <- data.table::fread(
  FOLD_METRIC_FILE,
  data.table = FALSE
)

aggregated_predictions <- data.table::fread(
  AGGREGATED_PREDICTION_FILE,
  data.table = FALSE
)

model_status <- data.table::fread(
  MODEL_STATUS_FILE,
  data.table = FALSE
)

for (data_name in c(
  "repeat_metrics",
  "fold_metrics",
  "aggregated_predictions",
  "model_status"
)) {
  data_object <- get(data_name)

  if (!("Model" %in% names(data_object))) {
    stop(data_name, " does not contain a Model column.")
  }

  data_object$Model <- factor(
    as.character(data_object$Model),
    levels = MODEL_ORDER
  )

  data_object <- data_object[
    !is.na(data_object$Model),
    ,
    drop = FALSE
  ]

  assign(data_name, data_object)
}

repeat_metrics$Repeat <- as.integer(
  repeat_metrics$Repeat
)

fold_metrics$Repeat <- as.integer(
  fold_metrics$Repeat
)

fold_metrics$Fold <- as.integer(
  fold_metrics$Fold
)

aggregated_predictions$Y_true <- as.integer(
  aggregated_predictions$Y_true
)

# ------------------------------------------------
# 5. Model-status and completeness audit
# ------------------------------------------------

if (!("Degenerate_prediction" %in% names(model_status))) {
  model_status$Degenerate_prediction <- FALSE
}

if (!("Probability_sd_across_patients" %in% names(model_status))) {
  model_status$Probability_sd_across_patients <- NA_real_
}

model_status$Degenerate_prediction <- as.logical(
  model_status$Degenerate_prediction
)

status_audit <- model_status |>
  dplyr::group_by(Model) |>
  dplyr::summarise(
    Total_fits = dplyr::n(),
    OK_fits = sum(Status == "OK", na.rm = TRUE),
    Failed_fits = sum(Status == "FAILED", na.rm = TRUE),
    Degenerate_fits = sum(
      Status == "DEGENERATE" |
        Degenerate_prediction,
      na.rm = TRUE
    ),
    Minimum_prediction_SD = {
      finite_sd <- Probability_sd_across_patients[
        is.finite(Probability_sd_across_patients)
      ]

      if (length(finite_sd) == 0L) {
        NA_real_
      } else {
        min(finite_sd)
      }
    },
    Maximum_prediction_SD = {
      finite_sd <- Probability_sd_across_patients[
        is.finite(Probability_sd_across_patients)
      ]

      if (length(finite_sd) == 0L) {
        NA_real_
      } else {
        max(finite_sd)
      }
    },
    .groups = "drop"
  )

safe_fwrite(
  status_audit,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_model_status_audit.csv"
  )
)

status_problem <- any(
  model_status$Status != "OK",
  na.rm = TRUE
)

if (
  "Probability_sd_across_patients" %in%
    names(model_status)
) {
  status_problem <- status_problem ||
    any(
      model_status$
        Probability_sd_across_patients <=
        DEGENERATE_PREDICTION_TOL,
      na.rm = TRUE
    )
}

if (
  status_problem &&
  STOP_IF_MODEL_STATUS_PROBLEM
) {
  stop(
    "The model-status audit found failed or degenerate fits. ",
    "Inspect rpa_standardized_model_status_audit.csv before inference."
  )
}

fold_completeness <- fold_metrics |>
  dplyr::count(
    Repeat,
    Model,
    name = "Observed_folds"
  ) |>
  dplyr::mutate(
    Expected_folds = EXPECTED_FOLDS,
    Complete = (
      Observed_folds ==
        Expected_folds
    )
  )

safe_fwrite(
  fold_completeness,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_fold_completeness.csv"
  )
)

repeat_completeness <- repeat_metrics |>
  dplyr::count(
    Repeat,
    Model,
    name = "Observed_rows"
  ) |>
  dplyr::mutate(
    Expected_rows = 1L,
    Complete = (
      Observed_rows == 1L
    )
  )

safe_fwrite(
  repeat_completeness,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_repeat_completeness.csv"
  )
)

if (
  any(!fold_completeness$Complete) ||
  any(!repeat_completeness$Complete)
) {
  stop(
    "Fold or repeat completeness failed. ",
    "Inspect the completeness files."
  )
}

detected_repeats <- dplyr::n_distinct(
  repeat_metrics$Repeat
)

detected_folds <- dplyr::n_distinct(
  fold_metrics$Fold
)

if (detected_repeats != EXPECTED_REPEATS) {
  warning(
    "Expected ",
    EXPECTED_REPEATS,
    " repeats but detected ",
    detected_repeats,
    "."
  )
}

if (detected_folds != EXPECTED_FOLDS) {
  warning(
    "Expected ",
    EXPECTED_FOLDS,
    " folds but detected ",
    detected_folds,
    "."
  )
}

# ------------------------------------------------
# 6. Metric eligibility and missingness
# ------------------------------------------------

available_inferential_metrics <- intersect(
  INFERENTIAL_METRICS,
  names(repeat_metrics)
)

available_descriptive_metrics <- intersect(
  DESCRIPTIVE_ONLY_METRICS,
  names(repeat_metrics)
)

all_available_metrics <- unique(
  c(
    available_inferential_metrics,
    available_descriptive_metrics
  )
)

if (length(available_inferential_metrics) == 0L) {
  stop("No inferential metrics were found.")
}

repeat_value_count_audit <- dplyr::bind_rows(
  lapply(
    all_available_metrics,
    function(metric_name) {
      dplyr::bind_rows(
        lapply(
          MODEL_ORDER,
          function(model_name) {
            model_index <- as.character(
              repeat_metrics$Model
            ) == model_name

            values <- repeat_metrics[[metric_name]][
              model_index
            ]

            tibble::tibble(
              Metric = metric_name,
              Model = model_name,
              Expected_rows = EXPECTED_REPEATS,
              Observed_rows = length(values),
              Finite_rows = sum(
                is.finite(
                  suppressWarnings(
                    as.numeric(values)
                  )
                )
              ),
              Complete_row_count = (
                length(values) == EXPECTED_REPEATS
              )
            )
          }
        )
      )
    }
  )
)

safe_fwrite(
  repeat_value_count_audit,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_repeat_value_count_audit.csv"
  )
)

if (any(!repeat_value_count_audit$Complete_row_count)) {
  stop(
    "The repeat-value count audit failed. Inspect ",
    "rpa_standardized_repeat_value_count_audit.csv."
  )
}

missing_metric_diagnostics <- dplyr::bind_rows(
  lapply(
    all_available_metrics,
    function(metric_name) {
      repeat_metrics |>
        dplyr::group_by(Model) |>
        dplyr::summarise(
          Total_repeats = dplyr::n(),
          Finite_repeats = sum(
            is.finite(.data[[metric_name]])
          ),
          Missing_or_nonfinite = sum(
            !is.finite(.data[[metric_name]])
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          Metric = metric_name,
          .before = 1
        )
    }
  )
)

safe_fwrite(
  missing_metric_diagnostics,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_missing_metric_diagnostics.csv"
  )
)

metric_eligibility <- dplyr::bind_rows(
  lapply(
    available_inferential_metrics,
    function(metric_name) {
      model_counts <- repeat_metrics |>
        dplyr::group_by(Model) |>
        dplyr::summarise(
          Finite_repeats = sum(
            is.finite(.data[[metric_name]])
          ),
          .groups = "drop"
        )

      complete_repeat_count <- repeat_metrics |>
        dplyr::select(
          Repeat,
          Model,
          dplyr::all_of(metric_name)
        ) |>
        tidyr::pivot_wider(
          names_from = Model,
          values_from =
            dplyr::all_of(metric_name)
        ) |>
        dplyr::filter(
          dplyr::if_all(
            dplyr::all_of(MODEL_ORDER),
            is.finite
          )
        ) |>
        nrow()

      tibble::tibble(
        Metric = metric_name,
        Minimum_finite_repeats =
          min(model_counts$Finite_repeats),
        Complete_paired_repeats =
          complete_repeat_count,
        Eligible = (
          min(model_counts$Finite_repeats) >=
            MIN_COMPLETE_REPEATS &&
            complete_repeat_count >=
            MIN_COMPLETE_REPEATS
        )
      )
    }
  )
)

safe_fwrite(
  metric_eligibility,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_metric_eligibility.csv"
  )
)

eligible_metrics <- metric_eligibility |>
  dplyr::filter(Eligible) |>
  dplyr::pull(Metric)

# ------------------------------------------------
# 7. Repeat-level descriptive statistics
# ------------------------------------------------

repeat_descriptive_statistics <- dplyr::bind_rows(
  lapply(
    all_available_metrics,
    function(metric_name) {
      dplyr::bind_rows(
        lapply(
          MODEL_ORDER,
          function(model_name) {
            values <- extract_metric_values(
              data = repeat_metrics,
              metric_name = metric_name,
              model_name = model_name,
              expected_repeats = EXPECTED_REPEATS
            )

            safe_descriptive_summary(
              values
            ) |>
              dplyr::mutate(
                Metric = metric_name,
                Model = model_name,
                Direction =
                  metric_direction(
                    metric_name
                  ),
                .before = 1
              )
          }
        )
      )
    }
  )
)

safe_fwrite(
  repeat_metrics,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_repeat_level_metrics.csv"
  )
)

safe_fwrite(
  repeat_descriptive_statistics,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_repeat_level_descriptive_statistics.csv"
  )
)

# ------------------------------------------------
# 8. Metric rankings
# ------------------------------------------------

ranking_rows <- list()

for (metric_name in eligible_metrics) {
  metric_means <- repeat_metrics |>
    dplyr::group_by(Model) |>
    dplyr::summarise(
      Mean = mean(
        .data[[metric_name]],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  if (metric_name %in% HIGHER_IS_BETTER) {
    metric_means <- metric_means |>
      dplyr::arrange(
        dplyr::desc(Mean)
      )
  } else if (metric_name %in% LOWER_IS_BETTER) {
    metric_means <- metric_means |>
      dplyr::arrange(Mean)
  } else {
    next
  }

  ranking_rows[[length(ranking_rows) + 1L]] <-
    metric_means |>
    dplyr::mutate(
      Metric = metric_name,
      Rank = dplyr::row_number(),
      Direction =
        metric_direction(metric_name),
      .before = 1
    )
}

model_rankings <- dplyr::bind_rows(
  ranking_rows
)

safe_fwrite(
  model_rankings,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_model_ranking_by_metric.csv"
  )
)

# ------------------------------------------------
# 9. Shapiro-Wilk diagnostics
# ------------------------------------------------

shapiro_by_model <- dplyr::bind_rows(
  lapply(
    eligible_metrics,
    function(metric_name) {
      dplyr::bind_rows(
        lapply(
          MODEL_ORDER,
          function(model_name) {
            values <- extract_metric_values(
              data = repeat_metrics,
              metric_name = metric_name,
              model_name = model_name,
              expected_repeats = EXPECTED_REPEATS
            )

            safe_shapiro(values) |>
              dplyr::mutate(
                Metric = metric_name,
                Model = model_name,
                .before = 1
              )
          }
        )
      )
    }
  )
)

safe_fwrite(
  shapiro_by_model,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_shapiro_wilk_by_model.csv"
  )
)

model_pairs <- combn(
  MODEL_ORDER,
  2L,
  simplify = FALSE
)

paired_difference_values <- list()
shapiro_paired_differences <- list()
paired_difference_completeness <- list()

for (metric_name in eligible_metrics) {
  wide_data <- repeat_metrics |>
    dplyr::select(
      Repeat,
      Model,
      dplyr::all_of(metric_name)
    ) |>
    tidyr::pivot_wider(
      names_from = Model,
      values_from =
        dplyr::all_of(metric_name)
    )

  for (model_pair in model_pairs) {
    model_1 <- model_pair[1L]
    model_2 <- model_pair[2L]

    valid <- (
      is.finite(wide_data[[model_1]]) &
        is.finite(wide_data[[model_2]])
    )

    difference <- (
      wide_data[[model_1]][valid] -
        wide_data[[model_2]][valid]
    )

    paired_difference_values[[
      length(paired_difference_values) + 1L
    ]] <- tibble::tibble(
      Metric = metric_name,
      Comparison = paste0(
        model_1,
        " - ",
        model_2
      ),
      Model_1 = model_1,
      Model_2 = model_2,
      Repeat = wide_data$Repeat[valid],
      Difference = difference
    )

    paired_difference_completeness[[
      length(paired_difference_completeness) + 1L
    ]] <- tibble::tibble(
      Metric = metric_name,
      Comparison = paste0(
        model_1,
        " - ",
        model_2
      ),
      Complete_pairs = sum(valid),
      Missing_pairs = sum(!valid)
    )

    shapiro_paired_differences[[
      length(shapiro_paired_differences) + 1L
    ]] <- safe_shapiro(
      difference
    ) |>
      dplyr::mutate(
        Metric = metric_name,
        Comparison = paste0(
          model_1,
          " - ",
          model_2
        ),
        Model_1 = model_1,
        Model_2 = model_2,
        .before = 1
      )
  }
}

paired_difference_values <-
  dplyr::bind_rows(
    paired_difference_values
  )

paired_difference_completeness <-
  dplyr::bind_rows(
    paired_difference_completeness
  )

shapiro_paired_differences <-
  dplyr::bind_rows(
    shapiro_paired_differences
  )

safe_fwrite(
  paired_difference_values,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_paired_difference_values.csv"
  )
)

safe_fwrite(
  paired_difference_completeness,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_paired_difference_completeness.csv"
  )
)

safe_fwrite(
  shapiro_paired_differences,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_shapiro_wilk_paired_differences.csv"
  )
)

test_family_recommendation <-
  shapiro_paired_differences |>
  dplyr::group_by(Metric) |>
  dplyr::summarise(
    Number_of_pairwise_tests = dplyr::n(),
    Number_rejecting_normality = sum(
      is.finite(p_value) &
        p_value <= ALPHA
    ),
    Number_indeterminate = sum(
      !is.finite(p_value)
    ),
    Recommended_family = if (
      all(
        is.finite(p_value) &
          p_value > ALPHA
      )
    ) {
      "parametric"
    } else {
      "nonparametric"
    },
    .groups = "drop"
  )

safe_fwrite(
  test_family_recommendation,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_test_family_recommendation.csv"
  )
)

# ------------------------------------------------
# 10. Repeat-level omnibus and pairwise tests
# ------------------------------------------------

omnibus_results <- list()
pairwise_results <- list()

if (RUN_REPEAT_LEVEL_TESTS) {
  for (metric_name in eligible_metrics) {
    wide_data <- repeat_metrics |>
      dplyr::select(
        Repeat,
        Model,
        dplyr::all_of(metric_name)
      ) |>
      tidyr::pivot_wider(
        names_from = Model,
        values_from =
          dplyr::all_of(metric_name)
      ) |>
      dplyr::filter(
        dplyr::if_all(
          dplyr::all_of(MODEL_ORDER),
          is.finite
        )
      )

    long_data <- wide_data |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(MODEL_ORDER),
        names_to = "Model",
        values_to = "Value"
      ) |>
      dplyr::mutate(
        Model = factor(
          Model,
          levels = MODEL_ORDER
        )
      )

    recommended_family <-
      test_family_recommendation |>
      dplyr::filter(
        Metric == metric_name
      ) |>
      dplyr::pull(
        Recommended_family
      )

    selected_family <- if (
      REPEAT_TEST_FAMILY == "auto"
    ) {
      recommended_family
    } else {
      REPEAT_TEST_FAMILY
    }

    if (selected_family == "parametric") {
      block_model <- stats::lm(
        Value ~ Model +
          factor(Repeat),
        data = long_data
      )

      anova_table <- stats::anova(
        block_model
      )

      model_row <- which(
        rownames(anova_table) == "Model"
      )

      omnibus_results[[
        length(omnibus_results) + 1L
      ]] <- tibble::tibble(
        Metric = metric_name,
        Test = "Randomized-block ANOVA",
        Test_family = "parametric",
        Statistic = anova_table$
          `F value`[model_row],
        Degrees_of_freedom =
          anova_table$Df[model_row],
        p_value = anova_table$
          `Pr(>F)`[model_row],
        Complete_repeats =
          nrow(wide_data)
      )
    } else {
      friedman_result <- stats::friedman.test(
        Value ~ Model | Repeat,
        data = long_data
      )

      omnibus_results[[
        length(omnibus_results) + 1L
      ]] <- tibble::tibble(
        Metric = metric_name,
        Test = "Friedman test",
        Test_family = "nonparametric",
        Statistic = unname(
          friedman_result$statistic
        ),
        Degrees_of_freedom = unname(
          friedman_result$parameter
        ),
        p_value = unname(
          friedman_result$p.value
        ),
        Complete_repeats =
          nrow(wide_data)
      )
    }

    metric_pairwise_rows <- list()

    for (model_pair in model_pairs) {
      model_1 <- model_pair[1L]
      model_2 <- model_pair[2L]

      x <- wide_data[[model_1]]
      y <- wide_data[[model_2]]
      difference <- x - y

      if (selected_family == "parametric") {
        test_result <- stats::t.test(
          x,
          y,
          paired = TRUE,
          conf.level = CONFIDENCE_LEVEL
        )

        confidence_interval <-
          unname(test_result$conf.int)

        metric_pairwise_rows[[
          length(metric_pairwise_rows) + 1L
        ]] <- tibble::tibble(
          Metric = metric_name,
          Comparison = paste0(
            model_1,
            " - ",
            model_2
          ),
          Model_1 = model_1,
          Model_2 = model_2,
          Test = "Paired t-test",
          Test_family = "parametric",
          n_pairs = length(difference),
          Mean_model_1 = mean(x),
          Mean_model_2 = mean(y),
          Mean_difference =
            mean(difference),
          Median_difference =
            stats::median(difference),
          CI_lower =
            confidence_interval[1L],
          CI_upper =
            confidence_interval[2L],
          Statistic = unname(
            test_result$statistic
          ),
          p_value_raw = unname(
            test_result$p.value
          ),
          Effect_size = cohens_dz(
            difference
          ),
          Effect_size_name = "Cohen dz",
          Direction =
            metric_direction(
              metric_name
            ),
          Favored_model =
            favored_model(
              metric_name,
              model_1,
              model_2,
              mean(difference)
            )
        )
      } else {
        test_result <- suppressWarnings(
          stats::wilcox.test(
            x,
            y,
            paired = TRUE,
            exact = FALSE,
            correct = TRUE,
            conf.int = TRUE,
            conf.level =
              CONFIDENCE_LEVEL
          )
        )

        confidence_interval <- if (
          is.null(test_result$conf.int)
        ) {
          c(NA_real_, NA_real_)
        } else {
          unname(test_result$conf.int)
        }

        metric_pairwise_rows[[
          length(metric_pairwise_rows) + 1L
        ]] <- tibble::tibble(
          Metric = metric_name,
          Comparison = paste0(
            model_1,
            " - ",
            model_2
          ),
          Model_1 = model_1,
          Model_2 = model_2,
          Test =
            "Paired Wilcoxon signed-rank test",
          Test_family =
            "nonparametric",
          n_pairs = length(difference),
          Mean_model_1 = mean(x),
          Mean_model_2 = mean(y),
          Mean_difference =
            mean(difference),
          Median_difference =
            stats::median(difference),
          CI_lower =
            confidence_interval[1L],
          CI_upper =
            confidence_interval[2L],
          Statistic = unname(
            test_result$statistic
          ),
          p_value_raw = unname(
            test_result$p.value
          ),
          Effect_size =
            rank_biserial_paired(
              difference
            ),
          Effect_size_name =
            "Matched-pairs rank-biserial correlation",
          Direction =
            metric_direction(
              metric_name
            ),
          Favored_model =
            favored_model(
              metric_name,
              model_1,
              model_2,
              stats::median(
                difference
              )
            )
        )
      }
    }

    metric_pairwise_data <- dplyr::bind_rows(
      metric_pairwise_rows
    ) |>
      dplyr::mutate(
        p_value_adjusted =
          stats::p.adjust(
            p_value_raw,
            method =
              P_ADJUST_METHOD
          ),
        Significant_adjusted =
          p_value_adjusted <
          ALPHA
      )

    pairwise_results[[
      length(pairwise_results) + 1L
    ]] <- metric_pairwise_data
  }
}

omnibus_results <- dplyr::bind_rows(
  omnibus_results
) |>
  dplyr::mutate(
    Significant =
      p_value < ALPHA,
    Interpretation = paste0(
      "Stability analysis across repeated CV; ",
      "repeats reuse the same patients."
    )
  )

pairwise_results <- dplyr::bind_rows(
  pairwise_results
) |>
  dplyr::mutate(
    Interpretation = paste0(
      "Stability analysis across repeated CV; ",
      "repeats reuse the same patients."
    )
  )

safe_fwrite(
  omnibus_results,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_omnibus_model_comparisons.csv"
  )
)

safe_fwrite(
  pairwise_results,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_pairwise_model_comparisons.csv"
  )
)

# ------------------------------------------------
# 11. Paired patient bootstrap
# ------------------------------------------------

bootstrap_model_estimates <- tibble::tibble()
bootstrap_pairwise_results <- tibble::tibble()
bootstrap_draws_long <- tibble::tibble()

if (RUN_PATIENT_BOOTSTRAP) {
  required_prediction_columns <- c(
    "Row_ID",
    "Model",
    "Y_true",
    "p_mean",
    "p_sd",
    "p_lwr",
    "p_upr"
  )

  missing_prediction_columns <- setdiff(
    required_prediction_columns,
    names(aggregated_predictions)
  )

  if (length(missing_prediction_columns) > 0L) {
    stop(
      "Aggregated prediction file is missing: ",
      paste(
        missing_prediction_columns,
        collapse = ", "
      )
    )
  }

  prediction_counts <- aggregated_predictions |>
    dplyr::count(
      Row_ID,
      name = "Model_count"
    )

  if (
    any(
      prediction_counts$Model_count !=
        length(MODEL_ORDER)
    )
  ) {
    stop(
      "Each patient must have one aggregated prediction ",
      "for every model."
    )
  }

  prediction_by_model <- lapply(
    MODEL_ORDER,
    function(model_name) {
      aggregated_predictions |>
        dplyr::filter(
          Model == model_name
        ) |>
        dplyr::arrange(Row_ID)
    }
  )

  names(prediction_by_model) <- MODEL_ORDER

  reference_row_ids <-
    prediction_by_model[[1L]]$Row_ID

  reference_y <-
    prediction_by_model[[1L]]$Y_true

  for (model_name in MODEL_ORDER) {
    model_data <- prediction_by_model[[model_name]]

    if (!identical(
      model_data$Row_ID,
      reference_row_ids
    )) {
      stop(
        "Patient ordering differs across models."
      )
    }

    if (!identical(
      as.integer(model_data$Y_true),
      as.integer(reference_y)
    )) {
      stop(
        "Outcome values differ across model prediction files."
      )
    }
  }

  bootstrap_metric_names <- unique(
    c(
      INFERENTIAL_METRICS,
      "mean_uncertainty",
      "mean_interval_width",
      "selective_coverage"
    )
  )

  point_estimate_rows <- list()

  for (model_name in MODEL_ORDER) {
    model_data <- prediction_by_model[[model_name]]

    values <- calculate_metric_vector(
      y = model_data$Y_true,
      p_mean = model_data$p_mean,
      p_sd = model_data$p_sd,
      p_lwr = model_data$p_lwr,
      p_upr = model_data$p_upr
    )

    point_estimate_rows[[
      length(point_estimate_rows) + 1L
    ]] <- tibble::tibble(
      Model = model_name,
      Metric = names(values),
      Point_estimate = as.numeric(values)
    )
  }

  point_estimates <- dplyr::bind_rows(
    point_estimate_rows
  ) |>
    dplyr::filter(
      Metric %in%
        bootstrap_metric_names
    )

  set.seed(BOOTSTRAP_SEED)

  bootstrap_rows <- vector(
    mode = "list",
    length = BOOTSTRAP_REPS *
      length(MODEL_ORDER)
  )

  bootstrap_index <- 0L
  patient_count <- length(reference_row_ids)

  for (bootstrap_id in seq_len(BOOTSTRAP_REPS)) {
    sampled_index <- sample.int(
      patient_count,
      size = patient_count,
      replace = TRUE
    )

    for (model_name in MODEL_ORDER) {
      bootstrap_index <- bootstrap_index + 1L
      model_data <- prediction_by_model[[model_name]]

      values <- calculate_metric_vector(
        y = model_data$Y_true[
          sampled_index
        ],
        p_mean = model_data$p_mean[
          sampled_index
        ],
        p_sd = model_data$p_sd[
          sampled_index
        ],
        p_lwr = model_data$p_lwr[
          sampled_index
        ],
        p_upr = model_data$p_upr[
          sampled_index
        ]
      )

      bootstrap_rows[[bootstrap_index]] <-
        tibble::as_tibble_row(
          c(
            list(
              Bootstrap = bootstrap_id,
              Model = model_name
            ),
            as.list(values)
          )
        )
    }

    if (
      bootstrap_id %% 100L == 0L
    ) {
      cat(
        sprintf(
          "\rPatient bootstrap %d/%d",
          bootstrap_id,
          BOOTSTRAP_REPS
        )
      )
    }
  }

  cat("\n")

  bootstrap_draws <- dplyr::bind_rows(
    bootstrap_rows
  )

  bootstrap_draws_long <- bootstrap_draws |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(
        bootstrap_metric_names
      ),
      names_to = "Metric",
      values_to = "Value"
    )

  bootstrap_model_estimates <-
    bootstrap_draws_long |>
    dplyr::filter(
      is.finite(Value)
    ) |>
    dplyr::group_by(
      Model,
      Metric
    ) |>
    dplyr::summarise(
      Bootstrap_finite_replicates =
        dplyr::n(),
      Bootstrap_mean =
        mean(Value),
      Bootstrap_standard_error =
        stats::sd(Value),
      CI_lower =
        as.numeric(
          stats::quantile(
            Value,
            probs =
              (1 - CONFIDENCE_LEVEL) / 2,
            names = FALSE
          )
        ),
      CI_upper =
        as.numeric(
          stats::quantile(
            Value,
            probs =
              1 -
              (1 - CONFIDENCE_LEVEL) / 2,
            names = FALSE
          )
        ),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      point_estimates,
      by = c(
        "Model",
        "Metric"
      )
    ) |>
    dplyr::mutate(
      Direction =
        vapply(
          Metric,
          metric_direction,
          character(1)
        )
    )

  bootstrap_pairwise_rows <- list()

  for (metric_name in bootstrap_metric_names) {
    metric_draws <- bootstrap_draws_long |>
      dplyr::filter(
        Metric == metric_name
      ) |>
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
      dplyr::filter(
        Metric == metric_name
      )

    metric_pair_rows <- list()

    for (model_pair in model_pairs) {
      model_1 <- model_pair[1L]
      model_2 <- model_pair[2L]

      difference <- (
        metric_draws[[model_1]] -
          metric_draws[[model_2]]
      )

      difference <- difference[
        is.finite(difference)
      ]

      point_1 <- metric_points |>
        dplyr::filter(
          Model == model_1
        ) |>
        dplyr::pull(
          Point_estimate
        )

      point_2 <- metric_points |>
        dplyr::filter(
          Model == model_2
        ) |>
        dplyr::pull(
          Point_estimate
        )

      point_difference <-
        point_1 - point_2

      if (length(difference) == 0L) {
        ci_lower <- NA_real_
        ci_upper <- NA_real_
        p_value_raw <- NA_real_
      } else {
        ci_lower <- as.numeric(
          stats::quantile(
            difference,
            probs =
              (1 - CONFIDENCE_LEVEL) / 2,
            names = FALSE
          )
        )

        ci_upper <- as.numeric(
          stats::quantile(
            difference,
            probs =
              1 -
              (1 - CONFIDENCE_LEVEL) / 2,
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

      metric_pair_rows[[
        length(metric_pair_rows) + 1L
      ]] <- tibble::tibble(
        Metric = metric_name,
        Comparison = paste0(
          model_1,
          " - ",
          model_2
        ),
        Model_1 = model_1,
        Model_2 = model_2,
        Model_1_point_estimate =
          point_1,
        Model_2_point_estimate =
          point_2,
        Point_difference =
          point_difference,
        Bootstrap_finite_replicates =
          length(difference),
        CI_lower = ci_lower,
        CI_upper = ci_upper,
        p_value_raw = p_value_raw,
        Direction =
          metric_direction(
            metric_name
          ),
        Favored_model =
          favored_model(
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
        p_value_adjusted =
          stats::p.adjust(
            p_value_raw,
            method =
              P_ADJUST_METHOD
          ),
        CI_excludes_zero = (
          is.finite(CI_lower) &
            is.finite(CI_upper) &
            (
              CI_lower > 0 |
                CI_upper < 0
            )
        ),
        Significant_adjusted = (
          is.finite(
            p_value_adjusted
          ) &
            p_value_adjusted <
            ALPHA
        ),
        Interpretation = paste0(
          "Paired patient bootstrap using fixed aggregated ",
          "out-of-fold predictions."
        )
      )

    bootstrap_pairwise_rows[[
      length(bootstrap_pairwise_rows) + 1L
    ]] <- metric_pair_data
  }

  bootstrap_pairwise_results <-
    dplyr::bind_rows(
      bootstrap_pairwise_rows
    )

  safe_fwrite(
    bootstrap_model_estimates,
    file.path(
      OUTPUT_DIR,
      "rpa_standardized_patient_bootstrap_model_estimates.csv"
    )
  )

  safe_fwrite(
    bootstrap_pairwise_results,
    file.path(
      OUTPUT_DIR,
      "rpa_standardized_patient_bootstrap_pairwise_comparisons.csv"
    )
  )

  # Bootstrap draws are useful for auditing but can be omitted from a
  # manuscript supplement if file size is a concern.
  safe_fwrite(
    bootstrap_draws_long,
    file.path(
      OUTPUT_DIR,
      "rpa_standardized_patient_bootstrap_draws_long.csv"
    )
  )
}

# ------------------------------------------------
# 12. Figures
# ------------------------------------------------

repeat_plot_data <- repeat_metrics |>
  dplyr::select(
    Repeat,
    Model,
    dplyr::all_of(
      all_available_metrics
    )
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      all_available_metrics
    ),
    names_to = "Metric",
    values_to = "Value"
  )

repeat_metric_plot <- ggplot(
  repeat_plot_data,
  aes(
    x = Model,
    y = Value,
    fill = Model
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.35
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(values = model_colors) +
  labs(
    title = "RPA Standardized: Repeat-Level Metric Distributions",
    subtitle = paste0(
      "Each point distribution is based on pooled out-of-fold ",
      "predictions within a CV repeat"
    ),
    x = "Model",
    y = "Metric value",
    fill = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title =
      element_text(face = "bold"),
    strip.text =
      element_text(face = "bold"),
    panel.grid.minor =
      element_blank()
  )

save_plot_pdf_png(
  repeat_metric_plot,
  file.path(
    FIGURE_DIR,
    "rpa_standardized_repeat_level_metric_distributions"
  ),
  width = 14,
  height = 12
)

print(repeat_metric_plot)

if (
  RUN_PATIENT_BOOTSTRAP &&
  nrow(bootstrap_pairwise_results) > 0L
) {
  bootstrap_forest_plot <- ggplot(
    bootstrap_pairwise_results |>
      dplyr::filter(
        Metric %in%
          INFERENTIAL_METRICS
      ),
    aes(
      x = Point_difference,
      y = Comparison
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray45"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_lower,
        xmax = CI_upper
      ),
      height = 0.2
    ) +
    geom_point(size = 2.4) +
    facet_wrap(
      ~ Metric,
      scales = "free_x",
      ncol = 3
    ) +
    labs(
      title = "RPA Standardized: Paired Patient-Bootstrap Differences",
      subtitle = paste0(
        "Difference = first model minus second model; ",
        round(
          CONFIDENCE_LEVEL * 100
        ),
        "% percentile intervals"
      ),
      x = "Metric difference",
      y = "Model comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title =
        element_text(face = "bold"),
      strip.text =
        element_text(face = "bold"),
      panel.grid.minor =
        element_blank()
    )

  save_plot_pdf_png(
    bootstrap_forest_plot,
    file.path(
      FIGURE_DIR,
      "rpa_standardized_patient_bootstrap_pairwise_forest"
    ),
    width = 13,
    height = 11
  )

  print(bootstrap_forest_plot)
}

# ------------------------------------------------
# 13. Analysis manifest
# ------------------------------------------------

analysis_manifest <- tibble::tibble(
  Item = c(
    "Repeat-level input",
    "Aggregated-prediction input",
    "Model-status input",
    "Expected repeats",
    "Detected repeats",
    "Expected folds",
    "Detected folds",
    "Models",
    "Repeat-level metric unit",
    "Repeat-value extraction",
    "Repeat test family mode",
    "Repeat-level tests executed",
    "Patient bootstrap executed",
    "Patient bootstrap replicates",
    "Bootstrap seed",
    "Alpha",
    "Confidence level",
    "P-value adjustment",
    "Eligible inferential metrics",
    "Descriptive-only metrics",
    "Primary inference note",
    "Repeat-level inference note"
  ),
  Value = c(
    REPEAT_METRIC_FILE,
    AGGREGATED_PREDICTION_FILE,
    MODEL_STATUS_FILE,
    EXPECTED_REPEATS,
    detected_repeats,
    EXPECTED_FOLDS,
    detected_folds,
    paste(MODEL_ORDER, collapse = ", "),
    paste0(
      "Metrics calculated from all pooled out-of-fold ",
      "predictions within each repeat"
    ),
    paste0(
      "Complete vectors extracted with data[[metric_name]][model_index]; ",
      "exactly EXPECTED_REPEATS values required per model and metric"
    ),
    REPEAT_TEST_FAMILY,
    RUN_REPEAT_LEVEL_TESTS,
    RUN_PATIENT_BOOTSTRAP,
    BOOTSTRAP_REPS,
    BOOTSTRAP_SEED,
    ALPHA,
    CONFIDENCE_LEVEL,
    P_ADJUST_METHOD,
    paste(
      eligible_metrics,
      collapse = ", "
    ),
    paste(
      available_descriptive_metrics,
      collapse = ", "
    ),
    paste0(
      "Paired patient bootstrap conditions on final aggregated ",
      "out-of-fold predictions"
    ),
    paste0(
      "Repeated-CV tests are stability analyses because the same ",
      "patients are reused across repeats"
    )
  )
)

safe_fwrite(
  analysis_manifest,
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_statistical_analysis_manifest.csv"
  )
)

saveRDS(
  list(
    status_audit =
      status_audit,
    fold_completeness =
      fold_completeness,
    repeat_completeness =
      repeat_completeness,
    repeat_value_count_audit =
      repeat_value_count_audit,
    missing_metric_diagnostics =
      missing_metric_diagnostics,
    metric_eligibility =
      metric_eligibility,
    repeat_descriptive_statistics =
      repeat_descriptive_statistics,
    model_rankings =
      model_rankings,
    shapiro_by_model =
      shapiro_by_model,
    shapiro_paired_differences =
      shapiro_paired_differences,
    test_family_recommendation =
      test_family_recommendation,
    omnibus_results =
      omnibus_results,
    pairwise_results =
      pairwise_results,
    bootstrap_model_estimates =
      bootstrap_model_estimates,
    bootstrap_pairwise_results =
      bootstrap_pairwise_results,
    analysis_manifest =
      analysis_manifest
  ),
  file.path(
    OUTPUT_DIR,
    "rpa_standardized_updated_statistical_analysis_v2.rds"
  ),
  compress = "xz"
)

# ------------------------------------------------
# 14. Final summary
# ------------------------------------------------

cat("\n================================================\n")
cat("B-RISA STANDARDIZED RPA STATISTICAL ANALYSIS COMPLETE\n")
cat("================================================\n")

cat("\nModel-status audit:\n")
print(status_audit)

cat("\nEligible inferential metrics:\n")
print(eligible_metrics)

cat("\nRepeat-level omnibus results:\n")
print(omnibus_results)

if (RUN_PATIENT_BOOTSTRAP) {
  cat("\nPatient-bootstrap model estimates:\n")
  print(bootstrap_model_estimates)

  cat("\nPatient-bootstrap pairwise comparisons:\n")
  print(bootstrap_pairwise_results)
}

cat("\nTables saved in:\n")
cat(OUTPUT_DIR, "\n")

cat("\nFigures saved in:\n")
cat(FIGURE_DIR, "\n")

cat("\nPrimary files for interpretation:\n")
cat("  rpa_standardized_model_status_audit.csv\n")
cat("  rpa_standardized_repeat_value_count_audit.csv\n")
cat("  rpa_standardized_repeat_level_descriptive_statistics.csv\n")
cat("  rpa_standardized_patient_bootstrap_model_estimates.csv\n")
cat("  rpa_standardized_patient_bootstrap_pairwise_comparisons.csv\n")
cat("  rpa_standardized_omnibus_model_comparisons.csv\n")
cat("  rpa_standardized_pairwise_model_comparisons.csv\n")
cat("================================================\n")

