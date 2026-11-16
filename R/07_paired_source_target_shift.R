# ================================================================
# B-RISA Paired Source-to-Target Shift Experiment
#
# Generates 500 paired RPA-like source/target replicates and evaluates frozen
# BLR, BLDA, and BART fits under no shift, prevalence shift, neutrophil-
# missingness shift, and combined shift. The script also computes the Monte
# Carlo confidence intervals reported in the manuscript.
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

PARAM_DIR <- file.path(
  PROJECT_ROOT,
  "data",
  "derived"
)

RPA_MU_FILE <- file.path(
  PARAM_DIR,
  "rpa_standardized_corrected_mean_vectors.csv"
)
RPA_SIGMA0_FILE <- file.path(
  PARAM_DIR,
  "rpa_standardized_corrected_Sigma_Y0_no_appendicitis.csv"
)
RPA_SIGMA1_FILE <- file.path(
  PARAM_DIR,
  "rpa_standardized_corrected_Sigma_Y1_appendicitis.csv"
)

# The Bayesian model implementations are embedded below.
# Bayesian model implementations are embedded below; no external model script is sourced.

needed_parameter_files <- c(
  RPA_MU_FILE,
  RPA_SIGMA0_FILE,
  RPA_SIGMA1_FILE
)

if (any(!file.exists(needed_parameter_files))) {
  stop(
    "Missing RPA empirical parameter file(s):\n",
    paste(
      needed_parameter_files[!file.exists(needed_parameter_files)],
      collapse = "\n"
    )
  )
}

OUTPUT_ROOT <- file.path(
  PROJECT_ROOT,
  "analysis_outputs",
  "paired_source_target_shift"
)

TABLE_DIR <- file.path(OUTPUT_ROOT, "tables")
FIGURE_DIR <- file.path(OUTPUT_ROOT, "figures")
CHECKPOINT_DIR <- file.path(OUTPUT_ROOT, "checkpoints")
OBJECT_DIR <- file.path(OUTPUT_ROOT, "r_objects")

for (d in c(
  OUTPUT_ROOT,
  TABLE_DIR,
  FIGURE_DIR,
  CHECKPOINT_DIR,
  OBJECT_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------
# 1. Targeted simulation settings
# ------------------------------------------------

BASE_SEED <- 20261101L

N_REPS <- 500L
# For a pilot, temporarily set:
# N_REPS <- 20L

CHECKPOINT_EVERY <- 10L
RESUME_ANALYSIS <- TRUE

N_SOURCE <- 780L
N_TARGET <- 301L
PREDICTOR_COUNT <- 3L

SOURCE_PREVALENCE <- 463 / 780
TARGET_PREVALENCE <- 230 / 301

SOURCE_NEUTROPHIL_MISSING <- 0.12949
TARGET_NEUTROPHIL_MISSING <- 0.73090

MODEL_SET <- c("BLR", "BLDA", "BART")

CLASSIFICATION_THRESHOLD <- 0.50
PRIMARY_COVERAGE <- 0.80

# Bayesian settings used throughout the manuscript.
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

BLDA_MEAN_PRIOR_KAPPA <- 0.01
BLDA_IW_DF_OFFSET <- 2L
BLDA_PRIOR_COVARIANCE_SCALE <- 1.0
BLDA_CLASS_PRIOR_ALPHA <- 0.5
BLDA_CLASS_PRIOR_BETA <- 0.5
BLDA_PRECISION_EIGEN_FLOOR <- 1e-10

EPS <- 1e-6
DEGENERATE_PREDICTION_TOL <- 1e-12

# ------------------------------------------------
# 2. Packages
# ------------------------------------------------

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble",
  "pROC",
  "PRROC",
  "MASS",
  "Matrix",
  "BART"
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
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(pROC)
  library(PRROC)
  library(MASS)
  library(Matrix)
})

if (!requireNamespace("BART", quietly = TRUE)) {
  stop("Package 'BART' is required.")
}

# ------------------------------------------------
# 3. Self-contained preprocessing and Bayesian model functions
# ------------------------------------------------
#
# These functions reproduce the corresponding V7 implementations used
# in the current B-RISA simulation analysis. Embedding them here avoids
# dependence on another model-script filename or location.

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
      "At least one predictor is entirely missing in the source training data."
    )
  }

  imputed <- train_data[
    ,
    feature_columns,
    drop = FALSE
  ]

  for (feature_name in feature_columns) {
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
    missing_index <- is.na(x[[feature_name]])

    if (any(missing_index)) {
      x[[feature_name]][missing_index] <-
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

softplus_stable <- function(x) {
  pmax(x, 0) +
    log1p(exp(-abs(x)))
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

  repaired <-
    decomposition$vectors %*%
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
  n_draws = N_POST_DRAWS_BLR
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  if (length(unique(y_train)) < 2L) {
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

  prior_precision <- 1 / prior_sd^2

  negative_log_posterior <- function(beta) {
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
    gr = gradient_negative_log_posterior,
    method = "BFGS",
    control = list(
      maxit = BLR_OPTIM_MAXIT,
      reltol = BLR_OPTIM_RELTOL
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
    stop(
      "BLR Hessian has non-finite eigenvalues."
    )
  }

  posterior_covariance <- tryCatch(
    solve(hessian),
    error = function(e) NULL
  )

  if (
    is.null(posterior_covariance) ||
    any(!is.finite(posterior_covariance))
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
    Sigma = posterior_covariance
  )

  beta_draws <- matrix(
    beta_draws,
    nrow = n_draws,
    ncol = parameter_count
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

  probability_draws <-
    validate_prediction_draws(
      draws = probability_draws,
      expected_draws = n_draws,
      expected_test_size = nrow(x_test),
      model_name = "BLR"
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
    draws = probability_draws,
    diagnostics = diagnostics
  )
}

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

  nu_prior <-
    predictor_count +
    BLDA_IW_DF_OFFSET

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

  posterior_scale <-
    nearest_positive_definite(
      posterior_scale,
      eigen_floor =
        BLDA_PRECISION_EIGEN_FLOOR
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

  posterior_precision_scale <-
    nearest_positive_definite(
      posterior_precision_scale,
      eigen_floor =
        BLDA_PRECISION_EIGEN_FLOOR
    )

  precision_draw_array <- stats::rWishart(
    n = n_draws,
    df = nu_post,
    Sigma = posterior_precision_scale
  )

  class_1_prior_draws <- stats::rbeta(
    n_draws,
    shape1 =
      BLDA_CLASS_PRIOR_ALPHA + n_1,
    shape2 =
      BLDA_CLASS_PRIOR_BETA + n_0
  )

  probability_draws <- matrix(
    NA_real_,
    nrow = n_draws,
    ncol = test_count
  )

  sampled_precision_min_eigen <-
    numeric(n_draws)
  sampled_precision_condition <-
    numeric(n_draws)

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

    precision_chol <- chol(
      precision_draw
    )

    mean_0_draw <-
      posterior_mean_0 +
      as.vector(
        backsolve(
          precision_chol,
          stats::rnorm(
            predictor_count
          )
        )
      ) /
      sqrt(kappa_post_0)

    mean_1_draw <-
      posterior_mean_1 +
      as.vector(
        backsolve(
          precision_chol,
          stats::rnorm(
            predictor_count
          )
        )
      ) /
      sqrt(kappa_post_1)

    class_1_prior <-
      class_1_prior_draws[draw_id]
    class_0_prior <-
      1 - class_1_prior

    precision_mean_0 <- as.vector(
      precision_draw %*%
        mean_0_draw
    )

    precision_mean_1 <- as.vector(
      precision_draw %*%
        mean_1_draw
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

    probability_draws[
      draw_id,
    ] <- plogis(
      discriminant_1 -
        discriminant_0
    )
  }

  probability_draws <-
    clip_prob(
      probability_draws
    )

  probability_draws <-
    validate_prediction_draws(
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
      mean(
        sampled_precision_min_eigen
      ),
    BLDA_max_sampled_precision_condition =
      max(
        sampled_precision_condition
      ),
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

predict_bart <- function(
  x_train,
  y_train,
  x_test,
  n_draws = N_POST_DRAWS_BART
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  if (length(unique(y_train)) < 2L) {
    stop(
      "BART training labels contain only one class."
    )
  }

  fit <- BART::pbart(
    x.train = x_train,
    y.train = y_train,
    x.test = x_test,
    ntree = BART_NTREE,
    ndpost = n_draws,
    nskip = N_BURN_BART,
    nkeeptest = n_draws,
    printevery = 100000L
  )

  if (is.null(fit$yhat.test)) {
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

  reversed_dimension <-
    rev(
      expected_dimension
    )

  orientation_action <-
    "accepted_documented_draw_by_patient_orientation"

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
      draws = probability_draws,
      expected_draws = n_draws,
      expected_test_size = nrow(x_test),
      model_name = "BART"
    )

  list(
    draws = probability_draws,
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
    BLR = N_POST_DRAWS_BLR,
    BLDA = N_POST_DRAWS_BLDA,
    BART = N_POST_DRAWS_BART,
    stop(
      "Unknown model: ",
      model_name
    )
  )
}

# ------------------------------------------------
# 4. Read empirical RPA class parameters
# ------------------------------------------------

read_mean_vectors <- function(path) {
  df <- data.table::fread(
    path,
    data.table = FALSE
  )

  required <- c(
    "Feature",
    "mu_Y0_no_appendicitis",
    "mu_Y1_appendicitis"
  )

  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(
      "Mean-vector file missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  list(
    features = as.character(df$Feature),
    mu0 = as.numeric(df$mu_Y0_no_appendicitis),
    mu1 = as.numeric(df$mu_Y1_appendicitis)
  )
}

read_matrix_csv <- function(path) {
  df <- data.table::fread(
    path,
    data.table = FALSE
  )

  rn <- as.character(df[[1L]])
  mat <- as.matrix(df[, -1L, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- rn
  colnames(mat) <- names(df)[-1L]
  mat
}

nearest_pd <- function(Sigma, ridge = 1e-8) {
  Sigma <- (as.matrix(Sigma) + t(as.matrix(Sigma))) / 2
  eig <- eigen(
    Sigma,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (all(is.finite(eig)) && min(eig) > 1e-8) {
    return(Sigma)
  }

  out <- as.matrix(
    Matrix::nearPD(
      Sigma,
      corr = FALSE
    )$mat
  )

  out + diag(ridge, nrow(out))
}

mu <- read_mean_vectors(RPA_MU_FILE)
Sigma0 <- nearest_pd(
  read_matrix_csv(RPA_SIGMA0_FILE)
)
Sigma1 <- nearest_pd(
  read_matrix_csv(RPA_SIGMA1_FILE)
)

if (length(mu$mu0) != PREDICTOR_COUNT ||
    length(mu$mu1) != PREDICTOR_COUNT) {
  stop("Expected exactly three RPA predictors.")
}

FEATURES <- paste0("X", seq_len(PREDICTOR_COUNT))

# ------------------------------------------------
# 5. Data generation
# ------------------------------------------------

make_exact_binary_labels <- function(
  n,
  prevalence,
  seed
) {
  set.seed(seed)
  n1 <- round(n * prevalence)
  n0 <- n - n1
  sample(
    c(rep(0L, n0), rep(1L, n1)),
    size = n,
    replace = FALSE
  )
}

generate_class_conditional_sample <- function(
  n,
  prevalence,
  seed
) {
  y <- make_exact_binary_labels(
    n,
    prevalence,
    seed
  )

  X <- matrix(
    NA_real_,
    nrow = n,
    ncol = PREDICTOR_COUNT
  )

  idx0 <- which(y == 0L)
  idx1 <- which(y == 1L)

  if (length(idx0) > 0L) {
    X[idx0, ] <- MASS::mvrnorm(
      n = length(idx0),
      mu = mu$mu0,
      Sigma = Sigma0
    )
  }

  if (length(idx1) > 0L) {
    X[idx1, ] <- MASS::mvrnorm(
      n = length(idx1),
      mu = mu$mu1,
      Sigma = Sigma1
    )
  }

  colnames(X) <- FEATURES

  list(
    X = X,
    y = y
  )
}

apply_neutrophil_mask <- function(
  X,
  mask
) {
  out <- X
  out[mask, 3L] <- NA_real_
  out
}

make_source_sample <- function(rep_id) {
  seed <- BASE_SEED + rep_id * 10000L + 101L

  generated <- generate_class_conditional_sample(
    N_SOURCE,
    SOURCE_PREVALENCE,
    seed
  )

  set.seed(
    BASE_SEED +
      rep_id * 10000L +
      102L
  )

  source_missing <- runif(N_SOURCE) <
    SOURCE_NEUTROPHIL_MISSING

  X <- apply_neutrophil_mask(
    generated$X,
    source_missing
  )

  df <- as.data.frame(X)
  df$Y_true <- as.integer(generated$y)
  df$Y_observed <- as.integer(generated$y)

  list(
    data = df,
    missing_rate = mean(source_missing)
  )
}

# Generate common class-specific pools large enough for both target prevalences.
make_target_conditions <- function(rep_id) {
  n1_source_prev <- round(N_TARGET * SOURCE_PREVALENCE)
  n0_source_prev <- N_TARGET - n1_source_prev

  n1_target_prev <- round(N_TARGET * TARGET_PREVALENCE)
  n0_target_prev <- N_TARGET - n1_target_prev

  n1_pool <- max(n1_source_prev, n1_target_prev)
  n0_pool <- max(n0_source_prev, n0_target_prev)

  seed0 <- BASE_SEED + rep_id * 10000L + 201L
  seed1 <- BASE_SEED + rep_id * 10000L + 202L

  set.seed(seed0)
  X0 <- MASS::mvrnorm(
    n = n0_pool,
    mu = mu$mu0,
    Sigma = Sigma0
  )

  set.seed(seed1)
  X1 <- MASS::mvrnorm(
    n = n1_pool,
    mu = mu$mu1,
    Sigma = Sigma1
  )

  make_target <- function(n0, n1, prevalence_label) {
    X <- rbind(
      X0[seq_len(n0), , drop = FALSE],
      X1[seq_len(n1), , drop = FALSE]
    )

    y <- c(
      rep(0L, n0),
      rep(1L, n1)
    )

    # Deterministic common shuffle for this prevalence composition.
    set.seed(
      BASE_SEED +
        rep_id * 10000L +
        if (prevalence_label == "source") 203L else 204L
    )

    ord <- sample.int(
      length(y),
      size = length(y),
      replace = FALSE
    )

    X <- X[ord, , drop = FALSE]
    y <- y[ord]

    colnames(X) <- FEATURES

    # One common uniform creates nested low/high missingness masks.
    set.seed(
      BASE_SEED +
        rep_id * 10000L +
        if (prevalence_label == "source") 205L else 206L
    )

    u <- runif(length(y))

    low_mask <- u < SOURCE_NEUTROPHIL_MISSING
    high_mask <- u < TARGET_NEUTROPHIL_MISSING

    list(
      low = list(
        X = apply_neutrophil_mask(X, low_mask),
        y = y,
        missing_rate = mean(low_mask)
      ),
      high = list(
        X = apply_neutrophil_mask(X, high_mask),
        y = y,
        missing_rate = mean(high_mask)
      )
    )
  }

  source_prev <- make_target(
    n0_source_prev,
    n1_source_prev,
    "source"
  )

  target_prev <- make_target(
    n0_target_prev,
    n1_target_prev,
    "target"
  )

  list(
    no_shift = source_prev$low,
    prevalence_shift = target_prev$low,
    missingness_shift = source_prev$high,
    combined_shift = target_prev$high
  )
}

# ------------------------------------------------
# 6. Metrics
# ------------------------------------------------

safe_auc <- function(y, p) {
  if (length(unique(y)) < 2L) return(NA_real_)
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
  )
}

calibration_parameters <- function(y, p) {
  p <- pmin(pmax(p, EPS), 1 - EPS)
  lp <- qlogis(p)

  if (length(unique(y)) < 2L ||
      stats::sd(lp) <= 1e-12) {
    return(c(
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  fit <- tryCatch(
    stats::glm(
      y ~ lp,
      family = stats::binomial()
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(c(
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  cf <- stats::coef(fit)

  c(
    intercept = unname(cf[[1L]]),
    slope = unname(cf[[2L]])
  )
}

selective_risk <- function(
  y,
  p_mean,
  p_sd,
  rule,
  coverage = PRIMARY_COVERAGE
) {
  n <- length(y)
  keep_n <- max(
    1L,
    min(n, ceiling(coverage * n))
  )

  pred <- as.integer(
    p_mean >= CLASSIFICATION_THRESHOLD
  )

  error <- as.integer(pred != y)

  if (rule == "posterior_sd") {
    ord <- order(p_sd, seq_len(n))
  } else if (rule == "confidence_margin") {
    margin <- abs(
      p_mean -
        CLASSIFICATION_THRESHOLD
    )
    ord <- order(-margin, seq_len(n))
  } else {
    stop("Unknown selective-prediction rule.")
  }

  mean(error[ord[seq_len(keep_n)]])
}

calc_target_metrics <- function(
  y,
  p_mean,
  p_sd
) {
  p_mean <- pmin(
    pmax(p_mean, EPS),
    1 - EPS
  )

  pred <- as.integer(
    p_mean >= CLASSIFICATION_THRESHOLD
  )

  tp <- sum(pred == 1L & y == 1L)
  tn <- sum(pred == 0L & y == 0L)
  fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L)

  sens <- if ((tp + fn) > 0L) {
    tp / (tp + fn)
  } else {
    NA_real_
  }

  spec <- if ((tn + fp) > 0L) {
    tn / (tn + fp)
  } else {
    NA_real_
  }

  cal <- calibration_parameters(
    y,
    p_mean
  )

  risk_sd <- selective_risk(
    y,
    p_mean,
    p_sd,
    "posterior_sd"
  )

  risk_margin <- selective_risk(
    y,
    p_mean,
    p_sd,
    "confidence_margin"
  )

  tibble::tibble(
    Accuracy = mean(pred == y),
    Sensitivity = sens,
    Specificity = spec,
    AUROC = safe_auc(y, p_mean),
    Brier = mean((p_mean - y)^2),
    Calibration_Intercept = cal[["intercept"]],
    Calibration_Slope = cal[["slope"]],
    Posterior_SD_Mean = mean(p_sd),
    Risk80_Posterior_SD = risk_sd,
    Risk80_Confidence_Margin = risk_margin,
    Risk80_SD_Minus_Margin =
      risk_sd - risk_margin
  )
}

# ------------------------------------------------
# 7. Analyze one paired replicate
# ------------------------------------------------

analyze_one_rep <- function(rep_id) {
  source <- make_source_sample(rep_id)
  targets <- make_target_conditions(rep_id)

  source_df <- source$data

  preprocessing <- fit_preprocess(
    source_df,
    FEATURES
  )

  x_source <- apply_preprocess(
    source_df,
    FEATURES,
    preprocessing
  )

  y_source <- as.integer(
    source_df$Y_observed
  )

  target_order <- c(
    "no_shift",
    "prevalence_shift",
    "missingness_shift",
    "combined_shift"
  )

  target_matrices <- vector(
    "list",
    length(target_order)
  )
  target_y <- vector(
    "list",
    length(target_order)
  )
  target_manifest <- vector(
    "list",
    length(target_order)
  )

  names(target_matrices) <- target_order
  names(target_y) <- target_order
  names(target_manifest) <- target_order

  for (target_name in target_order) {
    tdat <- targets[[target_name]]

    df <- as.data.frame(tdat$X)
    names(df) <- FEATURES

    target_matrices[[target_name]] <-
      apply_preprocess(
        df,
        FEATURES,
        preprocessing
      )

    target_y[[target_name]] <-
      as.integer(tdat$y)

    target_manifest[[target_name]] <-
      tibble::tibble(
        Replicate = rep_id,
        Target = target_name,
        N = length(tdat$y),
        True_Prevalence = mean(tdat$y),
        Neutrophil_Missing_Rate =
          tdat$missing_rate
      )
  }

  x_target_all <- do.call(
    rbind,
    target_matrices
  )

  target_sizes <- vapply(
    target_matrices,
    nrow,
    integer(1)
  )

  target_end <- cumsum(target_sizes)
  target_start <- c(
    1L,
    head(target_end, -1L) + 1L
  )

  metrics_list <- list()
  fit_status_list <- list()
  row_counter <- 1L

  for (model_index in seq_along(MODEL_SET)) {
    model_name <- MODEL_SET[[model_index]]

    fit_seed <- (
      BASE_SEED +
        rep_id * 10000L +
        model_index * 100L +
        301L
    )

    set.seed(fit_seed)

    start_time <- Sys.time()

    fit_result <- tryCatch(
      fit_model(
        model_name,
        x_source,
        y_source,
        x_target_all
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
      fit_status_list[[model_index]] <-
        tibble::tibble(
          Replicate = rep_id,
          Model = model_name,
          Status = "FAILED",
          Failure_Message =
            fit_result$message,
          Runtime_Seconds =
            runtime_seconds
        )
      next
    }

    draws <- validate_prediction_draws(
      fit_result$draws,
      expected_draws =
        expected_draw_count(
          model_name
        ),
      expected_test_size =
        nrow(x_target_all),
      model_name =
        model_name
    )

    fit_status_list[[model_index]] <-
      tibble::tibble(
        Replicate = rep_id,
        Model = model_name,
        Status = "OK",
        Failure_Message = NA_character_,
        Runtime_Seconds = runtime_seconds
      )

    for (target_index in seq_along(target_order)) {
      target_name <- target_order[[target_index]]
      columns <- target_start[[target_index]]:
        target_end[[target_index]]

      target_draws <- draws[
        ,
        columns,
        drop = FALSE
      ]

      summary <- summarise_draws(
        target_draws
      )

      met <- calc_target_metrics(
        y = target_y[[target_name]],
        p_mean = summary$p_mean,
        p_sd = summary$p_sd
      ) |>
        dplyr::mutate(
          Replicate = rep_id,
          Model = model_name,
          Target = target_name,
          Source_Observed_Prevalence =
            mean(y_source),
          Source_Neutrophil_Missing_Rate =
            source$missing_rate,
          Target_True_Prevalence =
            mean(target_y[[target_name]]),
          Target_Neutrophil_Missing_Rate =
            targets[[target_name]]$missing_rate,
          .before = 1
        )

      metrics_list[[row_counter]] <- met
      row_counter <- row_counter + 1L
    }

    rm(fit_result, draws)
    gc(FALSE)
  }

  list(
    metrics = dplyr::bind_rows(metrics_list),
    fit_status =
      dplyr::bind_rows(fit_status_list),
    target_manifest =
      dplyr::bind_rows(
        target_manifest
      )
  )
}

# ------------------------------------------------
# 8. Run / resume
# ------------------------------------------------

checkpoint_file <- file.path(
  CHECKPOINT_DIR,
  "paired_shift_checkpoint.rds"
)

if (
  RESUME_ANALYSIS &&
  file.exists(checkpoint_file)
) {
  checkpoint <- readRDS(
    checkpoint_file
  )

  all_metrics <- checkpoint$all_metrics
  all_status <- checkpoint$all_status
  all_target_manifest <-
    checkpoint$all_target_manifest
  completed_reps <-
    checkpoint$completed_reps

  cat(
    "Resuming after replicate ",
    completed_reps,
    ".\n",
    sep = ""
  )
} else {
  all_metrics <- list()
  all_status <- list()
  all_target_manifest <- list()
  completed_reps <- 0L
}

if (completed_reps < N_REPS) {
  for (rep_id in seq.int(
    completed_reps + 1L,
    N_REPS
  )) {
    cat(
      "\nReplicate ",
      rep_id,
      " of ",
      N_REPS,
      "\n",
      sep = ""
    )

    result <- analyze_one_rep(
      rep_id
    )

    all_metrics[[rep_id]] <-
      result$metrics
    all_status[[rep_id]] <-
      result$fit_status
    all_target_manifest[[rep_id]] <-
      result$target_manifest

    if (
      rep_id %% CHECKPOINT_EVERY == 0L ||
      rep_id == N_REPS
    ) {
      saveRDS(
        list(
          all_metrics =
            all_metrics,
          all_status =
            all_status,
          all_target_manifest =
            all_target_manifest,
          completed_reps =
            rep_id
        ),
        checkpoint_file,
        compress = "xz"
      )
    }
  }
}

metrics <- dplyr::bind_rows(
  all_metrics
)
fit_status <- dplyr::bind_rows(
  all_status
)
target_manifest <- dplyr::bind_rows(
  all_target_manifest
)

# Stop if any requested model fit failed.
failed <- fit_status |>
  dplyr::filter(Status != "OK")

if (nrow(failed) > 0L) {
  safe_path <- file.path(
    TABLE_DIR,
    "model_fit_failures.csv"
  )
  data.table::fwrite(
    failed,
    safe_path
  )

  warning(
    "At least one model fit failed. ",
    "See: ",
    safe_path
  )
}

# ------------------------------------------------
# 9. Paired within-replicate shift contrasts
# ------------------------------------------------

metric_names <- c(
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "AUROC",
  "Brier",
  "Calibration_Intercept",
  "Calibration_Slope",
  "Posterior_SD_Mean",
  "Risk80_Posterior_SD",
  "Risk80_Confidence_Margin",
  "Risk80_SD_Minus_Margin"
)

long_metrics <- metrics |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(
      metric_names
    ),
    names_to = "Metric",
    values_to = "Value"
  )

baseline <- long_metrics |>
  dplyr::filter(
    Target == "no_shift"
  ) |>
  dplyr::select(
    Replicate,
    Model,
    Metric,
    Baseline_Value = Value
  )

paired_contrasts <- long_metrics |>
  dplyr::filter(
    Target != "no_shift"
  ) |>
  dplyr::left_join(
    baseline,
    by = c(
      "Replicate",
      "Model",
      "Metric"
    )
  ) |>
  dplyr::mutate(
    Delta_From_No_Shift =
      Value - Baseline_Value
  )

summarise_vector <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(
      tibble::tibble(
        N = 0L,
        Mean = NA_real_,
        SD = NA_real_,
        Median = NA_real_,
        Q025 = NA_real_,
        Q975 = NA_real_
      )
    )
  }

  tibble::tibble(
    N = length(x),
    Mean = mean(x),
    SD = if (length(x) > 1L) sd(x) else NA_real_,
    Median = median(x),
    Q025 = as.numeric(
      quantile(
        x,
        0.025,
        names = FALSE
      )
    ),
    Q975 = as.numeric(
      quantile(
        x,
        0.975,
        names = FALSE
      )
    )
  )
}

target_summary <- long_metrics |>
  dplyr::group_by(
    Model,
    Target,
    Metric
  ) |>
  dplyr::group_modify(
    ~ summarise_vector(.x$Value)
  ) |>
  dplyr::ungroup()

contrast_summary <- paired_contrasts |>
  dplyr::group_by(
    Model,
    Target,
    Metric
  ) |>
  dplyr::group_modify(
    ~ summarise_vector(
      .x$Delta_From_No_Shift
    )
  ) |>
  dplyr::ungroup()

# ------------------------------------------------
# 10. Figures
# ------------------------------------------------

target_labels <- c(
  no_shift = "No shift",
  prevalence_shift = "Prevalence shift",
  missingness_shift = "Missingness shift",
  combined_shift = "Combined shift"
)

plot_data <- metrics |>
  dplyr::mutate(
    Target_Label =
      factor(
        target_labels[Target],
        levels = unname(target_labels)
      )
  )

risk_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = Target_Label,
    y = Risk80_SD_Minus_Margin,
    group = Model,
    linetype = Model
  )
) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    size = 2
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  ggplot2::labs(
    title = "Paired source-to-target shift simulation",
    subtitle = "Positive values mean posterior-SD deferral has higher error than confidence-margin deferral",
    x = NULL,
    y = "Risk@80%: posterior SD minus confidence margin"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      ),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(
    FIGURE_DIR,
    "paired_shift_selective_risk_difference.pdf"
  ),
  risk_plot,
  width = 9,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  file.path(
    FIGURE_DIR,
    "paired_shift_selective_risk_difference.png"
  ),
  risk_plot,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

cal_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = Target_Label,
    y = Calibration_Intercept,
    group = Model,
    linetype = Model
  )
) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    size = 2
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  ggplot2::labs(
    title = "Calibration intercept under paired target shifts",
    x = NULL,
    y = "Calibration intercept"
  ) +
  ggplot2::theme_minimal(
    base_size = 12
  ) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    axis.text.x =
      ggplot2::element_text(
        angle = 20,
        hjust = 1
      ),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(
    FIGURE_DIR,
    "paired_shift_calibration_intercept.pdf"
  ),
  cal_plot,
  width = 9,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  file.path(
    FIGURE_DIR,
    "paired_shift_calibration_intercept.png"
  ),
  cal_plot,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------
# 11. Save outputs
# ------------------------------------------------

data.table::fwrite(
  metrics,
  file.path(
    TABLE_DIR,
    "01_paired_shift_replication_metrics.csv"
  )
)

data.table::fwrite(
  target_manifest,
  file.path(
    TABLE_DIR,
    "02_paired_shift_target_manifest.csv"
  )
)

data.table::fwrite(
  fit_status,
  file.path(
    TABLE_DIR,
    "03_paired_shift_model_fit_status.csv"
  )
)

data.table::fwrite(
  target_summary,
  file.path(
    TABLE_DIR,
    "04_paired_shift_target_summary.csv"
  )
)

data.table::fwrite(
  paired_contrasts,
  file.path(
    TABLE_DIR,
    "05_paired_shift_within_replication_contrasts.csv"
  )
)

data.table::fwrite(
  contrast_summary,
  file.path(
    TABLE_DIR,
    "06_paired_shift_contrast_summary.csv"
  )
)

method_manifest <- tibble::tibble(
  Item = c(
    "Bayesian model implementation",
    "RPA mean-vector source",
    "RPA class-0 covariance source",
    "RPA class-1 covariance source",
    "Replicates",
    "Source sample size",
    "Target sample size",
    "Source prevalence",
    "Target prevalence",
    "Source neutrophil missingness",
    "Target neutrophil missingness",
    "Target conditions",
    "Pairing strategy",
    "Preprocessing",
    "Models",
    "Posterior draws",
    "Primary selective coverage"
  ),
  Value = c(
    "Embedded manuscript-equivalent BLR, genuine BLDA, and BART functions",
    RPA_MU_FILE,
    RPA_SIGMA0_FILE,
    RPA_SIGMA1_FILE,
    N_REPS,
    N_SOURCE,
    N_TARGET,
    SOURCE_PREVALENCE,
    TARGET_PREVALENCE,
    SOURCE_NEUTROPHIL_MISSING,
    TARGET_NEUTROPHIL_MISSING,
    "No shift; prevalence shift; neutrophil-missingness shift; combined shift",
    "Common class-specific target pools; nested missingness uniforms within prevalence composition",
    "Fit on source only; source medians/centering/scaling applied unchanged to all targets",
    paste(MODEL_SET, collapse = ", "),
    "300 per model",
    PRIMARY_COVERAGE
  )
)

data.table::fwrite(
  method_manifest,
  file.path(
    TABLE_DIR,
    "07_method_manifest.csv"
  )
)

saveRDS(
  list(
    metrics = metrics,
    target_manifest = target_manifest,
    fit_status = fit_status,
    target_summary = target_summary,
    paired_contrasts = paired_contrasts,
    contrast_summary = contrast_summary,
    method_manifest = method_manifest
  ),
  file.path(
    OBJECT_DIR,
    "paired_source_target_shift_complete_analysis.rds"
  ),
  compress = "xz"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    OUTPUT_ROOT,
    "sessionInfo.txt"
  )
)

cat("\n============================================================\n")
cat("B-RISA paired source-to-target shift analysis completed\n")
cat("============================================================\n\n")

cat("Target summaries for selective-risk difference:\n")
print(
  target_summary |>
    dplyr::filter(
      Metric == "Risk80_SD_Minus_Margin"
    )
)

cat("\nPaired deltas from no-shift condition:\n")
print(
  contrast_summary |>
    dplyr::filter(
      Metric %in% c(
        "AUROC",
        "Brier",
        "Calibration_Intercept",
        "Risk80_SD_Minus_Margin"
      )
    )
)

cat("\nOutputs written to:\n", OUTPUT_ROOT, "\n")


############
# ================================================================
# 95% MONTE CARLO CONFIDENCE INTERVALS
# ================================================================

# ------------------------------------------------
# Helper: mean + 95% Monte Carlo CI
# ------------------------------------------------

mean_mc_ci <- function(x, conf_level = 0.95) {
  
  x <- x[is.finite(x)]
  
  n <- length(x)
  
  if (n == 0L) {
    return(
      tibble::tibble(
        N = 0L,
        Mean = NA_real_,
        SD = NA_real_,
        SE = NA_real_,
        CI_Lower = NA_real_,
        CI_Upper = NA_real_
      )
    )
  }
  
  if (n == 1L) {
    return(
      tibble::tibble(
        N = 1L,
        Mean = x,
        SD = NA_real_,
        SE = NA_real_,
        CI_Lower = NA_real_,
        CI_Upper = NA_real_
      )
    )
  }
  
  mean_x <- mean(x)
  sd_x <- stats::sd(x)
  se_x <- sd_x / sqrt(n)
  
  alpha <- 1 - conf_level
  
  t_critical <- stats::qt(
    1 - alpha / 2,
    df = n - 1L
  )
  
  tibble::tibble(
    N = n,
    Mean = mean_x,
    SD = sd_x,
    SE = se_x,
    CI_Lower = mean_x - t_critical * se_x,
    CI_Upper = mean_x + t_critical * se_x
  )
}

# ================================================================
# A. 95% CI FOR MEAN PERFORMANCE WITHIN EACH TARGET CONDITION
# ================================================================

metrics_for_ci <- c(
  "Accuracy",
  "Sensitivity",
  "Specificity",
  "AUROC",
  "Brier",
  "Calibration_Intercept",
  "Calibration_Slope",
  "Posterior_SD_Mean",
  "Risk80_Posterior_SD",
  "Risk80_Confidence_Margin",
  "Risk80_SD_Minus_Margin"
)

target_mean_95ci <- metrics |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(metrics_for_ci),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  dplyr::group_by(
    Model,
    Target,
    Metric
  ) |>
  dplyr::group_modify(
    ~ mean_mc_ci(.x$Value)
  ) |>
  dplyr::ungroup()

data.table::fwrite(
  target_mean_95ci,
  file.path(
    TABLE_DIR,
    "08_paired_shift_mean_monte_carlo_95CI.csv"
  )
)

# ================================================================
# B. 95% CI FOR PAIRED CHANGE FROM NO-SHIFT CONDITION
# ================================================================

paired_delta_95ci <- paired_contrasts |>
  dplyr::group_by(
    Model,
    Target,
    Metric
  ) |>
  dplyr::group_modify(
    ~ mean_mc_ci(
      .x$Delta_From_No_Shift
    )
  ) |>
  dplyr::ungroup()

data.table::fwrite(
  paired_delta_95ci,
  file.path(
    TABLE_DIR,
    "09_paired_shift_delta_from_no_shift_monte_carlo_95CI.csv"
  )
)

# ================================================================
# C. PRIMARY RESULT:
# POSTERIOR-SD MINUS CONFIDENCE-MARGIN SELECTIVE RISK
# ================================================================

selective_risk_95ci <- target_mean_95ci |>
  dplyr::filter(
    Metric == "Risk80_SD_Minus_Margin"
  ) |>
  dplyr::mutate(
    Interpretation = dplyr::case_when(
      CI_Lower > 0 ~
        "Posterior-SD deferral has higher mean selective risk than confidence margin",
      CI_Upper < 0 ~
        "Posterior-SD deferral has lower mean selective risk than confidence margin",
      TRUE ~
        "Monte Carlo 95% CI includes zero"
    ),
    Estimate_95CI = sprintf(
      "%.4f (95%% Monte Carlo CI %.4f to %.4f)",
      Mean,
      CI_Lower,
      CI_Upper
    )
  ) |>
  dplyr::arrange(
    Target,
    Model
  )

data.table::fwrite(
  selective_risk_95ci,
  file.path(
    TABLE_DIR,
    "10_selective_risk_SD_minus_margin_monte_carlo_95CI.csv"
  )
)

# ================================================================
# D. KEY PAIRED SHIFT EFFECTS FOR THE MANUSCRIPT
# ================================================================

key_shift_95ci <- paired_delta_95ci |>
  dplyr::filter(
    Metric %in% c(
      "AUROC",
      "Brier",
      "Calibration_Intercept",
      "Calibration_Slope",
      "Risk80_Posterior_SD",
      "Risk80_Confidence_Margin",
      "Risk80_SD_Minus_Margin"
    )
  ) |>
  dplyr::mutate(
    Estimate_95CI = sprintf(
      "%.4f (95%% Monte Carlo CI %.4f to %.4f)",
      Mean,
      CI_Lower,
      CI_Upper
    )
  ) |>
  dplyr::arrange(
    Target,
    Metric,
    Model
  )

data.table::fwrite(
  key_shift_95ci,
  file.path(
    TABLE_DIR,
    "11_key_paired_shift_effects_monte_carlo_95CI.csv"
  )
)

# ================================================================
# E. OPTIONAL PUBLICATION-READY TABLE
# ================================================================

publication_table <- selective_risk_95ci |>
  dplyr::select(
    Model,
    Target,
    N,
    Mean,
    SD,
    SE,
    CI_Lower,
    CI_Upper,
    Estimate_95CI,
    Interpretation
  )

data.table::fwrite(
  publication_table,
  file.path(
    TABLE_DIR,
    "12_publication_ready_selective_risk_95CI.csv"
  )
)

# ================================================================
# F. PRINT IMPORTANT RESULTS TO CONSOLE
# ================================================================

cat("\n")
cat("============================================================\n")
cat("95% MONTE CARLO CONFIDENCE INTERVALS COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Primary result: Risk@80% posterior SD minus confidence margin\n\n"
)

print(
  selective_risk_95ci |>
    dplyr::select(
      Model,
      Target,
      N,
      Mean,
      SE,
      CI_Lower,
      CI_Upper,
      Estimate_95CI
    )
)

cat("\n")
cat(
  "Key paired changes relative to no-shift condition:\n\n"
)

print(
  key_shift_95ci |>
    dplyr::select(
      Model,
      Target,
      Metric,
      N,
      Mean,
      SE,
      CI_Lower,
      CI_Upper,
      Estimate_95CI
    )
)

cat("\nSaved confidence-interval tables to:\n")
cat(TABLE_DIR, "\n")

