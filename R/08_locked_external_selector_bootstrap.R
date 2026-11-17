# ================================================================
# B-RISA LOCKED EXTERNAL SELECTOR BOOTSTRAP
# Repository-ready version
#
# PURPOSE
# -------
# Reproduce the exploratory/post-hoc Düsseldorf selective-prediction
# comparisons reported in Table 4 of the manuscript.
#
# IMPORTANT
# ---------
# - This script does NOT replace the fully nested RPA evaluation.
# - The primary RPA incremental-value analysis is reproduced by the
#   fully nested selector script.
# - The final external selector is fitted using RPA only.
# - No Düsseldorf outcome is used to fit, standardize, tune, or select
#   the external selector.
# - Düsseldorf outcomes are used only for exploratory evaluation after
#   the selector has been frozen.
#
# This version uses repository-relative paths only.
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

# ----------------------------------------------------------------
# 0. Repository root and relative paths
# ----------------------------------------------------------------

find_project_root <- function() {
  env_root <- Sys.getenv("BRISA_ROOT", unset = "")

  if (nzchar(env_root)) {
    root <- normalizePath(
      env_root,
      winslash = "/",
      mustWork = FALSE
    )

    if (
      dir.exists(file.path(root, "R")) &&
      dir.exists(file.path(root, "data")) &&
      dir.exists(file.path(root, "results"))
    ) {
      return(root)
    }
  }

  candidates <- unique(
    normalizePath(
      c(
        getwd(),
        dirname(getwd())
      ),
      winslash = "/",
      mustWork = FALSE
    )
  )

  for (candidate in candidates) {
    if (
      dir.exists(file.path(candidate, "R")) &&
      dir.exists(file.path(candidate, "data")) &&
      dir.exists(file.path(candidate, "results"))
    ) {
      return(candidate)
    }
  }

  stop(
    "Cannot identify the B-RISA repository root. ",
    "Run this script from the repository root or from the R/ directory, ",
    "or set environment variable BRISA_ROOT to the repository root."
  )
}

PROJECT_ROOT <- find_project_root()

RPA_ALL_CV_FILE <- file.path(
  PROJECT_ROOT,
  "results",
  "rpa",
  "rpa_standardized_all_cv_predictions.csv"
)

DUSSELDORF_LOCKED_FILE <- file.path(
  PROJECT_ROOT,
  "results",
  "dusseldorf",
  "dusseldorf_external_locked_predictions.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "analysis_outputs",
  "locked_external_selector"
)

if (!file.exists(RPA_ALL_CV_FILE)) {
  stop(
    "Missing required RPA prediction file: ",
    RPA_ALL_CV_FILE,
    "\nRun R/04_rpa_clinical_analysis.R first, or restore the committed file."
  )
}

if (!file.exists(DUSSELDORF_LOCKED_FILE)) {
  stop(
    "Missing required locked Düsseldorf prediction file: ",
    DUSSELDORF_LOCKED_FILE,
    "\nRun R/05_dusseldorf_external_validation.R first, or restore the committed file."
  )
}

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n============================================================\n")
cat("B-RISA locked external selector bootstrap\n")
cat("============================================================\n")
cat("Repository root:\n  ", PROJECT_ROOT, "\n")
cat("RPA prediction input:\n  ", RPA_ALL_CV_FILE, "\n")
cat("Düsseldorf prediction input:\n  ", DUSSELDORF_LOCKED_FILE, "\n")
cat("Output directory:\n  ", OUTPUT_DIR, "\n\n")

# ----------------------------------------------------------------
# 1. Analysis settings
# ----------------------------------------------------------------

MODEL_ORDER <- c(
  "BLR",
  "BLDA",
  "BART"
)

CLASSIFICATION_THRESHOLD <- 0.50
SELECTOR_RIDGE <- 1e-4

COVERAGE_GRID <- seq(
  0.50,
  1.00,
  by = 0.05
)

PRIMARY_COVERAGE <- 0.80

BOOTSTRAP_REPS <- 2000L
BOOTSTRAP_SEED <- 20261031L
ALPHA <- 0.05

EPS <- 1e-8

RULE_ORDER <- c(
  "Posterior SD",
  "Confidence margin",
  "Combined selector"
)

# In the analysis that generated the manuscript's Düsseldorf bootstrap,
# cohort order was c("RPA", "Dusseldorf"), so Düsseldorf had index 2.
# Keeping this offset reproduces the same deterministic bootstrap stream.
DUSSELDORF_COHORT_INDEX_FOR_SEED <- 2L

# ----------------------------------------------------------------
# 2. Packages
# ----------------------------------------------------------------

required_packages <- c(
  "data.table",
  "dplyr",
  "tibble"
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
    "Missing required R packages: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ".\nInstall them first with:\ninstall.packages(c(",
    paste0(
      '"',
      missing_packages,
      '"',
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

# ----------------------------------------------------------------
# 3. General helpers
# ----------------------------------------------------------------

save_csv <- function(x, filename) {
  path <- file.path(
    OUTPUT_DIR,
    filename
  )

  data.table::fwrite(
    as.data.frame(x),
    path,
    na = "NA"
  )

  invisible(path)
}

clip_prob <- function(
  p,
  eps = EPS
) {
  pmin(
    pmax(
      as.numeric(p),
      eps
    ),
    1 - eps
  )
}

safe_sd <- function(x) {
  s <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(s) ||
    s <= 1e-12
  ) {
    1
  } else {
    s
  }
}

softplus_stable <- function(x) {
  pmax(
    x,
    0
  ) +
    log1p(
      exp(
        -abs(x)
      )
    )
}

# Ridge-penalized logistic regression.
# Intercept is not penalized.
fit_penalized_logistic <- function(
  X,
  y,
  lambda = SELECTOR_RIDGE
) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  y <- as.numeric(y)

  if (
    nrow(X) != length(y)
  ) {
    stop(
      "X/y dimension mismatch in selector fit."
    )
  }

  if (
    length(
      unique(y)
    ) < 2L
  ) {
    stop(
      "Selector training data contain only one error class."
    )
  }

  Xd <- cbind(
    "(Intercept)" = 1,
    X
  )

  p <- ncol(Xd)

  objective <- function(beta) {
    eta <- as.vector(
      Xd %*% beta
    )

    nll <- sum(
      softplus_stable(eta) -
        y * eta
    )

    penalty <- 0.5 *
      lambda *
      sum(
        beta[-1L]^2
      )

    nll + penalty
  }

  gradient <- function(beta) {
    eta <- as.vector(
      Xd %*% beta
    )

    prob <- plogis(eta)

    g <- as.vector(
      crossprod(
        Xd,
        prob - y
      )
    )

    g[-1L] <- g[-1L] +
      lambda *
      beta[-1L]

    g
  }

  fit <- stats::optim(
    par = rep(
      0,
      p
    ),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      maxit = 2000L,
      reltol = 1e-10
    ),
    hessian = TRUE
  )

  if (
    fit$convergence != 0L
  ) {
    warning(
      "Selector optimizer convergence code: ",
      fit$convergence
    )
  }

  list(
    coefficients = setNames(
      fit$par,
      colnames(Xd)
    ),
    convergence = fit$convergence,
    value = fit$value,
    lambda = lambda
  )
}

predict_penalized_logistic <- function(
  fit,
  X
) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  Xd <- cbind(
    "(Intercept)" = 1,
    X
  )

  beta <- fit$coefficients

  if (
    !identical(
      colnames(Xd),
      names(beta)
    )
  ) {
    stop(
      "Selector prediction columns do not match fitted coefficients."
    )
  }

  clip_prob(
    plogis(
      as.vector(
        Xd %*% beta
      )
    )
  )
}

order_for_rule <- function(
  data,
  rule
) {
  if (
    rule == "Posterior SD"
  ) {
    order(
      data$posterior_sd,
      data$Row_ID
    )
  } else if (
    rule == "Confidence margin"
  ) {
    order(
      -data$confidence_margin,
      data$Row_ID
    )
  } else if (
    rule == "Combined selector"
  ) {
    order(
      data$combined_error_risk,
      data$Row_ID
    )
  } else {
    stop(
      "Unknown selection rule: ",
      rule
    )
  }
}

risk_at_coverage <- function(
  data,
  rule,
  coverage
) {
  n <- nrow(data)

  if (
    n == 0L
  ) {
    return(
      NA_real_
    )
  }

  keep_n <- max(
    1L,
    min(
      n,
      ceiling(
        coverage * n
      )
    )
  )

  ord <- order_for_rule(
    data,
    rule
  )

  retained <- data[
    ord[
      seq_len(
        keep_n
      )
    ],
    ,
    drop = FALSE
  ]

  mean(
    retained$error
  )
}

selection_summary_at_coverage <- function(
  data,
  rule,
  coverage
) {
  n <- nrow(data)

  keep_n <- max(
    1L,
    min(
      n,
      ceiling(
        coverage * n
      )
    )
  )

  ord <- order_for_rule(
    data,
    rule
  )

  keep_idx <- ord[
    seq_len(
      keep_n
    )
  ]

  defer_idx <- setdiff(
    seq_len(n),
    keep_idx
  )

  retained <- data[
    keep_idx,
    ,
    drop = FALSE
  ]

  deferred <- data[
    defer_idx,
    ,
    drop = FALSE
  ]

  total_fn <- sum(
    data$Y_true == 1L &
      data$pred_class == 0L
  )

  total_fp <- sum(
    data$Y_true == 0L &
      data$pred_class == 1L
  )

  fn_deferred <- sum(
    deferred$Y_true == 1L &
      deferred$pred_class == 0L
  )

  fp_deferred <- sum(
    deferred$Y_true == 0L &
      deferred$pred_class == 1L
  )

  tibble::tibble(
    Rule = rule,
    Coverage =
      keep_n /
      n,
    N_Total = n,
    N_Retained =
      nrow(retained),
    N_Deferred =
      nrow(deferred),
    Retained_Errors =
      sum(
        retained$error
      ),
    Deferred_Errors =
      sum(
        deferred$error
      ),
    Selective_Risk =
      mean(
        retained$error
      ),
    Baseline_Error_Risk =
      mean(
        data$error
      ),
    Total_FN =
      total_fn,
    FN_Deferred =
      fn_deferred,
    FN_Remaining_Retained =
      total_fn -
      fn_deferred,
    FN_Capture_Among_Deferred =
      if (
        total_fn > 0L
      ) {
        fn_deferred /
          total_fn
      } else {
        NA_real_
      },
    Total_FP =
      total_fp,
    FP_Deferred =
      fp_deferred,
    FP_Remaining_Retained =
      total_fp -
      fp_deferred
  )
}

make_risk_curve <- function(data) {
  dplyr::bind_rows(
    lapply(
      RULE_ORDER,
      function(rule) {
        dplyr::bind_rows(
          lapply(
            COVERAGE_GRID,
            function(cov) {
              tibble::tibble(
                Rule = rule,
                Coverage = cov,
                Selective_Risk =
                  risk_at_coverage(
                    data,
                    rule,
                    cov
                  )
              )
            }
          )
        )
      }
    )
  )
}

trapezoid_auc <- function(
  x,
  y
) {
  ord <- order(x)

  x <- x[ord]
  y <- y[ord]

  if (
    length(x) < 2L
  ) {
    return(
      NA_real_
    )
  }

  sum(
    diff(x) *
      (
        head(
          y,
          -1L
        ) +
          tail(
            y,
            -1L
          )
      ) /
      2
  )
}

summarise_aurc <- function(
  curve_data
) {
  curve_data |>
    dplyr::group_by(
      Rule
    ) |>
    dplyr::summarise(
      AURC_50_100 =
        trapezoid_auc(
          Coverage,
          Selective_Risk
        ),
      nAURC_50_100 =
        AURC_50_100 /
        (
          max(Coverage) -
            min(Coverage)
        ),
      .groups = "drop"
    )
}

# ----------------------------------------------------------------
# 4. Read and validate RPA repeated held-out predictions
# ----------------------------------------------------------------

rpa_cv <- data.table::fread(
  RPA_ALL_CV_FILE,
  data.table = FALSE
)

required_rpa <- c(
  "Repeat",
  "Fold",
  "Row_ID",
  "Model",
  "Y_true",
  "p_mean",
  "p_sd"
)

missing_rpa <- setdiff(
  required_rpa,
  names(rpa_cv)
)

if (
  length(missing_rpa) > 0L
) {
  stop(
    "RPA prediction file is missing columns: ",
    paste(
      missing_rpa,
      collapse = ", "
    )
  )
}

rpa_cv <- rpa_cv |>
  dplyr::mutate(
    Model =
      as.character(Model),
    Row_ID =
      as.character(Row_ID),
    Y_true =
      as.integer(Y_true),
    p_mean =
      as.numeric(p_mean),
    p_sd =
      as.numeric(p_sd)
  ) |>
  dplyr::filter(
    Model %in%
      MODEL_ORDER
  )

repeat_audit <- rpa_cv |>
  dplyr::count(
    Model,
    Row_ID,
    Repeat,
    name = "n"
  )

if (
  min(
    repeat_audit$n
  ) != 1L ||
  max(
    repeat_audit$n
  ) != 1L
) {
  stop(
    "RPA repeated predictions are not one row ",
    "per patient/model/repetition."
  )
}

# Transport-compatible patient-level RPA summary.
# This excludes between-repeat variation:
# posterior_sd = sqrt(mean(within-fit posterior SD^2)).
rpa_patient <- rpa_cv |>
  dplyr::group_by(
    Row_ID,
    Model,
    Y_true
  ) |>
  dplyr::summarise(
    N_Heldout_Repetitions =
      dplyr::n(),
    p_mean =
      mean(p_mean),
    posterior_sd =
      sqrt(
        mean(
          p_sd^2
        )
      ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    pred_class =
      as.integer(
        p_mean >=
          CLASSIFICATION_THRESHOLD
      ),
    error =
      as.integer(
        pred_class !=
          Y_true
      ),
    confidence_margin =
      abs(
        p_mean -
          CLASSIFICATION_THRESHOLD
      ),
    ambiguity =
      CLASSIFICATION_THRESHOLD -
      confidence_margin
  )

if (
  any(
    rpa_patient$N_Heldout_Repetitions !=
      50L
  )
) {
  stop(
    "At least one RPA patient/model does not have ",
    "exactly 50 held-out repetitions."
  )
}

cat(
  "RPA patient-level selector features created: ",
  nrow(rpa_patient),
  " rows\n"
)

# ----------------------------------------------------------------
# 5. Read and validate locked Düsseldorf predictions
# ----------------------------------------------------------------

duss <- data.table::fread(
  DUSSELDORF_LOCKED_FILE,
  data.table = FALSE
)

required_duss <- c(
  "Row_ID",
  "Model",
  "Y_true",
  "p_mean",
  "p_sd"
)

missing_duss <- setdiff(
  required_duss,
  names(duss)
)

if (
  length(missing_duss) > 0L
) {
  stop(
    "Düsseldorf prediction file is missing columns: ",
    paste(
      missing_duss,
      collapse = ", "
    )
  )
}

if (
  "Analysis_Set" %in%
    names(duss)
) {
  duss <- duss |>
    dplyr::filter(
      as.character(
        Analysis_Set
      ) ==
        "primary_imputed_all"
    )
}

duss <- duss |>
  dplyr::mutate(
    Model =
      as.character(Model),
    Row_ID =
      as.character(Row_ID),
    Y_true =
      as.integer(Y_true),
    p_mean =
      as.numeric(p_mean),
    posterior_sd =
      as.numeric(p_sd)
  ) |>
  dplyr::filter(
    Model %in%
      MODEL_ORDER
  ) |>
  dplyr::distinct(
    Model,
    Row_ID,
    .keep_all = TRUE
  ) |>
  dplyr::mutate(
    pred_class =
      as.integer(
        p_mean >=
          CLASSIFICATION_THRESHOLD
      ),
    error =
      as.integer(
        pred_class !=
          Y_true
      ),
    confidence_margin =
      abs(
        p_mean -
          CLASSIFICATION_THRESHOLD
      ),
    ambiguity =
      CLASSIFICATION_THRESHOLD -
      confidence_margin
  )

duss_audit <- duss |>
  dplyr::count(
    Row_ID,
    name = "Model_Count"
  )

if (
  any(
    duss_audit$Model_Count !=
      length(MODEL_ORDER)
  )
) {
  stop(
    "Each Düsseldorf patient must have exactly ",
    "one prediction per model."
  )
}

cat(
  "Locked Düsseldorf predictions loaded: ",
  nrow(duss),
  " rows\n\n"
)

# ----------------------------------------------------------------
# 6. Fit final RPA-only selectors
# ----------------------------------------------------------------

fit_final_selector_one_model <- function(
  model_data,
  model_name
) {
  model_data <- model_data |>
    dplyr::arrange(
      Row_ID
    )

  mu_amb <- mean(
    model_data$ambiguity
  )

  sd_amb <- safe_sd(
    model_data$ambiguity
  )

  mu_sd <- mean(
    model_data$posterior_sd
  )

  sd_sd <- safe_sd(
    model_data$posterior_sd
  )

  X_final <- cbind(
    ambiguity_z =
      (
        model_data$ambiguity -
          mu_amb
      ) /
      sd_amb,
    posterior_sd_z =
      (
        model_data$posterior_sd -
          mu_sd
      ) /
      sd_sd
  )

  fit <- fit_penalized_logistic(
    X_final,
    model_data$error
  )

  list(
    model = model_name,
    data = model_data,
    fit = fit,
    scaler = list(
      ambiguity_mean =
        mu_amb,
      ambiguity_sd =
        sd_amb,
      posterior_sd_mean =
        mu_sd,
      posterior_sd_sd =
        sd_sd
    )
  )
}

final_selectors <- lapply(
  MODEL_ORDER,
  function(model_name) {
    fit_final_selector_one_model(
      rpa_patient |>
        dplyr::filter(
          Model ==
            model_name
        ),
      model_name
    )
  }
)

names(
  final_selectors
) <- MODEL_ORDER

selector_manifest <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(model_name) {
      obj <-
        final_selectors[[model_name]]

      coef <-
        obj$fit$coefficients

      tibble::tibble(
        Model =
          model_name,
        Selector_Training_Cohort =
          "RPA only",
        Selector_Training_N =
          nrow(
            obj$data
          ),
        Ambiguity_Mean_RPA =
          obj$scaler$ambiguity_mean,
        Ambiguity_SD_RPA =
          obj$scaler$ambiguity_sd,
        Posterior_SD_Mean_RPA =
          obj$scaler$posterior_sd_mean,
        Posterior_SD_SD_RPA =
          obj$scaler$posterior_sd_sd,
        Intercept =
          coef[["(Intercept)"]],
        Ambiguity_Coefficient =
          coef[["ambiguity_z"]],
        Posterior_SD_Coefficient =
          coef[["posterior_sd_z"]],
        Ridge_Penalty =
          SELECTOR_RIDGE,
        Optimizer_Convergence =
          obj$fit$convergence
      )
    }
  )
)

if (
  any(
    selector_manifest$Optimizer_Convergence !=
      0L
  )
) {
  warning(
    "At least one final RPA selector did not ",
    "return convergence code 0."
  )
}

# ----------------------------------------------------------------
# 7. Apply frozen RPA selectors to Düsseldorf
# ----------------------------------------------------------------

duss_selector_predictions <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(model_name) {
      obj <-
        final_selectors[[model_name]]

      d <- duss |>
        dplyr::filter(
          Model ==
            model_name
        ) |>
        dplyr::arrange(
          Row_ID
        )

      X_external <- cbind(
        ambiguity_z =
          (
            d$ambiguity -
              obj$scaler$ambiguity_mean
          ) /
          obj$scaler$ambiguity_sd,
        posterior_sd_z =
          (
            d$posterior_sd -
              obj$scaler$posterior_sd_mean
          ) /
          obj$scaler$posterior_sd_sd
      )

      d$combined_error_risk <-
        predict_penalized_logistic(
          obj$fit,
          X_external
        )

      d
    }
  )
)

# ----------------------------------------------------------------
# 8. Düsseldorf selective-risk summaries
# ----------------------------------------------------------------

primary_coverage_table <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(model_name) {
      d <-
        duss_selector_predictions |>
        dplyr::filter(
          Model ==
            model_name
        )

      dplyr::bind_rows(
        lapply(
          RULE_ORDER,
          function(rule) {
            selection_summary_at_coverage(
              d,
              rule,
              PRIMARY_COVERAGE
            )
          }
        )
      ) |>
        dplyr::mutate(
          Cohort =
            "Dusseldorf",
          Model =
            model_name,
          .before = 1
        )
    }
  )
)

risk_curve_table <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(model_name) {
      d <-
        duss_selector_predictions |>
        dplyr::filter(
          Model ==
            model_name
        )

      make_risk_curve(
        d
      ) |>
        dplyr::mutate(
          Cohort =
            "Dusseldorf",
          Model =
            model_name,
          .before = 1
        )
    }
  )
)

aurc_table <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(model_name) {
      summarise_aurc(
        risk_curve_table |>
          dplyr::filter(
            Model ==
              model_name
          )
      ) |>
        dplyr::mutate(
          Cohort =
            "Dusseldorf",
          Model =
            model_name,
          .before = 1
        )
    }
  )
)

# ----------------------------------------------------------------
# 9. Paired Düsseldorf patient bootstrap
# ----------------------------------------------------------------

bootstrap_one_model <- function(
  data,
  model_name,
  B = BOOTSTRAP_REPS
) {
  model_index <- match(
    model_name,
    MODEL_ORDER
  )

  # Exact deterministic seed stream used for Düsseldorf in the
  # analysis that generated the manuscript values.
  set.seed(
    BOOTSTRAP_SEED +
      DUSSELDORF_COHORT_INDEX_FOR_SEED *
        10000L +
      model_index *
        1000L
  )

  observed <- c(
    SD_minus_Margin =
      risk_at_coverage(
        data,
        "Posterior SD",
        PRIMARY_COVERAGE
      ) -
      risk_at_coverage(
        data,
        "Confidence margin",
        PRIMARY_COVERAGE
      ),

    Combined_minus_Margin =
      risk_at_coverage(
        data,
        "Combined selector",
        PRIMARY_COVERAGE
      ) -
      risk_at_coverage(
        data,
        "Confidence margin",
        PRIMARY_COVERAGE
      ),

    SD_minus_Combined =
      risk_at_coverage(
        data,
        "Posterior SD",
        PRIMARY_COVERAGE
      ) -
      risk_at_coverage(
        data,
        "Combined selector",
        PRIMARY_COVERAGE
      )
  )

  boot <- matrix(
    NA_real_,
    nrow = B,
    ncol =
      length(observed),
    dimnames = list(
      NULL,
      names(observed)
    )
  )

  n <- nrow(data)

  for (
    b in seq_len(B)
  ) {
    index <- sample.int(
      n,
      size = n,
      replace = TRUE
    )

    db <- data[
      index,
      ,
      drop = FALSE
    ]

    # Unique IDs are used only for deterministic tie breaking.
    db$Row_ID <- paste0(
      db$Row_ID,
      "_boot_",
      seq_len(
        nrow(db)
      )
    )

    r_sd <- risk_at_coverage(
      db,
      "Posterior SD",
      PRIMARY_COVERAGE
    )

    r_margin <- risk_at_coverage(
      db,
      "Confidence margin",
      PRIMARY_COVERAGE
    )

    r_combined <- risk_at_coverage(
      db,
      "Combined selector",
      PRIMARY_COVERAGE
    )

    boot[
      b,
    ] <- c(
      SD_minus_Margin =
        r_sd -
        r_margin,

      Combined_minus_Margin =
        r_combined -
        r_margin,

      SD_minus_Combined =
        r_sd -
        r_combined
    )
  }

  dplyr::bind_rows(
    lapply(
      names(observed),
      function(
        comparison_name
      ) {
        values <-
          boot[
            ,
            comparison_name
          ]

        obs <-
          unname(
            observed[[comparison_name]]
          )

        ci <- stats::quantile(
          values,
          probs = c(
            ALPHA / 2,
            1 -
              ALPHA / 2
          ),
          names = FALSE,
          type = 8,
          na.rm = TRUE
        )

        centered <-
          values -
          mean(
            values,
            na.rm = TRUE
          )

        valid <-
          is.finite(
            centered
          )

        extreme_count <- sum(
          abs(
            centered[
              valid
            ]
          ) >=
            abs(obs)
        )

        p_value <- (
          extreme_count +
            1
        ) / (
          sum(valid) +
            1
        )

        tibble::tibble(
          Cohort =
            "Dusseldorf",
          Model =
            model_name,
          Comparison =
            comparison_name,
          Difference =
            obs,
          CI_Lower =
            ci[[1L]],
          CI_Upper =
            ci[[2L]],
          P_Unadjusted =
            p_value
        )
      }
    )
  )
}

bootstrap_summary <- dplyr::bind_rows(
  lapply(
    MODEL_ORDER,
    function(
      model_name
    ) {
      bootstrap_one_model(
        duss_selector_predictions |>
          dplyr::filter(
            Model ==
              model_name
          ),
        model_name
      )
    }
  )
) |>
  dplyr::group_by(
    Cohort,
    Model
  ) |>
  dplyr::mutate(
    P_Holm =
      p.adjust(
        P_Unadjusted,
        method = "holm"
      )
  ) |>
  dplyr::ungroup()

# ----------------------------------------------------------------
# 10. Manuscript-value audit
# ----------------------------------------------------------------

expected_primary <- tibble::tribble(
  ~Model, ~Risk_SD, ~Risk_Margin, ~Risk_Combined,
  "BLR",  0.2116, 0.1950, 0.1909,
  "BLDA", 0.2075, 0.1867, 0.1909,
  "BART", 0.2780, 0.2199, 0.2241
)

expected_bootstrap <- tibble::tribble(
  ~Model, ~Comparison, ~Difference, ~CI_Lower, ~CI_Upper, ~P_Holm,

  "BLR",
  "SD_minus_Margin",
  0.0166,
  -0.0166,
  0.0539,
  0.756,

  "BLR",
  "Combined_minus_Margin",
  -0.0041,
  -0.0083,
  0.0124,
  0.756,

  "BLDA",
  "SD_minus_Margin",
  0.0207,
  -0.0083,
  0.0664,
  0.801,

  "BLDA",
  "Combined_minus_Margin",
  0.0041,
  -0.0124,
  0.0166,
  0.801,

  "BART",
  "SD_minus_Margin",
  0.0581,
  0.0083,
  0.0954,
  0.0210,

  "BART",
  "Combined_minus_Margin",
  0.0041,
  -0.0124,
  0.0193,
  0.652
)

get_primary_risk <- function(
  model_name,
  rule_name
) {
  x <- primary_coverage_table |>
    dplyr::filter(
      Model ==
        model_name,
      Rule ==
        rule_name
    )

  if (
    nrow(x) != 1L
  ) {
    stop(
      "Unexpected number of primary-risk rows for ",
      model_name,
      " / ",
      rule_name
    )
  }

  x$Selective_Risk[[1L]]
}

primary_check <- tibble::tibble(
  Model =
    MODEL_ORDER,

  Risk_SD =
    vapply(
      MODEL_ORDER,
      function(m) {
        get_primary_risk(
          m,
          "Posterior SD"
        )
      },
      numeric(1)
    ),

  Risk_Margin =
    vapply(
      MODEL_ORDER,
      function(m) {
        get_primary_risk(
          m,
          "Confidence margin"
        )
      },
      numeric(1)
    ),

  Risk_Combined =
    vapply(
      MODEL_ORDER,
      function(m) {
        get_primary_risk(
          m,
          "Combined selector"
        )
      },
      numeric(1)
    )
)

round4 <- function(x) {
  round(
    as.numeric(x),
    4
  )
}

primary_pass <- all(
  round4(
    primary_check$Risk_SD
  ) ==
    expected_primary$Risk_SD
) &&
  all(
    round4(
      primary_check$Risk_Margin
    ) ==
      expected_primary$Risk_Margin
  ) &&
  all(
    round4(
      primary_check$Risk_Combined
    ) ==
      expected_primary$Risk_Combined
  )

bootstrap_for_check <-
  bootstrap_summary |>
  dplyr::filter(
    Comparison %in%
      c(
        "SD_minus_Margin",
        "Combined_minus_Margin"
      )
  ) |>
  dplyr::arrange(
    match(
      Model,
      MODEL_ORDER
    ),
    match(
      Comparison,
      c(
        "SD_minus_Margin",
        "Combined_minus_Margin"
      )
    )
  )

expected_bootstrap_ordered <-
  expected_bootstrap |>
  dplyr::arrange(
    match(
      Model,
      MODEL_ORDER
    ),
    match(
      Comparison,
      c(
        "SD_minus_Margin",
        "Combined_minus_Margin"
      )
    )
  )

# Manuscript reports differences and confidence limits to four decimals.
# Holm-adjusted p values are generally reported to three decimals
# (BART SD-minus-margin appears as 0.0210 in Table 4, which is the same
# value at the reported precision). Compare p values with a 0.0005
# tolerance so the audit reflects manuscript rounding rather than
# requiring the unrounded bootstrap value to equal the printed value.
bootstrap_pass <- (
  nrow(
    bootstrap_for_check
  ) ==
    nrow(
      expected_bootstrap_ordered
    )
) &&
  all(
    round4(
      bootstrap_for_check$Difference
    ) ==
      expected_bootstrap_ordered$Difference
  ) &&
  all(
    round4(
      bootstrap_for_check$CI_Lower
    ) ==
      expected_bootstrap_ordered$CI_Lower
  ) &&
  all(
    round4(
      bootstrap_for_check$CI_Upper
    ) ==
      expected_bootstrap_ordered$CI_Upper
  ) &&
  all(
    abs(
      bootstrap_for_check$P_Holm -
        expected_bootstrap_ordered$P_Holm
    ) <= 0.0005
  )

audit_table <- tibble::tibble(
  Check = c(
    "Dusseldorf primary Risk@80 values",
    "Dusseldorf paired bootstrap Table 4 values"
  ),
  Passed = c(
    primary_pass,
    bootstrap_pass
  )
)

# ----------------------------------------------------------------
# 11. Save outputs
# ----------------------------------------------------------------

save_csv(
  selector_manifest,
  "01_final_rpa_selector_for_external_application.csv"
)

save_csv(
  duss_selector_predictions,
  "02_dusseldorf_frozen_selector_predictions.csv"
)

save_csv(
  primary_coverage_table,
  "03_dusseldorf_selective_risk_primary_80pct.csv"
)

save_csv(
  risk_curve_table,
  "04_dusseldorf_selective_risk_coverage.csv"
)

save_csv(
  aurc_table,
  "05_dusseldorf_selective_risk_naurc.csv"
)

save_csv(
  bootstrap_summary,
  "06_dusseldorf_paired_bootstrap_rule_comparisons_80pct.csv"
)

method_manifest <- tibble::tibble(
  Item = c(
    "RPA input",
    "Dusseldorf input",
    "RPA uncertainty used by selector",
    "Selector model",
    "Selector ridge penalty",
    "Classification threshold",
    "Primary retained coverage",
    "Coverage grid",
    "Bootstrap replicates",
    "Bootstrap seed",
    "Dusseldorf cohort-index seed offset",
    "External fitting rule",
    "Interpretation"
  ),

  Value = c(
    RPA_ALL_CV_FILE,
    DUSSELDORF_LOCKED_FILE,

    paste0(
      "sqrt(mean(within-fit posterior SD^2)) across ",
      "50 held-out repetitions; excludes between-repeat variation"
    ),

    paste0(
      "Ridge-penalized logistic error model: ",
      "error ~ ambiguity + posterior SD"
    ),

    SELECTOR_RIDGE,
    CLASSIFICATION_THRESHOLD,
    PRIMARY_COVERAGE,

    paste(
      COVERAGE_GRID,
      collapse = ", "
    ),

    BOOTSTRAP_REPS,
    BOOTSTRAP_SEED,
    DUSSELDORF_COHORT_INDEX_FOR_SEED,

    paste0(
      "Final selector fitted using RPA only; ",
      "RPA-derived standardization and coefficients ",
      "applied unchanged in Dusseldorf"
    ),

    "Exploratory/post hoc locked external incremental-value analysis"
  )
)

save_csv(
  method_manifest,
  "07_method_manifest.csv"
)

save_csv(
  audit_table,
  "08_manuscript_value_audit.csv"
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTPUT_DIR,
    "sessionInfo.txt"
  )
)

# ----------------------------------------------------------------
# 12. Console summary
# ----------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "B-RISA locked external selector bootstrap completed\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Final RPA selector coefficients:\n"
)

print(
  selector_manifest
)

cat(
  "\nDüsseldorf selective risk at 80% retained coverage:\n"
)

print(
  primary_coverage_table |>
    dplyr::select(
      Model,
      Rule,
      Coverage,
      Selective_Risk,
      Retained_Errors,
      FN_Deferred,
      FN_Capture_Among_Deferred
    )
)

cat(
  "\nDüsseldorf paired bootstrap rule comparisons:\n"
)

print(
  bootstrap_summary
)

cat(
  "\nManuscript-value audit:\n"
)

print(
  audit_table
)

if (
  !all(
    audit_table$Passed
  )
) {
  warning(
    "At least one reproduced value does not match the manuscript ",
    "after rounding. Inspect the saved audit and input files."
  )
}

cat(
  "\nOutputs written to:\n  ",
  OUTPUT_DIR,
  "\n"
)
