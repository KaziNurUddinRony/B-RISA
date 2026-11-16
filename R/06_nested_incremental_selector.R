# ==============================================================================
# B-RISA Fully Nested Incremental Deferral Selector and Audit
#
# Reproduces the fully nested RPA selector evaluation, leakage audit, locked
# Düsseldorf selector application, WBC harmonization audit, paired bootstrap
# contrasts, false-negative capture, and illustrative asymmetric-cost analyses
# reported in the manuscript.
#
# Run from the repository root, from the R/ directory, or set BRISA_ROOT
# to the repository root before running.
# ==============================================================================

rm(list = ls())
gc()
graphics.off()

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 180,
  dplyr.summarise.inform = FALSE
)

# ------------------------------------------------------------------------------
# 0. USER CONFIGURATION
# ------------------------------------------------------------------------------

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

RPA_FILE <- file.path(
  PROJECT_ROOT, "data", "derived",
  "rpa_standardized_selected_cleaned_data.csv"
)

DUS_FILE <- file.path(
  PROJECT_ROOT, "data", "raw", "app_data_ext.csv"
)

OUT_DIR <- file.path(PROJECT_ROOT, "analysis_outputs", "nested_incremental_selector")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Final analysis settings.
SEED_NESTED <- 20261115L
SEED_BOOT <- 20261116L
SEED_EXTERNAL <- 20260730L

OUTER_FOLDS <- 10L

# Set to 50L for the submission analysis so the selector-training features use
# the same repeated-CV construction as the current manuscript.
INNER_REPEATS <- 50L
INNER_FOLDS <- 10L

# Optional fast debugging mode. Set TRUE only to test the script.
DEBUG_MODE <- FALSE
if (DEBUG_MODE) {
  INNER_REPEATS <- 2L
}

MODEL_SET <- c("BLR", "BLDA", "BART")
FEATURES <- c(
  "Body_Temperature",
  "WBC_Count",
  "Neutrophil_Percentage"
)

N_DRAWS <- 300L
BART_NTREE <- 50L
BART_BURN <- 150L

CLASS_THRESHOLD <- 0.50
PRIMARY_COVERAGE <- 0.80
COVERAGE_GRID <- seq(0.50, 1.00, by = 0.05)

RIDGE_LAMBDA <- 1e-4
N_BOOT <- 2000L

# Illustrative only; these are NOT clinically established appendicitis utilities.
FN_COST_RATIOS <- c(2, 5)

# The source authors' public external-validation preprocessing uses >1000
# followed by division by 1000 for WBCCount. The audit also verifies whether
# any values fall in a range that could make an alternate threshold matter.
WBC_PRIMARY_THRESHOLD <- 1000

# RPA upstream standardization constants already documented in B-RISA.
RPA_UPSTREAM <- list(
  mean = c(
    Body_Temperature = 37.397341,
    WBC_Count = 12.763885,
    Neutrophil_Percentage = 71.836041
  ),
  sd = c(
    Body_Temperature = 0.916253,
    WBC_Count = 5.344710,
    Neutrophil_Percentage = 14.459023
  )
)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "MASS",
  "pROC",
  "BART"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages before running the final analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(MASS)
  library(pROC)
})

# ------------------------------------------------------------------------------
# 2. GENERAL HELPERS
# ------------------------------------------------------------------------------

clip_prob <- function(p, eps = 1e-8) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

safe_auc <- function(y, score) {
  y <- as.integer(y)
  score <- as.numeric(score)
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]
  score <- score[ok]
  if (length(unique(y)) < 2L) return(NA_real_)
  tryCatch(
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = y,
          predictor = score,
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        )
      )
    ),
    error = function(e) NA_real_
  )
}

trapz_base <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (length(x) < 2L) return(NA_real_)
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

save_csv <- function(x, filename) {
  data.table::fwrite(
    as.data.frame(x),
    file.path(OUT_DIR, filename),
    na = "NA"
  )
}

save_rds <- function(x, filename) {
  saveRDS(x, file.path(OUT_DIR, filename))
}

# ------------------------------------------------------------------------------
# 3. DATA LOADING AND PREPROCESSING
# ------------------------------------------------------------------------------

load_rpa <- function(path) {
  if (!file.exists(path)) stop("RPA file not found: ", path)

  dat <- data.table::fread(path, data.table = FALSE)

  required <- c("Y", FEATURES)
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop("RPA file missing columns: ", paste(missing, collapse = ", "))
  }

  dat <- dat[, required, drop = FALSE]
  dat$Y <- as.integer(dat$Y)
  dat <- dat[!is.na(dat$Y) & dat$Y %in% c(0L, 1L), , drop = FALSE]

  for (nm in FEATURES) {
    dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  }

  dat$Row_ID <- seq_len(nrow(dat))
  dat
}

fit_preprocessor <- function(train_dat) {
  med <- vapply(
    train_dat[, FEATURES, drop = FALSE],
    function(z) stats::median(as.numeric(z), na.rm = TRUE),
    numeric(1)
  )

  if (any(!is.finite(med))) stop("Non-finite training median.")

  x_imp <- train_dat[, FEATURES, drop = FALSE]
  for (j in seq_along(FEATURES)) {
    z <- as.numeric(x_imp[[j]])
    z[is.na(z)] <- med[j]
    x_imp[[j]] <- z
  }

  mu <- vapply(x_imp, mean, numeric(1))
  sig <- vapply(x_imp, stats::sd, numeric(1))
  sig[!is.finite(sig) | sig <= 1e-12] <- 1

  list(median = med, mean = mu, sd = sig)
}

apply_preprocessor <- function(dat, prep) {
  x <- dat[, FEATURES, drop = FALSE]

  for (j in seq_along(FEATURES)) {
    z <- suppressWarnings(as.numeric(x[[j]]))
    z[is.na(z)] <- prep$median[j]
    z <- (z - prep$mean[j]) / prep$sd[j]
    x[[j]] <- z
  }

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  colnames(x) <- FEATURES

  if (any(!is.finite(x))) stop("Non-finite value after preprocessing.")
  x
}

# ------------------------------------------------------------------------------
# 4. STRATIFIED FOLDS
# ------------------------------------------------------------------------------

make_stratified_folds <- function(y, k = 10L, seed = 1L) {
  y <- as.integer(y)
  if (min(table(y)) < k) {
    stop("Not enough observations in a class for ", k, " folds.")
  }

  set.seed(seed)

  idx0 <- sample(which(y == 0L))
  idx1 <- sample(which(y == 1L))

  f0 <- split(idx0, rep(seq_len(k), length.out = length(idx0)))
  f1 <- split(idx1, rep(seq_len(k), length.out = length(idx1)))

  lapply(seq_len(k), function(j) {
    test <- sort(c(f0[[j]], f1[[j]]))
    train <- setdiff(seq_along(y), test)
    list(fold = j, train = train, test = test)
  })
}

make_repeated_stratified_folds <- function(y, k, repeats, seed) {
  out <- vector("list", repeats * k)
  z <- 0L

  for (r in seq_len(repeats)) {
    fs <- make_stratified_folds(y, k = k, seed = seed + r * 1009L)
    for (j in seq_len(k)) {
      z <- z + 1L
      out[[z]] <- list(
        repeat_id = r,
        fold = j,
        train = fs[[j]]$train,
        test = fs[[j]]$test
      )
    }
  }
  out
}

# ------------------------------------------------------------------------------
# 5. FINAL B-RISA MODEL IMPLEMENTATIONS
# ------------------------------------------------------------------------------

# ----- 5.1 BLR: MAP + Laplace --------------------------------------------------

softplus_stable <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

nearest_pd <- function(M, floor = 1e-8) {
  M <- (M + t(M)) / 2
  e <- eigen(M, symmetric = TRUE)
  e$values[e$values < floor] <- floor
  out <- e$vectors %*% diag(e$values, nrow = length(e$values)) %*% t(e$vectors)
  (out + t(out)) / 2
}

predict_blr <- function(x_train, y_train, x_test, n_draws = N_DRAWS, seed = 1L) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.numeric(y_train)

  Xtr <- cbind("(Intercept)" = 1, x_train)
  Xte <- cbind("(Intercept)" = 1, x_test)

  q <- ncol(Xtr)
  prev0 <- (sum(y_train) + 0.5) / (length(y_train) + 1)

  prior_mean <- c(qlogis(prev0), rep(0, q - 1L))
  prior_sd <- rep(2.5, q)
  prior_prec <- 1 / prior_sd^2

  nlp <- function(beta) {
    eta <- as.vector(Xtr %*% beta)
    sum(softplus_stable(eta) - y_train * eta) +
      0.5 * sum(prior_prec * (beta - prior_mean)^2)
  }

  grad <- function(beta) {
    eta <- as.vector(Xtr %*% beta)
    p <- plogis(eta)
    as.vector(crossprod(Xtr, p - y_train)) +
      prior_prec * (beta - prior_mean)
  }

  opt <- stats::optim(
    par = prior_mean,
    fn = nlp,
    gr = grad,
    method = "BFGS",
    control = list(maxit = 2000L, reltol = 1e-10)
  )

  if (opt$convergence != 0L && max(abs(grad(opt$par))) > 1e-5) {
    stop("BLR optimization failed.")
  }

  pfit <- plogis(as.vector(Xtr %*% opt$par))
  w <- pmax(pfit * (1 - pfit), 1e-8)

  H <- crossprod(Xtr, Xtr * w) + diag(prior_prec, nrow = q)
  V <- nearest_pd(solve((H + t(H)) / 2))

  set.seed(seed)
  b <- MASS::mvrnorm(n = n_draws, mu = opt$par, Sigma = V)
  b <- matrix(b, nrow = n_draws, ncol = q)

  draws <- plogis(b %*% t(Xte))

  list(
    draws = draws,
    diagnostics = list(
      convergence = opt$convergence,
      max_abs_gradient = max(abs(grad(opt$par))),
      training_prevalence = mean(y_train)
    )
  )
}

# ----- 5.2 BLDA: shared-covariance Normal-Inverse-Wishart ----------------------

rinvwishart <- function(df, S) {
  solve(
    stats::rWishart(
      1L,
      df = df,
      Sigma = solve(S)
    )[, , 1L]
  )
}

predict_blda <- function(
  x_train,
  y_train,
  x_test,
  n_draws = N_DRAWS,
  seed = 1L,
  kappa0 = 0.01,
  nu_extra = 2L,
  lambda0_scale = 1,
  beta_a = 0.5,
  beta_b = 0.5
) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  y_train <- as.integer(y_train)

  p <- ncol(x_train)
  n <- nrow(x_train)

  idx0 <- which(y_train == 0L)
  idx1 <- which(y_train == 1L)
  n0 <- length(idx0)
  n1 <- length(idx1)

  if (n0 < 2L || n1 < 2L) stop("BLDA needs both classes.")

  x0 <- x_train[idx0, , drop = FALSE]
  x1 <- x_train[idx1, , drop = FALSE]

  xb0 <- colMeans(x0)
  xb1 <- colMeans(x1)

  S0 <- crossprod(sweep(x0, 2, xb0, "-"))
  S1 <- crossprod(sweep(x1, 2, xb1, "-"))

  m00 <- rep(0, p)
  m01 <- rep(0, p)

  k0n <- kappa0 + n0
  k1n <- kappa0 + n1

  m0n <- (kappa0 * m00 + n0 * xb0) / k0n
  m1n <- (kappa0 * m01 + n1 * xb1) / k1n

  nu0 <- p + nu_extra
  nu_n <- nu0 + n

  Lambda0 <- diag(lambda0_scale, p)

  Lambda_n <- Lambda0 +
    S0 + S1 +
    (kappa0 * n0 / k0n) * tcrossprod(xb0 - m00) +
    (kappa0 * n1 / k1n) * tcrossprod(xb1 - m01)

  set.seed(seed)

  pi1_draw <- stats::rbeta(
    n_draws,
    beta_a + n1,
    beta_b + n0
  )

  draws <- matrix(NA_real_, nrow = n_draws, ncol = nrow(x_test))

  for (s in seq_len(n_draws)) {
    Sigma <- rinvwishart(nu_n, Lambda_n)

    mu0 <- as.numeric(
      MASS::mvrnorm(
        1L,
        mu = m0n,
        Sigma = Sigma / k0n
      )
    )

    mu1 <- as.numeric(
      MASS::mvrnorm(
        1L,
        mu = m1n,
        Sigma = Sigma / k1n
      )
    )

    invS <- solve(Sigma)

    d0 <- as.vector(x_test %*% invS %*% mu0) -
      0.5 * as.numeric(crossprod(mu0, invS %*% mu0)) +
      log1p(-pi1_draw[s])

    d1 <- as.vector(x_test %*% invS %*% mu1) -
      0.5 * as.numeric(crossprod(mu1, invS %*% mu1)) +
      log(pi1_draw[s])

    draws[s, ] <- plogis(d1 - d0)
  }

  list(
    draws = draws,
    diagnostics = list(
      training_prevalence = mean(y_train),
      kappa0 = kappa0,
      nu0 = nu0,
      lambda0_scale = lambda0_scale
    )
  )
}

# ----- 5.3 BART ----------------------------------------------------------------

predict_bart <- function(x_train, y_train, x_test, n_draws = N_DRAWS, seed = 1L) {
  set.seed(seed)

  fit <- BART::pbart(
    x.train = as.matrix(x_train),
    y.train = as.numeric(y_train),
    x.test = as.matrix(x_test),
    ntree = BART_NTREE,
    ndpost = n_draws,
    nskip = BART_BURN,
    printevery = 100000L
  )

  if (is.null(fit$yhat.test)) stop("BART returned no yhat.test.")

  latent <- as.matrix(fit$yhat.test)
  draws <- stats::pnorm(latent)

  if (nrow(draws) == nrow(x_test) && ncol(draws) == n_draws) {
    draws <- t(draws)
  }

  if (!all(dim(draws) == c(n_draws, nrow(x_test)))) {
    stop(
      "Unexpected BART draw dimensions: ",
      paste(dim(draws), collapse = " x ")
    )
  }

  list(
    draws = draws,
    diagnostics = list(
      training_prevalence = mean(y_train),
      ntree = BART_NTREE
    )
  )
}

fit_base_model <- function(model, x_train, y_train, x_test, seed) {
  if (model == "BLR") {
    return(predict_blr(x_train, y_train, x_test, N_DRAWS, seed))
  }

  if (model == "BLDA") {
    return(predict_blda(x_train, y_train, x_test, N_DRAWS, seed))
  }

  if (model == "BART") {
    return(predict_bart(x_train, y_train, x_test, N_DRAWS, seed))
  }

  stop("Unknown model: ", model)
}

summarize_draws <- function(draws) {
  draws <- as.matrix(draws)
  tibble(
    p_mean = colMeans(draws),
    p_sd = apply(draws, 2L, stats::sd),
    p_lwr = apply(draws, 2L, stats::quantile, probs = 0.025, names = FALSE),
    p_upr = apply(draws, 2L, stats::quantile, probs = 0.975, names = FALSE)
  )
}

model_seed <- function(model, outer = 0L, repeat_id = 0L, fold = 0L, offset = 0L) {
  code <- match(model, MODEL_SET)
  as.integer(
    SEED_NESTED +
      outer * 1000000L +
      repeat_id * 10000L +
      fold * 100L +
      code +
      offset
  )
}

# ------------------------------------------------------------------------------
# 6. BASE-MODEL OOF FEATURES INSIDE AN OUTER TRAINING SET
# ------------------------------------------------------------------------------

generate_repeated_oof_features <- function(
  dat,
  model,
  repeats = INNER_REPEATS,
  k = INNER_FOLDS,
  seed = SEED_NESTED,
  outer_id = 0L,
  manifest_prefix = "inner"
) {
  n <- nrow(dat)

  p_sum <- numeric(n)
  v_sum <- numeric(n)
  count <- integer(n)

  folds <- make_repeated_stratified_folds(
    dat$Y,
    k = k,
    repeats = repeats,
    seed = seed + outer_id * 100000L
  )

  manifest <- vector("list", length(folds))

  for (z in seq_along(folds)) {
    f <- folds[[z]]

    tr <- f$train
    te <- f$test

    prep <- fit_preprocessor(dat[tr, , drop = FALSE])
    xtr <- apply_preprocessor(dat[tr, , drop = FALSE], prep)
    xte <- apply_preprocessor(dat[te, , drop = FALSE], prep)

    fit <- fit_base_model(
      model,
      xtr,
      dat$Y[tr],
      xte,
      seed = model_seed(
        model,
        outer = outer_id,
        repeat_id = f$repeat_id,
        fold = f$fold
      )
    )

    sm <- summarize_draws(fit$draws)

    p_sum[te] <- p_sum[te] + sm$p_mean
    v_sum[te] <- v_sum[te] + sm$p_sd^2
    count[te] <- count[te] + 1L

    manifest[[z]] <- tibble(
      outer_fold = outer_id,
      inner_repeat = f$repeat_id,
      inner_fold = f$fold,
      model = model,
      n_inner_train = length(tr),
      n_inner_test = length(te),
      min_row_id_train = min(dat$Row_ID[tr]),
      max_row_id_train = max(dat$Row_ID[tr]),
      min_row_id_test = min(dat$Row_ID[te]),
      max_row_id_test = max(dat$Row_ID[te])
    )
  }

  if (any(count != repeats)) {
    stop("Each patient should have exactly one OOF prediction per repeat.")
  }

  features <- tibble(
    Row_ID = dat$Row_ID,
    Y = dat$Y,
    p_mean = p_sum / count,
    p_sd = sqrt(v_sum / count)
  )

  features$margin <- abs(features$p_mean - CLASS_THRESHOLD)
  features$ambiguity <- CLASS_THRESHOLD - features$margin
  features$error <- as.integer(
    as.integer(features$p_mean >= CLASS_THRESHOLD) != features$Y
  )

  list(
    features = features,
    manifest = bind_rows(manifest)
  )
}

# ------------------------------------------------------------------------------
# 7. RIDGE ERROR-RISK SELECTORS
# ------------------------------------------------------------------------------

standardizer <- function(X) {
  mu <- colMeans(X)
  sig <- apply(X, 2L, stats::sd)
  sig[!is.finite(sig) | sig <= 1e-12] <- 1
  list(mean = mu, sd = sig)
}

apply_standardizer <- function(X, st) {
  sweep(sweep(X, 2L, st$mean, "-"), 2L, st$sd, "/")
}

fit_ridge_logistic <- function(X, y, lambda = RIDGE_LAMBDA) {
  X <- as.matrix(X)
  y <- as.integer(y)

  if (length(unique(y)) < 2L) {
    stop("Error-risk selector training outcome has only one class.")
  }

  Xd <- cbind(1, X)
  q <- ncol(Xd)

  nll <- function(beta) {
    eta <- as.vector(Xd %*% beta)
    sum(softplus_stable(eta) - y * eta) +
      0.5 * lambda * sum(beta[-1L]^2)
  }

  grad <- function(beta) {
    eta <- as.vector(Xd %*% beta)
    p <- plogis(eta)
    g <- as.vector(crossprod(Xd, p - y))
    g[-1L] <- g[-1L] + lambda * beta[-1L]
    g
  }

  init <- c(qlogis((sum(y) + 0.5) / (length(y) + 1)), rep(0, q - 1L))

  opt <- stats::optim(
    init,
    nll,
    grad,
    method = "BFGS",
    control = list(maxit = 2000L, reltol = 1e-12)
  )

  if (opt$convergence != 0L) {
    warning("Selector optimizer convergence code = ", opt$convergence)
  }

  list(coef = opt$par, convergence = opt$convergence)
}

predict_ridge_logistic <- function(fit, X) {
  X <- as.matrix(X)
  Xd <- cbind(1, X)
  plogis(as.vector(Xd %*% fit$coef))
}

fit_selector_pair <- function(train_features) {
  X_margin <- cbind(ambiguity = train_features$ambiguity)
  X_combined <- cbind(
    ambiguity = train_features$ambiguity,
    p_sd = train_features$p_sd
  )

  st_margin <- standardizer(X_margin)
  st_combined <- standardizer(X_combined)

  fit_margin <- fit_ridge_logistic(
    apply_standardizer(X_margin, st_margin),
    train_features$error
  )

  fit_combined <- fit_ridge_logistic(
    apply_standardizer(X_combined, st_combined),
    train_features$error
  )

  list(
    margin = list(fit = fit_margin, st = st_margin),
    combined = list(fit = fit_combined, st = st_combined)
  )
}

apply_selector_pair <- function(selector_fit, features) {
  Xm <- cbind(ambiguity = features$ambiguity)
  Xc <- cbind(
    ambiguity = features$ambiguity,
    p_sd = features$p_sd
  )

  tibble(
    margin_error_prob = predict_ridge_logistic(
      selector_fit$margin$fit,
      apply_standardizer(Xm, selector_fit$margin$st)
    ),
    combined_error_prob = predict_ridge_logistic(
      selector_fit$combined$fit,
      apply_standardizer(Xc, selector_fit$combined$st)
    )
  )
}

# ------------------------------------------------------------------------------
# 8. FULLY NESTED OUTER-FOLD RPA EVALUATION
# ------------------------------------------------------------------------------

run_nested_rpa <- function(rpa) {
  outer <- make_stratified_folds(
    rpa$Y,
    k = OUTER_FOLDS,
    seed = SEED_NESTED
  )

  all_predictions <- list()
  all_dependencies <- list()
  all_coefficients <- list()

  z_pred <- 0L
  z_dep <- 0L
  z_coef <- 0L

  for (o in seq_along(outer)) {
    message("Outer fold ", o, "/", length(outer))

    outer_train_idx <- outer[[o]]$train
    outer_test_idx <- outer[[o]]$test

    outer_train <- rpa[outer_train_idx, , drop = FALSE]
    outer_test <- rpa[outer_test_idx, , drop = FALSE]

    for (model in MODEL_SET) {
      message("  ", model)

      # Crucial nesting:
      # Every selector-training feature is generated using ONLY patients in
      # outer_train. No patient in outer_test can enter any inner base-model fit.
      inner <- generate_repeated_oof_features(
        outer_train,
        model = model,
        repeats = INNER_REPEATS,
        k = INNER_FOLDS,
        seed = SEED_NESTED + o * 100000L,
        outer_id = o
      )

      selector_fit <- fit_selector_pair(inner$features)

      # Fit base model once on the full outer-training set and predict only the
      # untouched outer test fold.
      prep_outer <- fit_preprocessor(outer_train)
      xtr <- apply_preprocessor(outer_train, prep_outer)
      xte <- apply_preprocessor(outer_test, prep_outer)

      base_fit <- fit_base_model(
        model,
        xtr,
        outer_train$Y,
        xte,
        seed = model_seed(model, outer = o, offset = 90000000L)
      )

      sm <- summarize_draws(base_fit$draws)

      test_features <- tibble(
        Row_ID = outer_test$Row_ID,
        Y = outer_test$Y,
        p_mean = sm$p_mean,
        p_sd = sm$p_sd
      )

      test_features$margin <- abs(test_features$p_mean - CLASS_THRESHOLD)
      test_features$ambiguity <- CLASS_THRESHOLD - test_features$margin
      test_features$error <- as.integer(
        as.integer(test_features$p_mean >= CLASS_THRESHOLD) != test_features$Y
      )

      sel_score <- apply_selector_pair(selector_fit, test_features)

      z_pred <- z_pred + 1L
      all_predictions[[z_pred]] <- bind_cols(
        test_features,
        tibble(
          model = model,
          outer_fold = o
        ),
        sel_score
      )

      # Explicit dependency audit.
      # Primary incremental-value quantity.
      intersection_count <- length(
        intersect(
          outer_test$Row_ID,
          outer_train$Row_ID
        )
      )

      z_dep <- z_dep + 1L
      all_dependencies[[z_dep]] <- tibble(
        outer_fold = o,
        model = model,
        n_outer_train = nrow(outer_train),
        n_outer_test = nrow(outer_test),
        outer_test_in_selector_training_patients = intersection_count,
        leakage_pass = intersection_count == 0L
      )

      z_coef <- z_coef + 1L
      all_coefficients[[z_coef]] <- tibble(
        outer_fold = o,
        model = model,
        selector_name = c("margin", "combined", "combined"),
        term = c("ambiguity", "ambiguity", "posterior_sd"),
        coefficient = c(
          selector_fit$margin$fit$coef[2L],
          selector_fit$combined$fit$coef[2L],
          selector_fit$combined$fit$coef[3L]
        ),
        convergence = c(
          selector_fit$margin$fit$convergence,
          selector_fit$combined$fit$convergence,
          selector_fit$combined$fit$convergence
        )
      )
    }
  }

  pred <- bind_rows(all_predictions) |>
    arrange(model, Row_ID)

  dep <- bind_rows(all_dependencies)
  coef <- bind_rows(all_coefficients)

  if (any(!dep$leakage_pass)) {
    stop("Nested dependency audit failed.")
  }

  # Every patient should occur once per model in the OUTER evaluation.
  check <- pred |>
    count(model, Row_ID, name = "n") |>
    dplyr::filter(n != 1L)

  if (nrow(check) > 0L) {
    stop("Nested outer predictions are not one-row-per-patient-per-model.")
  }

  list(predictions = pred, dependency = dep, coefficients = coef)
}

# ------------------------------------------------------------------------------
# 9. SELECTIVE PERFORMANCE
# ------------------------------------------------------------------------------

retained_indices <- function(score, coverage, higher_is_better = TRUE) {
  n <- length(score)
  keep_n <- max(1L, min(n, as.integer(ceiling(coverage * n))))
  ord <- order(score, decreasing = higher_is_better, na.last = NA)
  ord[seq_len(keep_n)]
}

selective_summary <- function(
  dat,
  rule = c("sd", "margin", "combined"),
  coverage = PRIMARY_COVERAGE,
  fn_cost = 1,
  fp_cost = 1
) {
  rule <- match.arg(rule)

  if (rule == "sd") {
    idx <- retained_indices(dat$p_sd, coverage, higher_is_better = FALSE)
  } else if (rule == "margin") {
    idx <- retained_indices(dat$margin, coverage, higher_is_better = TRUE)
  } else {
    idx <- retained_indices(
      dat$combined_error_prob,
      coverage,
      higher_is_better = FALSE
    )
  }

  d <- dat[idx, , drop = FALSE]

  pred <- as.integer(d$p_mean >= CLASS_THRESHOLD)
  y <- as.integer(d$Y)

  tp <- sum(pred == 1L & y == 1L)
  tn <- sum(pred == 0L & y == 0L)
  fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L)

  tibble(
    rule = rule,
    coverage_target = coverage,
    n_retained = nrow(d),
    coverage_realized = nrow(d) / nrow(dat),
    errors = fp + fn,
    selective_risk = (fp + fn) / nrow(d),
    weighted_risk = (fp_cost * fp + fn_cost * fn) / nrow(d),
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn,
    sensitivity = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
    specificity = ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
  )
}

risk_curve <- function(dat, rule, coverages = COVERAGE_GRID) {
  bind_rows(
    lapply(
      coverages,
      function(cov) selective_summary(dat, rule, cov)
    )
  )
}

naurc <- function(curve) {
  area <- trapz_base(curve$coverage_realized, curve$selective_risk)
  denom <- max(curve$coverage_realized) - min(curve$coverage_realized)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  area / denom
}

error_prediction_metrics <- function(y_error, p_error) {
  p_error <- clip_prob(p_error)
  tibble(
    brier = mean((p_error - y_error)^2),
    logloss = -mean(
      y_error * log(p_error) +
        (1 - y_error) * log(1 - p_error)
    ),
    auc = safe_auc(y_error, p_error)
  )
}

# ------------------------------------------------------------------------------
# 10. PAIRED PATIENT BOOTSTRAP FOR RULE CONTRASTS
# ------------------------------------------------------------------------------

bootstrap_rule_contrasts <- function(dat, B = N_BOOT, seed = SEED_BOOT) {
  obs <- c(
    sd_margin =
      selective_summary(dat, "sd")$selective_risk -
      selective_summary(dat, "margin")$selective_risk,
    combined_margin =
      selective_summary(dat, "combined")$selective_risk -
      selective_summary(dat, "margin")$selective_risk,
    sd_combined =
      selective_summary(dat, "sd")$selective_risk -
      selective_summary(dat, "combined")$selective_risk
  )

  set.seed(seed)

  boot <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(obs),
    dimnames = list(NULL, names(obs))
  )

  n <- nrow(dat)

  for (b in seq_len(B)) {
    ii <- sample.int(n, size = n, replace = TRUE)
    d <- dat[ii, , drop = FALSE]

    rsd <- selective_summary(d, "sd")$selective_risk
    rma <- selective_summary(d, "margin")$selective_risk
    rco <- selective_summary(d, "combined")$selective_risk

    boot[b, ] <- c(
      sd_margin = rsd - rma,
      combined_margin = rco - rma,
      sd_combined = rsd - rco
    )
  }

  out <- lapply(seq_along(obs), function(j) {
    v <- boot[, j]
    v <- v[is.finite(v)]

    ci <- stats::quantile(
      v,
      probs = c(0.025, 0.975),
      names = FALSE,
      type = 8
    )

    centered <- v - obs[j]

    # Finite-bootstrap corrected two-sided centered-bootstrap p-value.
    p_raw <- (
      1 + sum(abs(centered) >= abs(obs[j]))
    ) / (
      length(centered) + 1
    )

    tibble(
      comparison = names(obs)[j],
      difference = unname(obs[j]),
      ci_low = ci[1],
      ci_high = ci[2],
      raw_p = p_raw
    )
  }) |>
    bind_rows()

  out$holm_p <- p.adjust(out$raw_p, method = "holm")

  list(summary = out, draws = as.data.frame(boot))
}

# ------------------------------------------------------------------------------
# 11. COMPILE NESTED RPA RESULTS
# ------------------------------------------------------------------------------

compile_nested_results <- function(pred) {
  primary <- list()
  curves <- list()
  error_metrics <- list()
  boot_summaries <- list()
  boot_draws <- list()

  zp <- zc <- ze <- zb <- 0L

  for (model in MODEL_SET) {
    d <- pred[pred$model == model, , drop = FALSE]

    # Primary 80% risks.
    for (rule in c("sd", "margin", "combined")) {
      zp <- zp + 1L
      primary[[zp]] <- selective_summary(d, rule) |>
        mutate(model = model)
    }

    # Risk-coverage and nAURC.
    for (rule in c("sd", "margin", "combined")) {
      rc <- risk_curve(d, rule)
      rc$model <- model
      rc$nAURC <- naurc(rc)
      zc <- zc + 1L
      curves[[zc]] <- rc
    }

    # Error-risk selector diagnostic metrics.
    em1 <- error_prediction_metrics(d$error, d$margin_error_prob) |>
      mutate(model = model, selector = "margin_error_model")

    em2 <- error_prediction_metrics(d$error, d$combined_error_prob) |>
      mutate(model = model, selector = "combined_error_model")

    ze <- ze + 1L
    error_metrics[[ze]] <- bind_rows(em1, em2)

    # Paired patient bootstrap.
    bb <- bootstrap_rule_contrasts(
      d,
      B = N_BOOT,
      seed = SEED_BOOT + match(model, MODEL_SET) * 1000L
    )

    zb <- zb + 1L
    boot_summaries[[zb]] <- bb$summary |>
      mutate(model = model)

    boot_draws[[zb]] <- bb$draws |>
      mutate(
        model = model,
        bootstrap = seq_len(nrow(bb$draws))
      )
  }

  primary <- bind_rows(primary)
  curves <- bind_rows(curves)
  error_metrics <- bind_rows(error_metrics)
  boot_summary <- bind_rows(boot_summaries)
  boot_draws <- bind_rows(boot_draws)

  # Compact primary-results table.
  naurc_table <- curves |>
    group_by(model, rule) |>
    summarise(nAURC = first(nAURC), .groups = "drop") |>
    pivot_wider(
      names_from = rule,
      values_from = nAURC,
      names_prefix = "nAURC_"
    )

  compact <- boot_summary |>
    dplyr::filter(comparison %in% c("sd_margin", "combined_margin")) |>
    dplyr::select(
      model,
      comparison,
      difference,
      ci_low,
      ci_high,
      raw_p,
      holm_p
    ) |>
    left_join(naurc_table, by = "model")

  list(
    primary = primary,
    curves = curves,
    error_metrics = error_metrics,
    bootstrap_summary = boot_summary,
    bootstrap_draws = boot_draws,
    compact = compact
  )
}

# ------------------------------------------------------------------------------
# 12. BUILD FINAL RPA SELECTORS FOR LOCKED EXTERNAL APPLICATION
# ------------------------------------------------------------------------------

fit_final_rpa_selectors <- function(rpa) {
  out <- list()
  features_all <- list()
  manifest_all <- list()

  for (model in MODEL_SET) {
    message("Final full-RPA OOF feature construction: ", model)

    oof <- generate_repeated_oof_features(
      rpa,
      model = model,
      repeats = INNER_REPEATS,
      k = INNER_FOLDS,
      seed = SEED_NESTED + 50000000L,
      outer_id = 99L
    )

    selector_fit <- fit_selector_pair(oof$features)

    out[[model]] <- selector_fit
    features_all[[model]] <- oof$features |>
      mutate(model = model)

    manifest_all[[model]] <- oof$manifest
  }

  list(
    selectors = out,
    features = bind_rows(features_all),
    manifest = bind_rows(manifest_all)
  )
}

# ------------------------------------------------------------------------------
# 13. DÜSSELDORF RAW DATA, WBC AUDIT, AND LOCKED PREDICTIONS
# ------------------------------------------------------------------------------

first_existing_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (!length(hit)) {
    stop("None of these columns found: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

map_dusseldorf_outcome <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("[[:space:]_\\-]", "", z)

  out <- rep(NA_integer_, length(z))
  out[z %in% c("appendicitis", "yes", "1")] <- 1L
  out[z %in% c("noappendicitis", "no", "0")] <- 0L
  out
}

load_dusseldorf_raw <- function(path) {
  if (!file.exists(path)) stop("Düsseldorf file not found: ", path)

  dat <- data.table::fread(path, data.table = FALSE, check.names = FALSE)

  # Accept either the public source names or janitor-cleaned local names.
  col_temp <- first_existing_col(
    dat,
    c("BodyTemp", "body_temp", "body_temperature")
  )
  col_wbc <- first_existing_col(
    dat,
    c("WBCCount", "wbc_count")
  )
  col_neut <- first_existing_col(
    dat,
    c("NeutrophilPerc", "neutrophil_perc", "neutrophil_percentage")
  )
  col_y <- first_existing_col(
    dat,
    c("DiagnosisByCriteria", "diagnosis_by_criteria", "appendicitis_status")
  )

  tibble(
    Row_ID = seq_len(nrow(dat)),
    Y = map_dusseldorf_outcome(dat[[col_y]]),
    temperature_raw = suppressWarnings(as.numeric(dat[[col_temp]])),
    wbc_raw = suppressWarnings(as.numeric(dat[[col_wbc]])),
    neutrophil_raw = suppressWarnings(as.numeric(dat[[col_neut]]))
  )
}

wbc_audit <- function(dus) {
  w <- dus$wbc_raw

  tibble(
    n_total = length(w),
    n_nonmissing = sum(is.finite(w)),
    n_gt_100 = sum(w > 100, na.rm = TRUE),
    n_gt_1000 = sum(w > 1000, na.rm = TRUE),
    n_between_100_and_1000 = sum(w > 100 & w <= 1000, na.rm = TRUE),
    min_gt_100 = ifelse(any(w > 100, na.rm = TRUE), min(w[w > 100], na.rm = TRUE), NA_real_),
    max_wbc = max(w, na.rm = TRUE),
    legacy_gt100_equals_source_gt1000 =
      all(
        ifelse(w > 100, w / 1000, w) ==
          ifelse(w > 1000, w / 1000, w),
        na.rm = TRUE
      )
  )
}

harmonize_wbc <- function(w, rule = c("source_gt1000", "legacy_gt100", "none")) {
  rule <- match.arg(rule)
  z <- as.numeric(w)

  if (rule == "source_gt1000") {
    z[is.finite(z) & z > 1000] <- z[is.finite(z) & z > 1000] / 1000
  }

  if (rule == "legacy_gt100") {
    z[is.finite(z) & z > 100] <- z[is.finite(z) & z > 100] / 1000
  }

  z
}

make_dusseldorf_standardized <- function(dus, wbc_rule) {
  wbc <- harmonize_wbc(dus$wbc_raw, wbc_rule)

  out <- tibble(
    Y = dus$Y,
    Body_Temperature =
      (dus$temperature_raw - RPA_UPSTREAM$mean["Body_Temperature"]) /
      RPA_UPSTREAM$sd["Body_Temperature"],
    WBC_Count =
      (wbc - RPA_UPSTREAM$mean["WBC_Count"]) /
      RPA_UPSTREAM$sd["WBC_Count"],
    Neutrophil_Percentage =
      (dus$neutrophil_raw - RPA_UPSTREAM$mean["Neutrophil_Percentage"]) /
      RPA_UPSTREAM$sd["Neutrophil_Percentage"],
    Row_ID = dus$Row_ID
  )

  out
}

apply_final_selector <- function(selector_fit, features) {
  score <- apply_selector_pair(selector_fit, features)
  bind_cols(features, score)
}

fit_frozen_external_predictions <- function(rpa, dus_raw, final_selectors) {
  prep <- fit_preprocessor(rpa)
  xtr <- apply_preprocessor(rpa, prep)

  variants <- c("source_gt1000", "legacy_gt100", "none")

  dus_std_list <- lapply(
    variants,
    function(v) make_dusseldorf_standardized(dus_raw, v)
  )
  names(dus_std_list) <- variants

  X_list <- lapply(dus_std_list, apply_preprocessor, prep = prep)
  X_all <- do.call(rbind, X_list)

  n_ext <- nrow(dus_raw)

  all_out <- list()

  for (model in MODEL_SET) {
    message("Frozen external prediction + WBC sensitivity: ", model)

    fit <- fit_base_model(
      model,
      xtr,
      rpa$Y,
      X_all,
      seed = SEED_EXTERNAL + match(model, MODEL_SET)
    )

    sm_all <- summarize_draws(fit$draws)

    for (v_id in seq_along(variants)) {
      v <- variants[v_id]
      idx <- ((v_id - 1L) * n_ext + 1L):(v_id * n_ext)

      feat <- tibble(
        Row_ID = dus_raw$Row_ID,
        Y = dus_raw$Y,
        p_mean = sm_all$p_mean[idx],
        p_sd = sm_all$p_sd[idx]
      )

      feat$margin <- abs(feat$p_mean - CLASS_THRESHOLD)
      feat$ambiguity <- CLASS_THRESHOLD - feat$margin
      feat$error <- as.integer(
        as.integer(feat$p_mean >= CLASS_THRESHOLD) != feat$Y
      )

      scored <- apply_final_selector(
        final_selectors[[model]],
        feat
      )

      all_out[[length(all_out) + 1L]] <- scored |>
        mutate(
          model = model,
          wbc_rule = v
        )
    }
  }

  bind_rows(all_out)
}

# ------------------------------------------------------------------------------
# 14. EXTERNAL WBC SENSITIVITY METRICS
# ------------------------------------------------------------------------------

calibration_intercept_slope <- function(y, p) {
  p <- clip_prob(p)
  lp <- qlogis(p)

  fit <- tryCatch(
    suppressWarnings(glm(y ~ lp, family = binomial())),
    error = function(e) NULL
  )

  if (is.null(fit) || length(coef(fit)) < 2L) {
    return(c(intercept = NA_real_, slope = NA_real_))
  }

  c(
    intercept = unname(coef(fit)[1L]),
    slope = unname(coef(fit)[2L])
  )
}

base_external_metrics <- function(dat) {
  y <- dat$Y
  p <- clip_prob(dat$p_mean)
  pred <- as.integer(p >= CLASS_THRESHOLD)

  cal <- calibration_intercept_slope(y, p)

  tibble(
    n = nrow(dat),
    prevalence = mean(y),
    accuracy = mean(pred == y),
    sensitivity = sum(pred == 1 & y == 1) / sum(y == 1),
    specificity = sum(pred == 0 & y == 0) / sum(y == 0),
    auroc = safe_auc(y, p),
    brier = mean((p - y)^2),
    logloss = -mean(y * log(p) + (1 - y) * log(1 - p)),
    calibration_intercept = cal["intercept"],
    calibration_slope = cal["slope"],
    risk80_sd = selective_summary(dat, "sd")$selective_risk,
    risk80_margin = selective_summary(dat, "margin")$selective_risk,
    risk80_combined = selective_summary(dat, "combined")$selective_risk
  )
}

compile_wbc_sensitivity <- function(ext_pred) {
  ext_pred |>
    group_by(model, wbc_rule) |>
    group_modify(~ base_external_metrics(.x)) |>
    ungroup()
}

# ------------------------------------------------------------------------------
# 15. ASYMMETRIC-COST SENSITIVITY
# ------------------------------------------------------------------------------

compile_cost_sensitivity <- function(pred, cohort_label) {
  rows <- list()
  z <- 0L

  for (model in MODEL_SET) {
    d <- pred[pred$model == model, , drop = FALSE]

    for (ratio in FN_COST_RATIOS) {
      for (rule in c("sd", "margin", "combined")) {
        z <- z + 1L

        s <- selective_summary(
          d,
          rule = rule,
          coverage = PRIMARY_COVERAGE,
          fn_cost = ratio,
          fp_cost = 1
        )

        rows[[z]] <- s |>
          mutate(
            cohort = cohort_label,
            model = model,
            fn_to_fp_cost_ratio = paste0(ratio, ":1"),
            .before = 1
          )
      }
    }
  }

  bind_rows(rows)
}

# ------------------------------------------------------------------------------
# 16. MAIN RUN
# ------------------------------------------------------------------------------

message("Loading data...")
rpa <- load_rpa(RPA_FILE)
dus_raw <- load_dusseldorf_raw(DUS_FILE)

# ----- WBC harmonization audit -------------------------------------------------

audit_wbc <- wbc_audit(dus_raw)
save_csv(audit_wbc, "01_wbc_unit_audit.csv")

affected <- tibble(
  Row_ID = dus_raw$Row_ID,
  WBC_raw = dus_raw$wbc_raw,
  WBC_legacy_gt100 = harmonize_wbc(dus_raw$wbc_raw, "legacy_gt100"),
  WBC_source_gt1000 = harmonize_wbc(dus_raw$wbc_raw, "source_gt1000"),
  WBC_no_conversion = harmonize_wbc(dus_raw$wbc_raw, "none")
) |>
  dplyr::filter(
    is.finite(WBC_raw) &
      (
        WBC_raw > 100 |
        WBC_legacy_gt100 != WBC_source_gt1000
      )
  )

save_csv(affected, "02_wbc_affected_records.csv")

message("WBC audit:")
print(audit_wbc)

# ----- Fully nested RPA selector evaluation -------------------------------------

message("Starting fully nested RPA selector evaluation...")
nested <- run_nested_rpa(rpa)

save_csv(
  nested$predictions,
  "03_nested_rpa_patient_predictions.csv"
)
save_csv(
  nested$dependency,
  "04_nested_dependency_audit.csv"
)
save_csv(
  nested$coefficients,
  "05_nested_selector_coefficients.csv"
)

nested_results <- compile_nested_results(nested$predictions)

save_csv(
  nested_results$primary,
  "06_nested_primary_risk80.csv"
)
save_csv(
  nested_results$curves,
  "07_nested_risk_coverage.csv"
)
save_csv(
  nested_results$error_metrics,
  "08_nested_error_prediction_metrics.csv"
)
save_csv(
  nested_results$bootstrap_summary,
  "09_nested_bootstrap_contrasts.csv"
)
save_csv(
  nested_results$bootstrap_draws,
  "10_nested_bootstrap_draws.csv"
)
save_csv(
  nested_results$compact,
  "11_primary_incremental_value_table.csv"
)

# ----- Fit final RPA selector using RPA only; no Düsseldorf outcomes -----------

message("Building final RPA selector objects for locked external application...")
final_rpa <- fit_final_rpa_selectors(rpa)

save_csv(
  final_rpa$features,
  "12_final_rpa_selector_training_features.csv"
)
save_csv(
  final_rpa$manifest,
  "13_final_rpa_selector_oof_manifest.csv"
)
save_rds(
  final_rpa$selectors,
  "14_final_rpa_selectors.rds"
)

# ----- WBC sensitivity / locked Düsseldorf prediction --------------------------

message("Running locked Düsseldorf WBC sensitivity...")
ext_pred <- fit_frozen_external_predictions(
  rpa,
  dus_raw,
  final_rpa$selectors
)

save_csv(
  ext_pred,
  "15_dusseldorf_predictions_wbc_sensitivity.csv"
)

wbc_metrics <- compile_wbc_sensitivity(ext_pred)
save_csv(
  wbc_metrics,
  "16_dusseldorf_wbc_sensitivity_metrics.csv"
)

# ----- Asymmetric-cost sensitivity ---------------------------------------------

rpa_cost <- compile_cost_sensitivity(
  nested$predictions,
  "RPA_nested"
)

ext_primary <- ext_pred |>
  dplyr::filter(wbc_rule == "source_gt1000")

dus_cost <- compile_cost_sensitivity(
  ext_primary,
  "Dusseldorf_locked"
)

cost_results <- bind_rows(rpa_cost, dus_cost)

save_csv(
  cost_results,
  "17_asymmetric_fn_cost_sensitivity.csv"
)

# ----- False-negative capture at 80% retained coverage -------------------------

fn_capture <- list()
z <- 0L

for (cohort_name in c("RPA_nested", "Dusseldorf_locked")) {
  dat0 <- if (cohort_name == "RPA_nested") {
    nested$predictions
  } else {
    ext_primary
  }

  for (model in MODEL_SET) {
    d <- dat0[dat0$model == model, , drop = FALSE]

    pred_full <- as.integer(d$p_mean >= CLASS_THRESHOLD)
    total_fn <- sum(pred_full == 0L & d$Y == 1L)

    for (rule in c("sd", "margin", "combined")) {
      if (rule == "sd") {
        kept <- retained_indices(d$p_sd, PRIMARY_COVERAGE, FALSE)
      } else if (rule == "margin") {
        kept <- retained_indices(d$margin, PRIMARY_COVERAGE, TRUE)
      } else {
        kept <- retained_indices(d$combined_error_prob, PRIMARY_COVERAGE, FALSE)
      }

      deferred <- setdiff(seq_len(nrow(d)), kept)

      fn_deferred <- sum(
        pred_full[deferred] == 0L &
          d$Y[deferred] == 1L
      )

      z <- z + 1L
      fn_capture[[z]] <- tibble(
        cohort = cohort_name,
        model = model,
        rule = rule,
        total_false_negatives = total_fn,
        deferred_false_negatives = fn_deferred,
        false_negative_capture =
          ifelse(total_fn > 0, fn_deferred / total_fn, NA_real_)
      )
    }
  }
}

fn_capture <- bind_rows(fn_capture)
save_csv(
  fn_capture,
  "18_false_negative_capture.csv"
)

# ----- Session information ------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "19_sessionInfo.txt")
)

# ----- Console summary ----------------------------------------------------------

cat("\n============================================================\n")
cat("B-RISA nested incremental-selector analysis completed.\n")
cat("============================================================\n\n")

cat("Nested dependency audit passed: ",
    all(nested$dependency$leakage_pass), "\n", sep = "")

cat("\nPrimary nested incremental-value table:\n")
print(nested_results$compact)

cat("\nWBC audit:\n")
print(audit_wbc)

cat("\nDüsseldorf WBC sensitivity metrics:\n")
print(wbc_metrics)

cat("\nIllustrative asymmetric-cost sensitivity:\n")
print(cost_results)

cat("\nFalse-negative capture:\n")
print(fn_capture)

cat("\nOutputs saved to:\n", OUT_DIR, "\n", sep = "")

# ==============================================================================
# END
# ==============================================================================
