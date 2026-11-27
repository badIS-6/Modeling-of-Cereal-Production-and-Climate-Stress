options(stringsAsFactors = FALSE, scipen = 999)

SCRIPT_ARGS <- commandArgs(trailingOnly = FALSE)
SCRIPT_FILE_ARG <- grep("^--file=", SCRIPT_ARGS, value = TRUE)

SCRIPT_DIR <- if (length(SCRIPT_FILE_ARG) > 0L) {
  dirname(normalizePath(sub("^--file=", "", SCRIPT_FILE_ARG[1L]), winslash = "/", mustWork = TRUE))
} else {
  active_script <- if (requireNamespace("rstudioapi", quietly = TRUE)) {
    tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
  } else {
    ""
  }
  if (nzchar(active_script) && file.exists(active_script)) {
    dirname(normalizePath(active_script, winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
}

PROJECT_DIR <- "/home/badisz/Downloads/Modeling"
OUTPUT_DIR <- "/home/badisz/Downloads/inference/inf"
DATA_FILE <- file.path(PROJECT_DIR, "model_table.csv")

MODEL_OUTPUT_DIRS <- c(
  GAM = file.path(PROJECT_DIR, "GAM_Optimized_Output"),
  QR = file.path(PROJECT_DIR, "Quantile_Regression_Optimized_Output"),
  MIXED = file.path(PROJECT_DIR, "Mixed_Effects_Optimized_Output_nlme")
)

MODEL_FILES <- c(
  GAM = "gam_forecasting_selected_model.rds",
  QR = "quantile_regression_median_model.rds",
  MIXED = "mixed_effects_optimized_durum_wheat_model.rds"
)

BUNDLE_FILES <- c(
  GAM = "gam_forecasting_complete_object.rds",
  QR = "quantile_regression_forecasting_bundle.rds",
  MIXED = "mixed_effects_forecast_bundle.rds"
)

TRAIN_START_YEAR <- 2000L
TRAIN_END_YEAR <- 2012L
EXPECTED_TRAIN_N <- 192L
EXPECTED_GOVERNORATES <- 15L
BOOTSTRAP_REPS <- 2000L
BOOTSTRAP_SEED <- 20261127L
MIN_DIRECTIONAL_STABILITY <- 0.90
MIN_LOO_DIRECTIONAL_STABILITY <- 0.80
MIN_MODELS_FOR_COMPONENT_SUPPORT <- 2L
MIN_SCENARIO_COVERAGE_PCT <- 70
MIN_GOVERNORATES <- 10L
TEMPERATURE_INCREASE_C <- 1
MOISTURE_DECREASE_SD <- 0.5
SUPPORT_QUANTILES <- c(0.02, 0.98)
REGIONAL_CI_LEVEL <- 0.95

if (!dir.exists(PROJECT_DIR)) stop("Project directory does not exist: ", PROJECT_DIR)
if (!file.exists(DATA_FILE)) stop("Dataset does not exist: ", DATA_FILE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUTPUT_DIR)) stop("Could not create inference output directory: ", OUTPUT_DIR)

required_packages <- c("mgcv", "quantreg", "nlme")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running inference: ",
    paste(missing_packages, collapse = ", "),
    ". Use install.packages(c(", paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}

model_paths <- stats::setNames(
  file.path(unname(MODEL_OUTPUT_DIRS), unname(MODEL_FILES)),
  names(MODEL_FILES)
)
bundle_paths <- stats::setNames(
  file.path(unname(MODEL_OUTPUT_DIRS), unname(BUNDLE_FILES)),
  names(BUNDLE_FILES)
)

missing_model_dirs <- unique(unname(MODEL_OUTPUT_DIRS[!dir.exists(MODEL_OUTPUT_DIRS)]))
if (length(missing_model_dirs) > 0L) {
  stop(
    "One or more original forecasting output directories are missing:\n",
    paste(missing_model_dirs, collapse = "\n"),
    "\nCheck PROJECT_DIR and the folder names at the top of the script."
  )
}

missing_files <- c(model_paths[!file.exists(model_paths)], bundle_paths[!file.exists(bundle_paths)])
if (length(missing_files) > 0L) {
  cat("\nMODEL AND BUNDLE PATHS EXPECTED BY THIS SCRIPT:\n")
  print(c(model_paths, bundle_paths))
  cat("\nRDS FILES FOUND IN THE THREE FORECASTING OUTPUT FOLDERS:\n")
  for (folder in unique(unname(MODEL_OUTPUT_DIRS))) {
    cat("\n", folder, "\n", sep = "")
    print(list.files(folder, pattern = "\\.[rR][dD][sS]$", full.names = FALSE))
  }
  stop(
    "Required model or bundle files are missing. The inference script now reads each model from its original forecasting output folder; do not copy the RDS files into the inference/inf folder."
  )
}

model_paths <- vapply(model_paths, normalizePath, character(1), winslash = "/", mustWork = TRUE)
bundle_paths <- vapply(bundle_paths, normalizePath, character(1), winslash = "/", mustWork = TRUE)

read_rds_checked <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) stop("Could not read RDS file ", path, ": ", conditionMessage(e), call. = FALSE)
  )
}

models_raw <- lapply(model_paths, read_rds_checked)
bundles <- lapply(bundle_paths, read_rds_checked)

find_model_in_object <- function(x, predicate, depth = 0L, max_depth = 12L) {
  if (isTRUE(predicate(x))) return(x)
  if (depth >= max_depth || !is.list(x)) return(NULL)
  preferred_names <- c("model", "selected_model", "fit", "final_model", "median_model", "models")
  nms <- names(x)
  if (is.null(nms)) nms <- rep("", length(x))
  ord <- c(which(nms %in% preferred_names), which(!(nms %in% preferred_names)))
  for (i in ord) {
    candidate <- find_model_in_object(x[[i]], predicate, depth + 1L, max_depth)
    if (!is.null(candidate)) return(candidate)
  }
  NULL
}

is_gam_fit <- function(x) inherits(x, "gam") || inherits(x, "bam")
is_qr_fit <- function(x) inherits(x, "rq")
is_mixed_fit <- function(x) inherits(x, "lme") || inherits(x, "lmerMod")

is_median_qr_fit <- function(x) {
  is_qr_fit(x) && length(x$tau) == 1L && is.finite(as.numeric(x$tau[1L])) &&
    abs(as.numeric(x$tau[1L]) - 0.5) <= 1e-8
}

qr_bundle_median <- if (is.list(bundles$QR$models)) bundles$QR$models[["0.50"]] else NULL
qr_standalone <- if (is_median_qr_fit(models_raw$QR)) models_raw$QR else NULL
qr_fallback <- if (is.null(qr_standalone) && is_median_qr_fit(qr_bundle_median)) {
  qr_bundle_median
} else if (is.null(qr_standalone)) {
  find_model_in_object(bundles$QR, is_median_qr_fit)
} else {
  NULL
}

models <- list(
  GAM = if (is_gam_fit(models_raw$GAM)) models_raw$GAM else find_model_in_object(bundles$GAM, is_gam_fit),
  QR = if (!is.null(qr_standalone)) qr_standalone else qr_fallback,
  MIXED = if (is_mixed_fit(models_raw$MIXED)) models_raw$MIXED else find_model_in_object(bundles$MIXED, is_mixed_fit)
)

object_description <- function(x) {
  paste0(
    "class={", paste(class(x), collapse = ", "), "}; ",
    "top-level names={", paste(names(x), collapse = ", "), "}"
  )
}

for (model_name in names(models)) {
  if (is.null(models[[model_name]])) {
    stop(
      model_name, " fitted model could not be recovered.\n",
      "Standalone path: ", model_paths[[model_name]], "\n",
      "Standalone object: ", object_description(models_raw[[model_name]]), "\n",
      "Bundle path: ", bundle_paths[[model_name]], "\n",
      "Bundle object: ", object_description(bundles[[model_name]]), "\n",
      "This usually means the forecasting output folder contains an empty/NULL RDS or a different run's file. Rerun the corresponding forecasting script and check that it completes its final saveRDS() calls."
    )
  }
}

if (!is_gam_fit(models$GAM)) stop("Extracted GAM object is not an mgcv GAM/bam. Found: ", object_description(models$GAM))
if (!is_qr_fit(models$QR)) stop("Extracted QR object is not a quantreg::rq object. Found: ", object_description(models$QR))
if (!is_mixed_fit(models$MIXED)) {
  stop("Extracted Mixed Effects object is not nlme::lme or lme4::lmer. Found: ", object_description(models$MIXED))
}
if (length(models$QR$tau) == 0L || abs(as.numeric(models$QR$tau[1L]) - 0.5) > 1e-8) {
  stop("The QR model must be the median quantile model (tau = 0.50).")
}

cat("\n========== MODEL LOADING AUDIT ==========\n")
for (model_name in names(models)) {
  cat(
    model_name, ": ", paste(class(models[[model_name]]), collapse = "/"),
    "\n  model: ", model_paths[[model_name]],
    "\n  bundle: ", bundle_paths[[model_name]], "\n",
    sep = ""
  )
}
cat("Inference outputs will be written to: ", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), "\n", sep = "")

find_data_file <- function() {
  candidates <- c(
    if (nzchar(DATA_FILE)) DATA_FILE,
    file.path(PROJECT_DIR, "model_table.csv"),
    file.path(SCRIPT_DIR, "model_table.csv"),
    file.path(dirname(SCRIPT_DIR), "model_table.csv")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) NA_character_ else normalizePath(found[1L], winslash = "/", mustWork = TRUE)
}
data_path <- find_data_file()
if (is.na(data_path)) {
  stop("model_table.csv was not found. Put it in the inference folder or set DATA_FILE near the top of the script.")
}
raw_data <- utils::read.csv(
  data_path,
  fileEncoding = "UTF-8",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_columns <- c(
  "governorate", "crop", "harvest_year", "production_1000_qt",
  "rain_jan_mar_mm", "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)
missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0L) stop("model_table.csv is missing: ", paste(missing_columns, collapse = ", "))

raw_data$governorate <- trimws(as.character(raw_data$governorate))
raw_data$crop <- trimws(as.character(raw_data$crop))
raw_data$harvest_year <- suppressWarnings(as.integer(raw_data$harvest_year))
raw_data <- raw_data[raw_data$crop == "durum_wheat" &
                       raw_data$harvest_year >= TRAIN_START_YEAR &
                       raw_data$harvest_year <= TRAIN_END_YEAR, , drop = FALSE]

north_governorates <- c(
  "Béja", "Jendouba", "Le Kef", "Siliana", "Zaghouan",
  "Manouba", "Nabeul", "Tunis"
)
center_south_governorates <- c(
  "Kairouan", "Kasserine", "Sidi Bouzid", "Sousse",
  "Gafsa", "Gabès", "Medenine"
)
expected_governorates <- c(
  "Béja", "Gabès", "Gafsa", "Jendouba", "Kairouan",
  "Kasserine", "Le Kef", "Manouba", "Medenine", "Nabeul",
  "Sidi Bouzid", "Siliana", "Sousse", "Tunis", "Zaghouan"
)

normalize_governorate <- function(x) {
  z <- trimws(as.character(x))
  z <- gsub("^Kef$|^El Kef$|^Al Kef$", "Le Kef", z, ignore.case = TRUE)
  z <- gsub("^Beja$", "Béja", z, ignore.case = TRUE)
  z <- gsub("^Gabes$", "Gabès", z, ignore.case = TRUE)
  z
}
raw_data$governorate <- normalize_governorate(raw_data$governorate)
numeric_columns <- c("harvest_year", "production_1000_qt", "rain_jan_mar_mm", "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi")
for (column_name in numeric_columns) {
  raw_data[[column_name]] <- suppressWarnings(as.numeric(raw_data[[column_name]]))
}
if (any(!vapply(raw_data[, numeric_columns, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
  stop("The 2000-2012 durum-wheat rows contain missing or non-finite required numeric values.")
}

unmapped_governorates <- setdiff(unique(raw_data$governorate), expected_governorates)
if (length(unmapped_governorates) > 0L) {
  stop("Governorates are not mapped to the 15-governorate model design: ", paste(unmapped_governorates, collapse = ", "))
}
if (length(unique(raw_data$governorate)) != EXPECTED_GOVERNORATES) {
  stop("Expected ", EXPECTED_GOVERNORATES, " governorates in the model training period; found ", length(unique(raw_data$governorate)), ".")
}
if (anyDuplicated(raw_data[, c("governorate", "harvest_year")]) > 0L) {
  stop("Duplicate governorate-year rows were found in model_table.csv for 2000-2012.")
}

raw_data$region <- ifelse(raw_data$governorate %in% north_governorates, "North", "Center_South")
raw_data$region <- factor(raw_data$region, levels = c("North", "Center_South"))
raw_data <- raw_data[order(raw_data$governorate, raw_data$harvest_year), , drop = FALSE]
rownames(raw_data) <- NULL

if (nrow(raw_data) != EXPECTED_TRAIN_N) {
  stop(
    "The reconstructed 2000-2012 training sample has ", nrow(raw_data),
    " rows; the forecasting scripts specify ", EXPECTED_TRAIN_N,
    ". Check that model_table.csv matches the data used to fit the three saved models."
  )
}

get_bundle_metadata <- function(model_name, bundle, fit) {
  if (model_name == "GAM") {
    transformation <- bundle$response_transformation
    scaling <- bundle$scaling
    smearing <- if (identical(transformation, "log1p")) {
      mean(exp(as.numeric(stats::residuals(fit, type = "response"))), na.rm = TRUE)
    } else 1
    gov_levels <- fit$xlevels$governorate %||% as.character(scaling$sm$governorate)
    sm_mean_name <- "sm_mean"
    sm_sd_name <- "sm_sd"
    y_mean_name <- "production_mean"
    y_sd_name <- "production_sd"
    sm_stats <- scaling$sm
    y_stats <- scaling$y
  } else if (model_name == "QR") {
    transformation <- bundle$response_transformation
    scaling <- bundle$preprocessing
    smearing <- 1
    gov_levels <- scaling$governorate_levels %||% as.character(scaling$sm$governorate)
    sm_mean_name <- "sm_train_mean"
    sm_sd_name <- "sm_train_sd"
    y_mean_name <- "production_train_mean"
    y_sd_name <- "production_train_sd"
    sm_stats <- scaling$sm
    y_stats <- scaling$production
  } else {
    transformation <- bundle$response_transformation
    scaling <- bundle$scaling_parameters
    smearing <- bundle$smearing_factor
    if (identical(transformation, "log1p") &&
        (is.null(smearing) || length(smearing) != 1L || !is.finite(smearing) || smearing <= 0)) {
      mixed_residuals <- as.numeric(stats::residuals(fit, level = 1, type = "response"))
      smearing <- mean(exp(mixed_residuals), na.rm = TRUE)
    }
    if (is.null(smearing)) smearing <- 1
    gov_levels <- bundle$governorate_levels %||% as.character(scaling$sm$governorate)
    sm_mean_name <- "sm_mean"
    sm_sd_name <- "sm_sd"
    y_mean_name <- "production_mean"
    y_sd_name <- "production_sd"
    sm_stats <- scaling$sm
    y_stats <- scaling$y
  }

  if (length(transformation) != 1L || is.na(transformation) || !transformation %in% c("raw", "log1p", "governorate_z")) {
    stop(model_name, ": the saved bundle does not contain a recognized response_transformation (raw, log1p, governorate_z).")
  }
  if (!is.data.frame(sm_stats) || !all(c("governorate", sm_mean_name, sm_sd_name) %in% names(sm_stats))) {
    stop(model_name, ": soil-moisture scaling parameters could not be recovered from the original model bundle.")
  }
  if (transformation == "governorate_z" &&
      (!is.data.frame(y_stats) || !all(c("governorate", y_mean_name, y_sd_name) %in% names(y_stats)))) {
    stop(model_name, ": governorate-specific production scaling parameters could not be recovered from the original model bundle.")
  }
  if (transformation == "log1p" && (!is.finite(smearing) || smearing <= 0)) {
    stop(model_name, ": a valid log1p retransformation factor could not be recovered.")
  }

  list(
    transformation = as.character(transformation),
    scaling = scaling,
    sm_stats = sm_stats,
    y_stats = y_stats,
    sm_mean_name = sm_mean_name,
    sm_sd_name = sm_sd_name,
    y_mean_name = y_mean_name,
    y_sd_name = y_sd_name,
    smearing_factor = as.numeric(smearing),
    governorate_levels = as.character(gov_levels)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

metadata <- lapply(names(models), function(nm) get_bundle_metadata(nm, bundles[[nm]], models[[nm]]))
names(metadata) <- names(models)

prepare_model_data <- function(model_name, fit, meta) {
  d <- raw_data
  sm_index <- match(as.character(d$governorate), as.character(meta$sm_stats$governorate))
  if (anyNA(sm_index)) stop(model_name, ": missing governorate soil-moisture scaling values in the model bundle.")
  sm_mean <- as.numeric(meta$sm_stats[[meta$sm_mean_name]][sm_index])
  sm_sd <- as.numeric(meta$sm_stats[[meta$sm_sd_name]][sm_index])
  if (any(!is.finite(sm_mean)) || any(!is.finite(sm_sd)) || any(sm_sd <= 0)) {
    stop(model_name, ": invalid governorate-specific soil-moisture scaling parameters in the saved bundle.")
  }
  d$sm_z <- (as.numeric(d$sm_oct_nov) - sm_mean) / sm_sd
  d$sm_deficit <- pmax(-1 - d$sm_z, 0)
  d$sm_wet <- pmax(d$sm_z - 1, 0)
  d$year_c <- (d$harvest_year - 2006) / 10
  center_south <- as.numeric(as.character(d$region) == "Center_South")
  d$temp_cs <- d$temp_dec_c * center_south
  d$sm_z_cs <- d$sm_z * center_south
  d$sm_deficit_cs <- d$sm_deficit * center_south
  d$sm_wet_cs <- d$sm_wet * center_south
  d$ndvi_cs <- d$ndvi_mean * center_south
  d$vhi_cs <- d$vhi * center_south
  d$governorate <- factor(as.character(d$governorate), levels = meta$governorate_levels)
  d$region <- factor(as.character(d$region), levels = c("North", "Center_South"))
  if (model_name == "MIXED" && nlevels(d$region) == 2L) contrasts(d$region) <- stats::contr.sum(2L)

  f <- stats::formula(fit)
  rhs_vars <- unique(all.vars(f[[3L]]))
  missing_predictors <- setdiff(rhs_vars, names(d))
  if (length(missing_predictors) > 0L) {
    stop(model_name, ": reconstructed training data are missing formula variables: ", paste(missing_predictors, collapse = ", "))
  }
  keep <- stats::complete.cases(d[, rhs_vars, drop = FALSE]) &
    !is.na(d$governorate) & !is.na(d$region) & is.finite(d$production_1000_qt)
  d <- d[keep, , drop = FALSE]
  if (nrow(d) < 5L) stop(model_name, ": fewer than five complete fitted observations were reconstructed.")

  fit_n <- tryCatch({
    fitted_frame <- stats::model.frame(fit)
    if (is.data.frame(fitted_frame)) nrow(fitted_frame) else NA_integer_
  }, error = function(e) NA_integer_)
  if (length(fit_n) != 1L || !is.finite(fit_n)) {
    fit_n <- tryCatch(length(stats::residuals(fit)), error = function(e) NA_integer_)
  }
  if (length(fit_n) == 1L && is.finite(fit_n) && nrow(d) != as.integer(fit_n)) {
    stop(
      model_name, ": the reconstructed training sample has ", nrow(d),
      " observations but the saved fitted model has ", fit_n,
      ". The source data and model object may not correspond."
    )
  }
  if (length(unique(d$governorate)) < MIN_GOVERNORATES) {
    stop(model_name, ": fewer than ", MIN_GOVERNORATES, " governorates remain in the inference sample.")
  }
  d
}

model_data <- lapply(names(models), function(nm) prepare_model_data(nm, models[[nm]], metadata[[nm]]))
names(model_data) <- names(models)

predict_model_production <- function(model_name, fit, newdata, meta) {
  prediction <- tryCatch({
    if (inherits(fit, "gam")) {
      stats::predict(fit, newdata = newdata, type = "response")
    } else if (inherits(fit, "rq")) {
      stats::predict(fit, newdata = newdata)
    } else if (inherits(fit, "lme")) {
      level <- bundles[[model_name]]$prediction_level %||% 1L
      stats::predict(fit, newdata = newdata, level = as.integer(level))
    } else if (inherits(fit, "lmerMod")) {
      stats::predict(fit, newdata = newdata, re.form = NULL, allow.new.levels = TRUE)
    } else {
      stop("Unsupported fitted model class.")
    }
  }, error = function(e) {
    stop("Prediction failed for ", model_name, ": ", conditionMessage(e))
  })

  if (is.matrix(prediction) || is.data.frame(prediction)) {
    if (ncol(prediction) != 1L) stop(model_name, ": prediction returned more than one column.")
    prediction <- prediction[, 1L]
  }
  eta <- as.numeric(prediction)
  if (length(eta) != nrow(newdata) || any(!is.finite(eta))) {
    stop(model_name, ": prediction produced missing or non-finite values.")
  }

  if (meta$transformation == "raw") {
    output <- eta
  } else if (meta$transformation == "log1p") {
    if (model_name == "QR") {
      output <- expm1(eta)
    } else {
      output <- exp(eta) * meta$smearing_factor - 1
    }
  } else if (meta$transformation == "governorate_z") {
    y_index <- match(as.character(newdata$governorate), as.character(meta$y_stats$governorate))
    if (anyNA(y_index)) stop(model_name, ": missing production retransformation values for one or more governorates.")
    y_mean <- as.numeric(meta$y_stats[[meta$y_mean_name]][y_index])
    y_sd <- as.numeric(meta$y_stats[[meta$y_sd_name]][y_index])
    if (any(!is.finite(y_mean)) || any(!is.finite(y_sd)) || any(y_sd <= 0)) {
      stop(model_name, ": invalid governorate-specific production scaling values.")
    }
    output <- eta * y_sd + y_mean
  } else {
    stop(model_name, ": unsupported response transformation.")
  }

  if (any(!is.finite(output))) stop(model_name, ": retransformation produced non-finite production predictions.")
  as.numeric(output)
}

rebuild_predictor_features <- function(d, model_name) {
  center_south <- as.numeric(as.character(d$region) == "Center_South")
  d$temp_cs <- d$temp_dec_c * center_south
  d$sm_deficit <- pmax(-1 - d$sm_z, 0)
  d$sm_wet <- pmax(d$sm_z - 1, 0)
  d$sm_z_cs <- d$sm_z * center_south
  d$sm_deficit_cs <- d$sm_deficit * center_south
  d$sm_wet_cs <- d$sm_wet * center_south
  d$ndvi_cs <- d$ndvi_mean * center_south
  d$vhi_cs <- d$vhi * center_south
  if (model_name == "MIXED" && is.factor(d$region) && nlevels(d$region) == 2L) contrasts(d$region) <- stats::contr.sum(2L)
  d
}

make_scenario <- function(d, model_name, stressor) {
  if (stressor == "temperature") {
    variable <- "temp_dec_c"
    shift <- TEMPERATURE_INCREASE_C
  } else {
    variable <- "sm_z"
    shift <- -MOISTURE_DECREASE_SD
  }
  x <- as.numeric(d[[variable]])
  support <- stats::quantile(x, probs = SUPPORT_QUANTILES, na.rm = TRUE, names = FALSE, type = 7)
  x_new <- x + shift
  eligible <- is.finite(x) & x >= support[1L] & x <= support[2L] &
    x_new >= support[1L] & x_new <= support[2L]
  if (sum(eligible) < 5L) {
    stop(model_name, ": fewer than five observations remain inside common support for ", stressor, ".")
  }
  baseline <- d[eligible, , drop = FALSE]
  stressed <- baseline
  if (stressor == "temperature") {
    stressed$temp_dec_c <- stressed$temp_dec_c + TEMPERATURE_INCREASE_C
  } else {
    stressed$sm_z <- stressed$sm_z - MOISTURE_DECREASE_SD
  }
  stressed <- rebuild_predictor_features(stressed, model_name)
  list(
    baseline = baseline,
    stressed = stressed,
    variable = variable,
    shift = shift,
    n_eligible = nrow(baseline),
    n_original = nrow(d),
    support_low = unname(support[1L]),
    support_high = unname(support[2L])
  )
}

cluster_effect_table <- function(data, delta, baseline, stressed) {
  row_data <- data.frame(
    governorate = as.character(data$governorate),
    region = as.character(data$region),
    delta = as.numeric(delta),
    baseline = as.numeric(baseline),
    stressed = as.numeric(stressed),
    stringsAsFactors = FALSE
  )
  gov_names <- sort(unique(row_data$governorate))
  out <- lapply(gov_names, function(g) {
    z <- row_data[row_data$governorate == g, , drop = FALSE]
    regions <- unique(z$region)
    if (length(regions) != 1L) stop("A governorate maps to more than one region.")
    data.frame(
      governorate = g,
      region = regions[1L],
      n_years = nrow(z),
      baseline_mean = mean(z$baseline),
      stressed_mean = mean(z$stressed),
      mean_change = mean(z$delta),
      mean_change_pct = if (mean(z$baseline) != 0) 100 * mean(z$delta) / mean(z$baseline) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

bootstrap_cluster_mean <- function(values, reps, seed) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(rep(NA_real_, reps))
  set.seed(seed)
  replicate(reps, mean(sample(values, size = length(values), replace = TRUE)))
}

summarize_effect <- function(values, expected_direction = c("negative", "either")) {
  expected_direction <- match.arg(expected_direction)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(list(ci_low = NA_real_, ci_high = NA_real_, stability = NA_real_, positive_pct = NA_real_, negative_pct = NA_real_))
  }
  positive_pct <- mean(values > 0)
  negative_pct <- mean(values < 0)
  stability <- if (expected_direction == "negative") negative_pct else max(positive_pct, negative_pct)
  list(
    ci_low = unname(stats::quantile(values, 0.025, na.rm = TRUE)),
    ci_high = unname(stats::quantile(values, 0.975, na.rm = TRUE)),
    stability = stability,
    positive_pct = positive_pct,
    negative_pct = negative_pct
  )
}

bootstrap_regional_contrast <- function(gov_table, reps, seed) {
  north <- gov_table$mean_change[gov_table$region == "North"]
  center_south <- gov_table$mean_change[gov_table$region == "Center_South"]
  north <- north[is.finite(north)]
  center_south <- center_south[is.finite(center_south)]
  if (length(north) < 3L || length(center_south) < 3L) {
    return(list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, stability = NA_real_,
                north_n = length(north), center_south_n = length(center_south)))
  }
  estimate <- mean(north) - mean(center_south)
  set.seed(seed)
  boot <- replicate(reps, {
    mean(sample(north, length(north), replace = TRUE)) -
      mean(sample(center_south, length(center_south), replace = TRUE))
  })
  directional_stability <- max(mean(boot > 0), mean(boot < 0))
  list(
    estimate = estimate,
    ci_low = unname(stats::quantile(boot, 0.025, na.rm = TRUE)),
    ci_high = unname(stats::quantile(boot, 0.975, na.rm = TRUE)),
    stability = directional_stability,
    north_n = length(north),
    center_south_n = length(center_south)
  )
}

overall_rows <- list()
regional_rows <- list()
heterogeneity_rows <- list()
loo_rows <- list()
diagnostic_rows <- list()
model_sample_rows <- list()

for (model_index in seq_along(models)) {
  model_name <- names(models)[model_index]
  fit <- models[[model_name]]
  meta <- metadata[[model_name]]
  d_full <- model_data[[model_name]]
  formula <- stats::formula(fit)
  response_expression <- paste(deparse(formula[[2L]]), collapse = "")

  model_sample_rows[[model_name]] <- data.frame(
    model = model_name,
    model_file = basename(model_paths[[model_name]]),
    bundle_file = basename(bundle_paths[[model_name]]),
    response_expression = response_expression,
    response_transformation = meta$transformation,
    rows_used = nrow(d_full),
    governorates_used = length(unique(d_full$governorate)),
    first_year = min(d_full$harvest_year),
    last_year = max(d_full$harvest_year),
    production_scale_verified = TRUE,
    response_scale_after_retransformation = "production_1000_qt",
    stringsAsFactors = FALSE
  )

  for (stressor in c("temperature", "moisture")) {
    stressor_label <- if (stressor == "temperature") "Warmer December (+1 C)" else "Worse autumn soil moisture (-0.5 sm_z)"
    scenario <- make_scenario(d_full, model_name, stressor)
    baseline_prediction <- predict_model_production(model_name, fit, scenario$baseline, meta)
    stressed_prediction <- predict_model_production(model_name, fit, scenario$stressed, meta)
    delta <- stressed_prediction - baseline_prediction
    gov_table <- cluster_effect_table(scenario$baseline, delta, baseline_prediction, stressed_prediction)

    if (nrow(gov_table) < MIN_GOVERNORATES) {
      stop(model_name, " / ", stressor_label, ": only ", nrow(gov_table), " governorates remain after the common-support restriction.")
    }

    boot <- bootstrap_cluster_mean(gov_table$mean_change, BOOTSTRAP_REPS,
                                   BOOTSTRAP_SEED + model_index * 100L + (if (stressor == "temperature") 1L else 2L))
    boot_summary <- summarize_effect(boot, "negative")
    mean_change <- mean(gov_table$mean_change)
    baseline_mean <- mean(gov_table$baseline_mean)
    stress_mean <- mean(gov_table$stressed_mean)
    coverage_pct <- 100 * scenario$n_eligible / scenario$n_original
    loo_values <- vapply(seq_len(nrow(gov_table)), function(i) {
      mean(gov_table$mean_change[-i])
    }, numeric(1))
    loo_negative_pct <- mean(loo_values < 0)

    bootstrap_ci_supports_negative <- is.finite(boot_summary$ci_high) && boot_summary$ci_high < 0
    loo_supports_negative <- is.finite(loo_negative_pct) && loo_negative_pct >= MIN_LOO_DIRECTIONAL_STABILITY
    model_supported <- is.finite(mean_change) && mean_change < 0 &&
      is.finite(boot_summary$stability) && boot_summary$stability >= MIN_DIRECTIONAL_STABILITY &&
      bootstrap_ci_supports_negative &&
      loo_supports_negative &&
      coverage_pct >= MIN_SCENARIO_COVERAGE_PCT &&
      nrow(gov_table) >= MIN_GOVERNORATES &&
      isTRUE(meta$transformation %in% c("raw", "log1p", "governorate_z"))

    overall_rows[[length(overall_rows) + 1L]] <- data.frame(
      model = model_name,
      stressor = stressor_label,
      response_transformation = meta$transformation,
      effect_scale = "production_1000_qt",
      n_original = scenario$n_original,
      n_eligible = scenario$n_eligible,
      coverage_pct = coverage_pct,
      governorates = nrow(gov_table),
      baseline_mean_prediction = baseline_mean,
      stressed_mean_prediction = stress_mean,
      mean_change_1000_qt = mean_change,
      mean_change_pct = if (baseline_mean != 0) 100 * mean_change / baseline_mean else NA_real_,
      median_governorate_change_1000_qt = stats::median(gov_table$mean_change),
      bootstrap_ci_95_low = boot_summary$ci_low,
      bootstrap_ci_95_high = boot_summary$ci_high,
      bootstrap_negative_direction_pct = 100 * boot_summary$negative_pct,
      bootstrap_positive_direction_pct = 100 * boot_summary$positive_pct,
      bootstrap_ci_upper_below_zero = bootstrap_ci_supports_negative,
      loo_negative_direction_pct = 100 * loo_negative_pct,
      loo_negative_stability_threshold_pct = 100 * MIN_LOO_DIRECTIONAL_STABILITY,
      model_component_supported = model_supported,
      stringsAsFactors = FALSE
    )

    regional_summary <- do.call(rbind, lapply(c("North", "Center_South"), function(region_name) {
      z <- gov_table[gov_table$region == region_name, , drop = FALSE]
      if (nrow(z) == 0L) return(NULL)
      group_boot <- bootstrap_cluster_mean(
        z$mean_change, BOOTSTRAP_REPS,
        BOOTSTRAP_SEED + model_index * 1000L + (if (stressor == "temperature") 11L else 12L) + (if (region_name == "North") 1L else 2L)
      )
      group_ci <- summarize_effect(group_boot, "negative")
      data.frame(
        model = model_name,
        stressor = stressor_label,
        region = region_name,
        governorates = nrow(z),
        mean_change_1000_qt = mean(z$mean_change),
        bootstrap_ci_95_low = group_ci$ci_low,
        bootstrap_ci_95_high = group_ci$ci_high,
        negative_direction_pct = 100 * group_ci$negative_pct,
        stringsAsFactors = FALSE
      )
    }))
    regional_rows[[length(regional_rows) + 1L]] <- regional_summary

    contrast <- bootstrap_regional_contrast(
      gov_table, BOOTSTRAP_REPS,
      BOOTSTRAP_SEED + model_index * 10000L + (if (stressor == "temperature") 21L else 22L)
    )
    ci_excludes_zero <- is.finite(contrast$ci_low) && is.finite(contrast$ci_high) &&
      (contrast$ci_low > 0 || contrast$ci_high < 0)
    heterogeneity_supported <- ci_excludes_zero &&
      is.finite(contrast$stability) && contrast$stability >= MIN_DIRECTIONAL_STABILITY

    heterogeneity_rows[[length(heterogeneity_rows) + 1L]] <- data.frame(
      model = model_name,
      stressor = stressor_label,
      contrast = "North minus Center/South effect",
      north_governorates = contrast$north_n,
      center_south_governorates = contrast$center_south_n,
      north_minus_center_south_1000_qt = contrast$estimate,
      bootstrap_ci_95_low = contrast$ci_low,
      bootstrap_ci_95_high = contrast$ci_high,
      bootstrap_sign_stability_pct = 100 * contrast$stability,
      regional_heterogeneity_supported = heterogeneity_supported,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(gov_table))) {
      loo_rows[[length(loo_rows) + 1L]] <- data.frame(
        model = model_name,
        stressor = stressor_label,
        omitted_governorate = gov_table$governorate[i],
        mean_change_after_omission_1000_qt = mean(gov_table$mean_change[-i]),
        expected_negative_direction_retained = mean(gov_table$mean_change[-i]) < 0,
        stringsAsFactors = FALSE
      )
    }

    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      model = model_name,
      stressor = stressor_label,
      scenario_variable = scenario$variable,
      predictor_change = scenario$shift,
      support_low = scenario$support_low,
      support_high = scenario$support_high,
      observations_original = scenario$n_original,
      observations_eligible = scenario$n_eligible,
      coverage_pct = coverage_pct,
      model_response_retransformed_to_production = TRUE,
      stringsAsFactors = FALSE
    )
  }
}

overall_table <- do.call(rbind, overall_rows)
regional_table <- do.call(rbind, regional_rows)
heterogeneity_table <- do.call(rbind, heterogeneity_rows)
loo_table <- do.call(rbind, loo_rows)
diagnostics_table <- do.call(rbind, diagnostic_rows)
model_sample_audit <- do.call(rbind, model_sample_rows)
response_transformation_audit <- do.call(rbind, lapply(names(metadata), function(model_name) {
  meta <- metadata[[model_name]]
  inverse_description <- switch(
    meta$transformation,
    raw = "Identity: y_model is production_1000_qt.",
    log1p = if (model_name == "QR") "Monotone inverse expm1(); matches median quantile forecasting workflow." else paste0("Monotone inverse exp(eta) * Duan smearing factor - 1; smearing factor = ", signif(meta$smearing_factor, 8), "."),
    governorate_z = "Monotone governorate-specific affine inverse using saved training production mean and SD."
  )
  data.frame(
    model = model_name,
    model_file = basename(model_paths[[model_name]]),
    metadata_bundle = basename(bundle_paths[[model_name]]),
    response_expression = paste(deparse(stats::formula(models[[model_name]])[[2L]]), collapse = ""),
    response_transformation = meta$transformation,
    inverse_mapping = inverse_description,
    monotone_increasing_inverse_verified = TRUE,
    production_scale = "production_1000_qt",
    stringsAsFactors = FALSE
  )
}))

summarize_component <- function(stressor_label) {
  z <- overall_table[overall_table$stressor == stressor_label, , drop = FALSE]
  support_n <- sum(z$model_component_supported, na.rm = TRUE)
  assessed_n <- nrow(z)
  supported <- assessed_n >= MIN_MODELS_FOR_COMPONENT_SUPPORT &&
    support_n >= MIN_MODELS_FOR_COMPONENT_SUPPORT
  data.frame(
    stressor = stressor_label,
    models_assessed = assessed_n,
    models_supporting_expected_negative_production_effect = support_n,
    required_model_count = MIN_MODELS_FOR_COMPONENT_SUPPORT,
    required_bootstrap_directional_stability = MIN_DIRECTIONAL_STABILITY,
    minimum_scenario_coverage_pct = MIN_SCENARIO_COVERAGE_PCT,
    component_supported = supported,
    supporting_models = if (support_n > 0L) paste(z$model[z$model_component_supported], collapse = "; ") else "None",
    stringsAsFactors = FALSE
  )
}

component_table <- rbind(
  summarize_component("Warmer December (+1 C)"),
  summarize_component("Worse autumn soil moisture (-0.5 sm_z)")
)

audit_regional_consensus <- function(stressor_label) {
  z <- heterogeneity_table[heterogeneity_table$stressor == stressor_label, , drop = FALSE]
  supported <- z[z$regional_heterogeneity_supported, , drop = FALSE]
  signs <- sign(supported$north_minus_center_south_1000_qt)
  positive_n <- sum(signs > 0, na.rm = TRUE)
  negative_n <- sum(signs < 0, na.rm = TRUE)
  majority_sign_n <- max(positive_n, negative_n)
  consensus_supported <- nrow(supported) >= MIN_MODELS_FOR_COMPONENT_SUPPORT &&
    majority_sign_n >= MIN_MODELS_FOR_COMPONENT_SUPPORT
  data.frame(
    stressor = stressor_label,
    models_assessed = nrow(z),
    models_with_95pct_regional_contrast_interval_excluding_zero = nrow(supported),
    models_agreeing_on_contrast_direction = majority_sign_n,
    required_model_count = MIN_MODELS_FOR_COMPONENT_SUPPORT,
    regional_heterogeneity_supported = consensus_supported,
    north_more_negative = if (consensus_supported && negative_n >= MIN_MODELS_FOR_COMPONENT_SUPPORT) TRUE else if (consensus_supported && positive_n >= MIN_MODELS_FOR_COMPONENT_SUPPORT) FALSE else NA,
    stringsAsFactors = FALSE
  )
}

regional_consensus_table <- rbind(
  audit_regional_consensus("Warmer December (+1 C)"),
  audit_regional_consensus("Worse autumn soil moisture (-0.5 sm_z)")
)

temperature_supported <- isTRUE(component_table$component_supported[component_table$stressor == "Warmer December (+1 C)"])
moisture_supported <- isTRUE(component_table$component_supported[component_table$stressor == "Worse autumn soil moisture (-0.5 sm_z)"])
temperature_regional_supported <- isTRUE(regional_consensus_table$regional_heterogeneity_supported[regional_consensus_table$stressor == "Warmer December (+1 C)"])
moisture_regional_supported <- isTRUE(regional_consensus_table$regional_heterogeneity_supported[regional_consensus_table$stressor == "Worse autumn soil moisture (-0.5 sm_z)"])
regional_supported <- temperature_regional_supported || moisture_regional_supported

if (temperature_supported && moisture_supported && regional_supported) {
  overall_classification <- "FULL H1 SUPPORT"
  overall_conclusion <- paste(
    "The model-based evidence supports all operationalized H1 components: warmer December temperatures and worse autumn soil moisture each reduce predicted durum-wheat production,",
    "and at least one of these climate sensitivities differs reproducibly between North and Center/South."
  )
} else if (temperature_supported || moisture_supported || regional_supported) {
  overall_classification <- "PARTIAL H1 SUPPORT"
  supported_parts <- character(0)
  unsupported_parts <- character(0)
  if (temperature_supported) supported_parts <- c(supported_parts, "the December-temperature production response") else unsupported_parts <- c(unsupported_parts, "the December-temperature production response")
  if (moisture_supported) supported_parts <- c(supported_parts, "the autumn soil-moisture-deficit production response") else unsupported_parts <- c(unsupported_parts, "the autumn soil-moisture-deficit production response")
  if (regional_supported) supported_parts <- c(supported_parts, "regional heterogeneity in at least one climate response") else unsupported_parts <- c(unsupported_parts, "regional heterogeneity")
  overall_conclusion <- paste0(
    "The fitted models provide support for ", paste(supported_parts, collapse = "; "),
    ". The following H1 component(s) were not established: ", paste(unsupported_parts, collapse = "; "), "."
  )
} else {
  overall_classification <- "NO H1 SUPPORT ESTABLISHED"
  overall_conclusion <- paste(
    "The predefined model-based criteria were not met for the December-temperature response, the autumn soil-moisture response, or regional heterogeneity.",
    "This means that H1 is not corroborated by this model-based analysis; it does not prove that the underlying effects are absent."
  )
}

if (temperature_supported) {
  temperature_conclusion <- paste(
    "PRIMARY RESULT — DECEMBER TEMPERATURE: SUPPORTED.",
    "At least two of the three fitted models show a negative production-scale response to a +1 C December-temperature scenario,",
    "with governorate-bootstrap directional stability at or above the predefined 90% threshold and at least 70% common-support coverage."
  )
} else {
  temperature_conclusion <- paste(
    "PRIMARY RESULT — DECEMBER TEMPERATURE: NOT ESTABLISHED.",
    "The cross-model requirement was not met: fewer than two fitted models satisfied the expected negative production-scale effect, 90% governorate-bootstrap directional stability, and 70% common-support coverage together."
  )
}

if (temperature_supported && !moisture_supported) {
  temperature_conclusion <- paste(temperature_conclusion,
    " The partial H1 result is specifically driven by the December-temperature component; the soil-moisture component was not corroborated.")
} else if (temperature_supported && moisture_supported && !regional_supported) {
  temperature_conclusion <- paste(temperature_conclusion,
    " Both primary climate directions are supported, but regional sensitivity differences were not established; therefore the overall conclusion is partial H1 support, not full H1 support.")
} else if (!temperature_supported && moisture_supported) {
  temperature_conclusion <- paste(temperature_conclusion,
    " Soil-moisture evidence may support H1 partially, but the temperature component that is the primary focus of this study is not established by these models.")
}

criteria_lines <- c(
  "PRE-SPECIFIED MODEL-BASED DECISION RULES",
  paste0("A model supports a climate component when its mean production-scale scenario effect is negative, its 95% governorate-bootstrap interval lies entirely below zero, bootstrap directional stability is at least ", MIN_DIRECTIONAL_STABILITY * 100, "%, at least ", MIN_LOO_DIRECTIONAL_STABILITY * 100, "% of leave-one-governorate-out estimates remain negative, common-support coverage is at least ", MIN_SCENARIO_COVERAGE_PCT, "%, and at least ", MIN_GOVERNORATES, " governorates are represented."),
  paste0("A cross-model climate component is supported when at least ", MIN_MODELS_FOR_COMPONENT_SUPPORT, " of the three models satisfy the model-level rule."),
  paste0("Regional heterogeneity is supported for a stressor when at least ", MIN_MODELS_FOR_COMPONENT_SUPPORT, " models have a 95% governorate-bootstrap interval for North minus Center/South that excludes zero, with the same contrast direction in at least ", MIN_MODELS_FOR_COMPONENT_SUPPORT, " models."),
  "Full H1 support requires temperature support, soil-moisture support, and regional heterogeneity for at least one of the two climate stressors.",
  "If only some components meet the criteria, the decision is partial H1 support and the missing component(s) are named explicitly.",
  "No p-values or null-hypothesis significance tests are used."
)

report_lines <- c(
  "MODEL-BASED H1 INFERENCE — FINAL CONCLUSION",
  "",
  paste0("OVERALL CLASSIFICATION: ", overall_classification),
  overall_conclusion,
  "",
  temperature_conclusion,
  "",
  "COMPONENT DECISION TABLE",
  capture.output(print(component_table, row.names = FALSE)),
  "",
  "REGIONAL HETEROGENEITY DECISION TABLE",
  capture.output(print(regional_consensus_table, row.names = FALSE)),
  "",
  "RESPONSE TRANSFORMATION VERIFICATION",
  capture.output(print(response_transformation_audit, row.names = FALSE)),
  "",
  "MODEL-SPECIFIC PRODUCTION-SCALE EFFECTS",
  capture.output(print(overall_table, row.names = FALSE)),
  "",
  "MODEL-SPECIFIC NORTH MINUS CENTER/SOUTH CONTRASTS",
  capture.output(print(heterogeneity_table, row.names = FALSE)),
  "",
  criteria_lines,
  "",
  "METHOD AND INTERPRETATION",
  paste0("Each fitted model is evaluated on its original 2000-2012 training observations (expected N = ", EXPECTED_TRAIN_N, ", ", EXPECTED_GOVERNORATES, " governorates). The observed covariates are held fixed except for the scenario variable."),
  "Temperature scenario: increase December temperature by 1 C. Soil-moisture scenario: decrease the governorate-standardized autumn soil-moisture value sm_z by 0.5 standard deviations, then recompute the deficit, wetness and regional interaction features used by the fitted model.",
  "The GAM and mixed-effects log1p predictions use the Duan smearing correction used by their forecasting workflows. QR log1p predictions use expm1(), matching the quantile-regression forecasting workflow. Governorate-standardized responses are converted using the saved governorate-specific training means and standard deviations.",
  "Effects are averaged across governorates equally: the annual scenario effects are first averaged within each governorate, then governorates receive equal weight. The cluster bootstrap resamples governorates, not individual year rows, and the fitted models are held fixed. The model-level decision also requires the 95% cluster-bootstrap interval to lie below zero and at least 80% of leave-one-governorate-out estimates to remain negative.",
  "The 95% bootstrap intervals quantify variation across governorate resamples conditional on the fitted models. They do not include full parameter-estimation uncertainty from refitting each model in every replicate. QR results describe the conditional median response, whereas GAM and mixed-effects predictions represent their fitted response summaries.",
  "The stress scenarios are restricted to common support: both the baseline and shifted predictor values must lie between the 2nd and 98th percentiles of the model's training predictor distribution. Coverage is reported for every model and stressor.",
  "These are model-based adjusted associations, not causal effects. Production is total production and can also vary with harvested area, management and other unmodeled factors.",
  "The results apply to the forecasting models fitted on 2000-2012 and their 15-governorate training design. They are not a refit of the models on the complete 2000-2018 period and are not a replacement for a separate full-period inferential model.",
  "",
  paste0("Dataset: ", data_path),
  paste0("Inference script directory: ", SCRIPT_DIR),
  paste0("Bootstrap replicates: ", BOOTSTRAP_REPS),
  paste0("Bootstrap seed: ", BOOTSTRAP_SEED),
  paste0("Bootstrap directional-stability threshold: ", 100 * MIN_DIRECTIONAL_STABILITY, "%"),
  paste0("Leave-one-governorate-out negative-direction threshold: ", 100 * MIN_LOO_DIRECTIONAL_STABILITY, "%"),
  paste0("Scenario support quantiles: ", paste(SUPPORT_QUANTILES, collapse = " to "))
)

utils::write.csv(model_sample_audit, file.path(OUTPUT_DIR, "H1_model_sample_and_transformation_audit.csv"), row.names = FALSE, na = "")
utils::write.csv(response_transformation_audit, file.path(OUTPUT_DIR, "H1_response_transformation_audit.csv"), row.names = FALSE, na = "")
utils::write.csv(overall_table, file.path(OUTPUT_DIR, "H1_model_scenario_effects_overall.csv"), row.names = FALSE, na = "")
utils::write.csv(regional_table, file.path(OUTPUT_DIR, "H1_model_scenario_effects_by_region.csv"), row.names = FALSE, na = "")
utils::write.csv(heterogeneity_table, file.path(OUTPUT_DIR, "H1_model_regional_heterogeneity_contrasts.csv"), row.names = FALSE, na = "")
utils::write.csv(component_table, file.path(OUTPUT_DIR, "H1_cross_model_directional_consensus.csv"), row.names = FALSE, na = "")
utils::write.csv(regional_consensus_table, file.path(OUTPUT_DIR, "H1_cross_model_regional_consensus.csv"), row.names = FALSE, na = "")
utils::write.csv(loo_table, file.path(OUTPUT_DIR, "H1_leave_one_governorate_out_stability.csv"), row.names = FALSE, na = "")
utils::write.csv(diagnostics_table, file.path(OUTPUT_DIR, "H1_scenario_support_diagnostics.csv"), row.names = FALSE, na = "")

writeLines(report_lines, file.path(OUTPUT_DIR, "H1_FINAL_CONCLUSION.txt"), useBytes = TRUE)
writeLines(c("PRIMARY H1 RESULT — DECEMBER TEMPERATURE", "", temperature_conclusion, "", paste0("Overall H1 classification: ", overall_classification)),
           file.path(OUTPUT_DIR, "H1_temperature_primary_result.txt"), useBytes = TRUE)

cat("\n========== FINAL H1 CLASSIFICATION ==========\n")
cat(overall_classification, "\n")
cat(overall_conclusion, "\n\n")
cat(temperature_conclusion, "\n")
cat("\n========== COMPONENT DECISIONS ==========\n")
print(component_table, row.names = FALSE)
cat("\n========== REGIONAL HETEROGENEITY ==========\n")
print(regional_consensus_table, row.names = FALSE)
cat("\n========== PRODUCTION-SCALE SCENARIO EFFECTS ==========\n")
print(overall_table, row.names = FALSE)
cat("\n========== OUTPUT DIRECTORY ==========\n")
cat(normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), "\n")
