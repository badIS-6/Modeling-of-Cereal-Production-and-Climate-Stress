options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_DIR <- "/home/badisz/Downloads/Modeling"
DATA_PATH <- file.path(PROJECT_DIR, "model_table.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Mixed_Effects_Optimized_Output_nlme")
TRAIN_START_YEAR <- 2000L
TRAIN_END_YEAR <- 2012L
TEST_START_YEAR <- 2013L
TEST_END_YEAR <- 2018L
EXPECTED_N <- 281L
EXPECTED_TRAIN_N <- 192L
EXPECTED_TEST_N <- 89L
EXPECTED_GOVERNORATES <- 15L
EXPECTED_ZERO_PRODUCTION <- 11L
SEED <- 20261126L
BOUNDARY_SD_RATIO <- 0.001

if (!requireNamespace("nlme", quietly = TRUE)) {
  stop("The nlme package is required. Install it with install.packages('nlme').")
}

if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory not found: ", PROJECT_DIR)
}
if (!file.exists(DATA_PATH)) {
  stop("Dataset not found: ", DATA_PATH)
}
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

write_csv <- function(x, filename) {
  utils::write.csv(
    x,
    file = file.path(OUTPUT_DIR, filename),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

capture_text <- function(object, filename) {
  capture.output(object, file = file.path(OUTPUT_DIR, filename))
}

data <- utils::read.csv(
  DATA_PATH,
  fileEncoding = "UTF-8",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_cols <- c(
  "governorate", "crop", "harvest_year", "production_1000_qt",
  "rain_jan_mar_mm", "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)

missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (nrow(data) != EXPECTED_N) {
  stop("Expected ", EXPECTED_N, " observations; found ", nrow(data), ".")
}
if (anyNA(data[, required_cols, drop = FALSE])) {
  stop("Missing values were found in required columns.")
}

numeric_cols <- c(
  "harvest_year", "production_1000_qt", "rain_jan_mar_mm",
  "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)

if (!all(vapply(data[, numeric_cols, drop = FALSE], is.numeric, logical(1)))) {
  stop("One or more required numeric columns were imported as non-numeric.")
}
if (!all(vapply(data[, numeric_cols, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
  stop("Non-finite values were found in required numeric columns.")
}
if (!all(data$crop == "durum_wheat")) {
  stop("The dataset contains crop categories other than durum_wheat.")
}
if (any(data$production_1000_qt < 0)) {
  stop("Negative production values were found; inspect the target before fitting.")
}
if (sum(data$production_1000_qt == 0) != EXPECTED_ZERO_PRODUCTION) {
  stop(
    "Expected ", EXPECTED_ZERO_PRODUCTION,
    " zero-production observations; found ", sum(data$production_1000_qt == 0), "."
  )
}
if (any(data$harvest_year != as.integer(data$harvest_year))) {
  stop("harvest_year contains non-integer values.")
}
data$harvest_year <- as.integer(data$harvest_year)

expected_governorates <- c(
  "Béja", "Gabès", "Gafsa", "Jendouba", "Kairouan",
  "Kasserine", "Le Kef", "Manouba", "Medenine", "Nabeul",
  "Sidi Bouzid", "Siliana", "Sousse", "Tunis", "Zaghouan"
)

north_governorates <- c(
  "Béja", "Jendouba", "Le Kef", "Siliana", "Zaghouan",
  "Manouba", "Nabeul", "Tunis"
)

center_south_governorates <- c(
  "Kairouan", "Kasserine", "Sidi Bouzid", "Sousse",
  "Gafsa", "Gabès", "Medenine"
)

observed_governorates <- sort(unique(trimws(data$governorate)))
if (anyNA(data$governorate) || any(!nzchar(trimws(data$governorate)))) {
  stop("Empty or missing governorate names were found.")
}
if (length(observed_governorates) != EXPECTED_GOVERNORATES) {
  stop(
    "Expected ", EXPECTED_GOVERNORATES,
    " governorates; found ", length(observed_governorates), "."
  )
}
if (!setequal(observed_governorates, expected_governorates)) {
  stop(
    "Governorate names do not match the expected set. Unexpected: ",
    paste(setdiff(observed_governorates, expected_governorates), collapse = ", "),
    ". Missing: ",
    paste(setdiff(expected_governorates, observed_governorates), collapse = ", ")
  )
}

data$governorate <- trimws(data$governorate)
duplicate_keys <- duplicated(data[, c("governorate", "harvest_year")])
if (any(duplicate_keys)) {
  duplicate_rows <- data[duplicate_keys, c("governorate", "harvest_year"), drop = FALSE]
  write_csv(duplicate_rows, "duplicate_governorate_year_records.csv")
  stop("Duplicate governorate-year records were found. See duplicate_governorate_year_records.csv.")
}

if (any(data$harvest_year < TRAIN_START_YEAR | data$harvest_year > TEST_END_YEAR)) {
  stop("Observations outside the 2000-2018 study period were found.")
}

if (!setequal(sort(unique(data$harvest_year)), seq.int(TRAIN_START_YEAR, TEST_END_YEAR))) {
  stop("The dataset does not contain every harvest year from 2000 through 2018.")
}



governorate_levels <- expected_governorates
data$governorate <- factor(data$governorate, levels = governorate_levels)
data$region <- factor(
  ifelse(as.character(data$governorate) %in% north_governorates, "North", "Center_South"),
  levels = c("North", "Center_South")
)
contrasts(data$region) <- stats::contr.sum(2)
data <- data[order(data$harvest_year, data$governorate), , drop = FALSE]
rownames(data) <- NULL

coverage_matrix <- table(
  factor(data$governorate, levels = governorate_levels),
  factor(data$harvest_year, levels = seq.int(TRAIN_START_YEAR, TEST_END_YEAR))
)
coverage_table <- data.frame(
  governorate = rownames(coverage_matrix),
  Production_Years = rowSums(coverage_matrix > 0),
  Missing_Years = apply(coverage_matrix, 1, function(x) {
    missing <- colnames(coverage_matrix)[x == 0]
    if (length(missing)) paste(missing, collapse = ";") else ""
  }),
  Training_Years = rowSums(coverage_matrix[, as.character(TRAIN_START_YEAR:TRAIN_END_YEAR), drop = FALSE] > 0),
  Test_Years = rowSums(coverage_matrix[, as.character(TEST_START_YEAR:TEST_END_YEAR), drop = FALSE] > 0),
  row.names = NULL
)
write_csv(coverage_table, "governorate_production_coverage.csv")

train_data_raw <- data[
  data$harvest_year >= TRAIN_START_YEAR & data$harvest_year <= TRAIN_END_YEAR,
  , drop = FALSE
]
test_data_raw <- data[
  data$harvest_year >= TEST_START_YEAR & data$harvest_year <= TEST_END_YEAR,
  , drop = FALSE
]

if (nrow(train_data_raw) != EXPECTED_TRAIN_N || nrow(test_data_raw) != EXPECTED_TEST_N) {
  stop(
    "Expected temporal split of ", EXPECTED_TRAIN_N, "/", EXPECTED_TEST_N,
    "; found ", nrow(train_data_raw), "/", nrow(test_data_raw), "."
  )
}
if (!setequal(as.character(train_data_raw$governorate), as.character(test_data_raw$governorate))) {
  stop("Training and test periods do not contain the same governorates.")
}
if (any(table(train_data_raw$governorate) == 0) || any(table(test_data_raw$governorate) == 0)) {
  stop("At least one governorate is absent from the training or test period.")
}

split_audit <- data.frame(
  Dataset = c("Full data", "Training", "Test"),
  First_Year = c(min(data$harvest_year), min(train_data_raw$harvest_year), min(test_data_raw$harvest_year)),
  Last_Year = c(max(data$harvest_year), max(train_data_raw$harvest_year), max(test_data_raw$harvest_year)),
  Observations = c(nrow(data), nrow(train_data_raw), nrow(test_data_raw)),
  Governorates = c(
    length(unique(data$governorate)),
    length(unique(train_data_raw$governorate)),
    length(unique(test_data_raw$governorate))
  ),
  Zero_Production = c(
    sum(data$production_1000_qt == 0),
    sum(train_data_raw$production_1000_qt == 0),
    sum(test_data_raw$production_1000_qt == 0)
  )
)
write_csv(split_audit, "temporal_split_audit.csv")

fit_scaling_stats <- function(train) {
  sm_groups <- split(train$sm_oct_nov, as.character(train$governorate), drop = TRUE)
  sm_stats <- do.call(rbind, lapply(names(sm_groups), function(governorate_name) {
    values <- sm_groups[[governorate_name]]
    data.frame(
      governorate = governorate_name,
      sm_mean = mean(values),
      sm_sd = stats::sd(values),
      stringsAsFactors = FALSE
    )
  }))

  y_groups <- split(train$production_1000_qt, as.character(train$governorate), drop = TRUE)
  y_stats <- do.call(rbind, lapply(names(y_groups), function(governorate_name) {
    values <- y_groups[[governorate_name]]
    data.frame(
      governorate = governorate_name,
      production_mean = mean(values),
      production_sd = stats::sd(values),
      stringsAsFactors = FALSE
    )
  }))

  if (any(!is.finite(sm_stats$sm_mean)) || any(!is.finite(sm_stats$sm_sd)) || any(sm_stats$sm_sd <= 0)) {
    stop("Invalid training-period soil-moisture scaling parameters.")
  }
  if (any(!is.finite(y_stats$production_mean))) {
    stop("Invalid training-period governorate production means.")
  }

  list(sm = sm_stats, y = y_stats)
}

apply_features <- function(newdata, stats_obj, response_mode, include_response = TRUE, year_origin = 2006L) {
  d <- newdata
  sm_index <- match(as.character(d$governorate), stats_obj$sm$governorate)
  if (anyNA(sm_index)) {
    stop("A governorate has no training-period soil-moisture scaling parameters.")
  }

  d$sm_z <- (d$sm_oct_nov - stats_obj$sm$sm_mean[sm_index]) / stats_obj$sm$sm_sd[sm_index]
  d$sm_deficit <- pmax(-1 - d$sm_z, 0)
  d$sm_wet <- pmax(d$sm_z - 1, 0)
  d$year_c <- (d$harvest_year - year_origin) / 10
  d$governorate <- factor(as.character(d$governorate), levels = governorate_levels)
  d$region <- factor(as.character(d$region), levels = c("North", "Center_South"))
  contrasts(d$region) <- stats::contr.sum(2)

  if (include_response) {
    if (response_mode == "raw") {
      d$y_model <- d$production_1000_qt
    } else if (response_mode == "log1p") {
      d$y_model <- log1p(d$production_1000_qt)
    } else if (response_mode == "governorate_z") {
      y_index <- match(as.character(d$governorate), stats_obj$y$governorate)
      if (anyNA(y_index)) {
        stop("A governorate has no training-period production scaling parameters.")
      }
      y_sd <- stats_obj$y$production_sd[y_index]
      if (any(!is.finite(y_sd)) || any(y_sd <= 0)) {
        stop("Governorate-specific production standardization is undefined in this training fold.")
      }
      d$y_model <- (d$production_1000_qt - stats_obj$y$production_mean[y_index]) / y_sd
    } else {
      stop("Unknown response transformation: ", response_mode)
    }
  } else {
    d$y_model <- NA_real_
  }

  if (any(!is.finite(d$sm_z)) || any(!is.finite(d$sm_deficit)) || any(!is.finite(d$sm_wet))) {
    stop("Non-finite soil-moisture features were generated.")
  }
  d
}

fixed_specs <- list(
  baseline_climate = paste(
    "region * (temp_dec_c + sm_z + ndvi_mean) +",
    "I(sm_z^2) + rain_jan_mar_mm + vhi + year_c"
  ),
  threshold_climate = paste(
    "region * (temp_dec_c + sm_deficit + sm_wet + ndvi_mean + vhi) +",
    "rain_jan_mar_mm + year_c"
  ),
  parsimonious_threshold = paste(
    "region * (temp_dec_c + sm_deficit + sm_wet) +",
    "ndvi_mean + vhi + rain_jan_mar_mm + year_c"
  )
)

random_specs <- c("random_intercept", "random_year_slope", "random_moisture_slope")
correlation_specs <- c("independent", "AR1")
response_modes <- c("raw", "log1p", "governorate_z")

candidates <- expand.grid(
  Fixed_Specification = names(fixed_specs),
  Random_Structure = random_specs,
  Residual_Correlation = correlation_specs,
  Response_Transformation = response_modes,
  stringsAsFactors = FALSE
)
candidates$Candidate_ID <- sprintf("M%02d", seq_len(nrow(candidates)))
candidates$Model_Name <- paste(
  candidates$Fixed_Specification,
  candidates$Random_Structure,
  candidates$Residual_Correlation,
  candidates$Response_Transformation,
  sep = "__"
)

make_fixed_formula <- function(candidate) {
  stats::as.formula(paste("y_model ~", fixed_specs[[candidate$Fixed_Specification]]))
}

make_random_structure <- function(random_name) {
  switch(
    random_name,
    random_intercept = stats::as.formula("~ 1 | governorate"),
    random_year_slope = list(governorate = nlme::pdDiag(~ 1 + year_c)),
    random_moisture_slope = list(governorate = nlme::pdDiag(~ 1 + sm_z)),
    stop("Unknown random structure: ", random_name)
  )
}

make_correlation_structure <- function(correlation_name) {
  switch(
    correlation_name,
    independent = NULL,
    AR1 = nlme::corAR1(value = 0.2, form = ~ harvest_year | governorate),
    stop("Unknown residual correlation structure: ", correlation_name)
  )
}

lme_control <- nlme::lmeControl(
  opt = "optim",
  maxIter = 200,
  msMaxIter = 500,
  niterEM = 50,
  tolerance = 1e-7,
  msTol = 1e-7,
  returnObject = FALSE
)

fit_mixed_model <- function(fixed_formula, random_name, correlation_name, model_data, method = "ML") {
  model_data <- model_data[order(model_data$governorate, model_data$harvest_year), , drop = FALSE]
  fixed_terms <- stats::delete.response(stats::terms(fixed_formula))
  design_matrix <- stats::model.matrix(fixed_terms, data = model_data)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("The fixed-effects design matrix is rank deficient.")
  }

  fit_warnings <- character(0)
  fit_args <- list(
    fixed = fixed_formula,
    random = make_random_structure(random_name),
    data = model_data,
    method = method,
    control = lme_control,
    na.action = stats::na.fail
  )
  correlation_object <- make_correlation_structure(correlation_name)
  if (!is.null(correlation_object)) {
    fit_args$correlation <- correlation_object
  }

  model <- withCallingHandlers(
    do.call(nlme::lme, fit_args),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  critical_warning <- grepl(
    "converg|iteration limit|false convergence|singular|backsolve|not positive definite|non-positive definite|failed to converge",
    fit_warnings,
    ignore.case = TRUE
  )
  if (any(critical_warning)) {
    stop("Critical nlme warning: ", paste(unique(fit_warnings[critical_warning]), collapse = " | "))
  }

  fixed_coefficients <- nlme::fixed.effects(model)
  if (any(!is.finite(fixed_coefficients))) {
    stop("The fitted model contains non-finite fixed-effect coefficients.")
  }

  ap_var <- model$apVar
  if (is.character(ap_var) && grepl("positive definite|singular|failed", ap_var, ignore.case = TRUE)) {
    stop("The approximate variance-covariance matrix is not reliable: ", ap_var)
  }
  if (is.matrix(ap_var) && any(!is.finite(ap_var))) {
    stop("The approximate variance-covariance matrix contains non-finite values.")
  }

  attr(model, "fit_warnings") <- fit_warnings
  list(model = model, warnings = fit_warnings)
}

random_effect_diagnostics <- function(model, boundary_ratio = BOUNDARY_SD_RATIO) {
  vc <- as.matrix(nlme::VarCorr(model))
  if (is.null(rownames(vc)) || !("StdDev" %in% colnames(vc))) {
    stop("Could not extract random-effect standard deviations from VarCorr().")
  }

  sd_values <- suppressWarnings(as.numeric(vc[, "StdDev"]))
  names(sd_values) <- rownames(vc)
  residual_index <- grep("Residual", names(sd_values), ignore.case = TRUE)
  if (length(residual_index) != 1L || !is.finite(sd_values[residual_index]) || sd_values[residual_index] <= 0) {
    stop("Could not identify a valid residual standard deviation.")
  }

  random_values <- sd_values[-residual_index]
  residual_sd <- unname(sd_values[residual_index])
  ratios <- random_values / residual_sd
  boundary <- !length(ratios) || any(!is.finite(ratios)) || any(ratios < boundary_ratio)

  output <- data.frame(
    Random_Effect = names(random_values),
    Random_SD = unname(random_values),
    Residual_SD = residual_sd,
    Random_to_Residual_SD_Ratio = unname(ratios),
    Boundary_Flag = if (length(ratios)) unname(ratios < boundary_ratio) else logical(0),
    stringsAsFactors = FALSE
  )
  attr(output, "boundary") <- boundary
  output
}

get_smearing_factor <- function(model, response_mode) {
  if (response_mode != "log1p") {
    return(1)
  }
  model_residuals <- as.numeric(stats::residuals(model, level = 1, type = "response"))
  smearing <- mean(exp(model_residuals), na.rm = TRUE)
  if (!is.finite(smearing) || smearing <= 0) {
    stop("Could not calculate a finite Duan smearing factor for log1p predictions.")
  }
  smearing
}

back_transform_predictions <- function(model, newdata, response_mode, target_stats) {
  eta <- as.numeric(stats::predict(model, newdata = newdata, level = 1))
  if (length(eta) != nrow(newdata) || any(!is.finite(eta))) {
    stop("The model returned missing or non-finite predictions.")
  }

  if (response_mode == "raw") {
    return(eta)
  }
  if (response_mode == "log1p") {
    smearing <- get_smearing_factor(model, response_mode)
    prediction <- exp(eta) * smearing - 1
    if (any(!is.finite(prediction))) {
      stop("The log1p retransformation generated non-finite predictions.")
    }
    return(prediction)
  }
  if (response_mode == "governorate_z") {
    y_index <- match(as.character(newdata$governorate), target_stats$governorate)
    if (anyNA(y_index)) {
      stop("Could not match governorates to training-period target-scaling parameters.")
    }
    prediction <- eta * target_stats$production_sd[y_index] + target_stats$production_mean[y_index]
    if (any(!is.finite(prediction))) {
      stop("The governorate-standardized retransformation generated non-finite predictions.")
    }
    return(prediction)
  }
  stop("Unknown response transformation: ", response_mode)
}

calculate_metrics <- function(actual, predicted, label) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors have different lengths.")
  }
  if (any(!is.finite(actual)) || any(!is.finite(predicted))) {
    stop("Metrics cannot be calculated from non-finite values.")
  }
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  data.frame(
    Prediction = label,
    N = length(actual),
    R2 = if (sst > 0) 1 - sse / sst else NA_real_,
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted)),
    Mean_Error = mean(predicted - actual),
    row.names = NULL
  )
}

complexity_fixed <- c(
  baseline_climate = 3L,
  threshold_climate = 2L,
  parsimonious_threshold = 1L
)
complexity_random <- c(
  random_intercept = 0L,
  random_year_slope = 1L,
  random_moisture_slope = 1L
)
complexity_correlation <- c(independent = 0L, AR1 = 1L)
complexity_response <- c(raw = 0L, log1p = 1L, governorate_z = 1L)

folds <- list(
  list(name = "2000-2006_to_2007-2008", train_end = 2006L, valid_start = 2007L, valid_end = 2008L),
  list(name = "2000-2008_to_2009-2010", train_end = 2008L, valid_start = 2009L, valid_end = 2010L),
  list(name = "2000-2010_to_2011-2012", train_end = 2010L, valid_start = 2011L, valid_end = 2012L)
)

cv_predictions_by_candidate <- vector("list", nrow(candidates))
candidate_summary_rows <- vector("list", nrow(candidates))
fold_diagnostic_rows <- list()
baseline_cv_predictions <- lapply(folds, function(fold) {
  tr_raw <- train_data_raw[train_data_raw$harvest_year <= fold$train_end, , drop = FALSE]
  va_raw <- train_data_raw[
    train_data_raw$harvest_year >= fold$valid_start &
      train_data_raw$harvest_year <= fold$valid_end,
    , drop = FALSE
  ]
  baseline_means <- tapply(tr_raw$production_1000_qt, as.character(tr_raw$governorate), mean)
  baseline_pred <- as.numeric(baseline_means[as.character(va_raw$governorate)])
  if (!nrow(tr_raw) || !nrow(va_raw) || any(!is.finite(baseline_pred))) {
    stop("Could not calculate the fold-specific governorate-mean baseline.")
  }
  data.frame(
    Fold = fold$name,
    harvest_year = va_raw$harvest_year,
    governorate = as.character(va_raw$governorate),
    actual = va_raw$production_1000_qt,
    predicted = baseline_pred,
    stringsAsFactors = FALSE
  )
})

for (i in seq_len(nrow(candidates))) {
  candidate <- candidates[i, , drop = FALSE]
  fixed_formula_i <- make_fixed_formula(candidate)
  candidate_predictions <- list()
  candidate_failed <- FALSE
  candidate_failure <- NA_character_
  candidate_boundary <- FALSE
  candidate_boundary_folds <- 0L
  candidate_warnings <- character(0)

  for (j in seq_along(folds)) {
    fold <- folds[[j]]
    tr_raw <- train_data_raw[
      train_data_raw$harvest_year <= fold$train_end,
      , drop = FALSE
    ]
    va_raw <- train_data_raw[
      train_data_raw$harvest_year >= fold$valid_start &
        train_data_raw$harvest_year <= fold$valid_end,
      , drop = FALSE
    ]

    fold_result <- tryCatch({
      if (!nrow(tr_raw) || !nrow(va_raw)) {
        stop("The fold has an empty training or validation partition.")
      }
      if (!setequal(as.character(tr_raw$governorate), as.character(va_raw$governorate))) {
        stop("A governorate is absent from the fold training or validation partition.")
      }

      scaling_stats <- fit_scaling_stats(tr_raw)
      tr <- apply_features(tr_raw, scaling_stats, candidate$Response_Transformation, include_response = TRUE)
      va <- apply_features(va_raw, scaling_stats, candidate$Response_Transformation, include_response = FALSE)
      fit_info <- fit_mixed_model(
        fixed_formula = fixed_formula_i,
        random_name = candidate$Random_Structure,
        correlation_name = candidate$Residual_Correlation,
        model_data = tr,
        method = "ML"
      )
      boundary_info <- random_effect_diagnostics(fit_info$model)
      boundary_flag <- isTRUE(attr(boundary_info, "boundary"))
      predicted <- back_transform_predictions(
        fit_info$model,
        va,
        candidate$Response_Transformation,
        scaling_stats$y
      )

      list(
        predictions = data.frame(
          Candidate_ID = candidate$Candidate_ID,
          Model_Name = candidate$Model_Name,
          Fold = fold$name,
          harvest_year = va$harvest_year,
          governorate = as.character(va$governorate),
          actual = va$production_1000_qt,
          predicted = predicted,
          stringsAsFactors = FALSE
        ),
        boundary = boundary_flag,
        boundary_details = boundary_info,
        warnings = fit_info$warnings
      )
    }, error = function(e) {
      list(error = conditionMessage(e))
    })

    if (!is.null(fold_result$error)) {
      candidate_failed <- TRUE
      candidate_failure <- paste0(fold$name, ": ", fold_result$error)
      fold_diagnostic_rows[[length(fold_diagnostic_rows) + 1L]] <- data.frame(
        Candidate_ID = candidate$Candidate_ID,
        Model_Name = candidate$Model_Name,
        Fold = fold$name,
        Fit_Succeeded = FALSE,
        Boundary_Fit = NA,
        Warning_Count = NA_integer_,
        Warning_Details = NA_character_,
        Failure_Reason = fold_result$error,
        stringsAsFactors = FALSE
      )
      break
    }

    candidate_predictions[[j]] <- fold_result$predictions
    candidate_boundary <- candidate_boundary || fold_result$boundary
    candidate_boundary_folds <- candidate_boundary_folds + as.integer(fold_result$boundary)
    candidate_warnings <- c(candidate_warnings, fold_result$warnings)

    fold_diagnostic_rows[[length(fold_diagnostic_rows) + 1L]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fold = fold$name,
      Fit_Succeeded = TRUE,
      Boundary_Fit = fold_result$boundary,
      Warning_Count = length(fold_result$warnings),
      Warning_Details = if (length(fold_result$warnings)) paste(unique(fold_result$warnings), collapse = " | ") else "",
      Failure_Reason = "",
      stringsAsFactors = FALSE
    )
  }

  complete_folds <- length(candidate_predictions) == length(folds) && !candidate_failed
  candidate_cv <- NULL
  if (complete_folds) {
    candidate_cv <- do.call(rbind, candidate_predictions)
    cv_predictions_by_candidate[[i]] <- candidate_cv
    cv_metric <- calculate_metrics(
      candidate_cv$actual,
      candidate_cv$predicted,
      candidate$Model_Name
    )
    candidate_summary_rows[[i]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fixed_Specification = candidate$Fixed_Specification,
      Random_Structure = candidate$Random_Structure,
      Residual_Correlation = candidate$Residual_Correlation,
      Response_Transformation = candidate$Response_Transformation,
      CV_N = cv_metric$N,
      CV_R2 = cv_metric$R2,
      CV_RMSE = cv_metric$RMSE,
      CV_MAE = cv_metric$MAE,
      CV_Mean_Error = cv_metric$Mean_Error,
      Boundary_Folds = candidate_boundary_folds,
      Boundary_Observed = candidate_boundary,
      Successful_Folds = length(candidate_predictions),
      Warning_Count = length(candidate_warnings),
      Warning_Details = if (length(candidate_warnings)) paste(unique(candidate_warnings), collapse = " | ") else "",
      CV_Eligible = !candidate_boundary,
      Failure_Reason = "",
      Complexity_Rank = unname(complexity_fixed[candidate$Fixed_Specification]) +
        unname(complexity_random[candidate$Random_Structure]) +
        unname(complexity_correlation[candidate$Residual_Correlation]) +
        unname(complexity_response[candidate$Response_Transformation]),
      stringsAsFactors = FALSE
    )
  } else {
    candidate_summary_rows[[i]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fixed_Specification = candidate$Fixed_Specification,
      Random_Structure = candidate$Random_Structure,
      Residual_Correlation = candidate$Residual_Correlation,
      Response_Transformation = candidate$Response_Transformation,
      CV_N = NA_real_,
      CV_R2 = NA_real_,
      CV_RMSE = NA_real_,
      CV_MAE = NA_real_,
      CV_Mean_Error = NA_real_,
      Boundary_Folds = candidate_boundary_folds,
      Boundary_Observed = candidate_boundary,
      Successful_Folds = length(candidate_predictions),
      Warning_Count = length(candidate_warnings),
      Warning_Details = if (length(candidate_warnings)) paste(unique(candidate_warnings), collapse = " | ") else "",
      CV_Eligible = FALSE,
      Failure_Reason = candidate_failure,
      Complexity_Rank = unname(complexity_fixed[candidate$Fixed_Specification]) +
        unname(complexity_random[candidate$Random_Structure]) +
        unname(complexity_correlation[candidate$Residual_Correlation]) +
        unname(complexity_response[candidate$Response_Transformation]),
      stringsAsFactors = FALSE
    )
  }

  cat(
    "Candidate", i, "of", nrow(candidates), ":", candidate$Model_Name,
    if (candidate_failed) paste("FAILED:", candidate_failure) else if (candidate_boundary) "completed with boundary fit(s)" else "completed",
    "\n"
  )
}

candidate_cv_summary <- do.call(rbind, candidate_summary_rows)
fold_diagnostics <- if (length(fold_diagnostic_rows)) do.call(rbind, fold_diagnostic_rows) else data.frame()
write_csv(candidate_cv_summary, "mixed_effects_all_candidate_cv_results.csv")
write_csv(fold_diagnostics, "mixed_effects_fold_fit_diagnostics.csv")

eligible_cv <- candidate_cv_summary[candidate_cv_summary$CV_Eligible, , drop = FALSE]
if (!nrow(eligible_cv)) {
  stop("No candidate completed all validation folds without a boundary fit. Inspect the saved candidate and fold diagnostics.")
}

if (length(unique(eligible_cv$CV_N)) != 1L) {
  stop("Eligible candidates were evaluated on different numbers of validation observations.")
}

eligible_cv$Rank_RMSE <- rank(eligible_cv$CV_RMSE, ties.method = "average")
eligible_cv$Rank_MAE <- rank(eligible_cv$CV_MAE, ties.method = "average")
eligible_cv$Average_Rank <- rowMeans(eligible_cv[, c("Rank_RMSE", "Rank_MAE")])
eligible_cv <- eligible_cv[
  order(
    eligible_cv$Average_Rank,
    eligible_cv$CV_RMSE,
    eligible_cv$CV_MAE,
    eligible_cv$Warning_Count,
    eligible_cv$Complexity_Rank
  ),
  , drop = FALSE
]
eligible_cv$CV_Rank <- seq_len(nrow(eligible_cv))
rownames(eligible_cv) <- NULL
write_csv(eligible_cv, "mixed_effects_eligible_candidate_cv_rankings.csv")

top_three <- utils::head(eligible_cv, min(3L, nrow(eligible_cv)))
write_csv(top_three, "mixed_effects_top_three_cv_candidates.csv")

cv_prediction_export <- do.call(
  rbind,
  cv_predictions_by_candidate[!vapply(cv_predictions_by_candidate, is.null, logical(1))]
)
write_csv(cv_prediction_export, "mixed_effects_rolling_origin_cv_predictions.csv")

baseline_cv_export <- do.call(rbind, baseline_cv_predictions)
baseline_cv_metrics <- calculate_metrics(
  baseline_cv_export$actual,
  baseline_cv_export$predicted,
  "Governorate-specific training mean baseline"
)
write_csv(baseline_cv_export, "mixed_effects_rolling_origin_baseline_predictions.csv")
write_csv(baseline_cv_metrics, "mixed_effects_rolling_origin_baseline_metrics.csv")

selected_model <- NULL
selected_candidate <- NULL
selected_formula <- NULL
selected_scaling_stats <- NULL
selected_train <- NULL
selected_test <- NULL
selected_fit_info <- NULL
selected_boundary_info <- NULL
selected_test_predictions <- NULL
selected_smearing_factor <- NULL
final_refit_attempts <- list()

for (rank_index in seq_len(nrow(eligible_cv))) {
  ranked_row <- eligible_cv[rank_index, , drop = FALSE]
  candidate <- candidates[candidates$Candidate_ID == ranked_row$Candidate_ID, , drop = FALSE]
  fixed_formula <- make_fixed_formula(candidate)
  final_scaling <- tryCatch(fit_scaling_stats(train_data_raw), error = function(e) e)

  attempt <- tryCatch({
    if (inherits(final_scaling, "error")) stop(conditionMessage(final_scaling))
    final_train_data <- apply_features(
      train_data_raw,
      final_scaling,
      candidate$Response_Transformation,
      include_response = TRUE
    )
    final_train_data <- final_train_data[
      order(final_train_data$governorate, final_train_data$harvest_year),
      , drop = FALSE
    ]
    final_test_data <- apply_features(
      test_data_raw,
      final_scaling,
      candidate$Response_Transformation,
      include_response = FALSE
    )
    final_fit <- fit_mixed_model(
      fixed_formula = fixed_formula,
      random_name = candidate$Random_Structure,
      correlation_name = candidate$Residual_Correlation,
      model_data = final_train_data,
      method = "REML"
    )
    final_boundary <- random_effect_diagnostics(final_fit$model)
    if (isTRUE(attr(final_boundary, "boundary"))) {
      stop("The final REML refit produced a boundary random-effect variance.")
    }
    final_smearing <- get_smearing_factor(final_fit$model, candidate$Response_Transformation)
    final_test_predictions <- back_transform_predictions(
      final_fit$model,
      final_test_data,
      candidate$Response_Transformation,
      final_scaling$y
    )
    list(
      success = TRUE,
      model = final_fit$model,
      fit_info = final_fit,
      boundary_info = final_boundary,
      formula = fixed_formula,
      candidate = candidate,
      scaling = final_scaling,
      train = final_train_data,
      test = final_test_data,
      predictions = final_test_predictions,
      smearing_factor = final_smearing,
      error = ""
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })

  final_refit_attempts[[length(final_refit_attempts) + 1L]] <- data.frame(
    CV_Rank = ranked_row$CV_Rank,
    Candidate_ID = ranked_row$Candidate_ID,
    Model_Name = ranked_row$Model_Name,
    Final_REML_Success = attempt$success,
    Failure_Reason = attempt$error,
    stringsAsFactors = FALSE
  )

  if (isTRUE(attempt$success)) {
    selected_model <- attempt$model
    selected_fit_info <- attempt$fit_info
    selected_boundary_info <- attempt$boundary_info
    selected_test_predictions <- attempt$predictions
    selected_smearing_factor <- attempt$smearing_factor
    selected_formula <- attempt$formula
    selected_candidate <- attempt$candidate
    selected_scaling_stats <- attempt$scaling
    selected_train <- attempt$train
    selected_test <- attempt$test
    break
  }
}

final_refit_attempts <- do.call(rbind, final_refit_attempts)
write_csv(final_refit_attempts, "mixed_effects_final_refit_attempts.csv")

if (is.null(selected_model)) {
  stop("No CV-eligible candidate could be refitted with REML without a boundary fit. Inspect mixed_effects_final_refit_attempts.csv.")
}

selected_rank_row <- eligible_cv[eligible_cv$Candidate_ID == selected_candidate$Candidate_ID, , drop = FALSE]
predicted_unclipped <- selected_test_predictions
if (length(predicted_unclipped) != nrow(test_data_raw) || any(!is.finite(predicted_unclipped))) {
  stop("The selected model did not produce valid predictions for every test observation.")
}
predicted_clipped <- pmax(predicted_unclipped, 0)
actual <- test_data_raw$production_1000_qt

training_means <- selected_scaling_stats$y
baseline_index <- match(as.character(test_data_raw$governorate), training_means$governorate)
if (anyNA(baseline_index)) {
  stop("Could not match the test governorates to their training-period production means.")
}
baseline_prediction <- training_means$production_mean[baseline_index]

test_performance <- rbind(
  calculate_metrics(actual, predicted_unclipped, "Selected mixed-effects model - unclipped"),
  calculate_metrics(actual, predicted_clipped, "Selected mixed-effects model - clipped at zero"),
  calculate_metrics(actual, baseline_prediction, "Governorate-specific training-period mean baseline")
)
write_csv(test_performance, "mixed_effects_final_test_performance.csv")

predictions <- data.frame(
  harvest_year = test_data_raw$harvest_year,
  governorate = as.character(test_data_raw$governorate),
  region = as.character(test_data_raw$region),
  actual_production_1000_qt = actual,
  predicted_unclipped = predicted_unclipped,
  predicted_clipped = predicted_clipped,
  baseline_governorate_mean = baseline_prediction,
  residual_unclipped = actual - predicted_unclipped,
  residual_clipped = actual - predicted_clipped,
  absolute_error_unclipped = abs(actual - predicted_unclipped),
  absolute_error_clipped = abs(actual - predicted_clipped),
  stringsAsFactors = FALSE
)
write_csv(predictions, "mixed_effects_final_test_predictions.csv")

per_governorate_metrics <- do.call(
  rbind,
  lapply(split(seq_len(nrow(predictions)), predictions$governorate), function(index) {
    z <- predictions[index, , drop = FALSE]
    sst <- sum((z$actual_production_1000_qt - mean(z$actual_production_1000_qt))^2)
    data.frame(
      governorate = z$governorate[1],
      region = z$region[1],
      N_Test = nrow(z),
      R2_Unclipped = if (sst > 0) 1 - sum(z$residual_unclipped^2) / sst else NA_real_,
      RMSE_Unclipped = sqrt(mean(z$residual_unclipped^2)),
      MAE_Unclipped = mean(abs(z$residual_unclipped)),
      Bias_Unclipped = mean(predicted_unclipped[index] - actual[index]),
      RMSE_Clipped = sqrt(mean(z$residual_clipped^2)),
      MAE_Clipped = mean(abs(z$residual_clipped)),
      Baseline_RMSE = sqrt(mean((z$actual_production_1000_qt - z$baseline_governorate_mean)^2)),
      stringsAsFactors = FALSE
    )
  })
)
per_governorate_metrics <- per_governorate_metrics[
  order(per_governorate_metrics$RMSE_Unclipped),
  , drop = FALSE
]
write_csv(per_governorate_metrics, "mixed_effects_per_governorate_test_metrics.csv")

coefficient_matrix <- summary(selected_model)$tTable
coefficient_export <- data.frame(
  Term = rownames(coefficient_matrix),
  Estimate = coefficient_matrix[, "Value"],
  Std_Error = coefficient_matrix[, "Std.Error"],
  DF = coefficient_matrix[, "DF"],
  T_Value = coefficient_matrix[, "t-value"],
  P_Value = coefficient_matrix[, "p-value"],
  Response_Transformation = selected_candidate$Response_Transformation,
  Inference_Use = "Forecasting output only; not the final H1 inference procedure",
  stringsAsFactors = FALSE,
  row.names = NULL
)
write_csv(coefficient_export, "selected_model_fixed_effects.csv")

random_effect_export <- selected_boundary_info
write_csv(random_effect_export, "selected_model_random_effect_diagnostics.csv")
write_csv(selected_scaling_stats$sm, "soil_moisture_training_scaling.csv")
write_csv(selected_scaling_stats$y, "governorate_production_training_scaling.csv")

random_effects <- nlme::ranef(selected_model)
random_effects_export <- data.frame(
  governorate = rownames(random_effects),
  random_effects,
  check.names = FALSE,
  row.names = NULL
)
write_csv(random_effects_export, "selected_model_governorate_random_effects.csv")

correlation_coefficients <- if (is.null(selected_model$modelStruct$corStruct)) {
  data.frame(Structure = "Independent", Parameter = NA_character_, Estimate = NA_real_)
} else {
  correlation_values <- stats::coef(selected_model$modelStruct$corStruct, unconstrained = FALSE)
  correlation_names <- names(correlation_values)
  if (is.null(correlation_names) || any(!nzchar(correlation_names))) {
    correlation_names <- paste0("AR1_parameter_", seq_along(correlation_values))
  }
  data.frame(
    Structure = "AR1",
    Parameter = correlation_names,
    Estimate = as.numeric(correlation_values),
    row.names = NULL
  )
}
write_csv(correlation_coefficients, "selected_model_residual_correlation.csv")

smearing_factor <- selected_smearing_factor
write_csv(
  data.frame(
    Response_Transformation = selected_candidate$Response_Transformation,
    Duan_Smearing_Factor = smearing_factor,
    Prediction_Level = "Governorate-specific conditional prediction (level 1)",
    Primary_Test_Metric_Uses_Unclipped_Predictions = TRUE
  ),
  "selected_model_prediction_settings.csv"
)

train_fitted <- as.numeric(stats::fitted(selected_model, level = 1))
train_normalized_residuals <- as.numeric(stats::residuals(selected_model, level = 1, type = "normalized"))
if (length(train_fitted) != nrow(selected_train) || length(train_normalized_residuals) != nrow(selected_train)) {
  stop("Could not align final-model fitted values and residuals with the training observations.")
}
train_diagnostics <- data.frame(
  harvest_year = selected_train$harvest_year,
  governorate = as.character(selected_train$governorate),
  region = as.character(selected_train$region),
  observed_model_scale = selected_train$y_model,
  fitted_model_scale = train_fitted,
  normalized_residual = train_normalized_residuals,
  stringsAsFactors = FALSE
)
write_csv(train_diagnostics, "selected_model_training_residual_diagnostics.csv")

lag1_residual_diagnostics <- do.call(
  rbind,
  lapply(split(train_diagnostics, train_diagnostics$governorate), function(z) {
    z <- z[order(z$harvest_year), , drop = FALSE]
    if (nrow(z) < 2L) {
      return(data.frame(
        governorate = z$governorate[1],
        Consecutive_Year_Pairs = 0L,
        Residual_Lag1_Correlation = NA_real_
      ))
    }
    consecutive <- diff(z$harvest_year) == 1L
    left <- z$normalized_residual[-nrow(z)][consecutive]
    right <- z$normalized_residual[-1L][consecutive]
    correlation <- if (length(left) >= 3L && stats::sd(left) > 0 && stats::sd(right) > 0) {
      stats::cor(left, right)
    } else {
      NA_real_
    }
    data.frame(
      governorate = z$governorate[1],
      Consecutive_Year_Pairs = length(left),
      Residual_Lag1_Correlation = correlation
    )
  })
)
write_csv(lag1_residual_diagnostics, "selected_model_residual_lag1_diagnostics.csv")

capture_text(summary(selected_model), "selected_model_summary.txt")
capture_text(nlme::VarCorr(selected_model), "selected_model_variance_components.txt")
capture_text(selected_formula, "selected_model_formula.txt")
capture_text(selected_fit_info$warnings, "selected_model_fit_warnings.txt")
acf_result <- tryCatch(
  nlme::ACF(selected_model, resType = "normalized"),
  error = function(e) e
)
if (inherits(acf_result, "error")) {
  writeLines(
    paste("ACF diagnostic could not be calculated:", conditionMessage(acf_result)),
    con = file.path(OUTPUT_DIR, "selected_model_residual_acf.txt"),
    useBytes = TRUE
  )
} else {
  write_csv(as.data.frame(acf_result), "selected_model_residual_acf.csv")
  capture_text(acf_result, "selected_model_residual_acf.txt")
}
saveRDS(selected_model, file.path(OUTPUT_DIR, "mixed_effects_optimized_durum_wheat_model.rds"))
saveRDS(selected_model, file.path(OUTPUT_DIR, "mixed_effects_selected_forecasting_model.rds"))
saveRDS(
  list(
    model = selected_model,
    formula = selected_formula,
    candidate = selected_candidate,
    response_transformation = selected_candidate$Response_Transformation,
    scaling_parameters = selected_scaling_stats,
    governorate_levels = governorate_levels,
    region_levels = c("North", "Center_South"),
    north_governorates = north_governorates,
    center_south_governorates = center_south_governorates,
    year_origin = 2006L,
    prediction_level = 1L,
    smearing_factor = smearing_factor,
    selection_metrics = selected_rank_row,
    test_metrics = test_performance
  ),
  file.path(OUTPUT_DIR, "mixed_effects_forecast_bundle.rds")
)

png(file.path(OUTPUT_DIR, "test_actual_vs_predicted.png"), width = 1300, height = 1100, res = 160)
plot(
  actual,
  predicted_unclipped,
  xlab = "Observed production (1000 quintals)",
  ylab = "Predicted production (1000 quintals)",
  main = "Selected mixed-effects model: observed vs predicted",
  pch = 19
)
abline(a = 0, b = 1, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "test_residuals_by_year.png"), width = 1300, height = 900, res = 160)
plot(
  predictions$harvest_year,
  predictions$residual_unclipped,
  xlab = "Harvest year",
  ylab = "Observed minus predicted production",
  main = "Test residuals by harvest year",
  pch = 19
)
abline(h = 0, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "training_normalized_residuals_diagnostics.png"), width = 1500, height = 1100, res = 160)
old_par <- par(mfrow = c(2, 2))
plot(
  train_diagnostics$fitted_model_scale,
  train_diagnostics$normalized_residual,
  xlab = "Fitted values (model scale)",
  ylab = "Normalized residuals",
  main = "Residuals versus fitted values",
  pch = 19
)
abline(h = 0, lty = 2)
qqnorm(train_diagnostics$normalized_residual, main = "Normal Q-Q plot")
qqline(train_diagnostics$normalized_residual)
plot(
  train_diagnostics$harvest_year,
  train_diagnostics$normalized_residual,
  xlab = "Harvest year",
  ylab = "Normalized residuals",
  main = "Normalized residuals by year",
  pch = 19
)
abline(h = 0, lty = 2)
boxplot(
  normalized_residual ~ governorate,
  data = train_diagnostics,
  las = 2,
  cex.axis = 0.65,
  xlab = "Governorate",
  ylab = "Normalized residuals",
  main = "Normalized residuals by governorate"
)
par(old_par)
dev.off()

selected_model_summary <- data.frame(
  Selected_Candidate_ID = selected_candidate$Candidate_ID,
  Selected_Model_Name = selected_candidate$Model_Name,
  CV_Rank = selected_rank_row$CV_Rank,
  CV_R2 = selected_rank_row$CV_R2,
  CV_RMSE = selected_rank_row$CV_RMSE,
  CV_MAE = selected_rank_row$CV_MAE,
  Fixed_Specification = selected_candidate$Fixed_Specification,
  Random_Structure = selected_candidate$Random_Structure,
  Residual_Correlation = selected_candidate$Residual_Correlation,
  Response_Transformation = selected_candidate$Response_Transformation,
  Final_Fit_Method = "REML",
  Final_Model_AIC = stats::AIC(selected_model),
  Final_Model_BIC = stats::BIC(selected_model),
  Final_Model_LogLik = as.numeric(stats::logLik(selected_model)),
  Final_Fit_Warning_Count = length(selected_fit_info$warnings),
  Boundary_Fit = FALSE,
  Boundary_SD_Ratio_Threshold = BOUNDARY_SD_RATIO,
  Training_N = nrow(train_data_raw),
  Test_N = nrow(test_data_raw),
  Test_R2_Unclipped = test_performance$R2[test_performance$Prediction == "Selected mixed-effects model - unclipped"],
  Test_RMSE_Unclipped = test_performance$RMSE[test_performance$Prediction == "Selected mixed-effects model - unclipped"],
  Test_MAE_Unclipped = test_performance$MAE[test_performance$Prediction == "Selected mixed-effects model - unclipped"],
  stringsAsFactors = FALSE
)
write_csv(selected_model_summary, "selected_model_selection_summary.csv")

cat("\n========== REVISED NLME MIXED-EFFECTS FORECASTING COMPLETE ==========\n")
cat("Training observations:", nrow(train_data_raw), "\n")
cat("Test observations:", nrow(test_data_raw), "\n")
cat("Rolling-origin validation observations:", nrow(baseline_cv_export), "\n")
cat("Eligible candidate models:", nrow(eligible_cv), "of", nrow(candidates), "\n")
cat("Selected candidate:", selected_candidate$Model_Name, "\n")
cat("Selected fixed-effects formula:\n")
print(selected_formula)
cat("Random structure:", selected_candidate$Random_Structure, "\n")
cat("Residual correlation:", selected_candidate$Residual_Correlation, "\n")
cat("Response transformation:", selected_candidate$Response_Transformation, "\n")
cat("\nTop rolling-origin CV candidates:\n")
print(
  top_three[, c(
    "CV_Rank", "Model_Name", "CV_R2", "CV_RMSE", "CV_MAE",
    "Average_Rank", "Boundary_Observed", "Warning_Count"
  )],
  row.names = FALSE
)
cat("\nFinal 2013-2018 test performance:\n")
print(test_performance, row.names = FALSE)
cat("\nResidual lag-1 diagnostics by governorate:\n")
print(lag1_residual_diagnostics, row.names = FALSE)
cat("\nOutputs saved to:\n", OUTPUT_DIR, "\n")
