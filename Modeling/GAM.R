options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_DIR <- "/home/badisz/Downloads/Modeling"
DATA_PATH <- file.path(PROJECT_DIR, "model_table.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "GAM_Optimized_Output")
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
GAMMA <- 1.4
SMOOTH_K <- 4L
BOUNDARY_SD_RATIO <- 0.001

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("The mgcv package is required. Install it with install.packages('mgcv').")
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

read_data <- utils::read.csv(
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

missing_cols <- setdiff(required_cols, names(read_data))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
if (nrow(read_data) != EXPECTED_N) {
  stop("Expected ", EXPECTED_N, " observations; found ", nrow(read_data), ".")
}
if (anyNA(read_data[, required_cols, drop = FALSE])) {
  stop("Missing values were found in required columns.")
}

numeric_cols <- c(
  "harvest_year", "production_1000_qt", "rain_jan_mar_mm",
  "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)

if (!all(vapply(read_data[, numeric_cols, drop = FALSE], is.numeric, logical(1)))) {
  stop("One or more required numeric columns were imported as non-numeric.")
}
if (!all(vapply(read_data[, numeric_cols, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
  stop("Non-finite values were found in required numeric columns.")
}
if (!all(read_data$crop == "durum_wheat")) {
  stop("The dataset contains crop categories other than durum_wheat.")
}
if (any(read_data$production_1000_qt < 0)) {
  stop("Negative production values were found; inspect the target before fitting.")
}
if (sum(read_data$production_1000_qt == 0) != EXPECTED_ZERO_PRODUCTION) {
  stop(
    "Expected ", EXPECTED_ZERO_PRODUCTION,
    " zero-production observations; found ", sum(read_data$production_1000_qt == 0), "."
  )
}
if (any(read_data$harvest_year != as.integer(read_data$harvest_year))) {
  stop("harvest_year contains non-integer values.")
}
read_data$harvest_year <- as.integer(read_data$harvest_year)

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

if (anyNA(read_data$governorate) || any(!nzchar(trimws(read_data$governorate)))) {
  stop("Empty or missing governorate names were found.")
}
read_data$governorate <- trimws(read_data$governorate)
observed_governorates <- sort(unique(read_data$governorate))
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

if (any(duplicated(read_data[, c("governorate", "harvest_year")])) ) {
  duplicate_rows <- read_data[
    duplicated(read_data[, c("governorate", "harvest_year")]) |
      duplicated(read_data[, c("governorate", "harvest_year")], fromLast = TRUE),
    c("governorate", "harvest_year"),
    drop = FALSE
  ]
  write_csv(duplicate_rows, "duplicate_governorate_year_records.csv")
  stop("Duplicate governorate-year records were found. Inspect duplicate_governorate_year_records.csv.")
}

if (any(read_data$harvest_year < TRAIN_START_YEAR | read_data$harvest_year > TEST_END_YEAR)) {
  stop("Observations outside the 2000-2018 study period were found.")
}
if (!setequal(sort(unique(read_data$harvest_year)), seq.int(TRAIN_START_YEAR, TEST_END_YEAR))) {
  stop("The dataset does not contain every harvest year from 2000 through 2018.")
}

missing_region_map <- setdiff(observed_governorates, c(north_governorates, center_south_governorates))
if (length(missing_region_map)) {
  stop("Unmapped governorates: ", paste(missing_region_map, collapse = ", "))
}

governorate_levels <- expected_governorates
read_data$governorate <- factor(read_data$governorate, levels = governorate_levels)
read_data$region <- factor(
  ifelse(as.character(read_data$governorate) %in% north_governorates, "North", "Center_South"),
  levels = c("North", "Center_South")
)
read_data <- read_data[order(read_data$governorate, read_data$harvest_year), , drop = FALSE]
rownames(read_data) <- NULL

coverage_matrix <- table(
  factor(read_data$governorate, levels = governorate_levels),
  factor(read_data$harvest_year, levels = seq.int(TRAIN_START_YEAR, TEST_END_YEAR))
)
coverage_table <- data.frame(
  governorate = rownames(coverage_matrix),
  Production_Years = rowSums(coverage_matrix > 0),
  Missing_Years = apply(coverage_matrix, 1, function(x) {
    missing_years <- colnames(coverage_matrix)[x == 0]
    if (length(missing_years)) paste(missing_years, collapse = ";") else ""
  }),
  Training_Years = rowSums(coverage_matrix[, as.character(TRAIN_START_YEAR:TRAIN_END_YEAR), drop = FALSE] > 0),
  Test_Years = rowSums(coverage_matrix[, as.character(TEST_START_YEAR:TEST_END_YEAR), drop = FALSE] > 0),
  row.names = NULL
)
write_csv(coverage_table, "gam_governorate_production_coverage.csv")

train_data_raw <- read_data[
  read_data$harvest_year >= TRAIN_START_YEAR & read_data$harvest_year <= TRAIN_END_YEAR,
  , drop = FALSE
]
test_data_raw <- read_data[
  read_data$harvest_year >= TEST_START_YEAR & read_data$harvest_year <= TEST_END_YEAR,
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
  First_Year = c(min(read_data$harvest_year), min(train_data_raw$harvest_year), min(test_data_raw$harvest_year)),
  Last_Year = c(max(read_data$harvest_year), max(train_data_raw$harvest_year), max(test_data_raw$harvest_year)),
  Observations = c(nrow(read_data), nrow(train_data_raw), nrow(test_data_raw)),
  Governorates = c(
    length(unique(read_data$governorate)),
    length(unique(train_data_raw$governorate)),
    length(unique(test_data_raw$governorate))
  ),
  Zero_Production = c(
    sum(read_data$production_1000_qt == 0),
    sum(train_data_raw$production_1000_qt == 0),
    sum(test_data_raw$production_1000_qt == 0)
  )
)
write_csv(split_audit, "gam_temporal_split_audit.csv")

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
  if (any(!is.finite(y_stats$production_mean)) || any(!is.finite(y_stats$production_sd)) || any(y_stats$production_sd <= 0)) {
    stop("Invalid training-period production scaling parameters.")
  }
  list(sm = sm_stats, y = y_stats)
}

apply_features <- function(newdata, stats_obj, response_mode, include_response = TRUE) {
  d <- newdata
  sm_index <- match(as.character(d$governorate), stats_obj$sm$governorate)
  if (anyNA(sm_index)) {
    stop("A governorate has no training-period soil-moisture scaling parameters.")
  }

  d$sm_z <- (d$sm_oct_nov - stats_obj$sm$sm_mean[sm_index]) / stats_obj$sm$sm_sd[sm_index]
  d$year_c <- (d$harvest_year - 2006) / 10
  d$governorate <- factor(as.character(d$governorate), levels = governorate_levels)
  d$region <- factor(as.character(d$region), levels = c("North", "Center_South"))

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
      d$y_model <- (d$production_1000_qt - stats_obj$y$production_mean[y_index]) /
        stats_obj$y$production_sd[y_index]
    } else {
      stop("Unknown response transformation: ", response_mode)
    }
  } else {
    d$y_model <- NA_real_
  }

  if (any(!is.finite(d$sm_z)) || any(!is.finite(d$year_c))) {
    stop("Non-finite derived predictor values were generated.")
  }
  if (anyNA(d$governorate) || anyNA(d$region)) {
    stop("Factor conversion generated missing governorate or region values.")
  }
  d
}

make_ar_start <- function(model_data) {
  model_data <- model_data[order(model_data$governorate, model_data$harvest_year), , drop = FALSE]
  starts <- rep(FALSE, nrow(model_data))
  if (nrow(model_data)) {
    starts[1] <- TRUE
    if (nrow(model_data) > 1L) {
      same_governorate <- as.character(model_data$governorate[-1]) ==
        as.character(model_data$governorate[-nrow(model_data)])
      consecutive_year <- diff(model_data$harvest_year) == 1L
      starts[-1] <- !same_governorate | !consecutive_year
    }
  }
  starts
}

estimate_pooled_ar1 <- function(model, model_data) {
  residual_values <- as.numeric(stats::residuals(model, type = "response"))
  if (length(residual_values) != nrow(model_data) || any(!is.finite(residual_values))) {
    stop("Could not estimate the AR(1) parameter from finite training residuals.")
  }

  residual_data <- data.frame(
    governorate = as.character(model_data$governorate),
    harvest_year = model_data$harvest_year,
    residual = residual_values
  )
  residual_data <- residual_data[order(residual_data$governorate, residual_data$harvest_year), , drop = FALSE]
  groups <- split(residual_data, residual_data$governorate, drop = TRUE)
  previous_values <- numeric(0)
  current_values <- numeric(0)

  for (group_data in groups) {
    if (nrow(group_data) < 2L) next
    group_data <- group_data[order(group_data$harvest_year), , drop = FALSE]
    pair_index <- which(diff(group_data$harvest_year) == 1L)
    if (length(pair_index)) {
      previous_values <- c(previous_values, group_data$residual[pair_index])
      current_values <- c(current_values, group_data$residual[pair_index + 1L])
    }
  }

  if (length(previous_values) < 10L || stats::sd(previous_values) == 0 || stats::sd(current_values) == 0) {
    return(list(rho = 0, pairs = length(previous_values), estimated = FALSE))
  }

  rho <- sum(previous_values * current_values) / sum(previous_values^2)
  if (!is.finite(rho)) {
    return(list(rho = 0, pairs = length(previous_values), estimated = FALSE))
  }
  rho <- max(min(rho, 0.80), -0.80)
  list(rho = rho, pairs = length(previous_values), estimated = TRUE)
}

fixed_specs <- list(
  regional_full_smooth = paste(
    "region + s(governorate, bs = 're') +",
    "s(temp_dec_c, by = region, k = 4, bs = 'tp', id = 1) +",
    "s(sm_z, by = region, k = 4, bs = 'tp', id = 2) +",
    "s(rain_jan_mar_mm, k = 5, bs = 'tp') +",
    "s(ndvi_mean, by = region, k = 4, bs = 'tp', id = 3) +",
    "s(vhi, by = region, k = 4, bs = 'tp', id = 4) +",
    "s(year_c, k = 4, bs = 'tp')"
  ),
  regional_climate_shared_vegetation = paste(
    "region + s(governorate, bs = 're') +",
    "s(temp_dec_c, by = region, k = 4, bs = 'tp', id = 1) +",
    "s(sm_z, by = region, k = 4, bs = 'tp', id = 2) +",
    "s(rain_jan_mar_mm, k = 5, bs = 'tp') +",
    "s(ndvi_mean, k = 5, bs = 'tp') +",
    "s(vhi, k = 5, bs = 'tp') +",
    "s(year_c, k = 4, bs = 'tp')"
  ),
  regional_climate_linear_year = paste(
    "region + s(governorate, bs = 're') +",
    "s(temp_dec_c, by = region, k = 4, bs = 'tp', id = 1) +",
    "s(sm_z, by = region, k = 4, bs = 'tp', id = 2) +",
    "s(rain_jan_mar_mm, k = 5, bs = 'tp') +",
    "s(ndvi_mean, k = 5, bs = 'tp') +",
    "s(vhi, k = 5, bs = 'tp') + year_c"
  )
)

correlation_specs <- c("independent", "AR1")
response_modes <- c("raw", "log1p", "governorate_z")

candidates <- expand.grid(
  Fixed_Specification = names(fixed_specs),
  Residual_Correlation = correlation_specs,
  Response_Transformation = response_modes,
  stringsAsFactors = FALSE
)
candidates$Candidate_ID <- sprintf("G%02d", seq_len(nrow(candidates)))
candidates$Model_Name <- paste(
  candidates$Fixed_Specification,
  candidates$Residual_Correlation,
  candidates$Response_Transformation,
  sep = "__"
)

make_gam_formula <- function(candidate) {
  stats::as.formula(paste("y_model ~", fixed_specs[[candidate$Fixed_Specification]]))
}

fit_gam_candidate <- function(formula, model_data, correlation_name) {
  model_data <- model_data[order(model_data$governorate, model_data$harvest_year), , drop = FALSE]
  model_data$governorate <- factor(as.character(model_data$governorate), levels = governorate_levels)
  model_data$region <- factor(as.character(model_data$region), levels = c("North", "Center_South"))
  if (anyNA(model_data$governorate) || anyNA(model_data$region)) {
    stop("The GAM training data contain invalid factor levels.")
  }

  warning_messages <- character(0)
  capture_warnings <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }

  fit_args <- list(
    formula = formula,
    data = model_data,
    family = stats::gaussian(link = "identity"),
    method = "fREML",
    select = TRUE,
    gamma = GAMMA,
    na.action = stats::na.fail,
    control = mgcv::gam.control(maxit = 200, epsilon = 1e-7)
  )

  preliminary_model <- capture_warnings(do.call(mgcv::bam, fit_args))
  rho_info <- list(rho = 0, pairs = 0L, estimated = FALSE)
  final_model <- preliminary_model

  if (correlation_name == "AR1") {
    rho_info <- estimate_pooled_ar1(preliminary_model, model_data)
    fit_args$rho <- rho_info$rho
    fit_args$AR.start <- make_ar_start(model_data)
    final_model <- capture_warnings(do.call(mgcv::bam, fit_args))
  } else if (correlation_name != "independent") {
    stop("Unknown residual correlation structure: ", correlation_name)
  }

  critical_warning <- grepl(
    "converg|iteration limit|step failure|non-finite|singular fit|rank deficient|not positive definite|false convergence|failed to converge",
    warning_messages,
    ignore.case = TRUE
  )
  if (any(critical_warning)) {
    stop("Critical GAM warning: ", paste(unique(warning_messages[critical_warning]), collapse = " | "))
  }

  coefficients <- stats::coef(final_model)
  if (!length(coefficients) || any(!is.finite(coefficients))) {
    stop("The fitted GAM contains missing or non-finite coefficients.")
  }
  if (!is.null(final_model$converged) && !isTRUE(final_model$converged)) {
    stop("The GAM optimizer did not report convergence.")
  }

  list(
    model = final_model,
    warnings = unique(warning_messages),
    rho = rho_info$rho,
    rho_pairs = rho_info$pairs,
    rho_estimated = rho_info$estimated
  )
}

get_gov_random_effect_diagnostics <- function(model) {
  vc_object <- tryCatch(
    suppressWarnings(mgcv::gam.vcomp(model, rescale = TRUE)),
    error = function(e) NULL
  )
  vc_table <- NULL
  if (is.list(vc_object) && !is.null(vc_object$vc)) {
    vc_table <- vc_object$vc
  } else if (!is.null(vc_object)) {
    vc_table <- vc_object
  }

  random_sd <- NA_real_
  if (!is.null(dim(vc_table)) && !is.null(rownames(vc_table))) {
    random_index <- grep("governorate", rownames(vc_table), ignore.case = TRUE)
    if (length(random_index)) {
      random_sd <- suppressWarnings(as.numeric(vc_table[random_index[1], 1]))
    }
  } else if (is.numeric(vc_table) && !is.null(names(vc_table))) {
    random_index <- grep("governorate", names(vc_table), ignore.case = TRUE)
    if (length(random_index)) {
      random_sd <- as.numeric(vc_table[random_index[1]])
    }
  }

  residual_sd <- if (!is.null(model$sig2) && is.finite(model$sig2) && model$sig2 > 0) sqrt(model$sig2) else NA_real_
  ratio <- if (is.finite(random_sd) && is.finite(residual_sd) && residual_sd > 0) random_sd / residual_sd else NA_real_
  boundary <- if (is.finite(ratio)) ratio < BOUNDARY_SD_RATIO else NA

  data.frame(
    Governorate_Random_Effect_SD = random_sd,
    Residual_SD = residual_sd,
    Random_to_Residual_SD_Ratio = ratio,
    Random_Effect_Boundary_Flag = boundary
  )
}

get_smearing_factor <- function(model, response_mode) {
  if (response_mode != "log1p") return(1)
  residual_values <- as.numeric(stats::residuals(model, type = "response"))
  smearing <- mean(exp(residual_values), na.rm = TRUE)
  if (!is.finite(smearing) || smearing <= 0) {
    stop("Could not calculate a finite Duan smearing factor for log1p predictions.")
  }
  smearing
}

back_transform_predictions <- function(model, newdata, response_mode, target_stats) {
  eta <- as.numeric(stats::predict(model, newdata = newdata, type = "response"))
  if (length(eta) != nrow(newdata) || any(!is.finite(eta))) {
    stop("The GAM returned missing or non-finite predictions.")
  }

  if (response_mode == "raw") return(eta)
  if (response_mode == "log1p") {
    prediction <- exp(eta) * get_smearing_factor(model, response_mode) - 1
    if (any(!is.finite(prediction))) stop("The log1p retransformation generated non-finite predictions.")
    return(prediction)
  }
  if (response_mode == "governorate_z") {
    y_index <- match(as.character(newdata$governorate), target_stats$governorate)
    if (anyNA(y_index)) stop("Could not match governorates to training-period production scaling parameters.")
    prediction <- eta * target_stats$production_sd[y_index] + target_stats$production_mean[y_index]
    if (any(!is.finite(prediction))) stop("The governorate-standardized retransformation generated non-finite predictions.")
    return(prediction)
  }
  stop("Unknown response transformation: ", response_mode)
}

calculate_metrics <- function(actual, predicted, label) {
  if (length(actual) != length(predicted) || !length(actual)) {
    stop("Actual and predicted vectors must have the same positive length.")
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
    stringsAsFactors = FALSE
  )
}

folds <- list(
  list(name = "2000-2006_to_2007-2008", train_end = 2006L, valid_start = 2007L, valid_end = 2008L),
  list(name = "2000-2008_to_2009-2010", train_end = 2008L, valid_start = 2009L, valid_end = 2010L),
  list(name = "2000-2010_to_2011-2012", train_end = 2010L, valid_start = 2011L, valid_end = 2012L)
)

cv_predictions_by_candidate <- vector("list", nrow(candidates))
candidate_summary_rows <- vector("list", nrow(candidates))
fold_diagnostic_rows <- list()
baseline_cv_predictions <- list()

for (fold_index in seq_along(folds)) {
  fold <- folds[[fold_index]]
  fold_train <- train_data_raw[train_data_raw$harvest_year <= fold$train_end, , drop = FALSE]
  fold_valid <- train_data_raw[
    train_data_raw$harvest_year >= fold$valid_start &
      train_data_raw$harvest_year <= fold$valid_end,
    , drop = FALSE
  ]
  if (!nrow(fold_train) || !nrow(fold_valid)) {
    stop("A rolling-origin fold contains no training or validation observations: ", fold$name)
  }
  if (!setequal(as.character(fold_train$governorate), as.character(fold_valid$governorate))) {
    stop("A rolling-origin fold does not contain all governorates in both periods: ", fold$name)
  }
  fold_stats <- fit_scaling_stats(fold_train)
  fold_valid_features <- apply_features(fold_valid, fold_stats, "raw", include_response = FALSE)
  fold_baseline_index <- match(as.character(fold_valid$governorate), fold_stats$y$governorate)
  if (anyNA(fold_baseline_index)) stop("Missing governorate training mean in fold: ", fold$name)
  baseline_cv_predictions[[fold_index]] <- data.frame(
    Fold = fold$name,
    harvest_year = fold_valid$harvest_year,
    governorate = as.character(fold_valid$governorate),
    actual = fold_valid$production_1000_qt,
    predicted = fold_stats$y$production_mean[fold_baseline_index]
  )
}

for (candidate_index in seq_len(nrow(candidates))) {
  candidate <- candidates[candidate_index, , drop = FALSE]
  formula_i <- make_gam_formula(candidate)
  candidate_predictions <- list()
  candidate_failed <- FALSE
  failure_reason <- ""
  candidate_warnings <- character(0)
  candidate_boundary <- FALSE
  boundary_unknown <- FALSE

  for (fold_index in seq_along(folds)) {
    fold <- folds[[fold_index]]
    fold_train_raw <- train_data_raw[train_data_raw$harvest_year <= fold$train_end, , drop = FALSE]
    fold_valid_raw <- train_data_raw[
      train_data_raw$harvest_year >= fold$valid_start &
        train_data_raw$harvest_year <= fold$valid_end,
      , drop = FALSE
    ]

    fold_result <- tryCatch({
      fold_stats <- fit_scaling_stats(fold_train_raw)
      fold_train <- apply_features(
        fold_train_raw,
        fold_stats,
        candidate$Response_Transformation,
        include_response = TRUE
      )
      fold_valid <- apply_features(
        fold_valid_raw,
        fold_stats,
        candidate$Response_Transformation,
        include_response = FALSE
      )
      fit_info <- fit_gam_candidate(formula_i, fold_train, candidate$Residual_Correlation)
      random_diagnostics <- get_gov_random_effect_diagnostics(fit_info$model)
      boundary_flag <- isTRUE(random_diagnostics$Random_Effect_Boundary_Flag[1])
      if (is.na(random_diagnostics$Random_Effect_Boundary_Flag[1])) boundary_unknown <- TRUE
      prediction <- back_transform_predictions(
        fit_info$model,
        fold_valid,
        candidate$Response_Transformation,
        fold_stats$y
      )
      if (length(prediction) != nrow(fold_valid_raw) || any(!is.finite(prediction))) {
        stop("The model did not produce finite predictions for every validation observation.")
      }
      list(
        predictions = data.frame(
          Candidate_ID = candidate$Candidate_ID,
          Model_Name = candidate$Model_Name,
          Fold = fold$name,
          harvest_year = fold_valid_raw$harvest_year,
          governorate = as.character(fold_valid_raw$governorate),
          region = as.character(fold_valid_raw$region),
          actual = fold_valid_raw$production_1000_qt,
          predicted = prediction
        ),
        warnings = fit_info$warnings,
        rho = fit_info$rho,
        rho_pairs = fit_info$rho_pairs,
        rho_estimated = fit_info$rho_estimated,
        random_diagnostics = random_diagnostics,
        boundary = boundary_flag,
        error = ""
      )
    }, error = function(e) {
      list(error = conditionMessage(e))
    })

    if (nzchar(fold_result$error)) {
      candidate_failed <- TRUE
      failure_reason <- paste0(fold$name, ": ", fold_result$error)
      fold_diagnostic_rows[[length(fold_diagnostic_rows) + 1L]] <- data.frame(
        Candidate_ID = candidate$Candidate_ID,
        Model_Name = candidate$Model_Name,
        Fold = fold$name,
        Fit_Success = FALSE,
        AR1_Rho = NA_real_,
        AR1_Pairs = NA_integer_,
        AR1_Estimated = NA,
        Governorate_RE_SD = NA_real_,
        Residual_SD = NA_real_,
        RE_to_Residual_SD_Ratio = NA_real_,
        RE_Boundary = NA,
        Warning_Count = NA_integer_,
        Warning_Details = "",
        Failure_Reason = fold_result$error,
        stringsAsFactors = FALSE
      )
      break
    }

    candidate_predictions[[fold_index]] <- fold_result$predictions
    candidate_warnings <- c(candidate_warnings, fold_result$warnings)
    candidate_boundary <- candidate_boundary || fold_result$boundary
    fold_diagnostic_rows[[length(fold_diagnostic_rows) + 1L]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fold = fold$name,
      Fit_Success = TRUE,
      AR1_Rho = fold_result$rho,
      AR1_Pairs = fold_result$rho_pairs,
      AR1_Estimated = fold_result$rho_estimated,
      Governorate_RE_SD = fold_result$random_diagnostics$Governorate_Random_Effect_SD[1],
      Residual_SD = fold_result$random_diagnostics$Residual_SD[1],
      RE_to_Residual_SD_Ratio = fold_result$random_diagnostics$Random_to_Residual_SD_Ratio[1],
      RE_Boundary = fold_result$random_diagnostics$Random_Effect_Boundary_Flag[1],
      Warning_Count = length(fold_result$warnings),
      Warning_Details = if (length(fold_result$warnings)) paste(unique(fold_result$warnings), collapse = " | ") else "",
      Failure_Reason = "",
      stringsAsFactors = FALSE
    )
  }

  complete_folds <- !candidate_failed && length(candidate_predictions) == length(folds)
  if (complete_folds) {
    candidate_cv <- do.call(rbind, candidate_predictions)
    cv_predictions_by_candidate[[candidate_index]] <- candidate_cv
    cv_metric <- calculate_metrics(
      candidate_cv$actual,
      candidate_cv$predicted,
      candidate$Model_Name
    )
    candidate_summary_rows[[candidate_index]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fixed_Specification = candidate$Fixed_Specification,
      Residual_Correlation = candidate$Residual_Correlation,
      Response_Transformation = candidate$Response_Transformation,
      CV_N = cv_metric$N,
      CV_R2 = cv_metric$R2,
      CV_RMSE = cv_metric$RMSE,
      CV_MAE = cv_metric$MAE,
      CV_Mean_Error = cv_metric$Mean_Error,
      Successful_Folds = length(candidate_predictions),
      Random_Effect_Boundary_Observed = candidate_boundary,
      Random_Effect_Diagnostic_Unknown = boundary_unknown,
      Warning_Count = length(unique(candidate_warnings)),
      Warning_Details = if (length(candidate_warnings)) paste(unique(candidate_warnings), collapse = " | ") else "",
      CV_Eligible = !candidate_boundary,
      Failure_Reason = "",
      Complexity_Rank = match(candidate$Fixed_Specification, names(fixed_specs)) +
        match(candidate$Residual_Correlation, correlation_specs) - 1L +
        match(candidate$Response_Transformation, response_modes) - 1L,
      stringsAsFactors = FALSE
    )
  } else {
    candidate_summary_rows[[candidate_index]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID,
      Model_Name = candidate$Model_Name,
      Fixed_Specification = candidate$Fixed_Specification,
      Residual_Correlation = candidate$Residual_Correlation,
      Response_Transformation = candidate$Response_Transformation,
      CV_N = NA_real_,
      CV_R2 = NA_real_,
      CV_RMSE = NA_real_,
      CV_MAE = NA_real_,
      CV_Mean_Error = NA_real_,
      Successful_Folds = length(candidate_predictions),
      Random_Effect_Boundary_Observed = candidate_boundary,
      Random_Effect_Diagnostic_Unknown = boundary_unknown,
      Warning_Count = length(unique(candidate_warnings)),
      Warning_Details = if (length(candidate_warnings)) paste(unique(candidate_warnings), collapse = " | ") else "",
      CV_Eligible = FALSE,
      Failure_Reason = failure_reason,
      Complexity_Rank = match(candidate$Fixed_Specification, names(fixed_specs)) +
        match(candidate$Residual_Correlation, correlation_specs) - 1L +
        match(candidate$Response_Transformation, response_modes) - 1L,
      stringsAsFactors = FALSE
    )
  }

  cat(
    "Candidate", candidate_index, "of", nrow(candidates), ":", candidate$Model_Name,
    if (candidate_failed) paste("FAILED:", failure_reason) else if (candidate_boundary) "completed with a boundary governorate effect" else "completed",
    "\n"
  )
}

candidate_cv_summary <- do.call(rbind, candidate_summary_rows)
fold_diagnostics <- if (length(fold_diagnostic_rows)) do.call(rbind, fold_diagnostic_rows) else data.frame()
write_csv(candidate_cv_summary, "gam_all_candidate_cv_results.csv")
write_csv(fold_diagnostics, "gam_fold_fit_diagnostics.csv")

eligible_cv <- candidate_cv_summary[candidate_cv_summary$CV_Eligible, , drop = FALSE]
if (!nrow(eligible_cv)) {
  stop("No candidate completed all validation folds without a boundary governorate random effect. Inspect gam_fold_fit_diagnostics.csv and gam_all_candidate_cv_results.csv.")
}
if (length(unique(eligible_cv$CV_N)) != 1L) {
  stop("Eligible GAM candidates were evaluated on different numbers of validation observations.")
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
write_csv(eligible_cv, "gam_eligible_candidate_cv_rankings.csv")
write_csv(utils::head(eligible_cv, min(3L, nrow(eligible_cv))), "gam_top_three_cv_candidates.csv")

cv_prediction_export <- do.call(
  rbind,
  cv_predictions_by_candidate[!vapply(cv_predictions_by_candidate, is.null, logical(1))]
)
write_csv(cv_prediction_export, "gam_rolling_origin_cv_predictions.csv")

baseline_cv_export <- do.call(rbind, baseline_cv_predictions)
baseline_cv_metrics <- calculate_metrics(
  baseline_cv_export$actual,
  baseline_cv_export$predicted,
  "Governorate-specific training-period mean baseline"
)
write_csv(baseline_cv_export, "gam_rolling_origin_baseline_predictions.csv")
write_csv(baseline_cv_metrics, "gam_rolling_origin_baseline_metrics.csv")

selected_model <- NULL
selected_fit_info <- NULL
selected_candidate <- NULL
selected_formula <- NULL
selected_scaling <- NULL
selected_train <- NULL
selected_test <- NULL
selected_test_predictions <- NULL
selected_random_diagnostics <- NULL
final_refit_attempts <- list()

for (rank_index in seq_len(nrow(eligible_cv))) {
  ranked_row <- eligible_cv[rank_index, , drop = FALSE]
  candidate <- candidates[candidates$Candidate_ID == ranked_row$Candidate_ID, , drop = FALSE]
  formula_i <- make_gam_formula(candidate)

  attempt <- tryCatch({
    scaling <- fit_scaling_stats(train_data_raw)
    final_train_data <- apply_features(
      train_data_raw,
      scaling,
      candidate$Response_Transformation,
      include_response = TRUE
    )
    final_test_data <- apply_features(
      test_data_raw,
      scaling,
      candidate$Response_Transformation,
      include_response = FALSE
    )
    fit_info <- fit_gam_candidate(formula_i, final_train_data, candidate$Residual_Correlation)
    random_diagnostics <- get_gov_random_effect_diagnostics(fit_info$model)
    if (isTRUE(random_diagnostics$Random_Effect_Boundary_Flag[1])) {
      stop("The final refit produced a boundary governorate random-effect variance.")
    }
    predictions <- back_transform_predictions(
      fit_info$model,
      final_test_data,
      candidate$Response_Transformation,
      scaling$y
    )
    if (length(predictions) != nrow(test_data_raw) || any(!is.finite(predictions))) {
      stop("The final GAM did not predict every test observation with a finite value.")
    }
    list(
      success = TRUE,
      model = fit_info$model,
      fit_info = fit_info,
      candidate = candidate,
      formula = formula_i,
      scaling = scaling,
      train = final_train_data,
      test = final_test_data,
      predictions = predictions,
      random_diagnostics = random_diagnostics,
      error = ""
    )
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })

  final_refit_attempts[[length(final_refit_attempts) + 1L]] <- data.frame(
    CV_Rank = ranked_row$CV_Rank,
    Candidate_ID = candidate$Candidate_ID,
    Model_Name = candidate$Model_Name,
    Final_Refit_Success = attempt$success,
    Failure_Reason = attempt$error,
    stringsAsFactors = FALSE
  )

  if (isTRUE(attempt$success)) {
    selected_model <- attempt$model
    selected_fit_info <- attempt$fit_info
    selected_candidate <- attempt$candidate
    selected_formula <- attempt$formula
    selected_scaling <- attempt$scaling
    selected_train <- attempt$train
    selected_test <- attempt$test
    selected_test_predictions <- attempt$predictions
    selected_random_diagnostics <- attempt$random_diagnostics
    break
  }
}

final_refit_attempts <- do.call(rbind, final_refit_attempts)
write_csv(final_refit_attempts, "gam_final_refit_attempts.csv")
if (is.null(selected_model)) {
  stop("No CV-eligible GAM candidate could be refitted successfully on the full training period. Inspect gam_final_refit_attempts.csv.")
}

selected_rank_row <- eligible_cv[eligible_cv$Candidate_ID == selected_candidate$Candidate_ID, , drop = FALSE]
predicted_unclipped <- selected_test_predictions
predicted_clipped <- pmax(predicted_unclipped, 0)
actual <- test_data_raw$production_1000_qt

baseline_index <- match(as.character(test_data_raw$governorate), selected_scaling$y$governorate)
if (anyNA(baseline_index)) stop("Could not match test governorates to training-period production means.")
baseline_prediction <- selected_scaling$y$production_mean[baseline_index]

test_performance <- rbind(
  calculate_metrics(actual, predicted_unclipped, "Selected GAM - unclipped"),
  calculate_metrics(actual, predicted_clipped, "Selected GAM - clipped at zero"),
  calculate_metrics(actual, baseline_prediction, "Governorate-specific training-period mean baseline")
)
write_csv(test_performance, "gam_final_test_performance.csv")

test_predictions <- data.frame(
  harvest_year = test_data_raw$harvest_year,
  governorate = as.character(test_data_raw$governorate),
  region = as.character(test_data_raw$region),
  actual = actual,
  predicted_unclipped = predicted_unclipped,
  predicted_clipped = predicted_clipped,
  residual_unclipped = actual - predicted_unclipped,
  residual_clipped = actual - predicted_clipped,
  absolute_error_unclipped = abs(actual - predicted_unclipped),
  stringsAsFactors = FALSE
)
write_csv(test_predictions, "gam_final_test_predictions.csv")

per_governorate_metrics <- do.call(rbind, lapply(
  split(seq_len(nrow(test_predictions)), test_predictions$governorate),
  function(index) {
    selected <- test_predictions[index, , drop = FALSE]
    data.frame(
      Governorate = selected$governorate[1],
      Region = selected$region[1],
      N_Test = nrow(selected),
      RMSE_Unclipped = sqrt(mean((selected$actual - selected$predicted_unclipped)^2)),
      MAE_Unclipped = mean(abs(selected$actual - selected$predicted_unclipped)),
      RMSE_Clipped = sqrt(mean((selected$actual - selected$predicted_clipped)^2)),
      MAE_Clipped = mean(abs(selected$actual - selected$predicted_clipped)),
      Mean_Error_Unclipped = mean(selected$predicted_unclipped - selected$actual),
      stringsAsFactors = FALSE
    )
  }
))
per_governorate_metrics <- per_governorate_metrics[order(per_governorate_metrics$RMSE_Unclipped), , drop = FALSE]
write_csv(per_governorate_metrics, "gam_per_governorate_test_metrics.csv")

summary_object <- summary(selected_model)
parametric_matrix <- summary_object$p.table
if (!is.null(parametric_matrix) && nrow(parametric_matrix)) {
  parametric_inference <- data.frame(
    Term = rownames(parametric_matrix),
    Estimate = parametric_matrix[, 1],
    Std_Error = parametric_matrix[, 2],
    Statistic = parametric_matrix[, 3],
    P_value = parametric_matrix[, 4],
    CI_Lower_Normal_Approx = parametric_matrix[, 1] - stats::qnorm(0.975) * parametric_matrix[, 2],
    CI_Upper_Normal_Approx = parametric_matrix[, 1] + stats::qnorm(0.975) * parametric_matrix[, 2],
    row.names = NULL
  )
} else {
  parametric_inference <- data.frame()
}
write_csv(parametric_inference, "gam_selected_parametric_summary.csv")

smooth_matrix <- summary_object$s.table
if (!is.null(smooth_matrix) && nrow(smooth_matrix)) {
  smooth_inference <- data.frame(
    Term = rownames(smooth_matrix),
    smooth_matrix,
    check.names = FALSE,
    row.names = NULL
  )
} else {
  smooth_inference <- data.frame()
}
write_csv(smooth_inference, "gam_selected_smooth_summary.csv")

concurvity_matrix <- mgcv::concurvity(selected_model, full = TRUE)
concurvity_output <- data.frame(
  Term = rownames(t(concurvity_matrix)),
  t(concurvity_matrix),
  check.names = FALSE,
  row.names = NULL
)
write_csv(concurvity_output, "gam_selected_concurvity.csv")

variance_component_object <- tryCatch(
  suppressWarnings(mgcv::gam.vcomp(selected_model, rescale = TRUE)),
  error = function(e) e
)
if (inherits(variance_component_object, "error")) {
  variance_component_output <- data.frame(
    Output = paste("gam.vcomp could not be calculated:", conditionMessage(variance_component_object))
  )
} else if (is.list(variance_component_object) && !is.null(variance_component_object$all)) {
  all_components <- as.numeric(variance_component_object$all)
  component_names <- names(variance_component_object$all)
  if (is.null(component_names) || length(component_names) != length(all_components)) {
    component_names <- paste0("Variance_Component_", seq_along(all_components))
  }
  variance_component_output <- data.frame(
    Term = component_names,
    Variance_Component = all_components,
    row.names = NULL
  )
} else if (is.list(variance_component_object) && !is.null(variance_component_object$vc)) {
  component_matrix <- as.matrix(variance_component_object$vc)
  if (is.null(rownames(component_matrix))) {
    rownames(component_matrix) <- paste0("Variance_Component_", seq_len(nrow(component_matrix)))
  }
  variance_component_output <- data.frame(
    Term = rownames(component_matrix),
    component_matrix,
    check.names = FALSE,
    row.names = NULL
  )
} else if (is.list(variance_component_object)) {
  variance_component_output <- data.frame(Output = "gam.vcomp returned a list without a recognized variance-component field.")
} else if (!is.null(dim(variance_component_object))) {
  component_matrix <- as.matrix(variance_component_object)
  if (is.null(rownames(component_matrix))) {
    rownames(component_matrix) <- paste0("Variance_Component_", seq_len(nrow(component_matrix)))
  }
  variance_component_output <- data.frame(
    Term = rownames(component_matrix),
    component_matrix,
    check.names = FALSE,
    row.names = NULL
  )
} else {
  all_components <- as.numeric(variance_component_object)
  component_names <- names(variance_component_object)
  if (is.null(component_names) || length(component_names) != length(all_components)) {
    component_names <- paste0("Variance_Component_", seq_along(all_components))
  }
  variance_component_output <- data.frame(
    Term = component_names,
    Variance_Component = all_components,
    row.names = NULL
  )
}
write_csv(variance_component_output, "gam_selected_variance_components.csv")
write_csv(selected_random_diagnostics, "gam_selected_governorate_random_effect_diagnostics.csv")
write_csv(selected_scaling$sm, "gam_training_soil_moisture_scaling.csv")
write_csv(selected_scaling$y, "gam_training_production_scaling.csv")

train_residuals <- as.numeric(stats::residuals(selected_model, type = "response"))
train_fitted <- as.numeric(stats::fitted(selected_model))
if (length(train_residuals) != nrow(selected_train) || length(train_fitted) != nrow(selected_train)) {
  stop("The selected GAM returned residual or fitted-value vectors with unexpected lengths.")
}

pearson_residuals <- as.numeric(stats::residuals(selected_model, type = "pearson"))
ar1_standardized_residuals <- if (!is.null(selected_model$std.rsd)) {
  as.numeric(selected_model$std.rsd)
} else {
  rep(NA_real_, nrow(selected_train))
}
if (length(pearson_residuals) != nrow(selected_train)) {
  pearson_residuals <- rep(NA_real_, nrow(selected_train))
}
if (length(ar1_standardized_residuals) != nrow(selected_train)) {
  ar1_standardized_residuals <- rep(NA_real_, nrow(selected_train))
}

residual_diagnostics <- data.frame(
  governorate = as.character(selected_train$governorate),
  region = as.character(selected_train$region),
  harvest_year = selected_train$harvest_year,
  fitted = train_fitted,
  residual = train_residuals,
  pearson_residual = pearson_residuals,
  AR1_standardized_residual = ar1_standardized_residuals,
  stringsAsFactors = FALSE
)
residual_diagnostics <- residual_diagnostics[order(residual_diagnostics$governorate, residual_diagnostics$harvest_year), , drop = FALSE]
write_csv(residual_diagnostics, "gam_training_residual_diagnostics.csv")

lag_pairs <- list()
residual_groups <- split(residual_diagnostics, residual_diagnostics$governorate, drop = TRUE)
for (governorate_name in names(residual_groups)) {
  group_data <- residual_groups[[governorate_name]]
  group_data <- group_data[order(group_data$harvest_year), , drop = FALSE]
  pair_index <- which(diff(group_data$harvest_year) == 1L)
  if (length(pair_index)) {
    lag_pairs[[length(lag_pairs) + 1L]] <- data.frame(
      governorate = governorate_name,
      year_previous = group_data$harvest_year[pair_index],
      year_current = group_data$harvest_year[pair_index + 1L],
      residual_previous = group_data$residual[pair_index],
      residual_current = group_data$residual[pair_index + 1L],
      stringsAsFactors = FALSE
    )
  }
}
if (length(lag_pairs)) {
  lag_pair_data <- do.call(rbind, lag_pairs)
  residual_lag1_correlation <- if (nrow(lag_pair_data) >= 3L) {
    stats::cor(lag_pair_data$residual_previous, lag_pair_data$residual_current)
  } else {
    NA_real_
  }
  residual_lag1_slope <- if (nrow(lag_pair_data) >= 3L && sum(lag_pair_data$residual_previous^2) > 0) {
    sum(lag_pair_data$residual_previous * lag_pair_data$residual_current) /
      sum(lag_pair_data$residual_previous^2)
  } else {
    NA_real_
  }
} else {
  lag_pair_data <- data.frame()
  residual_lag1_correlation <- NA_real_
  residual_lag1_slope <- NA_real_
}
write_csv(lag_pair_data, "gam_training_residual_lag1_pairs.csv")
std_resid_values <- residual_diagnostics$AR1_standardized_residual
std_resid_lag1_correlation <- NA_real_
if (length(std_resid_values) == nrow(residual_diagnostics) && all(is.finite(std_resid_values))) {
  std_residual_pairs <- list()
  std_groups <- split(residual_diagnostics, residual_diagnostics$governorate, drop = TRUE)
  for (governorate_name in names(std_groups)) {
    group_data <- std_groups[[governorate_name]]
    group_data <- group_data[order(group_data$harvest_year), , drop = FALSE]
    pair_index <- which(diff(group_data$harvest_year) == 1L)
    if (length(pair_index)) {
      std_residual_pairs[[length(std_residual_pairs) + 1L]] <- data.frame(
        previous = group_data$AR1_standardized_residual[pair_index],
        current = group_data$AR1_standardized_residual[pair_index + 1L]
      )
    }
  }
  if (length(std_residual_pairs)) {
    std_residual_pairs <- do.call(rbind, std_residual_pairs)
    if (nrow(std_residual_pairs) >= 3L && stats::sd(std_residual_pairs$previous) > 0 && stats::sd(std_residual_pairs$current) > 0) {
      std_resid_lag1_correlation <- stats::cor(std_residual_pairs$previous, std_residual_pairs$current)
    }
  }
}

residual_autocorrelation_summary <- data.frame(
  Selected_Residual_Correlation = selected_candidate$Residual_Correlation,
  Selected_AR1_Rho = selected_fit_info$rho,
  AR1_Estimation_Pairs = selected_fit_info$rho_pairs,
  AR1_Estimated_From_Training = selected_fit_info$rho_estimated,
  Raw_Residual_Lag1_Correlation = residual_lag1_correlation,
  Raw_Residual_Lag1_Regression_Slope = residual_lag1_slope,
  AR1_Standardized_Residual_Lag1_Correlation = std_resid_lag1_correlation,
  Consecutive_Year_Pairs = nrow(lag_pair_data)
)
write_csv(residual_autocorrelation_summary, "gam_training_residual_autocorrelation_summary.csv")

k_check_result <- tryCatch(
  mgcv::k.check(selected_model),
  error = function(e) e
)
if (inherits(k_check_result, "error")) {
  k_check_output <- data.frame(Output = conditionMessage(k_check_result))
} else {
  k_check_output <- data.frame(
    Smooth = rownames(k_check_result),
    k_check_result,
    check.names = FALSE,
    row.names = NULL
  )
}
write_csv(k_check_output, "gam_selected_k_check.csv")

selected_model_record <- data.frame(
  Candidate_ID = selected_candidate$Candidate_ID,
  Model_Name = selected_candidate$Model_Name,
  Fixed_Specification = selected_candidate$Fixed_Specification,
  Residual_Correlation = selected_candidate$Residual_Correlation,
  Response_Transformation = selected_candidate$Response_Transformation,
  CV_Rank = selected_rank_row$CV_Rank,
  CV_RMSE = selected_rank_row$CV_RMSE,
  CV_MAE = selected_rank_row$CV_MAE,
  CV_R2 = selected_rank_row$CV_R2,
  AR1_Rho = selected_fit_info$rho,
  AR1_Pairs = selected_fit_info$rho_pairs,
  AR1_Estimated = selected_fit_info$rho_estimated,
  stringsAsFactors = FALSE
)
write_csv(selected_model_record, "gam_selected_model_record.csv")

capture_text(summary(selected_model), "gam_selected_model_summary.txt")
capture_text(selected_formula, "gam_selected_model_formula.txt")
capture_text(selected_model, "gam_selected_model_object.txt")
capture_text(selected_fit_info$warnings, "gam_selected_model_warnings.txt")
saveRDS(selected_model, file.path(OUTPUT_DIR, "gam_forecasting_selected_model.rds"))
saveRDS(
  list(
    model = selected_model,
    formula = selected_formula,
    candidate = selected_candidate,
    scaling = selected_scaling,
    ar1_rho = selected_fit_info$rho,
    ar1_pairs = selected_fit_info$rho_pairs,
    response_transformation = selected_candidate$Response_Transformation
  ),
  file.path(OUTPUT_DIR, "gam_forecasting_complete_object.rds")
)

writeLines(
  c(
    "This is a chronological forecasting model fitted on 2000-2012 and evaluated on 2013-2018.",
    "All candidate formulas estimate region-specific smooth responses for December temperature and standardized soil moisture.",
    "Rolling-origin validation is used only within the training period for candidate selection.",
    "The test-period predictors are observed values; the results are retrospective conditional predictions, not forecasts of future climate predictors.",
    "The parametric and smooth p-values exported from mgcv are model summaries, not governorate-clustered H1 inference.",
    "Use the separate full-period inference workflow for H1 tests, confidence intervals, regional contrasts and multiple-testing correction."
  ),
  file.path(OUTPUT_DIR, "gam_forecasting_scope_and_inference_limits.txt"),
  useBytes = TRUE
)

png(file.path(OUTPUT_DIR, "gam_selected_smooth_effects.png"), width = 2000, height = 1600, res = 180)
par(mfrow = c(2, 2))
plot(selected_model, pages = 1, shade = TRUE, residuals = FALSE, scale = 0)
dev.off()

png(file.path(OUTPUT_DIR, "gam_test_actual_vs_predicted.png"), width = 1300, height = 1100, res = 160)
plot(
  actual,
  predicted_unclipped,
  xlab = "Observed production (1000 quintals)",
  ylab = "Predicted production (1000 quintals)",
  main = "Selected GAM: observed vs predicted production",
  pch = 19
)
abline(a = 0, b = 1, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "gam_test_residuals_by_year.png"), width = 1400, height = 950, res = 160)
plot(
  test_predictions$harvest_year,
  test_predictions$residual_unclipped,
  xlab = "Harvest year",
  ylab = "Observed minus predicted production",
  main = "GAM test residuals by harvest year",
  pch = 19
)
abline(h = 0, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "gam_training_residual_diagnostics.png"), width = 1600, height = 1200, res = 160)
par(mfrow = c(2, 2))
plot(
  train_fitted,
  train_residuals,
  xlab = "Fitted production",
  ylab = "Residual",
  main = "Residuals vs fitted",
  pch = 19
)
abline(h = 0, lty = 2)
qqnorm(train_residuals, main = "Normal Q-Q plot")
qqline(train_residuals)
hist(train_residuals, main = "Training residual distribution", xlab = "Residual")
plot(
  selected_train$harvest_year,
  train_residuals,
  xlab = "Harvest year",
  ylab = "Residual",
  main = "Residuals over training years",
  pch = 19
)
abline(h = 0, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "gam_selected_gam_check.png"), width = 1600, height = 1200, res = 160)
par(mfrow = c(2, 2))
gam_check_text <- capture.output(mgcv::gam.check(selected_model))
dev.off()
writeLines(gam_check_text, file.path(OUTPUT_DIR, "gam_selected_gam_check.txt"), useBytes = TRUE)

cat("\n========== GAM FORECASTING OPTIMIZATION COMPLETE ==========\n")
cat("Training observations:", nrow(train_data_raw), "\n")
cat("Test observations:", nrow(test_data_raw), "\n")
cat("Selected candidate:", selected_candidate$Model_Name, "\n")
cat("Selected formula:\n")
print(selected_formula)
cat("Residual correlation:", selected_candidate$Residual_Correlation, "\n")
cat("Estimated AR(1) rho:", selected_fit_info$rho, "\n")
cat("\nTop rolling-origin CV candidates:\n")
print(
  utils::head(
    eligible_cv[, c("CV_Rank", "Model_Name", "CV_R2", "CV_RMSE", "CV_MAE", "Average_Rank", "Warning_Count")],
    3L
  ),
  row.names = FALSE
)
cat("\nFinal 2013-2018 test performance:\n")
print(test_performance, row.names = FALSE)
cat("\nSelected governorate random-effect diagnostics:\n")
print(selected_random_diagnostics, row.names = FALSE)
cat("\nResidual autocorrelation diagnostics:\n")
print(residual_autocorrelation_summary, row.names = FALSE)
cat("\nOutputs saved to:\n", OUTPUT_DIR, "\n")
