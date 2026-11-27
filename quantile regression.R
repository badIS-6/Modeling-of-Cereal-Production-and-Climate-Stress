options(stringsAsFactors = FALSE, scipen = 999)

PROJECT_DIR <- "/home/badisz/Downloads/Modeling"
DATA_PATH <- file.path(PROJECT_DIR, "model_table.csv")
OUTPUT_DIR <- file.path(PROJECT_DIR, "Quantile_Regression_Optimized_Output")
TRAIN_START_YEAR <- 2000L
TRAIN_END_YEAR <- 2012L
TEST_START_YEAR <- 2013L
TEST_END_YEAR <- 2018L
EXPECTED_N <- 281L
EXPECTED_TRAIN_N <- 192L
EXPECTED_TEST_N <- 89L
EXPECTED_GOVERNORATES <- 15L
EXPECTED_ZERO_PRODUCTION <- 11L
YEAR_ORIGIN <- 2006L
SEED <- 20261126L
TAU_GRID <- c(0.10, 0.25, 0.50, 0.75, 0.90)
PRIMARY_TAU <- 0.50
ALLOW_NONUNIQUE_WARNING <- TRUE

if (!requireNamespace("quantreg", quietly = TRUE)) {
  stop("The quantreg package is required. Install it with install.packages('quantreg').")
}
if (!dir.exists(PROJECT_DIR)) stop("Project directory not found: ", PROJECT_DIR)
if (!file.exists(DATA_PATH)) stop("Dataset not found: ", DATA_PATH)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

write_csv <- function(x, filename) {
  utils::write.csv(
    x, file = file.path(OUTPUT_DIR, filename), row.names = FALSE,
    fileEncoding = "UTF-8", na = ""
  )
}

required_cols <- c(
  "governorate", "crop", "harvest_year", "production_1000_qt",
  "rain_jan_mar_mm", "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)
numeric_cols <- c(
  "harvest_year", "production_1000_qt", "rain_jan_mar_mm",
  "temp_dec_c", "sm_oct_nov", "ndvi_mean", "vhi"
)

data <- utils::read.csv(
  DATA_PATH, fileEncoding = "UTF-8", stringsAsFactors = FALSE,
  check.names = FALSE, na.strings = c("", "NA", "NaN")
)

missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
if (nrow(data) != EXPECTED_N) stop("Expected ", EXPECTED_N, " observations; found ", nrow(data), ".")
if (anyNA(data[, required_cols, drop = FALSE])) stop("Missing values were found in required columns.")
if (!all(vapply(data[, numeric_cols, drop = FALSE], is.numeric, logical(1)))) {
  stop("One or more required numeric columns were imported as non-numeric.")
}
if (!all(vapply(data[, numeric_cols, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
  stop("Non-finite values were found in required numeric columns.")
}
if (!all(data$crop == "durum_wheat")) stop("The dataset contains crop categories other than durum_wheat.")
if (any(data$production_1000_qt < 0)) stop("Negative production values were found.")
if (sum(data$production_1000_qt == 0) != EXPECTED_ZERO_PRODUCTION) {
  stop("Expected ", EXPECTED_ZERO_PRODUCTION, " zero-production observations; found ",
       sum(data$production_1000_qt == 0), ".")
}
if (length(unique(data$governorate)) != EXPECTED_GOVERNORATES) {
  stop("Expected ", EXPECTED_GOVERNORATES, " governorates; found ",
       length(unique(data$governorate)), ".")
}
if (any(data$harvest_year < TRAIN_START_YEAR | data$harvest_year > TEST_END_YEAR)) {
  stop("The dataset contains harvest years outside 2000-2018.")
}
if (anyDuplicated(data[, c("governorate", "harvest_year")]) > 0L) {
  stop("Duplicate governorate-year records were found.")
}

north_governorates <- c(
  "Béja", "Jendouba", "Le Kef", "Siliana",
  "Zaghouan", "Manouba", "Nabeul", "Tunis"
)
center_south_governorates <- c(
  "Kairouan", "Kasserine", "Sidi Bouzid", "Sousse",
  "Gafsa", "Gabès", "Medenine"
)
unmapped <- setdiff(unique(data$governorate), c(north_governorates, center_south_governorates))
if (length(unmapped)) stop("Unmapped governorates: ", paste(unmapped, collapse = ", "))

data$harvest_year <- as.integer(data$harvest_year)
data$region <- factor(
  ifelse(data$governorate %in% north_governorates, "North", "Center_South"),
  levels = c("North", "Center_South")
)
data <- data[order(data$harvest_year, data$governorate), , drop = FALSE]
rownames(data) <- NULL

coverage_rows <- do.call(rbind, lapply(sort(unique(data$governorate)), function(g) {
  years <- sort(unique(data$harvest_year[data$governorate == g]))
  data.frame(
    Governorate = g,
    Observed_Production_Years = length(years),
    Expected_Production_Years = TRAIN_END_YEAR - TRAIN_START_YEAR + 1L +
      TEST_END_YEAR - TEST_START_YEAR + 1L,
    Missing_Production_Years = 19L - length(years),
    Coverage_Percent = 100 * length(years) / 19,
    First_Observed_Year = min(years),
    Last_Observed_Year = max(years)
  )
}))

train_data_raw <- data[data$harvest_year >= TRAIN_START_YEAR &
                         data$harvest_year <= TRAIN_END_YEAR, , drop = FALSE]
test_data_raw <- data[data$harvest_year >= TEST_START_YEAR &
                        data$harvest_year <= TEST_END_YEAR, , drop = FALSE]
if (nrow(train_data_raw) != EXPECTED_TRAIN_N || nrow(test_data_raw) != EXPECTED_TEST_N) {
  stop("Expected temporal split ", EXPECTED_TRAIN_N, "/", EXPECTED_TEST_N,
       "; found ", nrow(train_data_raw), "/", nrow(test_data_raw), ".")
}
if (!setequal(unique(train_data_raw$governorate), unique(test_data_raw$governorate))) {
  stop("Training and test periods do not contain the same governorate set.")
}

fit_preprocessor <- function(train) {
  sm_groups <- split(train$sm_oct_nov, as.character(train$governorate))
  sm_stats <- do.call(rbind, lapply(names(sm_groups), function(g) {
    x <- sm_groups[[g]]
    data.frame(governorate = g, sm_train_mean = mean(x),
               sm_train_sd = stats::sd(x), N_Training = length(x))
  }))
  y_groups <- split(train$production_1000_qt, as.character(train$governorate))
  y_stats <- do.call(rbind, lapply(names(y_groups), function(g) {
    x <- y_groups[[g]]
    data.frame(governorate = g, production_train_mean = mean(x),
               production_train_sd = stats::sd(x),
               production_train_median = stats::median(x), N_Training = length(x))
  }))
  if (any(!is.finite(sm_stats$sm_train_sd)) || any(sm_stats$sm_train_sd <= 0)) {
    stop("Invalid within-governorate soil-moisture SD in this training fold.")
  }
  if (any(!is.finite(y_stats$production_train_sd)) || any(y_stats$production_train_sd <= 0)) {
    stop("Invalid within-governorate production SD in this training fold.")
  }
  list(sm = sm_stats, production = y_stats,
       governorate_levels = sort(unique(as.character(train$governorate))))
}

apply_features <- function(newdata, stats_obj, response_transformation, include_response = FALSE) {
  d <- newdata
  gov <- as.character(d$governorate)
  unknown <- setdiff(unique(gov), stats_obj$governorate_levels)
  if (length(unknown)) stop("Governorates absent from training fold: ", paste(unknown, collapse = ", "))
  si <- match(gov, stats_obj$sm$governorate)
  yi <- match(gov, stats_obj$production$governorate)
  if (anyNA(si) || anyNA(yi)) stop("Missing training scaling parameters.")
  d$governorate <- factor(gov, levels = stats_obj$governorate_levels)
  d$region <- factor(as.character(d$region), levels = c("North", "Center_South"))
  if (anyNA(d$region)) stop("Invalid region classification.")
  d$sm_z <- (d$sm_oct_nov - stats_obj$sm$sm_train_mean[si]) / stats_obj$sm$sm_train_sd[si]
  d$sm_deficit <- pmax(-1 - d$sm_z, 0)
  d$sm_wet <- pmax(d$sm_z - 1, 0)
  d$year_c <- (d$harvest_year - YEAR_ORIGIN) / 10
  cs <- as.numeric(d$region == "Center_South")
  d$temp_cs <- d$temp_dec_c * cs
  d$sm_z_cs <- d$sm_z * cs
  d$sm_deficit_cs <- d$sm_deficit * cs
  d$sm_wet_cs <- d$sm_wet * cs
  d$ndvi_cs <- d$ndvi_mean * cs
  d$vhi_cs <- d$vhi * cs
  transformed <- c("sm_z", "sm_deficit", "sm_wet", "year_c")
  if (any(!vapply(d[, transformed, drop = FALSE], function(x) all(is.finite(x)), logical(1)))) {
    stop("Non-finite transformed predictors.")
  }
  if (include_response) {
    if (response_transformation == "raw") {
      d$y_model <- d$production_1000_qt
    } else if (response_transformation == "log1p") {
      d$y_model <- log1p(d$production_1000_qt)
    } else if (response_transformation == "governorate_z") {
      d$y_model <- (d$production_1000_qt -
                      stats_obj$production$production_train_mean[yi]) /
        stats_obj$production$production_train_sd[yi]
    } else stop("Unknown response transformation.")
  } else d$y_model <- NA_real_
  d
}

back_transform_predictions <- function(predictions, newdata, transformation, stats_obj) {
  p <- as.numeric(predictions)
  if (transformation == "raw") return(p)
  if (transformation == "log1p") return(expm1(p))
  if (transformation == "governorate_z") {
    i <- match(as.character(newdata$governorate), stats_obj$production$governorate)
    if (anyNA(i)) stop("Missing production scaling parameters for back-transformation.")
    return(p * stats_obj$production$production_train_sd[i] +
             stats_obj$production$production_train_mean[i])
  }
  stop("Unknown response transformation.")
}

fixed_specs <- list(
  baseline_quadratic = paste(
    "y_model ~ governorate + temp_dec_c + temp_cs + sm_z + sm_z_cs +",
    "I(sm_z^2) + rain_jan_mar_mm + ndvi_mean + vhi + year_c"
  ),
  threshold_full = paste(
    "y_model ~ governorate + temp_dec_c + temp_cs +",
    "sm_deficit + sm_deficit_cs + sm_wet + sm_wet_cs +",
    "rain_jan_mar_mm + ndvi_mean + ndvi_cs + vhi + vhi_cs + year_c"
  ),
  threshold_parsimonious = paste(
    "y_model ~ governorate + temp_dec_c + temp_cs +",
    "sm_deficit + sm_deficit_cs + sm_wet +",
    "rain_jan_mar_mm + ndvi_mean + vhi + year_c"
  )
)
complexity_rank <- c(baseline_quadratic = 2L, threshold_full = 3L,
                     threshold_parsimonious = 1L)
response_transformations <- c("raw", "log1p", "governorate_z")
candidates <- expand.grid(
  Fixed_Specification = names(fixed_specs),
  Response_Transformation = response_transformations,
  stringsAsFactors = FALSE
)
candidates$Candidate_ID <- sprintf("Q%02d", seq_len(nrow(candidates)))
candidates$Model_Name <- paste(candidates$Fixed_Specification,
                               candidates$Response_Transformation, sep = "__")
candidates$Complexity_Rank <- complexity_rank[candidates$Fixed_Specification]

make_formula <- function(specification_name) {
  stats::as.formula(fixed_specs[[specification_name]])
}

fit_rq_checked <- function(formula, model_data, tau_value, model_label) {
  mf <- stats::model.frame(formula, data = model_data, na.action = stats::na.fail)
  X <- stats::model.matrix(formula, data = mf)
  y <- stats::model.response(mf)
  if (nrow(X) != nrow(model_data) || any(!is.finite(X)) || any(!is.finite(y))) {
    stop("Invalid or non-finite model matrix/response: ", model_label)
  }
  if (nrow(X) < ncol(X) + 2L) {
    stop("Insufficient observations relative to design columns: ", model_label,
         " (n=", nrow(X), ", p=", ncol(X), ").")
  }
  rank_x <- qr(X)$rank
  if (rank_x < ncol(X)) {
    stop("Rank-deficient design matrix: ", model_label,
         " (rank=", rank_x, ", columns=", ncol(X), ").")
  }

  warning_messages <- character()
  fit_error <- NULL
  model <- tryCatch(
    withCallingHandlers(
      quantreg::rq(
        formula = formula, tau = tau_value, data = mf,
        method = "br", na.action = stats::na.fail, model = TRUE
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(model)) stop("Quantile-regression fitting error for ", model_label, ": ", fit_error)

  b <- stats::coef(model)
  if (length(b) != ncol(X) || any(!is.finite(b))) {
    stop("Incomplete or non-finite coefficients: ", model_label)
  }

  nonunique_warning <- grepl(
    "solution may be non.?unique|non.?unique solution",
    warning_messages, ignore.case = TRUE
  )
  other_warning <- warning_messages[!nonunique_warning]
  fatal_warning <- grepl(
    "singular|rank.?deficient|failed|iteration limit|non.?finite|infinite|did not converge",
    other_warning, ignore.case = TRUE
  )
  if (any(fatal_warning)) {
    stop("Serious numerical warning for ", model_label, ": ",
         paste(unique(other_warning[fatal_warning]), collapse = " | "))
  }
  if (any(nonunique_warning) && !ALLOW_NONUNIQUE_WARNING) {
    stop("Nonunique solution warning not permitted for ", model_label)
  }

  list(
    model = model,
    warnings = unique(warning_messages),
    nonunique_warning = any(nonunique_warning),
    design_rank = rank_x,
    design_columns = ncol(X),
    n = nrow(X),
    tau = tau_value
  )
}

metric_vector <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted)
  a <- actual[ok]
  p <- predicted[ok]
  if (!length(a)) return(c(N = 0, R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_))
  sse <- sum((a - p)^2)
  sst <- sum((a - mean(a))^2)
  c(N = length(a), R2 = if (sst > 0) 1 - sse / sst else NA_real_,
    RMSE = sqrt(mean((a - p)^2)), MAE = mean(abs(a - p)))
}
pinball_loss <- function(actual, predicted, tau) {
  ok <- is.finite(actual) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  r <- actual[ok] - predicted[ok]
  mean((tau - as.numeric(r < 0)) * r)
}

folds <- list(
  list(name = "2000-2006_to_2007-2008", train_end = 2006L, valid_start = 2007L, valid_end = 2008L),
  list(name = "2000-2008_to_2009-2010", train_end = 2008L, valid_start = 2009L, valid_end = 2010L),
  list(name = "2000-2010_to_2011-2012", train_end = 2010L, valid_start = 2011L, valid_end = 2012L)
)

cv_rows <- list()
cv_predictions <- list()
fold_metric_rows <- list()
candidate_diagnostic_rows <- list()

for (i in seq_len(nrow(candidates))) {
  candidate <- candidates[i, , drop = FALSE]
  formula_i <- make_formula(candidate$Fixed_Specification)
  fold_predictions <- list()
  fold_diagnostics <- list()

  for (j in seq_along(folds)) {
    fold <- folds[[j]]
    fold_train_raw <- train_data_raw[train_data_raw$harvest_year <= fold$train_end, , drop = FALSE]
    fold_valid_raw <- train_data_raw[
      train_data_raw$harvest_year >= fold$valid_start &
        train_data_raw$harvest_year <= fold$valid_end, , drop = FALSE
    ]

    outcome <- tryCatch({
      prep <- fit_preprocessor(fold_train_raw)
      tr <- apply_features(fold_train_raw, prep,
                           candidate$Response_Transformation, TRUE)
      va <- apply_features(fold_valid_raw, prep,
                           candidate$Response_Transformation, FALSE)
      fit <- fit_rq_checked(
        formula_i, tr, PRIMARY_TAU,
        paste(candidate$Model_Name, fold$name, sep = "__")
      )
      pred_model <- as.numeric(stats::predict(fit$model, newdata = va))
      pred <- back_transform_predictions(
        pred_model, va, candidate$Response_Transformation, prep
      )
      if (length(pred) != nrow(va) || any(!is.finite(pred))) {
        stop("Validation predictions are missing, non-finite or incorrectly sized.")
      }
      list(
        status = "SUCCESS",
        prediction = data.frame(
          Candidate_ID = candidate$Candidate_ID, Model_Name = candidate$Model_Name,
          Fold = fold$name, harvest_year = va$harvest_year,
          governorate = as.character(va$governorate), region = as.character(va$region),
          actual = va$production_1000_qt, predicted_median = pred
        ),
        diagnostic = data.frame(
          Candidate_ID = candidate$Candidate_ID, Model_Name = candidate$Model_Name,
          Fold = fold$name, Status = "SUCCESS",
          N_Training = nrow(tr), N_Validation = nrow(va),
          Design_Rank = fit$design_rank, Design_Columns = fit$design_columns,
          Nonunique_Warning = fit$nonunique_warning,
          Warning_Count = length(fit$warnings),
          Warning_Details = if (length(fit$warnings)) paste(fit$warnings, collapse = " | ") else NA_character_,
          Failure_Reason = NA_character_
        )
      )
    }, error = function(e) {
      list(
        status = "FAILED",
        diagnostic = data.frame(
          Candidate_ID = candidate$Candidate_ID, Model_Name = candidate$Model_Name,
          Fold = fold$name, Status = "FAILED",
          N_Training = nrow(fold_train_raw), N_Validation = nrow(fold_valid_raw),
          Design_Rank = NA_integer_, Design_Columns = NA_integer_,
          Nonunique_Warning = NA, Warning_Count = NA_integer_,
          Warning_Details = NA_character_,
          Failure_Reason = conditionMessage(e)
        )
      )
    })

    fold_diagnostics[[j]] <- outcome$diagnostic
    if (identical(outcome$status, "SUCCESS")) {
      fold_predictions[[j]] <- outcome$prediction
    }
  }

  candidate_diagnostic_rows[[i]] <- do.call(rbind, fold_diagnostics)
  successful_folds <- length(fold_predictions)

  if (successful_folds == length(folds)) {
    pred_data <- do.call(rbind, fold_predictions)
    m <- metric_vector(pred_data$actual, pred_data$predicted_median)
    cv_rows[[length(cv_rows) + 1L]] <- data.frame(
      Candidate_ID = candidate$Candidate_ID, Model_Name = candidate$Model_Name,
      Fixed_Specification = candidate$Fixed_Specification,
      Response_Transformation = candidate$Response_Transformation,
      CV_N = unname(m["N"]), CV_R2 = unname(m["R2"]),
      CV_RMSE = unname(m["RMSE"]), CV_MAE = unname(m["MAE"]),
      CV_Median_Pinball_Loss = pinball_loss(pred_data$actual, pred_data$predicted_median, PRIMARY_TAU),
      Complexity_Rank = candidate$Complexity_Rank,
      Successful_Folds = successful_folds,
      Folds_With_Nonunique_Warning = sum(vapply(fold_diagnostics, function(x) {
        isTRUE(x$Nonunique_Warning[1])
      }, logical(1))),
      stringsAsFactors = FALSE
    )
    cv_predictions[[candidate$Candidate_ID]] <- pred_data
    fold_metric_rows[[length(fold_metric_rows) + 1L]] <- do.call(rbind,
      lapply(split(pred_data, pred_data$Fold), function(z) {
        mm <- metric_vector(z$actual, z$predicted_median)
        data.frame(
          Candidate_ID = candidate$Candidate_ID, Model_Name = candidate$Model_Name,
          Fold = z$Fold[1], N = unname(mm["N"]),
          RMSE = unname(mm["RMSE"]), MAE = unname(mm["MAE"]),
          Median_Pinball_Loss = pinball_loss(z$actual, z$predicted_median, PRIMARY_TAU)
        )
      })
    )
  }

  cat("Candidate", candidate$Candidate_ID, "of", nrow(candidates), ":",
      candidate$Model_Name, "-",
      successful_folds, "of", length(folds), "folds completed\n")
}

all_candidate_diagnostics <- do.call(rbind, candidate_diagnostic_rows)
write_csv(all_candidate_diagnostics, "quantile_regression_candidate_fit_diagnostics.csv")

if (!length(cv_rows)) {
  failed_summary <- aggregate(
    Status ~ Candidate_ID + Model_Name,
    data = all_candidate_diagnostics,
    FUN = function(x) paste(sum(x == "SUCCESS"), "of", length(x), "folds succeeded")
  )
  write_csv(failed_summary, "quantile_regression_candidate_failure_summary.csv")
  stop(
    "No candidate completed all three rolling-origin folds. ",
    "See quantile_regression_candidate_fit_diagnostics.csv and ",
    "quantile_regression_candidate_failure_summary.csv. ",
    "Nonunique-solution warnings are now recorded rather than automatically rejected."
  )
}

cv_results <- do.call(rbind, cv_rows)
cv_results$Rank_RMSE <- rank(cv_results$CV_RMSE, ties.method = "average")
cv_results$Rank_MAE <- rank(cv_results$CV_MAE, ties.method = "average")
cv_results$Average_Rank <- rowMeans(cv_results[, c("Rank_RMSE", "Rank_MAE")])
cv_results <- cv_results[order(
  cv_results$Average_Rank, cv_results$CV_MAE,
  cv_results$CV_RMSE, cv_results$Complexity_Rank
), , drop = FALSE]
cv_results$CV_Rank <- seq_len(nrow(cv_results))
rownames(cv_results) <- NULL

top_three <- utils::head(cv_results, min(3L, nrow(cv_results)))
selected <- cv_results[1, , drop = FALSE]
selected_candidate <- candidates[candidates$Candidate_ID == selected$Candidate_ID, , drop = FALSE]
selected_formula <- make_formula(selected_candidate$Fixed_Specification)

final_preprocessing <- fit_preprocessor(train_data_raw)
final_train <- apply_features(train_data_raw, final_preprocessing,
                              selected_candidate$Response_Transformation, TRUE)
final_test <- apply_features(test_data_raw, final_preprocessing,
                             selected_candidate$Response_Transformation, FALSE)

final_models <- list()
final_fit_diagnostics <- list()
for (tau_value in TAU_GRID) {
  tau_name <- format(tau_value, nsmall = 2, trim = TRUE)
  fit <- tryCatch(
    fit_rq_checked(
      selected_formula, final_train, tau_value,
      paste0(selected$Model_Name, "__tau_", tau_name)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    final_fit_diagnostics[[tau_name]] <- data.frame(
      Tau = tau_value, Status = "FAILED", N_Training = nrow(final_train),
      Design_Rank = NA_integer_, Design_Columns = NA_integer_,
      Nonunique_Warning = NA, Warning_Count = NA_integer_,
      Warning_Details = NA_character_, Failure_Reason = conditionMessage(fit)
    )
  } else {
    final_models[[tau_name]] <- fit$model
    final_fit_diagnostics[[tau_name]] <- data.frame(
      Tau = tau_value, Status = "SUCCESS", N_Training = nrow(final_train),
      Design_Rank = fit$design_rank, Design_Columns = fit$design_columns,
      Nonunique_Warning = fit$nonunique_warning,
      Warning_Count = length(fit$warnings),
      Warning_Details = if (length(fit$warnings)) paste(fit$warnings, collapse = " | ") else NA_character_,
      Failure_Reason = NA_character_
    )
  }
}
final_fit_diagnostics_table <- do.call(rbind, final_fit_diagnostics)
primary_tau_name <- format(PRIMARY_TAU, nsmall = 2, trim = TRUE)
if (is.null(final_models[[primary_tau_name]])) {
  write_csv(final_fit_diagnostics_table, "quantile_regression_final_fit_diagnostics.csv")
  stop("The selected final median model failed. Inspect final-fit diagnostics.")
}

prediction_matrix <- matrix(
  NA_real_, nrow = nrow(final_test), ncol = length(TAU_GRID),
  dimnames = list(NULL, paste0("tau_", gsub("\\.", "_", format(TAU_GRID, nsmall = 2))))
)
for (tau_value in TAU_GRID) {
  tau_name <- format(tau_value, nsmall = 2, trim = TRUE)
  model <- final_models[[tau_name]]
  if (is.null(model)) next
  p_model <- as.numeric(stats::predict(model, newdata = final_test))
  prediction_matrix[, paste0("tau_", gsub("\\.", "_", tau_name))] <-
    back_transform_predictions(p_model, final_test,
                               selected_candidate$Response_Transformation,
                               final_preprocessing)
}

actual <- test_data_raw$production_1000_qt
median_column <- paste0("tau_", gsub("\\.", "_", primary_tau_name))
median_predictions <- as.numeric(prediction_matrix[, median_column])
if (any(!is.finite(median_predictions))) stop("Final median predictions are non-finite.")
median_predictions_clipped <- pmax(median_predictions, 0)
all_quantiles_available <- all(is.finite(prediction_matrix))

rearranged_matrix <- matrix(NA_real_, nrow(prediction_matrix), ncol(prediction_matrix),
                            dimnames = dimnames(prediction_matrix))
if (all_quantiles_available) {
  quantile_crossing_by_observation <- apply(prediction_matrix, 1, function(x) any(diff(x) < 0))
  quantile_crossing_pair_count <- sum(apply(prediction_matrix, 1, function(x) sum(diff(x) < 0)))
  rearranged_matrix <- t(apply(prediction_matrix, 1, sort))
  colnames(rearranged_matrix) <- colnames(prediction_matrix)
} else {
  quantile_crossing_by_observation <- rep(NA, nrow(prediction_matrix))
  quantile_crossing_pair_count <- NA_integer_
}
rearranged_median_predictions <- if (all_quantiles_available) {
  as.numeric(rearranged_matrix[, median_column])
} else rep(NA_real_, nrow(final_test))

baseline_index <- match(as.character(test_data_raw$governorate),
                        final_preprocessing$production$governorate)
baseline_governorate_median <- final_preprocessing$production$production_train_median[baseline_index]

calculate_performance_row <- function(a, p, label) {
  m <- metric_vector(a, p)
  data.frame(
    Model = selected$Model_Name, Prediction = label,
    N = unname(m["N"]), R2 = unname(m["R2"]),
    RMSE = unname(m["RMSE"]), MAE = unname(m["MAE"]),
    Median_Pinball_Loss = pinball_loss(a, p, PRIMARY_TAU)
  )
}
performance_rows <- list(
  calculate_performance_row(actual, median_predictions, "Median quantile - unclipped; primary"),
  calculate_performance_row(actual, median_predictions_clipped, "Median quantile - clipped at zero; secondary"),
  calculate_performance_row(actual, baseline_governorate_median, "Governorate training-period median baseline")
)
if (all_quantiles_available) {
  performance_rows[[length(performance_rows) + 1L]] <- calculate_performance_row(
    actual, rearranged_median_predictions, "Rearranged median - unclipped; supplementary"
  )
  performance_rows[[length(performance_rows) + 1L]] <- calculate_performance_row(
    actual, pmax(rearranged_median_predictions, 0), "Rearranged median - clipped; supplementary"
  )
}
performance_comparison <- do.call(rbind, performance_rows)
performance_comparison <- performance_comparison[order(performance_comparison$RMSE), , drop = FALSE]
rownames(performance_comparison) <- NULL

prediction_output <- data.frame(
  harvest_year = test_data_raw$harvest_year,
  governorate = as.character(test_data_raw$governorate),
  region = as.character(test_data_raw$region),
  actual_production_1000_qt = actual,
  predicted_median_unclipped = median_predictions,
  predicted_median_clipped = median_predictions_clipped,
  baseline_governorate_training_median = baseline_governorate_median,
  residual_median = actual - median_predictions,
  absolute_error_median = abs(actual - median_predictions),
  quantile_crossing_detected = quantile_crossing_by_observation
)
for (j in seq_along(TAU_GRID)) {
  suffix <- gsub("\\.", "_", format(TAU_GRID[j], nsmall = 2))
  prediction_output[[paste0("predicted_tau_", suffix, "_raw")]] <- prediction_matrix[, j]
  prediction_output[[paste0("predicted_tau_", suffix, "_rearranged")]] <- rearranged_matrix[, j]
}

quantile_evaluation <- do.call(rbind, lapply(seq_along(TAU_GRID), function(j) {
  qhat <- prediction_matrix[, j]
  ok <- is.finite(actual) & is.finite(qhat)
  data.frame(
    Tau = TAU_GRID[j], N = sum(ok),
    Pinball_Loss = pinball_loss(actual, qhat, TAU_GRID[j]),
    Empirical_Coverage = if (any(ok)) mean(actual[ok] <= qhat[ok]) else NA_real_,
    Coverage_Difference = if (any(ok)) mean(actual[ok] <= qhat[ok]) - TAU_GRID[j] else NA_real_,
    Mean_Predicted_Quantile = if (any(ok)) mean(qhat[ok]) else NA_real_
  )
}))

quantile_crossing_summary <- data.frame(
  All_Quantiles_Fitted = all_quantiles_available,
  Test_Observations = nrow(prediction_matrix),
  Observations_With_At_Least_One_Crossing =
    if (all_quantiles_available) sum(quantile_crossing_by_observation) else NA_integer_,
  Proportion_With_Crossing =
    if (all_quantiles_available) mean(quantile_crossing_by_observation) else NA_real_,
  Adjacent_Crossing_Pairs = quantile_crossing_pair_count,
  Rearrangement_Applied_To_Supplementary_Output = all_quantiles_available
)

per_governorate_metrics <- do.call(rbind, lapply(
  split(seq_len(nrow(test_data_raw)), test_data_raw$governorate), function(index) {
    a <- actual[index]
    p <- median_predictions[index]
    pc <- median_predictions_clipped[index]
    data.frame(
      Governorate = as.character(test_data_raw$governorate[index[1]]),
      Region = as.character(test_data_raw$region[index[1]]),
      N_Test = length(index),
      RMSE_Unclipped = sqrt(mean((a - p)^2)),
      MAE_Unclipped = mean(abs(a - p)),
      RMSE_Clipped = sqrt(mean((a - pc)^2)),
      MAE_Clipped = mean(abs(a - pc))
    )
  }
))
per_governorate_metrics <- per_governorate_metrics[
  order(per_governorate_metrics$RMSE_Unclipped), , drop = FALSE
]

coefficient_profile <- do.call(rbind, lapply(names(final_models), function(tau_name) {
  b <- stats::coef(final_models[[tau_name]])
  data.frame(
    Tau = as.numeric(tau_name), Term = names(b), Estimate = as.numeric(b),
    Response_Transformation = selected_candidate$Response_Transformation,
    Estimate_Scale = paste0("Model response scale: ", selected_candidate$Response_Transformation),
    row.names = NULL
  )
}))

median_coefficients <- stats::coef(final_models[[primary_tau_name]])
get_coefficient <- function(term) {
  if (term %in% names(median_coefficients)) as.numeric(median_coefficients[[term]]) else 0
}
regional_slope_row <- function(component, north, center_south) {
  data.frame(
    Component = component, North_Estimate = north,
    Center_South_Estimate = center_south,
    Estimate_Scale = paste0("Median quantile model response scale: ",
                            selected_candidate$Response_Transformation)
  )
}
north_temp <- get_coefficient("temp_dec_c")
south_temp <- north_temp + get_coefficient("temp_cs")
if (selected_candidate$Fixed_Specification == "baseline_quadratic") {
  q2 <- get_coefficient("I(sm_z^2)")
  north_sm <- get_coefficient("sm_z") + 2 * (-1.5) * q2
  south_sm <- get_coefficient("sm_z") + get_coefficient("sm_z_cs") + 2 * (-1.5) * q2
} else {
  north_sm <- -get_coefficient("sm_deficit")
  south_sm <- -(get_coefficient("sm_deficit") + get_coefficient("sm_deficit_cs"))
}
regional_slope_estimates <- do.call(rbind, list(
  regional_slope_row("December temperature slope", north_temp, south_temp),
  regional_slope_row("Soil-moisture derivative at sm_z = -1.5", north_sm, south_sm),
  regional_slope_row("NDVI slope", get_coefficient("ndvi_mean"),
                     get_coefficient("ndvi_mean") + get_coefficient("ndvi_cs")),
  regional_slope_row("VHI slope", get_coefficient("vhi"),
                     get_coefficient("vhi") + get_coefficient("vhi_cs"))
))

all_cv_predictions <- do.call(rbind, cv_predictions)
all_fold_metrics <- if (length(fold_metric_rows)) do.call(rbind, fold_metric_rows) else data.frame()
selection_summary <- data.frame(
  Selected_Candidate_ID = selected$Candidate_ID,
  Selected_Model_Name = selected$Model_Name,
  Fixed_Specification = selected$Fixed_Specification,
  Response_Transformation = selected$Response_Transformation,
  Primary_Quantile = PRIMARY_TAU,
  Selection_Method = "Three-fold chronological rolling-origin validation; average RMSE and MAE ranks",
  Training_Start_Year = TRAIN_START_YEAR, Training_End_Year = TRAIN_END_YEAR,
  Test_Start_Year = TEST_START_YEAR, Test_End_Year = TEST_END_YEAR,
  Training_N = nrow(train_data_raw), Test_N = nrow(test_data_raw),
  All_Final_Quantiles_Successful = all(final_fit_diagnostics_table$Status == "SUCCESS"),
  Quantile_Crossing_Rate_Available = all_quantiles_available,
  stringsAsFactors = FALSE
)

write_csv(cv_results, "quantile_regression_rolling_origin_cv_rankings.csv")
write_csv(top_three, "quantile_regression_top_three_cv_candidates.csv")
write_csv(all_candidate_diagnostics, "quantile_regression_candidate_fit_diagnostics.csv")
write_csv(all_cv_predictions, "quantile_regression_rolling_origin_cv_predictions.csv")
write_csv(all_fold_metrics, "quantile_regression_fold_metrics.csv")
write_csv(final_fit_diagnostics_table, "quantile_regression_final_fit_diagnostics.csv")
write_csv(performance_comparison, "quantile_regression_test_performance.csv")
write_csv(prediction_output, "quantile_regression_test_predictions.csv")
write_csv(quantile_evaluation, "quantile_regression_quantile_calibration.csv")
write_csv(quantile_crossing_summary, "quantile_regression_quantile_crossing_summary.csv")
write_csv(per_governorate_metrics, "quantile_regression_governorate_metrics.csv")
write_csv(coverage_rows, "quantile_regression_governorate_coverage_audit.csv")
write_csv(final_preprocessing$sm, "quantile_regression_training_soil_moisture_scaling.csv")
write_csv(final_preprocessing$production, "quantile_regression_training_production_scaling.csv")
write_csv(coefficient_profile, "quantile_regression_coefficient_profile.csv")
write_csv(regional_slope_estimates, "quantile_regression_regional_slope_estimates.csv")
write_csv(selection_summary, "quantile_regression_selection_summary.csv")

capture.output(selected_formula, file = file.path(OUTPUT_DIR, "quantile_regression_selected_formula.txt"))
capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "quantile_regression_session_info.txt"))
saveRDS(final_models, file = file.path(OUTPUT_DIR, "quantile_regression_models.rds"))
saveRDS(final_models[[primary_tau_name]],
        file = file.path(OUTPUT_DIR, "quantile_regression_median_model.rds"))
saveRDS(list(
  models = final_models, selected_candidate = selected,
  selected_formula = selected_formula,
  response_transformation = selected_candidate$Response_Transformation,
  preprocessing = final_preprocessing, tau_grid = TAU_GRID,
  primary_tau = PRIMARY_TAU, test_quantile_predictions = prediction_matrix,
  test_quantile_predictions_rearranged = rearranged_matrix,
  selection_summary = selection_summary
), file = file.path(OUTPUT_DIR, "quantile_regression_forecasting_bundle.rds"))

png(file.path(OUTPUT_DIR, "quantile_regression_test_actual_vs_predicted.png"),
    width = 1300, height = 1100, res = 160)
plot(actual, median_predictions,
     xlab = "Observed production (1000 quintals)",
     ylab = "Predicted median production (1000 quintals)",
     main = "Quantile regression: observed vs predicted production", pch = 19)
abline(a = 0, b = 1, lty = 2)
dev.off()

png(file.path(OUTPUT_DIR, "quantile_regression_test_residuals_by_year.png"),
    width = 1300, height = 900, res = 160)
plot(test_data_raw$harvest_year, actual - median_predictions,
     xlab = "Harvest year", ylab = "Observed minus predicted production",
     main = "Median quantile-regression test residuals by year", pch = 19)
abline(h = 0, lty = 2)
dev.off()

if (all(is.finite(quantile_evaluation$Empirical_Coverage))) {
  png(file.path(OUTPUT_DIR, "quantile_regression_test_quantile_calibration.png"),
      width = 1100, height = 900, res = 160)
  plot(quantile_evaluation$Tau, quantile_evaluation$Empirical_Coverage,
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "Nominal quantile",
       ylab = "Empirical test coverage: proportion(actual <= predicted quantile)",
       main = "Test-set quantile calibration", pch = 19)
  abline(a = 0, b = 1, lty = 2)
  dev.off()
}

cat("\n========== QUANTILE REGRESSION FORECASTING COMPLETE ==========\n")
cat("Training observations:", nrow(train_data_raw), "\n")
cat("Test observations:", nrow(test_data_raw), "\n")
cat("Selected candidate:", selected$Model_Name, "\n")
cat("Selected formula:\n")
print(selected_formula)
cat("Response transformation:", selected_candidate$Response_Transformation, "\n")
cat("Final quantiles fitted:", paste(names(final_models), collapse = ", "), "\n")
cat("All final quantiles successful:", all(final_fit_diagnostics_table$Status == "SUCCESS"), "\n")
cat("Quantile crossing rate available:", all_quantiles_available, "\n")
if (all_quantiles_available) {
  cat("Test observations with at least one quantile crossing:",
      sum(quantile_crossing_by_observation), "of", length(quantile_crossing_by_observation), "\n")
}
cat("\nTop rolling-origin validation candidates:\n")
print(top_three[, c("CV_Rank", "Model_Name", "CV_R2", "CV_RMSE", "CV_MAE",
                    "CV_Median_Pinball_Loss", "Average_Rank")], row.names = FALSE)
cat("\nFinal test performance:\n")
print(performance_comparison, row.names = FALSE)
cat("\nFinal quantile-fit diagnostics:\n")
print(final_fit_diagnostics_table, row.names = FALSE)
cat("\nRegional median-quantile slope estimates (descriptive; no inferential p-values):\n")
print(regional_slope_estimates, row.names = FALSE)
cat("\nOutputs saved to:\n", OUTPUT_DIR, "\n")
