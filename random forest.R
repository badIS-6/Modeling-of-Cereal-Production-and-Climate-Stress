if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("Package 'ranger' is required. Install it with install.packages('ranger').")
}

set.seed(42)

data <- read.csv(
  "model_table.csv",
  fileEncoding = "UTF-8",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "governorate",
  "crop",
  "harvest_year",
  "production_1000_qt",
  "rain_jan_mar_mm",
  "temp_dec_c",
  "sm_oct_nov",
  "ndvi_mean",
  "vhi"
)

predictor_cols <- c(
  "temp_dec_c",
  "sm_oct_nov",
  "rain_jan_mar_mm",
  "ndvi_mean",
  "vhi"
)

model_cols <- c(
  "production_1000_qt",
  "governorate",
  predictor_cols
)

stopifnot(all(required_cols %in% names(data)))
stopifnot(nrow(data) == 281L)
stopifnot(length(unique(data$governorate)) == 15L)
stopifnot(sum(is.na(data)) == 0L)
stopifnot(sum(data$production_1000_qt == 0) == 11L)
stopifnot(all(vapply(data[, predictor_cols, drop = FALSE], is.numeric, logical(1))))
stopifnot(is.numeric(data$production_1000_qt))
stopifnot(is.numeric(data$harvest_year))
stopifnot(all(is.finite(as.matrix(data[, predictor_cols, drop = FALSE]))))

governorate_levels <- sort(unique(data$governorate))
data$governorate <- factor(data$governorate, levels = governorate_levels)

train_data <- data[
  data$harvest_year >= 2000 & data$harvest_year <= 2012,
  ,
  drop = FALSE
]

test_data <- data[
  data$harvest_year >= 2013 & data$harvest_year <= 2018,
  ,
  drop = FALSE
]

stopifnot(nrow(train_data) == 192L)
stopifnot(nrow(test_data) == 89L)
stopifnot(identical(levels(train_data$governorate), levels(test_data$governorate)))
stopifnot(setequal(unique(as.character(train_data$governorate)), governorate_levels))
stopifnot(setequal(unique(as.character(test_data$governorate)), governorate_levels))

rf_formula <- stats::reformulate(
  termlabels = c("governorate", predictor_cols),
  response = "production_1000_qt"
)

folds <- data.frame(
  fold = 1:3,
  train_end = c(2007L, 2009L, 2011L),
  validation_start = c(2008L, 2010L, 2012L),
  validation_end = c(2009L, 2011L, 2012L)
)

tune_grid <- expand.grid(
  mtry = seq_along(c("governorate", predictor_cols)),
  min.node.size = c(3L, 5L, 8L, 12L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

tuning_results <- vector("list", nrow(tune_grid))

for (grid_index in seq_len(nrow(tune_grid))) {
  configuration_predictions <- vector("list", nrow(folds))

  for (fold_index in seq_len(nrow(folds))) {
    fold_spec <- folds[fold_index, ]

    fold_train <- train_data[
      train_data$harvest_year <= fold_spec$train_end,
      model_cols,
      drop = FALSE
    ]

    fold_validation <- train_data[
      train_data$harvest_year >= fold_spec$validation_start &
        train_data$harvest_year <= fold_spec$validation_end,
      model_cols,
      drop = FALSE
    ]

    stopifnot(nrow(fold_train) > 0L, nrow(fold_validation) > 0L)
    stopifnot(setequal(unique(as.character(fold_train$governorate)), governorate_levels))

    set.seed(42L + fold_index)

    cv_model <- ranger::ranger(
      formula = rf_formula,
      data = fold_train,
      num.trees = 500L,
      mtry = as.integer(tune_grid$mtry[grid_index]),
      min.node.size = as.integer(tune_grid$min.node.size[grid_index]),
      importance = "none",
      respect.unordered.factors = "order",
      num.threads = 1L,
      seed = 42L + fold_index
    )

    cv_prediction <- stats::predict(
      cv_model,
      data = fold_validation
    )$predictions

    configuration_predictions[[fold_index]] <- data.frame(
      actual = fold_validation$production_1000_qt,
      predicted = as.numeric(cv_prediction)
    )
  }

  cv_predictions <- do.call(rbind, configuration_predictions)
  cv_sse <- sum((cv_predictions$actual - cv_predictions$predicted)^2)
  cv_sst <- sum((cv_predictions$actual - mean(cv_predictions$actual))^2)

  tuning_results[[grid_index]] <- data.frame(
    mtry = as.integer(tune_grid$mtry[grid_index]),
    min.node.size = as.integer(tune_grid$min.node.size[grid_index]),
    CV_R2 = 1 - cv_sse / cv_sst,
    CV_RMSE = sqrt(mean((cv_predictions$actual - cv_predictions$predicted)^2)),
    CV_MAE = mean(abs(cv_predictions$actual - cv_predictions$predicted))
  )
}

tuning_results <- do.call(rbind, tuning_results)
tuning_results <- tuning_results[
  order(tuning_results$CV_RMSE, tuning_results$CV_MAE),
]
rownames(tuning_results) <- NULL

best_parameters <- tuning_results[1, , drop = FALSE]

final_train_data <- train_data[, model_cols, drop = FALSE]
final_test_data <- test_data[, model_cols, drop = FALSE]

set.seed(42)

rf_final <- ranger::ranger(
  formula = rf_formula,
  data = final_train_data,
  num.trees = 1000L,
  mtry = as.integer(best_parameters$mtry),
  min.node.size = as.integer(best_parameters$min.node.size),
  importance = "permutation",
  scale.permutation.importance = TRUE,
  respect.unordered.factors = "order",
  num.threads = 1L,
  seed = 42L
)

predictions <- as.numeric(
  stats::predict(rf_final, data = final_test_data)$predictions
)

actual <- test_data$production_1000_qt
predictions_clipped <- pmax(predictions, 0)

calculate_metrics <- function(actual, predicted, label) {
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)

  data.frame(
    Prediction = label,
    R2 = 1 - sse / sst,
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted)),
    stringsAsFactors = FALSE
  )
}

performance <- rbind(
  calculate_metrics(actual, predictions, "Unclipped"),
  calculate_metrics(actual, predictions_clipped, "Clipped at zero")
)

performance <- performance[order(performance$RMSE, performance$MAE), ]
rownames(performance) <- NULL
performance$mtry <- as.integer(best_parameters$mtry)
performance$min.node.size <- as.integer(best_parameters$min.node.size)

test_results <- data.frame(
  governorate = test_data$governorate,
  harvest_year = test_data$harvest_year,
  actual_production = actual,
  predicted_production = predictions,
  predicted_clipped = predictions_clipped,
  error = actual - predictions,
  absolute_error = abs(actual - predictions)
)

per_governorate <- do.call(
  rbind,
  lapply(governorate_levels, function(governorate_name) {
    governorate_rows <- test_results[
      as.character(test_results$governorate) == governorate_name,
      ,
      drop = FALSE
    ]

    data.frame(
      governorate = governorate_name,
      n = nrow(governorate_rows),
      RMSE = sqrt(mean(governorate_rows$error^2)),
      MAE = mean(governorate_rows$absolute_error),
      RMSE_clipped = sqrt(mean((
        governorate_rows$actual_production -
          governorate_rows$predicted_clipped
      )^2)),
      MAE_clipped = mean(abs(
        governorate_rows$actual_production -
          governorate_rows$predicted_clipped
      )),
      stringsAsFactors = FALSE
    )
  })
)

variable_importance <- data.frame(
  variable = names(rf_final$variable.importance),
  permutation_importance = as.numeric(rf_final$variable.importance),
  stringsAsFactors = FALSE
)

variable_importance <- variable_importance[
  order(variable_importance$permutation_importance, decreasing = TRUE),
]
rownames(variable_importance) <- NULL

write.csv(
  tuning_results,
  "random_forest_tuning_results.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  performance,
  "random_forest_metrics.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  test_results,
  "random_forest_predictions.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  per_governorate,
  "random_forest_governorate_metrics.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  variable_importance,
  "random_forest_variable_importance.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

saveRDS(rf_final, "random_forest_model.rds")
writeLines(capture.output(sessionInfo()), "random_forest_session_info.txt")

cat("Training observations:", nrow(train_data), "\n")
cat("Testing observations:", nrow(test_data), "\n")
cat("Best mtry:", best_parameters$mtry, "\n")
cat("Best min.node.size:", best_parameters$min.node.size, "\n")
cat("\nTuning results:\n")
print(tuning_results, row.names = FALSE)
cat("\nTest-set performance:\n")
print(performance, row.names = FALSE)
cat("\nNegative predictions:", sum(predictions < 0), "\n")
cat("\nPermutation variable importance:\n")
print(variable_importance, row.names = FALSE)
cat("\nPer-governorate test metrics:\n")
print(per_governorate, row.names = FALSE)
