
suppressPackageStartupMessages(library(quantreg))

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

if (!all(required_cols %in% names(data))) {
  stop("The dataset does not contain all required columns.")
}

if (nrow(data) != 281L) {
  stop("Expected 281 observations; found ", nrow(data), ".")
}

if (anyNA(data[, required_cols])) {
  stop("Missing values were found in required columns.")
}

if (length(unique(data$governorate)) != 15L) {
  stop("Expected 15 governorates.")
}

if (sum(data$production_1000_qt == 0) != 11L) {
  stop("Expected 11 zero-production observations.")
}

if (!all(data$crop == "durum_wheat")) {
  stop("The dataset contains crop categories other than durum_wheat.")
}

governorate_levels <- c(
  "Béja", "Gabès", "Gafsa", "Jendouba", "Kairouan",
  "Kasserine", "Le Kef", "Manouba", "Medenine", "Nabeul",
  "Sidi Bouzid", "Siliana", "Sousse", "Tunis", "Zaghouan"
)

north_governorates <- c(
  "Béja", "Jendouba", "Le Kef", "Siliana",
  "Zaghouan", "Manouba", "Nabeul", "Tunis"
)

center_south_governorates <- c(
  "Kairouan", "Kasserine", "Sidi Bouzid", "Sousse",
  "Gafsa", "Gabès", "Medenine"
)

unmapped_governorates <- setdiff(
  unique(data$governorate),
  c(north_governorates, center_south_governorates)
)

if (length(unmapped_governorates) > 0) {
  stop("Unmapped governorates: ", paste(unmapped_governorates, collapse = ", "))
}

if (length(setdiff(unique(data$governorate), governorate_levels)) > 0) {
  stop("Unexpected governorate names were found.")
}

data$governorate <- factor(data$governorate, levels = governorate_levels)
data$region <- factor(
  ifelse(as.character(data$governorate) %in% north_governorates, "North", "Center_South"),
  levels = c("North", "Center_South")
)

data <- data[order(data$harvest_year, data$governorate), ]

train_data <- data[data$harvest_year >= 2000 & data$harvest_year <= 2012, ]
test_data <- data[data$harvest_year >= 2013 & data$harvest_year <= 2018, ]

if (nrow(train_data) != 192L || nrow(test_data) != 89L) {
  stop("The temporal split does not match the expected 192/89 observations.")
}

sm_by_governorate <- split(
  train_data$sm_oct_nov,
  as.character(train_data$governorate)
)

sm_stats_matrix <- do.call(
  rbind,
  lapply(sm_by_governorate, function(x) {
    c(train_mean = mean(x), train_sd = sd(x))
  })
)

sm_stats <- data.frame(
  governorate = rownames(sm_stats_matrix),
  sm_train_mean = sm_stats_matrix[, "train_mean"],
  sm_train_sd = sm_stats_matrix[, "train_sd"],
  row.names = NULL
)

if (any(!is.finite(sm_stats$sm_train_sd)) || any(sm_stats$sm_train_sd <= 0)) {
  stop("At least one governorate has an invalid training-period soil-moisture SD.")
}

train_sm_index <- match(as.character(train_data$governorate), sm_stats$governorate)
test_sm_index <- match(as.character(test_data$governorate), sm_stats$governorate)

if (anyNA(train_sm_index) || anyNA(test_sm_index)) {
  stop("A governorate in the training or test set has no training soil-moisture parameters.")
}

train_data$sm_z <- (train_data$sm_oct_nov - sm_stats$sm_train_mean[train_sm_index]) / sm_stats$sm_train_sd[train_sm_index]
test_data$sm_z <- (test_data$sm_oct_nov - sm_stats$sm_train_mean[test_sm_index]) / sm_stats$sm_train_sd[test_sm_index]

train_data$sm_deficit <- pmax(-1 - train_data$sm_z, 0)
train_data$sm_wet <- pmax(train_data$sm_z - 1, 0)
test_data$sm_deficit <- pmax(-1 - test_data$sm_z, 0)
test_data$sm_wet <- pmax(test_data$sm_z - 1, 0)

rq_formula <- production_1000_qt ~
  governorate +
  temp_dec_c +
  sm_deficit +
  sm_wet +
  rain_jan_mar_mm +
  ndvi_mean +
  vhi +
  region:temp_dec_c +
  region:ndvi_mean +
  region:vhi

rq_design <- model.matrix(rq_formula, data = train_data)
if (qr(rq_design)$rank < ncol(rq_design)) {
  stop("The quantile-regression design matrix is rank deficient.")
}

set.seed(42)

tau_grid <- c(0.10, 0.25, 0.50, 0.75, 0.90)

rq_models <- setNames(
  lapply(tau_grid, function(tau_value) {
    quantreg::rq(
      formula = rq_formula,
      tau = tau_value,
      data = train_data,
      method = "fn",
      na.action = na.fail
    )
  }),
  as.character(tau_grid)
)

rq_median <- rq_models[["0.5"]]

prediction_list <- lapply(
  rq_models,
  predict,
  newdata = test_data
)

prediction_table <- as.data.frame(prediction_list)
names(prediction_table) <- paste0(
  "predicted_tau_",
  gsub("\\.", "_", names(rq_models))
)

actual <- test_data$production_1000_qt
median_predictions <- as.numeric(prediction_list[["0.5"]])
median_predictions_clipped <- pmax(median_predictions, 0)

calculate_metrics <- function(actual_values, predicted_values, prediction_type) {
  sse <- sum((actual_values - predicted_values)^2)
  sst <- sum((actual_values - mean(actual_values))^2)
  data.frame(
    Prediction = prediction_type,
    R2 = 1 - sse / sst,
    RMSE = sqrt(mean((actual_values - predicted_values)^2)),
    MAE = mean(abs(actual_values - predicted_values))
  )
}

performance_comparison <- rbind(
  calculate_metrics(actual, median_predictions, "Median quantile - unclipped"),
  calculate_metrics(actual, median_predictions_clipped, "Median quantile - clipped at zero")
)

performance_comparison <- performance_comparison[
  order(performance_comparison$RMSE),
]
rownames(performance_comparison) <- NULL

quantile_profile <- do.call(
  rbind,
  lapply(names(rq_models), function(tau_name) {
    model_coefficients <- coef(rq_models[[tau_name]])
    data.frame(
      Tau = as.numeric(tau_name),
      Term = names(model_coefficients),
      Estimate = as.numeric(model_coefficients),
      row.names = NULL
    )
  })
)

quantile_inference <- do.call(
  rbind,
  lapply(names(rq_models), function(tau_name) {
    model_summary <- summary(rq_models[[tau_name]], se = "nid")
    coefficient_matrix <- as.matrix(model_summary$coefficients)
    data.frame(
      Tau = as.numeric(tau_name),
      Term = rownames(coefficient_matrix),
      Estimate = coefficient_matrix[, 1],
      Std_Error = coefficient_matrix[, 2],
      Statistic = coefficient_matrix[, 3],
      P_value = coefficient_matrix[, 4],
      CI_Lower = coefficient_matrix[, 1] - qnorm(0.975) * coefficient_matrix[, 2],
      CI_Upper = coefficient_matrix[, 1] + qnorm(0.975) * coefficient_matrix[, 2],
      row.names = NULL
    )
  })
)

median_summary <- summary(rq_median, se = "nid", covariance = TRUE)
median_coefficient_matrix <- as.matrix(median_summary$coefficients)

median_inference <- data.frame(
  Term = rownames(median_coefficient_matrix),
  Estimate = median_coefficient_matrix[, 1],
  Std_Error = median_coefficient_matrix[, 2],
  Statistic = median_coefficient_matrix[, 3],
  P_value = median_coefficient_matrix[, 4],
  CI_Lower = median_coefficient_matrix[, 1] - qnorm(0.975) * median_coefficient_matrix[, 2],
  CI_Upper = median_coefficient_matrix[, 1] + qnorm(0.975) * median_coefficient_matrix[, 2],
  row.names = NULL
)

beta <- coef(rq_median)
rq_vcov <- median_summary$cov

if (is.null(rq_vcov)) {
  stop("The quantreg summary did not return a covariance matrix.")
}

rq_vcov <- as.matrix(rq_vcov)

if (is.null(names(beta)) || any(names(beta) == "")) {
  stop("The median regression coefficients do not have valid names.")
}

if (!identical(dim(rq_vcov), c(length(beta), length(beta)))) {
  stop("The covariance matrix dimensions do not match the coefficient vector.")
}

if (is.null(rownames(rq_vcov))) {
  rownames(rq_vcov) <- names(beta)
}

if (is.null(colnames(rq_vcov))) {
  colnames(rq_vcov) <- names(beta)
}

rq_vcov <- rq_vcov[names(beta), names(beta), drop = FALSE]

find_region_interaction <- function(variable_name) {
  candidates <- c(
    paste0("regionCenter_South:", variable_name),
    paste0(variable_name, ":regionCenter_South")
  )
  matches <- candidates[candidates %in% names(beta)]
  if (length(matches) != 1L) {
    stop("Could not uniquely identify the Center_South interaction for ", variable_name)
  }
  matches
}

calculate_slope <- function(variable_name, region_name) {
  contrast <- setNames(rep(0, length(beta)), names(beta))
  if (!variable_name %in% names(beta)) {
    stop("Missing main slope coefficient: ", variable_name)
  }
  contrast[variable_name] <- 1
  if (region_name == "Center_South") {
    interaction_name <- find_region_interaction(variable_name)
    contrast[interaction_name] <- 1
  }
  estimate <- sum(contrast * beta)
  standard_error <- sqrt(drop(t(contrast) %*% rq_vcov %*% contrast))
  statistic <- estimate / standard_error
  data.frame(
    Variable = variable_name,
    Region = region_name,
    Estimate = estimate,
    Std_Error = standard_error,
    Statistic = statistic,
    P_value = 2 * pnorm(abs(statistic), lower.tail = FALSE),
    CI_Lower = estimate - qnorm(0.975) * standard_error,
    CI_Upper = estimate + qnorm(0.975) * standard_error
  )
}

regional_slopes <- do.call(
  rbind,
  lapply(
    c("temp_dec_c", "ndvi_mean", "vhi"),
    function(variable_name) {
      rbind(
        calculate_slope(variable_name, "North"),
        calculate_slope(variable_name, "Center_South")
      )
    }
  )
)

rq_slope_reduced_formula <- production_1000_qt ~
  governorate +
  temp_dec_c +
  sm_deficit +
  sm_wet +
  rain_jan_mar_mm +
  ndvi_mean +
  vhi

rq_slope_reduced <- quantreg::rq(
  formula = rq_slope_reduced_formula,
  tau = 0.50,
  data = train_data,
  method = "fn",
  na.action = na.fail
)

regional_slope_test <- anova(
  rq_slope_reduced,
  rq_median,
  test = "Wald",
  se = "nid"
)

rq_region_alt_formula <- production_1000_qt ~
  region +
  temp_dec_c +
  sm_deficit +
  sm_wet +
  rain_jan_mar_mm +
  ndvi_mean +
  vhi +
  region:temp_dec_c +
  region:ndvi_mean +
  region:vhi

rq_region_alt <- quantreg::rq(
  formula = rq_region_alt_formula,
  tau = 0.50,
  data = train_data,
  method = "fn",
  na.action = na.fail
)

rq_region_alt_summary <- summary(
  rq_region_alt,
  se = "boot",
  R = 999
)

region_coefficient_names <- grep(
  "^region",
  rownames(rq_region_alt_summary$coefficients),
  value = TRUE
)

if (length(region_coefficient_names) == 0) {
  stop("No region main-effect coefficient was found in the regional model.")
}

regional_baseline_test <- data.frame(
  Term = region_coefficient_names,
  rq_region_alt_summary$coefficients[region_coefficient_names, , drop = FALSE],
  check.names = FALSE,
  row.names = NULL
)

per_governorate_metrics <- do.call(
  rbind,
  lapply(split(seq_len(nrow(test_data)), test_data$governorate), function(index) {
    governorate_actual <- actual[index]
    governorate_predicted <- median_predictions[index]
    data.frame(
      Governorate = as.character(test_data$governorate[index[1]]),
      N = length(index),
      RMSE = sqrt(mean((governorate_actual - governorate_predicted)^2)),
      MAE = mean(abs(governorate_actual - governorate_predicted))
    )
  })
)

per_governorate_metrics <- per_governorate_metrics[
  order(per_governorate_metrics$RMSE),
]

median_predictions_output <- data.frame(
  governorate = as.character(test_data$governorate),
  region = as.character(test_data$region),
  harvest_year = test_data$harvest_year,
  actual_production_1000_qt = actual,
  prediction_table,
  predicted_median_clipped = median_predictions_clipped,
  residual_median = actual - median_predictions,
  absolute_error_median = abs(actual - median_predictions),
  check.names = FALSE
)

write.csv(performance_comparison, "quantile_regression_metrics.csv", row.names = FALSE)
write.csv(median_predictions_output, "quantile_regression_test_predictions.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(quantile_profile, "quantile_regression_coefficient_profile.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(quantile_inference, "quantile_regression_quantile_inference.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(median_inference, "quantile_regression_median_inference.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(regional_slopes, "quantile_regression_regional_slopes.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per_governorate_metrics, "quantile_regression_governorate_metrics.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sm_stats, "quantile_regression_training_soil_moisture_scaling.csv", row.names = FALSE, fileEncoding = "UTF-8")
regional_slope_test_output <- capture.output(print(regional_slope_test))
write.csv(data.frame(Output = regional_slope_test_output), "quantile_regression_regional_slope_test.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(as.data.frame(regional_baseline_test), "quantile_regression_regional_baseline_test.csv", row.names = FALSE)

saveRDS(rq_models, "quantile_regression_models.rds")
saveRDS(rq_median, "quantile_regression_median_model.rds")

cat("Training observations:", nrow(train_data), "\n")
cat("Testing observations:", nrow(test_data), "\n")
cat("\nMedian quantile-regression summary:\n")
print(median_summary)
cat("\nTest-set performance:\n")
print(performance_comparison, row.names = FALSE)
cat("\nRegional slopes at the median quantile:\n")
print(regional_slopes, row.names = FALSE)
cat("\nJoint regional slope test:\n")
print(regional_slope_test)
cat("\nRegional baseline-intercept test:\n")
print(regional_baseline_test)
cat("\nPer-governorate test metrics:\n")
print(per_governorate_metrics, row.names = FALSE)
