data_path <- "model_table.csv"

out_metrics <- "linear_baseline_test_metrics.csv"
out_predictions <- "linear_baseline_test_predictions.csv"
out_coefficients <- "linear_baseline_coefficients.csv"
out_diagnostics <- "linear_baseline_diagnostics.pdf"

predictors <- c(
  "temp_dec_c",
  "sm_oct_nov",
  "rain_jan_mar_mm",
  "ndvi_mean",
  "vhi"
)

target <- "production_1000_qt"
required_columns <- c(
  "governorate", "crop", "harvest_year", target, predictors
)

if (!file.exists(data_path)) {
  stop(
    "Cannot find '", data_path, "'. Set the working directory to the folder ",
    "containing model_table.csv or update data_path."
  )
}

data <- read.csv(
  data_path,
  fileEncoding = "UTF-8",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0L) {
  stop(
    "The input CSV is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

numeric_columns <- c("harvest_year", target, predictors)
for (column in numeric_columns) {
  data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
}

if (anyNA(data[required_columns])) {
  stop("Missing values were found in required columns. Resolve them explicitly; this script does not impute or delete rows.")
}

numeric_matrix <- as.matrix(data[numeric_columns])
if (any(!is.finite(numeric_matrix))) {
  stop("Non-finite numeric values (Inf or -Inf) were found in the required numeric columns.")
}

if (any(data$harvest_year != as.integer(data$harvest_year))) {
  stop("harvest_year contains non-integer values.")
}
data$harvest_year <- as.integer(data$harvest_year)

if (any(trimws(data$governorate) == "") || any(trimws(data$crop) == "")) {
  stop("Blank governorate or crop labels were found.")
}

observed_crops <- unique(trimws(data$crop))
if (length(observed_crops) != 1L || tolower(observed_crops) != "durum_wheat") {
  stop(
    "Expected only crop = 'durum_wheat'; found: ",
    paste(observed_crops, collapse = ", ")
  )
}

if (any(data$harvest_year < 2000L | data$harvest_year > 2018L)) {
  stop("The dataset contains harvest years outside the specified 2000-2018 study period.")
}

if (anyDuplicated(data[c("governorate", "harvest_year")]) > 0L) {
  stop("Duplicate governorate-harvest_year records were found.")
}

if (nrow(data) != 281L) {
  stop("Expected 281 observations but found ", nrow(data), ". Check the input file before modeling.")
}
if (length(unique(data$governorate)) != 15L) {
  stop("Expected 15 governorates but found ", length(unique(data$governorate)), ".")
}
if (sum(data[[target]] == 0) != 11L) {
  stop("Expected 11 zero-production observations. Check the input data; no zero values are removed or transformed.")
}

data <- data[order(data$harvest_year, data$governorate), ]

cat("Data validation passed.\n")
cat("Total observations retained:", nrow(data), "\n")
cat("Governorates:", length(unique(data$governorate)), "\n")
cat("Zero-production observations retained:", sum(data[[target]] == 0), "\n\n")

train_data <- data[
  data$harvest_year >= 2000L & data$harvest_year <= 2012L,
  , drop = FALSE
]

test_data <- data[
  data$harvest_year >= 2013L & data$harvest_year <= 2018L,
  , drop = FALSE
]

if (nrow(train_data) != 192L || nrow(test_data) != 89L) {
  stop(
    "Unexpected temporal split sizes. Expected 192 training and 89 test rows; found ",
    nrow(train_data), " training and ", nrow(test_data), " test rows."
  )
}

if (max(train_data$harvest_year) >= min(test_data$harvest_year)) {
  stop("Temporal leakage detected: training years must precede all test years.")
}

cat("Training period:", min(train_data$harvest_year), "-",
    max(train_data$harvest_year), "| n =", nrow(train_data), "\n")
cat("Testing period :", min(test_data$harvest_year), "-",
    max(test_data$harvest_year), "| n =", nrow(test_data), "\n\n")

baseline_formula <- production_1000_qt ~
  temp_dec_c +
  sm_oct_nov +
  rain_jan_mar_mm +
  ndvi_mean +
  vhi

baseline_lm <- lm(
  formula = baseline_formula,
  data = train_data,
  na.action = na.fail
)

baseline_design <- model.matrix(baseline_lm)
if (qr(baseline_design)$rank < ncol(baseline_design)) {
  print(alias(baseline_lm))
  stop("The OLS design matrix is rank deficient. Inspect predictor dependencies before proceeding.")
}

cat("================ OLS BASELINE SUMMARY ================\n")
print(summary(baseline_lm))

ols_summary <- summary(baseline_lm)
coefficient_table <- as.data.frame(ols_summary$coefficients)
coefficient_table$term <- rownames(coefficient_table)
rownames(coefficient_table) <- NULL

names(coefficient_table)[1:4] <- c(
  "estimate", "std_error", "t_statistic", "p_value"
)

confidence_intervals <- confint(baseline_lm, level = 0.95)
coefficient_table$conf_low_95 <- confidence_intervals[
  coefficient_table$term, 1
]
coefficient_table$conf_high_95 <- confidence_intervals[
  coefficient_table$term, 2
]
coefficient_table <- coefficient_table[
  , c("term", "estimate", "std_error", "t_statistic", "p_value",
      "conf_low_95", "conf_high_95")
]

cat("\n================ COEFFICIENT ESTIMATES ================\n")
print(coefficient_table, row.names = FALSE)

predicted_raw <- as.numeric(
  predict(baseline_lm, newdata = test_data)
)

if (anyNA(predicted_raw) || any(!is.finite(predicted_raw))) {
  stop("The baseline model generated missing or non-finite test predictions.")
}

predicted_clipped <- pmax(predicted_raw, 0)

actual_test <- test_data[[target]]

calculate_metrics <- function(actual, predicted, prediction_rule) {
  if (length(actual) != length(predicted)) {
    stop("Actual and predicted vectors have different lengths.")
  }
  if (anyNA(actual) || anyNA(predicted)) {
    stop("Missing values found while calculating test metrics.")
  }

  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)

  test_r2 <- if (sst > 0) 1 - sse / sst else NA_real_

  data.frame(
    model = "Pooled OLS baseline",
    prediction_rule = prediction_rule,
    n_test = length(actual),
    R2_test = test_r2,
    RMSE = sqrt(mean((actual - predicted)^2)),
    MAE = mean(abs(actual - predicted)),
    stringsAsFactors = FALSE
  )
}

metrics_table <- rbind(
  calculate_metrics(actual_test, predicted_raw, "Unclipped predictions"),
  calculate_metrics(actual_test, predicted_clipped, "Negative predictions clipped at zero")
)

metrics_table <- metrics_table[
  order(metrics_table$RMSE, metrics_table$MAE),
  , drop = FALSE
]
rownames(metrics_table) <- NULL

cat("\n================ TEST-SET PERFORMANCE ================\n")
cat("Negative raw predictions:", sum(predicted_raw < 0), "of", length(predicted_raw), "\n")
print(metrics_table, row.names = FALSE, digits = 5)

test_predictions <- data.frame(
  governorate = test_data$governorate,
  harvest_year = test_data$harvest_year,
  observed_production_1000_qt = actual_test,
  predicted_raw_1000_qt = predicted_raw,
  predicted_clipped_1000_qt = predicted_clipped,
  error_raw = actual_test - predicted_raw,
  absolute_error_raw = abs(actual_test - predicted_raw),
  error_clipped = actual_test - predicted_clipped,
  absolute_error_clipped = abs(actual_test - predicted_clipped),
  stringsAsFactors = FALSE
)

write.csv(
  metrics_table,
  file = out_metrics,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  test_predictions,
  file = out_predictions,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  coefficient_table,
  file = out_coefficients,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

pdf(out_diagnostics, width = 8, height = 8, onefile = TRUE)
par(mfrow = c(2, 2))
plot(baseline_lm)
dev.off()

cat("\nSaved outputs:\n")
cat(" -", out_metrics, "\n")
cat(" -", out_predictions, "\n")
cat(" -", out_coefficients, "\n")
cat(" -", out_diagnostics, "\n")
