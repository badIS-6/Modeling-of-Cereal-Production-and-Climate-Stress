
suppressPackageStartupMessages(library(mgcv))

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

if (!all(required_cols %in% names(data))) {
  stop("The dataset does not contain all required columns.")
}

if (nrow(data) != 281L) {
  stop("Expected 281 observations; found ", nrow(data), ".")
}

if (anyNA(data[, required_cols])) {
  stop("Missing values were found in required columns.")
}

numeric_cols <- c(
  "harvest_year",
  "production_1000_qt",
  "rain_jan_mar_mm",
  "temp_dec_c",
  "sm_oct_nov",
  "ndvi_mean",
  "vhi"
)

if (!all(vapply(data[, numeric_cols, drop = FALSE], is.numeric, logical(1)))) {
  stop("One or more required numeric columns were imported as non-numeric.")
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

sm_stats <- do.call(
  rbind,
  lapply(names(sm_by_governorate), function(governorate_name) {
    values <- sm_by_governorate[[governorate_name]]
    data.frame(
      governorate = governorate_name,
      sm_train_mean = mean(values),
      sm_train_sd = sd(values),
      stringsAsFactors = FALSE
    )
  })
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

if (any(!is.finite(train_data$sm_z)) || any(!is.finite(test_data$sm_z))) {
  stop("Non-finite standardized soil-moisture values were generated.")
}

if (!identical(levels(train_data$governorate), levels(test_data$governorate)) ||
    !identical(levels(train_data$region), levels(test_data$region))) {
  stop("Training and test factor levels are inconsistent.")
}

gam_formula <- production_1000_qt ~
  region +
  s(governorate, bs = "re") +
  s(temp_dec_c, k = 5, bs = "tp") +
  s(sm_z, k = 5, bs = "tp") +
  s(rain_jan_mar_mm, k = 5, bs = "tp") +
  s(ndvi_mean, by = region, k = 4, bs = "tp") +
  s(vhi, by = region, k = 4, bs = "tp")

gam_model <- mgcv::gam(
  formula = gam_formula,
  data = train_data,
  family = gaussian(link = "identity"),
  method = "REML",
  select = TRUE,
  gamma = 1.4,
  na.action = na.fail
)

if (!isTRUE(gam_model$converged)) {
  warning("The GAM fitting routine did not report convergence.")
}

predicted_unclipped <- as.numeric(predict(gam_model, newdata = test_data, type = "response"))
predicted_clipped <- pmax(predicted_unclipped, 0)
actual <- test_data$production_1000_qt

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
  calculate_metrics(actual, predicted_unclipped, "GAM - unclipped"),
  calculate_metrics(actual, predicted_clipped, "GAM - clipped at zero")
)

performance_comparison <- performance_comparison[order(performance_comparison$RMSE), ]
rownames(performance_comparison) <- NULL

test_predictions <- data.frame(
  harvest_year = test_data$harvest_year,
  governorate = as.character(test_data$governorate),
  region = as.character(test_data$region),
  actual = actual,
  predicted_unclipped = predicted_unclipped,
  predicted_clipped = predicted_clipped,
  residual_unclipped = actual - predicted_unclipped,
  residual_clipped = actual - predicted_clipped
)

governorate_indices <- split(seq_len(nrow(test_predictions)), test_predictions$governorate)

per_governorate_metrics <- do.call(
  rbind,
  lapply(names(governorate_indices), function(governorate_name) {
    selected <- test_predictions[governorate_indices[[governorate_name]], ]
    data.frame(
      governorate = governorate_name,
      n_test = nrow(selected),
      RMSE_unclipped = sqrt(mean((selected$actual - selected$predicted_unclipped)^2)),
      MAE_unclipped = mean(abs(selected$actual - selected$predicted_unclipped)),
      RMSE_clipped = sqrt(mean((selected$actual - selected$predicted_clipped)^2)),
      MAE_clipped = mean(abs(selected$actual - selected$predicted_clipped)),
      stringsAsFactors = FALSE
    )
  })
)

summary_object <- summary(gam_model)
parametric_matrix <- as.matrix(summary_object$p.table)
parametric_inference <- data.frame(
  Term = rownames(parametric_matrix),
  Estimate = parametric_matrix[, 1],
  Std_Error = parametric_matrix[, 2],
  Statistic = parametric_matrix[, 3],
  P_value = parametric_matrix[, 4],
  CI_Lower = parametric_matrix[, 1] - qnorm(0.975) * parametric_matrix[, 2],
  CI_Upper = parametric_matrix[, 1] + qnorm(0.975) * parametric_matrix[, 2],
  row.names = NULL
)

if (!is.null(summary_object$s.table)) {
  smooth_matrix <- as.matrix(summary_object$s.table)
  smooth_inference <- data.frame(
    Term = rownames(smooth_matrix),
    smooth_matrix,
    check.names = FALSE,
    row.names = NULL
  )
} else {
  smooth_inference <- data.frame()
}

concurvity_matrix <- mgcv::concurvity(gam_model, full = TRUE)
concurvity_output <- data.frame(
  Term = colnames(concurvity_matrix),
  t(concurvity_matrix),
  check.names = FALSE,
  row.names = NULL
)

write.csv(performance_comparison, "gam_test_performance.csv", row.names = FALSE)
write.csv(test_predictions, "gam_test_predictions.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per_governorate_metrics, "gam_per_governorate_metrics.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(parametric_inference, "gam_parametric_inference.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(smooth_inference, "gam_smooth_inference.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(concurvity_output, "gam_concurvity.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sm_stats, "gam_training_soil_moisture_scaling.csv", row.names = FALSE, fileEncoding = "UTF-8")

capture.output(summary(gam_model), file = "gam_model_summary.txt")
capture.output(gam_model, file = "gam_model_object.txt")
saveRDS(gam_model, file = "gam_durum_wheat_model.rds")

png("gam_smooth_effects.png", width = 1800, height = 1400, res = 180)
par(mfrow = c(2, 2))
plot(gam_model, pages = 1, shade = TRUE, residuals = FALSE, scale = 0)
dev.off()

png("gam_test_actual_vs_predicted.png", width = 1200, height = 1000, res = 160)
plot(
  actual,
  predicted_unclipped,
  xlab = "Observed production (1000 quintals)",
  ylab = "Predicted production (1000 quintals)",
  main = "GAM: observed vs predicted production",
  pch = 19
)
abline(a = 0, b = 1, lty = 2)
dev.off()

png("gam_diagnostics.png", width = 1600, height = 1200, res = 160)
par(mfrow = c(2, 2))
gam_check_text <- capture.output(gam.check(gam_model))
dev.off()
writeLines(gam_check_text, "gam_diagnostics.txt", useBytes = TRUE)

cat("GAM fitting completed.\n")
cat("Training observations:", nrow(train_data), "\n")
cat("Test observations:", nrow(test_data), "\n")
print(performance_comparison)
print(smooth_inference)
