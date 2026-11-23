
install.packages("lme4")
install.packages("lmerTest")

suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))

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

contrasts(data$region) <- contr.sum(2)

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

contrasts(train_data$region) <- contr.sum(2)
contrasts(test_data$region) <- contr.sum(2)

mixed_formula_random_intercept <- production_1000_qt ~
  region * (temp_dec_c + sm_z + ndvi_mean) +
  I(sm_z^2) +
  rain_jan_mar_mm +
  vhi +
  (1 | governorate)

mixed_formula_random_slope <- production_1000_qt ~
  region * (temp_dec_c + sm_z + ndvi_mean) +
  I(sm_z^2) +
  rain_jan_mar_mm +
  vhi +
  (1 + sm_z || governorate)

model_control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000),
  check.rankX = "message+drop.cols"
)

model_random_intercept_ml <- lmerTest::lmer(
  formula = mixed_formula_random_intercept,
  data = train_data,
  REML = FALSE,
  control = model_control,
  na.action = na.fail
)

model_random_slope_ml <- tryCatch(
  lmerTest::lmer(
    formula = mixed_formula_random_slope,
    data = train_data,
    REML = FALSE,
    control = model_control,
    na.action = na.fail
  ),
  error = function(e) NULL
)

random_intercept_singular <- lme4::isSingular(model_random_intercept_ml, tol = 1e-4)

model_converged <- function(model) {
  if (is.null(model)) return(FALSE)
  optimizer_code <- model@optinfo$conv$opt
  lme4_messages <- model@optinfo$conv$lme4$messages
  optimizer_ok <- is.null(optimizer_code) || isTRUE(as.numeric(optimizer_code)[1] == 0)
  optimizer_ok && is.null(lme4_messages)
}

random_intercept_converged <- model_converged(model_random_intercept_ml)

if (is.null(model_random_slope_ml)) {
  random_slope_singular <- NA
  random_slope_converged <- FALSE
  likelihood_ratio <- NA_real_
  boundary_p_value <- NA_real_
  random_slope_aic <- NA_real_
  random_slope_bic <- NA_real_
  random_slope_loglik <- NA_real_
  random_slope_selected <- FALSE
  selection_reason <- "Random-slope model could not be fitted; selected random-intercept model."
} else {
  random_slope_singular <- lme4::isSingular(model_random_slope_ml, tol = 1e-4)
  random_slope_converged <- model_converged(model_random_slope_ml)
  likelihood_ratio <- max(0, 2 * (as.numeric(logLik(model_random_slope_ml)) - as.numeric(logLik(model_random_intercept_ml))))
  boundary_p_value <- if (likelihood_ratio <= .Machine$double.eps) 1 else 0.5 * (1 - pchisq(likelihood_ratio, df = 1))
  random_slope_aic <- AIC(model_random_slope_ml)
  random_slope_bic <- BIC(model_random_slope_ml)
  random_slope_loglik <- as.numeric(logLik(model_random_slope_ml))
  random_slope_selected <- random_intercept_converged &&
    random_slope_converged &&
    !random_slope_singular &&
    is.finite(boundary_p_value) &&
    boundary_p_value < 0.05 &&
    AIC(model_random_slope_ml) < AIC(model_random_intercept_ml)
  selection_reason <- if (random_slope_selected) {
    "Selected the uncorrelated soil-moisture random-slope model because it converged, was non-singular, improved AIC, and passed the boundary-adjusted likelihood-ratio comparison."
  } else {
    "Selected the random-intercept model; the random-slope model did not meet all convergence, singularity, AIC, and likelihood-ratio criteria."
  }
}

model_comparison <- data.frame(
  Model = c("Random intercept", "Random intercept plus soil-moisture random slope"),
  LogLik = c(as.numeric(logLik(model_random_intercept_ml)), random_slope_loglik),
  AIC = c(AIC(model_random_intercept_ml), random_slope_aic),
  BIC = c(BIC(model_random_intercept_ml), random_slope_bic),
  Singular = c(random_intercept_singular, random_slope_singular),
  Converged = c(random_intercept_converged, random_slope_converged),
  LRT_vs_random_intercept = c(NA_real_, likelihood_ratio),
  Boundary_adjusted_p = c(NA_real_, boundary_p_value),
  Selected = c(!random_slope_selected, random_slope_selected),
  stringsAsFactors = FALSE
)

selected_formula <- if (random_slope_selected) mixed_formula_random_slope else mixed_formula_random_intercept
selected_model_name <- if (random_slope_selected) "Random intercept plus soil-moisture random slope" else "Random intercept"

mixed_model <- lmerTest::lmer(
  formula = selected_formula,
  data = train_data,
  REML = TRUE,
  control = model_control,
  na.action = na.fail
)

if (lme4::isSingular(mixed_model, tol = 1e-4)) {
  warning("The selected REML model is singular. Interpret random-effect variance estimates cautiously.")
}

prediction_conditional <- as.numeric(predict(
  mixed_model,
  newdata = test_data,
  re.form = NULL,
  allow.new.levels = FALSE
))

prediction_fixed_only <- as.numeric(predict(
  mixed_model,
  newdata = test_data,
  re.form = NA,
  allow.new.levels = FALSE
))

prediction_conditional_clipped <- pmax(prediction_conditional, 0)
prediction_fixed_only_clipped <- pmax(prediction_fixed_only, 0)
actual <- test_data$production_1000_qt

calculate_metrics <- function(actual_values, predicted_values, prediction_type) {
  sse <- sum((actual_values - predicted_values)^2)
  sst <- sum((actual_values - mean(actual_values))^2)
  data.frame(
    Prediction = prediction_type,
    R2 = if (sst > 0) 1 - sse / sst else NA_real_,
    RMSE = sqrt(mean((actual_values - predicted_values)^2)),
    MAE = mean(abs(actual_values - predicted_values))
  )
}

performance_comparison <- rbind(
  calculate_metrics(actual, prediction_conditional, "Conditional - unclipped"),
  calculate_metrics(actual, prediction_conditional_clipped, "Conditional - clipped at zero"),
  calculate_metrics(actual, prediction_fixed_only, "Fixed effects only - unclipped"),
  calculate_metrics(actual, prediction_fixed_only_clipped, "Fixed effects only - clipped at zero")
)

performance_comparison$Selected_model <- selected_model_name
performance_comparison <- performance_comparison[order(performance_comparison$RMSE), ]
rownames(performance_comparison) <- NULL

test_predictions <- data.frame(
  harvest_year = test_data$harvest_year,
  governorate = as.character(test_data$governorate),
  region = as.character(test_data$region),
  actual = actual,
  predicted_conditional = prediction_conditional,
  predicted_conditional_clipped = prediction_conditional_clipped,
  predicted_fixed_only = prediction_fixed_only,
  predicted_fixed_only_clipped = prediction_fixed_only_clipped,
  residual_conditional = actual - prediction_conditional,
  residual_fixed_only = actual - prediction_fixed_only
)

governorate_indices <- split(seq_len(nrow(test_predictions)), test_predictions$governorate)

per_governorate_metrics <- do.call(
  rbind,
  lapply(names(governorate_indices), function(governorate_name) {
    selected <- test_predictions[governorate_indices[[governorate_name]], ]
    data.frame(
      governorate = governorate_name,
      n_test = nrow(selected),
      RMSE_conditional = sqrt(mean((selected$actual - selected$predicted_conditional)^2)),
      MAE_conditional = mean(abs(selected$actual - selected$predicted_conditional)),
      RMSE_fixed_only = sqrt(mean((selected$actual - selected$predicted_fixed_only)^2)),
      MAE_fixed_only = mean(abs(selected$actual - selected$predicted_fixed_only)),
      stringsAsFactors = FALSE
    )
  })
)

coefficient_matrix <- as.matrix(coef(summary(mixed_model)))
coefficient_df <- data.frame(
  Term = rownames(coefficient_matrix),
  Estimate = coefficient_matrix[, "Estimate"],
  Std_Error = coefficient_matrix[, "Std. Error"],
  df = coefficient_matrix[, "df"],
  t_value = coefficient_matrix[, "t value"],
  p_value = coefficient_matrix[, "Pr(>|t|)"],
  CI_Lower = coefficient_matrix[, "Estimate"] - qt(0.975, df = coefficient_matrix[, "df"]) * coefficient_matrix[, "Std. Error"],
  CI_Upper = coefficient_matrix[, "Estimate"] + qt(0.975, df = coefficient_matrix[, "df"]) * coefficient_matrix[, "Std. Error"],
  row.names = NULL
)

fixed_effect_anova <- as.data.frame(stats::anova(mixed_model, type = 3, ddf = "Satterthwaite"))
fixed_effect_anova$Term <- rownames(fixed_effect_anova)
rownames(fixed_effect_anova) <- NULL
fixed_effect_anova <- fixed_effect_anova[, c("Term", setdiff(names(fixed_effect_anova), "Term")), drop = FALSE]

random_effects <- as.data.frame(lme4::ranef(mixed_model, condVar = TRUE)$governorate)
random_effects$governorate <- rownames(random_effects)
rownames(random_effects) <- NULL
random_effects <- random_effects[, c("governorate", setdiff(names(random_effects), "governorate")), drop = FALSE]

random_effect_variances <- as.data.frame(lme4::VarCorr(mixed_model))

training_residuals <- residuals(mixed_model)
training_fitted <- fitted(mixed_model)
training_diagnostics <- data.frame(
  fitted = training_fitted,
  residual = training_residuals,
  standardized_residual = as.numeric(scale(training_residuals))
)

write.csv(model_comparison, "mixed_effects_random_effect_structure_comparison.csv", row.names = FALSE)
write.csv(performance_comparison, "mixed_effects_test_performance.csv", row.names = FALSE)
write.csv(test_predictions, "mixed_effects_test_predictions.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per_governorate_metrics, "mixed_effects_per_governorate_metrics.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(coefficient_df, "mixed_effects_fixed_effect_inference.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fixed_effect_anova, "mixed_effects_type3_anova.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(random_effects, "mixed_effects_governorate_random_effects.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(random_effect_variances, "mixed_effects_random_effect_variances.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sm_stats, "mixed_effects_training_soil_moisture_scaling.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(training_diagnostics, "mixed_effects_training_diagnostics.csv", row.names = FALSE)

capture.output(summary(mixed_model), file = "mixed_effects_model_summary.txt")
capture.output(lme4::isSingular(mixed_model, tol = 1e-4), file = "mixed_effects_singularity_check.txt")
capture.output(selection_reason, file = "mixed_effects_random_slope_selection.txt")
saveRDS(mixed_model, file = "mixed_effects_durum_wheat_model.rds")

png("mixed_effects_test_actual_vs_predicted.png", width = 1200, height = 1000, res = 160)
plot(
  actual,
  prediction_conditional,
  xlab = "Observed production (1000 quintals)",
  ylab = "Predicted production (1000 quintals)",
  main = "Mixed-effects model: observed vs predicted production",
  pch = 19
)
abline(a = 0, b = 1, lty = 2)
dev.off()

png("mixed_effects_residual_diagnostics.png", width = 1600, height = 1200, res = 160)
par(mfrow = c(2, 2))
plot(training_fitted, training_residuals, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs fitted", pch = 19)
abline(h = 0, lty = 2)
qqnorm(training_residuals, main = "Normal Q-Q plot of residuals")
qqline(training_residuals, lty = 2)
hist(training_residuals, main = "Training residual distribution", xlab = "Residuals")
plot(train_data$harvest_year, training_residuals, xlab = "Harvest year", ylab = "Residuals", main = "Residuals over time", pch = 19)
abline(h = 0, lty = 2)
dev.off()

cat("Mixed-effects model fitting completed.\n")
cat("Selected model:", selected_model_name, "\n")
cat("Training observations:", nrow(train_data), "\n")
cat("Test observations:", nrow(test_data), "\n")
cat("Random slope selection:", selection_reason, "\n")
print(performance_comparison)
print(coefficient_df)
print(fixed_effect_anova)
