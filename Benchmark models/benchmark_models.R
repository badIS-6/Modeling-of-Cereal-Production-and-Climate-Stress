project_dir <- "/home/badisz/Downloads/Modeling"
output_dir <- file.path(project_dir, "model_benchmark_results")

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prediction_files <- list(
  "Mixed Effects" = file.path(project_dir, "mixed effect _output", "mixed_effects_test_predictions.csv"),
  "Quantile Regression" = file.path(project_dir, "quantile regression _output", "quantile_regression_test_predictions.csv"),
  "GAM" = file.path(project_dir, "GAM _output", "gam_test_predictions.csv")
)

required_columns <- list(
  "Mixed Effects" = c("harvest_year", "governorate", "actual", "predicted_conditional"),
  "Quantile Regression" = c("harvest_year", "governorate", "actual_production_1000_qt", "predicted_tau_0_5"),
  "GAM" = c("harvest_year", "governorate", "actual", "predicted_unclipped")
)

read_prediction_file <- function(model_name, path) {
  if (!file.exists(path)) {
    stop("Required prediction file not found for ", model_name, ": ", path)
  }

  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  missing_columns <- setdiff(required_columns[[model_name]], names(x))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing columns in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (model_name == "Mixed Effects") {
    actual_name <- "actual"
    prediction_name <- "predicted_conditional"
  } else if (model_name == "Quantile Regression") {
    actual_name <- "actual_production_1000_qt"
    prediction_name <- "predicted_tau_0_5"
  } else {
    actual_name <- "actual"
    prediction_name <- "predicted_unclipped"
  }

  out <- data.frame(
    harvest_year = as.integer(x$harvest_year),
    governorate = trimws(as.character(x$governorate)),
    actual = as.numeric(x[[actual_name]]),
    prediction = as.numeric(x[[prediction_name]]),
    stringsAsFactors = FALSE
  )

  out$key <- paste(out$governorate, out$harvest_year, sep = "__")

  if (anyNA(out$harvest_year) || anyNA(out$governorate) ||
      any(out$governorate == "") || anyDuplicated(out$key)) {
    stop("Invalid or duplicated governorate-year keys in: ", path)
  }

  if (any(!is.finite(out$actual)) || any(!is.finite(out$prediction))) {
    stop("Non-finite actual or prediction values found in: ", path)
  }

  out
}

prediction_data <- lapply(names(prediction_files), function(model_name) {
  read_prediction_file(model_name, prediction_files[[model_name]])
})
names(prediction_data) <- names(prediction_files)

key_sets <- lapply(prediction_data, function(x) sort(x$key))
same_keys <- all(vapply(key_sets[-1], identical, logical(1), key_sets[[1]]))

if (!same_keys) {
  key_diagnostics <- data.frame(
    Model = names(prediction_data),
    Rows = vapply(prediction_data, nrow, integer(1)),
    Unique_Governorate_Year_Keys = vapply(
      prediction_data, function(x) length(unique(x$key)), integer(1)
    ),
    First_Year = vapply(prediction_data, function(x) min(x$harvest_year), integer(1)),
    Last_Year = vapply(prediction_data, function(x) max(x$harvest_year), integer(1))
  )
  write.csv(
    key_diagnostics,
    file.path(output_dir, "prediction_file_diagnostics.csv"),
    row.names = FALSE
  )
  stop(
    "The three prediction files do not contain identical governorate-year test observations. ",
    "See prediction_file_diagnostics.csv. No ranking was produced."
  )
}

reference <- prediction_data[[1]][order(prediction_data[[1]]$key), ]

for (model_name in names(prediction_data)[-1]) {
  current <- prediction_data[[model_name]][
    match(reference$key, prediction_data[[model_name]]$key),
  ]
  if (!isTRUE(all.equal(reference$actual, current$actual, tolerance = 1e-8))) {
    stop("Actual production values differ across prediction files for ", model_name, ".")
  }
}

if (nrow(reference) != 89L) {
  warning(
    "The common test set contains ", nrow(reference),
    " observations rather than the expected 89."
  )
}

if (!all(reference$harvest_year >= 2013L & reference$harvest_year <= 2018L)) {
  stop("The common prediction set includes observations outside 2013-2018.")
}

metric_rows <- lapply(names(prediction_data), function(model_name) {
  x <- prediction_data[[model_name]][
    match(reference$key, prediction_data[[model_name]]$key),
  ]

  errors <- x$actual - x$prediction
  sse <- sum(errors^2)
  sst <- sum((x$actual - mean(x$actual))^2)

  r2_value <- if (sst > 0) 1 - sse / sst else NA_real_
  rmse_value <- sqrt(mean(errors^2))
  mae_value <- mean(abs(errors))

  data.frame(
    Model = model_name,
    Prediction_Rule = switch(
      model_name,
      "Mixed Effects" = "Conditional unclipped prediction",
      "Quantile Regression" = "Median tau=0.5 unclipped prediction",
      "GAM" = "Unclipped prediction"
    ),
    N_Test = length(errors),
    First_Test_Year = min(x$harvest_year),
    Last_Test_Year = max(x$harvest_year),
    R2 = r2_value,
    RMSE = rmse_value,
    MAE = mae_value,
    Metrics_Valid = all(is.finite(c(r2_value, rmse_value, mae_value))),
    stringsAsFactors = FALSE
  )
})

benchmark <- do.call(rbind, metric_rows)
benchmark$R2_Rank <- rank(-benchmark$R2, ties.method = "average", na.last = "keep")
benchmark$RMSE_Rank <- rank(benchmark$RMSE, ties.method = "average", na.last = "keep")
benchmark$MAE_Rank <- rank(benchmark$MAE, ties.method = "average", na.last = "keep")
benchmark$Combined_Mean_Rank <- rowMeans(
  benchmark[, c("R2_Rank", "RMSE_Rank", "MAE_Rank")],
  na.rm = FALSE
)
benchmark <- benchmark[order(benchmark$Combined_Mean_Rank), ]

benchmark$H1_Candidate <- TRUE
benchmark$H1_Role <- c(
  "Primary full-H1 inferential candidate",
  "Primary partial-H1 inferential candidate",
  "Nonlinear partial-H1 supporting candidate"
)[match(benchmark$Model, c("Mixed Effects", "Quantile Regression", "GAM"))]

specification_audit <- data.frame(
  Model = c("Mixed Effects", "Quantile Regression", "GAM",
            "Random Forest", "Linear Regression"),
  H1_Candidate = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  H1_Role = c(
    "Primary full-H1 inferential candidate",
    "Primary partial-H1 inferential candidate",
    "Nonlinear partial-H1 supporting candidate",
    "Predictive machine-learning comparator; not a coefficient-based H1 test model",
    "Pooled predictive baseline; deliberately excludes spatial heterogeneity"
  ),
  H1_Temperature_Coverage = c(
    "Region-specific linear slopes",
    "Region-specific linear slopes",
    "Common nonlinear smooth; no region-specific temperature smooth",
    "Potential nonlinear effects/interactions; no explicit regional slope contrast",
    "Common linear slope only"
  ),
  H1_Soil_Moisture_Coverage = c(
    "Quadratic response with region-specific linear component",
    "Piecewise deficit/wet response common across regions",
    "Common nonlinear smooth; no region-specific soil-moisture smooth",
    "Potential nonlinear response; no explicit drought-branch parameterization",
    "Common linear raw soil-moisture slope only"
  ),
  Regional_Heterogeneity_Coverage = c(
    "Temperature and soil-moisture regional differences estimable",
    "Temperature regional difference estimable; soil-moisture regional difference not estimable",
    "Primary temperature and soil-moisture regional differences not estimable; NDVI and VHI have region-specific smooths",
    "Can explore regional patterns using model-based explanations; no direct coefficient test",
    "Not estimable: model is intentionally pooled"
  ),
  Specification_Warnings = c(
    "No explicit harvest_year trend or year effect; residual temporal structure should be assessed in a separate specification sensitivity analysis. Soil moisture is standardized within governorate using training-period mean and SD.",
    "No explicit harvest_year trend or year effect. No region:sm_deficit or region:sm_wet interaction, so regional soil-moisture sensitivity cannot be tested. Soil moisture is standardized within governorate.",
    "No explicit harvest_year trend or year effect. Temperature and soil-moisture smooths are shared across regions, so spatial H1 heterogeneity for these stressors is not tested. Soil moisture is standardized within governorate.",
    "No explicit harvest_year predictor. Variable importance is not a signed causal effect or a conventional coefficient hypothesis test.",
    "No governorate or region effect by deliberate baseline design; it cannot validate spatial heterogeneity. No explicit harvest_year effect."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  benchmark,
  file.path(output_dir, "benchmark_H1_candidates.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  specification_audit,
  file.path(output_dir, "H1_model_specification_audit.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  do.call(rbind, lapply(names(prediction_data), function(model_name) {
    x <- prediction_data[[model_name]]
    data.frame(
      Model = model_name,
      File = prediction_files[[model_name]],
      Rows = nrow(x),
      Unique_Governorates = length(unique(x$governorate)),
      First_Year = min(x$harvest_year),
      Last_Year = max(x$harvest_year),
      Prediction_Column = switch(
        model_name,
        "Mixed Effects" = "predicted_conditional",
        "Quantile Regression" = "predicted_tau_0_5",
        "GAM" = "predicted_unclipped"
      ),
      Actual_Column = switch(
        model_name,
        "Mixed Effects" = "actual",
        "Quantile Regression" = "actual_production_1000_qt",
        "GAM" = "actual"
      ),
      stringsAsFactors = FALSE
    )
  })),
  file.path(output_dir, "prediction_file_diagnostics.csv"),
  row.names = FALSE
)

cat("\nH1 MODEL SHORTLIST AND SPECIFICATION AUDIT\n")
cat("Project directory:", project_dir, "\n")
cat("Common test observations:", nrow(reference), "\n")
cat("Test years:", min(reference$harvest_year), "-", max(reference$harvest_year), "\n")
cat("\nMETRICS COMPUTED DIRECTLY FROM UNCLIPPED PREDICTIONS\n")
print(benchmark, row.names = FALSE)

cat("\nMODEL-SPECIFICATION AUDIT\n")
print(specification_audit, row.names = FALSE)

cat("\nSaved files in:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
