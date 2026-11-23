options(stringsAsFactors = FALSE)

script_file <- ""
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_file <- sub("^--file=", "", script_arg[1])
}
if (!nzchar(script_file)) {
  frames <- sys.frames()
  for (fr in rev(frames)) {
    if (!is.null(fr$ofile) && nzchar(fr$ofile)) {
      script_file <- fr$ofile
      break
    }
  }
}

root_override <- Sys.getenv("MODELLING_DIR", unset = "")
if (nzchar(root_override)) {
  root_dir <- normalizePath(root_override, mustWork = TRUE)
} else if (nzchar(script_file)) {
  root_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
} else {
  root_dir <- normalizePath(getwd(), mustWork = TRUE)
}

normalise_key <- function(x) {
  x <- gsub("²", "2", as.character(x), fixed = TRUE)
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

metric_kind <- function(x) {
  if (length(x) == 0 || is.na(x[1])) return(NA_character_)
  z <- normalise_key(x[1])
  if (grepl("r2|rsquared|coefficientofdetermination", z)) return("R2")
  if (grepl("rmse|rootmeansquarederror", z)) return("RMSE")
  if (grepl("mae|meanabsoluteerror", z)) return("MAE")
  NA_character_
}

parse_number <- function(x, metric) {
  original <- trimws(as.character(x))
  had_percent <- grepl("%", original, fixed = TRUE)
  cleaned <- gsub("%", "", original, fixed = TRUE)
  cleaned <- gsub(",", "", cleaned, fixed = TRUE)
  value <- suppressWarnings(as.numeric(cleaned))
  if (length(value) == 0 || is.na(value[1])) return(NA_real_)
  value <- value[1]
  if (identical(metric, "R2") && had_percent) value <- value / 100
  value
}

read_model_metrics <- function(path) {
  dat <- tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) {
    return(list(R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_, Scope = "Unreadable or empty CSV"))
  }

  keys <- normalise_key(names(dat))
  metric_col <- which(grepl("metric|measure|statistic", keys))
  value_col <- which(grepl("^value$|score|metricvalue|result|estimate", keys))
  scope_text <- "Unspecified; assumed to be validation metrics"
  values <- c(R2 = NA_real_, RMSE = NA_real_, MAE = NA_real_)

  if (length(metric_col) > 0) {
    mc <- metric_col[1]
    if (length(value_col) > 0) {
      preferred_value_col <- value_col[grepl("test|valid|holdout", keys[value_col])]
      vc <- if (length(preferred_value_col) > 0) preferred_value_col[1] else value_col[1]
    } else {
      numeric_candidates <- which(vapply(dat, function(x) {
        any(is.finite(suppressWarnings(as.numeric(gsub(",", "", as.character(x))))))
      }, logical(1)))
      numeric_candidates <- setdiff(numeric_candidates, mc)
      preferred_numeric_col <- numeric_candidates[grepl("test|valid|holdout", keys[numeric_candidates])]
      vc <- if (length(preferred_numeric_col) > 0) preferred_numeric_col[1] else if (length(numeric_candidates) > 0) numeric_candidates[1] else NA_integer_
    }

    if (!is.na(vc)) {
      metric_labels <- as.character(dat[[mc]])
      split_cols <- which(grepl("split|set|subset|partition|sample|stage|dataset|type", keys))
      row_labels <- metric_labels
      if (length(split_cols) > 0) {
        row_labels <- paste(row_labels, apply(dat[, split_cols, drop = FALSE], 1, paste, collapse = " "))
      }
      if (any(grepl("test|valid|holdout", tolower(row_labels)))) {
        scope_text <- "Validation/Test rows detected"
      }
      for (m in names(values)) {
        idx <- which(vapply(metric_labels, function(v) identical(metric_kind(v), m), logical(1)))
        if (length(idx) == 0) {
          idx <- which(vapply(row_labels, function(v) identical(metric_kind(v), m), logical(1)))
        }
        if (length(idx) > 0) {
          preferred <- idx[grepl("test|valid|holdout", tolower(row_labels[idx]))]
          use_idx <- if (length(preferred) > 0) preferred[1] else idx[1]
          values[m] <- parse_number(dat[[vc]][use_idx], m)
        }
      }
    }
  }

  if (all(is.na(values))) {
    col_kinds <- vapply(names(dat), metric_kind, character(1))
    for (m in names(values)) {
      idx <- which(col_kinds == m)
      if (length(idx) > 0) {
        preferred <- idx[grepl("test|valid|holdout", keys[idx])]
        use_col <- if (length(preferred) > 0) preferred[1] else idx[1]
        row_idx <- 1L
        scope_cols <- which(grepl("split|set|subset|partition|sample|stage|dataset|type", keys))
        if (nrow(dat) > 1 && length(scope_cols) > 0) {
          matching_rows <- which(apply(dat[, scope_cols, drop = FALSE], 1, function(row) {
            grepl("test|valid|holdout", tolower(paste(row, collapse = " ")))
          }))
          if (length(matching_rows) > 0) row_idx <- matching_rows[1]
        }
        values[m] <- parse_number(dat[[use_col]][row_idx], m)
        if (grepl("test|valid|holdout", keys[use_col])) {
          scope_text <- "Validation/Test columns detected"
        } else if (nrow(dat) > 1 && length(scope_cols) > 0) {
          scope_text <- "Validation/Test row detected"
        }
      }
    }
  }

  if (all(is.na(values)) && ncol(dat) >= 2 && nrow(dat) >= 3) {
    labels <- as.character(dat[[1]])
    for (m in names(values)) {
      idx <- which(vapply(labels, function(v) identical(metric_kind(v), m), logical(1)))
      if (length(idx) > 0) values[m] <- parse_number(dat[[2]][idx[1]], m)
    }
  }

  list(R2 = unname(values["R2"]), RMSE = unname(values["RMSE"]), MAE = unname(values["MAE"]), Scope = scope_text)
}

classify_model <- function(folder_name) {
  z <- normalise_key(folder_name)
  if (grepl("mixedeffect|mixedmodel|nlme|lme", z)) return("Mixed Effects")
  if (grepl("quantile", z)) return("Quantile Regression")
  if (grepl("generalizedadditive|gam", z)) return("GAM")
  if (grepl("randomforest|ranger", z)) return("Random Forest")
  if (grepl("linearregression|ols|linear", z)) return("Linear Regression")
  folder_name
}

metric_files <- list.files(
  path = root_dir,
  pattern = "^R2_RMSE_MAE\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(metric_files) == 0) {
  stop(paste0("No R2_RMSE_MAE.csv files were found under: ", root_dir,
              "\nPlace this script in the modelling folder containing the model output folders, or set MODELLING_DIR."))
}

records <- lapply(metric_files, function(path) {
  folder <- basename(dirname(path))
  model_type <- classify_model(folder)
  metrics <- read_model_metrics(path)
  h1_primary <- model_type %in% c("Mixed Effects", "Quantile Regression")
  h1_role <- if (model_type == "Mixed Effects") {
    "Primary H1 candidate: hierarchical effects and regional climate slopes/interactions must be tested"
  } else if (model_type == "Quantile Regression") {
    "Primary H1 candidate: regional temperature and soil-moisture interactions must be tested"
  } else if (model_type == "GAM") {
    "Supporting nonlinear model; regional temperature and soil-moisture smooths are needed for full H1 testing"
  } else if (model_type == "Random Forest") {
    "Predictive benchmark; directional H1 inference requires separate explainability/inference analysis"
  } else if (model_type == "Linear Regression") {
    "Baseline only; does not address the stated spatial dependence without regional interactions"
  } else {
    "Model type not classified for primary H1 inference"
  }
  data.frame(
    Model = model_type,
    Output_Folder = folder,
    Metrics_File = normalizePath(path, winslash = "/", mustWork = TRUE),
    R2 = metrics$R2,
    RMSE = metrics$RMSE,
    MAE = metrics$MAE,
    Metric_Scope = metrics$Scope,
    H1_Primary_Candidate = h1_primary,
    H1_Role = h1_role,
    stringsAsFactors = FALSE
  )
})

benchmark <- do.call(rbind, records)
benchmark$Metrics_Readable <- with(benchmark, is.finite(R2) & is.finite(RMSE) & is.finite(MAE) & RMSE >= 0 & MAE >= 0)
benchmark$Validation_Metrics_Valid <- with(benchmark, Metrics_Readable & R2 > 0)
benchmark$R2_Rank <- NA_real_
benchmark$RMSE_Rank <- NA_real_
benchmark$MAE_Rank <- NA_real_
benchmark$Combined_Rank <- NA_real_

scorable <- which(benchmark$Metrics_Readable)
if (length(scorable) > 0) {
  benchmark$R2_Rank[scorable] <- rank(-benchmark$R2[scorable], ties.method = "min")
  benchmark$RMSE_Rank[scorable] <- rank(benchmark$RMSE[scorable], ties.method = "min")
  benchmark$MAE_Rank[scorable] <- rank(benchmark$MAE[scorable], ties.method = "min")
  benchmark$Combined_Rank[scorable] <- rowMeans(benchmark[scorable, c("R2_Rank", "RMSE_Rank", "MAE_Rank"), drop = FALSE])
}

benchmark <- benchmark[order(is.na(benchmark$Combined_Rank), benchmark$Combined_Rank, -ifelse(is.finite(benchmark$R2), benchmark$R2, -Inf), benchmark$RMSE, benchmark$MAE), ]
rownames(benchmark) <- NULL

candidate_pool <- benchmark[benchmark$H1_Primary_Candidate & benchmark$Validation_Metrics_Valid, , drop = FALSE]
if (nrow(candidate_pool) > 0) {
  candidate_pool <- candidate_pool[order(candidate_pool$Combined_Rank, -candidate_pool$R2, candidate_pool$RMSE, candidate_pool$MAE), , drop = FALSE]
  selected <- head(candidate_pool, 2)
} else {
  selected <- benchmark[0, , drop = FALSE]
}

output_dir <- file.path(root_dir, "model_benchmark_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(benchmark, file.path(output_dir, "benchmark_all_models.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")
write.csv(selected, file.path(output_dir, "selected_models_for_H1.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

if (nrow(selected) > 0) {
  checklist <- expand.grid(
    Model = selected$Model,
    H1_Component = c("Elevated December temperature", "Autumn soil-moisture deficit", "Regional sensitivity"),
    stringsAsFactors = FALSE
  )
  checklist$Expected_Evidence <- ifelse(
    checklist$H1_Component == "Elevated December temperature",
    "Negative temperature effect with p < 0.05",
    ifelse(
      checklist$H1_Component == "Autumn soil-moisture deficit",
      "Evidence that low autumn soil moisture reduces production; test a nonlinear low-moisture response where specified",
      "Significant region-by-temperature/soil-moisture interaction, regional slope, or regional contrast"
    )
  )
  checklist$Estimate_or_Test_Statistic <- NA_character_
  checklist$P_value <- NA_real_
  checklist$Evidence_File <- NA_character_
  checklist$Inference_Status <- "Pending: enter the selected model's coefficient/smooth/interaction test result"
} else {
  checklist <- data.frame(
    Model = character(0),
    H1_Component = character(0),
    Expected_Evidence = character(0),
    Estimate_or_Test_Statistic = character(0),
    P_value = numeric(0),
    Evidence_File = character(0),
    Inference_Status = character(0)
  )
}

write.csv(checklist, file.path(output_dir, "H1_inference_checklist.csv"), row.names = FALSE, na = "", fileEncoding = "UTF-8")

report_lines <- c(
  paste0("Model benchmarking directory: ", root_dir),
  paste0("Metric files read: ", length(metric_files)),
  "Ranking: equal-weight average rank of validation R2 (higher is better), RMSE (lower is better), and MAE (lower is better).",
  "A model is validation-eligible only when all three metrics are finite and R2 > 0.",
  "Primary H1 shortlist is restricted to Mixed Effects and Quantile Regression.",
  "At most two models are selected.",
  "Metric ranking does not itself validate H1; the H1 inference checklist requires coefficient, nonlinear-response, and regional-interaction evidence."
)

if (nrow(selected) > 0) {
  report_lines <- c(report_lines, "", "Selected H1 candidates:")
  for (i in seq_len(nrow(selected))) {
    report_lines <- c(report_lines, paste0(
      i, ". ", selected$Model[i],
      " | R2 = ", signif(selected$R2[i], 5),
      " | RMSE = ", signif(selected$RMSE[i], 5),
      " | MAE = ", signif(selected$MAE[i], 5),
      " | combined rank = ", signif(selected$Combined_Rank[i], 4)
    ))
  }
} else {
  report_lines <- c(report_lines, "", "No H1-primary model met the validation eligibility rule. Check the metric CSV column layout and whether Mixed Effects or Quantile Regression has R2 > 0.")
}

writeLines(report_lines, file.path(output_dir, "benchmark_summary.txt"), useBytes = TRUE)

cat(paste(report_lines, collapse = "\n"), "\n")
cat("\nFiles saved in: ", output_dir, "\n", sep = "")
