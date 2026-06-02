#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop("Usage: Rscript combine_rds_to_excel.R <input_directory> <output_excel_file>")
}

input_dir <- args[1]
output_excel <- args[2]

rds_files <- list.files(
  path = input_dir,
  pattern = "\\.RDS$",
  full.names = TRUE
)

if(length(rds_files) == 0){
  stop("No RDS files found in the specified directory.")
}

wb <- createWorkbook()

addUniqueWorksheet <- function(wb, sheet_name, data) {

  base_name <- sheet_name
  suffix <- 1

  while(base_name %in% names(wb)){
    base_name <- paste0(sheet_name, "_", suffix)
    suffix <- suffix + 1
  }

  if(nchar(base_name) > 31) {
    base_name <- substr(base_name, 1, 31)
  }

  addWorksheet(wb, base_name)

  writeDataTable(
    wb,
    sheet = base_name,
    x = data,
    tableStyle = "TableStyleLight9"
  )
}

aggregate_population_counts <- function(df) {

  required_cols <- c(
    "COSMIC_Sample_Count",
    "COSMIC_Targeted_Sample_Count"
  )

  missing_cols <- setdiff(required_cols, colnames(df))

  if(length(missing_cols) > 0){
    warning(
      paste(
        "Missing COSMIC columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
    return(df)
  }

  df$COSMIC_Sample_Count <- suppressWarnings(
    as.numeric(df$COSMIC_Sample_Count)
  )

  df$COSMIC_Targeted_Sample_Count <- suppressWarnings(
    as.numeric(df$COSMIC_Targeted_Sample_Count)
  )

  df$COSMIC_Sample_Count[
    is.na(df$COSMIC_Sample_Count)
  ] <- 0

  df$COSMIC_Targeted_Sample_Count[
    is.na(df$COSMIC_Targeted_Sample_Count)
  ] <- 0

  df$COSMIC_Total_Sample_Count <-
    df$COSMIC_Sample_Count +
    df$COSMIC_Targeted_Sample_Count
  
  if ("GENIE_Sample_Count" %in% colnames(df)) {
    df$GENIE_Sample_Count <- suppressWarnings(as.numeric(df$GENIE_Sample_Count))
    df$GENIE_Sample_Count[is.na(df$GENIE_Sample_Count)] <- 0
  }
  return(df)
}

all_variants <- list()

for(rds_file in rds_files){

  message("Processing file: ", rds_file)

  sample_data <- readRDS(rds_file)

  sample_data <- aggregate_population_counts(sample_data)

  sample_name <- tools::file_path_sans_ext(
    basename(rds_file)
  )

  if(nrow(sample_data) == 0){

    message(
      "  Warning: ",
      sample_name,
      " is empty. Creating placeholder sheet."
    )

    sample_data <- data.frame(
      Sample_ID = sample_name,
      Note = "No variants found",
      stringsAsFactors = FALSE
    )

  } else {

    if(!"Sample_ID" %in% colnames(sample_data)){
      sample_data$Sample_ID <- sample_name
    }
  }

  addUniqueWorksheet(
    wb,
    sample_name,
    sample_data
  )

  if(!"Note" %in% colnames(sample_data)){
    all_variants[[sample_name]] <- sample_data
  }
}

combined_variants <- bind_rows(all_variants)

addWorksheet(wb, "Summary")

writeDataTable(
  wb,
  sheet = "Summary",
  x = combined_variants,
  tableStyle = "TableStyleLight9"
)

saveWorkbook(
  wb,
  output_excel,
  overwrite = TRUE
)

cat(
  "Excel workbook created successfully at",
  output_excel,
  "\n"
)
