library(magrittr)

read_blast <- function(blast_file, include_empty){
  # Define columns
  column_names <- c("query_accession", "subject_accession", "percent_identity", "alignment_length", "mismatches", "gap_opens", "query_start", "query_end", "subject_start", "subject_end", "e_value", "bit_score")
  column_types <- "ccdiiiiiiidd"
  
  # Define sample name
  sample_file <- basename(path = blast_file)
  sample <- stringr::str_remove(string = sample_file, pattern = "\\.blast")

  # Read blast results
  blast <- readr::read_tsv(file = blast_file, col_names = column_names, col_types = column_types)
  
  # Check whether to add empty samples
  add_empty <- nrow(blast) == 0 & include_empty
  if(add_empty){
    
    # Generate empty row
    raw_row <- rep(NA, length(column_names))
    row_list <- as.list(raw_row)
    names(row_list) <- column_names
    
    row <- dplyr::as_tibble(row_list)
    
    # Implement empty row
    blast <- dplyr::add_row(blast, row)
  }
  
  # Add sample name to matches
  dplyr::mutate(blast, sample = sample)
}


read_index <- function(index_file){
  # Define columns
  column_names <- c("query", "length", "offset", "linebases", "linewidth")
  column_types <- "ciiii"
  
  # Read query index 
  readr::read_tsv(file = index_file, col_names = column_names, col_types = column_types)
}


calculate_coverage <- function(blast, index){
  # Merge Blast results and query indexes
  blast_indexed <- dplyr::left_join(x = blast, y = index, by = c("query_accession" = "query"))
  
  # Calculate query coverage
  dplyr::mutate(
    blast_indexed,
    percent_coverage = alignment_length / length * 100
  )
}


cleanup_blast_columns <- function(blast){
  # Wrangle NA values when applicable
  blast_wrangled <- dplyr::mutate(
    blast,
    percent_identity = dplyr::case_when(
      is.na(percent_identity) ~ 0,
      TRUE ~ percent_identity
    ),
    percent_coverage = dplyr::case_when(
      is.na(percent_coverage) ~ 0,
      TRUE ~ percent_coverage
    )
  )
  
  # Define column selection and order
  columns_selected <- c(
    "sample", "query_accession", "percent_identity", "percent_coverage",
    "alignment_length", "mismatches", "gap_opens", "bit_score", 
    "query_length" = "length", "query_start", "query_end",
    "subject_accession", "subject_start", "subject_end"
  )
  
  dplyr::select(blast_wrangled, dplyr::all_of(columns_selected))
}

import_blast <- function(blast_files, index_file, include_empty){
  blast <- purrr::map_dfr(.x = blast_files, .f = read_blast, include_empty = include_empty)
  index <- read_index(index_file)
  blast_coverage <- calculate_coverage(blast, index)
  
  cleanup_blast_columns(blast_coverage)
}


summarise_results <- function(blast, threshold, top_only){
  
  # Convert threshold to numeric
  if (isFALSE(threshold))
    threshold = 0
  
  # Group by samples and query
  blast_by_query <- dplyr::group_by(blast, sample, query_accession)
  
  # Check whether samples passes thresholds and are top results
  blast_checkup <- dplyr::filter(
    blast_by_query,
    
    # Filter by threshold
    percent_coverage >= threshold & percent_identity >= threshold,
    
    # Filter by top results (if applicable)
    dplyr::case_when(
      
      # For top results, select max coverage and identity
      top_only ~ percent_coverage == max(percent_coverage) & percent_identity == max(percent_identity),
      
      # Otherwise include all
      TRUE ~ TRUE
    )
  )
  
  # Remove groups
  dplyr::ungroup(blast_checkup)
}


write_junctions <- function(junctions, junctions_file){
  
  file_type <- tools::file_ext(junctions_file)
  
  success <- FALSE
  if (file_type == "xlsx"){
    writexl::write_xlsx(x = junctions, path = junctions_file)
    success <- TRUE
  } else if (file_type == "csv"){
    readr::write_csv(x = junctions, file = junctions_file)
    success <- TRUE
  } else if (file_type == "tsv"){
    readr::write_tsv(x = junctions, file = junctions_file)
    success <- TRUE
  }
  
  return(success)
}


summarize_junctions <- function(blast_files, index_file, threshold, junctions_file, include_empty, top_only){
  blast <- import_blast(blast_files, index_file, include_empty)
  junctions <- summarise_results(blast, threshold, top_only)
  
  write_successfull <- write_junctions(junctions, junctions_file)
  
  if (write_successfull){
    message(glue::glue("Junctions results successfully written to: {junctions_file}"))
  } else {
    stop(glue::glue("Junctions results file could not be written to: {junctions_file}"))
  }
}


blast_files <- snakemake@input$blast
index_file <- snakemake@input$index
threshold <- snakemake@params$threshold
top_only <- snakemake@params$top_only
include_empty <- snakemake@params$include_empty
debug <- snakemake@params$debug
junctions_file <- snakemake@output$junctions

no_blasts <- length(blast_files) == 0
if (no_blasts)
  stop("No blast files detected. Are your sure you pointed to the correct input directory?")

if (debug){
  debug_file <- file.path(dirname(junctions_file), "summarise_junctions.RData")
  message(glue::glue("Debug mode enabled, writing snakemake object to: {debug_file}"))
  saveRDS(snakemake, debug_file)
}

summarize_junctions(blast_files, index_file, threshold, junctions_file, include_empty, top_only)
message("Done!")