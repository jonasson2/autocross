#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

format_num <- function(x) {
  sprintf("%.2f", x)
}

estimate_path <- function(path) {
  data <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
  data <- data[order(data$year, decreasing = TRUE), ]
  result <- autocross::estimate_CI(2000 - data$year, data$SST, data$T)

  data.frame(
    station = tools::file_path_sans_ext(basename(path)),
    n = result$n,
    r = format_num(result$r),
    `95%-CI` = sprintf("[%s, %s]", format_num(result$ci[1]), format_num(result$ci[2])),
    taux = format_num(result$taux),
    tauy = format_num(result$tauy),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

if (length(args) == 0) {
  stop(
    "usage: Rscript pearsonT3-table.R data-file [data-file ...]\n",
    "example: Rscript pearsonT3-table.R data/Stórhöfði.txt\n",
    "example: Rscript pearsonT3-table.R data/*.txt",
    call. = FALSE
  )
}

is_data_file <- function(path) {
  if (!file.exists(path)) {
    stop("Data file not found: ", path, call. = FALSE)
  }
  header <- readLines(path, n = 1, warn = FALSE)
  fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
  all(c("year", "SST", "T") %in% fields)
}

data_files <- args[vapply(args, is_data_file, logical(1))]
if (length(data_files) == 0) {
  stop("No input files with columns year, SST, and T were found.", call. = FALSE)
}

estimate_file <- function(path) {
  if (!file.exists(path)) {
    stop("Data file not found: ", path, call. = FALSE)
  }
  suppressMessages(estimate_path(path))
}

columns <- c("station", "n", "r", "95%-CI", "taux", "tauy")
results <- do.call(rbind, lapply(data_files, estimate_file))
results <- results[, columns]

display_width <- function(x) {
  nchar(x, type = "width", allowNA = FALSE, keepNA = FALSE)
}

widths <- vapply(columns, function(column) {
  max(display_width(c(column, as.character(results[[column]]))))
}, numeric(1))

pad_left <- function(x, width) {
  paste0(strrep(" ", width - display_width(x)), x)
}

pad_right <- function(x, width) {
  paste0(x, strrep(" ", width - display_width(x)))
}

pad_center <- function(x, width) {
  left <- floor((width - display_width(x)) / 2)
  right <- width - display_width(x) - left
  paste0(strrep(" ", left), x, strrep(" ", right))
}

format_row <- function(values, header = FALSE) {
  values <- as.character(values)
  pieces <- character(length(values))
  for (i in seq_along(values)) {
    column <- columns[i]
    if (header && column %in% c("n", "r", "95%-CI")) {
      pieces[i] <- pad_center(values[i], widths[i])
    } else if (header || column == "station") {
      pieces[i] <- pad_right(values[i], widths[i])
    } else {
      pieces[i] <- pad_left(values[i], widths[i])
    }
  }
  paste(pieces, collapse = "   ")
}

cat(format_row(columns, header = TRUE), "\n", sep = "")
cat(strrep("–", sum(widths) + 3 * (length(widths) - 1)), "\n", sep = "")
for (i in seq_len(nrow(results))) {
  cat(format_row(results[i, columns]), "\n", sep = "")
}
