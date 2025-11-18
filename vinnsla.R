setwd("~/autocross")
library(ggplot2)
library(autocross)
library(pracma)

# All dataframes go backwards time and have a column yr

# READ TEMPERATUE DATA
debugSource('read_data.R')
t.res = get_temperature()
tasi = t.res$tasi
tasim = t.res$tasim
sty = t.res$sty
græn = t.res$græn

# READ ICE CORE DATA
dye3.res = get_dye3()
dye3.monthly = dye3.res$monthly
dye3.annual = dye3.res$annual
renland = get_renland()

# READ MARINE CORE DATA
sst = get_sst()

# READ HAK-HVT TEMPERATURE PROXY
lakes = get_lakes()
#lakes$tproxy = detrend(lakes$tproxy)

# PLOT MONTHLY AND ANNUAL MEANS
source('plot_d18_t.R')  # To do: búa til tvö plott.
plot_monthly(tasim, dye3.monthly)
plot_annual(dye3.annual, tasi, græn, sty)

pretty.print <- function(info, beta, M) {
  # info : a list of length n, where each info[[i]] is a list with
  #        elements $avg_int, $min_int, $max_int (all numeric).
  # beta : numeric vector of length p
  # M    : numeric matrix with n rows and p columns
  
  # 1) Extract and format the three interval‐fields from info[[i]]:
  avg_vals <- sapply(info, `[[`, "avg_int")
  min_vals <- sapply(info, `[[`, "min_int")
  max_vals <- sapply(info, `[[`, "max_int")
  
  avg_labels <- formatC(avg_vals, format = "f", digits = 3)
  min_labels <- formatC(min_vals, format = "f", digits = 3)
  max_labels <- formatC(max_vals, format = "f", digits = 3)
  
  # 2) Format the matrix M exactly as before (fixed, 3 decimals) and reshape:
  F_chr <- formatC(M, format = "f", digits = 3)
  F_mat <- matrix(F_chr, nrow = nrow(M), ncol = ncol(M))
  
  # 3) Format the beta labels (header for the numeric columns):
  beta_labels <- formatC(beta, format = "f", digits = 3)
  
  # 4) Figure out a uniform column width:
  #    We need to consider the longest of:
  #      - the literal column‐names "avg_int", "min_int", "max_int"
  #      - the β‐labels
  #      - the row‐labels (avg_labels, min_labels, max_labels)
  #      - every entry of F_mat
  all_text <- c(
    "avg_int", "min_int", "max_int",
    beta_labels,
    avg_labels, min_labels, max_labels,
    as.vector(F_mat)
  )
  col_width <- max(nchar(all_text)) + 2
  
  pad <- function(x) sprintf(paste0("%-", col_width, "s"), x)
  
  # 5) Print the header row:
  cat(
    pad("avg_int"),
    pad("min_int"),
    pad("max_int"),
    paste(sapply(beta_labels, pad), collapse = ""),
    "\n"
  )
  
  # 6) Print a separator line of dashes:
  total_cols <- 3 + length(beta_labels)
  cat(rep("-", col_width * total_cols), sep = "", "\n")
  
  # 7) Print each data row:
  for (i in seq_along(info)) {
    cat(
      pad(avg_labels[i]),
      pad(min_labels[i]),
      pad(max_labels[i]),
      paste(sapply(F_mat[i, ], pad), collapse = ""),
      "\n"
    )
  }
}

pretty.print.1 <- function(interval, beta, M) {
  F <- formatC(M, format = "f", digits = 3)
  F <- matrix(F, nrow = nrow(M), ncol = ncol(M))
  beta_labels <- formatC(beta, format = "f", digits = 3)
  interval_labels <- formatC(interval, format = "f", digits = 3)
  col_width <- max(nchar(c("interval", beta_labels, interval_labels, F))) + 2
  pad <- function(x) sprintf(paste0("%-", col_width, "s"), x)
  cat(pad("interval"), paste(sapply(beta_labels, pad), collapse = ""), "\n")
  total_cols <- 1 + length(beta_labels)
  cat(rep("-", col_width * total_cols), sep = "", "\n")
  for (i in seq_along(interval_labels)) {
    cat(pad(interval_labels[i]), paste(sapply(F[i, ], pad), collapse = ""), "\n")
  }
}

#s = find_CI(sty, 't', dye3.annual, 'd18', beta=1, n=10)
#s = find_CI(sty, 't', dye3.annual, 'd18', beta=1, n=10)

compute.table <- function(series1, var1, series2, var2, n.values) {
  cat("Computing tables for correlations")
  beta.values = seq(1, 2, by=0.2)
  M = length(n.values)
  N = length(beta.values)
  R = matrix(0, M, N)
  ci.width = matrix(0, M, N)
  info = vector('list', M)
  for (i in 1:M) {
    n = n.values[i]
    for (j in 1:N) {
      beta = beta.values[j]
      CI.result = find_CI(series1, var1, series2, var2, beta, n)
      sxy = CI.result$CI
      R[i,j] = sxy$r
      ci.width[i,j] = sxy$ci[2] - sxy$ci[1]
    }
    info[[i]] = CI.result$info
  }
  signif.indicator = ci.width/2/R
  cat('\nCORRELATIONS:\n')
  pretty.print(info, beta.values, R)
  cat('\nCONFIDENCE INTERVAL WIDHTS\n')
  pretty.print(info, beta.values, ci.width)
  cat('\nSIGNIFICANCE INDICATORS\n')
  pretty.print(info, beta.values, signif.indicator)
}

df = simulate_shifted(500, 0.4, 5)
n.values = c(20, 50, 100, 200, 400, 800)
compute.table(df, 'x', df, 'y', n.values)
# compute.table(dye3.annual, 'd18', lakes, 'tproxy', n.values)
# compute.table(lakes, 'tproxy', sst, 'SST')
# CIxy = find_CI(dye3.annual, 'd18', sst, 'SST', 0.8, 100)
# CIxy = find_CI(renland, 'd18', dye3.annual, 'd18', 0.9, 10)
# CI = CIxy$CI
# xy = CIxy$xy
# str(CI)
