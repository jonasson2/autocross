grimsey_ann <- read.table("grimsey_annual.txt", header = TRUE)

plot_grimsey_ann <- function() {
  yrange <- range(grimsey_ann$SST, grimsey_ann$T, na.rm = TRUE)

  plot(
    grimsey_ann$year, grimsey_ann$SST,
    type = "l",
    col = "steelblue",
    lwd = 2,
    ylim = yrange,
    xlab = "Year",
    ylab = "Temperature (deg C)",
    main = "Grimsey annual mean temperatures"
  )

  lines(grimsey_ann$year, grimsey_ann$T, col = "firebrick", lwd = 2)

  legend(
    "topleft",
    legend = c("SST", "T"),
    col = c("steelblue", "firebrick"),
    lwd = 2,
    bty = "n"
  )
}

png("grimsey_annual_plot.png", width = 1200, height = 800, res = 150)
plot_grimsey_ann()
dev.off()

if (interactive()) {
  plot_grimsey_ann()
} else if (Sys.info()[["sysname"]] == "Darwin") {
  system2("open", "grimsey_annual_plot.png", wait = FALSE)
}
