library('usethis')
create_packages('autocross')
usethis::use_mit_license("Kristín Björg Ólafsdóttir")
usethis::use_description(fields = list(
  Title = "Autocross",
  Description = paste(
    "Compute calibrated confidence interval for two auto and cross correlated",
    "time series with variable time step"),
  `Authors@R` = c(
    person("Kristín Björg", "Ólafsdóttir", email = "kbo@vedur.is", role = c("aut", "cre")),
    person("Kristján", "Jónasson", email = "jonasson@hi.is", role = "aut"),
    person("Manfred", "Mudelsee", email = "mudelsee@climate-risk-analysis.com", role = "aut")
  ),
  version='1.0.0',
))
