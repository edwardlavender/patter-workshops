# Extract the R code from selected .qmd files
# * This is useful for debugging
input <- "./doc/patter-v.1.0-1/patter-input.qmd"
output <- "./doc/patter-v.1.0-1/patter-input.R"
knitr::purl(input, output = output)
rstudioapi::navigateToFile(output)
