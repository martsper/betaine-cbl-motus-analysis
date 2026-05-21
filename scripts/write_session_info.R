################################################################################
## Write R session information
################################################################################

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(tidyr)
})

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
