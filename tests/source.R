cat("\n\nSource.R found\n\n")

library(MarkedPointProcess,
        lib=if (file.exists("/home/schlather/TMP/MarkedPointProcess"))
        "~/TMP")

.path <- "~/R/MPP/MarkedPointProcess/R/"
if (EXTENDED.TESTING <- file.exists(paste(.path, "mpp.R", sep=""))) {
#

  EXTENDED.TESTING <- FALSE
  Source <- function(x) {
    cat(x, "...\n")
    x <- paste(.path, x, sep="")
    source(x)
  }
#  source("/home/schlather/R/RF/RandomFields/R/rf.R")
  Source("mpp.R")
  Source("srd.R")
} else q()
