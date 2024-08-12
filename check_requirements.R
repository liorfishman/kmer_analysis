use <- function(package, version=0, ...) {
  package <- as.character(substitute(package))
  if(!require(package, ..., warn.conflicts = FALSE, character.only=TRUE)) {
    return(invisible(FALSE))
  }
  pver <- packageVersion(package)
  if (compareVersion(as.character(pver), as.character(version)) < 0) {
    cat("Version ", version, " of '", package,
         "' required, but only ", as.character(pver), " is available",
        file = stderr(), sep = "")
    return(invisible(FALSE))
  }
  return(invisible(TRUE))
}

use(dplyr, 1.0)
use(tidyr, 1.3)
use(seqinr, 4.2)
use(stringr, 1.5)
use(data.table, 1.14)
use(tidyverse, 1.3)
use(kebabs, 1.28)

.Internal(printDeferredWarnings())

if (!use(dplyr, 1.0)) install.packages('dplyr')
if (!use(tidyr, 1.3)) install.packages('tidyr')
if (!use(seqinr, 4.2)) install.packages('seqinr')
if (!use(stringr, 1.5)) install.packages('stringr')
if (!use(data.table, 1.14)) install.packages('data.table')
if (!use(tidyverse, 1.3)) install.packages('tidyverse')
if (!use(kebabs, 1.28)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("kebabs")
}
