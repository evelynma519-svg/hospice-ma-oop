# =============================================================================
# run_all.R  --  reproduce everything from the raw HRS files.
# Usage (from the repository root):  Rscript run_all.R
# =============================================================================
source(file.path("R", "00_setup.R"))
source(file.path("R", "01_build.R"))      # -> data/derived/analysis*.rds, sample-size diagnostics
source(file.path("R", "02_replicate.R"))  # -> original specifications in R
source(file.path("R", "03_improved.R"))   # -> own controls, common sample, two-part model
message("Done. Tables in output/tables/, figures in output/figures/.")
