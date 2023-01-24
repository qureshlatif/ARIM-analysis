# ARIM-analysis
Analysis of oil and gas development on sagebrush in the Atlantic Rim project area, BLM, Wyoming

# Data compilation
00-Data_wrangle.R - Compiles data for analysis to produce Data_compiled.RDATA (requires access to Bird Conservancy internal database; not necessary to run if this file has already been compiled).
01-Explore_covariates.R - Explores and tabulates summaries of covariates. Also includes a summary of elevations reported in the manuscript.
01-Tabulate_sample_sizes.R - Generates table of sample sizes (Table 1 in manuscript).

# Analysis
01-Analyze.R - Implements model fit using the R package nimble (not necessary to run if the R object file "mod_path" has already been generated).

# Results

# Sourced scripts
