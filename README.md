# ARIM-analysis
Analysis of oil and gas development on sagebrush in the Atlantic Rim project area, BLM, Wyoming

# Data compilation
00-Data_wrangle.R - Compiles data for analysis to produce Data_compiled.RDATA (requires access to Bird Conservancy internal database; not necessary to run if this file has already been compiled).
01-Explore_covariates.R - Explores and tabulates summaries of covariates. Also includes a summary of elevations reported in the manuscript.
01-Tabulate_sample_sizes.R - Generates table of sample sizes (Table 1 in manuscript).

# Analysis
02-Analyze.R - Implements model fit using the R package nimble (not necessary to run if the R object file "mod_path" has already been generated).

# Results
03-Summarise_development_community_effects.R - Tabulates summaries of development effects on species richness trends.
03-Summarise_development_spp_effects.R - Tabulates summaries of development effects on species occupancy trends.
03-Tabulate_parameters.R - Tabulates summaries of community occupancy and path model parameter estimates.
04-Plot_community_mechanisms.R - Generates plots of species richness in relation to mechanistic covariates (Figure 6).
04-Plot_community_trends.R - Generates plots of species richness trends (Figure 4).
04-Plot_spp_effects.R - Generates plots displaying regression coefficients for species covariates relationships (Figure 2, Appendices S5, S7).
04-Plot_spp_mechanisms.R - Generates plots displaying occupancy of select species in relation to mechanistic covariates (not included in manuscript).
04-Plot_spp_trends.R - Generates plots displaying species occupancy trends (Figure 3 and additional plot with more species not included in manuscript).
04-Plot_spp_trends_odds.R - Generates alternate version of plots displaying species trends in odds occupancy (Appendix S8).
04-Summarise_community_mechanisms.R - Tabulates summaries of evidence for mechanisms underlying development effects on richness trends (lower Table 3).
04-Summarise_spp_mechanisms.R - Tabulates summaries of evidence for mechanisms underlying development effects on species occupancy trends (upper Table 3).
04-Summarise_spp_mechanisms_odds.R - Tabulates summaries of evidence for mechanisms underlying development effects on species odds occupancy trends (Appendix S8).

# Sourced scripts
Calculate_mech_path_covariate_values.R - Used to derive stratum-specific covariate values for estimating mechanistic pathway contributions.
Data_processing.R - Final data processing implmented before analysis.
model_path.nimble - Model code defining community occupancy model and regression models for path analysis implemented in the R package nimble.