BEF Bayesian post-processing outputs

Generated: 2026-07-09 19:56:28.052477
Run ID: 20260709
Project root: /Users/hyli0001/wrd/b/Dynamic_allometrics
Output directory: processed_data/bef_bayes/20260709
Discovered model families: gamma, gaussian, lognormal, student
Model comparison family: gamma
Best comparison model: ftp_sp_k0_gamma

Scale notes:
- Figure datasets are on the original response scale: befa.st and befr.st.
- Student models saved on log responses are Jacobian-adjusted for LOO, but excluded from model comparison.
- Publication curve files use posterior expected responses, not full posterior predictive yrep draws.
- 06_publication_response_curve_draws_gamma.txt stores a compact draw sample for plotting, not all posterior draws.

Recommended publication inputs:
- Main fitted-curve figure: 06_publication_response_curve_summary_gamma.txt.
- Optional spaghetti curves: 06_publication_response_curve_draws_gamma.txt.
- Gamma model ranking: 05_loo_total_ranking.txt.
- Diagnostics and screening: 01_mcmc_diagnostics.txt and 03_ppcheck_observed_yrep_summary.txt.
- Column and file descriptions: 00_output_dictionary.txt.
