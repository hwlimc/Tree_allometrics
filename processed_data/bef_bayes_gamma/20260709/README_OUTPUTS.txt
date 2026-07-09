BEF Bayesian gamma-only post-processing outputs

Generated: 2026-07-09 10:46:50.601099
Run ID: 20260709
Project root: /home/hlim/wrd/b/Dynamic_allometrics
Output directory: processed_data/bef_bayes_gamma/20260709
Discovered gamma models: base_k0_gamma, PFT_k0_gamma, PFT_k1_gamma, PFT_sp_k0_gamma, PFT_sp_k1_gamma, PFT_sp_k2_gamma
Best gamma model: PFT_sp_k1_gamma

Scale notes:
- Only gamma model files are loaded by this workflow.
- Figure datasets are on the original response scale: befa.st and befr.st.
- LOO compares hierarchical structures within the gamma distribution.
- Publication curve files use posterior expected responses, not full posterior predictive yrep draws.
- 06_publication_response_curve_draws_gamma.txt stores a compact draw sample for plotting, not all posterior draws.

Recommended publication inputs:
- Main fitted-curve figure: 06_publication_response_curve_summary_gamma.txt.
- Optional spaghetti curves: 06_publication_response_curve_draws_gamma.txt.
- Gamma hierarchy ranking: 05_loo_total_ranking.txt.
- Diagnostics and screening: 01_mcmc_diagnostics.txt and 03_ppcheck_observed_yrep_summary.txt.
- Column and file descriptions: 00_output_dictionary.txt.
