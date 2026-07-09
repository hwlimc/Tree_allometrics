BEF Bayesian gamma-only post-processing outputs

Generated: 2026-07-09 21:55:29.455641
Run ID: 20260709
Project root: /Users/hyli0001/wrd/b/Dynamic_allometrics
Output directory: processed_data/bef_bayes_gamma/20260709
Discovered gamma models: base_k0_gamma, ftp_k0_gamma, ftp_k1_gamma, ftp_sp_k0_gamma, ftp_sp_k1_gamma, ftp_sp_k2_gamma, PFT_k0_gamma, PFT_k1_gamma, PFT_sp_k0_gamma, PFT_sp_k1_gamma, PFT_sp_k2_gamma, sp_k0_gamma, sp_k1_gamma
Best gamma model: ftp_sp_k0_gamma

Scale notes:
- Only gamma model files are loaded by this workflow.
- Figure datasets are on the original response scale: befa.st and befr.st.
- LOO compares hierarchical structures within the gamma distribution.
- Publication curve files use posterior expected responses, not full posterior predictive yrep draws.
- 06_publication_response_curve_draws_gamma.txt stores a compact draw sample for plotting, not all posterior draws.
- 05_response_vs_rsd_h1_h2.pdf uses conditional h1/h2 predictions where those grouping levels exist.
- Grouped h1/h2 prediction lines are clipped to the observed rsd range for each h1/h2 combination.

Recommended publication inputs:
- Main fitted-curve figure: 06_publication_response_curve_summary_gamma.txt.
- Optional spaghetti curves: 06_publication_response_curve_draws_gamma.txt.
- Grouped h1/h2 diagnostic figure: 05_response_vs_rsd_h1_h2.pdf.
- Gamma hierarchy ranking: 05_loo_total_ranking.txt.
- Diagnostics and screening: 01_mcmc_diagnostics.txt and 03_ppcheck_observed_yrep_summary.txt.
- Column and file descriptions: 00_output_dictionary.txt.
