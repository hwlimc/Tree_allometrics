BEF Bayesian distribution-comparison outputs

Generated: 2026-07-13 19:58:04.597606
Run ID: 20260713195804
Project root: /Users/hyli0001/wrd/b/Dynamic_allometrics
Output directory: processed_data/bef_bayes_dist_comparision
Discovered model structures: 13
Discovered families: gamma, lognormal, student

Step 1: MCMC diagnostics
- Use 01_mcmc_diagnostics.txt to screen divergences, max Rhat, and ESS.

Step 2: posterior predictive checks
- Use 02_ppcheck_observed_yrep_summary.txt to compare yrep behavior across families within a model structure.
- Values are back-transformed to the original BEF response scale.

Step 3: LOO distribution comparison
- Use 03_loo_distribution_comparison_by_structure.txt as the main ranking table.
- rank_within_structure = 1 is the best distribution assumption for that model structure by summed elpd_loo.
- elpd_diff_from_structure_best and looic_diff_from_structure_best are computed within each model structure.
- Student/log-response fits are Jacobian-adjusted before LOO comparison when their data use z1/z2.
- pareto_warning marks models with one or more Pareto k values above 0.7.
