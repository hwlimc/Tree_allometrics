BEF Bayesian Workflow
================

- [Overview](#overview)
- [Inputs](#inputs)
- [Diagnostics](#diagnostics)
- [Posterior Summaries](#posterior-summaries)
- [Posterior Predictive Checks](#posterior-predictive-checks)
- [Model Comparison](#model-comparison)
- [Outputs](#outputs)
- [Reproducibility](#reproducibility)

# Overview

This report is generated from the workflow script
[02_developing_BEF_Allometrics_Bayes.R](02_developing_BEF_Allometrics_Bayes.R).  
Edit the `.R` file only. The Markdown file is a GitHub-readable record
of that run.

``` r
cat("**Run ID:**", wf$run_id, "\n\n")
```

    ## **Run ID:** 20260708_220743

``` r
cat("**Run time:**", wf$run_generated_at, "\n\n")
```

    ## **Run time:** 2026-07-08 22:07:43.096461

``` r
cat("**Project root:**", project_root, "\n\n")
```

    ## **Project root:** /Users/hyli0001/wrd/b/Dynamic_allometrics

``` r
cat("**Output directory:**", wf$out_dir, "\n\n")
```

    ## **Output directory:** bayes_outputs/bef_bayes/20260708_220743

# Inputs

``` r
input_manifest <- wf$manifest
input_manifest$model_time <- format(as.POSIXct(input_manifest$model_time), "%Y-%m-%d %H:%M:%S")

knitr::kable(
  input_manifest[, c("model", "model_structure", "family", "hierarchy", "k_depth", "model_file", "model_time", "n_rows")],
  col.names = c("Model", "Structure", "Family", "Hierarchy", "k depth", "File", "File time", "Rows"),
  align = c("l", "l", "l", "l", "l", "l", "l", "r")
)
```

|     | Model               | Structure | Family    | Hierarchy   | k depth | File                                                         | File time           | Rows |
|:----|:--------------------|:----------|:----------|:------------|:--------|:-------------------------------------------------------------|:--------------------|-----:|
| 16  | base_k0_gamma       | base_k0   | gamma     | base        | k0      | bayes_outputs/xp_none_k0_gamma_rsd_4-4k-99-15.rds            | 2026-06-18 11:43:18 | 1132 |
| 1   | ftp_k0_gamma        | ftp_k0    | gamma     | ftp         | k0      | bayes_outputs/xp_ftp_k0_gamma_rsd_4-4k-99-15.rds             | 2026-06-18 13:38:13 | 1132 |
| 2   | ftp_k0_lognormal    | ftp_k0    | lognormal | ftp         | k0      | bayes_outputs/xp_ftp_k0_lognormal_rsd_4-4k-99-15.rds         | 2026-06-18 01:15:47 | 1132 |
| 3   | ftp_k0_student      | ftp_k0    | student   | ftp         | k0      | bayes_outputs/xp_ftp_k0_stud_rsd_4-4k-99-15.rds              | 2026-06-17 21:01:48 | 1132 |
| 4   | ftp_k1_gamma        | ftp_k1    | gamma     | ftp         | k1      | bayes_outputs/xp_ftp_k1_gamma_rsd_4-4k-99-15.rds             | 2026-06-18 14:34:29 | 1132 |
| 5   | ftp_k1_lognormal    | ftp_k1    | lognormal | ftp         | k1      | bayes_outputs/xp_ftp_k1_lognormal_rsd_4-4k-99-15.rds         | 2026-06-18 01:19:18 | 1132 |
| 6   | ftp_k1_student      | ftp_k1    | student   | ftp         | k1      | bayes_outputs/xp_ftp_k1_stud_rsd_4-4k-99-15.rds              | 2026-06-17 21:43:21 | 1132 |
| 7   | ftp_sp_k0_gamma     | ftp_sp_k0 | gamma     | ftp-sp_code | k0      | bayes_outputs/xp_ftp-sp_code_k0_gamma_rsd_4-4k-99-15.rds     | 2026-06-18 14:23:41 | 1132 |
| 8   | ftp_sp_k0_lognormal | ftp_sp_k0 | lognormal | ftp-sp_code | k0      | bayes_outputs/xp_ftp-sp_code_k0_lognormal_rsd_4-4k-99-15.rds | 2026-06-18 01:18:52 | 1132 |
| 9   | ftp_sp_k0_student   | ftp_sp_k0 | student   | ftp-sp_code | k0      | bayes_outputs/xp_ftp-sp_code_k0_stud_rsd_4-4k-99-15.rds      | 2026-06-17 21:22:56 | 1132 |
| 10  | ftp_sp_k1_gamma     | ftp_sp_k1 | gamma     | ftp-sp_code | k1      | bayes_outputs/xp_ftp-sp_code_k1_gamma_rsd_4-4k-99-15.rds     | 2026-06-18 14:31:54 | 1132 |
| 11  | ftp_sp_k1_lognormal | ftp_sp_k1 | lognormal | ftp-sp_code | k1      | bayes_outputs/xp_ftp-sp_code_k1_lognormal_rsd_4-4k-99-15.rds | 2026-06-18 01:22:29 | 1132 |
| 12  | ftp_sp_k1_student   | ftp_sp_k1 | student   | ftp-sp_code | k1      | bayes_outputs/xp_ftp-sp_code_k1_stud_rsd_4-4k-99-15.rds      | 2026-06-17 22:17:56 | 1132 |
| 13  | ftp_sp_k2_gamma     | ftp_sp_k2 | gamma     | ftp-sp_code | k2      | bayes_outputs/xp_ftp-sp_code_k2_gamma_rsd_4-4k-99-15.rds     | 2026-06-18 14:52:57 | 1132 |
| 14  | ftp_sp_k2_lognormal | ftp_sp_k2 | lognormal | ftp-sp_code | k2      | bayes_outputs/xp_ftp-sp_code_k2_lognormal_rsd_4-4k-99-15.rds | 2026-06-18 01:23:29 | 1132 |
| 15  | ftp_sp_k2_student   | ftp_sp_k2 | student   | ftp-sp_code | k2      | bayes_outputs/xp_ftp-sp_code_k2_stud_rsd_4-4k-99-15.rds      | 2026-06-17 22:19:08 | 1132 |
| 17  | PFT_k0_student      | PFT_k0    | student   | PFT         | k0      | bayes_outputs/xp_PFT_k0_stud_rsd_4-4k-99-15.rds              | 2026-06-17 21:14:35 | 1132 |
| 18  | PFT_k1_student      | PFT_k1    | student   | PFT         | k1      | bayes_outputs/xp_PFT_k1_stud_rsd_4-4k-99-15.rds              | 2026-06-17 22:18:46 | 1132 |
| 19  | PFT_sp_k0_student   | PFT_sp_k0 | student   | PFT-sp_code | k0      | bayes_outputs/xp_PFT-sp_code_k0_stud_rsd_4-4k-99-15.rds      | 2026-06-17 21:49:30 | 1132 |
| 20  | PFT_sp_k1_student   | PFT_sp_k1 | student   | PFT-sp_code | k1      | bayes_outputs/xp_PFT-sp_code_k1_stud_rsd_4-4k-99-15.rds      | 2026-06-17 22:14:17 | 1132 |
| 21  | PFT_sp_k2_student   | PFT_sp_k2 | student   | PFT-sp_code | k2      | bayes_outputs/xp_PFT-sp_code_k2_stud_rsd_4-4k-99-15.rds      | 2026-06-17 22:19:19 | 1132 |
| 22  | sp_k0_gamma         | sp_k0     | gamma     | sp_code     | k0      | bayes_outputs/xp_sp_code_k0_gamma_rsd_4-4k-99-15.rds         | 2026-06-18 12:04:21 | 1132 |
| 23  | sp_k0_lognormal     | sp_k0     | lognormal | sp_code     | k0      | bayes_outputs/xp_sp_code_k0_lognormal_rsd_4-4k-99-15.rds     | 2026-06-18 01:21:21 | 1132 |
| 24  | sp_k0_student       | sp_k0     | student   | sp_code     | k0      | bayes_outputs/xp_sp_code_k0_stud_rsd_4-4k-99-15.rds          | 2026-06-17 20:56:20 | 1132 |
| 25  | sp_k1_gamma         | sp_k1     | gamma     | sp_code     | k1      | bayes_outputs/xp_sp_code_k1_gamma_rsd_4-4k-99-15.rds         | 2026-06-18 12:56:34 | 1132 |
| 26  | sp_k1_lognormal     | sp_k1     | lognormal | sp_code     | k1      | bayes_outputs/xp_sp_code_k1_lognormal_rsd_4-4k-99-15.rds     | 2026-06-18 01:20:14 | 1132 |
| 27  | sp_k1_student       | sp_k1     | student   | sp_code     | k1      | bayes_outputs/xp_sp_code_k1_stud_rsd_4-4k-99-15.rds          | 2026-06-17 21:05:44 | 1132 |

# Diagnostics

``` r
diag_tbl <- wf$step1[, c(
  "model", "model_structure", "family", "divergences", "max_treedepth", "max_rhat",
  "max_rhat_parameter", "min_bulk_ess", "min_bulk_ess_parameter",
  "min_tail_ess", "min_tail_ess_parameter"
)]

knitr::kable(
  diag_tbl,
  digits = 3,
  align = c("l", "l", "l", "r", "r", "r", "l", "r", "l", "r", "l")
)
```

|                     | model               | model_structure | family    | divergences | max_treedepth | max_rhat | max_rhat_parameter                     | min_bulk_ess | min_bulk_ess_parameter                | min_tail_ess | min_tail_ess_parameter                  |
|:--------------------|:--------------------|:----------------|:----------|------------:|--------------:|---------:|:---------------------------------------|-------------:|:--------------------------------------|-------------:|:----------------------------------------|
| base_k0_gamma       | base_k0_gamma       | base_k0         | gamma     |           0 |             9 |    1.002 | lp\_\_                                 |     1583.585 | lp\_\_                                |     3083.463 | lp\_\_                                  |
| ftp_k0_gamma        | ftp_k0_gamma        | ftp_k0          | gamma     |           0 |            10 |    1.002 | lp\_\_                                 |     1587.842 | lp\_\_                                |     2853.027 | sd_h1\_\_y2_logL_Intercept              |
| ftp_k0_lognormal    | ftp_k0_lognormal    | ftp_k0          | lognormal |           0 |             7 |    1.003 | r_h1\_\_y1m1_logA\[mono_B,Intercept\]  |     2046.136 | r_h1\_\_y1m1_logA\[mono_N,Intercept\] |     1899.427 | sd_h1\_\_y2_logA_Intercept              |
| ftp_k0_student      | ftp_k0_student      | ftp_k0          | student   |           0 |             7 |    1.002 | b_z2_logL_Intercept                    |     1819.662 | lp\_\_                                |     1067.142 | r_h1\_\_z2_logA\[mono_B,Intercept\]     |
| ftp_k1_gamma        | ftp_k1_gamma        | ftp_k1          | gamma     |           0 |            11 |    1.002 | sd_h1\_\_y2_logk_Intercept             |     1857.799 | lp\_\_                                |     2897.967 | r_h1\_\_y1m1_logL\[mono_B,Intercept\]   |
| ftp_k1_lognormal    | ftp_k1_lognormal    | ftp_k1          | lognormal |          17 |             8 |    1.004 | r_h1\_\_y1m1_logA\[mono_N,Intercept\]  |     1941.052 | lp\_\_                                |     2351.604 | sd_h1\_\_y1m1_logk_Intercept            |
| ftp_k1_student      | ftp_k1_student      | ftp_k1          | student   |          28 |             8 |    1.001 | z_4\[1,1\]                             |     1868.054 | lp\_\_                                |     2169.409 | sd_h1\_\_z1_logk_Intercept              |
| ftp_sp_k0_gamma     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     |           0 |            12 |    1.005 | z_6\[1,9\]                             |      515.005 | b_y2_logk_Intercept                   |      783.671 | r_h2\_\_y2_logL\[mono_B.PxT,Intercept\] |
| ftp_sp_k0_lognormal | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal |           0 |             7 |    1.003 | sd_h1\_\_y2_logA_Intercept             |     1616.462 | sd_h1\_\_y1m1_logA_Intercept          |     1042.721 | r_h1\_\_y2_logA\[mono_N,Intercept\]     |
| ftp_sp_k0_student   | ftp_sp_k0_student   | ftp_sp_k0       | student   |           0 |             7 |    1.004 | b_z1_logk_Intercept                    |     1761.074 | lp\_\_                                |     1347.219 | r_h1\_\_z2_logA\[mono_B,Intercept\]     |
| ftp_sp_k1_gamma     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     |           0 |            10 |    1.008 | sd_h2\_\_y2_logA_Intercept             |      456.353 | sd_h2\_\_y2_logA_Intercept            |      767.130 | r_h2\_\_y2_logL\[mono_B.PxT,Intercept\] |
| ftp_sp_k1_lognormal | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal |          11 |             9 |    1.004 | lp\_\_                                 |     1766.882 | lp\_\_                                |     2101.443 | sd_h1\_\_y1m1_logk_Intercept            |
| ftp_sp_k1_student   | ftp_sp_k1_student   | ftp_sp_k1       | student   |          11 |             9 |    1.004 | lp\_\_                                 |     1805.797 | lp\_\_                                |     1863.482 | sd_h1\_\_z1_logk_Intercept              |
| ftp_sp_k2_gamma     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     |           0 |            11 |    1.013 | r_h2\_\_y2_logk\[mono_N.CO,Intercept\] |      391.653 | r_h1\_\_y2_logk\[mono_B,Intercept\]   |      294.480 | r_h1\_\_y2_logA\[mono_B,Intercept\]     |
| ftp_sp_k2_lognormal | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal |           7 |             9 |    1.006 | r_h1\_\_y1m1_logA\[mono_B,Intercept\]  |     1864.118 | lp\_\_                                |     2501.617 | sd_h1\_\_y1m1_logk_Intercept            |
| ftp_sp_k2_student   | ftp_sp_k2_student   | ftp_sp_k2       | student   |           9 |             9 |    1.003 | r_h1\_\_z2_logk\[mono_B,Intercept\]    |     1705.920 | lp\_\_                                |     2549.367 | r_h1\_\_z1_logk\[mono_N,Intercept\]     |
| PFT_k0_student      | PFT_k0_student      | PFT_k0          | student   |           0 |             7 |    1.003 | nu_z1                                  |     1822.275 | sd_h1\_\_z2_logA_Intercept            |     1745.575 | r_h1\_\_z2_logL\[ENF,Intercept\]        |
| PFT_k1_student      | PFT_k1_student      | PFT_k1          | student   |          32 |            10 |    1.004 | lp\_\_                                 |     1497.683 | lp\_\_                                |     2450.452 | sd_h1\_\_z2_logL_Intercept              |
| PFT_sp_k0_student   | PFT_sp_k0_student   | PFT_sp_k0       | student   |           0 |             8 |    1.003 | z_6\[1,20\]                            |     1914.167 | lp\_\_                                |     2852.686 | r_h1\_\_z2_logL\[ENF,Intercept\]        |
| PFT_sp_k1_student   | PFT_sp_k1_student   | PFT_sp_k1       | student   |          28 |             8 |    1.003 | z_2\[1,4\]                             |     1448.131 | lp\_\_                                |     1924.123 | sd_h1\_\_z2_logL_Intercept              |
| PFT_sp_k2_student   | PFT_sp_k2_student   | PFT_sp_k2       | student   |          49 |             9 |    1.003 | z_4\[1,12\]                            |     1705.237 | lp\_\_                                |     3072.493 | sd_h1\_\_z1_logk_Intercept              |
| sp_k0_gamma         | sp_k0_gamma         | sp_k0           | gamma     |           0 |             9 |    1.006 | b_y1m1_logL_Intercept                  |      590.006 | b_y2_logk_Intercept                   |      490.818 | r_h1\_\_y2_logL\[CJa,Intercept\]        |
| sp_k0_lognormal     | sp_k0_lognormal     | sp_k0           | lognormal |           0 |             7 |    1.003 | z_4\[1,10\]                            |     1683.657 | lp\_\_                                |     3161.633 | lp\_\_                                  |
| sp_k0_student       | sp_k0_student       | sp_k0           | student   |           0 |             7 |    1.004 | lp\_\_                                 |     1584.181 | lp\_\_                                |     3007.544 | lp\_\_                                  |
| sp_k1_gamma         | sp_k1_gamma         | sp_k1           | gamma     |           0 |             9 |    1.004 | sd_h1\_\_y2_logk_Intercept             |      722.826 | b_y2_logk_Intercept                   |     1025.062 | sd_h1\_\_y2_logL_Intercept              |
| sp_k1_lognormal     | sp_k1_lognormal     | sp_k1           | lognormal |           2 |             8 |    1.003 | r_h1\_\_y2_logL\[QG,Intercept\]        |     1897.562 | lp\_\_                                |     3497.029 | lp\_\_                                  |
| sp_k1_student       | sp_k1_student       | sp_k1           | student   |           0 |             7 |    1.003 | z_4\[1,9\]                             |     1542.606 | lp\_\_                                |     3172.626 | lp\_\_                                  |

# Posterior Summaries

``` r
post_tbl <- wf$step2
post_tbl <- post_tbl[grepl(
  paste(
    "^b_(y1m1|y2|z1|z2)_log(L|A|k)_Intercept$",
    "^sd_h[12]__(y1m1|y2|z1|z2)_log(L|A|k)_Intercept$",
    "^shape_(y1m1|y2)$",
    "^sigma_(y1m1|y2|z1|z2)$",
    "^nu_(z1|z2)$",
    sep = "|"
  ),
  post_tbl$variable
), ]

knitr::kable(
  post_tbl[, c("model", "model_structure", "family", "variable", "mean", "median", "sd", "q5", "q95")],
  digits = 4,
  align = c("l", "l", "l", "l", "r", "r", "r", "r", "r")
)
```

|                        | model               | model_structure | family    | variable                     |    mean |  median |      sd |      q5 |     q95 |
|:-----------------------|:--------------------|:----------------|:----------|:-----------------------------|--------:|--------:|--------:|--------:|--------:|
| base_k0_gamma.1        | base_k0_gamma       | base_k0         | gamma     | b_y1m1_logL_Intercept        | -1.0159 | -1.0152 |  0.0297 | -1.0659 | -0.9685 |
| base_k0_gamma.2        | base_k0_gamma       | base_k0         | gamma     | b_y1m1_logA_Intercept        |  0.2666 |  0.2640 |  0.1761 | -0.0195 |  0.5581 |
| base_k0_gamma.3        | base_k0_gamma       | base_k0         | gamma     | b_y1m1_logk_Intercept        |  2.1598 |  2.1622 |  0.1158 |  1.9718 |  2.3467 |
| base_k0_gamma.4        | base_k0_gamma       | base_k0         | gamma     | b_y2_logL_Intercept          | -1.7253 | -1.5407 |  0.6686 | -3.0499 | -0.9991 |
| base_k0_gamma.5        | base_k0_gamma       | base_k0         | gamma     | b_y2_logA_Intercept          | -1.0195 | -1.0011 |  0.2500 | -1.4673 | -0.6627 |
| base_k0_gamma.6        | base_k0_gamma       | base_k0         | gamma     | b_y2_logk_Intercept          |  0.0388 | -0.0803 |  0.6282 | -0.8174 |  1.1602 |
| base_k0_gamma.7        | base_k0_gamma       | base_k0         | gamma     | shape_y1m1                   |  2.3173 |  2.3160 |  0.0917 |  2.1673 |  2.4690 |
| base_k0_gamma.8        | base_k0_gamma       | base_k0         | gamma     | shape_y2                     |  3.7757 |  3.7695 |  0.2263 |  3.4106 |  4.1522 |
| ftp_k0_gamma.1         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y1m1_logL_Intercept        | -1.3748 | -1.2472 |  0.4976 | -2.3717 | -0.7892 |
| ftp_k0_gamma.2         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y1m1_logA_Intercept        | -0.3588 | -0.3048 |  0.5042 | -1.2692 |  0.3724 |
| ftp_k0_gamma.3         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y1m1_logk_Intercept        |  1.9721 |  1.9747 |  0.1326 |  1.7505 |  2.1829 |
| ftp_k0_gamma.4         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y2_logL_Intercept          | -2.0424 | -1.9614 |  0.6742 | -3.2536 | -1.0873 |
| ftp_k0_gamma.5         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y2_logA_Intercept          | -0.9042 | -0.8894 |  0.3300 | -1.4467 | -0.3848 |
| ftp_k0_gamma.6         | ftp_k0_gamma        | ftp_k0          | gamma     | b_y2_logk_Intercept          |  0.0351 | -0.0346 |  0.4386 | -0.5850 |  0.8403 |
| ftp_k0_gamma.7         | ftp_k0_gamma        | ftp_k0          | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.7005 |  0.5051 |  0.6117 |  0.1526 |  1.9135 |
| ftp_k0_gamma.8         | ftp_k0_gamma        | ftp_k0          | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.8691 |  0.7086 |  0.6093 |  0.2322 |  2.0612 |
| ftp_k0_gamma.9         | ftp_k0_gamma        | ftp_k0          | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.9390 |  0.7500 |  0.7399 |  0.1817 |  2.3312 |
| ftp_k0_gamma.10        | ftp_k0_gamma        | ftp_k0          | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.4109 |  0.2722 |  0.4511 |  0.0208 |  1.2812 |
| ftp_k0_gamma.11        | ftp_k0_gamma        | ftp_k0          | gamma     | shape_y1m1                   |  2.3877 |  2.3870 |  0.0945 |  2.2346 |  2.5457 |
| ftp_k0_gamma.12        | ftp_k0_gamma        | ftp_k0          | gamma     | shape_y2                     |  4.1636 |  4.1575 |  0.2475 |  3.7640 |  4.5748 |
| ftp_k0_lognormal.1     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y1m1_logL_Intercept        | -2.9259 | -2.9271 |  0.8389 | -4.3143 | -1.5562 |
| ftp_k0_lognormal.2     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y1m1_logA_Intercept        | -1.3740 | -1.3741 |  0.7165 | -2.5495 | -0.1977 |
| ftp_k0_lognormal.3     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y1m1_logk_Intercept        |  1.0059 |  0.9226 |  1.2251 | -0.8948 |  2.9346 |
| ftp_k0_lognormal.4     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y2_logL_Intercept          | -3.0097 | -3.0006 |  0.8729 | -4.4658 | -1.5865 |
| ftp_k0_lognormal.5     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y2_logA_Intercept          | -1.4256 | -1.4266 |  0.7367 | -2.6525 | -0.2140 |
| ftp_k0_lognormal.6     | ftp_k0_lognormal    | ftp_k0          | lognormal | b_y2_logk_Intercept          |  0.6772 |  0.6220 |  1.0813 | -1.0356 |  2.5144 |
| ftp_k0_lognormal.7     | ftp_k0_lognormal    | ftp_k0          | lognormal | sd_h1\_\_y1m1_logL_Intercept |  3.6992 |  3.4762 |  1.5591 |  1.5984 |  6.5156 |
| ftp_k0_lognormal.8     | ftp_k0_lognormal    | ftp_k0          | lognormal | sd_h1\_\_y1m1_logA_Intercept |  3.2420 |  3.1988 |  1.9751 |  0.2266 |  6.5846 |
| ftp_k0_lognormal.9     | ftp_k0_lognormal    | ftp_k0          | lognormal | sd_h1\_\_y2_logL_Intercept   |  3.2427 |  3.0461 |  1.5771 |  1.0850 |  6.1052 |
| ftp_k0_lognormal.10    | ftp_k0_lognormal    | ftp_k0          | lognormal | sd_h1\_\_y2_logA_Intercept   |  3.4045 |  3.2777 |  1.7822 |  0.5452 |  6.5262 |
| ftp_k0_lognormal.11    | ftp_k0_lognormal    | ftp_k0          | lognormal | sigma_y1m1                   |  1.3229 |  1.3226 |  0.0275 |  1.2779 |  1.3688 |
| ftp_k0_lognormal.12    | ftp_k0_lognormal    | ftp_k0          | lognormal | sigma_y2                     |  1.1019 |  1.1016 |  0.0334 |  1.0477 |  1.1573 |
| ftp_k0_student.1       | ftp_k0_student      | ftp_k0          | student   | b_z1_logL_Intercept          | -2.1791 | -2.1742 |  0.8332 | -3.5635 | -0.8204 |
| ftp_k0_student.2       | ftp_k0_student      | ftp_k0          | student   | b_z1_logA_Intercept          | -0.9637 | -0.9631 |  0.7347 | -2.1613 |  0.2499 |
| ftp_k0_student.3       | ftp_k0_student      | ftp_k0          | student   | b_z1_logk_Intercept          |  0.9698 |  0.8445 |  1.2562 | -0.9782 |  2.9921 |
| ftp_k0_student.4       | ftp_k0_student      | ftp_k0          | student   | b_z2_logL_Intercept          | -2.2282 | -2.2142 |  0.8644 | -3.6741 | -0.8396 |
| ftp_k0_student.5       | ftp_k0_student      | ftp_k0          | student   | b_z2_logA_Intercept          | -0.8213 | -0.8277 |  0.7113 | -2.0162 |  0.3465 |
| ftp_k0_student.6       | ftp_k0_student      | ftp_k0          | student   | b_z2_logk_Intercept          |  0.5591 |  0.4838 |  1.0423 | -1.0536 |  2.4262 |
| ftp_k0_student.7       | ftp_k0_student      | ftp_k0          | student   | sd_h1\_\_z1_logL_Intercept   |  4.1819 |  3.9439 |  1.6034 |  2.0429 |  7.1523 |
| ftp_k0_student.8       | ftp_k0_student      | ftp_k0          | student   | sd_h1\_\_z1_logA_Intercept   |  3.5478 |  3.5335 |  2.0187 |  0.3124 |  7.0112 |
| ftp_k0_student.9       | ftp_k0_student      | ftp_k0          | student   | sd_h1\_\_z2_logL_Intercept   |  3.8163 |  3.5781 |  1.5872 |  1.6918 |  6.7665 |
| ftp_k0_student.10      | ftp_k0_student      | ftp_k0          | student   | sd_h1\_\_z2_logA_Intercept   |  3.9045 |  3.7443 |  1.7869 |  1.1841 |  7.0422 |
| ftp_k0_student.11      | ftp_k0_student      | ftp_k0          | student   | sigma_z1                     |  1.3095 |  1.3087 |  0.0284 |  1.2633 |  1.3564 |
| ftp_k0_student.12      | ftp_k0_student      | ftp_k0          | student   | sigma_z2                     |  1.0908 |  1.0904 |  0.0348 |  1.0341 |  1.1493 |
| ftp_k0_student.13      | ftp_k0_student      | ftp_k0          | student   | nu_z1                        | 62.4445 | 59.8054 | 19.2583 | 36.3245 | 98.5353 |
| ftp_k0_student.14      | ftp_k0_student      | ftp_k0          | student   | nu_z2                        | 50.9394 | 48.4014 | 17.6646 | 27.0821 | 83.8744 |
| ftp_k1_gamma.1         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y1m1_logL_Intercept        | -1.3692 | -1.2573 |  0.4901 | -2.3444 | -0.7635 |
| ftp_k1_gamma.2         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y1m1_logA_Intercept        | -0.1858 | -0.0866 |  0.4971 | -1.1381 |  0.4527 |
| ftp_k1_gamma.3         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y1m1_logk_Intercept        |  1.3197 |  1.5090 |  0.7485 | -0.1388 |  2.1979 |
| ftp_k1_gamma.4         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y2_logL_Intercept          | -2.1177 | -2.0377 |  0.6945 | -3.3594 | -1.1244 |
| ftp_k1_gamma.5         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y2_logA_Intercept          | -0.8588 | -0.8383 |  0.3147 | -1.3842 | -0.3916 |
| ftp_k1_gamma.6         | ftp_k1_gamma        | ftp_k1          | gamma     | b_y2_logk_Intercept          | -0.0304 | -0.0734 |  0.5316 | -0.8328 |  0.8965 |
| ftp_k1_gamma.7         | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.7043 |  0.5183 |  0.5897 |  0.1649 |  1.8442 |
| ftp_k1_gamma.8         | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.6895 |  0.5231 |  0.6127 |  0.0545 |  1.8693 |
| ftp_k1_gamma.9         | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y1m1_logk_Intercept |  0.9888 |  0.7749 |  0.8585 |  0.0643 |  2.6974 |
| ftp_k1_gamma.10        | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.8229 |  0.6318 |  0.7298 |  0.0737 |  2.2480 |
| ftp_k1_gamma.11        | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.3870 |  0.2449 |  0.4341 |  0.0209 |  1.2072 |
| ftp_k1_gamma.12        | ftp_k1_gamma        | ftp_k1          | gamma     | sd_h1\_\_y2_logk_Intercept   |  0.6241 |  0.4663 |  0.5764 |  0.0411 |  1.7573 |
| ftp_k1_gamma.13        | ftp_k1_gamma        | ftp_k1          | gamma     | shape_y1m1                   |  2.3889 |  2.3876 |  0.0934 |  2.2354 |  2.5457 |
| ftp_k1_gamma.14        | ftp_k1_gamma        | ftp_k1          | gamma     | shape_y2                     |  4.1624 |  4.1564 |  0.2416 |  3.7657 |  4.5728 |
| ftp_k1_lognormal.1     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y1m1_logL_Intercept        | -2.9230 | -2.9091 |  0.8424 | -4.3419 | -1.5744 |
| ftp_k1_lognormal.2     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y1m1_logA_Intercept        | -1.1340 | -1.1230 |  0.7104 | -2.3096 |  0.0102 |
| ftp_k1_lognormal.3     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y1m1_logk_Intercept        |  0.7054 |  0.6799 |  0.8877 | -0.7308 |  2.2240 |
| ftp_k1_lognormal.4     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y2_logL_Intercept          | -3.0023 | -3.0010 |  0.8733 | -4.4592 | -1.5872 |
| ftp_k1_lognormal.5     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y2_logA_Intercept          | -1.2724 | -1.2681 |  0.7081 | -2.4390 | -0.0928 |
| ftp_k1_lognormal.6     | ftp_k1_lognormal    | ftp_k1          | lognormal | b_y2_logk_Intercept          |  0.5960 |  0.6026 |  0.8745 | -0.8329 |  2.0305 |
| ftp_k1_lognormal.7     | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y1m1_logL_Intercept |  3.7062 |  3.4602 |  1.5816 |  1.5576 |  6.6560 |
| ftp_k1_lognormal.8     | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y1m1_logA_Intercept |  1.4754 |  0.9098 |  1.5578 |  0.0544 |  4.7352 |
| ftp_k1_lognormal.9     | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y1m1_logk_Intercept |  2.6151 |  2.4397 |  1.5838 |  0.3079 |  5.4105 |
| ftp_k1_lognormal.10    | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y2_logL_Intercept   |  3.2533 |  3.0279 |  1.6115 |  1.0283 |  6.1724 |
| ftp_k1_lognormal.11    | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y2_logA_Intercept   |  1.8847 |  1.3891 |  1.6832 |  0.0916 |  5.1407 |
| ftp_k1_lognormal.12    | ftp_k1_lognormal    | ftp_k1          | lognormal | sd_h1\_\_y2_logk_Intercept   |  2.4959 |  2.3116 |  1.6991 |  0.1872 |  5.5773 |
| ftp_k1_lognormal.13    | ftp_k1_lognormal    | ftp_k1          | lognormal | sigma_y1m1                   |  1.3225 |  1.3222 |  0.0278 |  1.2774 |  1.3684 |
| ftp_k1_lognormal.14    | ftp_k1_lognormal    | ftp_k1          | lognormal | sigma_y2                     |  1.1008 |  1.1006 |  0.0335 |  1.0463 |  1.1564 |
| ftp_k1_student.1       | ftp_k1_student      | ftp_k1          | student   | b_z1_logL_Intercept          | -2.1564 | -2.1446 |  0.8231 | -3.5242 | -0.8267 |
| ftp_k1_student.2       | ftp_k1_student      | ftp_k1          | student   | b_z1_logA_Intercept          | -0.7318 | -0.7193 |  0.6799 | -1.8812 |  0.3715 |
| ftp_k1_student.3       | ftp_k1_student      | ftp_k1          | student   | b_z1_logk_Intercept          |  0.6777 |  0.6804 |  0.8741 | -0.7520 |  2.1569 |
| ftp_k1_student.4       | ftp_k1_student      | ftp_k1          | student   | b_z2_logL_Intercept          | -2.2236 | -2.2117 |  0.8466 | -3.6487 | -0.8582 |
| ftp_k1_student.5       | ftp_k1_student      | ftp_k1          | student   | b_z2_logA_Intercept          | -0.6938 | -0.6929 |  0.7061 | -1.8606 |  0.4748 |
| ftp_k1_student.6       | ftp_k1_student      | ftp_k1          | student   | b_z2_logk_Intercept          |  0.5775 |  0.5799 |  0.8766 | -0.8638 |  2.0344 |
| ftp_k1_student.7       | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z1_logL_Intercept   |  4.1844 |  3.9573 |  1.5771 |  2.0372 |  7.1016 |
| ftp_k1_student.8       | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z1_logA_Intercept   |  1.4277 |  0.8430 |  1.5517 |  0.0627 |  4.7330 |
| ftp_k1_student.9       | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z1_logk_Intercept   |  2.7960 |  2.6122 |  1.5622 |  0.4467 |  5.5843 |
| ftp_k1_student.10      | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z2_logL_Intercept   |  3.8164 |  3.5988 |  1.5736 |  1.7105 |  6.7574 |
| ftp_k1_student.11      | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z2_logA_Intercept   |  1.8864 |  1.3084 |  1.7651 |  0.0800 |  5.3674 |
| ftp_k1_student.12      | ftp_k1_student      | ftp_k1          | student   | sd_h1\_\_z2_logk_Intercept   |  2.7678 |  2.6541 |  1.7466 |  0.2184 |  5.8407 |
| ftp_k1_student.13      | ftp_k1_student      | ftp_k1          | student   | sigma_z1                     |  1.3088 |  1.3081 |  0.0283 |  1.2633 |  1.3567 |
| ftp_k1_student.14      | ftp_k1_student      | ftp_k1          | student   | sigma_z2                     |  1.0902 |  1.0893 |  0.0345 |  1.0345 |  1.1481 |
| ftp_k1_student.15      | ftp_k1_student      | ftp_k1          | student   | nu_z1                        | 62.5010 | 59.8480 | 19.4812 | 35.8109 | 98.8292 |
| ftp_k1_student.16      | ftp_k1_student      | ftp_k1          | student   | nu_z2                        | 51.1375 | 48.4399 | 17.6349 | 27.4618 | 84.7803 |
| ftp_sp_k0_gamma.1      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y1m1_logL_Intercept        | -1.3809 | -1.2695 |  0.4617 | -2.3027 | -0.8380 |
| ftp_sp_k0_gamma.2      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y1m1_logA_Intercept        | -0.1451 | -0.0589 |  0.4468 | -1.0157 |  0.4317 |
| ftp_sp_k0_gamma.3      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y1m1_logk_Intercept        |  1.9281 |  1.9302 |  0.1152 |  1.7346 |  2.1146 |
| ftp_sp_k0_gamma.4      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y2_logL_Intercept          | -1.6949 | -1.5629 |  0.6324 | -2.9276 | -0.9270 |
| ftp_sp_k0_gamma.5      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y2_logA_Intercept          | -0.9962 | -0.9988 |  0.3569 | -1.5869 | -0.4220 |
| ftp_sp_k0_gamma.6      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | b_y2_logk_Intercept          |  0.5798 |  0.5583 |  0.5825 | -0.3200 |  1.5275 |
| ftp_sp_k0_gamma.7      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.5854 |  0.4047 |  0.5809 |  0.0460 |  1.7371 |
| ftp_sp_k0_gamma.8      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h2\_\_y1m1_logL_Intercept |  0.5506 |  0.5374 |  0.0969 |  0.4156 |  0.7269 |
| ftp_sp_k0_gamma.9      | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.5749 |  0.3983 |  0.5800 |  0.0347 |  1.7230 |
| ftp_sp_k0_gamma.10     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h2\_\_y1m1_logA_Intercept |  0.3508 |  0.3337 |  0.1191 |  0.1902 |  0.5725 |
| ftp_sp_k0_gamma.11     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.6269 |  0.4409 |  0.6254 |  0.0510 |  1.8257 |
| ftp_sp_k0_gamma.12     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h2\_\_y2_logL_Intercept   |  0.5094 |  0.4493 |  0.2347 |  0.2697 |  0.9574 |
| ftp_sp_k0_gamma.13     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.4166 |  0.2708 |  0.4653 |  0.0224 |  1.2872 |
| ftp_sp_k0_gamma.14     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | sd_h2\_\_y2_logA_Intercept   |  0.4676 |  0.4183 |  0.3103 |  0.0570 |  1.0328 |
| ftp_sp_k0_gamma.15     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | shape_y1m1                   |  3.3935 |  3.3899 |  0.1403 |  3.1658 |  3.6267 |
| ftp_sp_k0_gamma.16     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | shape_y2                     |  5.5334 |  5.5296 |  0.3493 |  4.9680 |  6.1128 |
| ftp_sp_k0_lognormal.1  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y1m1_logL_Intercept        | -2.9053 | -2.8978 |  0.8319 | -4.3215 | -1.5625 |
| ftp_sp_k0_lognormal.2  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y1m1_logA_Intercept        | -1.3686 | -1.3607 |  0.7114 | -2.5493 | -0.2155 |
| ftp_sp_k0_lognormal.3  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y1m1_logk_Intercept        |  0.9413 |  0.8204 |  1.2320 | -0.9578 |  2.9316 |
| ftp_sp_k0_lognormal.4  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y2_logL_Intercept          | -2.9712 | -2.9648 |  0.8711 | -4.4422 | -1.5543 |
| ftp_sp_k0_lognormal.5  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y2_logA_Intercept          | -1.4100 | -1.4061 |  0.7214 | -2.5997 | -0.2301 |
| ftp_sp_k0_lognormal.6  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | b_y2_logk_Intercept          |  0.6186 |  0.5597 |  1.0630 | -1.0468 |  2.4518 |
| ftp_sp_k0_lognormal.7  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h1\_\_y1m1_logL_Intercept |  3.8453 |  3.6415 |  1.5746 |  1.6715 |  6.7427 |
| ftp_sp_k0_lognormal.8  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h2\_\_y1m1_logL_Intercept |  0.6302 |  0.4631 |  0.5872 |  0.0355 |  1.7926 |
| ftp_sp_k0_lognormal.9  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h1\_\_y1m1_logA_Intercept |  3.4750 |  3.4210 |  2.0403 |  0.2400 |  7.1145 |
| ftp_sp_k0_lognormal.10 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h2\_\_y1m1_logA_Intercept |  0.6170 |  0.4537 |  0.5747 |  0.0338 |  1.7532 |
| ftp_sp_k0_lognormal.11 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h1\_\_y2_logL_Intercept   |  3.4250 |  3.2221 |  1.5957 |  1.2261 |  6.3487 |
| ftp_sp_k0_lognormal.12 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h2\_\_y2_logL_Intercept   |  0.6182 |  0.4555 |  0.5795 |  0.0335 |  1.7529 |
| ftp_sp_k0_lognormal.13 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h1\_\_y2_logA_Intercept   |  3.6466 |  3.4864 |  1.7636 |  0.8958 |  6.7996 |
| ftp_sp_k0_lognormal.14 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sd_h2\_\_y2_logA_Intercept   |  0.6208 |  0.4426 |  0.5796 |  0.0414 |  1.7776 |
| ftp_sp_k0_lognormal.15 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sigma_y1m1                   |  1.3232 |  1.3230 |  0.0277 |  1.2789 |  1.3699 |
| ftp_sp_k0_lognormal.16 | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | sigma_y2                     |  1.1016 |  1.1005 |  0.0349 |  1.0465 |  1.1601 |
| ftp_sp_k0_student.1    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z1_logL_Intercept          | -2.1476 | -2.1515 |  0.8388 | -3.5526 | -0.7733 |
| ftp_sp_k0_student.2    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z1_logA_Intercept          | -0.9518 | -0.9532 |  0.7150 | -2.1127 |  0.2433 |
| ftp_sp_k0_student.3    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z1_logk_Intercept          |  0.8609 |  0.7139 |  1.2238 | -1.0028 |  2.9388 |
| ftp_sp_k0_student.4    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z2_logL_Intercept          | -2.2033 | -2.2100 |  0.8354 | -3.5887 | -0.8286 |
| ftp_sp_k0_student.5    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z2_logA_Intercept          | -0.8125 | -0.8248 |  0.7168 | -1.9841 |  0.3796 |
| ftp_sp_k0_student.6    | ftp_sp_k0_student   | ftp_sp_k0       | student   | b_z2_logk_Intercept          |  0.5734 |  0.5149 |  1.0527 | -1.0628 |  2.4285 |
| ftp_sp_k0_student.7    | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h1\_\_z1_logL_Intercept   |  4.3493 |  4.0853 |  1.6567 |  2.1397 |  7.4726 |
| ftp_sp_k0_student.8    | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6500 |  0.4731 |  0.6147 |  0.0364 |  1.8585 |
| ftp_sp_k0_student.9    | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h1\_\_z1_logA_Intercept   |  3.8372 |  3.7987 |  2.0045 |  0.3971 |  7.2228 |
| ftp_sp_k0_student.10   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h2\_\_z1_logA_Intercept   |  0.6336 |  0.4625 |  0.6039 |  0.0382 |  1.8018 |
| ftp_sp_k0_student.11   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h1\_\_z2_logL_Intercept   |  3.9384 |  3.7319 |  1.5710 |  1.7737 |  6.7485 |
| ftp_sp_k0_student.12   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6332 |  0.4619 |  0.6005 |  0.0345 |  1.8229 |
| ftp_sp_k0_student.13   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h1\_\_z2_logA_Intercept   |  4.0617 |  3.9137 |  1.7778 |  1.3592 |  7.2109 |
| ftp_sp_k0_student.14   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sd_h2\_\_z2_logA_Intercept   |  0.6538 |  0.4839 |  0.6082 |  0.0379 |  1.8516 |
| ftp_sp_k0_student.15   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sigma_z1                     |  1.3095 |  1.3092 |  0.0285 |  1.2637 |  1.3569 |
| ftp_sp_k0_student.16   | ftp_sp_k0_student   | ftp_sp_k0       | student   | sigma_z2                     |  1.0897 |  1.0889 |  0.0349 |  1.0336 |  1.1493 |
| ftp_sp_k0_student.17   | ftp_sp_k0_student   | ftp_sp_k0       | student   | nu_z1                        | 62.7087 | 59.8478 | 19.8363 | 35.6672 | 98.8670 |
| ftp_sp_k0_student.18   | ftp_sp_k0_student   | ftp_sp_k0       | student   | nu_z2                        | 51.1613 | 48.2480 | 17.6875 | 27.7130 | 84.0348 |
| ftp_sp_k1_gamma.1      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y1m1_logL_Intercept        | -1.3919 | -1.2812 |  0.4604 | -2.3129 | -0.8563 |
| ftp_sp_k1_gamma.2      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y1m1_logA_Intercept        | -0.1457 | -0.0575 |  0.4574 | -1.0440 |  0.4393 |
| ftp_sp_k1_gamma.3      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y1m1_logk_Intercept        |  1.4397 |  1.6781 |  0.6636 |  0.0510 |  2.1234 |
| ftp_sp_k1_gamma.4      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y2_logL_Intercept          | -1.6396 | -1.5048 |  0.6228 | -2.8540 | -0.8753 |
| ftp_sp_k1_gamma.5      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y2_logA_Intercept          | -0.9445 | -0.9472 |  0.3958 | -1.5713 | -0.3097 |
| ftp_sp_k1_gamma.6      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | b_y2_logk_Intercept          |  0.5501 |  0.5524 |  0.6638 | -0.5205 |  1.6191 |
| ftp_sp_k1_gamma.7      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.5884 |  0.4104 |  0.5735 |  0.0386 |  1.6915 |
| ftp_sp_k1_gamma.8      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h2\_\_y1m1_logL_Intercept |  0.5476 |  0.5341 |  0.0997 |  0.4100 |  0.7306 |
| ftp_sp_k1_gamma.9      | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.5884 |  0.4120 |  0.5841 |  0.0351 |  1.7360 |
| ftp_sp_k1_gamma.10     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h2\_\_y1m1_logA_Intercept |  0.3476 |  0.3318 |  0.1164 |  0.1892 |  0.5615 |
| ftp_sp_k1_gamma.11     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y1m1_logk_Intercept |  0.7191 |  0.4571 |  0.7504 |  0.0289 |  2.2377 |
| ftp_sp_k1_gamma.12     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.6518 |  0.4674 |  0.6210 |  0.0492 |  1.9026 |
| ftp_sp_k1_gamma.13     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h2\_\_y2_logL_Intercept   |  0.4703 |  0.4113 |  0.2340 |  0.2579 |  0.8680 |
| ftp_sp_k1_gamma.14     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.4939 |  0.3375 |  0.5227 |  0.0288 |  1.4737 |
| ftp_sp_k1_gamma.15     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h2\_\_y2_logA_Intercept   |  0.6659 |  0.5348 |  0.4973 |  0.0685 |  1.6428 |
| ftp_sp_k1_gamma.16     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | sd_h1\_\_y2_logk_Intercept   |  0.7302 |  0.5390 |  0.6791 |  0.0432 |  2.0341 |
| ftp_sp_k1_gamma.17     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | shape_y1m1                   |  3.3905 |  3.3870 |  0.1367 |  3.1710 |  3.6202 |
| ftp_sp_k1_gamma.18     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | shape_y2                     |  5.5824 |  5.5687 |  0.3710 |  4.9944 |  6.2165 |
| ftp_sp_k1_lognormal.1  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y1m1_logL_Intercept        | -2.9021 | -2.9016 |  0.8457 | -4.3032 | -1.5280 |
| ftp_sp_k1_lognormal.2  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y1m1_logA_Intercept        | -1.1410 | -1.1338 |  0.6959 | -2.3021 |  0.0046 |
| ftp_sp_k1_lognormal.3  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y1m1_logk_Intercept        |  0.6575 |  0.6411 |  0.8872 | -0.7925 |  2.1298 |
| ftp_sp_k1_lognormal.4  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y2_logL_Intercept          | -2.9579 | -2.9512 |  0.8603 | -4.3826 | -1.5608 |
| ftp_sp_k1_lognormal.5  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y2_logA_Intercept          | -1.2491 | -1.2422 |  0.7220 | -2.4310 | -0.0617 |
| ftp_sp_k1_lognormal.6  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | b_y2_logk_Intercept          |  0.5893 |  0.5913 |  0.8693 | -0.8344 |  2.0363 |
| ftp_sp_k1_lognormal.7  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y1m1_logL_Intercept |  3.8506 |  3.5948 |  1.6010 |  1.7089 |  6.8056 |
| ftp_sp_k1_lognormal.8  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h2\_\_y1m1_logL_Intercept |  0.6427 |  0.4625 |  0.6107 |  0.0361 |  1.8547 |
| ftp_sp_k1_lognormal.9  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y1m1_logA_Intercept |  1.4679 |  0.9066 |  1.5394 |  0.0617 |  4.6749 |
| ftp_sp_k1_lognormal.10 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h2\_\_y1m1_logA_Intercept |  0.7303 |  0.5177 |  0.7304 |  0.0402 |  2.1296 |
| ftp_sp_k1_lognormal.11 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y1m1_logk_Intercept |  2.7530 |  2.6036 |  1.5870 |  0.3767 |  5.5966 |
| ftp_sp_k1_lognormal.12 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y2_logL_Intercept   |  3.4463 |  3.2161 |  1.6039 |  1.2409 |  6.3898 |
| ftp_sp_k1_lognormal.13 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h2\_\_y2_logL_Intercept   |  0.6247 |  0.4592 |  0.5851 |  0.0355 |  1.7478 |
| ftp_sp_k1_lognormal.14 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y2_logA_Intercept   |  1.9254 |  1.3940 |  1.7621 |  0.0824 |  5.3560 |
| ftp_sp_k1_lognormal.15 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h2\_\_y2_logA_Intercept   |  0.7438 |  0.5285 |  0.7439 |  0.0386 |  2.2180 |
| ftp_sp_k1_lognormal.16 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sd_h1\_\_y2_logk_Intercept   |  2.5471 |  2.3685 |  1.7252 |  0.1854 |  5.6567 |
| ftp_sp_k1_lognormal.17 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sigma_y1m1                   |  1.3224 |  1.3217 |  0.0276 |  1.2775 |  1.3687 |
| ftp_sp_k1_lognormal.18 | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | sigma_y2                     |  1.1022 |  1.1008 |  0.0343 |  1.0483 |  1.1606 |
| ftp_sp_k1_student.1    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z1_logL_Intercept          | -2.1388 | -2.1462 |  0.8072 | -3.4485 | -0.8163 |
| ftp_sp_k1_student.2    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z1_logA_Intercept          | -0.7430 | -0.7343 |  0.6938 | -1.8982 |  0.3917 |
| ftp_sp_k1_student.3    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z1_logk_Intercept          |  0.6722 |  0.6554 |  0.8935 | -0.7802 |  2.1875 |
| ftp_sp_k1_student.4    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z2_logL_Intercept          | -2.1952 | -2.2036 |  0.8515 | -3.6146 | -0.8047 |
| ftp_sp_k1_student.5    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z2_logA_Intercept          | -0.6851 | -0.6914 |  0.7146 | -1.8580 |  0.4827 |
| ftp_sp_k1_student.6    | ftp_sp_k1_student   | ftp_sp_k1       | student   | b_z2_logk_Intercept          |  0.5696 |  0.5781 |  0.8615 | -0.8343 |  2.0037 |
| ftp_sp_k1_student.7    | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z1_logL_Intercept   |  4.3322 |  4.0769 |  1.6259 |  2.1930 |  7.3557 |
| ftp_sp_k1_student.8    | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6481 |  0.4823 |  0.5923 |  0.0381 |  1.8261 |
| ftp_sp_k1_student.9    | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z1_logA_Intercept   |  1.4412 |  0.8599 |  1.5646 |  0.0648 |  4.7791 |
| ftp_sp_k1_student.10   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h2\_\_z1_logA_Intercept   |  0.7316 |  0.5240 |  0.7187 |  0.0375 |  2.1636 |
| ftp_sp_k1_student.11   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z1_logk_Intercept   |  2.8920 |  2.7332 |  1.6012 |  0.4447 |  5.8059 |
| ftp_sp_k1_student.12   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z2_logL_Intercept   |  3.9703 |  3.7371 |  1.5976 |  1.8212 |  6.8787 |
| ftp_sp_k1_student.13   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6358 |  0.4734 |  0.6011 |  0.0371 |  1.8188 |
| ftp_sp_k1_student.14   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z2_logA_Intercept   |  1.9928 |  1.3809 |  1.8792 |  0.0776 |  5.6499 |
| ftp_sp_k1_student.15   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h2\_\_z2_logA_Intercept   |  0.7700 |  0.5414 |  0.7737 |  0.0397 |  2.3022 |
| ftp_sp_k1_student.16   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sd_h1\_\_z2_logk_Intercept   |  2.8043 |  2.6774 |  1.8046 |  0.2169 |  5.9940 |
| ftp_sp_k1_student.17   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sigma_z1                     |  1.3090 |  1.3086 |  0.0283 |  1.2632 |  1.3562 |
| ftp_sp_k1_student.18   | ftp_sp_k1_student   | ftp_sp_k1       | student   | sigma_z2                     |  1.0899 |  1.0885 |  0.0348 |  1.0334 |  1.1488 |
| ftp_sp_k1_student.19   | ftp_sp_k1_student   | ftp_sp_k1       | student   | nu_z1                        | 62.4257 | 59.7565 | 19.5889 | 35.6704 | 98.6340 |
| ftp_sp_k1_student.20   | ftp_sp_k1_student   | ftp_sp_k1       | student   | nu_z2                        | 50.9123 | 48.2179 | 17.7722 | 26.9099 | 83.8738 |
| ftp_sp_k2_gamma.1      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y1m1_logL_Intercept        | -1.4041 | -1.2965 |  0.4567 | -2.3166 | -0.8695 |
| ftp_sp_k2_gamma.2      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y1m1_logA_Intercept        | -0.1786 | -0.0872 |  0.4868 | -1.1201 |  0.4537 |
| ftp_sp_k2_gamma.3      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y1m1_logk_Intercept        |  1.3768 |  1.5866 |  0.6819 |  0.0043 |  2.1335 |
| ftp_sp_k2_gamma.4      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y2_logL_Intercept          | -1.8945 | -1.8068 |  0.6797 | -3.1356 | -0.9671 |
| ftp_sp_k2_gamma.5      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y2_logA_Intercept          | -0.8663 | -0.8691 |  0.3905 | -1.4789 | -0.2142 |
| ftp_sp_k2_gamma.6      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | b_y2_logk_Intercept          |  0.2561 |  0.1909 |  0.6819 | -0.7471 |  1.4659 |
| ftp_sp_k2_gamma.7      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.5812 |  0.4035 |  0.5864 |  0.0436 |  1.7099 |
| ftp_sp_k2_gamma.8      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y1m1_logL_Intercept |  0.5426 |  0.5306 |  0.0996 |  0.4040 |  0.7224 |
| ftp_sp_k2_gamma.9      | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.6667 |  0.5043 |  0.6064 |  0.0546 |  1.8471 |
| ftp_sp_k2_gamma.10     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y1m1_logA_Intercept |  0.2963 |  0.2833 |  0.1473 |  0.0689 |  0.5562 |
| ftp_sp_k2_gamma.11     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y1m1_logk_Intercept |  0.7757 |  0.5416 |  0.7633 |  0.0338 |  2.3038 |
| ftp_sp_k2_gamma.12     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y1m1_logk_Intercept |  0.2705 |  0.2545 |  0.1596 |  0.0393 |  0.5570 |
| ftp_sp_k2_gamma.13     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.7281 |  0.5227 |  0.7048 |  0.0568 |  2.0969 |
| ftp_sp_k2_gamma.14     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y2_logL_Intercept   |  0.4693 |  0.4085 |  0.2721 |  0.1247 |  1.0024 |
| ftp_sp_k2_gamma.15     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.5415 |  0.3127 |  0.6419 |  0.0233 |  1.8783 |
| ftp_sp_k2_gamma.16     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y2_logA_Intercept   |  0.4427 |  0.3471 |  0.3642 |  0.0487 |  1.1917 |
| ftp_sp_k2_gamma.17     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h1\_\_y2_logk_Intercept   |  0.7236 |  0.4961 |  0.7395 |  0.0367 |  2.2474 |
| ftp_sp_k2_gamma.18     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | sd_h2\_\_y2_logk_Intercept   |  0.4852 |  0.4450 |  0.3259 |  0.0484 |  1.0570 |
| ftp_sp_k2_gamma.19     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | shape_y1m1                   |  3.4051 |  3.4050 |  0.1381 |  3.1820 |  3.6368 |
| ftp_sp_k2_gamma.20     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | shape_y2                     |  5.5742 |  5.5574 |  0.3652 |  5.0029 |  6.2140 |
| ftp_sp_k2_lognormal.1  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y1m1_logL_Intercept        | -2.8928 | -2.9002 |  0.8225 | -4.2465 | -1.5322 |
| ftp_sp_k2_lognormal.2  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y1m1_logA_Intercept        | -1.1486 | -1.1404 |  0.6958 | -2.2800 | -0.0212 |
| ftp_sp_k2_lognormal.3  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y1m1_logk_Intercept        |  0.6173 |  0.6133 |  0.8644 | -0.7777 |  2.0687 |
| ftp_sp_k2_lognormal.4  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y2_logL_Intercept          | -2.9704 | -2.9498 |  0.8617 | -4.4263 | -1.5805 |
| ftp_sp_k2_lognormal.5  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y2_logA_Intercept          | -1.2613 | -1.2697 |  0.7217 | -2.4635 | -0.0672 |
| ftp_sp_k2_lognormal.6  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | b_y2_logk_Intercept          |  0.5159 |  0.5059 |  0.8551 | -0.8803 |  1.9524 |
| ftp_sp_k2_lognormal.7  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y1m1_logL_Intercept |  3.8794 |  3.6494 |  1.6132 |  1.7041 |  6.8088 |
| ftp_sp_k2_lognormal.8  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y1m1_logL_Intercept |  0.6310 |  0.4621 |  0.5947 |  0.0351 |  1.7903 |
| ftp_sp_k2_lognormal.9  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y1m1_logA_Intercept |  1.6940 |  1.0633 |  1.7119 |  0.0688 |  5.2200 |
| ftp_sp_k2_lognormal.10 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y1m1_logA_Intercept |  0.7629 |  0.5361 |  0.7524 |  0.0435 |  2.2305 |
| ftp_sp_k2_lognormal.11 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y1m1_logk_Intercept |  2.7829 |  2.6262 |  1.6686 |  0.2659 |  5.7637 |
| ftp_sp_k2_lognormal.12 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y1m1_logk_Intercept |  0.6276 |  0.4121 |  0.6990 |  0.0299 |  1.9718 |
| ftp_sp_k2_lognormal.13 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y2_logL_Intercept   |  3.4284 |  3.2219 |  1.6184 |  1.1964 |  6.3962 |
| ftp_sp_k2_lognormal.14 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y2_logL_Intercept   |  0.6099 |  0.4460 |  0.5738 |  0.0349 |  1.7329 |
| ftp_sp_k2_lognormal.15 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y2_logA_Intercept   |  2.1211 |  1.6519 |  1.8478 |  0.1013 |  5.6569 |
| ftp_sp_k2_lognormal.16 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y2_logA_Intercept   |  0.7650 |  0.5438 |  0.7603 |  0.0398 |  2.2632 |
| ftp_sp_k2_lognormal.17 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h1\_\_y2_logk_Intercept   |  2.5734 |  2.3751 |  1.8303 |  0.1621 |  5.9094 |
| ftp_sp_k2_lognormal.18 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sd_h2\_\_y2_logk_Intercept   |  0.7302 |  0.4968 |  0.7631 |  0.0364 |  2.1857 |
| ftp_sp_k2_lognormal.19 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sigma_y1m1                   |  1.3224 |  1.3222 |  0.0277 |  1.2782 |  1.3685 |
| ftp_sp_k2_lognormal.20 | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | sigma_y2                     |  1.1002 |  1.0988 |  0.0337 |  1.0470 |  1.1575 |
| ftp_sp_k2_student.1    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z1_logL_Intercept          | -2.1458 | -2.1515 |  0.8249 | -3.4962 | -0.7705 |
| ftp_sp_k2_student.2    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z1_logA_Intercept          | -0.7561 | -0.7484 |  0.6967 | -1.9095 |  0.3811 |
| ftp_sp_k2_student.3    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z1_logk_Intercept          |  0.5996 |  0.5953 |  0.8877 | -0.8505 |  2.0787 |
| ftp_sp_k2_student.4    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z2_logL_Intercept          | -2.1953 | -2.1843 |  0.8401 | -3.5751 | -0.8388 |
| ftp_sp_k2_student.5    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z2_logA_Intercept          | -0.6811 | -0.6851 |  0.6995 | -1.8245 |  0.4647 |
| ftp_sp_k2_student.6    | ftp_sp_k2_student   | ftp_sp_k2       | student   | b_z2_logk_Intercept          |  0.5054 |  0.4916 |  0.8738 | -0.9371 |  1.9898 |
| ftp_sp_k2_student.7    | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z1_logL_Intercept   |  4.3014 |  4.0554 |  1.6283 |  2.1023 |  7.2713 |
| ftp_sp_k2_student.8    | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6442 |  0.4719 |  0.6008 |  0.0401 |  1.8398 |
| ftp_sp_k2_student.9    | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z1_logA_Intercept   |  1.7135 |  1.0111 |  1.7838 |  0.0668 |  5.4642 |
| ftp_sp_k2_student.10   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z1_logA_Intercept   |  0.7522 |  0.5346 |  0.7414 |  0.0429 |  2.2122 |
| ftp_sp_k2_student.11   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z1_logk_Intercept   |  2.9236 |  2.7937 |  1.7156 |  0.3267 |  5.9030 |
| ftp_sp_k2_student.12   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z1_logk_Intercept   |  0.6131 |  0.4271 |  0.6414 |  0.0339 |  1.8307 |
| ftp_sp_k2_student.13   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z2_logL_Intercept   |  3.9600 |  3.7195 |  1.5750 |  1.8498 |  6.9039 |
| ftp_sp_k2_student.14   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6270 |  0.4680 |  0.5824 |  0.0382 |  1.7858 |
| ftp_sp_k2_student.15   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z2_logA_Intercept   |  2.2664 |  1.7347 |  1.9832 |  0.1038 |  5.9901 |
| ftp_sp_k2_student.16   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z2_logA_Intercept   |  0.7661 |  0.5414 |  0.7638 |  0.0381 |  2.2765 |
| ftp_sp_k2_student.17   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h1\_\_z2_logk_Intercept   |  2.7217 |  2.5840 |  1.8931 |  0.1439 |  6.1179 |
| ftp_sp_k2_student.18   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sd_h2\_\_z2_logk_Intercept   |  0.7335 |  0.5116 |  0.7644 |  0.0371 |  2.2063 |
| ftp_sp_k2_student.19   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sigma_z1                     |  1.3092 |  1.3086 |  0.0283 |  1.2636 |  1.3566 |
| ftp_sp_k2_student.20   | ftp_sp_k2_student   | ftp_sp_k2       | student   | sigma_z2                     |  1.0896 |  1.0891 |  0.0345 |  1.0344 |  1.1468 |
| ftp_sp_k2_student.21   | ftp_sp_k2_student   | ftp_sp_k2       | student   | nu_z1                        | 62.6368 | 59.7379 | 19.3060 | 36.6983 | 97.6588 |
| ftp_sp_k2_student.22   | ftp_sp_k2_student   | ftp_sp_k2       | student   | nu_z2                        | 51.2173 | 48.3944 | 18.0682 | 27.1728 | 84.6200 |
| PFT_k0_student.1       | PFT_k0_student      | PFT_k0          | student   | b_z1_logL_Intercept          | -2.4668 | -2.4547 |  0.8761 | -3.9348 | -1.0592 |
| PFT_k0_student.2       | PFT_k0_student      | PFT_k0          | student   | b_z1_logA_Intercept          | -1.1717 | -1.1589 |  0.6751 | -2.3064 | -0.1045 |
| PFT_k0_student.3       | PFT_k0_student      | PFT_k0          | student   | b_z1_logk_Intercept          |  2.1886 |  2.5565 |  1.0791 | -0.1478 |  3.3265 |
| PFT_k0_student.4       | PFT_k0_student      | PFT_k0          | student   | b_z2_logL_Intercept          | -2.5639 | -2.5491 |  0.8863 | -4.0738 | -1.1303 |
| PFT_k0_student.5       | PFT_k0_student      | PFT_k0          | student   | b_z2_logA_Intercept          | -1.1376 | -1.1409 |  0.7405 | -2.3446 |  0.0913 |
| PFT_k0_student.6       | PFT_k0_student      | PFT_k0          | student   | b_z2_logk_Intercept          |  1.2922 |  1.3633 |  1.1913 | -0.7028 |  3.0182 |
| PFT_k0_student.7       | PFT_k0_student      | PFT_k0          | student   | sd_h1\_\_z1_logL_Intercept   |  4.4497 |  4.1977 |  1.7684 |  2.0124 |  7.7259 |
| PFT_k0_student.8       | PFT_k0_student      | PFT_k0          | student   | sd_h1\_\_z1_logA_Intercept   |  1.8864 |  0.9869 |  2.0600 |  0.0567 |  6.1238 |
| PFT_k0_student.9       | PFT_k0_student      | PFT_k0          | student   | sd_h1\_\_z2_logL_Intercept   |  3.8332 |  3.6471 |  1.7566 |  1.3074 |  6.9895 |
| PFT_k0_student.10      | PFT_k0_student      | PFT_k0          | student   | sd_h1\_\_z2_logA_Intercept   |  3.1854 |  3.1091 |  2.1995 |  0.1679 |  7.0431 |
| PFT_k0_student.11      | PFT_k0_student      | PFT_k0          | student   | sigma_z1                     |  1.3104 |  1.3097 |  0.0287 |  1.2649 |  1.3582 |
| PFT_k0_student.12      | PFT_k0_student      | PFT_k0          | student   | sigma_z2                     |  1.0920 |  1.0913 |  0.0346 |  1.0362 |  1.1497 |
| PFT_k0_student.13      | PFT_k0_student      | PFT_k0          | student   | nu_z1                        | 62.7883 | 60.2879 | 19.2929 | 36.2454 | 98.6583 |
| PFT_k0_student.14      | PFT_k0_student      | PFT_k0          | student   | nu_z2                        | 51.6230 | 48.6571 | 18.4485 | 27.2583 | 86.5070 |
| PFT_k1_student.1       | PFT_k1_student      | PFT_k1          | student   | b_z1_logL_Intercept          | -2.4645 | -2.4369 |  0.8455 | -3.8886 | -1.1185 |
| PFT_k1_student.2       | PFT_k1_student      | PFT_k1          | student   | b_z1_logA_Intercept          | -0.9144 | -0.9185 |  0.6999 | -2.0561 |  0.2339 |
| PFT_k1_student.3       | PFT_k1_student      | PFT_k1          | student   | b_z1_logk_Intercept          |  1.1751 |  1.1461 |  0.9944 | -0.4288 |  2.8252 |
| PFT_k1_student.4       | PFT_k1_student      | PFT_k1          | student   | b_z2_logL_Intercept          | -2.5738 | -2.5733 |  0.8941 | -4.0638 | -1.1115 |
| PFT_k1_student.5       | PFT_k1_student      | PFT_k1          | student   | b_z2_logA_Intercept          | -0.9173 | -0.9063 |  0.7367 | -2.1509 |  0.2884 |
| PFT_k1_student.6       | PFT_k1_student      | PFT_k1          | student   | b_z2_logk_Intercept          |  0.9259 |  0.9244 |  0.9371 | -0.6204 |  2.4903 |
| PFT_k1_student.7       | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z1_logL_Intercept   |  4.4462 |  4.2233 |  1.7405 |  2.0286 |  7.6168 |
| PFT_k1_student.8       | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z1_logA_Intercept   |  1.4595 |  0.9008 |  1.5529 |  0.0569 |  4.7657 |
| PFT_k1_student.9       | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z1_logk_Intercept   |  2.6524 |  2.4474 |  1.7780 |  0.1944 |  5.9125 |
| PFT_k1_student.10      | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z2_logL_Intercept   |  3.8271 |  3.6185 |  1.7479 |  1.3258 |  7.0003 |
| PFT_k1_student.11      | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z2_logA_Intercept   |  2.0657 |  1.5065 |  1.8651 |  0.1032 |  5.7568 |
| PFT_k1_student.12      | PFT_k1_student      | PFT_k1          | student   | sd_h1\_\_z2_logk_Intercept   |  2.6926 |  2.4952 |  1.8919 |  0.1826 |  6.1288 |
| PFT_k1_student.13      | PFT_k1_student      | PFT_k1          | student   | sigma_z1                     |  1.3096 |  1.3091 |  0.0289 |  1.2629 |  1.3571 |
| PFT_k1_student.14      | PFT_k1_student      | PFT_k1          | student   | sigma_z2                     |  1.0927 |  1.0923 |  0.0354 |  1.0361 |  1.1523 |
| PFT_k1_student.15      | PFT_k1_student      | PFT_k1          | student   | nu_z1                        | 62.5221 | 59.5742 | 19.4956 | 35.9466 | 98.4416 |
| PFT_k1_student.16      | PFT_k1_student      | PFT_k1          | student   | nu_z2                        | 51.3066 | 48.5986 | 18.2920 | 27.2108 | 84.9494 |
| PFT_sp_k0_student.1    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z1_logL_Intercept          | -2.4423 | -2.4271 |  0.8444 | -3.8402 | -1.0846 |
| PFT_sp_k0_student.2    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z1_logA_Intercept          | -1.1821 | -1.1816 |  0.6939 | -2.3489 | -0.0514 |
| PFT_sp_k0_student.3    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z1_logk_Intercept          |  2.1182 |  2.5160 |  1.1210 | -0.2451 |  3.3335 |
| PFT_sp_k0_student.4    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z2_logL_Intercept          | -2.5390 | -2.5366 |  0.8790 | -4.0227 | -1.1036 |
| PFT_sp_k0_student.5    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z2_logA_Intercept          | -1.0991 | -1.1130 |  0.7423 | -2.3083 |  0.1252 |
| PFT_sp_k0_student.6    | PFT_sp_k0_student   | PFT_sp_k0       | student   | b_z2_logk_Intercept          |  1.2128 |  1.2652 |  1.2079 | -0.7975 |  2.9970 |
| PFT_sp_k0_student.7    | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h1\_\_z1_logL_Intercept   |  4.6019 |  4.3770 |  1.7604 |  2.1666 |  7.7729 |
| PFT_sp_k0_student.8    | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6473 |  0.4795 |  0.5986 |  0.0367 |  1.8266 |
| PFT_sp_k0_student.9    | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h1\_\_z1_logA_Intercept   |  2.0721 |  1.1611 |  2.1475 |  0.0671 |  6.3827 |
| PFT_sp_k0_student.10   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h2\_\_z1_logA_Intercept   |  0.5663 |  0.4249 |  0.5211 |  0.0323 |  1.6210 |
| PFT_sp_k0_student.11   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h1\_\_z2_logL_Intercept   |  4.0024 |  3.8060 |  1.7705 |  1.4404 |  7.1792 |
| PFT_sp_k0_student.12   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6185 |  0.4634 |  0.5745 |  0.0384 |  1.7418 |
| PFT_sp_k0_student.13   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h1\_\_z2_logA_Intercept   |  3.4214 |  3.3816 |  2.1780 |  0.2346 |  7.1991 |
| PFT_sp_k0_student.14   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sd_h2\_\_z2_logA_Intercept   |  0.6176 |  0.4501 |  0.5742 |  0.0350 |  1.7658 |
| PFT_sp_k0_student.15   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sigma_z1                     |  1.3104 |  1.3098 |  0.0287 |  1.2641 |  1.3588 |
| PFT_sp_k0_student.16   | PFT_sp_k0_student   | PFT_sp_k0       | student   | sigma_z2                     |  1.0926 |  1.0919 |  0.0346 |  1.0367 |  1.1503 |
| PFT_sp_k0_student.17   | PFT_sp_k0_student   | PFT_sp_k0       | student   | nu_z1                        | 62.4946 | 60.1821 | 19.1834 | 35.9279 | 97.6561 |
| PFT_sp_k0_student.18   | PFT_sp_k0_student   | PFT_sp_k0       | student   | nu_z2                        | 51.0765 | 48.3065 | 17.7106 | 27.4437 | 82.9073 |
| PFT_sp_k1_student.1    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z1_logL_Intercept          | -2.4274 | -2.4102 |  0.8599 | -3.8760 | -1.0578 |
| PFT_sp_k1_student.2    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z1_logA_Intercept          | -0.9062 | -0.9003 |  0.7113 | -2.0946 |  0.2440 |
| PFT_sp_k1_student.3    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z1_logk_Intercept          |  1.1392 |  1.1080 |  0.9624 | -0.4154 |  2.7757 |
| PFT_sp_k1_student.4    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z2_logL_Intercept          | -2.5437 | -2.5296 |  0.8903 | -4.0403 | -1.0990 |
| PFT_sp_k1_student.5    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z2_logA_Intercept          | -0.8757 | -0.8786 |  0.7268 | -2.0647 |  0.3275 |
| PFT_sp_k1_student.6    | PFT_sp_k1_student   | PFT_sp_k1       | student   | b_z2_logk_Intercept          |  0.9172 |  0.8875 |  0.9517 | -0.6230 |  2.5427 |
| PFT_sp_k1_student.7    | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z1_logL_Intercept   |  4.6068 |  4.3749 |  1.7556 |  2.1696 |  7.8045 |
| PFT_sp_k1_student.8    | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6438 |  0.4807 |  0.5913 |  0.0400 |  1.8244 |
| PFT_sp_k1_student.9    | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z1_logA_Intercept   |  1.4297 |  0.8859 |  1.5496 |  0.0589 |  4.6548 |
| PFT_sp_k1_student.10   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h2\_\_z1_logA_Intercept   |  0.7342 |  0.5117 |  0.7501 |  0.0383 |  2.1960 |
| PFT_sp_k1_student.11   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z1_logk_Intercept   |  2.8306 |  2.6714 |  1.7829 |  0.2908 |  5.9610 |
| PFT_sp_k1_student.12   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z2_logL_Intercept   |  3.9557 |  3.7696 |  1.7281 |  1.4485 |  7.0895 |
| PFT_sp_k1_student.13   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6191 |  0.4509 |  0.5876 |  0.0321 |  1.7698 |
| PFT_sp_k1_student.14   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z2_logA_Intercept   |  1.9544 |  1.3670 |  1.8356 |  0.0908 |  5.6015 |
| PFT_sp_k1_student.15   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h2\_\_z2_logA_Intercept   |  0.7885 |  0.5608 |  0.7657 |  0.0446 |  2.2927 |
| PFT_sp_k1_student.16   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sd_h1\_\_z2_logk_Intercept   |  2.8938 |  2.6923 |  1.9244 |  0.2396 |  6.4140 |
| PFT_sp_k1_student.17   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sigma_z1                     |  1.3099 |  1.3094 |  0.0288 |  1.2640 |  1.3575 |
| PFT_sp_k1_student.18   | PFT_sp_k1_student   | PFT_sp_k1       | student   | sigma_z2                     |  1.0921 |  1.0910 |  0.0352 |  1.0356 |  1.1522 |
| PFT_sp_k1_student.19   | PFT_sp_k1_student   | PFT_sp_k1       | student   | nu_z1                        | 62.4537 | 59.6943 | 19.1734 | 36.3135 | 97.8293 |
| PFT_sp_k1_student.20   | PFT_sp_k1_student   | PFT_sp_k1       | student   | nu_z2                        | 50.9519 | 48.2440 | 17.8940 | 27.0798 | 84.1973 |
| PFT_sp_k2_student.1    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z1_logL_Intercept          | -2.4226 | -2.4215 |  0.8424 | -3.7997 | -1.0308 |
| PFT_sp_k2_student.2    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z1_logA_Intercept          | -0.9099 | -0.8990 |  0.6938 | -2.0574 |  0.2137 |
| PFT_sp_k2_student.3    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z1_logk_Intercept          |  1.0087 |  0.9731 |  0.9684 | -0.5101 |  2.6991 |
| PFT_sp_k2_student.4    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z2_logL_Intercept          | -2.5173 | -2.5021 |  0.8751 | -4.0063 | -1.0919 |
| PFT_sp_k2_student.5    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z2_logA_Intercept          | -0.8549 | -0.8508 |  0.7402 | -2.0802 |  0.3485 |
| PFT_sp_k2_student.6    | PFT_sp_k2_student   | PFT_sp_k2       | student   | b_z2_logk_Intercept          |  0.8039 |  0.7885 |  0.9282 | -0.6935 |  2.3541 |
| PFT_sp_k2_student.7    | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z1_logL_Intercept   |  4.6035 |  4.3872 |  1.7280 |  2.2249 |  7.8156 |
| PFT_sp_k2_student.8    | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z1_logL_Intercept   |  0.6443 |  0.4808 |  0.5990 |  0.0393 |  1.8328 |
| PFT_sp_k2_student.9    | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z1_logA_Intercept   |  1.7025 |  1.0329 |  1.7865 |  0.0721 |  5.4707 |
| PFT_sp_k2_student.10   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z1_logA_Intercept   |  0.7764 |  0.5435 |  0.7843 |  0.0408 |  2.3347 |
| PFT_sp_k2_student.11   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z1_logk_Intercept   |  2.9332 |  2.8295 |  1.8448 |  0.2367 |  6.1526 |
| PFT_sp_k2_student.12   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z1_logk_Intercept   |  0.6039 |  0.4136 |  0.6395 |  0.0334 |  1.8287 |
| PFT_sp_k2_student.13   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z2_logL_Intercept   |  4.0227 |  3.8434 |  1.7233 |  1.5235 |  7.1094 |
| PFT_sp_k2_student.14   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z2_logL_Intercept   |  0.6181 |  0.4532 |  0.5812 |  0.0384 |  1.7543 |
| PFT_sp_k2_student.15   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z2_logA_Intercept   |  2.2097 |  1.6683 |  1.9686 |  0.0913 |  6.0004 |
| PFT_sp_k2_student.16   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z2_logA_Intercept   |  0.7977 |  0.5635 |  0.7861 |  0.0481 |  2.3820 |
| PFT_sp_k2_student.17   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h1\_\_z2_logk_Intercept   |  2.8825 |  2.6946 |  1.9743 |  0.1875 |  6.4823 |
| PFT_sp_k2_student.18   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sd_h2\_\_z2_logk_Intercept   |  0.7443 |  0.5248 |  0.7483 |  0.0423 |  2.2524 |
| PFT_sp_k2_student.19   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sigma_z1                     |  1.3097 |  1.3094 |  0.0286 |  1.2628 |  1.3575 |
| PFT_sp_k2_student.20   | PFT_sp_k2_student   | PFT_sp_k2       | student   | sigma_z2                     |  1.0915 |  1.0901 |  0.0352 |  1.0353 |  1.1511 |
| PFT_sp_k2_student.21   | PFT_sp_k2_student   | PFT_sp_k2       | student   | nu_z1                        | 62.7608 | 59.8006 | 19.6335 | 36.0763 | 98.7409 |
| PFT_sp_k2_student.22   | PFT_sp_k2_student   | PFT_sp_k2       | student   | nu_z2                        | 51.5190 | 48.8869 | 17.8488 | 27.2386 | 84.9948 |
| sp_k0_gamma.1          | sp_k0_gamma         | sp_k0           | gamma     | b_y1m1_logL_Intercept        | -1.1258 | -1.1253 |  0.1254 | -1.3329 | -0.9249 |
| sp_k0_gamma.2          | sp_k0_gamma         | sp_k0           | gamma     | b_y1m1_logA_Intercept        |  0.1441 |  0.1448 |  0.1756 | -0.1460 |  0.4295 |
| sp_k0_gamma.3          | sp_k0_gamma         | sp_k0           | gamma     | b_y1m1_logk_Intercept        |  1.9556 |  1.9596 |  0.1111 |  1.7679 |  2.1317 |
| sp_k0_gamma.4          | sp_k0_gamma         | sp_k0           | gamma     | b_y2_logL_Intercept          | -1.4475 | -1.3019 |  0.4803 | -2.4554 | -0.9664 |
| sp_k0_gamma.5          | sp_k0_gamma         | sp_k0           | gamma     | b_y2_logA_Intercept          | -1.0139 | -1.0079 |  0.2360 | -1.4100 | -0.6446 |
| sp_k0_gamma.6          | sp_k0_gamma         | sp_k0           | gamma     | b_y2_logk_Intercept          |  0.6459 |  0.6253 |  0.5696 | -0.2334 |  1.5785 |
| sp_k0_gamma.7          | sp_k0_gamma         | sp_k0           | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.5603 |  0.5477 |  0.0972 |  0.4262 |  0.7381 |
| sp_k0_gamma.8          | sp_k0_gamma         | sp_k0           | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.3543 |  0.3389 |  0.1130 |  0.1994 |  0.5568 |
| sp_k0_gamma.9          | sp_k0_gamma         | sp_k0           | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.5255 |  0.4667 |  0.2337 |  0.2770 |  0.9813 |
| sp_k0_gamma.10         | sp_k0_gamma         | sp_k0           | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.4856 |  0.4346 |  0.3311 |  0.0467 |  1.0986 |
| sp_k0_gamma.11         | sp_k0_gamma         | sp_k0           | gamma     | shape_y1m1                   |  3.3933 |  3.3929 |  0.1354 |  3.1762 |  3.6157 |
| sp_k0_gamma.12         | sp_k0_gamma         | sp_k0           | gamma     | shape_y2                     |  5.5456 |  5.5359 |  0.3502 |  4.9833 |  6.1382 |
| sp_k0_lognormal.1      | sp_k0_lognormal     | sp_k0           | lognormal | b_y1m1_logL_Intercept        | -5.1735 | -5.1528 |  0.4330 | -5.9011 | -4.4954 |
| sp_k0_lognormal.2      | sp_k0_lognormal     | sp_k0           | lognormal | b_y1m1_logA_Intercept        | -1.6070 | -1.5931 |  0.6652 | -2.7211 | -0.5313 |
| sp_k0_lognormal.3      | sp_k0_lognormal     | sp_k0           | lognormal | b_y1m1_logk_Intercept        |  2.7070 |  2.7076 |  0.3818 |  2.0841 |  3.3348 |
| sp_k0_lognormal.4      | sp_k0_lognormal     | sp_k0           | lognormal | b_y2_logL_Intercept          | -4.7955 | -4.7786 |  0.4491 | -5.5678 | -4.0887 |
| sp_k0_lognormal.5      | sp_k0_lognormal     | sp_k0           | lognormal | b_y2_logA_Intercept          | -2.0082 | -1.9928 |  0.6399 | -3.0684 | -0.9584 |
| sp_k0_lognormal.6      | sp_k0_lognormal     | sp_k0           | lognormal | b_y2_logk_Intercept          |  2.4530 |  2.4656 |  0.4867 |  1.6430 |  3.2365 |
| sp_k0_lognormal.7      | sp_k0_lognormal     | sp_k0           | lognormal | sd_h1\_\_y1m1_logL_Intercept |  0.4122 |  0.3140 |  0.3724 |  0.0266 |  1.1212 |
| sp_k0_lognormal.8      | sp_k0_lognormal     | sp_k0           | lognormal | sd_h1\_\_y1m1_logA_Intercept |  0.5109 |  0.3920 |  0.4537 |  0.0337 |  1.4073 |
| sp_k0_lognormal.9      | sp_k0_lognormal     | sp_k0           | lognormal | sd_h1\_\_y2_logL_Intercept   |  0.4020 |  0.3014 |  0.3648 |  0.0240 |  1.1157 |
| sp_k0_lognormal.10     | sp_k0_lognormal     | sp_k0           | lognormal | sd_h1\_\_y2_logA_Intercept   |  0.5062 |  0.3755 |  0.4713 |  0.0286 |  1.4271 |
| sp_k0_lognormal.11     | sp_k0_lognormal     | sp_k0           | lognormal | sigma_y1m1                   |  1.3278 |  1.3273 |  0.0280 |  1.2826 |  1.3743 |
| sp_k0_lognormal.12     | sp_k0_lognormal     | sp_k0           | lognormal | sigma_y2                     |  1.1101 |  1.1094 |  0.0345 |  1.0548 |  1.1675 |
| sp_k0_student.1        | sp_k0_student       | sp_k0           | student   | b_z1_logL_Intercept          | -4.9845 | -4.9703 |  0.4104 | -5.6833 | -4.3464 |
| sp_k0_student.2        | sp_k0_student       | sp_k0           | student   | b_z1_logA_Intercept          | -1.2473 | -1.2296 |  0.6407 | -2.3307 | -0.2310 |
| sp_k0_student.3        | sp_k0_student       | sp_k0           | student   | b_z1_logk_Intercept          |  2.8266 |  2.8179 |  0.3675 |  2.2367 |  3.4448 |
| sp_k0_student.4        | sp_k0_student       | sp_k0           | student   | b_z2_logL_Intercept          | -4.5981 | -4.5814 |  0.4268 | -5.3306 | -3.9379 |
| sp_k0_student.5        | sp_k0_student       | sp_k0           | student   | b_z2_logA_Intercept          | -1.5556 | -1.5450 |  0.6312 | -2.6119 | -0.5300 |
| sp_k0_student.6        | sp_k0_student       | sp_k0           | student   | b_z2_logk_Intercept          |  2.6660 |  2.6585 |  0.4673 |  1.9207 |  3.4378 |
| sp_k0_student.7        | sp_k0_student       | sp_k0           | student   | sd_h1\_\_z1_logL_Intercept   |  0.4101 |  0.3053 |  0.4050 |  0.0259 |  1.1152 |
| sp_k0_student.8        | sp_k0_student       | sp_k0           | student   | sd_h1\_\_z1_logA_Intercept   |  0.4912 |  0.3612 |  0.4539 |  0.0323 |  1.3515 |
| sp_k0_student.9        | sp_k0_student       | sp_k0           | student   | sd_h1\_\_z2_logL_Intercept   |  0.4051 |  0.3027 |  0.3895 |  0.0258 |  1.1158 |
| sp_k0_student.10       | sp_k0_student       | sp_k0           | student   | sd_h1\_\_z2_logA_Intercept   |  0.5254 |  0.3807 |  0.4973 |  0.0297 |  1.5136 |
| sp_k0_student.11       | sp_k0_student       | sp_k0           | student   | sigma_z1                     |  1.3155 |  1.3151 |  0.0284 |  1.2695 |  1.3624 |
| sp_k0_student.12       | sp_k0_student       | sp_k0           | student   | sigma_z2                     |  1.1006 |  1.0999 |  0.0356 |  1.0432 |  1.1598 |
| sp_k0_student.13       | sp_k0_student       | sp_k0           | student   | nu_z1                        | 62.9835 | 60.1913 | 19.7465 | 35.9969 | 99.8172 |
| sp_k0_student.14       | sp_k0_student       | sp_k0           | student   | nu_z2                        | 51.6092 | 48.4438 | 18.1511 | 27.3412 | 86.0654 |
| sp_k1_gamma.1          | sp_k1_gamma         | sp_k1           | gamma     | b_y1m1_logL_Intercept        | -1.1349 | -1.1336 |  0.1293 | -1.3476 | -0.9247 |
| sp_k1_gamma.2          | sp_k1_gamma         | sp_k1           | gamma     | b_y1m1_logA_Intercept        |  0.1756 |  0.1801 |  0.1787 | -0.1224 |  0.4613 |
| sp_k1_gamma.3          | sp_k1_gamma         | sp_k1           | gamma     | b_y1m1_logk_Intercept        |  1.9541 |  1.9594 |  0.1361 |  1.7228 |  2.1684 |
| sp_k1_gamma.4          | sp_k1_gamma         | sp_k1           | gamma     | b_y2_logL_Intercept          | -1.7717 | -1.6352 |  0.6282 | -3.0060 | -1.0240 |
| sp_k1_gamma.5          | sp_k1_gamma         | sp_k1           | gamma     | b_y2_logA_Intercept          | -0.9131 | -0.8924 |  0.2265 | -1.3135 | -0.5773 |
| sp_k1_gamma.6          | sp_k1_gamma         | sp_k1           | gamma     | b_y2_logk_Intercept          |  0.2466 |  0.1348 |  0.6239 | -0.6172 |  1.4040 |
| sp_k1_gamma.7          | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y1m1_logL_Intercept |  0.5604 |  0.5495 |  0.0949 |  0.4271 |  0.7275 |
| sp_k1_gamma.8          | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y1m1_logA_Intercept |  0.3160 |  0.3042 |  0.1417 |  0.0973 |  0.5630 |
| sp_k1_gamma.9          | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y1m1_logk_Intercept |  0.2160 |  0.1994 |  0.1417 |  0.0236 |  0.4724 |
| sp_k1_gamma.10         | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y2_logL_Intercept   |  0.5495 |  0.4788 |  0.3186 |  0.1433 |  1.1708 |
| sp_k1_gamma.11         | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y2_logA_Intercept   |  0.4102 |  0.3624 |  0.2822 |  0.0472 |  0.9585 |
| sp_k1_gamma.12         | sp_k1_gamma         | sp_k1           | gamma     | sd_h1\_\_y2_logk_Intercept   |  0.4666 |  0.4296 |  0.3149 |  0.0471 |  1.0388 |
| sp_k1_gamma.13         | sp_k1_gamma         | sp_k1           | gamma     | shape_y1m1                   |  3.4027 |  3.4011 |  0.1367 |  3.1831 |  3.6325 |
| sp_k1_gamma.14         | sp_k1_gamma         | sp_k1           | gamma     | shape_y2                     |  5.5515 |  5.5434 |  0.3478 |  4.9900 |  6.1335 |
| sp_k1_lognormal.1      | sp_k1_lognormal     | sp_k1           | lognormal | b_y1m1_logL_Intercept        | -5.1679 | -5.1511 |  0.4296 | -5.9007 | -4.4998 |
| sp_k1_lognormal.2      | sp_k1_lognormal     | sp_k1           | lognormal | b_y1m1_logA_Intercept        | -1.6654 | -1.6438 |  0.6970 | -2.8400 | -0.5342 |
| sp_k1_lognormal.3      | sp_k1_lognormal     | sp_k1           | lognormal | b_y1m1_logk_Intercept        |  2.7189 |  2.7203 |  0.4284 |  2.0311 |  3.4158 |
| sp_k1_lognormal.4      | sp_k1_lognormal     | sp_k1           | lognormal | b_y2_logL_Intercept          | -4.7906 | -4.7679 |  0.4386 | -5.5452 | -4.0979 |
| sp_k1_lognormal.5      | sp_k1_lognormal     | sp_k1           | lognormal | b_y2_logA_Intercept          | -2.0834 | -2.0693 |  0.6692 | -3.2068 | -0.9981 |
| sp_k1_lognormal.6      | sp_k1_lognormal     | sp_k1           | lognormal | b_y2_logk_Intercept          |  2.4169 |  2.4465 |  0.5677 |  1.4751 |  3.2876 |
| sp_k1_lognormal.7      | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y1m1_logL_Intercept |  0.4094 |  0.3075 |  0.3830 |  0.0260 |  1.1208 |
| sp_k1_lognormal.8      | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y1m1_logA_Intercept |  0.5308 |  0.3940 |  0.4975 |  0.0328 |  1.4774 |
| sp_k1_lognormal.9      | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y1m1_logk_Intercept |  0.3055 |  0.2219 |  0.3170 |  0.0181 |  0.8553 |
| sp_k1_lognormal.10     | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y2_logL_Intercept   |  0.4136 |  0.3199 |  0.3739 |  0.0269 |  1.1331 |
| sp_k1_lognormal.11     | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y2_logA_Intercept   |  0.5393 |  0.4015 |  0.5180 |  0.0320 |  1.5008 |
| sp_k1_lognormal.12     | sp_k1_lognormal     | sp_k1           | lognormal | sd_h1\_\_y2_logk_Intercept   |  0.4292 |  0.2980 |  0.4768 |  0.0231 |  1.2570 |
| sp_k1_lognormal.13     | sp_k1_lognormal     | sp_k1           | lognormal | sigma_y1m1                   |  1.3275 |  1.3270 |  0.0283 |  1.2822 |  1.3750 |
| sp_k1_lognormal.14     | sp_k1_lognormal     | sp_k1           | lognormal | sigma_y2                     |  1.1103 |  1.1099 |  0.0342 |  1.0552 |  1.1678 |
| sp_k1_student.1        | sp_k1_student       | sp_k1           | student   | b_z1_logL_Intercept          | -4.9893 | -4.9770 |  0.4134 | -5.6923 | -4.3344 |
| sp_k1_student.2        | sp_k1_student       | sp_k1           | student   | b_z1_logA_Intercept          | -1.2828 | -1.2697 |  0.6675 | -2.4315 | -0.2274 |
| sp_k1_student.3        | sp_k1_student       | sp_k1           | student   | b_z1_logk_Intercept          |  2.8538 |  2.8525 |  0.3900 |  2.2269 |  3.5011 |
| sp_k1_student.4        | sp_k1_student       | sp_k1           | student   | b_z2_logL_Intercept          | -4.5953 | -4.5741 |  0.4280 | -5.3298 | -3.9198 |
| sp_k1_student.5        | sp_k1_student       | sp_k1           | student   | b_z2_logA_Intercept          | -1.5849 | -1.5635 |  0.6420 | -2.6485 | -0.5449 |
| sp_k1_student.6        | sp_k1_student       | sp_k1           | student   | b_z2_logk_Intercept          |  2.6623 |  2.6717 |  0.5102 |  1.8200 |  3.4874 |
| sp_k1_student.7        | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z1_logL_Intercept   |  0.4017 |  0.3043 |  0.3833 |  0.0246 |  1.0955 |
| sp_k1_student.8        | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z1_logA_Intercept   |  0.5248 |  0.3901 |  0.4857 |  0.0310 |  1.4786 |
| sp_k1_student.9        | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z1_logk_Intercept   |  0.2845 |  0.2157 |  0.2716 |  0.0157 |  0.7857 |
| sp_k1_student.10       | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z2_logL_Intercept   |  0.4137 |  0.3101 |  0.3869 |  0.0260 |  1.1376 |
| sp_k1_student.11       | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z2_logA_Intercept   |  0.5465 |  0.4007 |  0.5338 |  0.0309 |  1.5379 |
| sp_k1_student.12       | sp_k1_student       | sp_k1           | student   | sd_h1\_\_z2_logk_Intercept   |  0.3925 |  0.2817 |  0.4225 |  0.0232 |  1.0867 |
| sp_k1_student.13       | sp_k1_student       | sp_k1           | student   | sigma_z1                     |  1.3153 |  1.3150 |  0.0280 |  1.2699 |  1.3624 |
| sp_k1_student.14       | sp_k1_student       | sp_k1           | student   | sigma_z2                     |  1.1010 |  1.1004 |  0.0353 |  1.0441 |  1.1606 |
| sp_k1_student.15       | sp_k1_student       | sp_k1           | student   | nu_z1                        | 62.7812 | 59.9859 | 19.5252 | 35.9584 | 99.0548 |
| sp_k1_student.16       | sp_k1_student       | sp_k1           | student   | nu_z2                        | 51.5792 | 49.0681 | 17.7236 | 27.4972 | 84.7478 |

# Posterior Predictive Checks

The next table compares observed values with posterior predictive draws
on the original scale.

``` r
pp_tbl <- wf$step3_summary
knitr::kable(
  pp_tbl[, c(
    "model", "model_structure", "family", "response", "n", "observed_mean", "yrep_median_mean",
    "mean_residual", "rmse", "mean_abs_std_residual",
    "prop_outside_yrep90", "prop_pp_twotail_lt_0.10", "prop_pp_twotail_lt_0.05"
  )],
  digits = 4,
  align = c("l", "l", "l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r")
)
```

|                             | model               | model_structure | family    | response |    n | observed_mean | yrep_median_mean | mean_residual |   rmse | mean_abs_std_residual | prop_outside_yrep90 | prop_pp_twotail_lt_0.10 | prop_pp_twotail_lt_0.05 |
|:----------------------------|:--------------------|:----------------|:----------|:---------|-----:|--------------:|-----------------:|--------------:|-------:|----------------------:|--------------------:|------------------------:|------------------------:|
| base_k0_gamma.befa.st       | base_k0_gamma       | base_k0         | gamma     | befa.st  | 1132 |        1.4411 |           1.3780 |        0.0631 | 0.3230 |                0.7550 |              0.0866 |                  0.0848 |                  0.0512 |
| ftp_k0_gamma.befa.st        | ftp_k0_gamma        | ftp_k0          | gamma     | befa.st  | 1132 |        1.4411 |           1.3805 |        0.0605 | 0.3195 |                0.7590 |              0.0901 |                  0.0875 |                  0.0530 |
| ftp_k0_lognormal.befa.st    | ftp_k0_lognormal    | ftp_k0          | lognormal | befa.st  | 1132 |        1.4411 |           2.0052 |       -0.5642 | 0.6694 |                0.1291 |              0.0733 |                  0.0680 |                  0.0203 |
| ftp_k0_student.befa.st      | ftp_k0_student      | ftp_k0          | student   | befa.st  | 1132 |        1.4411 |           2.0053 |       -0.5643 | 0.6683 |                0.1151 |              0.0716 |                  0.0707 |                  0.0133 |
| ftp_k1_gamma.befa.st        | ftp_k1_gamma        | ftp_k1          | gamma     | befa.st  | 1132 |        1.4411 |           1.3804 |        0.0607 | 0.3188 |                0.7551 |              0.0883 |                  0.0875 |                  0.0504 |
| ftp_k1_lognormal.befa.st    | ftp_k1_lognormal    | ftp_k1          | lognormal | befa.st  | 1132 |        1.4411 |           2.0028 |       -0.5617 | 0.6643 |                0.1277 |              0.0769 |                  0.0760 |                  0.0168 |
| ftp_k1_student.befa.st      | ftp_k1_student      | ftp_k1          | student   | befa.st  | 1132 |        1.4411 |           2.0058 |       -0.5647 | 0.6668 |                0.1152 |              0.0716 |                  0.0707 |                  0.0168 |
| ftp_sp_k0_gamma.befa.st     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | befa.st  | 1132 |        1.4411 |           1.3958 |        0.0453 | 0.2640 |                0.7201 |              0.0963 |                  0.0954 |                  0.0495 |
| ftp_sp_k0_lognormal.befa.st | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | befa.st  | 1132 |        1.4411 |           2.0050 |       -0.5639 | 0.6677 |                0.1286 |              0.0777 |                  0.0751 |                  0.0186 |
| ftp_sp_k0_student.befa.st   | ftp_sp_k0_student   | ftp_sp_k0       | student   | befa.st  | 1132 |        1.4411 |           2.0034 |       -0.5623 | 0.6658 |                0.1148 |              0.0671 |                  0.0645 |                  0.0168 |
| ftp_sp_k1_gamma.befa.st     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | befa.st  | 1132 |        1.4411 |           1.3953 |        0.0457 | 0.2640 |                0.7172 |              0.0945 |                  0.0945 |                  0.0486 |
| ftp_sp_k1_lognormal.befa.st | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | befa.st  | 1132 |        1.4411 |           2.0041 |       -0.5630 | 0.6642 |                0.1277 |              0.0716 |                  0.0663 |                  0.0168 |
| ftp_sp_k1_student.befa.st   | ftp_sp_k1_student   | ftp_sp_k1       | student   | befa.st  | 1132 |        1.4411 |           2.0044 |       -0.5633 | 0.6659 |                0.1155 |              0.0707 |                  0.0698 |                  0.0141 |
| ftp_sp_k2_gamma.befa.st     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | befa.st  | 1132 |        1.4411 |           1.3951 |        0.0460 | 0.2635 |                0.7121 |              0.0936 |                  0.0928 |                  0.0451 |
| ftp_sp_k2_lognormal.befa.st | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | befa.st  | 1132 |        1.4411 |           2.0044 |       -0.5633 | 0.6674 |                0.1272 |              0.0733 |                  0.0716 |                  0.0194 |
| ftp_sp_k2_student.befa.st   | ftp_sp_k2_student   | ftp_sp_k2       | student   | befa.st  | 1132 |        1.4411 |           2.0043 |       -0.5632 | 0.6659 |                0.1170 |              0.0742 |                  0.0716 |                  0.0150 |
| PFT_k0_student.befa.st      | PFT_k0_student      | PFT_k0          | student   | befa.st  | 1132 |        1.4411 |           2.0075 |       -0.5664 | 0.6667 |                0.1148 |              0.0742 |                  0.0707 |                  0.0168 |
| PFT_k1_student.befa.st      | PFT_k1_student      | PFT_k1          | student   | befa.st  | 1132 |        1.4411 |           2.0064 |       -0.5654 | 0.6670 |                0.1155 |              0.0680 |                  0.0636 |                  0.0177 |
| PFT_sp_k0_student.befa.st   | PFT_sp_k0_student   | PFT_sp_k0       | student   | befa.st  | 1132 |        1.4411 |           2.0084 |       -0.5673 | 0.6674 |                0.1136 |              0.0698 |                  0.0689 |                  0.0168 |
| PFT_sp_k1_student.befa.st   | PFT_sp_k1_student   | PFT_sp_k1       | student   | befa.st  | 1132 |        1.4411 |           2.0053 |       -0.5643 | 0.6663 |                0.1163 |              0.0724 |                  0.0716 |                  0.0159 |
| PFT_sp_k2_student.befa.st   | PFT_sp_k2_student   | PFT_sp_k2       | student   | befa.st  | 1132 |        1.4411 |           2.0070 |       -0.5659 | 0.6679 |                0.1164 |              0.0751 |                  0.0716 |                  0.0177 |
| sp_k0_gamma.befa.st         | sp_k0_gamma         | sp_k0           | gamma     | befa.st  | 1132 |        1.4411 |           1.3955 |        0.0455 | 0.2640 |                0.7175 |              0.0989 |                  0.0989 |                  0.0504 |
| sp_k0_lognormal.befa.st     | sp_k0_lognormal     | sp_k0           | lognormal | befa.st  | 1132 |        1.4411 |           2.0137 |       -0.5726 | 0.6691 |                0.1272 |              0.0777 |                  0.0733 |                  0.0168 |
| sp_k0_student.befa.st       | sp_k0_student       | sp_k0           | student   | befa.st  | 1132 |        1.4411 |           2.0159 |       -0.5748 | 0.6703 |                0.1138 |              0.0689 |                  0.0671 |                  0.0150 |
| sp_k1_gamma.befa.st         | sp_k1_gamma         | sp_k1           | gamma     | befa.st  | 1132 |        1.4411 |           1.3946 |        0.0464 | 0.2629 |                0.7113 |              0.0981 |                  0.0963 |                  0.0468 |
| sp_k1_lognormal.befa.st     | sp_k1_lognormal     | sp_k1           | lognormal | befa.st  | 1132 |        1.4411 |           2.0144 |       -0.5733 | 0.6729 |                0.1262 |              0.0760 |                  0.0742 |                  0.0177 |
| sp_k1_student.befa.st       | sp_k1_student       | sp_k1           | student   | befa.st  | 1132 |        1.4411 |           2.0155 |       -0.5744 | 0.6724 |                0.1151 |              0.0733 |                  0.0707 |                  0.0159 |
| base_k0_gamma.befr.st       | base_k0_gamma       | base_k0         | gamma     | befr.st  |  524 |        0.4422 |           0.4038 |        0.0385 | 0.2696 |                0.7418 |              0.0916 |                  0.0897 |                  0.0534 |
| ftp_k0_gamma.befr.st        | ftp_k0_gamma        | ftp_k0          | gamma     | befr.st  |  524 |        0.4422 |           0.4070 |        0.0352 | 0.2594 |                0.7462 |              0.0878 |                  0.0878 |                  0.0515 |
| ftp_k0_lognormal.befr.st    | ftp_k0_lognormal    | ftp_k0          | lognormal | befr.st  |  524 |        0.4422 |           1.0084 |       -0.5662 | 0.6294 |                0.2166 |              0.0515 |                  0.0496 |                  0.0191 |
| ftp_k0_student.befr.st      | ftp_k0_student      | ftp_k0          | student   | befr.st  |  524 |        0.4422 |           1.0039 |       -0.5617 | 0.6242 |                0.1911 |              0.0420 |                  0.0401 |                  0.0191 |
| ftp_k1_gamma.befr.st        | ftp_k1_gamma        | ftp_k1          | gamma     | befr.st  |  524 |        0.4422 |           0.4074 |        0.0349 | 0.2587 |                0.7477 |              0.0878 |                  0.0878 |                  0.0515 |
| ftp_k1_lognormal.befr.st    | ftp_k1_lognormal    | ftp_k1          | lognormal | befr.st  |  524 |        0.4422 |           1.0049 |       -0.5627 | 0.6264 |                0.2160 |              0.0515 |                  0.0458 |                  0.0210 |
| ftp_k1_student.befr.st      | ftp_k1_student      | ftp_k1          | student   | befr.st  |  524 |        0.4422 |           1.0057 |       -0.5635 | 0.6278 |                0.1923 |              0.0496 |                  0.0477 |                  0.0191 |
| ftp_sp_k0_gamma.befr.st     | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     | befr.st  |  524 |        0.4422 |           0.4130 |        0.0292 | 0.2262 |                0.7073 |              0.0706 |                  0.0706 |                  0.0324 |
| ftp_sp_k0_lognormal.befr.st | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | befr.st  |  524 |        0.4422 |           1.0078 |       -0.5656 | 0.6294 |                0.2172 |              0.0515 |                  0.0515 |                  0.0229 |
| ftp_sp_k0_student.befr.st   | ftp_sp_k0_student   | ftp_sp_k0       | student   | befr.st  |  524 |        0.4422 |           1.0062 |       -0.5639 | 0.6279 |                0.1949 |              0.0515 |                  0.0515 |                  0.0229 |
| ftp_sp_k1_gamma.befr.st     | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     | befr.st  |  524 |        0.4422 |           0.4128 |        0.0294 | 0.2218 |                0.7033 |              0.0687 |                  0.0687 |                  0.0363 |
| ftp_sp_k1_lognormal.befr.st | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | befr.st  |  524 |        0.4422 |           1.0062 |       -0.5640 | 0.6258 |                0.2149 |              0.0534 |                  0.0458 |                  0.0210 |
| ftp_sp_k1_student.befr.st   | ftp_sp_k1_student   | ftp_sp_k1       | student   | befr.st  |  524 |        0.4422 |           1.0057 |       -0.5635 | 0.6262 |                0.1956 |              0.0534 |                  0.0515 |                  0.0210 |
| ftp_sp_k2_gamma.befr.st     | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     | befr.st  |  524 |        0.4422 |           0.4130 |        0.0293 | 0.2242 |                0.7035 |              0.0725 |                  0.0687 |                  0.0363 |
| ftp_sp_k2_lognormal.befr.st | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | befr.st  |  524 |        0.4422 |           1.0041 |       -0.5619 | 0.6248 |                0.2181 |              0.0496 |                  0.0496 |                  0.0229 |
| ftp_sp_k2_student.befr.st   | ftp_sp_k2_student   | ftp_sp_k2       | student   | befr.st  |  524 |        0.4422 |           1.0033 |       -0.5611 | 0.6243 |                0.1954 |              0.0458 |                  0.0420 |                  0.0191 |
| PFT_k0_student.befr.st      | PFT_k0_student      | PFT_k0          | student   | befr.st  |  524 |        0.4422 |           1.0083 |       -0.5661 | 0.6299 |                0.1937 |              0.0477 |                  0.0439 |                  0.0229 |
| PFT_k1_student.befr.st      | PFT_k1_student      | PFT_k1          | student   | befr.st  |  524 |        0.4422 |           1.0080 |       -0.5658 | 0.6283 |                0.1910 |              0.0496 |                  0.0477 |                  0.0210 |
| PFT_sp_k0_student.befr.st   | PFT_sp_k0_student   | PFT_sp_k0       | student   | befr.st  |  524 |        0.4422 |           1.0079 |       -0.5657 | 0.6273 |                0.1931 |              0.0496 |                  0.0496 |                  0.0191 |
| PFT_sp_k1_student.befr.st   | PFT_sp_k1_student   | PFT_sp_k1       | student   | befr.st  |  524 |        0.4422 |           1.0058 |       -0.5636 | 0.6269 |                0.1956 |              0.0496 |                  0.0477 |                  0.0229 |
| PFT_sp_k2_student.befr.st   | PFT_sp_k2_student   | PFT_sp_k2       | student   | befr.st  |  524 |        0.4422 |           1.0062 |       -0.5640 | 0.6268 |                0.1923 |              0.0458 |                  0.0420 |                  0.0172 |
| sp_k0_gamma.befr.st         | sp_k0_gamma         | sp_k0           | gamma     | befr.st  |  524 |        0.4422 |           0.4134 |        0.0289 | 0.2249 |                0.7081 |              0.0706 |                  0.0687 |                  0.0344 |
| sp_k0_lognormal.befr.st     | sp_k0_lognormal     | sp_k0           | lognormal | befr.st  |  524 |        0.4422 |           1.0151 |       -0.5729 | 0.6337 |                0.2140 |              0.0477 |                  0.0477 |                  0.0191 |
| sp_k0_student.befr.st       | sp_k0_student       | sp_k0           | student   | befr.st  |  524 |        0.4422 |           1.0193 |       -0.5771 | 0.6376 |                0.1923 |              0.0458 |                  0.0458 |                  0.0191 |
| sp_k1_gamma.befr.st         | sp_k1_gamma         | sp_k1           | gamma     | befr.st  |  524 |        0.4422 |           0.4129 |        0.0293 | 0.2268 |                0.7126 |              0.0706 |                  0.0706 |                  0.0382 |
| sp_k1_lognormal.befr.st     | sp_k1_lognormal     | sp_k1           | lognormal | befr.st  |  524 |        0.4422 |           1.0161 |       -0.5739 | 0.6357 |                0.2099 |              0.0477 |                  0.0458 |                  0.0210 |
| sp_k1_student.befr.st       | sp_k1_student       | sp_k1           | student   | befr.st  |  524 |        0.4422 |           1.0200 |       -0.5778 | 0.6385 |                0.1883 |              0.0439 |                  0.0382 |                  0.0191 |

# Model Comparison

Gamma and lognormal fits use the positive response scale directly.
Student fits saved with `z1`/`z2` use log-transformed responses, so
their pointwise log-likelihoods are Jacobian-adjusted before LOO
comparison.

``` r
loo_tbl <- wf$step5_total
knitr::kable(
  loo_tbl[, c(
    "model", "model_structure", "family", "elpd_loo", "p_loo", "looic", "max_pareto_k",
    "n_pareto_k_gt_0.7", "n_pareto_k_gt_1.0", "elpd_diff_from_best",
    "looic_diff_from_best", "elpd_diff_from_structure_best",
    "looic_diff_from_structure_best"
  )],
  digits = 3,
  align = c("l", "l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r")
)
```

|     | model               | model_structure | family    | elpd_loo |  p_loo |    looic | max_pareto_k | n_pareto_k_gt_0.7 | n_pareto_k_gt_1.0 | elpd_diff_from_best | looic_diff_from_best | elpd_diff_from_structure_best | looic_diff_from_structure_best |
|:----|:--------------------|:----------------|:----------|---------:|-------:|---------:|-------------:|------------------:|------------------:|--------------------:|---------------------:|------------------------------:|-------------------------------:|
| 4   | ftp_sp_k0_gamma     | ftp_sp_k0       | gamma     |  426.691 | 61.854 | -853.382 |    0.5978466 |                 0 |                 0 |               0.000 |                0.000 |                         0.000 |                          0.000 |
| 8   | sp_k1_gamma         | sp_k1           | gamma     |  426.657 | 66.001 | -853.315 |    0.7654125 |                 1 |                 0 |              -0.033 |                0.067 |                         0.000 |                          0.000 |
| 7   | sp_k0_gamma         | sp_k0           | gamma     |  426.418 | 62.896 | -852.836 |    0.7079919 |                 1 |                 0 |              -0.273 |                0.546 |                         0.000 |                          0.000 |
| 5   | ftp_sp_k1_gamma     | ftp_sp_k1       | gamma     |  425.886 | 67.521 | -851.771 |    0.6509906 |                 0 |                 0 |              -0.805 |                1.611 |                         0.000 |                          0.000 |
| 6   | ftp_sp_k2_gamma     | ftp_sp_k2       | gamma     |  425.367 | 71.679 | -850.734 |    0.6093178 |                 0 |                 0 |              -1.324 |                2.647 |                         0.000 |                          0.000 |
| 2   | ftp_k0_gamma        | ftp_k0          | gamma     |  149.178 | 11.897 | -298.356 |    0.3357957 |                 0 |                 0 |            -277.513 |              555.025 |                         0.000 |                          0.000 |
| 3   | ftp_k1_gamma        | ftp_k1          | gamma     |  149.125 | 12.807 | -298.249 |    0.3417675 |                 0 |                 0 |            -277.566 |              555.132 |                         0.000 |                          0.000 |
| 1   | base_k0_gamma       | base_k0         | gamma     |  103.961 |  8.824 | -207.923 |    0.2246412 |                 0 |                 0 |            -322.729 |              645.459 |                         0.000 |                          0.000 |
| 12  | ftp_sp_k1_lognormal | ftp_sp_k1       | lognormal | -982.121 |  1.159 | 1964.241 |    0.9936590 |                 3 |                 0 |           -1408.811 |             2817.623 |                     -1408.006 |                       2816.012 |
| 13  | ftp_sp_k2_lognormal | ftp_sp_k2       | lognormal | -982.233 |  1.145 | 1964.465 |    0.9431435 |                 8 |                 0 |           -1408.923 |             2817.847 |                     -1407.600 |                       2815.200 |
| 10  | ftp_k1_lognormal    | ftp_k1          | lognormal | -982.294 |  1.133 | 1964.588 |    0.8137229 |                 3 |                 0 |           -1408.985 |             2817.970 |                     -1131.419 |                       2262.838 |
| 11  | ftp_sp_k0_lognormal | ftp_sp_k0       | lognormal | -983.129 |  1.141 | 1966.257 |    0.6948225 |                 0 |                 0 |           -1409.819 |             2819.639 |                     -1409.819 |                       2819.639 |
| 9   | ftp_k0_lognormal    | ftp_k0          | lognormal | -983.238 |  1.082 | 1966.476 |    0.4873743 |                 0 |                 0 |           -1409.929 |             2819.857 |                     -1132.416 |                       2264.832 |
| 19  | ftp_sp_k1_student   | ftp_sp_k1       | student   | -989.010 |  1.150 | 1978.020 |    0.8674983 |                 4 |                 0 |           -1415.701 |             2831.402 |                     -1414.896 |                       2829.791 |
| 17  | ftp_k1_student      | ftp_k1          | student   | -989.072 |  1.139 | 1978.144 |    0.7376755 |                 1 |                 0 |           -1415.763 |             2831.525 |                     -1138.196 |                       2276.393 |
| 20  | ftp_sp_k2_student   | ftp_sp_k2       | student   | -989.091 |  1.135 | 1978.183 |    0.9562807 |                 3 |                 0 |           -1415.782 |             2831.564 |                     -1414.458 |                       2828.917 |
| 18  | ftp_sp_k0_student   | ftp_sp_k0       | student   | -989.930 |  1.114 | 1979.861 |    0.7048694 |                 1 |                 0 |           -1416.621 |             2833.242 |                     -1416.621 |                       2833.242 |
| 16  | ftp_k0_student      | ftp_k0          | student   | -990.089 |  1.093 | 1980.178 |    0.7176416 |                 1 |                 0 |           -1416.780 |             2833.559 |                     -1139.267 |                       2278.534 |
| 25  | PFT_sp_k2_student   | PFT_sp_k2       | student   | -990.725 |  1.243 | 1981.451 |    1.0347743 |                 3 |                 1 |           -1417.416 |             2834.832 |                         0.000 |                          0.000 |
| 24  | PFT_sp_k1_student   | PFT_sp_k1       | student   | -990.792 |  1.264 | 1981.584 |    0.9953054 |                 4 |                 0 |           -1417.483 |             2834.966 |                         0.000 |                          0.000 |
| 14  | sp_k0_lognormal     | sp_k0           | lognormal | -990.877 |  1.191 | 1981.753 |    0.8333682 |                 3 |                 0 |           -1417.567 |             2835.135 |                     -1417.294 |                       2834.589 |
| 22  | PFT_k1_student      | PFT_k1          | student   | -990.960 |  1.237 | 1981.920 |    0.9262235 |                 3 |                 0 |           -1417.651 |             2835.302 |                         0.000 |                          0.000 |
| 15  | sp_k1_lognormal     | sp_k1           | lognormal | -991.208 |  1.213 | 1982.415 |    0.7735805 |                 3 |                 0 |           -1417.898 |             2835.797 |                     -1417.865 |                       2835.730 |
| 23  | PFT_sp_k0_student   | PFT_sp_k0       | student   | -991.720 |  1.210 | 1983.439 |    0.7878087 |                 3 |                 0 |           -1418.411 |             2836.821 |                         0.000 |                          0.000 |
| 21  | PFT_k0_student      | PFT_k0          | student   | -991.859 |  1.193 | 1983.718 |    0.6712409 |                 0 |                 0 |           -1418.550 |             2837.099 |                         0.000 |                          0.000 |
| 26  | sp_k0_student       | sp_k0           | student   | -999.712 |  1.220 | 1999.424 |    0.8889120 |                 3 |                 0 |           -1426.403 |             2852.805 |                     -1426.130 |                       2852.259 |
| 27  | sp_k1_student       | sp_k1           | student   | -999.930 |  1.225 | 1999.859 |    0.9260131 |                 3 |                 0 |           -1426.620 |             2853.241 |                     -1426.587 |                       2853.174 |

# Outputs

``` r
out_files <- sort(list.files(wf$out_dir))
out_tbl <- data.frame(
  file = out_files,
  path = file.path(wf$out_dir, out_files),
  stringsAsFactors = FALSE
)

knitr::kable(
  out_tbl,
  col.names = c("File", "Path"),
  escape = FALSE,
  align = c("l", "l")
)
```

| File | Path |
|:-----|:-----|

# Reproducibility

``` r
cat("**R version:**", R.version.string, "\n\n")
```

    ## **R version:** R version 4.6.0 (2026-04-24)

``` r
cat("**Loaded packages:** `brms`, `posterior`, `loo`\n\n")
```

    ## **Loaded packages:** `brms`, `posterior`, `loo`

``` r
cat("**Project root:**", project_root, "\n\n")
```

    ## **Project root:** /Users/hyli0001/wrd/b/Dynamic_allometrics

``` r
cat("**Workflow script:**", file.path(script_dir, "02_developing_BEF_Allometrics_Bayes.R"), "\n\n")
```

    ## **Workflow script:** /Users/hyli0001/wrd/b/Dynamic_allometrics/scripts/02_developing_BEF_Allometrics_Bayes.R

``` r
cat("```text\n")
```

``` text

``` r
cat(paste(capture.output(sessionInfo()), collapse = "\n"))
```

R version 4.6.0 (2026-04-24) Platform: aarch64-apple-darwin23 Running
under: macOS Tahoe 26.5.1

Matrix products: default BLAS:
/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib
LAPACK:
/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;
LAPACK version 3.12.1

locale: \[1\] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Europe/Stockholm tzcode source: internal

attached base packages: \[1\] stats graphics grDevices utils datasets
methods base

other attached packages: \[1\] loo_2.9.0 posterior_1.7.0 brms_2.23.0
Rcpp_1.1.1-1.1

loaded via a namespace (and not attached): \[1\] bridgesampling_1.2-1
tensorA_0.36.2.1 generics_0.1.4  
\[4\] stringi_1.8.7 lattice_0.22-9 digest_0.6.39  
\[7\] magrittr_2.0.5 evaluate_1.0.5 grid_4.6.0  
\[10\] RColorBrewer_1.1-3 mvtnorm_1.4-0 fastmap_1.2.0  
\[13\] plyr_1.8.9 Matrix_1.7-5 processx_3.9.0  
\[16\] pkgbuild_1.4.8 backports_1.5.1 gridExtra_2.3  
\[19\] Brobdingnag_1.2-9 QuickJSR_1.10.0 scales_1.4.0  
\[22\] codetools_0.2-20 abind_1.4-8 cli_3.6.6  
\[25\] rlang_1.2.0 cmdstanr_0.8.0 yaml_2.3.12  
\[28\] StanHeaders_2.32.10 inline_0.3.21 rstan_2.32.7  
\[31\] tools_4.6.0 parallel_4.6.0 reshape2_1.4.5  
\[34\] rstantools_2.6.0 checkmate_2.3.4 coda_0.19-4.1  
\[37\] dplyr_1.2.1 ggplot2_4.0.3 vctrs_0.7.3  
\[40\] R6_2.6.1 stats4_4.6.0 matrixStats_1.5.0  
\[43\] lifecycle_1.0.5 stringr_1.6.0 pkgconfig_2.0.3  
\[46\] RcppParallel_5.1.11-2 pillar_1.11.1 gtable_0.3.6  
\[49\] glue_1.8.1 xfun_0.59 tibble_3.3.1  
\[52\] tidyselect_1.2.1 knitr_1.51 farver_2.1.2  
\[55\] bayesplot_1.15.0 htmltools_0.5.9 nlme_3.1-169  
\[58\] rmarkdown_2.31 compiler_4.6.0 S7_0.2.2  
\[61\] distributional_0.7.0

``` r
cat("\n```\n")
```


    <details>
    <summary>Workflow script</summary>

    ```r

    # ---- setup ----
    setwd("/Users/hyli0001/wrd/b/Dynamic_allometrics/")
    project_root <- getwd()
    library(brms)
    library(posterior)
    library(loo)

    run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
    run_generated_at <- as.character(Sys.time())
    out_dir <- file.path("bayes_outputs", "bef_bayes", run_id)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # ---- model-list ----
    model_dir <- "bayes_outputs"
    model_file_pattern <- "^xp_(.+)_(gamma|lognormal|stud|student)_rsd_4-4k-99-15[.]rds$"

    family_label <- function(family) {
      family <- tolower(family)
      ifelse(family %in% c("stud", "student"), "student", family)
    }

    model_structure_label <- function(structure_key) {
      out <- gsub("-", "_", structure_key)
      out <- gsub("sp_code", "sp", out)
      out[out == "none_k0"] <- "base_k0"
      out
    }

    discover_model_info <- function(model_dir) {
      files <- list.files(model_dir, pattern = "[.]rds$", full.names = TRUE)
      file_names <- basename(files)
      match_data <- regmatches(file_names, regexec(model_file_pattern, file_names))
      keep <- lengths(match_data) > 0

      if (!any(keep)) {
        stop("No Bayesian model .rds files found in ", model_dir)
      }

      match_data <- match_data[keep]
      structure_key <- vapply(match_data, `[[`, character(1), 2)
      family <- family_label(vapply(match_data, `[[`, character(1), 3))
      model_structure <- model_structure_label(structure_key)
      k_depth <- sub("^.*_(k[0-9]+)$", "\\1", structure_key)
      hierarchy <- sub("_[kK][0-9]+$", "", structure_key)
      hierarchy[hierarchy == "none"] <- "base"

      out <- data.frame(
        model = paste(model_structure, family, sep = "_"),
        model_structure = model_structure,
        family = family,
        hierarchy = hierarchy,
        k_depth = k_depth,
        model_file = files[keep],
        stringsAsFactors = FALSE
      )

      out[order(out$model_structure, out$family), ]
    }

    model_info <- discover_model_info(model_dir)

    missing_files <- model_info$model_file[!file.exists(model_info$model_file)]
    if (length(missing_files) > 0) {
      stop("Missing model files:\n", paste(missing_files, collapse = "\n"))
    }

    model_info$model_time <- as.character(file.info(model_info$model_file)$mtime)
    model_info$model_file_size <- file.info(model_info$model_file)$size

    fits <- setNames(lapply(model_info$model_file, readRDS), model_info$model)

    # ---- helper-functions ----
    write_txt <- function(x, file) {
      write.table(
        x,
        file.path(out_dir, file),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }

    add_model <- function(df, model) {
      i <- match(model, model_info$model)
      data.frame(
        model = model,
        model_structure = model_info$model_structure[i],
        family = model_info$family[i],
        df,
        row.names = NULL
      )
    }

    response_info_for_fit <- function(fit) {
      data_names <- names(fit$data)

      if (all(c("y1m1", "y2") %in% data_names)) {
        return(data.frame(
          resp = c("y1m1", "y2"),
          response = c("befa.st", "befr.st"),
          transform = c("shift_y1m1", "identity"),
          jacobian_log_response = c(FALSE, FALSE),
          stringsAsFactors = FALSE
        ))
      }

      if (all(c("z1", "z2") %in% data_names)) {
        return(data.frame(
          resp = c("z1", "z2"),
          response = c("befa.st", "befr.st"),
          transform = c("log_y1m1", "log_y2"),
          jacobian_log_response = c(TRUE, TRUE),
          stringsAsFactors = FALSE
        ))
      }

      stop("Unknown response columns in fit$data: ", paste(data_names, collapse = ", "))
    }

    resp_label <- function(resp_info) {
      resp_info$response[[1]]
    }

    backtransform <- function(resp_info, x) {
      switch(
        resp_info$transform[[1]],
        shift_y1m1 = x + 1,
        identity = x,
        log_y1m1 = exp(x) + 1,
        log_y2 = exp(x),
        stop("Unknown response transform: ", resp_info$transform[[1]])
      )
    }

    safe_max <- function(x) {
      if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }

    safe_min <- function(x) {
      if (length(x) == 0 || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    }

    safe_which_max <- function(x) {
      if (length(x) == 0 || all(is.na(x))) NA_integer_ else which.max(replace(x, is.na(x), -Inf))
    }

    safe_which_min <- function(x) {
      if (length(x) == 0 || all(is.na(x))) NA_integer_ else which.min(replace(x, is.na(x), Inf))
    }

    run_by_model <- function(fun) {
      do.call(rbind, Map(fun, fits, model_info$model))
    }

    run_by_model_response <- function(fun) {
      out <- list()
      k <- 1L

      for (i in seq_along(fits)) {
        resp_info <- response_info_for_fit(fits[[i]])
        for (j in seq_len(nrow(resp_info))) {
          out[[k]] <- fun(fits[[i]], model_info$model[[i]], resp_info[j, , drop = FALSE])
          k <- k + 1L
        }
      }

      do.call(rbind, out)
    }

    plot_grid <- function(n) {
      ncol <- min(4, n)
      nrow <- ceiling(n / ncol)
      c(nrow = nrow, ncol = ncol)
    }

    chain_id_for_fit <- function(fit, ndraws) {
      nchains <- NULL

      if (inherits(fit$fit, "stanfit")) {
        nchains <- fit$fit@sim$chains
      } else if (is.function(fit$fit$num_chains)) {
        nchains <- fit$fit$num_chains()
      }

      if (is.null(nchains) || ndraws %% nchains != 0) {
        return(NULL)
      }

      rep(seq_len(nchains), each = ndraws / nchains)
    }

    compute_loo <- function(fit, resp_info, newdata) {
      resp <- resp_info$resp[[1]]
      log_lik_matrix <- log_lik(fit, resp = resp, newdata = newdata)

      if (anyNA(log_lik_matrix)) {
        stop("NA values found in log-likelihood for response ", resp)
      }

      if (isTRUE(resp_info$jacobian_log_response[[1]])) {
        log_lik_matrix <- sweep(log_lik_matrix, 2, newdata[[resp]], "-")
      }

      chain_id <- chain_id_for_fit(fit, nrow(log_lik_matrix))
      r_eff <- if (is.null(chain_id)) {
        NULL
      } else {
        relative_eff(exp(log_lik_matrix), chain_id = chain_id)
      }

      loo(log_lik_matrix, r_eff = r_eff)
    }

    # ---- manifest ----
    manifest <- data.frame(
      run_id = run_id,
      run_generated_at = run_generated_at,
      model_info,
      n_rows = sapply(fits, function(x) nrow(x$data)),
      stringsAsFactors = FALSE
    )

    write_txt(manifest, "00_model_manifest.txt")
    manifest

    # ---- step-1-mcmc ----
    get_mcmc <- function(fit, model) {
      np <- nuts_params(fit)
      s <- summarise_draws(as_draws_df(fit))
      s <- s[!grepl("^Ymi_", s$variable), ]

      i_rhat <- safe_which_max(s$rhat)
      i_bulk <- safe_which_min(s$ess_bulk)
      i_tail <- safe_which_min(s$ess_tail)
      treedepth <- np$Value[np$Parameter == "treedepth__"]

      out <- data.frame(
        divergences = sum(np$Parameter == "divergent__" & np$Value == 1),
        max_treedepth = safe_max(treedepth),
        max_rhat = safe_max(s$rhat),
        max_rhat_parameter = if (is.na(i_rhat)) NA_character_ else s$variable[i_rhat],
        min_bulk_ess = safe_min(s$ess_bulk),
        min_bulk_ess_parameter = if (is.na(i_bulk)) NA_character_ else s$variable[i_bulk],
        min_tail_ess = safe_min(s$ess_tail),
        min_tail_ess_parameter = if (is.na(i_tail)) NA_character_ else s$variable[i_tail],
        stringsAsFactors = FALSE
      )

      add_model(out, model)
    }

    step1 <- run_by_model(get_mcmc)
    write_txt(step1, "01_mcmc_diagnostics.txt")
    step1

    # ---- step-2-posterior ----
    get_posteriors <- function(fit, model) {
      s <- summarise_draws(as_draws_df(fit))
      s <- s[grepl("^(b_|sd_|shape_|sigma_|nu_)", s$variable), ]

      keep <- intersect(
        c("variable", "mean", "median", "sd", "mad", "q5", "q95"),
        names(s)
      )

      add_model(s[, keep, drop = FALSE], model)
    }

    step2 <- run_by_model(get_posteriors)
    write_txt(step2, "02_posterior_parameter_summary.txt")
    head(step2)

    # ---- step-3-ppcheck-data ----
    get_ppc_yrep <- function(fit, model, resp_info, ndraws = 1000) {
      resp <- resp_info$resp[[1]]
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      obs <- backtransform(resp_info, dat[[resp]])

      yrep <- posterior_predict(
        fit,
        newdata = dat,
        resp = resp,
        ndraws = ndraws
      )
      yrep <- backtransform(resp_info, yrep)

      yrep_median <- apply(yrep, 2, median)
      yrep_q05 <- apply(yrep, 2, quantile, 0.05)
      yrep_q95 <- apply(yrep, 2, quantile, 0.95)
      yrep_sd <- apply(yrep, 2, sd)

      pp_p_upper <- colMeans(sweep(yrep, 2, obs, `>=`))
      pp_p_two_tail <- 2 * pmin(pp_p_upper, 1 - pp_p_upper)

      out <- data.frame(
        response = resp_label(resp_info),
        row_id = seq_len(nrow(dat)),
        x = dat$x,
        observed = obs,
        yrep_median = yrep_median,
        yrep_q05 = yrep_q05,
        yrep_q95 = yrep_q95,
        residual_obs_minus_yrep_median = obs - yrep_median,
        std_residual = (obs - yrep_median) / yrep_sd,
        pp_p_upper = pp_p_upper,
        pp_p_two_tail = pp_p_two_tail,
        outside_yrep90 = obs < yrep_q05 | obs > yrep_q95,
        stringsAsFactors = FALSE
      )

      add_model(out, model)
    }

    step3_yrep <- run_by_model_response(get_ppc_yrep)
    write_txt(step3_yrep, "03_ppcheck_observed_yrep_pointwise.txt")
    head(step3_yrep)

    # ---- step-3-ppcheck-summary ----
    summarize_ppcheck <- function(d) {
      data.frame(
        model = d$model[1],
        model_structure = d$model_structure[1],
        family = d$family[1],
        response = d$response[1],
        n = nrow(d),
        observed_mean = mean(d$observed, na.rm = TRUE),
        yrep_median_mean = mean(d$yrep_median, na.rm = TRUE),
        mean_residual = mean(d$residual_obs_minus_yrep_median, na.rm = TRUE),
        rmse = sqrt(mean(d$residual_obs_minus_yrep_median^2, na.rm = TRUE)),
        mean_abs_std_residual = mean(abs(d$std_residual), na.rm = TRUE),
        prop_outside_yrep90 = mean(d$outside_yrep90, na.rm = TRUE),
        prop_pp_twotail_lt_0.10 = mean(d$pp_p_two_tail < 0.10, na.rm = TRUE),
        prop_pp_twotail_lt_0.05 = mean(d$pp_p_two_tail < 0.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }

    step3_summary <- do.call(
      rbind,
      lapply(split(step3_yrep, list(step3_yrep$model, step3_yrep$response), drop = TRUE), summarize_ppcheck)
    )

    step3_summary <- step3_summary[order(step3_summary$response, step3_summary$model), ]
    write_txt(step3_summary, "03_ppcheck_observed_yrep_summary.txt")
    step3_summary

    # ---- step-3-ppcheck-density-plot ----
    plot_ppc_density_one <- function(fit, model, resp_info, ndraws = 50) {
      resp <- resp_info$resp[[1]]
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      obs <- backtransform(resp_info, dat[[resp]])

      yrep <- posterior_predict(
        fit,
        newdata = dat,
        resp = resp,
        ndraws = ndraws
      )
      yrep <- backtransform(resp_info, yrep)

      obs_den <- density(obs, na.rm = TRUE)
      yrep_den <- lapply(seq_len(nrow(yrep)), function(i) density(yrep[i, ], na.rm = TRUE))
      y_lim <- range(c(obs_den$y, unlist(lapply(yrep_den, `[[`, "y"))), na.rm = TRUE)

      plot(
        obs_den,
        lwd = 2,
        col = "black",
        ylim = y_lim,
        xlab = resp_label(resp_info),
        main = paste(resp_label(resp_info), model, sep = " | ")
      )

      for (i in seq_along(yrep_den)) {
        lines(yrep_den[[i]], col = rgb(1, 0, 0, 0.15))
      }

      lines(obs_den, lwd = 2, col = "black")
      legend(
        "topright",
        legend = c("Observed", "yrep"),
        col = c("black", "red"),
        lwd = 2,
        bty = "n"
      )
    }

    pdf(file.path(out_dir, "03_ppcheck_density_observed_yrep.pdf"), width = 12, height = 7)
    old_par <- par(no.readonly = TRUE)
    layout_dim <- plot_grid(length(fits))

    for (response_name in c("befa.st", "befr.st")) {
      par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (i in seq_along(fits)) {
        resp_info <- response_info_for_fit(fits[[i]])
        resp_info <- resp_info[resp_info$response == response_name, , drop = FALSE]
        plot_ppc_density_one(fits[[i]], model_info$model[[i]], resp_info)
      }

      mtext(paste("Density PPC:", response_name), outer = TRUE, cex = 1.2)
    }

    par(old_par)
    dev.off()

    # ---- step-4-observed-predicted-data ----
    get_obs_pred <- function(fit, model, resp_info) {
      resp <- resp_info$resp[[1]]
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      pred <- fitted(fit, newdata = dat, resp = resp)[, "Estimate"]

      out <- data.frame(
        response = resp_label(resp_info),
        row_id = seq_len(nrow(dat)),
        x = dat$x,
        observed = backtransform(resp_info, dat[[resp]]),
        predicted = backtransform(resp_info, pred),
        stringsAsFactors = FALSE
      )

      out$residual <- out$observed - out$predicted
      add_model(out, model)
    }

    step4 <- run_by_model_response(get_obs_pred)
    write_txt(step4, "04_observed_vs_predicted.txt")
    head(step4)

    # ---- step-4-observed-predicted-plot ----
    pdf(file.path(out_dir, "04_observed_vs_predicted.pdf"), width = 12, height = 7)
    old_par <- par(no.readonly = TRUE)
    layout_dim <- plot_grid(length(fits))

    for (resp in c("befa.st", "befr.st")) {
      par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (mod in names(fits)) {
        d <- step4[step4$response == resp & step4$model == mod, ]
        lim <- range(c(d$observed, d$predicted), na.rm = TRUE)

        plot(
          d$observed,
          d$predicted,
          xlim = lim,
          ylim = lim,
          pch = 16,
          col = rgb(0, 0, 0, 0.35),
          xlab = "Observed",
          ylab = "Predicted",
          main = paste(resp, mod, sep = " | ")
        )

        abline(0, 1, col = "red", lwd = 2)
      }

      mtext(paste("Observed vs predicted:", resp), outer = TRUE, cex = 1.2)
    }

    par(old_par)
    dev.off()

    # ---- step-5-loo ----
    get_loo <- function(fit, model, resp_info) {
      resp <- resp_info$resp[[1]]
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      lo <- compute_loo(fit, resp_info = resp_info, newdata = dat)
      est <- lo$estimates
      pk <- pareto_k_values(lo)

      out <- data.frame(
        response = resp_label(resp_info),
        n = nrow(dat),
        elpd_loo = est["elpd_loo", "Estimate"],
        elpd_loo_se = est["elpd_loo", "SE"],
        p_loo = est["p_loo", "Estimate"],
        looic = est["looic", "Estimate"],
        looic_se = est["looic", "SE"],
        max_pareto_k = max(pk, na.rm = TRUE),
        n_pareto_k_gt_0.7 = sum(pk > 0.7, na.rm = TRUE),
        n_pareto_k_gt_1.0 = sum(pk > 1.0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )

      add_model(out, model)
    }

    step5 <- run_by_model_response(get_loo)
    step5 <- step5[order(step5$response, -step5$elpd_loo), ]
    write_txt(step5, "05_loo_metrics_by_response.txt")
    step5

    # ---- step-5-total-ranking ----
    step5_total <- aggregate(
      cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~
        model + model_structure + family,
      data = step5,
      FUN = sum
    )

    step5_total$max_pareto_k <- tapply(step5$max_pareto_k, step5$model, max, na.rm = TRUE)[step5_total$model]
    step5_total$elpd_diff_from_best <- step5_total$elpd_loo - max(step5_total$elpd_loo)
    step5_total$looic_diff_from_best <- step5_total$looic - min(step5_total$looic)
    step5_total$elpd_diff_from_structure_best <- ave(
      step5_total$elpd_loo,
      step5_total$model_structure,
      FUN = function(x) x - max(x)
    )
    step5_total$looic_diff_from_structure_best <- ave(
      step5_total$looic,
      step5_total$model_structure,
      FUN = function(x) x - min(x)
    )
    step5_total <- step5_total[order(-step5_total$elpd_loo), ]

    write_txt(step5_total, "05_loo_total_ranking.txt")
    step5_total

    # ---- step-5-response-rsd-data ----
    get_response_curve <- function(fit, model, resp_info, n_grid = 100) {
      resp <- resp_info$resp[[1]]
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      x_grid <- seq(min(dat$x, na.rm = TRUE), max(dat$x, na.rm = TRUE), length.out = n_grid)
      newdat <- data.frame(x = x_grid)

      if ("h1" %in% names(fit$data)) {
        newdat$h1 <- fit$data$h1[1]
      }
      if ("h2" %in% names(fit$data)) {
        newdat$h2 <- fit$data$h2[1]
      }

      pred <- fitted(
        fit,
        newdata = newdat,
        resp = resp,
        re_formula = NA,
        probs = c(0.05, 0.95)
      )

      obs <- add_model(data.frame(
        response = resp_label(resp_info),
        row_id = seq_len(nrow(dat)),
        x = dat$x,
        observed = backtransform(resp_info, dat[[resp]]),
        stringsAsFactors = FALSE
      ), model)

      curve <- add_model(data.frame(
        response = resp_label(resp_info),
        x = x_grid,
        predicted = backtransform(resp_info, pred[, "Estimate"]),
        pred_q05 = backtransform(resp_info, pred[, "Q5"]),
        pred_q95 = backtransform(resp_info, pred[, "Q95"]),
        stringsAsFactors = FALSE
      ), model)

      list(obs = obs, curve = curve)
    }

    curve_list <- unlist(
      Map(function(fit, model) {
        resp_info <- response_info_for_fit(fit)
        lapply(seq_len(nrow(resp_info)), function(i) {
          get_response_curve(fit, model, resp_info[i, , drop = FALSE])
        })
      }, fits, model_info$model),
      recursive = FALSE
    )

    step5_obs <- do.call(rbind, lapply(curve_list, `[[`, "obs"))
    step5_curve <- do.call(rbind, lapply(curve_list, `[[`, "curve"))

    write_txt(step5_obs, "05_response_vs_rsd_observed.txt")
    write_txt(step5_curve, "05_response_vs_rsd_predicted.txt")
    head(step5_curve)

    # ---- step-5-response-rsd-plot ----
    pdf(file.path(out_dir, "05_response_vs_rsd.pdf"), width = 12, height = 7)
    old_par <- par(no.readonly = TRUE)
    layout_dim <- plot_grid(length(fits))

    for (resp in c("befa.st", "befr.st")) {
      par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (mod in names(fits)) {
        obs_d <- step5_obs[step5_obs$response == resp & step5_obs$model == mod, ]
        cur_d <- step5_curve[step5_curve$response == resp & step5_curve$model == mod, ]
        ord <- order(cur_d$x)
        ylim <- range(c(obs_d$observed, cur_d$pred_q05, cur_d$pred_q95), na.rm = TRUE)

        plot(
          obs_d$x,
          obs_d$observed,
          pch = 16,
          col = rgb(0, 0, 0, 0.35),
          xlab = "Relative stand density (rsd)",
          ylab = resp,
          main = paste(resp, mod, sep = " | "),
          ylim = ylim
        )

        lines(cur_d$x[ord], cur_d$predicted[ord], col = "red", lwd = 2)
        lines(cur_d$x[ord], cur_d$pred_q05[ord], col = "red", lty = 2)
        lines(cur_d$x[ord], cur_d$pred_q95[ord], col = "red", lty = 2)
      }

      mtext(paste("Response vs rsd:", resp), outer = TRUE, cex = 1.2)
    }

    par(old_par)
    dev.off()

    # ---- run-summary ----
    writeLines(
      c(
        paste("Run ID:", run_id),
        paste("Run generated at:", run_generated_at),
        paste("Project root:", project_root),
        paste("Output directory:", out_dir),
        "",
        "Files:",
        list.files(out_dir)
      ),
      file.path(out_dir, "RUN_SUMMARY.txt")
    )

    out_dir

</details>
