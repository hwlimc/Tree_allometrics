# BEF Bayesian Workflow

# 1 Overview

This report is generated from the workflow script [02_developing_BEF_Allometrics_Bayes.R](02_developing_BEF_Allometrics_Bayes.R).  
Edit the `.R` file only. This Markdown file is a GitHub-readable record of that run.

``` r
cat("**Run ID:**", wf$run_id, "\n\n")
```

    ## **Run ID:** 20260708_163138

``` r
cat("**Run time:**", wf$run_generated_at, "\n\n")
```

    ## **Run time:** 2026-07-08 16:31:38.275159

``` r
cat("**Project root:**", project_root, "\n\n")
```

    ## **Project root:** /Users/hyli0001/wrd/b/Dynamic_allometrics

``` r
cat("**Output directory:**", wf$out_dir, "\n\n")
```

    ## **Output directory:** bayes_outputs/bef_bayes/20260708_163138

# 2 Inputs

``` r
input_manifest <- wf$manifest
input_manifest$model_time <- format(as.POSIXct(input_manifest$model_time), "%Y-%m-%d %H:%M:%S")

knitr::kable(
  input_manifest[, c("model", "model_file", "model_time", "n_rows")],
  col.names = c("Model", "File", "File time", "Rows"),
  align = c("l", "l", "l", "r")
)
```

|           | Model     | File                                                     | File time           | Rows |
|:----------|:----------|:---------------------------------------------------------|:--------------------|-----:|
| base_k0   | base_k0   | bayes_outputs/xp_none_k0_gamma_rsd_4-4k-99-15.rds        | 2026-06-18 11:43:18 | 1132 |
| ftp_k0    | ftp_k0    | bayes_outputs/xp_ftp_k0_gamma_rsd_4-4k-99-15.rds         | 2026-06-18 13:38:13 | 1132 |
| ftp_k1    | ftp_k1    | bayes_outputs/xp_ftp_k1_gamma_rsd_4-4k-99-15.rds         | 2026-06-18 14:34:29 | 1132 |
| sp_k0     | sp_k0     | bayes_outputs/xp_sp_code_k0_gamma_rsd_4-4k-99-15.rds     | 2026-06-18 12:04:21 | 1132 |
| sp_k1     | sp_k1     | bayes_outputs/xp_sp_code_k1_gamma_rsd_4-4k-99-15.rds     | 2026-06-18 12:56:34 | 1132 |
| ftp_sp_k0 | ftp_sp_k0 | bayes_outputs/xp_ftp-sp_code_k0_gamma_rsd_4-4k-99-15.rds | 2026-06-18 14:23:41 | 1132 |
| ftp_sp_k1 | ftp_sp_k1 | bayes_outputs/xp_ftp-sp_code_k1_gamma_rsd_4-4k-99-15.rds | 2026-06-18 14:31:54 | 1132 |
| ftp_sp_k2 | ftp_sp_k2 | bayes_outputs/xp_ftp-sp_code_k2_gamma_rsd_4-4k-99-15.rds | 2026-06-18 14:52:57 | 1132 |

# 3 Diagnostics

``` r
diag_tbl <- wf$step1[, c(
  "model", "divergences", "max_treedepth", "max_rhat",
  "max_rhat_parameter", "min_bulk_ess", "min_bulk_ess_parameter",
  "min_tail_ess", "min_tail_ess_parameter"
)]

knitr::kable(
  diag_tbl,
  digits = 3,
  align = c("l", "r", "r", "r", "l", "r", "l", "r", "l")
)
```

|           | model     | divergences | max_treedepth | max_rhat | max_rhat_parameter                     | min_bulk_ess | min_bulk_ess_parameter              | min_tail_ess | min_tail_ess_parameter                  |
|:----------|:----------|------------:|--------------:|---------:|:---------------------------------------|-------------:|:------------------------------------|-------------:|:----------------------------------------|
| base_k0   | base_k0   |           0 |             9 |    1.002 | lp\_\_                                 |     1583.585 | lp\_\_                              |     3083.463 | lp\_\_                                  |
| ftp_k0    | ftp_k0    |           0 |            10 |    1.002 | lp\_\_                                 |     1587.842 | lp\_\_                              |     2853.027 | sd_h1\_\_y2_logL_Intercept              |
| ftp_k1    | ftp_k1    |           0 |            11 |    1.002 | sd_h1\_\_y2_logk_Intercept             |     1857.799 | lp\_\_                              |     2897.967 | r_h1\_\_y1m1_logL\[mono_B,Intercept\]   |
| sp_k0     | sp_k0     |           0 |             9 |    1.006 | b_y1m1_logL_Intercept                  |      590.006 | b_y2_logk_Intercept                 |      490.818 | r_h1\_\_y2_logL\[CJa,Intercept\]        |
| sp_k1     | sp_k1     |           0 |             9 |    1.004 | sd_h1\_\_y2_logk_Intercept             |      722.826 | b_y2_logk_Intercept                 |     1025.062 | sd_h1\_\_y2_logL_Intercept              |
| ftp_sp_k0 | ftp_sp_k0 |           0 |            12 |    1.005 | z_6\[1,9\]                             |      515.005 | b_y2_logk_Intercept                 |      783.671 | r_h2\_\_y2_logL\[mono_B.PxT,Intercept\] |
| ftp_sp_k1 | ftp_sp_k1 |           0 |            10 |    1.008 | sd_h2\_\_y2_logA_Intercept             |      456.353 | sd_h2\_\_y2_logA_Intercept          |      767.130 | r_h2\_\_y2_logL\[mono_B.PxT,Intercept\] |
| ftp_sp_k2 | ftp_sp_k2 |           0 |            11 |    1.013 | r_h2\_\_y2_logk\[mono_N.CO,Intercept\] |      391.653 | r_h1\_\_y2_logk\[mono_B,Intercept\] |      294.480 | r_h1\_\_y2_logA\[mono_B,Intercept\]     |

# 4 Posterior Summaries

``` r
post_tbl <- wf$step2
post_tbl <- post_tbl[post_tbl$variable %in% c(
  "b_y1m1_logL_Intercept", "b_y1m1_logA_Intercept", "b_y1m1_logk_Intercept",
  "b_y2_logL_Intercept", "b_y2_logA_Intercept", "b_y2_logk_Intercept",
  "sd_h1__y1m1_logL_Intercept", "sd_h1__y1m1_logA_Intercept",
  "sd_h1__y2_logL_Intercept", "sd_h1__y2_logA_Intercept"
), ]

knitr::kable(
  post_tbl[, c("model", "variable", "mean", "median", "sd", "q5", "q95")],
  digits = 4,
  align = c("l", "l", "r", "r", "r", "r", "r")
)
```

|              | model     | variable                     |    mean |  median |     sd |      q5 |     q95 |
|:-------------|:----------|:-----------------------------|--------:|--------:|-------:|--------:|--------:|
| base_k0.1    | base_k0   | b_y1m1_logL_Intercept        | -1.0159 | -1.0152 | 0.0297 | -1.0659 | -0.9685 |
| base_k0.2    | base_k0   | b_y1m1_logA_Intercept        |  0.2666 |  0.2640 | 0.1761 | -0.0195 |  0.5581 |
| base_k0.3    | base_k0   | b_y1m1_logk_Intercept        |  2.1598 |  2.1622 | 0.1158 |  1.9718 |  2.3467 |
| base_k0.4    | base_k0   | b_y2_logL_Intercept          | -1.7253 | -1.5407 | 0.6686 | -3.0499 | -0.9991 |
| base_k0.5    | base_k0   | b_y2_logA_Intercept          | -1.0195 | -1.0011 | 0.2500 | -1.4673 | -0.6627 |
| base_k0.6    | base_k0   | b_y2_logk_Intercept          |  0.0388 | -0.0803 | 0.6282 | -0.8174 |  1.1602 |
| ftp_k0.1     | ftp_k0    | b_y1m1_logL_Intercept        | -1.3748 | -1.2472 | 0.4976 | -2.3717 | -0.7892 |
| ftp_k0.2     | ftp_k0    | b_y1m1_logA_Intercept        | -0.3588 | -0.3048 | 0.5042 | -1.2692 |  0.3724 |
| ftp_k0.3     | ftp_k0    | b_y1m1_logk_Intercept        |  1.9721 |  1.9747 | 0.1326 |  1.7505 |  2.1829 |
| ftp_k0.4     | ftp_k0    | b_y2_logL_Intercept          | -2.0424 | -1.9614 | 0.6742 | -3.2536 | -1.0873 |
| ftp_k0.5     | ftp_k0    | b_y2_logA_Intercept          | -0.9042 | -0.8894 | 0.3300 | -1.4467 | -0.3848 |
| ftp_k0.6     | ftp_k0    | b_y2_logk_Intercept          |  0.0351 | -0.0346 | 0.4386 | -0.5850 |  0.8403 |
| ftp_k0.7     | ftp_k0    | sd_h1\_\_y1m1_logL_Intercept |  0.7005 |  0.5051 | 0.6117 |  0.1526 |  1.9135 |
| ftp_k0.8     | ftp_k0    | sd_h1\_\_y1m1_logA_Intercept |  0.8691 |  0.7086 | 0.6093 |  0.2322 |  2.0612 |
| ftp_k0.9     | ftp_k0    | sd_h1\_\_y2_logL_Intercept   |  0.9390 |  0.7500 | 0.7399 |  0.1817 |  2.3312 |
| ftp_k0.10    | ftp_k0    | sd_h1\_\_y2_logA_Intercept   |  0.4109 |  0.2722 | 0.4511 |  0.0208 |  1.2812 |
| ftp_k1.1     | ftp_k1    | b_y1m1_logL_Intercept        | -1.3692 | -1.2573 | 0.4901 | -2.3444 | -0.7635 |
| ftp_k1.2     | ftp_k1    | b_y1m1_logA_Intercept        | -0.1858 | -0.0866 | 0.4971 | -1.1381 |  0.4527 |
| ftp_k1.3     | ftp_k1    | b_y1m1_logk_Intercept        |  1.3197 |  1.5090 | 0.7485 | -0.1388 |  2.1979 |
| ftp_k1.4     | ftp_k1    | b_y2_logL_Intercept          | -2.1177 | -2.0377 | 0.6945 | -3.3594 | -1.1244 |
| ftp_k1.5     | ftp_k1    | b_y2_logA_Intercept          | -0.8588 | -0.8383 | 0.3147 | -1.3842 | -0.3916 |
| ftp_k1.6     | ftp_k1    | b_y2_logk_Intercept          | -0.0304 | -0.0734 | 0.5316 | -0.8328 |  0.8965 |
| ftp_k1.7     | ftp_k1    | sd_h1\_\_y1m1_logL_Intercept |  0.7043 |  0.5183 | 0.5897 |  0.1649 |  1.8442 |
| ftp_k1.8     | ftp_k1    | sd_h1\_\_y1m1_logA_Intercept |  0.6895 |  0.5231 | 0.6127 |  0.0545 |  1.8693 |
| ftp_k1.10    | ftp_k1    | sd_h1\_\_y2_logL_Intercept   |  0.8229 |  0.6318 | 0.7298 |  0.0737 |  2.2480 |
| ftp_k1.11    | ftp_k1    | sd_h1\_\_y2_logA_Intercept   |  0.3870 |  0.2449 | 0.4341 |  0.0209 |  1.2072 |
| sp_k0.1      | sp_k0     | b_y1m1_logL_Intercept        | -1.1258 | -1.1253 | 0.1254 | -1.3329 | -0.9249 |
| sp_k0.2      | sp_k0     | b_y1m1_logA_Intercept        |  0.1441 |  0.1448 | 0.1756 | -0.1460 |  0.4295 |
| sp_k0.3      | sp_k0     | b_y1m1_logk_Intercept        |  1.9556 |  1.9596 | 0.1111 |  1.7679 |  2.1317 |
| sp_k0.4      | sp_k0     | b_y2_logL_Intercept          | -1.4475 | -1.3019 | 0.4803 | -2.4554 | -0.9664 |
| sp_k0.5      | sp_k0     | b_y2_logA_Intercept          | -1.0139 | -1.0079 | 0.2360 | -1.4100 | -0.6446 |
| sp_k0.6      | sp_k0     | b_y2_logk_Intercept          |  0.6459 |  0.6253 | 0.5696 | -0.2334 |  1.5785 |
| sp_k0.7      | sp_k0     | sd_h1\_\_y1m1_logL_Intercept |  0.5603 |  0.5477 | 0.0972 |  0.4262 |  0.7381 |
| sp_k0.8      | sp_k0     | sd_h1\_\_y1m1_logA_Intercept |  0.3543 |  0.3389 | 0.1130 |  0.1994 |  0.5568 |
| sp_k0.9      | sp_k0     | sd_h1\_\_y2_logL_Intercept   |  0.5255 |  0.4667 | 0.2337 |  0.2770 |  0.9813 |
| sp_k0.10     | sp_k0     | sd_h1\_\_y2_logA_Intercept   |  0.4856 |  0.4346 | 0.3311 |  0.0467 |  1.0986 |
| sp_k1.1      | sp_k1     | b_y1m1_logL_Intercept        | -1.1349 | -1.1336 | 0.1293 | -1.3476 | -0.9247 |
| sp_k1.2      | sp_k1     | b_y1m1_logA_Intercept        |  0.1756 |  0.1801 | 0.1787 | -0.1224 |  0.4613 |
| sp_k1.3      | sp_k1     | b_y1m1_logk_Intercept        |  1.9541 |  1.9594 | 0.1361 |  1.7228 |  2.1684 |
| sp_k1.4      | sp_k1     | b_y2_logL_Intercept          | -1.7717 | -1.6352 | 0.6282 | -3.0060 | -1.0240 |
| sp_k1.5      | sp_k1     | b_y2_logA_Intercept          | -0.9131 | -0.8924 | 0.2265 | -1.3135 | -0.5773 |
| sp_k1.6      | sp_k1     | b_y2_logk_Intercept          |  0.2466 |  0.1348 | 0.6239 | -0.6172 |  1.4040 |
| sp_k1.7      | sp_k1     | sd_h1\_\_y1m1_logL_Intercept |  0.5604 |  0.5495 | 0.0949 |  0.4271 |  0.7275 |
| sp_k1.8      | sp_k1     | sd_h1\_\_y1m1_logA_Intercept |  0.3160 |  0.3042 | 0.1417 |  0.0973 |  0.5630 |
| sp_k1.10     | sp_k1     | sd_h1\_\_y2_logL_Intercept   |  0.5495 |  0.4788 | 0.3186 |  0.1433 |  1.1708 |
| sp_k1.11     | sp_k1     | sd_h1\_\_y2_logA_Intercept   |  0.4102 |  0.3624 | 0.2822 |  0.0472 |  0.9585 |
| ftp_sp_k0.1  | ftp_sp_k0 | b_y1m1_logL_Intercept        | -1.3809 | -1.2695 | 0.4617 | -2.3027 | -0.8380 |
| ftp_sp_k0.2  | ftp_sp_k0 | b_y1m1_logA_Intercept        | -0.1451 | -0.0589 | 0.4468 | -1.0157 |  0.4317 |
| ftp_sp_k0.3  | ftp_sp_k0 | b_y1m1_logk_Intercept        |  1.9281 |  1.9302 | 0.1152 |  1.7346 |  2.1146 |
| ftp_sp_k0.4  | ftp_sp_k0 | b_y2_logL_Intercept          | -1.6949 | -1.5629 | 0.6324 | -2.9276 | -0.9270 |
| ftp_sp_k0.5  | ftp_sp_k0 | b_y2_logA_Intercept          | -0.9962 | -0.9988 | 0.3569 | -1.5869 | -0.4220 |
| ftp_sp_k0.6  | ftp_sp_k0 | b_y2_logk_Intercept          |  0.5798 |  0.5583 | 0.5825 | -0.3200 |  1.5275 |
| ftp_sp_k0.7  | ftp_sp_k0 | sd_h1\_\_y1m1_logL_Intercept |  0.5854 |  0.4047 | 0.5809 |  0.0460 |  1.7371 |
| ftp_sp_k0.9  | ftp_sp_k0 | sd_h1\_\_y1m1_logA_Intercept |  0.5749 |  0.3983 | 0.5800 |  0.0347 |  1.7230 |
| ftp_sp_k0.11 | ftp_sp_k0 | sd_h1\_\_y2_logL_Intercept   |  0.6269 |  0.4409 | 0.6254 |  0.0510 |  1.8257 |
| ftp_sp_k0.13 | ftp_sp_k0 | sd_h1\_\_y2_logA_Intercept   |  0.4166 |  0.2708 | 0.4653 |  0.0224 |  1.2872 |
| ftp_sp_k1.1  | ftp_sp_k1 | b_y1m1_logL_Intercept        | -1.3919 | -1.2812 | 0.4604 | -2.3129 | -0.8563 |
| ftp_sp_k1.2  | ftp_sp_k1 | b_y1m1_logA_Intercept        | -0.1457 | -0.0575 | 0.4574 | -1.0440 |  0.4393 |
| ftp_sp_k1.3  | ftp_sp_k1 | b_y1m1_logk_Intercept        |  1.4397 |  1.6781 | 0.6636 |  0.0510 |  2.1234 |
| ftp_sp_k1.4  | ftp_sp_k1 | b_y2_logL_Intercept          | -1.6396 | -1.5048 | 0.6228 | -2.8540 | -0.8753 |
| ftp_sp_k1.5  | ftp_sp_k1 | b_y2_logA_Intercept          | -0.9445 | -0.9472 | 0.3958 | -1.5713 | -0.3097 |
| ftp_sp_k1.6  | ftp_sp_k1 | b_y2_logk_Intercept          |  0.5501 |  0.5524 | 0.6638 | -0.5205 |  1.6191 |
| ftp_sp_k1.7  | ftp_sp_k1 | sd_h1\_\_y1m1_logL_Intercept |  0.5884 |  0.4104 | 0.5735 |  0.0386 |  1.6915 |
| ftp_sp_k1.9  | ftp_sp_k1 | sd_h1\_\_y1m1_logA_Intercept |  0.5884 |  0.4120 | 0.5841 |  0.0351 |  1.7360 |
| ftp_sp_k1.12 | ftp_sp_k1 | sd_h1\_\_y2_logL_Intercept   |  0.6518 |  0.4674 | 0.6210 |  0.0492 |  1.9026 |
| ftp_sp_k1.14 | ftp_sp_k1 | sd_h1\_\_y2_logA_Intercept   |  0.4939 |  0.3375 | 0.5227 |  0.0288 |  1.4737 |
| ftp_sp_k2.1  | ftp_sp_k2 | b_y1m1_logL_Intercept        | -1.4041 | -1.2965 | 0.4567 | -2.3166 | -0.8695 |
| ftp_sp_k2.2  | ftp_sp_k2 | b_y1m1_logA_Intercept        | -0.1786 | -0.0872 | 0.4868 | -1.1201 |  0.4537 |
| ftp_sp_k2.3  | ftp_sp_k2 | b_y1m1_logk_Intercept        |  1.3768 |  1.5866 | 0.6819 |  0.0043 |  2.1335 |
| ftp_sp_k2.4  | ftp_sp_k2 | b_y2_logL_Intercept          | -1.8945 | -1.8068 | 0.6797 | -3.1356 | -0.9671 |
| ftp_sp_k2.5  | ftp_sp_k2 | b_y2_logA_Intercept          | -0.8663 | -0.8691 | 0.3905 | -1.4789 | -0.2142 |
| ftp_sp_k2.6  | ftp_sp_k2 | b_y2_logk_Intercept          |  0.2561 |  0.1909 | 0.6819 | -0.7471 |  1.4659 |
| ftp_sp_k2.7  | ftp_sp_k2 | sd_h1\_\_y1m1_logL_Intercept |  0.5812 |  0.4035 | 0.5864 |  0.0436 |  1.7099 |
| ftp_sp_k2.9  | ftp_sp_k2 | sd_h1\_\_y1m1_logA_Intercept |  0.6667 |  0.5043 | 0.6064 |  0.0546 |  1.8471 |
| ftp_sp_k2.13 | ftp_sp_k2 | sd_h1\_\_y2_logL_Intercept   |  0.7281 |  0.5227 | 0.7048 |  0.0568 |  2.0969 |
| ftp_sp_k2.15 | ftp_sp_k2 | sd_h1\_\_y2_logA_Intercept   |  0.5415 |  0.3127 | 0.6419 |  0.0233 |  1.8783 |

# 5 Posterior Predictive Checks

The next table compares observed values with posterior predictive draws on the original scale.

``` r
pp_tbl <- wf$step3_summary
knitr::kable(
  pp_tbl[, c(
    "model", "response", "n", "observed_mean", "yrep_median_mean",
    "mean_residual", "rmse", "mean_abs_std_residual",
    "prop_outside_yrep90", "prop_pp_twotail_lt_0.10", "prop_pp_twotail_lt_0.05"
  )],
  digits = 4,
  align = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r")
)
```

|                   | model     | response |    n | observed_mean | yrep_median_mean | mean_residual |   rmse | mean_abs_std_residual | prop_outside_yrep90 | prop_pp_twotail_lt_0.10 | prop_pp_twotail_lt_0.05 |
|:------------------|:----------|:---------|-----:|--------------:|-----------------:|--------------:|-------:|----------------------:|--------------------:|------------------------:|------------------------:|
| base_k0.befa.st   | base_k0   | befa.st  | 1132 |        1.4411 |           1.3776 |        0.0635 | 0.3229 |                0.7598 |              0.0875 |                  0.0857 |                  0.0504 |
| ftp_k0.befa.st    | ftp_k0    | befa.st  | 1132 |        1.4411 |           1.3805 |        0.0605 | 0.3188 |                0.7558 |              0.0875 |                  0.0839 |                  0.0512 |
| ftp_k1.befa.st    | ftp_k1    | befa.st  | 1132 |        1.4411 |           1.3798 |        0.0612 | 0.3188 |                0.7568 |              0.0901 |                  0.0857 |                  0.0477 |
| ftp_sp_k0.befa.st | ftp_sp_k0 | befa.st  | 1132 |        1.4411 |           1.3958 |        0.0452 | 0.2639 |                0.7144 |              0.0989 |                  0.0963 |                  0.0442 |
| ftp_sp_k1.befa.st | ftp_sp_k1 | befa.st  | 1132 |        1.4411 |           1.3953 |        0.0458 | 0.2637 |                0.7155 |              0.0945 |                  0.0936 |                  0.0477 |
| ftp_sp_k2.befa.st | ftp_sp_k2 | befa.st  | 1132 |        1.4411 |           1.3957 |        0.0454 | 0.2614 |                0.7074 |              0.0936 |                  0.0928 |                  0.0495 |
| sp_k0.befa.st     | sp_k0     | befa.st  | 1132 |        1.4411 |           1.3955 |        0.0456 | 0.2634 |                0.7145 |              0.0945 |                  0.0928 |                  0.0512 |
| sp_k1.befa.st     | sp_k1     | befa.st  | 1132 |        1.4411 |           1.3953 |        0.0458 | 0.2630 |                0.7108 |              0.0919 |                  0.0910 |                  0.0451 |
| base_k0.befr.st   | base_k0   | befr.st  |  524 |        0.4422 |           0.4039 |        0.0383 | 0.2685 |                0.7399 |              0.0840 |                  0.0802 |                  0.0553 |
| ftp_k0.befr.st    | ftp_k0    | befr.st  |  524 |        0.4422 |           0.4071 |        0.0351 | 0.2592 |                0.7434 |              0.0916 |                  0.0878 |                  0.0553 |
| ftp_k1.befr.st    | ftp_k1    | befr.st  |  524 |        0.4422 |           0.4072 |        0.0350 | 0.2594 |                0.7466 |              0.0935 |                  0.0916 |                  0.0477 |
| ftp_sp_k0.befr.st | ftp_sp_k0 | befr.st  |  524 |        0.4422 |           0.4134 |        0.0288 | 0.2259 |                0.7128 |              0.0706 |                  0.0706 |                  0.0382 |
| ftp_sp_k1.befr.st | ftp_sp_k1 | befr.st  |  524 |        0.4422 |           0.4129 |        0.0293 | 0.2226 |                0.7037 |              0.0706 |                  0.0668 |                  0.0382 |
| ftp_sp_k2.befr.st | ftp_sp_k2 | befr.st  |  524 |        0.4422 |           0.4135 |        0.0288 | 0.2242 |                0.7047 |              0.0706 |                  0.0706 |                  0.0344 |
| sp_k0.befr.st     | sp_k0     | befr.st  |  524 |        0.4422 |           0.4130 |        0.0292 | 0.2260 |                0.7130 |              0.0649 |                  0.0649 |                  0.0344 |
| sp_k1.befr.st     | sp_k1     | befr.st  |  524 |        0.4422 |           0.4128 |        0.0294 | 0.2261 |                0.7076 |              0.0706 |                  0.0687 |                  0.0363 |

# 6 Model Comparison

``` r
loo_tbl <- wf$step5_total
knitr::kable(
  loo_tbl[, c(
    "model", "elpd_loo", "p_loo", "looic", "max_pareto_k",
    "n_pareto_k_gt_0.7", "n_pareto_k_gt_1.0", "elpd_diff_from_best",
    "looic_diff_from_best"
  )],
  digits = 3,
  align = c("l", "r", "r", "r", "r", "r", "r", "r", "r")
)
```

|     | model     | elpd_loo |  p_loo |    looic | max_pareto_k | n_pareto_k_gt_0.7 | n_pareto_k_gt_1.0 | elpd_diff_from_best | looic_diff_from_best |
|:----|:----------|---------:|-------:|---------:|-------------:|------------------:|------------------:|--------------------:|---------------------:|
| 4   | ftp_sp_k0 |  426.691 | 61.854 | -853.382 |    0.5978466 |                 0 |                 0 |               0.000 |                0.000 |
| 8   | sp_k1     |  426.657 | 66.001 | -853.315 |    0.7654125 |                 1 |                 0 |              -0.033 |                0.067 |
| 7   | sp_k0     |  426.418 | 62.896 | -852.836 |    0.7079919 |                 1 |                 0 |              -0.273 |                0.546 |
| 5   | ftp_sp_k1 |  425.886 | 67.521 | -851.771 |    0.6509906 |                 0 |                 0 |              -0.805 |                1.611 |
| 6   | ftp_sp_k2 |  425.367 | 71.679 | -850.734 |    0.6093178 |                 0 |                 0 |              -1.324 |                2.647 |
| 2   | ftp_k0    |  149.178 | 11.897 | -298.356 |    0.3357957 |                 0 |                 0 |            -277.513 |              555.025 |
| 3   | ftp_k1    |  149.125 | 12.807 | -298.249 |    0.3417675 |                 0 |                 0 |            -277.566 |              555.132 |
| 1   | base_k0   |  103.961 |  8.824 | -207.923 |    0.2246412 |                 0 |                 0 |            -322.729 |              645.459 |

# 7 Outputs

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
|------|------|

# 8 Reproducibility

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

R version 4.6.0 (2026-04-24) Platform: aarch64-apple-darwin23 Running under: macOS Tahoe 26.5.1

Matrix products: default BLAS: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib; LAPACK version 3.12.1

locale: \[1\] C.UTF-8/UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Europe/Stockholm tzcode source: internal

attached base packages: \[1\] stats graphics grDevices utils datasets methods base

other attached packages: \[1\] loo_2.9.0 posterior_1.7.0 brms_2.23.0 Rcpp_1.1.1-1.1

loaded via a namespace (and not attached): \[1\] bridgesampling_1.2-1 tensorA_0.36.2.1 sass_0.4.10  
\[4\] generics_0.1.4 stringi_1.8.7 lattice_0.22-9  
\[7\] digest_0.6.39 magrittr_2.0.5 evaluate_1.0.5  
\[10\] grid_4.6.0 RColorBrewer_1.1-3 mvtnorm_1.4-0  
\[13\] fastmap_1.2.0 plyr_1.8.9 jsonlite_2.0.0  
\[16\] Matrix_1.7-5 processx_3.9.0 pkgbuild_1.4.8  
\[19\] backports_1.5.1 gridExtra_2.3 Brobdingnag_1.2-9  
\[22\] QuickJSR_1.10.0 scales_1.4.0 codetools_0.2-20  
\[25\] jquerylib_0.1.4 abind_1.4-8 cli_3.6.6  
\[28\] rlang_1.2.0 cmdstanr_0.8.0 cachem_1.1.0  
\[31\] yaml_2.3.12 StanHeaders_2.32.10 inline_0.3.21  
\[34\] rstan_2.32.7 tools_4.6.0 parallel_4.6.0  
\[37\] reshape2_1.4.5 rstantools_2.6.0 checkmate_2.3.4  
\[40\] coda_0.19-4.1 dplyr_1.2.1 ggplot2_4.0.3  
\[43\] vctrs_0.7.3 R6_2.6.1 stats4_4.6.0  
\[46\] matrixStats_1.5.0 lifecycle_1.0.5 stringr_1.6.0  
\[49\] pkgconfig_2.0.3 RcppParallel_5.1.11-2 pillar_1.11.1  
\[52\] bslib_0.11.0 gtable_0.3.6 glue_1.8.1  
\[55\] xfun_0.59 tibble_3.3.1 tidyselect_1.2.1  
\[58\] knitr_1.51 farver_2.1.2 bayesplot_1.15.0  
\[61\] htmltools_0.5.9 nlme_3.1-169 rmarkdown_2.31  
\[64\] compiler_4.6.0 S7_0.2.2 distributional_0.7.0

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
    model_info <- data.frame(
      model = c(
        "base_k0",
        "ftp_k0",
        "ftp_k1",
        "sp_k0",
        "sp_k1",
        "ftp_sp_k0",
        "ftp_sp_k1",
        "ftp_sp_k2"
      ),
      model_file = c(
        "bayes_outputs/xp_none_k0_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_ftp_k0_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_ftp_k1_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_sp_code_k0_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_sp_code_k1_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_ftp-sp_code_k0_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_ftp-sp_code_k1_gamma_rsd_4-4k-99-15.rds",
        "bayes_outputs/xp_ftp-sp_code_k2_gamma_rsd_4-4k-99-15.rds"
      ),
      stringsAsFactors = FALSE
    )

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
      data.frame(model = model, df, row.names = NULL)
    }

    resp_label <- function(resp) {
      if (identical(resp, "y1m1")) "befa.st" else "befr.st"
    }

    backtransform <- function(resp, x) {
      if (identical(resp, "y1m1")) x + 1 else x
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
      do.call(
        rbind,
        unlist(
          Map(function(fit, model) {
            list(
              fun(fit, model, "y1m1"),
              fun(fit, model, "y2")
            )
          }, fits, model_info$model),
          recursive = FALSE
        )
      )
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

    compute_loo <- function(fit, resp, newdata) {
      log_lik_matrix <- log_lik(fit, resp = resp, newdata = newdata)

      if (anyNA(log_lik_matrix)) {
        stop("NA values found in log-likelihood for response ", resp)
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
      s <- s[grepl("^(b_|sd_|shape_)", s$variable), ]

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
    get_ppc_yrep <- function(fit, model, resp, ndraws = 1000) {
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      obs <- backtransform(resp, dat[[resp]])

      yrep <- posterior_predict(
        fit,
        newdata = dat,
        resp = resp,
        ndraws = ndraws
      )
      yrep <- backtransform(resp, yrep)

      yrep_median <- apply(yrep, 2, median)
      yrep_q05 <- apply(yrep, 2, quantile, 0.05)
      yrep_q95 <- apply(yrep, 2, quantile, 0.95)
      yrep_sd <- apply(yrep, 2, sd)

      pp_p_upper <- colMeans(sweep(yrep, 2, obs, `>=`))
      pp_p_two_tail <- 2 * pmin(pp_p_upper, 1 - pp_p_upper)

      out <- data.frame(
        response = resp_label(resp),
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
    plot_ppc_density_one <- function(fit, model, resp, ndraws = 50) {
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      obs <- backtransform(resp, dat[[resp]])

      yrep <- posterior_predict(
        fit,
        newdata = dat,
        resp = resp,
        ndraws = ndraws
      )
      yrep <- backtransform(resp, yrep)

      obs_den <- density(obs, na.rm = TRUE)
      yrep_den <- lapply(seq_len(nrow(yrep)), function(i) density(yrep[i, ], na.rm = TRUE))
      y_lim <- range(c(obs_den$y, unlist(lapply(yrep_den, `[[`, "y"))), na.rm = TRUE)

      plot(
        obs_den,
        lwd = 2,
        col = "black",
        ylim = y_lim,
        xlab = resp_label(resp),
        main = paste(resp_label(resp), model, sep = " | ")
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

    for (resp in c("y1m1", "y2")) {
      par(
        mfrow = c(layout_dim["nrow"], layout_dim["ncol"]),
        mar = c(4, 4, 3, 1),
        oma = c(0, 0, 2, 0)
      )

      for (mod in names(fits)) {
        plot_ppc_density_one(fits[[mod]], mod, resp)
      }

      mtext(paste("Density PPC:", resp_label(resp)), outer = TRUE, cex = 1.2)
    }

    par(old_par)
    dev.off()

    # ---- step-4-observed-predicted-data ----
    get_obs_pred <- function(fit, model, resp) {
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      pred <- fitted(fit, newdata = dat, resp = resp)[, "Estimate"]

      out <- data.frame(
        response = resp_label(resp),
        row_id = seq_len(nrow(dat)),
        x = dat$x,
        observed = backtransform(resp, dat[[resp]]),
        predicted = backtransform(resp, pred),
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
    get_loo <- function(fit, model, resp) {
      dat <- fit$data[!is.na(fit$data[[resp]]), , drop = FALSE]
      lo <- compute_loo(fit, resp = resp, newdata = dat)
      est <- lo$estimates
      pk <- pareto_k_values(lo)

      out <- data.frame(
        response = resp_label(resp),
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
      cbind(elpd_loo, p_loo, looic, n, n_pareto_k_gt_0.7, n_pareto_k_gt_1.0) ~ model,
      data = step5,
      FUN = sum
    )

    step5_total$max_pareto_k <- tapply(step5$max_pareto_k, step5$model, max, na.rm = TRUE)[step5_total$model]
    step5_total$elpd_diff_from_best <- step5_total$elpd_loo - max(step5_total$elpd_loo)
    step5_total$looic_diff_from_best <- step5_total$looic - min(step5_total$looic)
    step5_total <- step5_total[order(-step5_total$elpd_loo), ]

    write_txt(step5_total, "05_loo_total_ranking.txt")
    step5_total

    # ---- step-5-response-rsd-data ----
    get_response_curve <- function(fit, model, resp, n_grid = 100) {
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
        response = resp_label(resp),
        row_id = seq_len(nrow(dat)),
        x = dat$x,
        observed = backtransform(resp, dat[[resp]]),
        stringsAsFactors = FALSE
      ), model)

      curve <- add_model(data.frame(
        response = resp_label(resp),
        x = x_grid,
        predicted = backtransform(resp, pred[, "Estimate"]),
        pred_q05 = backtransform(resp, pred[, "Q5"]),
        pred_q95 = backtransform(resp, pred[, "Q95"]),
        stringsAsFactors = FALSE
      ), model)

      list(obs = obs, curve = curve)
    }

    curve_list <- unlist(
      Map(function(fit, model) {
        list(
          get_response_curve(fit, model, "y1m1"),
          get_response_curve(fit, model, "y2")
        )
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
