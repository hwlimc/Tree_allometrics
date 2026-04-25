## Dynamic Grouping Message

The central question of this study is not only whether stand structure improves biomass expansion factors and allometric equations, but also whether those dynamic responses are similar enough across species to justify merging equations at broader biological levels. The updated analysis therefore evaluates dynamic relationships at `species`, `genus`, `family`, `PFT3`, and `PFT5` levels using the grouping information imported from `data/PFT.txt`.

The main message is that stand density matters, but it does not fully explain variation in biomass expansion factors and biomass allocation among plots, and the appropriate equation level depends on how consistently stand-structure effects are shared within a grouping. Relative stand density (`rsd`) remains a useful baseline crowding variable, but in many cases an additional descriptor of within-stand diameter distribution improves model performance. That means the study should be framed as a search for the coarsest grouping level at which dynamic responses to stand structure can still be represented without unacceptable loss of fit.

In practical terms, two stands with similar mean tree size and density can still differ in branch, foliage, or crown biomass if their size distributions differ. Likewise, two species within the same genus or PFT may still differ in which stand descriptor explains their BEF or biomass variation best. The correct grouping level therefore has to be tested explicitly rather than assumed.

The modeling implication is not to replace `rsd`, but to extend the baseline equation and then compare that extension across grouping levels. A defensible workflow is:

1. For allometric biomass components, first compare `log(y) ~ log(d)` and `log(y) ~ log(d^2 h)`.
2. Use the better size-only form as the base equation for that response and grouping.
3. Test one additional non-collinear stand descriptor such as `rsd`, `std`, `d.qm`, `d.sd`, `d.cv`, `d.iqr`, `d.p90_p10`, `d.skew`, or `d.gini`.
4. Compare the resulting dynamic responses at `species`, `genus`, `family`, `PFT3`, and `PFT5` levels.

Height-variation terms such as `h.sd` or `h.cv` can be presented as optional extensions. They should not be the default inputs for the main equations, but readers who already have height information can use them for potentially better projections in some groups. If no ecologically sensible stand term shows a stable relationship, then the appropriate conclusion is that no dynamic stand correction is needed there and a mean BEF or simpler equation should be used.

The present screening suggests that stem biomass often prefers `log(d^2 h)`, whereas branch and foliage biomass usually still prefer `log(d)`. The stand effect should therefore be interpreted as an additional correction to the best size-only allometric form, not as a substitute for choosing the right base equation.

This leads to the manuscript-level conclusion:

`Dynamic biomass expansion factors and allometric relationships are influenced by both stand crowding and within-stand size structure, but the strength and form of these effects are not necessarily identical among species. Therefore, the appropriate level for equation generalization should be tested explicitly from species to genus, family, and plant functional type.`

## Figures

Figure 1 shows correlation among the core stand-level variables and helps identify terms that should not be included together in the same final model.

[Figure 1. Correlation heatmap of candidate stand-level variables](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_1_stand_term_correlation_heatmap.pdf>)

Figure 2 shows the strongest improvements in model fit relative to the `dbh`-only baseline for the core operational set. This supports the argument that stand structure improves explanation, but that the best stand term differs across grouping levels and responses.

[Figure 2. Top stand terms ranked by improvement over the dbh-only model](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_2_best_stand_terms_delta_aic.pdf>)

Figures 3-7 show examples of dynamic stand-term relationships after pooling at broader grouping levels using the core variable set. They illustrate that some merged groups retain a strong and coherent stand-structure signal with simple operational predictors. Figure 8 shows an optional height-variance extension for readers who have those data available.

[Figure 3. Fagaceae family Bf vs d.qm](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_3_family_Fagaceae_Bf_d_qm.pdf>)

[Figure 4. Quercus genus Bf vs d.qm](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_4_genus_Quercus_Bf_d_qm.pdf>)

[Figure 5. Fagaceae family befa.st vs d.qm](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_5_family_Fagaceae_befa_st_d_qm.pdf>)

[Figure 6. Quercus genus befa.st vs d.qm](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_6_genus_Quercus_befa_st_d_qm.pdf>)

[Figure 7. PFT5 ENF befa.v vs rsd](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_7_pft5_ENF_befa_v_rsd.pdf>)

[Figure 8. Optional PFT3 DBF Bf vs h.sd](</Users/hyli0001/OneDrive/wrd/b/Dynamic_allometrics/manuscript/figs/Figure_8_optional_pft3_DBF_Bf_h_sd.pdf>)

## Suggested Short Version

If a shorter version is needed for the manuscript text:

`Allometric equations should first be built on the better size-only form, usually `log(d)` or `log(d^2 h)` depending on the biomass component. Stand crowding and diameter-distribution terms can then explain additional variation and may remain stable after merging to genus, family, or plant functional type. Height-variation terms may improve projection when available, but they should be treated as optional extensions rather than default inputs.`
