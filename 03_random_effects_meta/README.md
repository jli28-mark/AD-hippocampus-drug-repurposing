Primary analysis: four-study, intercept-only random-effects meta-analysis of log2FC.
Input effects are AD minus Control/NC, paired with the native audited model SE.
Each common gene contributes exactly four estimates. No study DEG filter is applied.
Sampling variance = SE squared; random-effects weight = 1/(SE squared + tau squared).
Tau squared: REML. Hypothesis testing and pooled 95% CIs: Standard Wald Z-test (standard normal distribution).
BH adjustment uses the entire planned common-gene family before effect-size/direction filters.
Failed fits remain in meta_all_genes.csv with NA P/FDR and count toward BH family size.
Meta_FDR_significant identifies statistical evidence for a nonzero average study effect.
Candidate signature: Meta_FDR < 0.05, abs(pooled log2FC) >= 0.58, at least 3/4 matching signs.
A zero study effect does not support either direction.
Single-study statistical significance is not required for directional support.
Study FDR values are preserved from each original differential analysis.
Opposite-direction study FDR evidence is flagged, not a primary candidate exclusion.
signature_no_observed_opposite_study_FDR.csv is an explicitly labelled additional subset.
The 4/4 subset is reported separately and is not required for the primary candidate list.
I2, tau2 and maximum study weight are reported; no cutoff on these is used.
Direction support and study weights do not prove that a result is free of single-study influence.
Forest study intervals use normal approximations from the input SE.
Forest pooled intervals use the model confidence interval.
Study estimates are treated as independent; this script cannot verify donor overlap.
Different tissue sampling and covariate models remain sources of between-study differences.
Documentation: https://wviechtb.github.io/metafor/reference/rma.uni.html

