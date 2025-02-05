# High-dim-mediation-analysis-with-survival-outcomes
This repository contains programming codes for the numerical studies in the article "Post-selection inference for high-dimensional
mediation analysis with survival outcomes" that is authored by Tzu-Jung Huang, Zhonghua Liu and Ian W. McKeague (https://arxiv.org/abs/2408.06517) and accepted by Scandinavian Journal of Statistics. Below is the layout of this repository.

1. The directory *standalone simulation codes* contains the main script *sim_maxMediSurv.R* that demonstrates how to conduct our proposed Stabilized One-Step Estimator ('SOSE'), along with its extended version and other competing methods (as mentioned in Competing Methods Section of the article) in one single Monte Carlo simulation run,
where our codes use auxiliary functions collected in the files *code_maxMediSurv.R*, *data_management.R* and *sim_data_acquisition.R* in the same directory.
  
2. The directory *real data application* contains the main script *TCGA_analysis.R* that applies Bonferroni One-Step Estimator ('BONF_OSE') and Stabilized One-Step Estimator ('SOSE') and their extended versions with one random ordering of data (that is, num_rdo = 1) to 
the TCGA lung cancer survival data, using auxiliary functions collected in the files *code_maxMediSurv.R*, *TCGA_dataacquisition.R* and *data_management.R* in the same directory. Further values of r depend on how many random orderings will be taken to generate the results for 'SOSE', 
say X-fold stabilized one-step estimator that uses num_rdo = X. In the real data application of this article, X = 100.
