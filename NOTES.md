Critical Functions
1. `compute_gwas_weights`
Location: R/gwas.R
Significant Inputs: The `geno` input (individual genotype dosage matrix), indiivdual phenotypes, and covariates (like age, sex, etc.).
Math: Linear regression for continuous phenotypes or logistic regression for binary traits. Does this independently for each SNP, adjusting for covariates
Without individual data? The `geno` input is missing in Route A because we only know the summary statistics, not specific individual genotypes. That makes it impossible for this function to actually calculate the basic linear regression line. Another missing piece is the individual phenotypes. This is essentially the dependent variable, again making it impossible to build the basic individualized data points for a regression.
Replacement: The compute_gwas_weights function is no longer needed because we now have the actual final calculated marginal effect sizes already for each SNP. Wherever the program requires output of this function, it can simply be replaced with the corresponding given summary statistics.
Here is the rubric for the second function, built exactly using your words:
2. `compute_pgs_matrix`
Location: R/gwas.R
Significant Inputs: Requires individual genotype dosages and SNP marginal effect sizes
Math: Executes matrix multiplication of the individual genotype dosages (0, 1, or 2) by the SNP marginal effect sizes. This outputs the Localized Polygenic Risk Score (LPRS) for ecah individual
Without individual data? The specific matrix multiplication fails when approached with summary statistics and renders this function unusable. This is because we do not have the individual genotype dosages (0, 1, or 2) to multiply by the SNP marginal effect sizes.
Replacement: A Linkage Disequilibrium (LD) correlation matrix can be used along with the marginal effect size values to calculate the approximate variance. This matrix represents the population-level probabilities that two SNPs are inherited together. By implementing this method, we no longer need the individual genotype data.
3. `screen_forward_max_region`
Location: R/screening.R
Significant Inputs: The Localized Polygenic Risk Scores (LPRS) for every individual and the actual phenotypic values of the individuals, and covariates for each individual
Math: Run joint multiple regression to find the window p-value. Independent variable is the LPRS (w/ covariates) and dependent variable is phenotype.
Without individual data? It is impossible to run this joint regression with only summary statistics. There are two key inputs that we do not explicitly know, with those being the Localized Polygenic Risk Scores (LPRS) for every individual and the actual phenotypic values of the individuals (both are individualistic thus are not given along with summary statistics).
Replacement: We can use population-level information to approximate the same joint p-value without the individual data. We can use a analytical quadratic form equation (e.g. generalized Wald test) to estimate the joint test stat and window p-value with just two key inputs: The independent effect sizes of the SNPs and the LD correlation matrix that represents how the SNPs interact with each other.
4. `ras_scan`
Location: R/pipeline.R
Significant Inputs: Requires the individual data for the 50/50 train/hold-out splits and num_rep loops.
Math: Loops through data multiple times to split individuals 50/50 for training/ testing to average out the results and prevent overfit
Without individual data? We do no longer need the 50/50 train/holdout splits because the weights are already given by the summary statistics, as explained earlier.
Replacement: It doesn't need to loop because it could just run the analytical equations directly through the SNP windows in a single pass to output the final Time Series sequence of Regional Association Scores (RAS) that can be used for Change Point Detection.