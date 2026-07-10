Critical Functions
1. `compute_gwas_weights`
Location: R/gwas.R
Significant Inputs: The `geno` input and the individual phenotypes.
Math: Calculate basic linear regression line
Without individual data? The `geno` input is missing in Route A because we only know the summary statistics, not specific individual genotypes. That makes it impossible for this function to actually calculate the basic linear regression line. Another missing piece is the individual phenotypes. This is essentially the dependent variable, again making it impossible to build the basic individualized data points for a regression.
Replacement: The compute_gwas_weights function is no longer needed because we now have the actual final calculated marginal effect sizes already for each SNP. Wherever the program requires output of this function, it can simply be replaced with the corresponding given summary statistics.
Here is the rubric for the second function, built exactly using your words:
2. `compute_pgs_matrix`
Location: man/compute_pgs_matrix.Rd and R/gwas.R
Significant Inputs: Requires individual genotype dosages and SNP marginal effect sizes
Math: Executes matrix multiplication of the individual genotype dosages (0, 1, or 2) by the SNP marginal effect sizes.
Without individual data? The specific matrix multiplication fails when approached with summary statistics and renders this function unusable. This is because we do not have the individual genotype dosages (0, 1, or 2) to multiply by the SNP marginal effect sizes.
Replacement: An Linkage Disequilibrium (LD) correlation matrix can be used along with the marginal effect size values to calculate the approximate variance. This matrix represents the population-level probabilities that two SNPs are inherited together. By implementing this method, we no longer need the individual genotype data.