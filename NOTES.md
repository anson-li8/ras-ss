Critical Functions
1. `compute_gwas_weights`
Location: R/gwas.R
Significant Inputs: The `geno` input and the individual phenotypes.
Math: Calculate basic linear regression line
Without individual data? The `geno` input is missing in Route A because we only know the summary statistics, not specific individual genotypes. That makes it impossible for this function to actually calculate the basic linear regression line. Another missing piece is the individual phenotypes. This is essentially the dependent variable, again making it impossible to build the basic individualized data points for a regression.
Replacement: The compute_gwas_weights function is no longer needed because we now have the actual final calculated marginal effect sizes already for each SNP. Wherever the program requires output of this function, it can simply be replaced with the corresponding given summary statistics.