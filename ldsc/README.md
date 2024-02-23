### SNP-based heritability and genetic correlations were estimated for each phenotype (combination) using Linkage Disequilibrium Score Regression (LDSC; Bulik-Sullivan et al., 2015).
> For installation and step-by-step user guides, please refer to the repository of the software on https://github.com/bulik/ldsc.

In short:
1. Munge your summary statistics with `munge.slurm.sh` to transform them into the standard LDSC format.
2. Use the munged summary statistics to run `h2.slurm.sh`, this will give you the SNP-based heritability estimates. Plot the results with `plot_ldsc_h2_results.R`
3. Use the munged summary statistics to run `rg.slurm.sh` and obtain all pairwise genetic correlations between the traits. Plot the results with `plot_ldsc_rg_results.ipynb`
