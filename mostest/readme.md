1. Create z-matrices with original and permuted GWAS z-scores across all phenotypes

2. Run MOSTest with ```submit_mostest.slurm```. The main input of MOSTest are
      - univariate GWAS summary statistics for each phenotype based on original genotypes (see plink folder)
      - univariate GWAS summary statistics for each phenotype based on permuted genotypes (see plink folder)
      - z-matrices created in 1.

3. With ```submit_pvals2csv.slurm```:
   
    - Create "FUMA-style" genome-wide significant loci, using https://github.com/precimed/python_convert

3. With ```submit_clump.slurm```:
   
    - Create "FUMA-style" genome-wide significant loci, using https://github.com/precimed/python_convert

4. With ```Manhattan_plots.R```:

    - Create Manhattan plot for each multivariate GWAS output
