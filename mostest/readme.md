1. Create z-matrices with original and permuted GWAS z-scores across all phenotypes with ```mostest_multimodal.ipynb```. 

      - Two NxM matrices are created: one for original genotypes and another for permuted genotypes, where N is the number of variants tested in the GWASs and M is the number of phenotypes. 
      - Each matrix is saved as three separate files with specified common prefix: ```<prefix>.cols``` file - a single column with variant IDs, ```<prefix>.rows``` file - a single column with phenotype names, ```<prefix>.dat``` - binary file with NxM matrix of z scores (in float32 format). 
      - NB! the order of variants and phenotypes should be the same in the original and corresponding permuted files. ```<prefix>.dat``` files are then used in the ```mostest_multimodal.m``` to estimate mostest test statistics and corresponding p-values.

2. Run MOSTest with ```submit_mostest.slurm```. The main input of MOSTest are
      - univariate GWAS summary statistics for each phenotype based on original genotypes (see plink folder)
      - univariate GWAS summary statistics for each phenotype based on permuted genotypes (see plink folder)
      - z-matrices created in 1.

3. With ```submit_pvals2csv.slurm```:
   
    - Convert matlab output to CSV file

3. With ```submit_clump.slurm```:
   
    - Create "FUMA-style" genome-wide significant loci, using https://github.com/precimed/python_convert

4. With ```Manhattan_plots.R```:

    - Create Manhattan plot for each multivariate GWAS output
