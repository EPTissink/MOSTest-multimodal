Assuming all relevant condFDR analyses are done. Do following steps:

1. With ```multimodal_pleioprs.ipynb```:

    - Extract subjects with relevant phenotypes from NoDa csv (save as pheno.txt file) and basic covariates for these subjects (covar.txt).

2. With ```extract_geno.sh```:

    - For subjects selected at step 1 extract genotypes.

3. With ```qc_geno.sh```:
   
    - Perform basic QC.

4. With ```multimodal_pleioprs.ipynb```:

    - Add FDR column to original sumstats keeping only variants where both original p-value and condFDR value are not NA.
    
    - Allign variant IDs between sumstats and plink genotype bim file, keeping only variants presenting both in sumstats and bim file.

5. With ```clump.sh```:

    - Clump sumstats produced in step 4 (with FDR column and alligned to bim file) first based on the original p-value then based on condFDR value using default PRSice2 clumping parameters.

6. With ```select_snp_subsets.sh```:

    - Select N = 10, 100, 1000, 10000, 50000, 100000, 15000 top independent variants both based on the original p-value and on condFDR-based clumping.

7. With ```run_prsice.sh```:

    - Run PRSice2 for each selected subset of clumped variants without additional clumping, using effect sizes from the original sumstats.

