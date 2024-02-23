Outline: run univariate GWAS for all phenotypes for both original (non-permuted) genotypes and for randomly permuted genotypes. Genotypes are split into chunks to allow parallel processing, providing considerable speedup when running on the HPC.

1. Produce dataset of unrelated individuals and remove rare (MAF<0.005) SNPs with ```remove_related.sh```. This dataset will be used for all subsequent alayses. 

2. Split list of all SNP IDs into chunks with 10000 SNPs (resulting in 907 chunks) with ```make_snp_chunks.sh```.

3. Use ```make_bfile_chunks.slurm.sh``` to
    * split genotypes into chunks using SNP chunks generated in step 2, i.e. 907 (non-permuted) genotype datasets are produced each containing all individuals and 10000 SNPs (the last chunk contains less SNPs).
    * make chunks of permuted genotypes randomly permuting each non-permuted genotype chunk (```permute_bed.py``` python script is used to generate permuted genotypes).

4. Use ```merge_chunks.sh``` to merge all chunks of permuted genotypes.

5. Run univariate GWAS for the original (non-permuted) genotypes (using ```gwas_chr.slurm.sh``` and genotypes produced in step 1) and randomly permuted genotypes (using ```gwas_perm_chr.slurm.sh``` and genotypes produced at step 4). For computational efficiency GWASs are run per-chromosome.

6. Use ```make_LD_ref.sh``` to make genotype files which are subsequently used as an LD reference for clumping od the MOSTest sumstats. ```LD_ref.txt``` file contains a list of 1000 individuals randomly selected from the dataset obtained in step 1.
