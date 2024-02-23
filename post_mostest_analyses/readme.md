After phenotye definition (`MOSTest-multimodal/phenotypes/`), running all univariate GWAS (`MOSTest-multimodal/plink/`), and combining them into multivariate summary statistics (`MOSTest-multimodal/mostest/`), we were interested in the following analyses:

1. Locus overlap between the three MOSTest summary statistics (sMRI, fMRI, dMRI) was defined as physically overlapping genome-wide significant loci. Herefore we used the GenomicRanges R-package, for more information see https://github.com/Bioconductor/GenomicRanges. In short, this compares the chromosome and start and end base pair positions of all loci between any pair of summary statistics as done in `Loci_overlap_GRanges.R`.

2. Also the genes that were found to be genome-wide significant in MAGMA [ran on https://fuma.ctglab.nl] were compared to provide a similar overview. The overlapping patterns were then plotted with the eulerr R-package, see `Gene_overlap.R`.

3. Functional consequences of lead SNPs were provided by ANNOVAR as implemented in FUMA [https://fuma.ctglab.nl] and then compared to Watanabe et al. 2019 (Nature Communications) to test for enrichment in `ANNOVAR_enrichment.R`

4. We used Gene Ontology gene-sets to test for enrichment of genes using hypergeometric testing as implemented in FUMA [https://fuma.ctglab.nl]. The exported significant Gene Ontolgoy terms were visualized as a graph using Cytoscape, EnrichmentMap and AutoAnnotate following the Nature Protocol by Reimand & Isserlin et al. and available in `Cytoscape_prep.R`. For more information on how to make Cytoscape network graphics, see: https://github.com/reimandlab/ActivePathways.

5. We visualized the temporal gene expression pattern across development using data from Kang et al [http://www.humanbraintranscriptome.org] using the code in `GeneExpression_simpleMean.R`.

6. We explored cell-type enrichment of genes by performing Fisher Exact Tests for each cell-type identified by Bhaduri et al. using single-cell RNA sequencing analyses of the fetal brain in `MOSTest_CellTypeEnrichment.R`, of which the results were then plotted in `plot_celltype_enrichment.R`.

7. For Supplementary Figure 9, we used FUMA [https://fuma.ctglab.nl] to test for the enrichment of disorder genes (positionally mapped from the original disorder GWAS loci vs conditional FDR loci) in tissue-specific gene-sets with upregulated expression and visualized the results with `Plot_tissue_specificity.R`.






