###############################################################################
##################### Testing for Cell Type Enrichment ########################
###############################################################################

library(tidyverse); library(data.table) 
library(org.Hs.eg.db); library(limma); library('biomaRt')

## read in MOSTest mapped genes

MOSTest <- read.table("crossmodality_genes_bonf.txt",
                      sep = "\t", header = T)
T1 <- read.table("T1_specific_genes_bonf.txt",
                 sep = "\t", header = T)
rsfMRI <- read.table("rsfmri_specific_genes_bonf.txt",
                     sep = "\t", header = T)
diffusion  <- read.table("diffusion_specific_genes_bonf.txt",
                         sep = "\t", header = T)

## ENSEMBL to Gene Symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

MOSTest_symbol <- getBM(filters = "ensembl_gene_id", 
                        attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        values = MOSTest$GENE, mart = mart)
table(duplicated(MOSTest_symbol$hgnc_symbol))
  ## missing gene symbol... manually remove
MOSTest_symbol <- filter(MOSTest_symbol, duplicated(hgnc_symbol) == F)
MOSTest_symbol <- MOSTest_symbol[!MOSTest_symbol$ensembl_gene_id == "ENSG00000100101",]

T1_symbol <- getBM(filters = "ensembl_gene_id", 
                   attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   values = T1$GENE, mart = mart)
table(duplicated(T1_symbol$hgnc_symbol))
  ## missing gene symbol... manually remove
T1_symbol <- filter(T1_symbol, duplicated(hgnc_symbol) == F)
T1_symbol <- T1_symbol[!T1_symbol$ensembl_gene_id == "ENSG00000125695",]

rsfMRI_symbol <- getBM(filters = "ensembl_gene_id", 
                       attributes = c("ensembl_gene_id", "hgnc_symbol"),
                       values = rsfMRI$GENE, mart = mart)
table(duplicated(rsfMRI_symbol$hgnc_symbol))

diffusion_symbol <- getBM(filters = "ensembl_gene_id", 
                          attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          values = diffusion$GENE, mart = mart)
table(duplicated(diffusion_symbol$hgnc_symbol))
diffusion_symbol <- filter(diffusion_symbol, duplicated(hgnc_symbol) == F)
  ## missing gene symbol... manually remove
diffusion_symbol <- diffusion_symbol[!diffusion_symbol$ensembl_gene_id == "ENSG00000173366",]


### Update/harmonize gene symbols
MOSTest_symbol$updated_GeneSymbol <- alias2SymbolTable(MOSTest_symbol$hgnc_symbol)
MOSTest_symbol$updated_GeneSymbol[is.na(MOSTest_symbol$updated_GeneSymbol)] <- 
  as.character(MOSTest_symbol$hgnc_symbol[is.na(MOSTest_symbol$updated_GeneSymbol)])

T1_symbol$updated_GeneSymbol <- alias2SymbolTable(T1_symbol$hgnc_symbol)
T1_symbol$updated_GeneSymbol[is.na(T1_symbol$updated_GeneSymbol)] <- 
  as.character(T1_symbol$hgnc_symbol[is.na(T1_symbol$updated_GeneSymbol)])

rsfMRI_symbol$updated_GeneSymbol <- alias2SymbolTable(rsfMRI_symbol$hgnc_symbol)
rsfMRI_symbol$updated_GeneSymbol[is.na(rsfMRI_symbol$updated_GeneSymbol)] <- 
  as.character(rsfMRI_symbol$hgnc_symbol[is.na(rsfMRI_symbol$updated_GeneSymbol)])

diffusion_symbol$updated_GeneSymbol <- alias2SymbolTable(diffusion_symbol$hgnc_symbol)
diffusion_symbol$updated_GeneSymbol[is.na(diffusion_symbol$updated_GeneSymbol)] <- 
  as.character(diffusion_symbol$hgnc_symbol[is.na(diffusion_symbol$updated_GeneSymbol)])



### Load background gene list
HGNC_AllGenes <- read.csv("HGNC_ALL_Genes_non_alt_loci_set.txt",
                          sep="\t")
HGNC_AllGenes$updated_GeneSymbol <- alias2SymbolTable(HGNC_AllGenes$symbol)
HGNC_AllGenes$updated_GeneSymbol[is.na(HGNC_AllGenes$updated_GeneSymbol)] <- 
  as.character(HGNC_AllGenes$symbol[is.na(HGNC_AllGenes$updated_GeneSymbol)])

table(duplicated(HGNC_AllGenes$updated_GeneSymbol))
HGNC_AllGenes <- HGNC_AllGenes[duplicated(HGNC_AllGenes$updated_GeneSymbol) == F,]


# Read in Cell Type File
Bhaduri_2021 <- read.table("Bhaduri_2021_WholeBrain.txt", 
                           sep="\t", header = T)

## remove cerebellum cell-genes + region specific cell type names + harmonize symbols
Bhaduri_2021 <- Bhaduri_2021[- grep("cerebellum", Bhaduri_2021$cluster),]
Bhaduri_2021$cluster <- gsub("_.*", "", Bhaduri_2021$cluster)

Bhaduri_2021$updated_GeneSymbol <- alias2SymbolTable(Bhaduri_2021$gene)
Bhaduri_2021$updated_GeneSymbol[is.na(Bhaduri_2021$updated_GeneSymbol)] <- 
  as.character(Bhaduri_2021$gene[is.na(Bhaduri_2021$updated_GeneSymbol)])

Bhaduri_HGNC_ALL <- Bhaduri_2021[Bhaduri_2021$updated_GeneSymbol %in%
                                   HGNC_AllGenes$updated_GeneSymbol,]
Bhaduri_HGNC_ALL <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>% 
  mutate(dup_gene = duplicated(updated_GeneSymbol)) %>%
  ungroup()
Bhaduri_HGNC_ALL <- Bhaduri_HGNC_ALL[Bhaduri_HGNC_ALL$dup_gene == FALSE,]


## MOSTest Genes
MOSTest_hyper_Bhaduri_OR <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(
    a= length(intersect(MOSTest_symbol$updated_GeneSymbol, updated_GeneSymbol)),
    b= length(updated_GeneSymbol) - a,
    c= length(MOSTest_symbol$updated_GeneSymbol) - a,
    d= length(HGNC_AllGenes$updated_GeneSymbol) - sum(a, b, c), 
    p.value = fisher.test(matrix(c(a,c,b,d), nrow=2))$p.value,
    OR = fisher.test(matrix(c(a,c,b,d), nrow=2))$estimate,
    OR_95CI_upper = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[2],
    OR_95CI_lower = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[1]
  )
print(MOSTest_hyper_Bhaduri_OR, n=nrow(MOSTest_hyper_Bhaduri_OR))


## T1 Genes
T1_hyper_Bhaduri_OR <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(
    a= length(intersect(T1_symbol$updated_GeneSymbol, updated_GeneSymbol)),
    b= length(updated_GeneSymbol) - a,
    c= length(T1_symbol$updated_GeneSymbol) - a,
    d= length(HGNC_AllGenes$updated_GeneSymbol) - sum(a, b, c), 
    p.value = fisher.test(matrix(c(a,c,b,d), nrow=2))$p.value,
    OR = fisher.test(matrix(c(a,c,b,d), nrow=2))$estimate,
    OR_95CI_upper = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[2],
    OR_95CI_lower = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[1]
  )
print(T1_hyper_Bhaduri_OR, n=nrow(T1_hyper_Bhaduri_OR))


## rsfMRI Genes
rsfMRI_hyper_Bhaduri_OR <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(
    a= length(intersect(rsfMRI_symbol$updated_GeneSymbol, updated_GeneSymbol)),
    b= length(updated_GeneSymbol) - a,
    c= length(rsfMRI_symbol$updated_GeneSymbol) - a,
    d= length(HGNC_AllGenes$updated_GeneSymbol) - sum(a, b, c), 
    p.value = fisher.test(matrix(c(a,c,b,d), nrow=2))$p.value,
    OR = fisher.test(matrix(c(a,c,b,d), nrow=2))$estimate,
    OR_95CI_upper = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[2],
    OR_95CI_lower = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[1]
  )
print(rsfMRI_hyper_Bhaduri_OR, n=nrow(rsfMRI_hyper_Bhaduri_OR))


## diffusion Genes
diffusion_hyper_Bhaduri_OR <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(
    a= length(intersect(diffusion_symbol$updated_GeneSymbol, updated_GeneSymbol)),
    b= length(updated_GeneSymbol) - a,
    c= length(diffusion_symbol$updated_GeneSymbol) - a,
    d= length(HGNC_AllGenes$updated_GeneSymbol) - sum(a, b, c), 
    p.value = fisher.test(matrix(c(a,c,b,d), nrow=2))$p.value,
    OR = fisher.test(matrix(c(a,c,b,d), nrow=2))$estimate,
    OR_95CI_upper = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[2],
    OR_95CI_lower = fisher.test(matrix(c(a,c,b,d), nrow=2))$conf.int[1]
  )
print(diffusion_hyper_Bhaduri_OR, n=nrow(diffusion_hyper_Bhaduri_OR))


########################################################################
## Generating Supplementary Tables
MOSTest_hyper_Bhaduri_ol <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(n_overlap=length(intersect(MOSTest_symbol$updated_GeneSymbol, 
                                       updated_GeneSymbol)),
            Genes=intersect(MOSTest_symbol$updated_GeneSymbol, 
                            updated_GeneSymbol)) %>% 
  summarise(n_overlap=n_overlap,
            Overlapping_Genes=toString(Genes)) %>% data.frame()

MOSTest_hyper_Bhaduri_ol_p <- merge(
  MOSTest_hyper_Bhaduri_OR, MOSTest_hyper_Bhaduri_ol, by="cluster", all = T)

MOSTest_hyper_Bhaduri_ol_p <- MOSTest_hyper_Bhaduri_ol_p[
  duplicated(MOSTest_hyper_Bhaduri_ol_p$cluster) ==F,
]


write.table(MOSTest_hyper_Bhaduri_ol_p, 
            sep="\t", quote = F, row.names = F,
            file="MOSTest_Bhaduri_FisherExact.txt",
)



T1_hyper_Bhaduri_ol <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(n_overlap=length(intersect(T1_symbol$updated_GeneSymbol, 
                                       updated_GeneSymbol)),
            Genes=intersect(T1_symbol$updated_GeneSymbol, 
                            updated_GeneSymbol)) %>% 
  summarise(n_overlap=n_overlap,
            Overlapping_Genes=toString(Genes)) %>% data.frame()

T1_hyper_Bhaduri_ol_p <- merge(
  T1_hyper_Bhaduri_OR, T1_hyper_Bhaduri_ol, by="cluster", all = T)

T1_hyper_Bhaduri_ol_p <- T1_hyper_Bhaduri_ol_p[
  duplicated(T1_hyper_Bhaduri_ol_p$cluster) ==F,
]

write.table(T1_hyper_Bhaduri_ol_p, 
            sep="\t", quote = F, row.names = F,
            file="T1_Bhaduri_FisherExact.txt",
)



rsfMRI_hyper_Bhaduri_ol <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(n_overlap=length(intersect(rsfMRI_symbol$updated_GeneSymbol, 
                                       updated_GeneSymbol)),
            Genes=intersect(rsfMRI_symbol$updated_GeneSymbol, 
                            updated_GeneSymbol)) %>% 
  summarise(n_overlap=n_overlap,
            Overlapping_Genes=toString(Genes)) %>% data.frame()

rsfMRI_hyper_Bhaduri_ol_p <- merge(
  rsfMRI_hyper_Bhaduri_OR, rsfMRI_hyper_Bhaduri_ol, by="cluster", all = T)

rsfMRI_hyper_Bhaduri_ol_p <- rsfMRI_hyper_Bhaduri_ol_p[
  duplicated(rsfMRI_hyper_Bhaduri_ol_p$cluster) ==F,
]

write.table(rsfMRI_hyper_Bhaduri_ol_p, sep="\t", quote = F, row.names = F,
            file="rsfMRI_Bhaduri_FisherExact.txt",
)



diffusion_hyper_Bhaduri_ol <- Bhaduri_HGNC_ALL %>% group_by(cluster) %>%
  summarise(n_overlap=length(intersect(diffusion_symbol$updated_GeneSymbol, 
                                       updated_GeneSymbol)),
            Genes=intersect(diffusion_symbol$updated_GeneSymbol, 
                            updated_GeneSymbol)) %>% 
  summarise(n_overlap=n_overlap,
            Overlapping_Genes=toString(Genes)) %>% data.frame()

diffusion_hyper_Bhaduri_ol_p <- merge(
  diffusion_hyper_Bhaduri_OR, diffusion_hyper_Bhaduri_ol, by="cluster", all = T)

diffusion_hyper_Bhaduri_ol_p <- diffusion_hyper_Bhaduri_ol_p[
  duplicated(diffusion_hyper_Bhaduri_ol_p$cluster) ==F,
]

write.table(diffusion_hyper_Bhaduri_ol_p, 
            sep="\t", quote = F, row.names = F,
            file="diffusion_Bhaduri_FisherExact.txt",
)

