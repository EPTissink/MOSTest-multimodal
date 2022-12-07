#! /usr/bin/env Rscript

rm(list=ls())
library(data.table)
library(ggplot2)
library(factoextra)
library(fastICA)
library(dplyr)
library(tidyr)
library(magick)

Path <- "/Users/Documents/"
setwd(Path)

### Read in Kang preproc files
ensgkey <- fread(paste0(Path,"ensg2symbol.txt"),data.table = F)
expr <- fread(paste0(Path,"GaP/Expression.csv"),data.table = F)
Rows <- fread(paste0(Path, "GaP/Rows.csv"),data.table = F)
Columns <- fread(paste0(Path, "GaP/columns_human.csv"),data.table = F)
Donors <- fread(paste0(Path, "GaP/Donors.csv"),data.table = F)
genes <- fread("/Users/Documents/MAGMA/NCBI37.3.gene.loc",data.table = F)
GenesAll <- Rows[["gene_symbol"]]

### Read in gene lists
magmaCross <- fread(paste0(Path,"crossmodality_genes_most.txt"), data.table = F)  %>% left_join(ensgkey,by=c("GENE"="GENE")) %>%
  left_join(genes,by=c("SYMBOL"="V6")) %>% filter(complete.cases(.))
magmaT1 <- fread(paste0(Path,"T1_specific_genes_most.txt"), data.table=F)  %>% left_join(ensgkey,by=c("GENE"="GENE")) %>%
  left_join(genes,by=c("SYMBOL"="V6")) %>% filter(complete.cases(.))
magmaDif <- fread(paste0(Path,"diffusion_specific_genes_most.txt"), data.table=F)  %>% left_join(ensgkey,by=c("GENE"="GENE")) %>%
  left_join(genes,by=c("SYMBOL"="V6")) %>% filter(complete.cases(.))
magmaRest <- fread(paste0(Path,"rsfmri_specific_genes_most.txt"), data.table=F)  %>% left_join(ensgkey,by=c("GENE"="GENE")) %>%
  left_join(genes,by=c("SYMBOL"="V6")) %>% filter(complete.cases(.))

### Select cortical regions
lobe="all"
#for(lobe in c("temporal","parietal","frontal","occipital")){
  if (lobe=="all") {
    CxCols <- Columns$id[Columns$structure_acronym %in% c("A1C","IPC","M1C","S1C","V1C","DFC","ITC","STC","MFC","OFC","VFC")]
  } else if (lobe=="temporal") {
    CxCols <- Columns$id[Columns$structure_acronym %in% c("A1C","ITC","STC")]
  } else if (lobe=="parietal") {
    CxCols <- Columns$id[Columns$structure_acronym %in% c("IPC","S1C")]
  } else if (lobe=="frontal") {
    CxCols <- Columns$id[Columns$structure_acronym %in% c("M1C","DFC","MFC","OFC","VFC")]
  } else if (lobe=="occipital") {
    CxCols <- Columns$id[Columns$structure_acronym %in% c("V1C")]
  } else {
    print("crap")
}

expr1 <- expr[,c("probe_num",CxCols)]
Columns1 <- Columns[Columns$id %in% CxCols,]

### Select probe with highest differential stability
Rows <- Rows[order(-Rows$DS),]
Rows <- Rows[match(unique(Rows$gene_symbol),Rows$gene_symbol),]

### make NCX df: mean over regions per donor
for(donor in unique(Columns1$donor_name)){
  colsdonor <- Columns1$id[Columns1$donor_name==donor]
   if (length(colsdonor)==1){
     expr1$tmp <- expr1[colsdonor]
     } else {expr1$tmp <-rowMeans(expr1[colsdonor])}
  colnames(expr1)[colnames(expr1)=="tmp"] <- donor
}
exprNCX <- expr1[,c("probe_num",unique(Columns1$donor_name))]

#### Select genes, and combine into plotting df
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

exprNCX <- as.data.frame(cbind(exprNCX[,1],t(apply(exprNCX[,-1],1,scale01))))
exprNCX[-1] <- exprNCX[-1] * 100
colnames(exprNCX)[1] <- "probe_num"

Prenatal <- unique(Columns1$donor_name[Columns1$DayAge<260])
Postnatal <- unique(Columns1$donor_name[Columns1$DayAge>260])

CrossDat <- exprNCX[exprNCX$probe_num %in% Rows$probe_num[Rows$gene_symbol %in% magmaCross$SYMBOL],]
CrossDat1 <- gather(CrossDat,Donor, Expression, 2:ncol(exprNCX)) %>% mutate(GeneSet="Cross-modality (n=1,275)") 
CrossDat1$Phase <- ifelse(CrossDat1$Donor %in% Prenatal,"Prenatal",ifelse(CrossDat1$Donor %in% Postnatal,"Postnatal",NA))

T1Dat <- exprNCX[exprNCX$probe_num %in% Rows$probe_num[Rows$gene_symbol %in% magmaT1$SYMBOL],]
T1Dat1 <- gather(T1Dat,Donor, Expression, 2:ncol(exprNCX)) %>% mutate(GeneSet="T1-weighted (n=515)") 
T1Dat1$Phase <- ifelse(T1Dat1$Donor %in% Prenatal,"Prenatal",ifelse(T1Dat1$Donor %in% Postnatal,"Postnatal",NA))

DifDat <- exprNCX[exprNCX$probe_num %in% Rows$probe_num[Rows$gene_symbol %in% magmaDif$SYMBOL],]
DifDat1 <- gather(DifDat,Donor, Expression, 2:ncol(exprNCX)) %>% mutate(GeneSet="Diffusion (n=380)") 
DifDat1$Phase <- ifelse(DifDat1$Donor %in% Prenatal,"Prenatal",ifelse(DifDat1$Donor %in% Postnatal,"Postnatal",NA))

RestDat <- exprNCX[exprNCX$probe_num %in% Rows$probe_num[Rows$gene_symbol %in% magmaRest$SYMBOL],]
RestDat1 <- gather(RestDat,Donor, Expression, 2:ncol(exprNCX)) %>% mutate(GeneSet="rs-fMRI (n=5)") 
RestDat1$Phase <- ifelse(RestDat1$Donor %in% Prenatal,"Prenatal",ifelse(RestDat1$Donor %in% Postnatal,"Postnatal",NA))

PlotDF <- rbind(CrossDat1,T1Dat1,DifDat1) #RestDat1
PlotDF <- left_join(PlotDF,Donors[,c("donor_name","DayAge")],by=c("Donor"="donor_name"))

##### Some plotting prettifying params
int <- list()
int[["developmental_stages"]]  <- c(12*7,24*7,40*7+19*30.44,40*7+6*365.25,40*7+12*365.25,40*7+20*365.25,40*7+65*365.25)
int[["developmental_labels"]]  <- c("2nd-Trimester","3rd-Trimester","Early childhood", "Late childhood", "Adolescence", "Adulthood","Old age")
int[["developmental_age"]]  <- c("12 weeks","24 weeks","19 month", "6 years", "12 years ", "20 years", "65 years")
int[["life_stages"]]  <-data.frame(stages=int$developmental_stages,labels=int$developmental_labels,age=int$developmental_age)
int[["xlim"]]  <- log10(c(70,30230.5))
int[["birth_date"]]  <- c(40*7)
int[["birth_label"]]  <- c("Birth")
FUN=log10
PlotDF$GeneSet <- as.factor(PlotDF$GeneSet)

ypos=39.5 # min ylim -2
library(ggformula)
############# line plot per gene set
BS <- ggplot(data=PlotDF,aes(x=FUN(DayAge),y=Expression,color=GeneSet)) + theme_bw() + 
  geom_smooth() +
  scale_color_manual(values=c('black',"#1B9E77","#D95F02"))+ #,"#7570B3"
  geom_vline(xintercept=FUN(int$developmental_stages), size=.7, color="lightgray") + 
  geom_vline(xintercept=FUN(int$birth_date), size=1.2, color="darkgray") + 
  ylab("Mean gene expression \n") + 
  coord_cartesian(ylim=c(40,56),expand=F, clip = "off") + 
  theme(plot.margin=unit(c(1,1,5,1.2),"cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 18),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 14),
        axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position=c(0.82,0.85),
        legend.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'))
ggsave(paste0(Path,"BrainSpan_SimpleMean_220727.png"), plot = BS, width = 10, height = 7,bg="transparent")
