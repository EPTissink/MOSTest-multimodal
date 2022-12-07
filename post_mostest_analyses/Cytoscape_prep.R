## Prepare FUMA GSEA results as input for EnrichmentMap in Cytoscape (inspired by code from ActivePathways)
library(ActivePathways)
library(data.table)
library(dplyr)

###################
## Load FUMA output
###################
crossmodality = fread("/path/crossmodality/FUMA_gene2func86752/GS.txt",data.table=F)
diffusion = fread("/path/diffusion_specific/FUMA_gene2func86754/GS.txt",data.table=F) 
t1 = fread("/path/Genes/T1_specific/FUMA_gene2func86756/GS.txt",data.table=F)

#Extract significant GO terms
GO_crossmodality = crossmodality[startsWith(crossmodality$Category,"GO"),] 
GO_diffusion = diffusion[startsWith(diffusion$Category,"GO"),]
GO_t1 = t1[startsWith(t1$Category,"GO"),]

#Create similar GO identifier as GMT file
GO_crossmodality$term.id = paste0(toupper(sub('_', '',GO_crossmodality$Category)),sub('GO_', '_', GO_crossmodality$GeneSet))
GO_diffusion$term.id = paste0(toupper(sub('_', '',GO_diffusion$Category)),sub('GO_', '_', GO_diffusion$GeneSet))
GO_t1$term.id = paste0(toupper(sub('_', '',GO_t1$Category)),sub('GO_', '_', GO_t1$GeneSet))

#Merge into one file
all_crossmodal_t1 = merge(GO_crossmodality,GO_t1,by=c("Category","GeneSet","N_genes","link","term.id"),all=T)
all_crossmodal_t1_diffusion = merge(all_crossmodal_t1,GO_diffusion,by=c("Category","GeneSet","N_genes","link","term.id"),all=T)
all = all_crossmodal_t1_diffusion
all$term.name = tools::toTitleCase(tolower(gsub("_"," ",substr(all$GeneSet, 4, nchar(all$GeneSet)))))

###################
# terms: A data.table object with the columns 'term.id', 'term.name', 'adjusted.p.val'.
###################
all$adjusted.p.val = apply(all[,c("adjP.x","adjP.y","adjP")],1,min,na.rm=T)
terms = all[,c("term.id","term.name","adjusted.p.val")]
write.table(terms,"/path/Cytoscape/terms.txt",quote = F,row.names = F,col.names = T,sep="\t")

###################  
# gmt: An abridged GMT object containing only the pathways that were found to be significant in the ActivePathways analysis.
###################
#GMT = read.GMT("/Users/elleke/Downloads/c5.go.v2022.1.Hs.symbols.gmt.txt")
GMT <- strsplit(readLines("/path/Downloads/c5.go.v2022.1.Hs.symbols.gmt.txt"), '\t')
names(GMT) <- sapply(GMT, `[`, 1)
gmt <- lapply(GMT, function(x) { list(id=x[1], 
                                      name=tools::toTitleCase(tolower(gsub("_"," ",gsub("GOBP_","",x[1])))), 
                                      genes=x[-c(1,2)]) })
class(gmt) <- 'GMT'
gmt = gmt[names(sapply(gmt,"[",c(1))) %in% paste0(terms$term.id,".id")] #subset gmt file to significant ones
write.GMT(gmt,filename = "/path/Cytoscape/GO_enriched.gmt")

###################
# col.significance: A data.table object with a column 'term.id' and a column for each type of evidence 
# indicating whether a term was also found to be signficiant or not when considering only the genes and 
# p-values in the corresponding column of the scores matrix. If term was not found, NA's are shown in columns, otherwise the relevant lists of genes are shown.
###################
all$Crossmodality <- with(all, ifelse(is.na(p.x), "0", '1'))
all$sMRI <- with(all, ifelse(is.na(p.y), "0", '1'))
all$dMRI <- with(all, ifelse(is.na(p), "0", '1'))
col.significance = all[,c("term.id","Crossmodality","sMRI","dMRI")]
col.colors = c("#FFEDAD","#CDEECC","#FDC0BE")
instruct.str <- paste('piechart:',
                      ' attributelist="', 
                      paste(colnames(col.significance)[-1], collapse=','),
                      '" colorlist="', 
                      paste(col.colors, collapse=','), 
                      '" showlabels=FALSE', sep='')
col.significance$instruct = instruct.str

write.table(col.significance,"/path/Cytoscape/subgroups.txt",quote = F,row.names = F,col.names = T,sep="\t")
