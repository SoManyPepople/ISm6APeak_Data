
options(stringsAsFactors = F)
setwd("/data/m6A_calling_strategy/Analysis")
.libPaths(c("/home/huangp/R/x86_64-pc-linux-gnu-library/4.2","/opt/R/4.2.0/lib/R/library", "/home/huangp/anaconda3/envs/m6A_seq/lib/R/library"))
options(scipen = 9)
for(p in c("data.table","tidyverse","foreach","doParallel","ggplot2","ggsci","ggpubr","plotgardener")){
  suppressWarnings(library(p,character.only = T))
}
options(width=300)
options(bedtools.path = "~/anaconda3/envs/m6A_seq/bin")
library(bedtoolsr)
library(mysterypackage)
source("/data/m6A_calling_strategy/Script/Wrapper_function_for_m6A_calling_using_MACS2SMePeak.R")

### function to enrichment #########
#self made overrepresentation fuction
OR_enrichment <- function(InputGenes=dt.Genelist$ENTREZID,Type=c("GO_BP","GO_MF","KEGG","RA")[1], Organism=c("human","mouse")[2]){
  if(Organism=="human"){
    library(org.Hs.eg.db)
    OrgDb <- org.Hs.eg.db
  }else{library(org.Mm.eg.db)
    OrgDb <- org.Mm.eg.db}
  if(Type == "GO_BP"){
    enrich_res.dt <- enrichGO(gene=InputGenes,OrgDb = OrgDb, ont="BP",pvalueCutoff = 0.05,readable = T)@result %>% as.data.table() %>% dplyr::filter(Count>=3 & pvalue<=0.05)
    if(nrow(enrich_res.dt)>0){enrich_res.dt <- enrich_res.dt %>% dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,geneName=geneID,Count) %>%
      mutate(Ontology=Type) %>% dplyr::arrange(p.adjust)}else{ enrich_res.dt <- NULL}
  }else if(Type=="GO_MF"){
    enrich_res.dt <- enrichGO(gene=InputGenes,OrgDb = OrgDb, ont="MF",pvalueCutoff = 0.05,readable = T)@result %>% as.data.table() %>% dplyr::filter(Count>=3 & pvalue<=0.05)
    if(nrow(enrich_res.dt)>0){enrich_res.dt <- enrich_res.dt %>% dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,geneName=geneID,Count) %>%
      mutate(Ontology=Type) %>% dplyr::arrange(p.adjust)}else{ enrich_res.dt <- NULL}
  }else if(Type == "KEGG"){
    enrich_res.dt <- setReadable(enrichKEGG(gene =InputGenes, use_internal_data = T,pvalueCutoff = 0.05,organism = ifelse(Organism=="human","hsa","mmu")),OrgDb = OrgDb, keyType = "ENTREZID")@result %>% as.data.table() %>% dplyr::filter(Count>=3 & pvalue<=0.05)
    all_path = as.data.frame(KEGG.db::KEGGPATHID2NAME)
    enrich_res.dt$Description = all_path$path_name[match(enrich_res.dt$ID, all_path$path_id)]
    if(nrow(enrich_res.dt)>0){enrich_res.dt <- enrich_res.dt %>% dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,geneName=geneID,Count) %>%
      mutate(Ontology=Type) %>% dplyr::arrange(p.adjust)}else{ enrich_res.dt <- NULL}
  }else if(Type == "RA"){
    enrich_res.dt <- ReactomePA::enrichPathway(InputGenes, pvalueCutoff = 0.05,organism = Organism,readable = T)@result %>% as.data.table() %>% dplyr::filter(Count>=3 & pvalue<=0.05)
    if(nrow(enrich_res.dt)>0){enrich_res.dt <- enrich_res.dt %>% dplyr::select(ID,Description,GeneRatio,BgRatio,pvalue,p.adjust,geneName=geneID,Count) %>% mutate(Ontology=Type) %>% dplyr::arrange(p.adjust)}else{ enrich_res.dt <- NULL}
  }
  return(enrich_res.dt)
}
#function to deredundant GO terms
TermDeredundant <- function(InputEnrichRes = dt.input, Overlapcutoff=0.8,nthread=10){
  if(nrow(InputEnrichRes)<=10){nthread<-1}
  #build list for all terms
  require(foreach)
  TermGeneList <- foreach(i=1:nrow(InputEnrichRes))%do%{
    str_split( InputEnrichRes$geneName[i],pattern = "/") %>% unlist()
  }
  names(TermGeneList) <- InputEnrichRes$ID
  #calculate the overlap cutoff
  require(doParallel)
  registerDoParallel(cl=nthread)
  dt.overlap.ratio <- foreach(j = 1:(nrow(InputEnrichRes)-1),.combine='rbind')%dopar%{
    foreach(k=(j+1):nrow(InputEnrichRes), .combine='rbind')%do%{
      nSize1 <- length(TermGeneList[[j]])
      nSize2 <- length(TermGeneList[[k]])
      nOverlap <- length(intersect(TermGeneList[[j]],TermGeneList[[k]]))
      data.table(Term1=names(TermGeneList)[j],nSize1, Term2=names(TermGeneList)[k], nSize2, nOverlap)
    }
  }
  stopImplicitCluster()
  dt.filter <- dt.overlap.ratio %>% dplyr::filter(nOverlap/nSize2>=Overlapcutoff)
  return(InputEnrichRes %>% dplyr::filter(!ID %in% dt.filter$Term2))
}
#function to transform gene ratio
mixedToFloat <- function(x){
  foreach(i=1:length(x), .combine = 'c')%do%{
    x1 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",1) %>% as.integer()
    x2 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",2) %>% as.integer()
    x1/x2*100
  }
}

##############  functional enrichment of Novel m6A+ gene detected via optimal parameters ###########################

library(KEGG.db)
library(genekitr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
# library(BiocParallel)
# register(SerialParam())#turn off parallel

#load Novel m6A+ gene
load("sm6APeak_Figure2.intermediate.results.RDS")
#Novel m6A+ gene in HEK293
dt.exomePeak2.BestCombo.specific.peak.benchmark.HEK293 %>% dplyr::filter(BestCombo.specific==TRUE)#
dt.sm6APeak.novel.shared.m6AGene <- dt.exomePeak2.BestCombo.specific.peak.benchmark.HEK293 %>% dplyr::filter(BestCombo.specific==TRUE) %>%
  separate_longer_delim(cols = "OverlappedGenes",delim = ",") %>% as.data.table() %>% group_by(OverlappedGenes) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>%
  distinct(OverlappedGenes,nSample) %>% dplyr::filter(nSample>=1)#1653 shared novel m6AGene
#function to perform multiple ontology
Genelist <- dt.sm6APeak.novel.shared.m6AGene[OverlappedGenes != "NA" & !is.na(OverlappedGenes),unique(OverlappedGenes)]
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)#1382 unique novel m6A+gene
Genelist[!Genelist %in% dt.Genelist$SYMBOL]
#start enrichment
dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA"), .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type,Organism = c("human","mouse")[1])
}
dt.enrich.res[,.N,by=.(Ontology)]
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
dt.enrich.res[p.adjust<=0.1,.N,by=.(Ontology)]

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.optimal.parameter.HEK293.RDS")

#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(pvalue<0.01) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table() %>% dplyr::filter(str_detect(Ontology,"GO"))
#visualization of functional enrichment 
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust),lgPvalue= -log10(pvalue), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgPvalue),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.OptimalParameter.HEK293 <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgPvalue, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.OptimalParameter.HEK293

#Novel m6A+ gene in mESC
dt.exomePeak2.BestCombo.specific.peak.benchmark.mESC %>% dplyr::filter(BestCombo.specific==TRUE)#
dt.sm6APeak.novel.shared.m6AGene <- dt.exomePeak2.BestCombo.specific.peak.benchmark.mESC %>% dplyr::filter(BestCombo.specific==TRUE) %>%
  separate_longer_delim(cols = "OverlappedGenes",delim = ",") %>% as.data.table() %>% group_by(OverlappedGenes) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>%
  distinct(OverlappedGenes,nSample) %>% dplyr::filter(nSample>=1)#1599 shared novel m6AGene
#function to perform multiple ontology
Genelist <- dt.sm6APeak.novel.shared.m6AGene[OverlappedGenes != "NA" & !is.na(OverlappedGenes),unique(OverlappedGenes)]
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Mm.eg.db) %>% as.data.table()#1517 unique novel m6A+gene
dt.Genelist
Genelist[!Genelist %in% dt.Genelist$SYMBOL]
#start enrichment
dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA")[1:2], .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type,Organism = c("human","mouse")[2])
}
dt.enrich.res[,.N,by=.(Ontology)]
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
dt.enrich.res[p.adjust<=0.1,.N,by=.(Ontology)]

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.optimal.parameter.mESC.RDS")

#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(pvalue<0.01 & p.adjust<=0.05) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table() %>% dplyr::filter(str_detect(Ontology,"GO"))
#visualization of functional enrichment 
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust),lgPvalue= -log10(pvalue), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgFDR),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.OptimalParameter.mESC <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgPvalue, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.OptimalParameter.mESC

library(plotgardener)
# pdf("FigureS_functional_enrichment_result_of_Novel_m6A_modified_genes_revealed_by_performance_improvement.pdf",width = 8.2677, height = 11.693)
#six subplots (three improvement, and one improvement in CRC)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.top.enriched.term.novel.m6AGene.OptimalParameter.mESC+ggtitle(label="Optimal parameter revealed novel m6A+ gene in mESC"), x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.novel.m6AGene.OptimalParameter.HEK293+ggtitle(label="Optimal parameter revealed novel m6A+ gene in HEK293"), x = 9.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 9, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)


dev.off()


##############  functional enrichment of Novel m6A+ gene detected via strand-specific intensity ###########################

library(KEGG.db)
library(genekitr)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
# library(BiocParallel)
# register(SerialParam())#turn off parallel

#load Novel m6A+ gene
load("sm6APeak_Figure3.intermediate.results.RDS")
#Novel m6A+ gene in HEK293
dt.MACS2.BestOption.specific.peak.benchmark.HEK293 %>% dplyr::filter(BestOption.specific==TRUE)#
dt.sm6APeak.novel.shared.m6AGene <- dt.MACS2.BestOption.specific.peak.benchmark.HEK293 %>% dplyr::filter(BestOption.specific==TRUE) %>%
  separate_longer_delim(cols = "OverlappedGenes",delim = ",") %>% as.data.table() %>% group_by(OverlappedGenes) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>%
  distinct(OverlappedGenes,nSample) %>% dplyr::filter(nSample>=1)#1766 shared novel m6AGene
#function to perform multiple ontology
Genelist <- dt.sm6APeak.novel.shared.m6AGene[OverlappedGenes != "NA" & !is.na(OverlappedGenes),unique(OverlappedGenes)]
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Hs.eg.db) %>% as.data.table()#1559 unique novel m6A+gene
dt.Genelist
Genelist[!Genelist %in% dt.Genelist$SYMBOL]
#start enrichment
dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA")[1:2], .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type,Organism = c("human","mouse")[1])
}
dt.enrich.res[,.N,by=.(Ontology)]
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
dt.enrich.res[p.adjust<=0.1,.N,by=.(Ontology)]

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.best.option.HEK293.RDS")

#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(pvalue<0.001) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table() %>% dplyr::filter(str_detect(Ontology,"GO"))
#visualization of functional enrichment 
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust),lgPvalue= -log10(pvalue), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgPvalue),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.BestOption.HEK293 <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgPvalue, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.BestOption.HEK293

#Novel m6A+ gene in mESC
dt.MACS2.BestOption.specific.peak.benchmark.mESC %>% dplyr::filter(BestOption.specific==TRUE)#
dt.sm6APeak.novel.shared.m6AGene <- dt.MACS2.BestOption.specific.peak.benchmark.mESC %>% dplyr::filter(BestOption.specific==TRUE) %>%
  separate_longer_delim(cols = "OverlappedGenes",delim = ",") %>% as.data.table() %>% group_by(OverlappedGenes) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>%
  distinct(OverlappedGenes,nSample) %>% dplyr::filter(nSample>=1)#2150 shared novel m6AGene
#function to perform multiple ontology
Genelist <- dt.sm6APeak.novel.shared.m6AGene[OverlappedGenes != "NA" & !is.na(OverlappedGenes),unique(OverlappedGenes)]
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Mm.eg.db) %>% as.data.table()#2065 unique novel m6A+gene
dt.Genelist
Genelist[!Genelist %in% dt.Genelist$SYMBOL]
#start enrichment
dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA")[1:2], .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type,Organism = c("human","mouse")[2])
}
dt.enrich.res[,.N,by=.(Ontology)]
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
dt.enrich.res[p.adjust<=0.1,.N,by=.(Ontology)]

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.best.option.mESC.RDS")

#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(pvalue<0.01 & p.adjust<=0.05) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table() %>% dplyr::filter(str_detect(Ontology,"GO"))
#visualization of functional enrichment 
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust),lgPvalue= -log10(pvalue), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgFDR),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.BestOption.mESC <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgPvalue, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.BestOption.mESC

library(plotgardener)
# pdf("FigureS_functional_enrichment_result_of_Novel_m6A_modified_genes_revealed_by_performance_improvement.pdf",width = 8.2677, height = 11.693)
#six subplots (three improvement, and one improvement in CRC)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.top.enriched.term.novel.m6AGene.BestOption.mESC+ggtitle(label="Strand-specific intensity revealed novel m6A+ gene in mESC"), x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.novel.m6AGene.BestOption.HEK293+ggtitle(label="Strand-specific intensity revealed novel m6A+ gene in HEK293"), x = 9.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 9, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)


dev.off()




##############  functional enrichment of Novel m6A+ gene identified by ISm6APeak in mouse blastocyst ###########################

library(KEGG.db)
library(genekitr)
library(org.Mm.eg.db)
library(clusterProfiler)
# library(BiocParallel)
# register(SerialParam())#turn off parallel

#load Novel m6A+ gene
load("sm6APeak_Figure4.intermediate.results.RDS")
dt.sm6APeak.novel.shared.m6AGene %>% dplyr::filter(IsNovelGene)#989 novel m6A+ gene
dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6AGene %>% dplyr::filter(IsNovelGene==TRUE) %>% separate_longer_delim(cols = "OverlappedGenes",delim = ",") %>% as.data.table()
#function to perform multiple ontology
Genelist <- dt.sm6APeak.novel.shared.m6AGene[OverlappedGenes != "NA" & !is.na(OverlappedGenes),unique(OverlappedGenes)]#866 unique novel m6A+gene
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Mm.eg.db)
Genelist[!Genelist %in% dt.Genelist$SYMBOL]

dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA"), .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type)
}
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
# 1:    GO_BP 159
# 2:    GO_MF  17 

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.ISm6APeak.in.mouse.blastocyst.RDS")


#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(p.adjust<=0.05) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table()
#visualization of functional enrichment 
mixedToFloat <- function(x){
  foreach(i=1:length(x), .combine = 'c')%do%{
    x1 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",1) %>% as.integer()
    x2 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",2) %>% as.integer()
    x1/x2*100
  }
}
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgFDR),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.blastocyst  <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgFDR, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.blastocyst

library(plotgardener)
# pdf("FigureS_functional_enrichment_result_of_Novel_m6A_modified_genes_revealed_by_performance_improvement.pdf",width = 8.2677, height = 11.693)
#six subplots (three improvement, and one improvement in CRC)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.top.enriched.term.novel.m6AGene.blastocyst+ggtitle(label = "ISm6APeak revealed novel m6A+ gene in mouse blastocyst"), x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)



dev.off()


##############  functional enrichment of Novel m6A+ gene identified by ISm6APeak in CRC ###########################

library(KEGG.db)
library(genekitr)
library(org.Hs.eg.db)
library(clusterProfiler)
# library(BiocParallel)
# register(SerialParam())#turn off parallel

#load Novel m6A+ gene
load("sm6APeak_Figure5_intermediate.results.RDS")
#Novel m6A+ gene in CRC
Reduce(x=strict.novel.concurrent.m6a.gene.list[1:3],intersect) %>% str()#N=636
#function to perform multiple ontology
Genelist <- Reduce(x=strict.novel.concurrent.m6a.gene.list[1:3],intersect) %>% unique()
dt.Genelist <- clusterProfiler::bitr(Genelist,fromType = c("SYMBOL"), toType = c("ENTREZID"),OrgDb = org.Hs.eg.db) %>% as.data.table()#498 unique novel m6A+gene
dt.Genelist
Genelist[!Genelist %in% dt.Genelist$SYMBOL]
#start enrichment
dt.enrich.res <- foreach(type=c("GO_BP","GO_MF","KEGG","RA")[1:2], .combine='rbind')%do%{
  message(paste0("Perform enrichment for ",type))
  OR_enrichment(InputGenes=dt.Genelist$ENTREZID %>% unique(),Type=type,Organism = c("human","mouse")[1])
}
dt.enrich.res[,.N,by=.(Ontology)]
dt.enrich.res[p.adjust<=0.05,.N,by=.(Ontology)]
# Ontology   N
dt.enrich.res[p.adjust<=0.1,.N,by=.(Ontology)]

saveRDS(dt.enrich.res, file="dt.enrich.res.novel.m6AGene.ISm6APeak.CRC.RDS")

#deredundant of enriched term
dt.input <- dt.enrich.res %>% dplyr::filter(pvalue<0.01) %>% dplyr::arrange(p.adjust)
deredundant_enrichment_res <- TermDeredundant(InputEnrichRes = dt.input,Overlapcutoff = 0.8,nthread = 10)
#select top5 terms for GOBP and GOMF
dt.top.enriched.term.strict.novel.m6agene <- deredundant_enrichment_res %>% group_by(Ontology) %>% dplyr::arrange(p.adjust) %>% slice_head(n=5) %>% as.data.table() %>% dplyr::filter(str_detect(Ontology,"GO"))
#visualization of functional enrichment 
mixedToFloat <- function(x){
  foreach(i=1:length(x), .combine = 'c')%do%{
    x1 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",1) %>% as.integer()
    x2 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",2) %>% as.integer()
    x1/x2*100
  }
}
dt.top.enriched.term <- dt.top.enriched.term.strict.novel.m6agene %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust),lgPvalue= -log10(pvalue), Description=stringr::str_wrap(Description,width = 80)) %>% as.data.table() %>%
  dplyr::arrange(desc(lgPvalue),Description) %>% mutate(Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.novel.m6AGene.ISm6APeak.CRC <- ggplot(dt.top.enriched.term , aes(x = EnrichFold, y = Description, color = lgPvalue, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(floor(min(dt.top.enriched.term$EnrichFold)),ceiling(max(dt.top.enriched.term$EnrichFold))),breaks = seq(0,max(dt.top.enriched.term$EnrichFold)*1.1,1))+
  guides(size=guide_legend(nrow = 1))+
  labs(x="Enrichment fold",y=NULL)+
  # facet_wrap( ~Group, scales = "free_y",ncol = 1) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.top.enriched.term.novel.m6AGene.ISm6APeak.CRC


library(plotgardener)
# pdf("FigureS_functional_enrichment_result_of_Novel_m6A_modified_genes_revealed_by_performance_improvement.pdf",width = 8.2677, height = 11.693)
#six subplots (three improvement, and one improvement in CRC)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.top.enriched.term.novel.m6AGene.ISm6APeak.CRC+ggtitle(label="ISm6APeak revealed novel m6A+ gene in CRC"), x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)


dev.off()


################# combine all enrichment plot together ##########
library(plotgardener)
pdf("FigureS_functional_enrichment_result_of_Novel_m6A_modified_genes_revealed_by_performance_improvement.pdf",width = 8.2677, height = 11.693)
#six subplots (three improvement, and one improvement in CRC)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.top.enriched.term.novel.m6AGene.OptimalParameter.mESC+ggtitle(label="Optimal parameter revealed novel m6A+ gene in mESC"), x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.novel.m6AGene.OptimalParameter.HEK293+ggtitle(label="Optimal parameter revealed novel m6A+ gene in HEK293"), x = 9.2, y=0.05, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 9, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)

plotGG(plot = p.top.enriched.term.novel.m6AGene.BestOption.mESC+ggtitle(label="Strand-specific intensity revealed novel m6A+ gene in mESC"), x = 0.2, y=8, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 8,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.novel.m6AGene.BestOption.HEK293+ggtitle(label="Strand-specific intensity revealed novel m6A+ gene in HEK293"), x = 9.2, y=8, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 9, y = 8,just = c("top","left"), default.units = "cm",draw=T)

plotGG(plot = p.top.enriched.term.novel.m6AGene.blastocyst+ggtitle(label = "ISm6APeak revealed novel m6A+ gene in mouse blastocyst"), x = 0.2, y=16, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05, y = 16,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.novel.m6AGene.ISm6APeak.CRC+ggtitle(label="ISm6APeak revealed novel m6A+ gene in CRC"), x = 9.2, y=16, default.units = "cm",width = 8.8, height = 7.8)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 9, y = 16,just = c("top","left"), default.units = "cm",draw=T)


dev.off()


