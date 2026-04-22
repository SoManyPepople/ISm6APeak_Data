options(stringsAsFactors = F)
setwd("/data/m6A_calling_strategy/Analysis")
.libPaths(c("/home/huangp/R/x86_64-pc-linux-gnu-library/4.2","/opt/R/4.2.0/lib/R/library", "/home/huangp/anaconda3/envs/m6A_seq/lib/R/library"))
options(scipen = 9)
for(p in c("data.table","tidyverse","foreach","doParallel","ggplot2","ggsci","ggpubr","plotgardener","bedtoolsr")){
  suppressWarnings(library(p,character.only = T))
}
options(width=300)
options(bedtools.path = "~/anaconda3/envs/m6A_seq/bin")
library(bedtoolsr)
library(mysterypackage)
library(BSgenome.Mmusculus.UCSC.mm39)
source("/data/m6A_calling_strategy/Script/Wrapper_function_for_m6A_calling_using_MACS2SMePeak.R")

### Part1 IVT m6A-seq experiment  the benchmark of m6A calling algorithm ###############
#figure 1a. The count of curated consistent GLORI and eTAMseq m6A sites with highly ratio of classic motif

#load HEK benchmark m6As
HEK293.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="HEK293",value=T)
names(HEK293.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="HEK293",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.HEK293 <-  foreach(f = 1:length(HEK293.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(HEK293.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}
#load mESC benchmark m6As
mESC.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="mESC",value=T)
names(mESC.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="mESC",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.mESC <-  foreach(f = 1:length(mESC.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(mESC.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}

dt.HEK293.benchmark.m6As <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/dt.HEK293.benchmark.m6As.RDS")
dt.mESC.benchmark.m6As <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/dt.mESC.benchmark.m6As.RDS")
dt.HeLa.benchmark.m6As <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/dt.HeLa.benchmark.m6As.RDS")
All_RRACH <- paste0(rep(paste0(c("AA","AG","GA","GG"), "AC"),3),rep(c("A","C","T"),each=4))
All_DRAC <-  paste0(rep(c("A","G","T"),each=2), rep(c("A","G"),3), "AC")
All_DRACH <- paste0(rep(All_DRAC, 3), rep(c("A", "C", "T"), each=6))
dt.count.RRACHratio.benchmark.m6As <- rbind(dt.HEK293.benchmark.m6As %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Seq),
                                            dt.mESC.benchmark.m6As %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Seq),
                                            dt.HeLa.benchmark.m6As %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Seq) ) %>% group_by(benchmark_name=Sample) %>%
  mutate(nSites=n_distinct(name), nSiteDRACH=n_distinct(name[Seq %in% All_DRACH])) %>% as.data.table() %>%
  distinct(benchmark_name,nSites,nSiteDRACH) %>% mutate(ratioDRACH=nSiteDRACH/nSites)

dt.count.RRACHratio.benchmark.m6As <- dt.count.RRACHratio.benchmark.m6As %>%
  arrange(desc(ratioDRACH)) %>%
  mutate(
    benchmark_name = factor(benchmark_name, levels = benchmark_name),
    label_pct = scales::percent(ratioDRACH, accuracy = 0.1),
    label_counts = paste0(nSiteDRACH, " / ", nSites)
  )
dt.count.RRACHratio.benchmark.m6As <- dt.count.RRACHratio.benchmark.m6As %>% dplyr::filter(benchmark_name %in% c("GLORI_HEK293","GLORI_mESC","eTAMseq_mESC","GLORI_HeLa","eTAMseq_HeLa")) %>%
  dplyr::arrange(desc(nSites)) %>% 
  mutate(benchmark_name = factor(benchmark_name, levels=unique(benchmark_name)))
p.count.RRACHratio.benchmark.m6As <- ggplot(dt.count.RRACHratio.benchmark.m6As,
                                            aes(y = benchmark_name, x = ratioDRACH, fill = benchmark_name)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  # counts inside the bar (centered horizontally)
  geom_text(aes(label = label_counts, x = ratioDRACH / 2),
            color = "white", size = 2, fontface = "plain", hjust = 0.5) +
  # percent label to the right of the bar
  geom_text(aes(label = label_pct, x = ratioDRACH - 0.03),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, 1.05), expand = c(0, 0)) +
  scale_fill_manual(values = c("#ffd47f","#7b93c6","#acd9ee","#c82406", "#f7c1cf"),breaks = c("GLORI_HEK293","GLORI_mESC","eTAMseq_mESC","GLORI_HeLa","eTAMseq_HeLa")) +
  labs(x = "DRACH ratio", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.RRACHratio.benchmark.m6As

pdf("new_Fig1A_count_and_RRACHratio_of_benchmark_m6As.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.count.RRACHratio.benchmark.m6As, x = 0.05, y=0.2, default.units = "cm",width = 5.8, height = 5.8)

dev.off()

#figure 1b. The overlap of GLORI, eTAMseq, and SACseq m6A sites and m6A+ genes
# dt.shared.benchmark.m6As.HEK293 <- dt.benchmark.m6A.HEK293 %>% group_by(name) %>% dplyr::filter(n_distinct(benchmark_name)==2) %>% as.data.table() %>%
#   dplyr::arrange(name,benchmark_name) %>% mutate(score=case_when(benchmark_name=="GLORI_HEK293" ~ score, benchmark_name=="SACseq_HEK293" ~ score/100)) %>%
#   pivot_wider(id_cols = c("name"), values_from = "score", names_from = "benchmark_name") %>% as.data.table()
# cor.test(dt.shared.benchmark.m6As.HEK293$GLORI_HEK293, dt.shared.benchmark.m6As.HEK293$SACseq_HEK293,method = "spearman")
dt.shared.benchmark.m6As.mESC <- dt.benchmark.m6A.mESC %>% group_by(name) %>% dplyr::filter(n_distinct(benchmark_name)==2) %>% as.data.table() %>%
  dplyr::arrange(name,benchmark_name) %>%
  pivot_wider(id_cols = c("name"), values_from = "score", names_from = "benchmark_name") %>% as.data.table()
Cor <- cor.test(dt.shared.benchmark.m6As.mESC$GLORI_mESC, dt.shared.benchmark.m6As.mESC$eTAMseq_mESC,method = "pearson")
library(ggpointdensity)
p.shared.benchmark.m6A.sites <- ggplot(dt.shared.benchmark.m6As.mESC, aes(x = GLORI_mESC, y = eTAMseq_mESC)) +
  ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d') +
  viridis::scale_color_viridis(name = "Point density", option = "magma", trans = "log") +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  labs(x = "GLORI m6A level in mESC", y = "eTAMseq m6A level in mESC", title = paste0(nrow(dt.shared.benchmark.m6As.mESC), " shared benchmark m6A sites"),
       subtitle = paste0("Cor=",round(Cor$estimate,2)," ",ifelse(Cor$p.value<2.2e-16, "p-value<2.2e-16", paste0("p-value=",signif(Cor$p.value))))) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = "none", plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.15)), plot.subtitle = element_text(hjust = 0.5))

ggsave(plot=p.shared.benchmark.m6A.sites, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 5.8,filename = "Figure1_p.shared.benchmark.m6A.sites.pdf")

#whether GLORI_HeLa and eTAMseq_HeLa display highly consistent m6A sites and m6A level?
dt.HeLa.benchmark.m6As[,.(N=n_distinct(name)),by=.(Sample)]
# 1: eTAMseq_HeLa  70635
# 2:   GLORI_HeLa 128654
# 3:  SACseq_HeLa  34989
unique(dt.HeLa.benchmark.m6As$name) %>% length()#165,759 sites
dt.shared.benchmark.m6As.HeLa <- dt.HeLa.benchmark.m6As %>% filter(Sample != "SACseq_HeLa") %>% #don't use SACseq results since its weak correlation with GLORI HeLa
  group_by(name) %>% dplyr::filter(n_distinct(Sample)>=2) %>% mutate(SharedSample=paste(unique(Sample),collapse = ",")) %>%
  as.data.table() %>% dplyr::arrange(name,Sample)#42,642 shared m6As at least in two datasets

dt.shared.benchmark.m6As.HeLa <- dt.shared.benchmark.m6As.HeLa %>% dplyr::filter(str_detect(SharedSample,"GLORI") & str_detect(SharedSample,"eTAM")) %>%
  filter(Sample %in% c("eTAMseq_HeLa","GLORI_HeLa")) %>% pivot_wider(id_cols = c("name"), values_from = "score", names_from = "Sample") %>% as.data.table()#42,642
Cor.HeLa <- cor.test(dt.shared.benchmark.m6As.HeLa$GLORI_HeLa, dt.shared.benchmark.m6As.HeLa$eTAMseq_HeLa,method = "pearson")
Cor.HeLa#cor=0.95

# dt.shared.benchmark.m6As.HeLa.GLORI_plus_SAC <- dt.shared.benchmark.m6As.HeLa %>% dplyr::filter(str_detect(SharedSample,"GLORI") & str_detect(SharedSample,"SAC")) %>%
#   filter(Sample %in% c("SACseq_HeLa","GLORI_HeLa")) %>% mutate(score=case_when(score>1 ~ score/100, .default = score)) %>%
#   pivot_wider(id_cols = c("name"), values_from = "score", names_from = "Sample",values_fn = mean) %>% as.data.table()#23,272
# dt.shared.benchmark.m6As.HeLa.GLORI_plus_SAC[is.na(GLORI_HeLa) | is.na(SACseq_HeLa),]
# Cor.GLORI_plus_SAC <- cor.test(dt.shared.benchmark.m6As.HeLa.GLORI_plus_SAC$GLORI_HeLa, dt.shared.benchmark.m6As.HeLa.GLORI_plus_SAC$SACseq_HeLa)
# Cor.GLORI_plus_SAC#0.39, give up strange result

library(ggpointdensity)
p.HeLa.shared.benchmark.m6A.sites <- ggplot(dt.shared.benchmark.m6As.HeLa, aes(x = GLORI_HeLa, y = eTAMseq_HeLa)) +
  ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d') +
  viridis::scale_color_viridis(name = "Point density", option = "magma", trans = "log") +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  labs(x = "GLORI m6A level in HeLa", y = "eTAMseq m6A level in HeLa", title = paste0(nrow(dt.shared.benchmark.m6As.HeLa), " shared benchmark m6A sites"),
       subtitle = paste0("Cor=",round(Cor$estimate,2)," ",ifelse(Cor$p.value<2.2e-16, "p-value<2.2e-16", paste0("p-value=",signif(Cor$p.value))))) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = "none", plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.15)), plot.subtitle = element_text(hjust = 0.5))

ggsave(plot=p.HeLa.shared.benchmark.m6A.sites, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 5.8,filename = "FigureS1_p.HeLa.shared.benchmark.m6A.sites.pdf")

#figure 1c. a proportion of confident benchmark m6As with further m6A pertubation validation
# dt.eTAM.mESC.confident <- dt.eTAM.mESC.consistent %>% left_join(x=.,y=eTAM.mESC.perturbation %>% dplyr::select(name,score.Mettl3cko=score),by=c("name")) %>%
#   mutate(PerturbEffect=score.Mettl3cko-score) %>% filter(nRep>1 | PerturbEffect <= -0.1) %>%
#   mutate(confidence =  dplyr::case_when(nRep>1 & PerturbEffect <= -0.5 ~ "high",
#                                         nRep >1 & PerturbEffect <= -0.1 ~ "medium", .default = "low") %>% factor(levels=c("low","medium","high"))) %>%
#   dplyr::arrange(name)
# dt.eTAM.mESC.confident[,.N,keyby=.(confidence)]
# #    confidence     N
# # 1:        low 20148
# # 2:     medium  5439
# # 3:       high   269
dt.eTAMseq.mESC.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/eTAMseq/dt.eTAMseq.mESC.confident.m6As.RDS")
# dt.GLORI.HEK293.confident <- dt.GLORI.HEK293.confident %>% filter(PerturbEffect.FTO <= -0.1 | PerturbEffect.METTL3i <= -0.1) %>%
#   mutate(Confidence=case_when(PerturbEffect.FTO <= -0.5 & PerturbEffect.METTL3i <= -0.5 ~ "high",
#                               PerturbEffect.FTO <= -0.1 & PerturbEffect.METTL3i <= -0.1 ~ "medium", .default = "low") %>% factor(levels=c("low","medium","high")))
# dt.GLORI.HEK293.confident[,.N,keyby=.(Confidence)]
# #    Confidence     N
# # 1:        low 34503
# # 2:     medium 12191
# # 3:       high  6452
dt.GLORI.HEK293.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/GLORI/dt.GLORI.HEK293.confident.m6As.RDS")

dt.benchmark.perturbation.validation <- rbind(dt.eTAMseq.mESC.confident %>% dplyr::filter(PerturbEffect<= -0.1) %>% dplyr::select(name,PerturbEffect) %>% mutate(PerturbGroup="eTAMseq_mESC.Mettl3.cKO"),
                                              dt.GLORI.HEK293.confident %>% dplyr::filter(PerturbEffect.FTO <= -0.1) %>% dplyr::select(name,PerturbEffect=PerturbEffect.FTO) %>% mutate(PerturbGroup="GLORI_HEK293.FTO.treated"),
                                              dt.GLORI.HEK293.confident %>% dplyr::filter(PerturbEffect.METTL3i <= -0.1) %>% dplyr::select(name,PerturbEffect=PerturbEffect.METTL3i) %>% mutate(PerturbGroup="GLORI_HEK293.METTL3.inhibitor")) %>%
  group_by(PerturbGroup) %>% mutate(Nsites=n_distinct(name)) %>% as.data.table() %>% mutate(PerturbationValidatedSite=paste0(PerturbGroup,"(Nsite=",Nsites,")"))

p.benchmark.perturbation.validation <- ggplot(dt.benchmark.perturbation.validation, aes(x = PerturbEffect, color=PerturbationValidatedSite)) +
  geom_density(alpha = 1, size = 0.5) +
  labs(x = "Delta of m6A level under perturbation",y = "Density")+
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  scale_color_manual(values = c("#a0d4ae","#0d8b43","#5562ab"), breaks = c("GLORI_HEK293.FTO.treated(Nsite=41566)", "GLORI_HEK293.METTL3.inhibitor(Nsite=30223)", "eTAMseq_mESC.Mettl3.cKO(Nsite=7870)"))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  guides(color=guide_legend(ncol = 1))+
  theme(
    legend.position = "top",legend.title = element_blank(),
    legend.key.height = unit(6,"pt"),legend.key.width = unit(6,"pt"),
    plot.title = element_text(face = "plain", size = rel(1.05), hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.9)),
    axis.text = element_text(size = rel(0.85)),
    panel.grid.major = element_line(colour = "grey92", size = 0.3),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = rel(0.7), hjust = 0)
  )

#figure 1d. Use antibody-independent m6A sites is better benchmarks than m6A peaks from IVT samples for performance estimation (IVT)
#MACS2 default combo (enrichment fold>=2) of WT and IVT of HEK_NEB, HEK_SYSY, and HEK_Abcam
dt.MACS2.default.peaks <- LoadPeakMethod(Method="MACS2",Method.peak.dir = "/data/m6A_calling_strategy/sm6APeak/MACS2/")
dt.MACS2.default.peaks <- dt.MACS2.default.peaks %>% dplyr::filter(seqnames %in% paste0("chr",1:22))
dt.MACS2.default.peaks[,.(NPeak=n_distinct(name), NPeak.FoldOver2=n_distinct(name[score>=2]),NPeak.FoldOver5=n_distinct(name[score>=5])),keyby=.(Sample)]
dt.MACS2.default.peaks <- dt.MACS2.default.peaks %>% dplyr::filter(score>=2) %>% mutate(SampleName=strsplit(Sample,split="/",fixed=T) %>% sapply("[",2)) %>%
  mutate(Condition=case_when(str_detect(Sample,"IVT")~"IVT",.default = "mRNA"), Antibody=case_when(str_detect(Sample,"SYSY") ~ "SYSY", str_detect(Sample,"Abcam") ~ "Abcam", .default = "NEB"),
         Cell=case_when(str_detect(Sample,"HEK") ~ "HEK", .default = "mESC"))
dt.sample.info <- dt.MACS2.default.peaks %>% distinct(Sample,SampleName,Condition,Antibody,Cell)
#obtain the IVT+ peaks (overlap >=50% in IVT1 or IVT2) and IVT- peaks for each sample
dt.MACS2.default.peaks.IVToverlap <- foreach(c = dt.sample.info$SampleName %>% grep(pattern="IVT",invert = T), .combine='rbind')%do%{
  foreach(ivt = dt.sample.info[Antibody==dt.sample.info$Antibody[c] & Cell==dt.sample.info$Cell[c] & Condition=="IVT",SampleName], .combine='rbind')%do%{
   bed.peak <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==dt.sample.info$SampleName[c]) %>% dplyr::select(seqnames,start,end,name,score,strand)
   bed.ivt <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==ivt) %>% dplyr::select(seqnames,start,end,name,score,strand)
   dt.overlap <- bt.intersect(a=bed.peak, b=bed.ivt, s=F, f=0.5,F=0.5, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
   dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
   dt.ivt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(IVToverlap=ifelse(sum.overlap/width>=0.5,"IVT+","IVT-"))
   bed.peak %>% left_join(x=.,y=dt.ivt.peak %>% distinct(name,IVToverlap) %>% mutate(IVTsample=ivt), by="name")
  } %>% mutate(SampleName=dt.sample.info$SampleName[c])
}
dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,IVTsample) %>%
  mutate(nPeak=n_distinct(name),IVT_pos=n_distinct(name[IVToverlap=="IVT+"]),IVT_neg=n_distinct(name[IVToverlap=="IVT-"])) %>% mutate(IVT_ratio=round(IVT_neg/IVT_pos,2)) %>%
  as.data.table() %>% distinct(SampleName,IVTsample,nPeak,IVT_pos,IVT_neg,IVT_ratio)
#         SampleName      IVTsample nPeak IVT_pos IVT_neg IVT_ratio
# 1: HEK_Abcam_mRNA1 HEK_Abcam_IVT1 13275     631   12644     20.04
# 2: HEK_Abcam_mRNA1 HEK_Abcam_IVT2 13275    3817    9458      2.48
# 3: HEK_Abcam_mRNA2 HEK_Abcam_IVT1 13072     263   12809     48.70
# 4: HEK_Abcam_mRNA2 HEK_Abcam_IVT2 13072    3710    9362      2.52
# 5:   HEK_NEB_mRNA1   HEK_NEB_IVT1 28358    5403   22955      4.25
# 6:   HEK_NEB_mRNA1   HEK_NEB_IVT2 28358    3776   24582      6.51
# 7:   HEK_NEB_mRNA2   HEK_NEB_IVT1 53547    6940   46607      6.72
# 8:   HEK_NEB_mRNA2   HEK_NEB_IVT2 53547    4703   48844     10.39
# 9:  HEK_SYSY_mRNA1  HEK_SYSY_IVT1 25231    4929   20302      4.12
# 10:  HEK_SYSY_mRNA1  HEK_SYSY_IVT2 25231    5650   19581      3.47
# 11:  HEK_SYSY_mRNA2  HEK_SYSY_IVT1 18443    3369   15074      4.47
# 12:  HEK_SYSY_mRNA2  HEK_SYSY_IVT2 18443    1183   17260     14.59
# 13:        mESC_WT1      mESC_IVT1 17212    1315   15897     12.09
# 14:        mESC_WT1      mESC_IVT2 17212    1274   15938     12.51
# 15:        mESC_WT2      mESC_IVT1 17690    1516   16174     10.67
# 16:        mESC_WT2      mESC_IVT2 17690    2021   15669      7.75
#obtain the overlap of benchmark sites for each peaks
#load HEK benchmark m6As
HEK293.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="HEK293",value=T)
names(HEK293.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="HEK293",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.HEK293 <-  foreach(f = 1:length(HEK293.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(HEK293.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}
#load confident HEK benchmark sites
dt.GLORI.HEK293.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/GLORI/dt.GLORI.HEK293.confident.m6As.RDS")
#load mESC benchmark m6As
mESC.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="mESC",value=T)
names(mESC.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="mESC",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.mESC <-  foreach(f = 1:length(mESC.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(mESC.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}
#load confident mESC benchmark sites
dt.eTAMseq.mESC.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/eTAMseq/dt.eTAMseq.mESC.confident.m6As.RDS")

dt.MACS2.default.peaks.Benchoverlap <- foreach(c = dt.sample.info$SampleName %>% grep(pattern="IVT",invert = T), .combine='rbind')%do%{
  bed.peak <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==dt.sample.info$SampleName[c]) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(dt.sample.info$SampleName[c],"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }else{
      bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(dt.sample.info$SampleName[c],"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(SampleName=dt.sample.info$SampleName[c])
}
#the count/ratio of benchmark+ peak among IVT+ and IVT- peaks
dt.MACS2.IVT.benchmark.association <- dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) #classify a peak as IVT- if not overlap with any IVT samples
dt.MACS2.IVT.benchmark.association <- dt.MACS2.IVT.benchmark.association %>% left_join(x=.,y=dt.MACS2.default.peaks.Benchoverlap,by=c("SampleName","name")) %>%
  group_by(SampleName,IVToverlapped) %>% mutate(TotalPeak=n_distinct(name)) %>%
  group_by(SampleName,IVToverlapped,TotalPeak,benchmark_name) %>% mutate(nPeak.bench=n_distinct(name[!is.na(nBenchSite) & nBenchSite>0]),SumBenchSites=sum(nBenchSite[!is.na(nBenchSite)])) %>%
  mutate(nPeak.bench.confident=n_distinct(name[!is.na(nConfidentBenchSite) & nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(SampleName,IVToverlapped,TotalPeak,benchmark_name,nPeak.bench,SumBenchSites,nPeak.bench.confident) %>% dplyr::filter(!is.na(benchmark_name))
# SampleName IVToverlapped TotalPeak benchmark_name nPeak.bench SumBenchSites nPeak.bench.confident
# 1: HEK_Abcam_mRNA1          IVT-      9300   GLORI_HEK293        7454         51754                  3626
# 2: HEK_Abcam_mRNA1          IVT+      3975   GLORI_HEK293        1816          9250                   505
# 3: HEK_Abcam_mRNA2          IVT-      9314   GLORI_HEK293        7390         55575                  3560
# 4: HEK_Abcam_mRNA2          IVT+      3758   GLORI_HEK293        1777          9354                   500
# 5:   HEK_NEB_mRNA1          IVT-     22231   GLORI_HEK293       12272         71170                  4592
# 6:   HEK_NEB_mRNA1          IVT+      6127   GLORI_HEK293        1898          6981                   471
# 7:   HEK_NEB_mRNA2          IVT-     45701   GLORI_HEK293       15525         85334                  5286
# 8:   HEK_NEB_mRNA2          IVT+      7846   GLORI_HEK293        2120          7350                   469
# 9:  HEK_SYSY_mRNA1          IVT-     18876   GLORI_HEK293       12455         80091                  4955
# 10:  HEK_SYSY_mRNA1          IVT+      6355   GLORI_HEK293        2137          9063                   615
# 11:  HEK_SYSY_mRNA2          IVT-     14909   GLORI_HEK293       10523         83192                  4478
# 12:  HEK_SYSY_mRNA2          IVT+      3534   GLORI_HEK293        1540          7043                   487
# 13:        mESC_WT1          IVT-     15352   eTAMseq_mESC        6766         31335                  2670
# 14:        mESC_WT1          IVT-     15352     GLORI_mESC        8570         47362                  2564
# 15:        mESC_WT1          IVT+      1860   eTAMseq_mESC         948          3478                   426
# 16:        mESC_WT1          IVT+      1860     GLORI_mESC        1161          4854                   421
# 17:        mESC_WT2          IVT-     15078   eTAMseq_mESC        6595         28741                  2612
# 18:        mESC_WT2          IVT-     15078     GLORI_mESC        8569         43965                  2506
# 19:        mESC_WT2          IVT+      2612     GLORI_mESC        1497          5603                   474
# 20:        mESC_WT2          IVT+      2612   eTAMseq_mESC        1226          4160                   480
#figure 1d the proportion of IVT+ (MACS2) peaks  among different cell x antibody (four sample)
dt.IVT.ratio.MACS2 <- dt.MACS2.IVT.benchmark.association %>% left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  distinct(SampleName,Antibody,Cell,IVToverlapped,TotalPeak) %>% pivot_wider(id_cols = c("SampleName","Antibody","Cell"), names_from = "IVToverlapped", values_from = "TotalPeak") %>%
  as.data.table() %>% mutate(PeakPerSample=`IVT-`+`IVT+`,IVTratio=`IVT+`/(`IVT-`+`IVT+`)) %>% group_by(Antibody,Cell) %>% mutate(avg.PeakPerSample=mean(PeakPerSample),avg.IVTratio=mean(IVTratio)) %>%
  distinct(Antibody,Cell,avg.PeakPerSample,avg.IVTratio) %>% mutate(GroupName=paste0(Antibody,"_",Cell)) %>% as.data.table()
dt.IVT.ratio.MACS2 <- dt.IVT.ratio.MACS2 %>% distinct(avg.PeakPerSample, avg.IVTratio, GroupName) %>%
  mutate(GroupName=str_replace(GroupName,pattern="HEK",replacement="HEK293"),avg.PeakPerSample=round(avg.PeakPerSample,0)) %>%
  mutate(label=paste0("NoPeak=",avg.PeakPerSample),label.pct=scales::percent(avg.IVTratio, accuracy = 0.1))
dt.IVT.ratio.MACS2 <- dt.IVT.ratio.MACS2 %>% dplyr::arrange(desc(avg.IVTratio)) %>% mutate(GroupName=factor(GroupName,levels=GroupName))
p.IVT.ratio.MACS2 <- ggplot(dt.IVT.ratio.MACS2,
                            aes(y = GroupName, x = avg.IVTratio, fill = GroupName)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  # percent label to the right of the bar
  geom_text(aes(label = label.pct, x =  avg.IVTratio + 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  #N peak label to the left of the bar
  geom_text(aes(label = label, x = 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  scale_fill_manual(values = c("#ffd47f","#f7c1cf","#7b93c6","#acd9ee"),breaks = c("Abcam_HEK293","SYSY_HEK293","NEB_HEK293","NEB_mESC")) +
  labs(x = "Proportion of IVT+ MACS2 peaks ", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

#compare the  IVT+ and IVT- peaks with random simulation peaks regarding containing GLORI sites
dt.MACS2.default.peaks.IVToverlap.type <- dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,score,SampleName,IVToverlapped)
dt.MACS2.default.peaks.IVToverlap.type[,.N,by=.(SampleName,IVToverlapped)]
#kept only genic peaks as input for simulation
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/m6A_seq/bin")
dt.MACS2.default.peaks.genic.annot <- rbind(annot_peak(dt.MACS2.default.peaks.IVToverlap %>% filter(!str_detect(name,"mESC")) %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                 strand = FALSE,fract = 0.5,gtf = "/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf", annot_type = c("gene", "feature", "both")[1]) %>%
                                              dplyr::select(seqnames,start,end,name,score,strand,OverlappedGeneIDs),
                                            annot_peak(dt.MACS2.default.peaks.IVToverlap %>% filter(str_detect(name,"mESC")) %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                       strand = FALSE,fract = 0.5,gtf = "/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf", annot_type = c("gene", "feature", "both")[1]) %>%
                                              dplyr::select(seqnames,start,end,name,score,strand,OverlappedGeneIDs)
                                            )

dt.MACS2.default.genic.peaks <- dt.MACS2.default.peaks.genic.annot %>% filter(!is.na(OverlappedGeneIDs))
dt.MACS2.default.genic.peaks <- dt.MACS2.default.genic.peaks %>% left_join(x=.,y=dt.MACS2.default.peaks.IVToverlap.type,by=c("name","score"))

#estimate cumulative curve for real IVT+ and IVT- peaks
dt.MACS2.default.genic.peaks <- dt.MACS2.default.genic.peaks %>% left_join(x=.,y=dt.MACS2.default.peaks.Benchoverlap %>% filter(str_detect(benchmark_name,"GLORI")) %>% distinct(name,benchmark_name,nBenchSite),by="name")
dt.MACS2.default.genic.peaks <- dt.MACS2.default.genic.peaks %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite))
dt.cum.data.MACS2.default.genic.peaks <- dt.MACS2.default.genic.peaks %>%
  # group_by(SampleName,IVToverlapped) %>%
  group_by(SampleName) %>%
  dplyr::arrange(desc(score)) %>%
  mutate(peak_width = end - start,
         cum_length = cumsum(peak_width),       # X-axis: Cumulative length
         cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
         ) %>% as.data.table() %>%
  mutate(IVToverlapped="all") %>% 
  distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites) %>% dplyr::arrange(SampleName,IVToverlapped,cum_length,cum_sites)

#write function to simulation and estimate cumulative curve
#use regionR package to generate random peaks for each samples IVT+ and IVT- peaks
library(regioneR)
library(GenomicRanges)
library(GenomicFeatures)

#simulation for HEK293 IVT+ and IVT- peaks
Txdb.hg38 <- makeTxDbFromGFF("~/genome_db/gencode.v44.annotation.gtf")
#extract all genic region to serve as the "allowed universe"
genes_universe <- genes(Txdb.hg38)
exons_universe <- exons(Txdb.hg38)
#use loop to generate 1000x for per sample per IVToverlapped condition and pick the 99%, 50% and 1% of simulation according to cumulative benchmark sites
dt.cum.data.MACS2.random.HEK.peaks <- foreach(S = unique(dt.MACS2.default.genic.peaks$SampleName) %>% grep(pattern="mESC",invert = T,value = T), .combine = 'rbind')%do%{
  # foreach(IVT = unique(dt.MACS2.default.genic.peaks$IVToverlapped), .combine = 'rbind')%do%{
  # real_peaks <- dt.MACS2.default.genic.peaks %>% filter(SampleName==S & IVToverlapped==IVT) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
  real_peaks <- dt.MACS2.default.genic.peaks %>% filter(SampleName==S) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
  real_peaks <- makeGRangesFromDataFrame(real_peaks,keep.extra.columns = T)
  #generate random peak
  # t1 <- Sys.time()
  nSimulation=100
  registerDoParallel(cl=10)
  dt.random.peaks <- foreach(R = 1:nSimulation, .combine = 'rbind')%dopar%{
    random_peaks <- randomizeRegions(real_peaks, genome = "hg38", allow.overlaps = FALSE, mask = NULL, 
                                     # universe = genes_universe)
                                     universe = exons_universe)
    #transform random peak into data table
    as.data.table(random_peaks) %>% mutate(name=paste("random",R,seqnames,start,end,strand,sep="_"),score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand)
  }
  stopImplicitCluster()
  # Sys.time() - t1
  #overlap with benchmark GLORI sites
  bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  dt.overlap <- bt.intersect(a=dt.random.peaks, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)#key parameter
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  # dt.overlap
  #benchmark sites with high confidence
  dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap <- dt.overlap %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2))
  #determine the 99%, 50% and 1% of simulation according to cumulative benchmark sites
  dt.overlap.cumsites <- dt.overlap %>% group_by(RandomID) %>% mutate(CumSites=sum(nBenchSite)) %>% as.data.table() %>%
    distinct(RandomID,CumSites) %>% dplyr::arrange(desc(CumSites)) 
  dt.randomid.selected <- dt.overlap.cumsites[c(0.01,0.5,0.99)*nSimulation,] %>% mutate(RandomRank=c("99%","50%","1%"))
  #output selected random peaks cum data for visualization
  dt.selected.overlap <- dt.overlap %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
  dt.selected.randompeaks <- dt.random.peaks %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2)) %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
  dt.selected.randompeaks <- dt.selected.randompeaks %>% left_join(x=.,y=dt.selected.overlap,by=c("name","RandomID")) %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite),
                                                                                                                                 nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))
  #cumulative output for selected random peaks
  dt.cum.data.random.peaks <- dt.selected.randompeaks %>% group_by(RandomID) %>% mutate(peak_width = end - start,
                                                                                        cum_length = cumsum(peak_width),       # X-axis: Cumulative length
                                                                                        cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
  ) %>% as.data.table() %>% distinct(RandomID,name,nBenchSite,cum_length,cum_sites) %>% left_join(x=.,y=dt.randomid.selected,by="RandomID")
  # message(paste0("Finished simulation of ", S, " ", IVT, " peaks"))
  message(paste0("Finished simulation of ", S, " ", " peaks"))
  rm(dt.random.peaks)
  rm(dt.overlap)
  # dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped=IVT)
  dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped="all")
  # }
}
dt.cum.data.MACS2.random.HEK.peaks <- dt.cum.data.MACS2.random.HEK.peaks %>% filter(IVToverlapped=="all") 
#simulation for mESC IVT+ and IVT- peaks
Txdb.mm39 <- makeTxDbFromGFF("~/genome_db/gencode.vM33.annotation.gtf")
#extract all genic region to serve as the "allowed universe"
genes_universe <- genes(Txdb.mm39)
exons_universe <- exons(Txdb.mm39)
#use loop to generate 1000x for per sample per IVToverlapped condition and pick the 99%, 50% and 1% of simulation according to cumulative benchmark sites
dt.cum.data.MACS2.random.mESC.peaks <- foreach(S = unique(dt.MACS2.default.genic.peaks$SampleName) %>% grep(pattern="mESC",invert = F,value = T), .combine = 'rbind')%do%{
  # foreach(IVT = unique(dt.MACS2.default.genic.peaks$IVToverlapped), .combine = 'rbind')%do%{
  # real_peaks <- dt.MACS2.default.genic.peaks %>% filter(SampleName==S & IVToverlapped==IVT) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
  real_peaks <- dt.MACS2.default.genic.peaks %>% filter(SampleName==S) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
  real_peaks <- makeGRangesFromDataFrame(real_peaks,keep.extra.columns = T)
  #generate random peak
  # t1 <- Sys.time()
  nSimulation=100
  registerDoParallel(cl=10)
  dt.random.peaks <- foreach(R = 1:nSimulation, .combine = 'rbind')%dopar%{
    random_peaks <- randomizeRegions(real_peaks, genome = "mm39", allow.overlaps = FALSE, mask = NULL, 
                                     # universe = genes_universe)
                                     universe = exons_universe)
    #transform random peak into data table
    as.data.table(random_peaks) %>% mutate(name=paste("random",R,seqnames,start,end,strand,sep="_"),score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand)
  }
  stopImplicitCluster()
  # Sys.time() - t1
  #overlap with benchmark GLORI sites
  bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::filter(benchmark_name=="GLORI_mESC") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  dt.overlap <- bt.intersect(a=dt.random.peaks, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)##key parameter 
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  dt.overlap
  #benchmark sites with high confidence
  dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap <- dt.overlap %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2))
  #determine the 99%, 50% and 1% of simulation according to cumulative benchmark sites
  dt.overlap.cumsites <- dt.overlap %>% group_by(RandomID) %>% mutate(CumSites=sum(nBenchSite)) %>% as.data.table() %>%
    distinct(RandomID,CumSites) %>% dplyr::arrange(desc(CumSites)) 
  dt.randomid.selected <- dt.overlap.cumsites[c(0.01,0.5,0.99)*nSimulation,] %>% mutate(RandomRank=c("99%","50%","1%"))
  #output selected random peaks cum data for visualization
  dt.selected.overlap <- dt.overlap %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
  dt.selected.randompeaks <- dt.random.peaks %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2)) %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
  dt.selected.randompeaks <- dt.selected.randompeaks %>% left_join(x=.,y=dt.selected.overlap,by=c("name","RandomID")) %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite),
                                                                                                                                 nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))
  #cumulative output for selected random peaks
  dt.cum.data.random.peaks <- dt.selected.randompeaks %>% group_by(RandomID) %>% mutate(peak_width = end - start,
                                                                                        cum_length = cumsum(peak_width),       # X-axis: Cumulative length
                                                                                        cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
  ) %>%
    as.data.table() %>% distinct(RandomID,name,nBenchSite,cum_length,cum_sites) %>% left_join(x=.,y=dt.randomid.selected,by="RandomID")
  # message(paste0("Finished simulation of ", S, " ", IVT, " peaks"))
  message(paste0("Finished simulation of ", S, " ", " peaks"))
  rm(dt.random.peaks)
  rm(dt.overlap)
  # dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped=IVT)
  dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped="all")
  # }
}
dt.cum.data.MACS2.random.mESC.peaks <- dt.cum.data.MACS2.random.mESC.peaks %>% filter(IVToverlapped=="all") 

#plot cumulative curve of MACS2 IVT+/ IVT- peaks and simulated peaks 
dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks <- rbind(dt.cum.data.MACS2.default.genic.peaks %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites) %>% mutate(PeakGroup="MACS2"),
                                                                 dt.cum.data.MACS2.random.HEK.peaks %>% mutate(PeakGroup=paste0("Random",RandomRank)) %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites,PeakGroup),
                                                                 dt.cum.data.MACS2.random.mESC.peaks %>% mutate(PeakGroup=paste0("Random",RandomRank)) %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites,PeakGroup)) %>%
  dplyr::arrange(SampleName,IVToverlapped,PeakGroup,cum_length,cum_sites)
dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>% slice_max(cum_length) %>% as.data.table() %>% dplyr::arrange(SampleName,IVToverlapped,PeakGroup)
#scaled cum_length and cum_sites for each Sample x IVToverlapped:use the median cum_sites as 100%
#pruned the extreme long cum_length x cum_sites (keep only 100 data points)
dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks[,.N,by=.(SampleName,IVToverlapped,PeakGroup)]
dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks <- foreach(S = unique(dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks$SampleName),.combine = 'rbind')%do%{
  foreach(IVT = unique(dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks$IVToverlapped),.combine = 'rbind')%do%{
    dt.cumdata <- dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% dplyr::filter(SampleName==S & IVToverlapped==IVT)
    foreach(G = unique(dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks$PeakGroup), .combine = 'rbind')%do%{
      dt.cumdata.subgroup <- dt.cumdata %>% dplyr::filter(PeakGroup==G) %>% dplyr::arrange(cum_length,cum_sites)
      dt.cumdata.subgroup[floor(nrow(dt.cumdata.subgroup)*seq(0,1,0.001)[-1]),] %>% mutate(scaled_cum_peaks=seq(0,1,0.001)[-1]*100)
    }
  }
}
#scaled the cum_sites, and cum_peaks
dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks <- dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% filter(PeakGroup %in% c("MACS2","Random50%","Random99%")) %>%
  group_by(SampleName,IVToverlapped,PeakGroup) %>% mutate(scaled_cum_sites=cum_sites/max(cum_sites)*100) %>% as.data.table()
#average of sample group
dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks <- dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% dplyr::filter(IVToverlapped=="all") %>%
  mutate(Sample=case_when(str_detect(SampleName,"HEK_Abcam") ~ "Abcam_HEK293",
                          str_detect(SampleName,"HEK_NEB") ~ "NEB_HEK293",
                          str_detect(SampleName,"HEK_SYSY") ~ "SYSY_HEK293",
                          str_detect(SampleName,"mESC_WT") ~ "NEB_mESC")) %>%
  group_by(Sample,PeakGroup,scaled_cum_peaks) %>% mutate(scaled_cum_peaks=mean(scaled_cum_peaks) %>% round(digits = 2),scaled_cum_sites=mean(scaled_cum_sites) %>% round(digits = 2)) %>%
  as.data.table() %>% distinct(Sample,IVToverlapped,PeakGroup,scaled_cum_peaks,scaled_cum_sites) %>%
  mutate(PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MACS2")))

p.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks <- ggplot(dt.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% dplyr::filter(IVToverlapped=="all"),
                                                                        aes(x = scaled_cum_peaks, y = scaled_cum_sites)) +
  geom_line(aes(linetype=PeakGroup,color=PeakGroup),linewidth=0.5) +
  geom_abline(intercept = 0,slope = 1,color="grey50",linetype="solid",linewidth=0.25)+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  scale_linetype_manual(values=c("solid","dashed","dotted"),breaks = c("MACS2","Random99%","Random50%"))+
  # scale_y_continuous(transform = "log10")+
  # geom_hline(yintercept = c(1000),color="grey50",linetype="solid",linewidth=0.25)+
  labs(
    title = NULL,
    x = "Scaled cumulative peak proportion (%)",
    y = "Scaled cumulative recovered GLORI sites (%)",
  ) +
  guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+
  # coord_cartesian(xlim=c(0,2.5))+
  facet_wrap(~Sample,ncol = 2)+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  guides(
    colour = guide_legend(
      override.aes = list(
        linetype = c("dotted","dashed","solid")   # ←←← 根据你实际的线型顺序修改
      ),
      title = "PeakGroup",
      nrow = 1,                    # 可选：横向排列
      keywidth = 1.5
    ),
    linetype = "none"              # 隐藏独立的 linetype legend
  )
p.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks

#the TPR  and FPR comparison of top25% v.s. tail25%
dt.top25vstail25.MACS2.default.genic.peaks <-  rbind(
  dt.MACS2.default.genic.peaks %>% group_by(SampleName) %>% dplyr::arrange(desc(score)) %>% mutate(totalsites=sum(nBenchSite)) %>% slice_head(prop = 0.25) %>%
    mutate(cum_sites_ratio=sum(nBenchSite)/totalsites*100,FP_peak_ratio=n_distinct(name[nBenchSite==0])/n_distinct(name)*100,Group="top25%") %>%
    as.data.table() %>% distinct(SampleName,Group,cum_sites_ratio,FP_peak_ratio),
  dt.MACS2.default.genic.peaks %>% group_by(SampleName) %>% dplyr::arrange(desc(score)) %>% mutate(totalsites=sum(nBenchSite)) %>% slice_tail(prop = 0.25) %>%
    mutate(cum_sites_ratio=sum(nBenchSite)/totalsites*100,FP_peak_ratio=n_distinct(name[nBenchSite==0])/n_distinct(name)*100,Group="tail25%") %>%
    as.data.table() %>% distinct(SampleName,Group,cum_sites_ratio,FP_peak_ratio)) %>%
  dplyr::arrange(SampleName,Group)
dt.top25vstail25.MACS2.default.genic.peaks 
dt.top25vstail25.MACS2.default.genic.peaks <- dt.top25vstail25.MACS2.default.genic.peaks %>% mutate(Group=factor(Group,levels=c("tail25%","top25%")))
#pair-wise baxplot visualize TPR  and FPR comparison of top25% v.s. tail25%
p.top25vstail25.MACS2.cum.sites.ratio <- ggplot(dt.top25vstail25.MACS2.default.genic.peaks, aes(x = Group, y = cum_sites_ratio, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values = c("#F3D28A", "#F39030"),breaks = c("tail25%","top25%")) +
  labs(y = str_wrap("Cumulative proportion (%) of benchmark m6As",width = 40), x = "MACS2 peak subset") +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.cum.sites.ratio <- dt.top25vstail25.MACS2.default.genic.peaks %>%
  group_by(Group) %>%
  summarise(mean_val = round(mean(cum_sites_ratio), 2))
p.top25vstail25.MACS2.cum.sites.ratio <- p.top25vstail25.MACS2.cum.sites.ratio +
  geom_text(data = mean_values.cum.sites.ratio,
            aes(x = Group, 
                y = max(dt.top25vstail25.MACS2.default.genic.peaks$cum_sites_ratio) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.top25vstail25.MACS2.cum.sites.ratio <- p.top25vstail25.MACS2.cum.sites.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.top25vstail25.MACS2.default.genic.peaks$cum_sites_ratio) * 1.05,
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.top25vstail25.MACS2.default.genic.peaks$cum_sites_ratio) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.top25vstail25.MACS2.cum.sites.ratio
p.top25vstail25.MACS2.FP.peak.ratio <- ggplot(dt.top25vstail25.MACS2.default.genic.peaks, aes(x = Group, y = FP_peak_ratio, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values = c("#F3D28A", "#F39030"),breaks = c("tail25%","top25%")) +
  labs(y = str_wrap("Proportion (%) of peaks contain no benchmark m6As",width=45), x = "MACS2 peak subset") +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.FP.peak.ratio <- dt.top25vstail25.MACS2.default.genic.peaks %>%
  group_by(Group) %>%
  summarise(mean_val = round(mean(FP_peak_ratio), 2))
p.top25vstail25.MACS2.FP.peak.ratio <- p.top25vstail25.MACS2.FP.peak.ratio +
  geom_text(data = mean_values.FP.peak.ratio,
            aes(x = Group, 
                y = max(dt.top25vstail25.MACS2.default.genic.peaks$FP_peak_ratio) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.top25vstail25.MACS2.FP.peak.ratio <- p.top25vstail25.MACS2.FP.peak.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.top25vstail25.MACS2.default.genic.peaks$FP_peak_ratio) * 1.05,
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.top25vstail25.MACS2.default.genic.peaks$FP_peak_ratio) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.top25vstail25.MACS2.FP.peak.ratio
cowplot::plot_grid(p.top25vstail25.MACS2.cum.sites.ratio,p.top25vstail25.MACS2.FP.peak.ratio) 
# save.image("updated_sm6APeak_Figure1.intermediate.results.RDS")
#barplot to display the max cumulative GLORI sites of MACS2 IVT+/IVT- peaks vs RandomPeaks (50% and 99%)
dt.density.GLORI.MACS2vsRandom <- dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>%
  mutate(TotalPeak=n(),CumWidth=max(cum_length)/1e6,CumSites=max(cum_sites)) %>% as.data.table() %>% distinct(SampleName,IVToverlapped,PeakGroup,TotalPeak,CumWidth,CumSites) %>%
  dplyr::arrange(SampleName,IVToverlapped,PeakGroup,TotalPeak,CumWidth,CumSites) %>% mutate(density.perMb=CumSites/CumWidth, density.per1000Peak=CumSites/TotalPeak*1000)
dt.density.GLORI.MACS2vsRandom <- dt.density.GLORI.MACS2vsRandom %>% dplyr::filter(PeakGroup !="Random1%") %>%
  mutate(SampleName=factor(SampleName),PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MACS2")))
#keep all peak 
p.density.GLORI.MACS2vsRandom.per1000Peak  <- ggplot(data=dt.density.GLORI.MACS2vsRandom%>% dplyr::filter(IVToverlapped=="all"),
                                                       aes(x=SampleName,y=density.per1000Peak, fill=PeakGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  labs(x=NULL,y="Density of GLORI m6A sites (per 1000 peak)")+
  guides(fill=guide_legend(title = "PeakGroup"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),transform = "log10")+
  # facet_wrap(~IVToverlapped,scales = "fixed")+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.density.GLORI.MACS2vsRandom.per1000Peak
#add enrich fold v.s Random50%
dt.density.fold.GLORI.MACS2vsRandom.per1000Peak <- dt.density.GLORI.MACS2vsRandom%>% dplyr::filter(IVToverlapped=="all") %>%
  group_by(SampleName) %>% mutate(DensityFold=density.per1000Peak[PeakGroup=="MACS2"]/density.per1000Peak[PeakGroup=="Random50%"]) %>%
  as.data.table() %>% distinct(SampleName,DensityFold)
p.density.GLORI.MACS2vsRandom.per1000Peak <- p.density.GLORI.MACS2vsRandom.per1000Peak+
  geom_text(data = dt.density.fold.GLORI.MACS2vsRandom.per1000Peak,
            inherit.aes = FALSE,
            aes(x = SampleName, 
                y = max(dt.density.GLORI.MACS2vsRandom$density.per1000Peak) * 1.05,
                label = paste(round(DensityFold),"x")),
            size = 2.5, fontface = "bold", color = "#EA9E58")
p.density.GLORI.MACS2vsRandom.per1000Peak

#compare the cumulative benchmark m6As of MACS2 vs. Random Peaks
dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks <- dt.cum.data.MACS2.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>%
  mutate(FPR=round(n_distinct(name[nBenchSite==0])/n_distinct(name)*100,2),TotalBenchSites=max(cum_sites)) %>%
  mutate(TPR=dplyr::case_when(str_detect(SampleName,"HEK") ~ round(TotalBenchSites/nrow(dt.benchmark.m6A.HEK293[benchmark_name=="GLORI_HEK293",])*100,2),
                              str_detect(SampleName,"mESC") ~ round(TotalBenchSites/nrow(dt.benchmark.m6A.mESC[benchmark_name=="GLORI_mESC",])*100,2))) %>% 
  as.data.table() %>% distinct(SampleName,IVToverlapped,PeakGroup,FPR,TPR)
dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks <- dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks %>% dplyr::filter(PeakGroup %in% c("MACS2","Random50%","Random99%")) %>%
  mutate(PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MACS2")))
p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- ggplot(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks, aes(x = PeakGroup, y = TPR, fill = PeakGroup)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color=PeakGroup),size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  scale_y_continuous(expand = expansion(add = c(0,0.3)),breaks = seq(0,max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$TPR),by=10))+
  labs(y = str_wrap("Cumulative proportion (%) of benchmark m6As",width = 40), x = NULL) +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.cum.sites.ratio <- dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks %>%
  group_by(PeakGroup) %>%
  summarise(mean_val = round(mean(TPR), 2))
p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio +
  geom_text(data = mean_values.cum.sites.ratio,
            aes(x = PeakGroup, 
                y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$TPR) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     comparisons = list(c("MACS2","Random50%")),
                     bracket.size = 0.1,
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$TPR) * 0.8,
                     label.y.npc = "bottom",
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$TPR) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio

#compare the ratio of peak contain benchmark m6As of MACS2 vs Random Peaks
p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- ggplot(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks, aes(x = PeakGroup, y = FPR, fill = PeakGroup)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color=PeakGroup),size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MACS2"))+
  scale_y_continuous(expand = expansion(add = c(0,0.3)),breaks = seq(0,max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$FPR),by=10))+
  labs(y = str_wrap("Proportion (%) of peaks contain no benchmark m6As",width = 40), x = NULL) +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.FP.peak.ratio<- dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks %>%
  group_by(PeakGroup) %>%
  summarise(mean_val = round(mean(FPR), 2))
p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio +
  geom_text(data = mean_values.FP.peak.ratio,
            aes(x = PeakGroup, 
                y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$FPR) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     comparisons = list(c("MACS2","Random50%")),
                     bracket.size = 0.1,
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$FPR) * 0.8,
                     label.y.npc = "bottom",
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.TPR.FPR.MACS2.default.genic.peaks.vs.random.peaks$FPR) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio
save.image("updated_sm6APeak_Figure1.intermediate.results.RDS")


#figure 1e the proportion of (MACS2) peaks contains benchmark sites among IVT+ and IVT- peak among different cellxantibody
dt.benchmark.ratio.MACS2 <- dt.MACS2.IVT.benchmark.association %>% mutate(benchmark_ratio=nPeak.bench/TotalPeak) %>%
  left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  group_by(Antibody,Cell,benchmark_name,IVToverlapped) %>% mutate(avg.benchmark_ratio=mean(benchmark_ratio)) %>% as.data.table() %>%
  distinct(Antibody,Cell,benchmark_name,IVToverlapped,avg.benchmark_ratio) %>% dplyr::arrange(Antibody,Cell,benchmark_name,IVToverlapped)
dt.benchmark.ratio.MACS2 <- dt.benchmark.ratio.MACS2 %>% mutate(GroupName=paste0(Antibody,"_",Cell," | ",strsplit(benchmark_name,split="_") %>% sapply("[",1))) %>%
  dplyr::arrange(IVToverlapped,desc(avg.benchmark_ratio)) %>% mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>% dplyr::arrange(IVToverlapped,GroupName)
ttest.p.benchmark.ratio.MACS2 <- signif(t.test(x=dt.benchmark.ratio.MACS2[IVToverlapped=="IVT+",avg.benchmark_ratio],y=dt.benchmark.ratio.MACS2[IVToverlapped=="IVT-",avg.benchmark_ratio],paired = T )$p.value,2)
p.benchmark.ratio.MACS2 <- ggplot(data=dt.benchmark.ratio.MACS2,aes(x=GroupName, y=avg.benchmark_ratio,fill=IVToverlapped))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#fcbb46","#ef756c"),breaks=c("IVT+", "IVT-"))+
  labs(x=NULL,y="% of peaks harbor benchmark m6A sites", subtitle = paste0("paired t-test p-value=", ttest.p.benchmark.ratio.MACS2))+
  guides(fill=guide_legend(title = "MACS2 peaks",nrow = 1))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.9,0.8),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.benchmark.ratio.MACS2
#figure 1f the proportion of (MACS2) peaks contains high confident benchmark sites (perturbation experiment validated) among IVT+ and IVT- peak among different cell x antibody
dt.confident.benchmark.ratio.MACS2 <- dt.MACS2.IVT.benchmark.association %>% mutate(confident.benchmark_ratio=nPeak.bench.confident/TotalPeak) %>%
  left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  group_by(Antibody,Cell,benchmark_name,IVToverlapped) %>% mutate(avg.confident.benchmark_ratio=mean(confident.benchmark_ratio)) %>% as.data.table() %>%
  distinct(Antibody,Cell,benchmark_name,IVToverlapped,avg.confident.benchmark_ratio) %>% dplyr::arrange(Antibody,Cell,benchmark_name,IVToverlapped)
dt.confident.benchmark.ratio.MACS2 <- dt.confident.benchmark.ratio.MACS2 %>% mutate(GroupName=paste0(Antibody,"_",Cell," | ",strsplit(benchmark_name,split="_") %>% sapply("[",1))) %>%
  dplyr::arrange(IVToverlapped,desc(avg.confident.benchmark_ratio)) %>% mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>% dplyr::arrange(IVToverlapped,GroupName)
ttest.p.confident.benchmark.ratio.MACS2 <- signif(t.test(x=dt.confident.benchmark.ratio.MACS2[IVToverlapped=="IVT+",avg.confident.benchmark_ratio],
                                                         y=dt.confident.benchmark.ratio.MACS2[IVToverlapped=="IVT-",avg.confident.benchmark_ratio],paired = T )$p.value,2)
p.confident.benchmark.ratio.MACS2 <- ggplot(data=dt.confident.benchmark.ratio.MACS2,aes(x=GroupName, y=avg.confident.benchmark_ratio,fill=IVToverlapped))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#fcbb46","#ef756c"),breaks=c("IVT+", "IVT-"))+
  labs(x=NULL,y=str_wrap("% of peaks harbor confident benchmark m6A sites",width = 40), subtitle = paste0("paired t-test p-value=", ttest.p.confident.benchmark.ratio.MACS2))+
  guides(fill=guide_legend(title = "MACS2 peaks",nrow = 1))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.9,0.8),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.confident.benchmark.ratio.MACS2
#figure 1e the count of m6A+ genes (MACS2) that could be missed if use IVT as negative control
#HEK293
#annotate the MACS2 peak of HEK293 to gene
dt.MACS2.default.peaks.IVToverlap.HEK293 <- dt.MACS2.default.peaks.IVToverlap %>% dplyr::filter(str_detect(SampleName,"HEK")) %>% distinct(seqnames, start,end,name,score, strand)
dt.MACS2.default.peaks.IVToverlap.HEK293.annotgene <- annot_peak(peak.bed=dt.MACS2.default.peaks.IVToverlap.HEK293,
                                                                 strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                 annot_type = "gene")
#determine the IVT+ specific m6A+ genes
dt.MACS2.default.peaks.IVToverlap.HEK293.annotgene <- dt.MACS2.default.peaks.IVToverlap.HEK293.annotgene %>%
  left_join(x=.,y=dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
              as.data.table() %>% distinct(name,SampleName,IVToverlapped), by=c("name"))
dt.MACS2.IVTpos.specific.m6Agene.HEK293 <-   dt.MACS2.default.peaks.IVToverlap.HEK293.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName,OverlappedGenes) %>% mutate(IVTpos.specific=all(IVToverlapped=="IVT+")) %>% as.data.table() %>%
  distinct(SampleName,OverlappedGenes,IVTpos.specific) %>% dplyr::filter(IVTpos.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)
dt.MACS2.IVTpos.specific.m6Agene.HEK293 <- dt.MACS2.IVTpos.specific.m6Agene.HEK293 %>%
  left_join(x=.,y=dt.MACS2.default.peaks.IVToverlap.HEK293.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ","), by=c("SampleName","OverlappedGenes"))
#obtain the IVT+ peaks overlapped with (high confident) benchmark sites
dt.benchmark.MACS2.peak.HEK293 <-  dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) %>% dplyr::filter(IVToverlapped=="IVT+") %>% dplyr::filter(str_detect(SampleName,"HEK")) %>%
  inner_join(x=.,y=dt.MACS2.default.peaks.Benchoverlap %>% dplyr::filter(nBenchSite>0),by=c("name","SampleName"))
#collect IVT+ specific m6A genes with benchmark m6A sites support
dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.HEK293 %>%
  inner_join(x=.,y=dt.benchmark.MACS2.peak.HEK293 %>% dplyr::select(name,benchmark_name,nBenchSite, nConfidentBenchSite),by="name")
#rank the IVT+ specific m6A gene with benchmark m6A sites support
dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.sample.info %>% inner_join(x=.,y=dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark,by="SampleName")
dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% group_by(Antibody,Cell,OverlappedGenes) %>%
  mutate(nSampleName=n_distinct(SampleName)) %>% as.data.table() %>% dplyr::arrange(Antibody,Cell,desc(nSampleName),desc(nConfidentBenchSite))
#obtain final IVT+ specific m6A+ gene list with benchmark sites
dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% distinct(Antibody,Cell,nSampleName,OverlappedGenes,benchmark_name,nBenchSite,nConfidentBenchSite)
#meanwhile,obtain final IVT+ specific m6A+ gene list with confident benchmark sites
dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% dplyr::filter(nConfidentBenchSite>0)

#mESC
#annotate the MACS2 peak of mESC to gene
dt.MACS2.default.peaks.IVToverlap.mESC <- dt.MACS2.default.peaks.IVToverlap %>% dplyr::filter(str_detect(SampleName,"mESC")) %>% distinct(seqnames, start,end,name,score, strand)
dt.MACS2.default.peaks.IVToverlap.mESC.annotgene <- annot_peak(peak.bed=dt.MACS2.default.peaks.IVToverlap.mESC,
                                                                 strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                                                 annot_type = "gene")
dt.MACS2.default.peaks.IVToverlap.mESC.annotgene <- dt.MACS2.default.peaks.IVToverlap.mESC.annotgene %>%
  left_join(x=.,y=dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
              as.data.table() %>% distinct(name,SampleName,IVToverlapped), by=c("name"))
#determine the IVT+ specific m6A+ genes
dt.MACS2.IVTpos.specific.m6Agene.mESC <-   dt.MACS2.default.peaks.IVToverlap.mESC.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName,OverlappedGenes) %>% mutate(IVTpos.specific=all(IVToverlapped=="IVT+")) %>% as.data.table() %>%
  distinct(SampleName,OverlappedGenes,IVTpos.specific) %>% dplyr::filter(IVTpos.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)
dt.MACS2.IVTpos.specific.m6Agene.mESC <- dt.MACS2.IVTpos.specific.m6Agene.mESC %>%
  left_join(x=.,y=dt.MACS2.default.peaks.IVToverlap.mESC.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ","), by=c("SampleName","OverlappedGenes"))
#obtain the IVT+ peaks overlapped with (high confident) benchmark sites
dt.benchmark.MACS2.peak.mESC <-  dt.MACS2.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) %>% dplyr::filter(IVToverlapped=="IVT+") %>% dplyr::filter(str_detect(SampleName,"mESC")) %>%
  inner_join(x=.,y=dt.MACS2.default.peaks.Benchoverlap %>% dplyr::filter(nBenchSite>0),by=c("name","SampleName"))
#collect IVT+ specific m6A genes with benchmark m6A sites support
dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.mESC %>%
  inner_join(x=.,y=dt.benchmark.MACS2.peak.mESC %>% dplyr::select(name,benchmark_name,nBenchSite, nConfidentBenchSite),by="name")
#rank the IVT+ specific m6A gene with benchmark m6A sites support
dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.sample.info %>% inner_join(x=.,y=dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark,by="SampleName")
dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark %>% group_by(Antibody,Cell,OverlappedGenes) %>%
  mutate(nSampleName=n_distinct(SampleName)) %>% as.data.table() %>% dplyr::arrange(Antibody,Cell,desc(nSampleName),desc(nConfidentBenchSite))
#obtain final IVT+ specific m6A+ gene list with benchmark sites
dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark %>% distinct(Antibody,Cell,nSampleName,OverlappedGenes,benchmark_name,nBenchSite,nConfidentBenchSite)
#meanwhile,obtain final IVT+ specific m6A+ gene list with confident benchmark sites
dt.MACS2.IVTpos.specific.m6Agene.mESC.with.confident.benchmark <- dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark %>% dplyr::filter(nConfidentBenchSite>0)

#calculate the count of IVTpos specific m6A+ gene in each AntibodyxCell
dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark <- rbind(rbind(dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark[nBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)],
                                                               dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark[nBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)]) %>% mutate(Group="Contain benchmark m6A sites"),
                                                               rbind(dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.benchmark[nConfidentBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)],
                                                                     dt.MACS2.IVTpos.specific.m6Agene.mESC.with.benchmark[nConfidentBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)]) %>% mutate(Group="Contain confident benchmark m6A sites")
                                                               )
dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark <- dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark %>% mutate(SampleLabel=paste0(Antibody,"_",str_replace(Cell,"HEK","HEK293"))) %>%
  dplyr::arrange(Group,desc(nGene)) %>% dplyr::mutate(SampleLabel=factor(SampleLabel,levels=unique(SampleLabel)),
                                                      Group=factor(Group, levels=c("Contain benchmark m6A sites","Contain confident benchmark m6A sites"), labels=c("Benchmark m6As","Confident benchmark m6As")))
dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark <- dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark %>% mutate(label.count=paste0(paste0("No.Gene=",nGene)))
p.count.MACS2.IVTpos.specific.m6Agene.with.benchmark  <- ggplot(data=dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark, aes(x=SampleLabel,y=nGene, fill=Group))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark %>% dplyr::filter(Group=="Benchmark m6As"),
            aes(label = label.count, y =  nGene/2, x=SampleLabel),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.MACS2.IVTpos.specific.m6Agene.with.benchmark %>% dplyr::filter(Group=="Confident benchmark m6As"),
            aes(label = label.count, y =  nGene/2, x=SampleLabel),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y="N of genes harbor only IVT+ MACS2 peaks")+
  guides(fill=guide_legend(title = "IVT+ peaks contain"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.9,0.8),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.MACS2.IVTpos.specific.m6Agene.with.benchmark
#figure 1f-g. the IGV plot of selected m6A+ genes harbor IVT+ MACS2 peaks that support by confident m6A sites in HEK293 and mESC
dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)


#HEK293 IVT+ specific m6A+ gene with high confident benchmark m6As
selected.m6Agene.SYSY.HEK <-  dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="SYSY" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=5)
# OverlappedGenes nSampleName nConfidentBenchSite
# 1:           RPRD2           2                  12
# 2:           NUMA1           2                   9
# 3:            RRP1           2                   7
# 4:         NECTIN2           2                   6
# 5:          RALGDS           2                   6
selected.m6Agene.NEB.HEK <-  dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="NEB" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=5)
# OverlappedGenes nSampleName nConfidentBenchSite
# 1:       LINC02976           2                   6
# 2:           NR2F2           2                   5
# 3:           HSH2D           2                   4
# 4:           RBM23           2                   4
# 5:      VPS9D1-AS1           2                   4
selected.m6Agene.NEB.Abcam <-  dt.MACS2.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="Abcam" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=5)
# OverlappedGenes nSampleName nConfidentBenchSite
# 1:           RPRD2           2                  12
# 2:            ARF6           2                   8
# 3:           USP10           2                   7
# 4:            PDP2           2                   7
# 5:           AKAP1           2                   6
#select well-studied genes
dt.selected.gene.HEK.SYSY.MACS2 <- dt.gene.hg38 %>% filter(gene_name=="NUMA1")
dt.selected.gene <- dt.selected.gene.HEK.SYSY.MACS2
dt.selected.gene.HEK.NEB.MACS2 <- dt.gene.hg38 %>% filter(gene_name=="NR2F2")
dt.selected.gene <- dt.selected.gene.HEK.NEB.MACS2
dt.selected.gene.HEK.Abcam.MACS2 <- dt.gene.hg38 %>% filter(gene_name=="ARF6")
dt.selected.gene <- dt.selected.gene.HEK.Abcam.MACS2
#track for bigwig files of study samples
# StudySamples <- c("HEK_SYSY_mRNA1", "HEK_SYSY_mRNA2","HEK_SYSY_IVT1","HEK_SYSY_IVT2")
# StudySamples <- c("HEK_NEB_mRNA1", "HEK_NEB_mRNA2","HEK_NEB_IVT1","HEK_NEB_IVT2")
StudySamples <- c("HEK_Abcam_mRNA1", "HEK_Abcam_mRNA2","HEK_Abcam_IVT1","HEK_Abcam_IVT2")
      ## Create a page (7.5*7.5cm)
      window.size=1#100kb
      pseudocount=1
      track.height=0.5
      # pdf("Fig1_SYSY_HEK_NUMA1_IVTpos_specific_peak.pdf")
      # pdf("Fig1_NEB_HEK_NR2F2_IVTpos_specific_peak.pdf")
      pdf("Fig1_Abcam_HEK_ARF6_IVTpos_specific_peak.pdf")

      pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
      region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
      #plot gene
      paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38", width=3.0)
      ## Plot genes small
      genesPlot <- plotGenes(params = paramsgene,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                             x = 0.5, y = 0, height = 1.0,width = 5,just = c("left", "top"), default.units = "cm",fontsize = 7)
      annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

      if(dt.selected.gene$strand=="-"){
        dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
          Input.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
          RIP.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
          Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
        }
        #determine the signal range for pos and neg strand
        # range_both <- ceiling(quantile(abs(dt.density$score),1))+1
        range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
          distinct(sample,range) %>% as.data.table()
        range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
          group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
        #plot
        for(s in 1:length(StudySamples)){
          signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                   params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                   just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                     params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          #add MACS2 peak in this region
          dt.peak.region <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
          if(nrow(dt.peak.region)>0){
            for(r in 1:nrow(dt.peak.region)){
              plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                          fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
            }
          }
          plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
          sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
          sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
        }
      }else{
        dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
          Input.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
          RIP.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
          Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
        }
        #determine the signal range for pos strand
        range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
          distinct(sample,range) %>% as.data.table()
        range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
          group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
        #plot
        for(s in 1:length(StudySamples)){
          signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                   params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                   just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                     params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          #add MACS2 peak in this region
          dt.peak.region <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
          if(nrow(dt.peak.region)>0){
            for(r in 1:nrow(dt.peak.region)){
              plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                          fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
            }
          }
          plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
          sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
          sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
        }
      }
      #track for benchmark m6A sites (perturbation effect)
      #FTO GLORI HEK293
      signal_confident.benchmark.FTO <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.FTO), params = region,
                                                   x = 0.5, y = 1.5+length(StudySamples)*track.height, width = 5, height = track.height,
                                                   just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
      plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
      sample_label <- plotText(label="FTO treatment", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
      sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
      #METTL3i GLORI HEK293
      signal_confident.benchmark.METTL3i <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.METTL3i), params = region,
                                                   x = 0.5, y = 1.5+(length(StudySamples)+1)*track.height, width = 5, height = track.height,
                                                   just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
      plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+1.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
      sample_label <- plotText(label="METTL3 inhibition", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+1.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
      sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

      #add legend for MeRIP
      legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                                     just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                     title = expression("log"["2"]~"(CPM)"),fondface="bold")
      legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in GLORI_HEK293"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(length(StudySamples)+2)*track.height,
                                                   width = 2.5, height = 1,
                                          just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fondface="bold")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")
      dev.off()

#mESC IVT+ specific m6A+ gene with high confident benchmark m6As
selected.m6Agene.mESC.NEB <-  dt.MACS2.IVTpos.specific.m6Agene.mESC.with.confident.benchmark %>% dplyr::filter(Antibody=="NEB" & Cell=="mESC" & nSampleName==2) %>%
                               distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"Gm")) %>%  slice_head(n=5)
#     OverlappedGenes nSampleName nConfidentBenchSite
# 1:          Trim13           2                   8
# 2:           Smyd5           2                   5
# 3:           Txnip           2                   5
# 4:           Zfp57           2                   5
# 5:            Ccnf           2                   4

#select well-studied genes
dt.selected.gene.mESC.NEB.MACS2 <- dt.gene.mm39 %>% filter(gene_name=="Trim13")
dt.selected.gene <- dt.selected.gene.mESC.NEB.MACS2
StudySamples <- c("mESC_WT1","mESC_WT2","mESC_IVT1", "mESC_IVT2")

      ## Create a page (7.5*7.5cm)
      window.size=5#100kb
      pseudocount=1
      track.height=0.5
      pdf("Fig1_NEB_mESC_Trim13_IVTpos_specific_peak.pdf")
      pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
      mm39.assembly <- assembly(BSgenome = BSgenome.Mmusculus.UCSC.mm39,Genome = "mm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",OrgDb = "org.Mm.eg.db")
      region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly)
      #plot gene
      paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly, width=3.0)
      ## Plot genes small
      genesPlot <- plotGenes(params = paramsgene,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                             x = 0.5, y = 0, height = 1.0,width = 5,just = c("left", "top"), default.units = "cm",fontsize = 7)
      annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
      #track for bigwig files of study samples

      if(dt.selected.gene$strand=="-"){
        dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
          Input.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
          RIP.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
          Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
        }
        #determine the signal range for pos and neg strand
        # range_both <- ceiling(quantile(abs(dt.density$score),1))+1
        range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
          distinct(sample,range) %>% as.data.table()
        range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
          group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
        #plot
        for(s in 1:length(StudySamples)){
          signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                   params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                   just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                     params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          #add MACS2 peak in this region
          dt.peak.region <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
          if(nrow(dt.peak.region)>0){
            for(r in 1:nrow(dt.peak.region)){
              plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                          fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
            }
          }
          plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
          sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
          sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
        }
      }else{
        dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
          Input.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
          RIP.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
          Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
          rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
        }
        #determine the signal range for pos strand
        range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
          distinct(sample,range) %>% as.data.table()
        range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
          group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
        #plot
        for(s in 1:length(StudySamples)){
          signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                   params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                   just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                                     params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
          #add MACS2 peak in this region
          dt.peak.region <- dt.MACS2.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
          if(nrow(dt.peak.region)>0){
            for(r in 1:nrow(dt.peak.region)){
              plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                          fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
            }
          }
          plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
          sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
          sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
        }
      }
      #track for benchmark m6A sites (perturbation effect)
      #METTL3cko
      signal_confident.benchmark.Mettl3cko <- plotSignal(data = dt.eTAMseq.mESC.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect), params = region,
                                                   x = 0.5, y = 1.5+length(StudySamples)*track.height, width = 5, height = track.height,
                                                   just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
      plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
      sample_label <- plotText(label="Mettl3_CKO", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
      sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

      #add legend for MeRIP
      legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+1)*track.height, width = 2.5, height = 1,
                                     just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                     title = expression("log"["2"]~"(CPM)"),fondface="bold")
      legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in eTAMseq_mESC"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(length(StudySamples)+1)*track.height,
                                                   width = 2.5, height = 1,
                                          just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fondface="bold")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(length(StudySamples)+2)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")
      dev.off()


#save intermediate variables for re-use
save.image("sm6APeak_Figure1.intermediate.results.RDS")

####################### combine all figure 1 together #############################
fig1.list <- list(p.count.RRACHratio.benchmark.m6As=p.count.RRACHratio.benchmark.m6As, p.shared.benchmark.m6A.sites=p.shared.benchmark.m6A.sites,
                  p.benchmark.perturbation.validation=p.benchmark.perturbation.validation, p.IVT.ratio.MACS2=p.IVT.ratio.MACS2,
                  p.benchmark.ratio.MACS2=p.benchmark.ratio.MACS2, p.confident.benchmark.ratio.MACS2=p.confident.benchmark.ratio.MACS2,
                  p.count.MACS2.IVTpos.specific.m6Agene.with.benchmark=p.count.MACS2.IVTpos.specific.m6Agene.with.benchmark)
saveRDS(fig1.list, file="Figure1.plot.list.RDS")

pdf("Figure1.IVT_m6Aseq_samples_is_not_perfect_negative_control_for_m6A_calling.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.count.RRACHratio.benchmark.m6As, x = 0.05, y=0.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.shared.benchmark.m6A.sites, x = 6, y=0.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 6, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.benchmark.perturbation.validation, x = 12, y=0.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 12, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.IVT.ratio.MACS2, x = 0.05, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.benchmark.ratio.MACS2, x = 6, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 6, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.confident.benchmark.ratio.MACS2, x = 12, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 12, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.MACS2.IVTpos.specific.m6Agene.with.benchmark, x = 0.05, y=12.2, default.units = "cm",width = 7.8, height = 5.8)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)
dev.off()




######### FigS1 IVT+ non-MACS2 (exomePeak2 and TRESS) peaks were compromised of noticeable true positive peak ##########

#figure s1a. Use antibody-independent m6A sites is better benchmarks than m6A peaks from IVT samples for performance estimation (IVT)
#MeTPeak default combo (enrichment fold>=2) of WT and IVT of HEK_NEB, HEK_SYSY, and HEK_Abcam
dt.MeTPeak.default.peaks <- LoadPeakMethod(Method="MeTPeak",Method.peak.dir = "/data/m6A_calling_strategy/sm6APeak/MeTPeak")
dt.MeTPeak.default.peaks <- dt.MeTPeak.default.peaks %>% dplyr::filter(seqnames %in% paste0("chr",1:22))
dt.MeTPeak.default.peaks[,.(NPeak=n_distinct(name), NPeak.FoldOver2=n_distinct(name[score>=2]),NPeak.FoldOver5=n_distinct(name[score>=5])),keyby=.(Sample)]
dt.MeTPeak.default.peaks <- dt.MeTPeak.default.peaks %>% dplyr::filter(score>=2) %>% mutate(SampleName=strsplit(Sample,split="/",fixed=T) %>% sapply("[",1)) %>%
  mutate(Condition=case_when(str_detect(Sample,"IVT")~"IVT",.default = "mRNA"), Antibody=case_when(str_detect(Sample,"SYSY") ~ "SYSY", str_detect(Sample,"Abcam") ~ "Abcam", .default = "NEB"),
         Cell=case_when(str_detect(Sample,"HEK") ~ "HEK", .default = "mESC"))
dt.sample.info <- dt.MeTPeak.default.peaks %>% distinct(Sample,SampleName,Condition,Antibody,Cell)
#obtain the IVT+ peaks (overlap >=50% in IVT1 or IVT2) and IVT- peaks for each sample
dt.MeTPeak.default.peaks.IVToverlap <- foreach(c = dt.sample.info$SampleName %>% grep(pattern="IVT",invert = T), .combine='rbind')%do%{
  foreach(ivt = dt.sample.info[Antibody==dt.sample.info$Antibody[c] & Cell==dt.sample.info$Cell[c] & Condition=="IVT",SampleName], .combine='rbind')%do%{
    bed.peak <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==dt.sample.info$SampleName[c]) %>% dplyr::select(seqnames,start,end,name,score,strand)
    bed.ivt <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==ivt) %>% dplyr::select(seqnames,start,end,name,score,strand)
    dt.overlap <- bt.intersect(a=bed.peak, b=bed.ivt, s=F, f=0.5,F=0.5, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
    dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
    dt.ivt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(IVToverlap=ifelse(sum.overlap/width>=0.5,"IVT+","IVT-"))
    bed.peak %>% left_join(x=.,y=dt.ivt.peak %>% distinct(name,IVToverlap) %>% mutate(IVTsample=ivt), by="name")
  } %>% mutate(SampleName=dt.sample.info$SampleName[c])
}
dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,IVTsample) %>%
  mutate(nPeak=n_distinct(name),IVT_pos=n_distinct(name[IVToverlap=="IVT+"]),IVT_neg=n_distinct(name[IVToverlap=="IVT-"])) %>% mutate(IVT_ratio=round(IVT_neg/IVT_pos,2)) %>%
  as.data.table() %>% distinct(SampleName,IVTsample,nPeak,IVT_pos,IVT_neg,IVT_ratio)
#         SampleName      IVTsample nPeak IVT_pos IVT_neg IVT_ratio
# 1: HEK_Abcam_mRNA1 HEK_Abcam_IVT1 13364    2233   11131      4.98
# 2: HEK_Abcam_mRNA1 HEK_Abcam_IVT2 13364    2181   11183      5.13
# 3: HEK_Abcam_mRNA2 HEK_Abcam_IVT1 13587    2344   11243      4.80
# 4: HEK_Abcam_mRNA2 HEK_Abcam_IVT2 13587    2263   11324      5.00
# 5:   HEK_NEB_mRNA1   HEK_NEB_IVT1 14742    1041   13701     13.16
# 6:   HEK_NEB_mRNA1   HEK_NEB_IVT2 14742     973   13769     14.15
# 7:   HEK_NEB_mRNA2   HEK_NEB_IVT1 18193    1133   17060     15.06
# 8:   HEK_NEB_mRNA2   HEK_NEB_IVT2 18193    1126   17067     15.16
# 9:  HEK_SYSY_mRNA1  HEK_SYSY_IVT1 19313    2495   16818      6.74
# 10:  HEK_SYSY_mRNA1  HEK_SYSY_IVT2 19313    2696   16617      6.16
# 11:  HEK_SYSY_mRNA2  HEK_SYSY_IVT1 19100    2499   16601      6.64
# 12:  HEK_SYSY_mRNA2  HEK_SYSY_IVT2 19100    2773   16327      5.89
# 13:        mESC_WT1      mESC_IVT1 12547     678   11869     17.51
# 14:        mESC_WT1      mESC_IVT2 12547     771   11776     15.27
# 15:        mESC_WT2      mESC_IVT1 13242     570   12672     22.23
# 16:        mESC_WT2      mESC_IVT2 13242     938   12304     13.12
#obtain the overlap of benchmark sites for each peaks
#load HEK benchmark m6As
HEK293.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="HEK293",value=T)
names(HEK293.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="HEK293",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.HEK293 <-  foreach(f = 1:length(HEK293.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(HEK293.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}
#load confident HEK benchmark sites
dt.GLORI.HEK293.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/GLORI/dt.GLORI.HEK293.confident.m6As.RDS")
#load mESC benchmark m6As
mESC.benchmark.m6As.files <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = T) %>% grep(pattern="mESC",value=T)
names(mESC.benchmark.m6As.files) <- list.files("/data/m6A_calling_strategy/RIPPeakS_resource/",pattern = "(benchmark).+(m6As.bed)$",full.names = F) %>%
  grep(pattern="mESC",value=T) %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.benchmark.m6A.mESC <-  foreach(f = 1:length(mESC.benchmark.m6As.files),.combine = 'rbind')%do%{
  dt.bench <- fread(mESC.benchmark.m6As.files[f]) %>% distinct()
  colnames(dt.bench) <- c("seqnames","start","end","name","score","strand","benchmark_name","nRep","m6AGene","Motif")
  dt.bench
}
#load confident mESC benchmark sites
dt.eTAMseq.mESC.confident <- readRDS(file="/data/m6A_calling_strategy/Benchmark_m6As/eTAMseq/dt.eTAMseq.mESC.confident.m6As.RDS")

dt.MeTPeak.default.peaks.Benchoverlap <- foreach(c = dt.sample.info$SampleName %>% grep(pattern="IVT",invert = T), .combine='rbind')%do%{
  bed.peak <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==dt.sample.info$SampleName[c]) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(dt.sample.info$SampleName[c],"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }else{
    bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(dt.sample.info$SampleName[c],"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(SampleName=dt.sample.info$SampleName[c])
}
#the count/ratio of benchmark+ peak among IVT+ and IVT- peaks
dt.MeTPeak.IVT.benchmark.association <- dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) #classify a peak as IVT- if not overlap with any IVT samples
dt.MeTPeak.IVT.benchmark.association <- dt.MeTPeak.IVT.benchmark.association %>% left_join(x=.,y=dt.MeTPeak.default.peaks.Benchoverlap,by=c("SampleName","name")) %>%
  group_by(SampleName,IVToverlapped) %>% mutate(TotalPeak=n_distinct(name)) %>%
  group_by(SampleName,IVToverlapped,TotalPeak,benchmark_name) %>% mutate(nPeak.bench=n_distinct(name[!is.na(nBenchSite) & nBenchSite>0]),SumBenchSites=sum(nBenchSite[!is.na(nBenchSite)])) %>%
  mutate(nPeak.bench.confident=n_distinct(name[!is.na(nConfidentBenchSite) & nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(SampleName,IVToverlapped,TotalPeak,benchmark_name,nPeak.bench,SumBenchSites,nPeak.bench.confident) %>% dplyr::filter(!is.na(benchmark_name))
#         SampleName IVToverlapped TotalPeak benchmark_name nPeak.bench SumBenchSites nPeak.bench.confident
# 1: HEK_Abcam_mRNA1          IVT-     10562   GLORI_HEK293        7750         48787                  3713
# 2: HEK_Abcam_mRNA1          IVT+      2802   GLORI_HEK293        1419          7448                   431
# 3: HEK_Abcam_mRNA2          IVT-     10656   GLORI_HEK293        7858         49602                  3752
# 4: HEK_Abcam_mRNA2          IVT+      2931   GLORI_HEK293        1484          7723                   412
# 5:   HEK_NEB_mRNA1          IVT-     13429   GLORI_HEK293       10064         81244                  4591
# 6:   HEK_NEB_mRNA1          IVT+      1313   GLORI_HEK293         866          6244                   344
# 7:   HEK_NEB_mRNA2          IVT-     16718   GLORI_HEK293       12129         93315                  5247
# 8:   HEK_NEB_mRNA2          IVT+      1475   GLORI_HEK293         960          6517                   357
# 9:  HEK_SYSY_mRNA1          IVT-     16069   GLORI_HEK293       11601         86038                  4936
# 10:  HEK_SYSY_mRNA1          IVT+      3244   GLORI_HEK293        1788         11089                   694
# 11:  HEK_SYSY_mRNA2          IVT-     15824   GLORI_HEK293       11406         82150                  4802
# 12:  HEK_SYSY_mRNA2          IVT+      3276   GLORI_HEK293        1763         10898                   686
# 13:        mESC_WT1          IVT-     11456   eTAMseq_mESC        5970         29331                  2488
# 14:        mESC_WT1          IVT-     11456     GLORI_mESC        7340         44289                  2394
# 15:        mESC_WT1          IVT+      1091   eTAMseq_mESC         667          3059                   334
# 16:        mESC_WT1          IVT+      1091     GLORI_mESC         808          4647                   328
# 17:        mESC_WT2          IVT-     12036   eTAMseq_mESC        5876         27451                  2408
# 18:        mESC_WT2          IVT-     12036     GLORI_mESC        7427         42220                  2308
# 19:        mESC_WT2          IVT+      1206   eTAMseq_mESC         713          3200                   338
# 20:        mESC_WT2          IVT+      1206     GLORI_mESC         865          4809                   335
#figure S1A the proportion of IVT+ (MeTPeak) peaks  among different cell x antibody (four sample)
dt.IVT.ratio.MeTPeak <- dt.MeTPeak.IVT.benchmark.association %>% left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  distinct(SampleName,Antibody,Cell,IVToverlapped,TotalPeak) %>% pivot_wider(id_cols = c("SampleName","Antibody","Cell"), names_from = "IVToverlapped", values_from = "TotalPeak") %>%
  as.data.table() %>% mutate(PeakPerSample=`IVT-`+`IVT+`,IVTratio=`IVT+`/(`IVT-`+`IVT+`)) %>% group_by(Antibody,Cell) %>% mutate(avg.PeakPerSample=mean(PeakPerSample),avg.IVTratio=mean(IVTratio)) %>%
  distinct(Antibody,Cell,avg.PeakPerSample,avg.IVTratio) %>% mutate(GroupName=paste0(Antibody,"_",Cell)) %>% as.data.table()
dt.IVT.ratio.MeTPeak <- dt.IVT.ratio.MeTPeak %>% distinct(avg.PeakPerSample, avg.IVTratio, GroupName) %>%
  mutate(GroupName=str_replace(GroupName,pattern="HEK",replacement="HEK293"),avg.PeakPerSample=round(avg.PeakPerSample,0)) %>%
  mutate(label=paste0("NoPeak=",avg.PeakPerSample),label.pct=scales::percent(avg.IVTratio, accuracy = 0.1))
dt.IVT.ratio.MeTPeak <- dt.IVT.ratio.MeTPeak %>% dplyr::arrange(desc(avg.IVTratio)) %>% mutate(GroupName=factor(GroupName,levels=GroupName))
p.IVT.ratio.MeTPeak <- ggplot(dt.IVT.ratio.MeTPeak,
                              aes(y = GroupName, x = avg.IVTratio, fill = GroupName)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  # percent label to the right of the bar
  geom_text(aes(label = label.pct, x =  avg.IVTratio + 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  #N peak label to the left of the bar
  geom_text(aes(label = label, x = 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  scale_fill_manual(values = c("#ffd47f","#f7c1cf","#7b93c6","#acd9ee"),breaks = c("Abcam_HEK293","SYSY_HEK293","NEB_HEK293","NEB_mESC")) +
  labs(x = "Proportion of IVT+ MeTPeak peaks ", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.IVT.ratio.MeTPeak
##figure S1B the proportion of (MeTPeak) peaks contains benchmark sites among IVT+ and IVT- peak among different cellxantibody
dt.benchmark.ratio.MeTPeak <- dt.MeTPeak.IVT.benchmark.association %>% mutate(benchmark_ratio=nPeak.bench/TotalPeak) %>%
  left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  group_by(Antibody,Cell,benchmark_name,IVToverlapped) %>% mutate(avg.benchmark_ratio=mean(benchmark_ratio)) %>% as.data.table() %>%
  distinct(Antibody,Cell,benchmark_name,IVToverlapped,avg.benchmark_ratio) %>% dplyr::arrange(Antibody,Cell,benchmark_name,IVToverlapped)
dt.benchmark.ratio.MeTPeak <- dt.benchmark.ratio.MeTPeak %>% mutate(GroupName=paste0(Antibody,"_",Cell," | ",strsplit(benchmark_name,split="_") %>% sapply("[",1))) %>%
  dplyr::arrange(IVToverlapped,desc(avg.benchmark_ratio)) %>% mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>% dplyr::arrange(IVToverlapped,GroupName)
ttest.p.benchmark.ratio.MeTPeak <- signif(t.test(x=dt.benchmark.ratio.MeTPeak[IVToverlapped=="IVT+",avg.benchmark_ratio],y=dt.benchmark.ratio.MeTPeak[IVToverlapped=="IVT-",avg.benchmark_ratio],paired = T )$p.value,2)
p.benchmark.ratio.MeTPeak <- ggplot(data=dt.benchmark.ratio.MeTPeak,aes(x=GroupName, y=avg.benchmark_ratio,fill=IVToverlapped))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#fcbb46","#ef756c"),breaks=c("IVT+", "IVT-"))+
  labs(x=NULL,y="% of peaks harbor benchmark m6A sites", subtitle = paste0("paired t-test p-value=", ttest.p.benchmark.ratio.MeTPeak))+
  guides(fill=guide_legend(title = "MeTPeak peaks",nrow=1))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.95,0.9),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.benchmark.ratio.MeTPeak
#figure S1C the proportion of (MeTPeak) peaks contains high confident benchmark sites (perturbation experiment validated) among IVT+ and IVT- peak among different cell x antibody
dt.confident.benchmark.ratio.MeTPeak <- dt.MeTPeak.IVT.benchmark.association %>% mutate(confident.benchmark_ratio=nPeak.bench.confident/TotalPeak) %>%
  left_join(x=.,y=dt.sample.info %>% dplyr::select(SampleName,Antibody,Cell),by="SampleName") %>%
  group_by(Antibody,Cell,benchmark_name,IVToverlapped) %>% mutate(avg.confident.benchmark_ratio=mean(confident.benchmark_ratio)) %>% as.data.table() %>%
  distinct(Antibody,Cell,benchmark_name,IVToverlapped,avg.confident.benchmark_ratio) %>% dplyr::arrange(Antibody,Cell,benchmark_name,IVToverlapped)
dt.confident.benchmark.ratio.MeTPeak <- dt.confident.benchmark.ratio.MeTPeak %>% mutate(GroupName=paste0(Antibody,"_",Cell," | ",strsplit(benchmark_name,split="_") %>% sapply("[",1))) %>%
  dplyr::arrange(IVToverlapped,desc(avg.confident.benchmark_ratio)) %>% mutate(GroupName=factor(GroupName,levels=unique(GroupName))) %>% dplyr::arrange(IVToverlapped,GroupName)
ttest.p.confident.benchmark.ratio.MeTPeak <- signif(t.test(x=dt.confident.benchmark.ratio.MeTPeak[IVToverlapped=="IVT+",avg.confident.benchmark_ratio],
                                                           y=dt.confident.benchmark.ratio.MeTPeak[IVToverlapped=="IVT-",avg.confident.benchmark_ratio],paired = T )$p.value,2)
p.confident.benchmark.ratio.MeTPeak <- ggplot(data=dt.confident.benchmark.ratio.MeTPeak,aes(x=GroupName, y=avg.confident.benchmark_ratio,fill=IVToverlapped))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#fcbb46","#ef756c"),breaks=c("IVT+", "IVT-"))+
  labs(x=NULL,y="% of peaks harbor confident benchmark m6A sites", subtitle = paste0("paired t-test p-value=", ttest.p.confident.benchmark.ratio.MeTPeak))+
  guides(fill=guide_legend(title = "MeTPeak peaks",nrow = 1))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.95,0.9),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.confident.benchmark.ratio.MeTPeak


#compare the  IVT+ and IVT- peaks with random simulation peaks regarding containing GLORI sites
dt.MeTPeak.default.peaks.IVToverlap.type <- dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,score,SampleName,IVToverlapped)
dt.MeTPeak.default.peaks.IVToverlap.type[,.N,by=.(SampleName,IVToverlapped)]
#kept only genic peaks as input for simulation
library(bedtoolsr)
options(bedtools.path = "~/anaconda3/envs/m6A_seq/bin")
dt.MeTPeak.default.peaks.genic.annot <- rbind(annot_peak(dt.MeTPeak.default.peaks.IVToverlap %>% filter(!str_detect(name,"mESC")) %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                       strand = TRUE,fract = 0.5,gtf = "/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf", annot_type = c("gene", "feature", "both")[1]) %>%
                                              dplyr::select(seqnames,start,end,name,score,strand,OverlappedGeneIDs),
                                            annot_peak(dt.MeTPeak.default.peaks.IVToverlap %>% filter(str_detect(name,"mESC")) %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                       strand = TRUE,fract = 0.5,gtf = "/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf", annot_type = c("gene", "feature", "both")[1]) %>%
                                              dplyr::select(seqnames,start,end,name,score,strand,OverlappedGeneIDs)
)

dt.MeTPeak.default.genic.peaks <- dt.MeTPeak.default.peaks.genic.annot %>% filter(!is.na(OverlappedGeneIDs))
dt.MeTPeak.default.genic.peaks <- dt.MeTPeak.default.genic.peaks %>% left_join(x=.,y=dt.MeTPeak.default.peaks.IVToverlap.type,by=c("name","score"))

#estimate cumulative curve for real IVT+ and IVT- peaks
dt.MeTPeak.default.genic.peaks <- dt.MeTPeak.default.genic.peaks %>% left_join(x=.,y=dt.MeTPeak.default.peaks.Benchoverlap %>% filter(str_detect(benchmark_name,"GLORI")) %>% distinct(name,benchmark_name,nBenchSite),by="name")
dt.MeTPeak.default.genic.peaks <- dt.MeTPeak.default.genic.peaks %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite))
dt.cum.data.MeTPeak.default.genic.peaks <- dt.MeTPeak.default.genic.peaks %>%
  #group_by(SampleName,IVToverlapped) %>%
  group_by(SampleName) %>%
  dplyr::arrange(desc(score)) %>%
  mutate(peak_width = end - start,
         cum_length = cumsum(peak_width),       # X-axis: Cumulative length
         cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
  ) %>% as.data.table() %>%
  mutate(IVToverlapped="all") %>% 
  distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites) %>% dplyr::arrange(SampleName,IVToverlapped,cum_length,cum_sites)

#write function to simulation and estimate cumulative curve
#use regionR package to generate random peaks for each samples IVT+ and IVT- peaks
library(regioneR)
library(GenomicRanges)
library(GenomicFeatures)

#simulation for HEK293 IVT+ and IVT- peaks
Txdb.hg38 <- makeTxDbFromGFF("~/genome_db/gencode.v44.annotation.gtf")
#extract all genic region to serve as the "allowed universe"
genes_universe <- genes(Txdb.hg38)
exons_universe <- exons(Txdb.hg38)
#use loop to generate 1000x for per sample per IVToverlapped condition and pick the 99%, 50% and 1% of simulation according to cumulative benchmark sites
dt.cum.data.MeTPeak.random.HEK.peaks <- foreach(S = unique(dt.MeTPeak.default.genic.peaks$SampleName) %>% grep(pattern="mESC",invert = T,value = T), .combine = 'rbind')%do%{
  # foreach(IVT = unique(dt.MeTPeak.default.genic.peaks$IVToverlapped), .combine = 'rbind')%do%{
    # real_peaks <- dt.MeTPeak.default.genic.peaks %>% filter(SampleName==S & IVToverlapped==IVT) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
   real_peaks <- dt.MeTPeak.default.genic.peaks %>% filter(SampleName==S) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
    real_peaks <- makeGRangesFromDataFrame(real_peaks,keep.extra.columns = T)
    #generate random peak
    # t1 <- Sys.time()
    nSimulation=100
    registerDoParallel(cl=10)
    dt.random.peaks <- foreach(R = 1:nSimulation, .combine = 'rbind')%dopar%{
      random_peaks <- randomizeRegions(real_peaks, genome = "hg38", allow.overlaps = FALSE, mask = NULL, 
                                       # universe = genes_universe)
                                       universe = exons_universe)
      #transform random peak into data table
      as.data.table(random_peaks) %>% mutate(name=paste("random",R,seqnames,start,end,strand,sep="_"),score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand)
    }
    stopImplicitCluster()
    # Sys.time() - t1
    #overlap with benchmark GLORI sites
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
    dt.overlap <- bt.intersect(a=dt.random.peaks, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)#key parameter
    dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
    # dt.overlap
    #benchmark sites with high confidence
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])
    dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
      as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
    dt.overlap <- dt.overlap %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2))
    #determine the 99%, 50% and 1% of simulation according to cumulative benchmark sites
    dt.overlap.cumsites <- dt.overlap %>% group_by(RandomID) %>% mutate(CumSites=sum(nBenchSite)) %>% as.data.table() %>%
      distinct(RandomID,CumSites) %>% dplyr::arrange(desc(CumSites)) 
    dt.randomid.selected <- dt.overlap.cumsites[c(0.01,0.5,0.99)*nSimulation,] %>% mutate(RandomRank=c("99%","50%","1%"))
    #output selected random peaks cum data for visualization
    dt.selected.overlap <- dt.overlap %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
    dt.selected.randompeaks <- dt.random.peaks %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2)) %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
    dt.selected.randompeaks <- dt.selected.randompeaks %>% left_join(x=.,y=dt.selected.overlap,by=c("name","RandomID")) %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite),
                                                                                                                                   nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))
    #cumulative output for selected random peaks
    dt.cum.data.random.peaks <- dt.selected.randompeaks %>% group_by(RandomID) %>% mutate(peak_width = end - start,
                                                                                          cum_length = cumsum(peak_width),       # X-axis: Cumulative length
                                                                                          cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
    ) %>% as.data.table() %>% distinct(RandomID,name,nBenchSite,cum_length,cum_sites) %>% left_join(x=.,y=dt.randomid.selected,by="RandomID")
    # message(paste0("Finished simulation of ", S, " ", IVT, " peaks"))
    message(paste0("Finished simulation of ", S, " ", " peaks"))
    rm(dt.random.peaks)
    rm(dt.overlap)
    # dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped=IVT)
    dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped="all")
  # }
}
dt.cum.data.MeTPeak.random.HEK.peaks <- dt.cum.data.MeTPeak.random.HEK.peaks %>% filter(IVToverlapped=="all") 
#simulation for mESC IVT+ and IVT- peaks
Txdb.mm39 <- makeTxDbFromGFF("~/genome_db/gencode.vM33.annotation.gtf")
#extract all genic region to serve as the "allowed universe"
genes_universe <- genes(Txdb.mm39)
exons_universe <- exons(Txdb.mm39)
#use loop to generate 1000x for per sample per IVToverlapped condition and pick the 99%, 50% and 1% of simulation according to cumulative benchmark sites
dt.cum.data.MeTPeak.random.mESC.peaks <- foreach(S = unique(dt.MeTPeak.default.genic.peaks$SampleName) %>% grep(pattern="mESC",invert = F,value = T), .combine = 'rbind')%do%{
  # foreach(IVT = unique(dt.MeTPeak.default.genic.peaks$IVToverlapped), .combine = 'rbind')%do%{
    # real_peaks <- dt.MeTPeak.default.genic.peaks %>% filter(SampleName==S & IVToverlapped==IVT) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
    real_peaks <- dt.MeTPeak.default.genic.peaks %>% filter(SampleName==S) %>% dplyr::distinct(seqnames,start,end,name,score,strand)
    real_peaks <- makeGRangesFromDataFrame(real_peaks,keep.extra.columns = T)
    #generate random peak
    # t1 <- Sys.time()
    nSimulation=100
    registerDoParallel(cl=10)
    dt.random.peaks <- foreach(R = 1:nSimulation, .combine = 'rbind')%dopar%{
      random_peaks <- randomizeRegions(real_peaks, genome = "mm39", allow.overlaps = FALSE, mask = NULL, 
                                       # universe = genes_universe)
                                       universe = exons_universe)
      #transform random peak into data table
      as.data.table(random_peaks) %>% mutate(name=paste("random",R,seqnames,start,end,strand,sep="_"),score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand)
    }
    stopImplicitCluster()
    # Sys.time() - t1
    #overlap with benchmark GLORI sites
    bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::filter(benchmark_name=="GLORI_mESC") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
    dt.overlap <- bt.intersect(a=dt.random.peaks, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)##key parameter 
    dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
    dt.overlap
    #benchmark sites with high confidence
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
    dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
      as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
    dt.overlap <- dt.overlap %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2))
    #determine the 99%, 50% and 1% of simulation according to cumulative benchmark sites
    dt.overlap.cumsites <- dt.overlap %>% group_by(RandomID) %>% mutate(CumSites=sum(nBenchSite)) %>% as.data.table() %>%
      distinct(RandomID,CumSites) %>% dplyr::arrange(desc(CumSites)) 
    dt.randomid.selected <- dt.overlap.cumsites[c(0.01,0.5,0.99)*nSimulation,] %>% mutate(RandomRank=c("99%","50%","1%"))
    #output selected random peaks cum data for visualization
    dt.selected.overlap <- dt.overlap %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
    dt.selected.randompeaks <- dt.random.peaks %>% mutate(RandomID=strsplit(name,split="_",fixed=T) %>% sapply("[",2)) %>% dplyr::filter(RandomID %in% unique(dt.randomid.selected$RandomID))
    dt.selected.randompeaks <- dt.selected.randompeaks %>% left_join(x=.,y=dt.selected.overlap,by=c("name","RandomID")) %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite),
                                                                                                                                   nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))
    #cumulative output for selected random peaks
    dt.cum.data.random.peaks <- dt.selected.randompeaks %>% group_by(RandomID) %>% mutate(peak_width = end - start,
                                                                                          cum_length = cumsum(peak_width),       # X-axis: Cumulative length
                                                                                          cum_sites = cumsum(nBenchSite)       # Y-axis: Cumulative GLORI sites recovered
                                                                                          ) %>%
      as.data.table() %>% distinct(RandomID,name,nBenchSite,cum_length,cum_sites) %>% left_join(x=.,y=dt.randomid.selected,by="RandomID")
    # message(paste0("Finished simulation of ", S, " ", IVT, " peaks"))
    message(paste0("Finished simulation of ", S, " ", " peaks"))
    rm(dt.random.peaks)
    rm(dt.overlap)
    # dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped=IVT)
    dt.cum.data.random.peaks %>% mutate(SampleName=S,IVToverlapped="all")
  # }
}
dt.cum.data.MeTPeak.random.mESC.peaks <- dt.cum.data.MeTPeak.random.mESC.peaks %>% filter(IVToverlapped=="all") 

#plot cumulative curve of MeTPeak IVT+/ IVT- peaks and simulated peaks 
dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks <- rbind(dt.cum.data.MeTPeak.default.genic.peaks %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites) %>% mutate(PeakGroup="MeTPeak"),
                                                               dt.cum.data.MeTPeak.random.HEK.peaks %>% mutate(PeakGroup=paste0("Random",RandomRank)) %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites,PeakGroup),
                                                               dt.cum.data.MeTPeak.random.mESC.peaks %>% mutate(PeakGroup=paste0("Random",RandomRank)) %>% distinct(SampleName,IVToverlapped,name,nBenchSite,cum_length,cum_sites,PeakGroup)) %>%
  dplyr::arrange(SampleName,IVToverlapped,PeakGroup,cum_length,cum_sites)
dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>% slice_max(cum_length) %>% as.data.table() %>% dplyr::arrange(SampleName,IVToverlapped,PeakGroup)
#scaled cum_length and cum_sites for each Sample x IVToverlapped:use the median cum_sites as 100%
#pruned the extreme long cum_length x cum_sites (keep only 100 data points)
dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks[,.N,by=.(SampleName,IVToverlapped,PeakGroup)]
dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks <- foreach(S = unique(dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks$SampleName),.combine = 'rbind')%do%{
  foreach(IVT = unique(dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks$IVToverlapped),.combine = 'rbind')%do%{
    dt.cumdata <- dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% dplyr::filter(SampleName==S & IVToverlapped==IVT)
    foreach(G = unique(dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks$PeakGroup), .combine = 'rbind')%do%{
      dt.cumdata.subgroup <- dt.cumdata %>% dplyr::filter(PeakGroup==G) %>% dplyr::arrange(cum_length,cum_sites)
      dt.cumdata.subgroup[floor(nrow(dt.cumdata.subgroup)*seq(0,1,0.001)[-1]),] %>% mutate(scaled_cum_peaks=seq(0,1,0.001)[-1]*100)
    }
  }
}
#scaled the cum_sites, and cum_peaks
dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks <- dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% filter(PeakGroup %in% c("MeTPeak","Random50%","Random99%")) %>%
  group_by(SampleName,IVToverlapped,PeakGroup) %>% mutate(scaled_cum_sites=cum_sites/max(cum_sites)*100) %>% as.data.table()
#average of sample group
dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks <- dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% dplyr::filter(IVToverlapped=="all") %>%
  mutate(Sample=case_when(str_detect(SampleName,"HEK_Abcam") ~ "Abcam_HEK293",
                          str_detect(SampleName,"HEK_NEB") ~ "NEB_HEK293",
                          str_detect(SampleName,"HEK_SYSY") ~ "SYSY_HEK293",
                          str_detect(SampleName,"mESC_WT") ~ "NEB_mESC")) %>%
  group_by(Sample,PeakGroup,scaled_cum_peaks) %>% mutate(scaled_cum_peaks=mean(scaled_cum_peaks) %>% round(digits = 2),scaled_cum_sites=mean(scaled_cum_sites) %>% round(digits = 2)) %>%
  as.data.table() %>% distinct(Sample,IVToverlapped,PeakGroup,scaled_cum_peaks,scaled_cum_sites) %>%
  mutate(PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MeTPeak")))

p.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks <- ggplot(dt.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% dplyr::filter(IVToverlapped=="all"),
                                                                      aes(x = scaled_cum_peaks, y = scaled_cum_sites)) +
  geom_line(aes(linetype=PeakGroup,color=PeakGroup),linewidth=0.5) +
  geom_abline(intercept = 0,slope = 1,color="grey50",linetype="solid",linewidth=0.25)+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  scale_linetype_manual(values=c("solid","dashed","dotted"),breaks = c("MeTPeak","Random99%","Random50%"))+
  # scale_y_continuous(transform = "log10")+
  # geom_hline(yintercept = c(1000),color="grey50",linetype="solid",linewidth=0.25)+
  labs(
    title = NULL,
    x = "Scaled cumulative peak proportion (%)",
    y = "Scaled cumulative recovered GLORI sites (%)",
  ) +
  guides(color=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+
  # coord_cartesian(xlim=c(0,2.5))+
  facet_wrap(~Sample,ncol = 2)+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  guides(
    colour = guide_legend(
      override.aes = list(
        linetype = c("dotted","dashed","solid")   # ←←← 根据你实际的线型顺序修改
      ),
      title = "PeakGroup",
      nrow = 1,                    # 可选：横向排列
      keywidth = 1.5
    ),
    linetype = "none"              # 隐藏独立的 linetype legend
  )
p.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks

#the TPR  and FPR comparison of top25% v.s. tail25%
dt.top25vstail25.MeTPeak.default.genic.peaks <-  rbind(
  dt.MeTPeak.default.genic.peaks %>% group_by(SampleName) %>% dplyr::arrange(desc(score)) %>% mutate(totalsites=sum(nBenchSite)) %>% slice_head(prop = 0.25) %>%
    mutate(cum_sites_ratio=sum(nBenchSite)/totalsites*100,FP_peak_ratio=n_distinct(name[nBenchSite==0])/n_distinct(name)*100,Group="top25%") %>%
    as.data.table() %>% distinct(SampleName,Group,cum_sites_ratio,FP_peak_ratio),
  dt.MeTPeak.default.genic.peaks %>% group_by(SampleName) %>% dplyr::arrange(desc(score)) %>% mutate(totalsites=sum(nBenchSite)) %>% slice_tail(prop = 0.25) %>%
    mutate(cum_sites_ratio=sum(nBenchSite)/totalsites*100,FP_peak_ratio=n_distinct(name[nBenchSite==0])/n_distinct(name)*100,Group="tail25%") %>%
    as.data.table() %>% distinct(SampleName,Group,cum_sites_ratio,FP_peak_ratio)) %>%
  dplyr::arrange(SampleName,Group)
dt.top25vstail25.MeTPeak.default.genic.peaks 
dt.top25vstail25.MeTPeak.default.genic.peaks <- dt.top25vstail25.MeTPeak.default.genic.peaks %>% mutate(Group=factor(Group,levels=c("tail25%","top25%")))
#pair-wise baxplot visualize TPR  and FPR comparison of top25% v.s. tail25%
p.top25vstail25.MeTPeak.cum.sites.ratio <- ggplot(dt.top25vstail25.MeTPeak.default.genic.peaks, aes(x = Group, y = cum_sites_ratio, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values = c("#F3D28A", "#F39030"),breaks = c("tail25%","top25%")) +
  labs(y = str_wrap("Cumulative proportion (%) of benchmark m6As",width = 40), x = "MeTPeak peak subset") +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.cum.sites.ratio <- dt.top25vstail25.MeTPeak.default.genic.peaks %>%
  group_by(Group) %>%
  summarise(mean_val = round(mean(cum_sites_ratio), 2))
p.top25vstail25.MeTPeak.cum.sites.ratio <- p.top25vstail25.MeTPeak.cum.sites.ratio +
  geom_text(data = mean_values.cum.sites.ratio,
            aes(x = Group, 
                y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$cum_sites_ratio) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.top25vstail25.MeTPeak.cum.sites.ratio <- p.top25vstail25.MeTPeak.cum.sites.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$cum_sites_ratio) * 1.05,
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$cum_sites_ratio) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.top25vstail25.MeTPeak.cum.sites.ratio
p.top25vstail25.MeTPeak.FP.peak.ratio <- ggplot(dt.top25vstail25.MeTPeak.default.genic.peaks, aes(x = Group, y = FP_peak_ratio, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values = c("#F3D28A", "#F39030"),breaks = c("tail25%","top25%")) +
  labs(y = str_wrap("Proportion (%) of peaks contain no benchmark m6As",width=45), x = "MeTPeak peak subset") +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.FP.peak.ratio <- dt.top25vstail25.MeTPeak.default.genic.peaks %>%
  group_by(Group) %>%
  summarise(mean_val = round(mean(FP_peak_ratio), 2))
p.top25vstail25.MeTPeak.FP.peak.ratio <- p.top25vstail25.MeTPeak.FP.peak.ratio +
  geom_text(data = mean_values.FP.peak.ratio,
            aes(x = Group, 
                y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$FP_peak_ratio) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.top25vstail25.MeTPeak.FP.peak.ratio <- p.top25vstail25.MeTPeak.FP.peak.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$FP_peak_ratio) * 1.05,
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.top25vstail25.MeTPeak.default.genic.peaks$FP_peak_ratio) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.top25vstail25.MeTPeak.FP.peak.ratio
cowplot::plot_grid(p.top25vstail25.MeTPeak.cum.sites.ratio,p.top25vstail25.MeTPeak.FP.peak.ratio) 
# save.image("updated_sm6APeak_FigureS1.intermediate.results.RDS")
#barplot to display the max cumulative GLORI sites of MeTPeak IVT+/IVT- peaks vs RandomPeaks (50% and 99%)
dt.density.GLORI.MeTPeakvsRandom <- dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>%
  mutate(TotalPeak=n(),CumWidth=max(cum_length)/1e6,CumSites=max(cum_sites)) %>% as.data.table() %>% distinct(SampleName,IVToverlapped,PeakGroup,TotalPeak,CumWidth,CumSites) %>%
  dplyr::arrange(SampleName,IVToverlapped,PeakGroup,TotalPeak,CumWidth,CumSites) %>% mutate(density.perMb=CumSites/CumWidth, density.per1000Peak=CumSites/TotalPeak*1000)
dt.density.GLORI.MeTPeakvsRandom <- dt.density.GLORI.MeTPeakvsRandom %>% dplyr::filter(PeakGroup !="Random1%") %>%
  mutate(SampleName=factor(SampleName),PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MeTPeak")))
#keep all peak 
p.density.GLORI.MeTPeakvsRandom.per1000Peak  <- ggplot(data=dt.density.GLORI.MeTPeakvsRandom%>% dplyr::filter(IVToverlapped=="all"),
                                                     aes(x=SampleName,y=density.per1000Peak, fill=PeakGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  labs(x=NULL,y="Density of GLORI m6A sites (per 1000 peak)")+
  guides(fill=guide_legend(title = "PeakGroup"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),transform = "log10")+
  # facet_wrap(~IVToverlapped,scales = "fixed")+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.density.GLORI.MeTPeakvsRandom.per1000Peak
#add enrich fold v.s Random50%
dt.density.fold.GLORI.MeTPeakvsRandom.per1000Peak <- dt.density.GLORI.MeTPeakvsRandom%>% dplyr::filter(IVToverlapped=="all") %>%
  group_by(SampleName) %>% mutate(DensityFold=density.per1000Peak[PeakGroup=="MeTPeak"]/density.per1000Peak[PeakGroup=="Random50%"]) %>%
  as.data.table() %>% distinct(SampleName,DensityFold)
p.density.GLORI.MeTPeakvsRandom.per1000Peak <- p.density.GLORI.MeTPeakvsRandom.per1000Peak+
  geom_text(data = dt.density.fold.GLORI.MeTPeakvsRandom.per1000Peak,
            inherit.aes = FALSE,
            aes(x = SampleName, 
                y = max(dt.density.GLORI.MeTPeakvsRandom$density.per1000Peak) * 1.05,
                label = paste(round(DensityFold),"x")),
            size = 2.5, fontface = "bold", color = "#EA9E58")
p.density.GLORI.MeTPeakvsRandom.per1000Peak

#compare the cumulative benchmark m6As of MeTPeak vs. Random Peaks
dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks <- dt.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks %>% group_by(SampleName,IVToverlapped,PeakGroup) %>%
  mutate(FPR=round(n_distinct(name[nBenchSite==0])/n_distinct(name)*100,2),TotalBenchSites=max(cum_sites)) %>%
  mutate(TPR=dplyr::case_when(str_detect(SampleName,"HEK") ~ round(TotalBenchSites/nrow(dt.benchmark.m6A.HEK293[benchmark_name=="GLORI_HEK293",])*100,2),
                              str_detect(SampleName,"mESC") ~ round(TotalBenchSites/nrow(dt.benchmark.m6A.mESC[benchmark_name=="GLORI_mESC",])*100,2))) %>% 
  as.data.table() %>% distinct(SampleName,IVToverlapped,PeakGroup,FPR,TPR)
dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks <- dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks %>% dplyr::filter(PeakGroup %in% c("MeTPeak","Random50%","Random99%")) %>%
  mutate(PeakGroup=factor(PeakGroup,levels=c("Random50%","Random99%","MeTPeak")))
p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- ggplot(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks, aes(x = PeakGroup, y = TPR, fill = PeakGroup)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color=PeakGroup),size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  scale_y_continuous(expand = expansion(add = c(0,0.3)),breaks = seq(0,max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$TPR),by=10))+
  labs(y = str_wrap("Cumulative proportion (%) of benchmark m6As",width = 40), x = NULL) +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.cum.sites.ratio <- dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks %>%
  group_by(PeakGroup) %>%
  summarise(mean_val = round(mean(TPR), 2))
p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio +
  geom_text(data = mean_values.cum.sites.ratio,
            aes(x = PeakGroup, 
                y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$TPR) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio <- p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     comparisons = list(c("MeTPeak","Random50%")),
                     bracket.size = 0.1,
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$TPR) * 0.8,
                     label.y.npc = "bottom",
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$TPR) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio

#compare the ratio of peak contain benchmark m6As of MeTPeak vs Random Peaks
p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- ggplot(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks, aes(x = PeakGroup, y = FPR, fill = PeakGroup)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(aes(color=PeakGroup),size = 1, alpha = 0.9,position = position_jitter(width = 0.2,seed=999)) +
  geom_line(aes(group = SampleName), color = "grey60", alpha = 0.6, linewidth = 0.6,position = position_jitter(width = 0.2,seed=999)) +   # 配对连线
  scale_fill_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  scale_color_manual(values=c("#7BC0CD","#4198AC","#EA9E58"),breaks=c("Random50%","Random99%","MeTPeak"))+
  scale_y_continuous(expand = expansion(add = c(0,0.3)),breaks = seq(0,max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$FPR),by=10))+
  labs(y = str_wrap("Proportion (%) of peaks contain no benchmark m6As",width = 40), x = NULL) +
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    legend.position = "none", # Adjust position inside the plot
    legend.background = element_rect(fill="transparent")
  )
# add mean value
mean_values.FP.peak.ratio<- dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks %>%
  group_by(PeakGroup) %>%
  summarise(mean_val = round(mean(FPR), 2))
p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio +
  geom_text(data = mean_values.FP.peak.ratio,
            aes(x = PeakGroup, 
                y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$FPR) * 1.01,
                label = paste0("avg =", mean_val,"%")),
            size = 2.5, fontface = "plain", color = "black")
#add pvalue
p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio <- p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio+
  stat_compare_means(method = "wilcox.test", 
                     paired = TRUE, 
                     comparisons = list(c("MeTPeak","Random50%")),
                     bracket.size = 0.1,
                     label = "p.format",
                     label.x = 1.35,
                     label.y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$FPR) * 0.8,
                     label.y.npc = "bottom",
                     size=3) +
  annotate("text", x = 1.15, y = max(dt.TPR.FPR.MeTPeak.default.genic.peaks.vs.random.peaks$FPR) * 1.06,
           label = "Wilcox", size = 2.5, hjust = 1)
p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio



save.image("updated_sm6APeak_FigureS1.intermediate.results.RDS")

#figure S1D the count of m6A+ genes (MeTPeak) that could be missed if use IVT as negative control
#HEK293
#annotate the MeTPeak peak of HEK293 to gene
dt.MeTPeak.default.peaks.IVToverlap.HEK293 <- dt.MeTPeak.default.peaks.IVToverlap %>% dplyr::filter(str_detect(SampleName,"HEK")) %>% distinct(seqnames, start,end,name,score, strand)
dt.MeTPeak.default.peaks.IVToverlap.HEK293.annotgene <- annot_peak(peak.bed=dt.MeTPeak.default.peaks.IVToverlap.HEK293,
                                                                   strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",#check the strandness
                                                                   annot_type = "gene")
#determine the IVT+ specific m6A+ genes
dt.MeTPeak.default.peaks.IVToverlap.HEK293.annotgene <- dt.MeTPeak.default.peaks.IVToverlap.HEK293.annotgene %>%
  left_join(x=.,y=dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
              as.data.table() %>% distinct(name,SampleName,IVToverlapped), by=c("name"))
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293 <-   dt.MeTPeak.default.peaks.IVToverlap.HEK293.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName,OverlappedGenes) %>% mutate(IVTpos.specific=all(IVToverlapped=="IVT+")) %>% as.data.table() %>%
  distinct(SampleName,OverlappedGenes,IVTpos.specific) %>% dplyr::filter(IVTpos.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293 <- dt.MeTPeak.IVTpos.specific.m6Agene.HEK293 %>%
  left_join(x=.,y=dt.MeTPeak.default.peaks.IVToverlap.HEK293.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ","), by=c("SampleName","OverlappedGenes"))
#obtain the IVT+ peaks overlapped with (high confident) benchmark sites
dt.benchmark.MeTPeak.peak.HEK293 <-  dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) %>% dplyr::filter(IVToverlapped=="IVT+") %>% dplyr::filter(str_detect(SampleName,"HEK")) %>%
  inner_join(x=.,y=dt.MeTPeak.default.peaks.Benchoverlap %>% dplyr::filter(nBenchSite>0),by=c("name","SampleName"))
#collect IVT+ specific m6A genes with benchmark m6A sites support
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.HEK293 %>%
  inner_join(x=.,y=dt.benchmark.MeTPeak.peak.HEK293 %>% dplyr::select(name,benchmark_name,nBenchSite, nConfidentBenchSite),by="name")
#rank the IVT+ specific m6A gene with benchmark m6A sites support
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.sample.info %>% inner_join(x=.,y=dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark,by="SampleName")
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% group_by(Antibody,Cell,OverlappedGenes) %>%
  mutate(nSampleName=n_distinct(SampleName)) %>% as.data.table() %>% dplyr::arrange(Antibody,Cell,desc(nSampleName),desc(nConfidentBenchSite))
#obtain final IVT+ specific m6A+ gene list with benchmark sites
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% distinct(Antibody,Cell,nSampleName,OverlappedGenes,benchmark_name,nBenchSite,nConfidentBenchSite)
#meanwhile,obtain final IVT+ specific m6A+ gene list with confident benchmark sites
dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark %>% dplyr::filter(nConfidentBenchSite>0)

#mESC
#annotate the MeTPeak peak of mESC to gene
dt.MeTPeak.default.peaks.IVToverlap.mESC <- dt.MeTPeak.default.peaks.IVToverlap %>% dplyr::filter(str_detect(SampleName,"mESC")) %>% distinct(seqnames, start,end,name,score, strand)
dt.MeTPeak.default.peaks.IVToverlap.mESC.annotgene <- annot_peak(peak.bed=dt.MeTPeak.default.peaks.IVToverlap.mESC,
                                                                 strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                                                 annot_type = "gene")
dt.MeTPeak.default.peaks.IVToverlap.mESC.annotgene <- dt.MeTPeak.default.peaks.IVToverlap.mESC.annotgene %>%
  left_join(x=.,y=dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
              as.data.table() %>% distinct(name,SampleName,IVToverlapped), by=c("name"))
#determine the IVT+ specific m6A+ genes
dt.MeTPeak.IVTpos.specific.m6Agene.mESC <-   dt.MeTPeak.default.peaks.IVToverlap.mESC.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName,OverlappedGenes) %>% mutate(IVTpos.specific=all(IVToverlapped=="IVT+")) %>% as.data.table() %>%
  distinct(SampleName,OverlappedGenes,IVTpos.specific) %>% dplyr::filter(IVTpos.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)
dt.MeTPeak.IVTpos.specific.m6Agene.mESC <- dt.MeTPeak.IVTpos.specific.m6Agene.mESC %>%
  left_join(x=.,y=dt.MeTPeak.default.peaks.IVToverlap.mESC.annotgene %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ","), by=c("SampleName","OverlappedGenes"))
#obtain the IVT+ peaks overlapped with (high confident) benchmark sites
dt.benchmark.MeTPeak.peak.mESC <-  dt.MeTPeak.default.peaks.IVToverlap %>% group_by(SampleName,name) %>%  mutate(IVToverlapped=case_when(all(IVToverlap=="IVT-") ~ "IVT-", .default = "IVT+")) %>%
  as.data.table() %>% distinct(name,SampleName,IVToverlapped) %>% dplyr::filter(IVToverlapped=="IVT+") %>% dplyr::filter(str_detect(SampleName,"mESC")) %>%
  inner_join(x=.,y=dt.MeTPeak.default.peaks.Benchoverlap %>% dplyr::filter(nBenchSite>0),by=c("name","SampleName"))
#collect IVT+ specific m6A genes with benchmark m6A sites support
dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.mESC %>%
  inner_join(x=.,y=dt.benchmark.MeTPeak.peak.mESC %>% dplyr::select(name,benchmark_name,nBenchSite, nConfidentBenchSite),by="name")
#rank the IVT+ specific m6A gene with benchmark m6A sites support
dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.sample.info %>% inner_join(x=.,y=dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark,by="SampleName")
dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark %>% group_by(Antibody,Cell,OverlappedGenes) %>%
  mutate(nSampleName=n_distinct(SampleName)) %>% as.data.table() %>% dplyr::arrange(Antibody,Cell,desc(nSampleName),desc(nConfidentBenchSite))
#obtain final IVT+ specific m6A+ gene list with benchmark sites
dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark %>% distinct(Antibody,Cell,nSampleName,OverlappedGenes,benchmark_name,nBenchSite,nConfidentBenchSite)
#meanwhile,obtain final IVT+ specific m6A+ gene list with confident benchmark sites
dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.confident.benchmark <- dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark %>% dplyr::filter(nConfidentBenchSite>0)

#calculate the count of IVTpos specific m6A+ gene in each AntibodyxCell
dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark <- rbind(rbind(dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark[nBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)],
                                                                       dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark[nBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)]) %>% mutate(Group="Contain benchmark m6A sites"),
                                                                 rbind(dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.benchmark[nConfidentBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)],
                                                                       dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.benchmark[nConfidentBenchSite>0,.(nGene=n_distinct(OverlappedGenes)),keyby=.(Antibody,Cell)]) %>% mutate(Group="Contain confident benchmark m6A sites")
)
dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark <- dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark %>% mutate(SampleLabel=paste0(Antibody,"_",str_replace(Cell,"HEK","HEK293"))) %>%
  dplyr::arrange(Group,desc(nGene)) %>% dplyr::mutate(SampleLabel=factor(SampleLabel,levels=unique(SampleLabel)),
                                                      Group=factor(Group, levels=c("Contain benchmark m6A sites","Contain confident benchmark m6A sites"), labels=c("Benchmark m6As","Confident benchmark m6As")))
dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark <- dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark %>% mutate(label.count=paste0(paste0("No.Gene=",nGene)))
p.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark  <- ggplot(data=dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark, aes(x=SampleLabel,y=nGene, fill=Group))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark %>% dplyr::filter(Group=="Benchmark m6As"),
            aes(label = label.count, y =  nGene/2, x=SampleLabel),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark %>% dplyr::filter(Group=="Confident benchmark m6As"),
            aes(label = label.count, y =  nGene/2, x=SampleLabel),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y="N of genes harbor only IVT+ MeTPeak peaks")+
  guides(fill=guide_legend(title = "IVT+ peaks contain"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.9,0.8),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark
#figure S1E-F. the IGV plot of selected m6A+ genes harbor IVT+ MeTPeak peaks that support by confident m6A sites in HEK293 and mESC
dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

#HEK293 IVT+ specific m6A+ gene with high confident benchmark m6As
selected.m6Agene.SYSY.HEK <-  dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="SYSY" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=5)
selected.m6Agene.SYSY.HEK
# OverlappedGenes nSampleName nConfidentBenchSite
# 1:           FOXC1           2                  15
# 2:         C6orf47           2                  10
# 3:     C6orf47-AS1           2                  10
# 4:           PRR11           2                  10
# 5:       PRR11-AS1           2                  10
selected.m6Agene.NEB.HEK <-  dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="NEB" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=5)
selected.m6Agene.NEB.HEK
# OverlappedGenes nSampleName nConfidentBenchSite
# 1:        ITPRIPL1           2                  14
# 2:        ITPRIPL1           2                  13
# 3:            RRP1           2                   7
# 4:        UBE2V2P3           2                   7
# 5:        C11orf24           2                   6
selected.m6Agene.NEB.Abcam <-  dt.MeTPeak.IVTpos.specific.m6Agene.HEK293.with.confident.benchmark %>% dplyr::filter(Antibody=="Abcam" & Cell=="HEK" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>%  slice_head(n=20)
selected.m6Agene.NEB.Abcam
# OverlappedGenes nSampleName nConfidentBenchSite

#select well-studied genes
# dt.selected.gene.HEK.SYSY.MeTPeak <- dt.gene.hg38 %>% filter(gene_name=="FOXC1")
# dt.selected.gene <- dt.selected.gene.HEK.SYSY.MeTPeak
# dt.selected.gene.HEK.NEB.MeTPeak <- dt.gene.hg38 %>% filter(gene_name=="C11orf24")
# dt.selected.gene <- dt.selected.gene.HEK.NEB.MeTPeak
dt.selected.gene.HEK.Abcam.MeTPeak <- dt.gene.hg38 %>% filter(gene_name=="ARF6")
dt.selected.gene <- dt.selected.gene.HEK.Abcam.MeTPeak
#track for bigwig files of study samples
# StudySamples <- c("HEK_SYSY_mRNA1", "HEK_SYSY_mRNA2","HEK_SYSY_IVT1","HEK_SYSY_IVT2")
# StudySamples <- c("HEK_NEB_mRNA1", "HEK_NEB_mRNA2","HEK_NEB_IVT1","HEK_NEB_IVT2")
StudySamples <- c("HEK_Abcam_mRNA1", "HEK_Abcam_mRNA2","HEK_Abcam_IVT1","HEK_Abcam_IVT2")
## Create a page (7.5*7.5cm)
window.size=1#100kb
pseudocount=1
track.height=0.5
# pdf("FigS1_SYSY_HEK_FOXC1_IVTpos_specific_peak.pdf")
# pdf("FigS1_NEB_HEK_C11orf24_IVTpos_specific_peak.pdf")
pdf("FigS1_Abcam_HEK_ARF6_IVTpos_specific_peak.pdf")

pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
#plot gene
paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38", width=3.0)
## Plot genes small
genesPlot <- plotGenes(params = paramsgene,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = 5,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos and neg strand
  # range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
    distinct(sample,range) %>% as.data.table()
  range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
    group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    #add MeTPeak peak in this region
    dt.peak.region <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
    if(nrow(dt.peak.region)>0){
      for(r in 1:nrow(dt.peak.region)){
        plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                    fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
      }
    }
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos strand
  range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
    distinct(sample,range) %>% as.data.table()
  range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
    group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    #add MeTPeak peak in this region
    dt.peak.region <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
    if(nrow(dt.peak.region)>0){
      for(r in 1:nrow(dt.peak.region)){
        plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                    fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
      }
    }
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
#track for benchmark m6A sites (perturbation effect)
#FTO GLORI HEK293
signal_confident.benchmark.FTO <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.FTO), params = region,
                                             x = 0.5, y = 1.5+length(StudySamples)*track.height, width = 5, height = track.height,
                                             just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="FTO treatment", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
#METTL3i GLORI HEK293
signal_confident.benchmark.METTL3i <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.METTL3i), params = region,
                                                 x = 0.5, y = 1.5+(length(StudySamples)+1)*track.height, width = 5, height = track.height,
                                                 just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+1.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="METTL3 inhibition", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+1.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in GLORI_HEK293"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(length(StudySamples)+2)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#mESC IVT+ specific m6A+ gene with high confident benchmark m6As
selected.m6Agene.mESC.NEB <-  dt.MeTPeak.IVTpos.specific.m6Agene.mESC.with.confident.benchmark %>% dplyr::filter(Antibody=="NEB" & Cell=="mESC" & nSampleName==2) %>%
  distinct(OverlappedGenes,nSampleName,nConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"Gm")) %>%  slice_head(n=10)
selected.m6Agene.mESC.NEB
#     OverlappedGenes nSampleName nConfidentBenchSite
# 6:          Tspyl1           2                   7
# 7:            Gcat           2                   6
# 8:            H1f0           2                   6
# 9:           Mcm10           2                   6

#select well-studied genes
dt.selected.gene.mESC.NEB.MeTPeak <- dt.gene.mm39 %>% filter(gene_name=="Tspyl1")
dt.selected.gene <- dt.selected.gene.mESC.NEB.MeTPeak
StudySamples <- c("mESC_WT1","mESC_WT2","mESC_IVT1", "mESC_IVT2")

## Create a page (7.5*7.5cm)
window.size=2#100kb
pseudocount=1
track.height=0.5
pdf("FigS1_NEB_mESC_Tspyl1_IVTpos_specific_peak.pdf")
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
mm39.assembly <- assembly(BSgenome = BSgenome.Mmusculus.UCSC.mm39,Genome = "mm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",OrgDb = "org.Mm.eg.db")
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly)
#plot gene
paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly, width=3.0)
## Plot genes small
genesPlot <- plotGenes(params = paramsgene,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = 5,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
#track for bigwig files of study samples

if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos and neg strand
  # range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
    distinct(sample,range) %>% as.data.table()
  range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
    group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    #add MeTPeak peak in this region
    dt.peak.region <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
    if(nrow(dt.peak.region)>0){
      for(r in 1:nrow(dt.peak.region)){
        plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                    fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
      }
    }
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data/m6A_calling_strategy/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos strand
  range_both <- dt.density %>% as.data.table() %>% group_by(sample) %>% mutate(range=ceiling(quantile(abs(score),1))+1) %>%
    distinct(sample,range) %>% as.data.table()
  range_both <- range_both %>%  mutate(study=strsplit(sample,"_") %>% sapply("[",1)) %>%
    group_by(study) %>% mutate(range=max(range)) %>% distinct(sample,study,range) %>% as.data.table()
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = 5, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    #add MeTPeak peak in this region
    dt.peak.region <- dt.MeTPeak.default.peaks %>% dplyr::filter(SampleName==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
    if(nrow(dt.peak.region)>0){
      for(r in 1:nrow(dt.peak.region)){
        plotgardener::annoHighlight(signal_Input, chrom = dt.peak.region$seqnames[r], chromstart = dt.peak.region$start[r], chromend = dt.peak.region$end[r],
                                    fill = "#eb4601", linecolor = "transparent",height = track.height, alpha = 0.35,y = 1.5+(s-1)*track.height,just = c("left", "top"), default.units = "cm")
      }
    }
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
#track for benchmark m6A sites (perturbation effect)
#METTL3cko
signal_confident.benchmark.Mettl3cko <- plotSignal(data = dt.eTAMseq.mESC.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect), params = region,
                                                   x = 0.5, y = 1.5+length(StudySamples)*track.height, width = 5, height = track.height,
                                                   just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="Mettl3_CKO", fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+1)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in eTAMseq_mESC"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(length(StudySamples)+1)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(length(StudySamples)+2)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#save intermediate variables for re-use
save.image("sm6APeak_FigureS1.intermediate.results.RDS")

####################### combine all figureS1 together
figs1.list <- list(p.IVT.ratio.MeTPeak=p.IVT.ratio.MeTPeak,
                   p.benchmark.ratio.MeTPeak=p.benchmark.ratio.MeTPeak,
                   p.confident.benchmark.ratio.MeTPeak=p.confident.benchmark.ratio.MeTPeak,
                   p.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark=p.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark)
saveRDS(figs1.list, file="FigureS1.plot.list.RDS")

pdf("FigureS1_MeTPeak_IVTpos_peaks were_composed_of_noticeable_true_positive_m6A_modification.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.IVT.ratio.MeTPeak, x = 0.05, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.benchmark.ratio.MeTPeak, x = 6, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 6, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.confident.benchmark.ratio.MeTPeak, x = 12, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 12, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.MeTPeak.IVTpos.specific.m6Agene.with.benchmark, x = 0.05, y=12.2, default.units = "cm",width = 7.8, height = 5.8)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)
dev.off()


##### FigS1 new IVT+ and IVT- peaks from mRNA both significantly enriched for benchmark m6As compared with random peak##########
#load the saved figures
load("updated_sm6APeak_Figure1.intermediate.results.RDS")
# pdf("FigureS1new_IVTpos_IVTneg_MACS2_m6A_peaks_from_mRNA_samples_enriched_for_benchmark_m6As.pdf",width = 8.2677, height = 11.693)
# pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
# plotGG(plot = p.density.GLORI.MACS2vsRandom.per1000Peak, x = 0.2, y=0.05, default.units = "cm",width = 8.8, height = 5.8)
# plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.density.GLORI.MACS2vsRandom.perMb, x = 0.2, y=6, default.units = "cm",width = 8.8, height = 5.8)
# plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05, y = 6,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks, x = 0.2, y=12, default.units = "cm",width = 8.8, height = 11.8)
# plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 12,just = c("top","left"), default.units = "cm",draw=T)
# dev.off()
#new edition
pdf("FigureS1new_MACS2_m6A_peaks_from_mRNA_samples_enriched_for_benchmark_m6As.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotGG(plot = p.MACS2.default.genic.peaks.vs.random.peaks.cum.sites.ratio, x = 0.2, y=0.05, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.MACS2.default.genic.peaks.vs.random.peaks.FP.peak.ratio, x = 4.6, y=0.05, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 4.45, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
#row2
plotGG(plot = p.density.GLORI.MACS2vsRandom.per1000Peak, x = 0.2, y=5.55, default.units = "cm",width = 8.8, height = 5.8)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 5.55,just = c("top","left"), default.units = "cm",draw=T)
#row3
plotGG(plot = p.pruned.cum.data.MACS2.default.genic.peaks.vs.random.peaks, x = 0.2, y=11.5, default.units = "cm",width = 8.8, height = 7.3)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 11.5,just = c("top","left"), default.units = "cm",draw=T)
#row4
plotGG(plot = p.top25vstail25.MACS2.cum.sites.ratio, x = 0.2, y=19, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05, y = 19,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top25vstail25.MACS2.FP.peak.ratio, x = 4.6, y=19, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 4.45, y = 19,just = c("top","left"), default.units = "cm",draw=T)
dev.off()


#load the saved figures
load("updated_sm6APeak_FigureS1.intermediate.results.RDS")
# pdf("FigureS1new_IVTpos_IVTneg_MeTPeak_m6A_peaks_from_mRNA_samples_enriched_for_benchmark_m6As.pdf",width = 8.2677, height = 11.693)
# pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
# plotGG(plot = p.density.GLORI.MeTPeakvsRandom.per1000Peak, x = 9.2, y=0.05, default.units = "cm",width = 8.8, height = 5.8)
# plotText(label = "D", fontsize = 8, fontface = "bold",x = 9.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.density.GLORI.MeTPeakvsRandom.perMb, x = 9.2, y=6, default.units = "cm",width = 8.8, height = 5.8)
# plotText(label = "E", fontsize = 8, fontface = "bold",x = 9.05, y = 6,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks, x = 9.2, y=12, default.units = "cm",width = 8.8, height = 11.8)
# plotText(label = "F", fontsize = 8, fontface = "bold",x = 9.05, y = 12,just = c("top","left"), default.units = "cm",draw=T)
# dev.off()

pdf("FigureS1new_MeTPeak_m6A_peaks_from_mRNA_samples_enriched_for_benchmark_m6As.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotGG(plot = p.MeTPeak.default.genic.peaks.vs.random.peaks.cum.sites.ratio, x = 9.2, y=0.05, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 9.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.MeTPeak.default.genic.peaks.vs.random.peaks.FP.peak.ratio, x = 13.6, y=0.05, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "H", fontsize = 8, fontface = "bold",x = 13.45, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
#row2
plotGG(plot = p.density.GLORI.MeTPeakvsRandom.per1000Peak, x = 9.2, y=5.55, default.units = "cm",width = 8.8, height = 5.8)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 9.05, y = 5.55,just = c("top","left"), default.units = "cm",draw=T)
#row3
plotGG(plot = p.pruned.cum.data.MeTPeak.default.genic.peaks.vs.random.peaks, x = 9.2, y=11.5, default.units = "cm",width = 8.8, height = 7.3)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 9.05, y = 11.5,just = c("top","left"), default.units = "cm",draw=T)
#row4
plotGG(plot = p.top25vstail25.MeTPeak.cum.sites.ratio, x = 9.2, y=19, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 9.05, y = 19,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top25vstail25.MeTPeak.FP.peak.ratio, x = 13.6, y=19, default.units = "cm",width = 4.2, height = 5.3)
plotText(label = "L", fontsize = 8, fontface = "bold",x = 13.45, y = 19,just = c("top","left"), default.units = "cm",draw=T)
dev.off()