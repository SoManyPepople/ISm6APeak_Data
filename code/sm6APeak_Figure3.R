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

### Part3 Incorporation of stranded peak intensity lead to higher performance of peak calling for stranded MeRIPseq data###############
#3a The stranded RIP/Input enrichment of false positive peaks of best combo (NEB_HEK293 and NEB_mESC)
#load the exomePeak best combo m6A peak in HEK293
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE)
dt.pruned.FPR.TPR.AUC.NEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HEK.all.combo.RDS")
dt.pruned.FPR.TPR.AUC.NEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.mESC.all.combo.RDS")
dt.GLORI.AUC.NEB <- rbind(dt.pruned.FPR.TPR.AUC.NEB.HEK %>% distinct(Method,ComboName,Benchmark,AUC.scaled, TPR20),
                          dt.pruned.FPR.TPR.AUC.NEB.mESC %>% distinct(Method,ComboName,Benchmark,AUC.scaled, TPR20)) %>%
  left_join(x=.,y=dt.method.combo %>% dplyr::select(Method,ComboName,IsDefaultCombo))

dt.BestCombo.NEB <- dt.GLORI.AUC.NEB %>% filter(str_detect(Benchmark,"HEK293")) %>% dplyr::arrange(Benchmark,Method) %>% group_by(Benchmark,Method) %>% dplyr::arrange(desc(AUC.scaled), desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()

dt.method.best.AUC.NEB <- rbind(dt.pruned.FPR.TPR.AUC.NEB.HEK %>% inner_join(x=.,y=dt.BestCombo.NEB %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")),
                                dt.pruned.FPR.TPR.AUC.NEB.mESC %>% inner_join(x=.,y=dt.BestCombo.NEB %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")))
dt.method.best.AUC.NEB %>% distinct(Method,ComboName,Benchmark)
#highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
dt.method.best.AUC.NEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
dt.method.best.AUC.NEB <- dt.method.best.AUC.NEB %>% mutate(Cell = case_when(str_detect(Benchmark,"HEK")~ "HEK293", .default = "mESC")) %>% mutate(Cell=factor(Cell,levels=c("HEK293","mESC"))) %>%
  dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))
dt.method.best.AUC.NEB %>% distinct(Method,ComboName,Benchmark) %>% filter(Method=="exomePeak") %>% distinct(Method,ComboName)

#obtain the corresponding cutoff and method combo for exomePeak best combo
dt.method.best.AUC.NEB %>% filter(Method=="exomePeak") %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
# FPR  TPR cutoff    Method ComboName                          Benchmark  AUC AUC.scaled maxTPRmFPR TPR20   Cell                    label
# 1: 0.2 0.39   9.42 exomePeak   Combo26 dt.benchmark.GLORI_HEK293.m6As.bed 0.18       0.36       0.21  0.39 HEK293 exomePeak:ScaledAUC=0.36
dt.selected.method.combo.cutoff.FPR20 <- dt.method.best.AUC.NEB %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>% filter(!is.na(cutoff) & TPR>0) %>%
                                                 dplyr::select(Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Best",Combo=paste0(Method,"_",ComboName))
# Method ComboName cutoff FPR  TPR Group              Combo
# 1: exomePeak2   Combo39   1.81 0.2 0.54  Best exomePeak2_Combo39
# 2:  exomePeak   Combo26   9.42 0.2 0.39  Best  exomePeak_Combo26
# 3:    MeTPeak    Combo2   6.17 0.2 0.38  Best     MeTPeak_Combo2
# 4:      TRESS   Combo27   1.28 0.2 0.36  Best      TRESS_Combo27
# 5: MeRIPtools   Combo13   3.08 0.2 0.26  Best MeRIPtools_Combo13
#load exomePeak2 (Default and Best Combo) for HEK_NEB and mESC_NEB
#HEK293
dt.Method.BestCombo.FPR20.m6A.HEK293 <- foreach(C = 1:nrow(dt.selected.method.combo.cutoff.FPR20),.combine = rbind)%do%{
  LoadPeakMethod(Method = dt.selected.method.combo.cutoff.FPR20$Method[C], Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB/", dt.selected.method.combo.cutoff.FPR20$Combo[C])) %>%
    mutate(Group=paste0(dt.selected.method.combo.cutoff.FPR20$Method[C],"Best")) %>% mutate(PeakOverCutoff=score>= dt.selected.method.combo.cutoff.FPR20$cutoff[C])
}
dt.Method.BestCombo.FPR20.m6A.HEK293 %>% group_by(Group) %>% slice_head(n=1)
dt.Method.BestCombo.FPR20.m6A.HEK293[,.N,by=.(Group,Sample,PeakOverCutoff)] %>% dplyr::arrange(Sample,Group,PeakOverCutoff)

#divided into false positive peak(FP peak), false negative peak (FP peak), and other ture positive peak
#overlap with benchmark sites to define FP, TP, and FN
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
dt.Method.BestCombo.FPR20.HEK293.peaks.Benchmarkoverlap <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.HEK293$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.HEK293 %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#combine benchmark overlap into peak
dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.HEK293 %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.HEK293.peaks.Benchmarkoverlap  %>% filter(benchmark_name=="GLORI_HEK293") %>% dplyr::distinct(name,Sample,nBenchSite,nConfidentBenchSite), by=c("name","Sample")) %>%
  mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite))
#FP and TP, and FN
dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap %>%
  mutate(PeakGroup=case_when(PeakOverCutoff==TRUE & nBenchSite>0 ~ "TP", PeakOverCutoff==FALSE & nBenchSite>0 ~ "FN", PeakOverCutoff==TRUE & nBenchSite==0 ~ "FP", .default = "TN"))
dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap[,.N,by=.(Sample,Group,PeakGroup)] %>% dplyr::arrange(Sample,Group,PeakGroup)


#mESC
dt.Method.BestCombo.FPR20.m6A.mESC <- foreach(C = 1:nrow(dt.selected.method.combo.cutoff.FPR20),.combine = rbind)%do%{
  LoadPeakMethod(Method = dt.selected.method.combo.cutoff.FPR20$Method[C], Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB_mESC/", dt.selected.method.combo.cutoff.FPR20$Combo[C])) %>%
    mutate(Group=paste0(dt.selected.method.combo.cutoff.FPR20$Method[C],"Best")) %>% mutate(PeakOverCutoff=score>= dt.selected.method.combo.cutoff.FPR20$cutoff[C])
}
dt.Method.BestCombo.FPR20.m6A.mESC %>% group_by(Group) %>% slice_head(n=1)
dt.Method.BestCombo.FPR20.m6A.mESC[,.N,by=.(Group,Sample,PeakOverCutoff)] %>% dplyr::arrange(Sample,Group,PeakOverCutoff)

#divided into false positive peak(FP peak), false negative peak (FP peak), and other ture positive peak
#overlap with benchmark sites to define FP, TP, and FN
#load HEK benchmark m6As
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
dt.Method.BestCombo.FPR20.mESC.peaks.Benchmarkoverlap <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.mESC$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.mESC %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"mESC")){
    bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#combine benchmark overlap into peak
dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.mESC %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.mESC.peaks.Benchmarkoverlap  %>% filter(benchmark_name=="GLORI_mESC") %>% dplyr::distinct(name,Sample,nBenchSite,nConfidentBenchSite), by=c("name","Sample")) %>%
  mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite))
#FP and TP, and FN
dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap %>%
  mutate(PeakGroup=case_when(PeakOverCutoff==TRUE & nBenchSite>0 ~ "TP", PeakOverCutoff==FALSE & nBenchSite>0 ~ "FN", PeakOverCutoff==TRUE & nBenchSite==0 ~ "FP", .default = "TN"))
dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap[,.N,by=.(Sample,Group,PeakGroup)] %>% dplyr::arrange(Sample,Group,PeakGroup)



#plot the count of FP and FN peaks of all method best combo at FPR20
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count  <- rbind(dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN"),
                                                 dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN")) %>%
                                           group_by(Sample,Method,Group,PeakGroup) %>% mutate(nPeak=n_distinct(name)) %>% as.data.table() %>% distinct(Sample,Method,Group,PeakGroup,nPeak)
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count %>% mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293",.default = "mESC"))
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count %>% mutate(Method=factor(Method,levels=rev(c("MeRIPtools","TRESS","MeTPeak","exomePeak","exomePeak2"))),
                                                                                                        PeakGroup=factor(PeakGroup,levels=c("TP","FP","FN")))
#add grey rect
# indices to highlight
pos <- seq_along(levels(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count$Method))
idx <- pos[pos %% 2!=0]
rects <- data.frame(xmin = pos[idx] - 0.5,xmax = pos[idx] + 0.5,ymin = -Inf,ymax = Inf)
p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count <- ggplot(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count, aes(x = Method, y = nPeak)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(aes(color=PeakGroup),fill="transparent",position = position_dodge2(padding = 0.05,width=0.6),width = 0.55,outlier.shape = NA,size = 0.4,alpha=0.7) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(aes(x = Method, y = nPeak, color=PeakGroup, shape=Cell), position = position_dodge2(padding = 0.05,width = 0.6),size=1.5, stroke = 0.6,alpha=0.7) +
  geom_rect(data = rects,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),inherit.aes = FALSE,fill = "grey92", alpha = 0.2)+
  # scales for shapes and fills
  scale_shape_manual(values = c(18,22),breaks=c("HEK293","mESC"), name = "NEB sample") +
  scale_color_manual(values = c("#e9abac","#b5d2e8","#cee4b4"), breaks = c("TP","FP","FN"), name = "PeakType") +
  labs(x = NULL, y = "N of m6A peak from method's best combo") +
  guides(shape=guide_legend(nrow=2),color=guide_legend(nrow=2))+
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )


#check the positive and negative RIP/Input enrichment of FP peak on positive and negative strand
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN  <- rbind(dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN"),
                                                 dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN"))
tmpdir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmpdir)){dir.create(tmpdir)}
#prepare method m6A
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN  <- rbind(dt.Method.BestCombo.FPR20.m6A.HEK293.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN"),
                                                       dt.Method.BestCombo.FPR20.m6A.mESC.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN"))
dt.Method.m6a.HEK <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN %>% dplyr::filter(str_detect(Sample,"HEK"))
#prepare the stranded.BAM.Depth
#calculate stranded bam depth
Samples <- unique(dt.Method.m6a.HEK$Sample)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/Figure3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID:mean_Ex_bin100_P1
selected.option <- "mean_Ex_bin100_P1"
dt.Method.m6a <- dt.Method.m6a.HEK %>% mutate(ComboName="Best") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
# org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
# org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite"
# org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
# org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed"
# org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/hg38_gencode.gene.bed"
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo="Best", .combine='rbind')%do%{
    foreach(s=Samples, .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      # selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",
                 intron.bed = "/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
   }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
# intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)

#mESC
dt.Method.m6a.mESC <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN %>% dplyr::filter(str_detect(Sample,"mESC"))
#prepare the stranded.BAM.Depth
#calculate stranded bam depth
Samples <- unique(dt.Method.m6a.mESC$Sample)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/Figure3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID:mean_Ex_bin100_P1
selected.option <- "mean_Ex_bin100_P1"
dt.Method.m6a <- dt.Method.m6a.mESC %>% mutate(ComboName="Best") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo="Best", .combine='rbind')%do%{
    foreach(s=Samples, .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      # selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_mm39.txt",
                 intron.bed = "/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.mm39.bed",
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
  }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
# intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)

#load calculate bedtools-based peak intensity in both pos strand and neg strand
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- foreach(S = unique(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN$Sample), .combine='rbind')%do%{
  foreach(M = as.character(unique(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN$Method)), .combine = 'rbind')%do%{
    dt.intensity <- fread(paste0(stranded.bam.dir,"/",S,"__",M,"__Best_",selected.option,"_bedtools_intensity.tsv"),header=T) %>%
      left_join(x=.,y=dt.Method.BestCombo.FPR20.m6A.TP.FP.FN %>% dplyr::select(name,strand,Sample,Method,PeakGroup), by="name")
    dt.intensity
  }
}
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity  %>% group_by(Sample,Method,PeakGroup) %>% mutate(nPeak=n_distinct(name), nPeak.samestrand=n_distinct(name[strand.new==strand])) %>%
  as.data.table() %>% distinct(Sample,Method,PeakGroup,nPeak,nPeak.samestrand) %>% dplyr::arrange(Sample,Method,PeakGroup)
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>%
  mutate(score.samestrand=case_when(strand=="+" ~ pos, strand=="-" ~ neg), score.reversestrand = case_when(strand=="-" ~ pos, strand=="+" ~ neg))

dt.bothstrand.Method.BestCombo.FPR20 <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% dplyr::filter(PeakGroup=="TP")
p.bothstrand.Method.BestCombo.FPR20.TP <-
  ggplot(dt.bothstrand.Method.BestCombo.FPR20, aes(x = score.samestrand, y = score.reversestrand)) +
  geom_hline(yintercept = 2,color = "grey60", linetype = "dashed")+
  geom_vline(xintercept = 2,color = "grey60", linetype = "dashed")+
  ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d') +
  viridis::scale_color_viridis(name = "Point density", option = "magma", trans = "log") +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  # facet_wrap(~PeakGroup,ncol=3)+
  labs(x = "Intensity on same strand of peak", y = "Intensity on reverse strand of peak", title="FP peaks") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = c(0.7,0.5), plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.15)), plot.subtitle = element_text(hjust = 0.5))
dt.bothstrand.Method.BestCombo.FPR20 <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% dplyr::filter(PeakGroup=="FP")
p.bothstrand.Method.BestCombo.FPR20.FP <-
  ggplot(dt.bothstrand.Method.BestCombo.FPR20, aes(x = score.samestrand, y = score.reversestrand)) +
  geom_hline(yintercept = 2,color = "grey60", linetype = "dashed")+
  geom_vline(xintercept = 2,color = "grey60", linetype = "dashed")+
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d'),dpi=300) +
  viridis::scale_color_viridis(name = "Density", option = "magma", trans = "log",  labels = function(x) scales::scientific(x,  digits = 2)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  # facet_wrap(~PeakGroup,ncol=3)+
  labs(x = "Intensity on same strand of peak", y = "Intensity on reverse strand of peak", title="FP peaks") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = c(0.8,0.5),legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
        plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1)),
        plot.subtitle = element_text(hjust = 0.5))
dt.bothstrand.Method.BestCombo.FPR20 <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% dplyr::filter(PeakGroup=="FN")
p.bothstrand.Method.BestCombo.FPR20.FN <-
  ggplot(dt.bothstrand.Method.BestCombo.FPR20, aes(x = score.samestrand, y = score.reversestrand)) +
  geom_hline(yintercept = 2,color = "grey60", linetype = "dashed")+
  geom_vline(xintercept = 2,color = "grey60", linetype = "dashed")+
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d'),dpi = 300) +
  viridis::scale_color_viridis(name = "Density", option = "magma", trans = "log",  labels = function(x) scales::scientific(x,  digits = 2)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  # facet_wrap(~PeakGroup,ncol=3)+
  labs(x = "Intensity on same strand of peak", y = "Intensity on reverse strand of peak", title="FN peaks") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = c(0.8,0.5),legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
        plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.15)), plot.subtitle = element_text(hjust = 0.5))

cowplot::plot_grid( p.bothstrand.Method.BestCombo.FPR20.FP, p.bothstrand.Method.BestCombo.FPR20.FN,nrow = 1)

# dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% group_by(PeakGroup,Method,Sample) %>%
#   mutate(nPeak=n_distinct(name), nPeak.samestrand.enrichment = n_distinct(name[score.samestrand>=2]), nPeak.reversestrand.enrichment = n_distinct(name[score.reversestrand>=2])) %>% as.data.table() %>%
#   dplyr::distinct(PeakGroup,Method,Sample,nPeak,nPeak.samestrand.enrichment,nPeak.reversestrand.enrichment) %>%
#   mutate(ratio.enrich.samestrand=nPeak.samestrand.enrichment/nPeak, ratio.enrich.reversestrand=nPeak.reversestrand.enrichment/nPeak) %>% dplyr::arrange(Sample,Method,PeakGroup)
# dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% group_by(PeakGroup,Method) %>%
#   mutate(avg.ratio.enrich.samestrand=mean(ratio.enrich.samestrand), avg.ratio.enrich.reversestrand=mean(ratio.enrich.reversestrand)) %>% as.data.table() %>%
#   arrange(Method,PeakGroup)
dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% group_by(PeakGroup,Method,Sample) %>%
  mutate(EnrichGroup=case_when(score.reversestrand>=2 & (score.samestrand<2 | is.na(score.samestrand)) ~ "ReverseStrandOnly",
                               score.samestrand>=2 & (score.reversestrand<2 | is.na(score.reversestrand)) ~ "SameStrandOnly",
                               score.samestrand>=2 & score.reversestrand>=2 ~ "BothStrand", .default = "NeitherStrand")) %>%
  mutate(nPeak=n_distinct(name)) %>% mutate(Ratio.sames.only = n_distinct(name[EnrichGroup=="SameStrandOnly"])/nPeak, Ratio.reverse.only = n_distinct(name[EnrichGroup=="ReverseStrandOnly"])/nPeak,
                                            Ratio.both = n_distinct(name[EnrichGroup=="BothStrand"])/nPeak, Ratio.neither=n_distinct(name[EnrichGroup=="NeitherStrand"])/nPeak) %>% as.data.table() %>%
  dplyr::distinct(PeakGroup,Method,Sample,nPeak,Ratio.sames.only,Ratio.reverse.only,Ratio.both,Ratio.neither) %>% dplyr::arrange(Sample,Method,PeakGroup)


dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity  <- dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>%
                                  mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293", .default = "mESC")) %>% mutate(PeakGroup=factor(PeakGroup, levels=c("TP","FP","FN"))) %>%
  pivot_longer(cols = c("Ratio.sames.only", "Ratio.reverse.only",   "Ratio.both", "Ratio.neither"), names_to = "EnrichedStrandnessGroup",names_prefix = "Ratio.", values_to = "Ratio") %>% as.data.table() %>%
  mutate(Ratio=Ratio*100) %>% mutate(EnrichedStrandnessGroup=factor(EnrichedStrandnessGroup,levels=c("sames.only", "both", "reverse.only", "neither")))
#add grey rect
# indices to highlight
pos <- seq_along(levels(dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity$Method))
idx <- pos[pos %% 2!=0]
rects <- data.frame(xmin = pos[idx] - 0.5,xmax = pos[idx] + 0.5,ymin = -Inf,ymax = Inf)
p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN <-
  ggplot(data=dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity, aes(x=Method, y=Ratio))+
  geom_boxplot(aes(color=EnrichedStrandnessGroup),fill="transparent",position = position_dodge2(padding = 0.01,width=0.6),width = 0.55,outlier.shape = NA,size = 0.4,alpha=0.7) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(aes(color=EnrichedStrandnessGroup, shape=Cell), position = position_dodge2(padding = 0.01,width = 0.6),size=1.5, stroke = 0.6,alpha=0.7) +
  geom_rect(data = rects,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),inherit.aes = FALSE,fill = "grey92", alpha = 0.35)+
  # scales for shapes and fills
  scale_shape_manual(values = c(18,22),breaks=c("HEK293","mESC"), name = "NEB sample") +
  # scale_color_manual(values = c("#d14e51","#5f93c4","#7da14e"), breaks = c("TP","FP","FN"), name = "PeakType") +
  scale_color_manual(values = c("#eaabac","#caafd4","#7bb0d5","#cfe4b6"), breaks = c("sames.only", "both", "reverse.only", "neither"), name = "EnrichedStrandnessGroup") +
  guides(shape=guide_legend(nrow = 2),color=guide_legend(nrow=2,title = "EnrichedType"))+
  labs(x = NULL, y = str_wrap("% of m6A peak display enrichment (intensity>2)",width=40)) +
  facet_wrap(~PeakGroup)+
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85)),
    axis.text.x = element_text(face = "plain", angle = 30,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )

cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity[PeakGroup=="TP",score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity[PeakGroup=="TP",score.samestrand],method = c("pearson","spearman")[1])

#whether benchmark sites m6A level correlate better with the new intensity?
dt.Method.BestCombo.FPR20.HEK293.peaks.BenchmarkScore <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.HEK293$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.HEK293 %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,site.score=V11,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".") %>% mutate(site.score=as.numeric(site.score))
  dt.overlap <- dt.overlap %>%  group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),max.BenchScore=max(site.score,na.rm=T), mean.BenchScore=mean(site.score,na.rm=T),median.BenchScore=median(site.score,na.rm=T)) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,max.BenchScore,mean.BenchScore,median.BenchScore)
  dt.overlap %>% mutate(Sample=c)
}
dt.Method.BestCombo.FPR20.mESC.peaks.BenchmarkScore <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.mESC$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.mESC %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,site.score=V11,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".") %>% mutate(site.score=as.numeric(site.score))
  dt.overlap <- dt.overlap %>%  group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),max.BenchScore=max(site.score,na.rm=T), mean.BenchScore=mean(site.score,na.rm=T),median.BenchScore=median(site.score,na.rm=T)) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,max.BenchScore,mean.BenchScore,median.BenchScore)
  dt.overlap %>% mutate(Sample=c)
}
#incorporate the benchmark score to the peak
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore <- rbind(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% filter(str_detect(Sample,"HEK")) %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.HEK293.peaks.BenchmarkScore %>% dplyr::select(name,contains("BenchScore")),by="name"),
  dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% filter(!str_detect(Sample,"HEK")) %>%
    left_join(x=.,y=dt.Method.BestCombo.FPR20.mESC.peaks.BenchmarkScore %>% filter(str_detect(benchmark_name,"GLORI")) %>% dplyr::select(name,contains("BenchScore")),by="name"))

#calculate the cor (spearman, pearson) of score.old and score.samestrand with benchmarkscore(max,mean,median) in TP peaks
dt.cor.with.benchmarkscore <- foreach(P=c("TP","FN"), .combine='rbind')%do%{
  foreach(M = c("pearson","spearman"), .combine='rbind')%do%{
    cor.maxscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,max.BenchScore],method = M)$estimate
    cor.maxscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,max.BenchScore],method = M)$estimate
    cor.meanscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,mean.BenchScore],method = M)$estimate
    cor.meanscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,mean.BenchScore],method = M)$estimate
    cor.medianscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,median.BenchScore],method = M)$estimate
    cor.medianscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore[PeakGroup==P,median.BenchScore],method = M)$estimate
    data.table(PeakGroup=P, CorMethod=M, cor.maxscore.old, cor.maxscore.new, cor.meanscore.old, cor.meanscore.new, cor.medianscore.old, cor.medianscore.new)
  }
}
# PeakGroup CorMethod cor.maxscore.old cor.maxscore.new cor.meanscore.old cor.meanscore.new cor.medianscore.old cor.medianscore.new
# 1:        TP   pearson        0.1367976        0.3253917        0.07789068        0.03349456          0.07407888          0.02350112
# 2:        TP  spearman        0.1769708        0.4145794        0.11363164        0.08177369          0.09958120          0.05978622
# 3:        FN   pearson        0.2555909        0.3752103        0.12374715        0.14249448          0.10351791          0.11195940
# 4:        FN  spearman        0.2287484        0.4249824        0.14522571        0.22679028          0.12215359          0.17886434

dt.cor <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore  %>% dplyr::select(name,Sample,Method,PeakGroup,score.old,score.samestrand,max.BenchScore) %>%
           pivot_longer(cols = c("score.samestrand","score.old"),values_to = "score.tocompare", names_to = "MetricType", names_prefix = "score.") %>% as.data.table()
plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark <- vector(mode="list",length=4)
names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark) <- paste(rep(c("TP","FN"),each=2), rep(c("old","samestrand"),2), sep="_")
for(i in 1:length(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark)){
  P <- strsplit(names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark)[i],split="_") %>% sapply("[",1)
  M <- strsplit(names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark)[i],split="_") %>% sapply("[",2)
  dt.cor.input <- dt.cor %>% dplyr::filter(MetricType== M & PeakGroup==P)
  p <- ggplot(dt.cor.input, aes(x = score.tocompare, y = max.BenchScore)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.1,method='kde2d'),dpi=300) +
     # geom_bin2d(inherit.aes = T,bins=200) +
    viridis::scale_color_viridis(name = "Point density", option = "magma", trans = "log") +
    # geom_density2d(aes(color=MetricType))+
    stat_cor(method = "spearman",label.x=quantile(dt.cor.input$score.tocompare,0.99,na.rm=T)*0.1, size=1.8)+
    # scale_color_manual(values =c("#f0766d","#839cd1"),breaks=c("samestrand","old"))+
    coord_cartesian(xlim = c(1,quantile(dt.cor.input$score.tocompare,0.99,na.rm=T)))+
    # facet_wrap(~MetricType,nrow=1)+
    labs(x = ifelse(M=="old","Original peak score", "New peak intensity"), y = "m6A level of bechmark sites in peak", title=paste0(P," peaks")) +
    theme_minimal(base_size = 6,base_family = "Helvetica") +
    theme(legend.position = "none", plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.)), plot.subtitle = element_text(hjust = 0.5))
  if(i %in% c(2,4) ){p <- p+theme(axis.title.y = element_blank())}
  if(i %in% c(1,2)){p <- p+theme(axis.title.x= element_blank())}
  plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark[[i]] <- p
}
p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark <- cowplot::plot_grid(plotlist = plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark, nrow = 2)

#the count of FP peaks (ratio of this type FP peaks with benchmark sites) due to enrichment in the reverse strand
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>%
  filter(PeakGroup=="FP" & score.reversestrand >=2 & score.samestrand<=2 & Method %in% c("MeTPeak","exomePeak","TRESS"))
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched[,.N,by=.(Sample,Method)] %>% dplyr::arrange(Sample,Method)
#overlap with benchmark
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.HEK <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched$Sample) %>% grep(pattern="HEK",value=T), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched %>% dplyr::filter(Sample==c) %>% mutate(seqnames=strsplit(name,split="_") %>% sapply("[",5),
                                                                                                       start=strsplit(name,split="_") %>% sapply("[",6),
                                                                                                       end=strsplit(name,split="_") %>% sapply("[",7),score=score.new,strand=strand.new) %>%
    dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.mESC <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched$Sample) %>% grep(pattern="mESC",value=T), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched %>% dplyr::filter(Sample==c) %>% mutate(seqnames=strsplit(name,split="_") %>% sapply("[",4),
                                                                                                       start=strsplit(name,split="_") %>% sapply("[",5),
                                                                                                       end=strsplit(name,split="_") %>% sapply("[",6),score=score.new,strand=strand.new) %>%
    dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"mESC")){
    bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.mESC <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.mESC %>% group_by(name) %>%
   dplyr::arrange(desc(nConfidentBenchSite),desc(nBenchSite)) %>% slice_head(n=1) %>% as.data.table()

dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap <- rbind(dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched %>% dplyr::filter() %>% dplyr::filter(str_detect(Sample,"mESC")) %>%
                                                                             left_join(x=.,y=dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.mESC %>% dplyr::select(name,nBenchSite,nConfidentBenchSite),by="name"),
                                                                           dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched %>% dplyr::filter(str_detect(Sample,"HEK")) %>%
                                                                             left_join(x=.,y=dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.HEK %>% dplyr::select(name,nBenchSite,nConfidentBenchSite),by="name")
                                                                           ) %>% mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite), nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))

#calculate the ratio of those FP reverse
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap %>% group_by(Sample,Method) %>% mutate(nPeak=n_distinct(name)) %>%
  mutate(ratio.Benchmark=n_distinct(name[nBenchSite>0])/nPeak, ratio.ConfidentBenchmark = n_distinct(name[nConfidentBenchSite>0])/nPeak) %>%
  as.data.table() %>% distinct(Sample,Method,nPeak,ratio.Benchmark,ratio.ConfidentBenchmark) %>% dplyr::arrange(Sample,Method)
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio %>% group_by(Method) %>%
  mutate(avg.ratio.Benchmark=mean(ratio.Benchmark), avg.ratio.ConfidentBenchmark=mean(ratio.ConfidentBenchmark)) %>% as.data.table() %>% distinct(Method,avg.ratio.Benchmark,avg.ratio.ConfidentBenchmark) %>%
  pivot_longer(cols = c("avg.ratio.Benchmark", "avg.ratio.ConfidentBenchmark"),names_to = "ContainSite",values_to = "Ratio",names_prefix = "avg.ratio.") %>% as.data.table() %>%
  mutate(ContainSite=factor(ContainSite,levels=c("Benchmark","ConfidentBenchmark"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>% mutate(label.pt=paste0(round(Ratio*100,1),"%")) %>%
  mutate(Method=factor(Method,levels=unique(Method)))

p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio <-
  ggplot(data=dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio,
         aes(x=Method,y=Ratio, fill=ContainSite))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio %>% dplyr::filter(ContainSite=="Benchmark m6As"),
            aes(label = label.pt, y =  Ratio+0.02, x=Method),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio %>% dplyr::filter(ContainSite=="Confident benchmark m6As"),
            aes(label = label.pt, y =  Ratio+0.02, x=Method),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y=str_wrap("% of FP peaks with intensity(>2) at reverse strand",width=40))+
  guides(fill=guide_legend(title = expression("Peaks\ncontain"),nrow = 2))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  # scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio
#3c The distribution of AUC of six method's  all combo with all alternative options (NEB_HEK293 and NEB_mESC)
#load all combo option AUC
dt.pruned.FPR.TPR.option.AUC.NEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HEK.RDS")
dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% distinct(Method,ComboName,optionID,AUC.scaled)
dt.pruned.FPR.TPR.option.AUC.NEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.mESC.RDS")
dt.pruned.FPR.TPR.option.AUC.NEB.mESC %>% distinct(Method,ComboName,optionID,AUC.scaled)
#load all combo AUC
dt.pruned.FPR.TPR.combo.AUC.NEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HEK.all.combo.RDS")
dt.pruned.FPR.TPR.combo.AUC.NEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.mESC.all.combo.RDS")
dt.BestCombo.NEB <- dt.pruned.FPR.TPR.combo.AUC.NEB.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.DefaultCombo.AUC.NEB <- rbind(dt.pruned.FPR.TPR.combo.AUC.NEB.HEK %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) %>% mutate(Cell="HEK293"),
                                dt.pruned.FPR.TPR.combo.AUC.NEB.mESC %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) %>% mutate(Cell="mESC"))
dt.BestCombo.AUC.NEB <- rbind(dt.pruned.FPR.TPR.combo.AUC.NEB.HEK %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.NEB %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) %>% mutate(Cell="HEK293"),
                              dt.pruned.FPR.TPR.combo.AUC.NEB.mESC %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.NEB %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) %>% mutate(Cell="mESC"))

dt.AUC.NEB.all.option <- rbind(dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>% mutate(Cell="HEK293"),
                              dt.pruned.FPR.TPR.option.AUC.NEB.mESC %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>% mutate(Cell="mESC"))

dt.BestComboOption.NEB <- dt.AUC.NEB.all.method.combo %>% filter(Cell=="HEK293") %>% group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>%
  dplyr::arrange(desc(AUC.scaled),desc(TPR20))

dt.AUC.NEB.all.option <- dt.AUC.NEB.all.option %>% left_join(x=.,y=dt.BestComboOption.NEB %>% mutate(IsBestComboOption=TRUE) %>% dplyr::select(Method,ComboName,optionID,IsBestComboOption),by=c("Method","ComboName","optionID")) %>%
                         mutate(ComboGroup=case_when(IsBestComboOption==TRUE ~ "BestComboOption", .default = "OthersComboOption")) %>% dplyr::select(-IsBestComboOption)
dt.AUC.NEB.all.option <- rbind(dt.AUC.NEB.all.option %>% mutate(Group=ComboGroup) %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group),
                               dt.DefaultCombo.AUC.NEB %>% mutate(optionID=NA, Group="DefaultCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group),
                               dt.BestCombo.AUC.NEB %>% mutate(optionID=NA, Group="BestCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group))
dt.AUC.NEB.all.option <- dt.AUC.NEB.all.option %>% mutate(Method=factor(Method,levels=unique(dt.BestComboOption.NEB$Method)),Group=factor(Group,levels=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo")))
p.AUC.NEB.all.option <- ggplot(dt.AUC.NEB.all.option, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.AUC.NEB.all.option %>% dplyr::filter(Group =="OthersComboOption"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=0.1,alpha=0.15, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.AUC.NEB.all.option %>% dplyr::filter(Group =="DefaultCombo"), aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=1.5)+
  geom_jitter(data = dt.AUC.NEB.all.option %>% dplyr::filter(Group =="BestCombo") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=1.5)+
  geom_point(data = dt.AUC.NEB.all.option %>% dplyr::filter(Group =="BestComboOption") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#eb4601",size=3)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled)", title = paste0("9120 (285 Combo x 32 Option) ComboOption combination")) +
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 30,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = c(0.92, 0.85),  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )
# legend-only plot for color: hide linetype legend
p_col_legend <- ggplot(data=dt.AUC.NEB.all.option,aes(x = Method, y = AUC.scaled)) +
  geom_jitter(aes(shape=Group,color=Group,size=Group), width = 0.18, height = 0,  stroke = 0.6,alpha=1) +
  scale_shape_manual(values = c(18,22,18,18),breaks=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  scale_size_manual(values = c(3,0.5,2,2), breaks = c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  scale_color_manual(values = c("#eb4601","grey70","#e9abac","#7bb0d5"), breaks = c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  theme_minimal() +
  guides(linetype = "none") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "right",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"), legend.title = element_blank(),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
# extract legend grobs
leg_col <- cowplot::get_legend(p_col_legend)
# assemble: main plot + legends at arbitrary positions
p.AUC.NEB.all.option  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.AUC.NEB.all.option , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1.15, y = 0.8, width = 0.22, height = 0.25, hjust = 1.2, vjust = 0.5)


p.AUC.NEB.all.option

pdf("new_Fig3G.AUC.NEB.all.option.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.AUC.NEB.all.option, x = 0.05, y=9.2, default.units = "cm",width = 10, height = 3.8)
dev.off()


#3d The AUC of six method's best combo and option (NEB_HEK293 and NEB_mESC)
dt.pruned.TPR.FPR.NEB.Method.BestOption <- rbind(dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% inner_join(x=.,y=dt.BestComboOption.NEB %>% dplyr::select(Method,ComboName,optionID),by=c("Method","ComboName","optionID")) %>%
                                        filter(str_detect(Benchmark,"GLORI")) %>% mutate(Cell="HEK293"),
                               dt.pruned.FPR.TPR.option.AUC.NEB.mESC %>% inner_join(x=.,y=dt.BestComboOption.NEB %>% dplyr::select(Method,ComboName,optionID),by=c("Method","ComboName","optionID")) %>%
                                 filter(str_detect(Benchmark,"GLORI")) %>%  mutate(Cell="mESC"))

dt.pruned.TPR.FPR.NEB.Method.BestOption.label <- dt.pruned.TPR.FPR.NEB.Method.BestOption %>% distinct(Method,ComboName,optionID,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>%
  mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.pruned.TPR.FPR.NEB.Method.BestOption$TPR),2)) %>% mutate(xpos=0.6)
p.pruned.TPR.FPR.NEB.Method.BestOption  <-  ggplot(data=dt.pruned.TPR.FPR.NEB.Method.BestOption ,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of m6A peak from best combo and option combination",y=str_wrap("Surrogate TPR of m6A peak from best combo and option combination",width = 40))+
  geom_text(data=dt.pruned.TPR.FPR.NEB.Method.BestOption.label, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  facet_wrap(~Cell)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.minor.ticks.x.bottom = element_blank(),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.pruned.TPR.FPR.NEB.Method.BestOption

########################## 3e the novel m6A+ gene in MACS2S vs exomePeak2Best at (FPR20 cutoff) in NEB_HEK and NEB_mESC ##############################

#calling  m6A peaks from MACS2S BestOption at FPR20 cutoff
dt.pruned.TPR.FPR.NEB.Method.BestOption %>% filter(Method=="MACS2" & FPR==0.2)
#    FPR  TPR cutoff Method ComboName         optionID                          Benchmark  AUC AUC.scaled maxTPRmFPR TPR20   Cell
# 1: 0.2 0.65   1.48  MACS2    Combo1 mean_Tx_bin50_P2 dt.benchmark.GLORI_HEK293.m6As.bed 0.31       0.62       0.46  0.65 HEK293
# 2: 0.2 0.64   1.86  MACS2    Combo1 mean_Tx_bin50_P2   dt.benchmark.GLORI_mESC.m6As.bed 0.30       0.60       0.46  0.64   mESC
dt.TPR20.MACS2.BestOption.cutoff <- dt.pruned.TPR.FPR.NEB.Method.BestOption %>% filter(Method=="MACS2" & FPR==0.2 & str_detect(Benchmark,"GLORI_HEK293")) %>% dplyr::select(Method,ComboName,optionID,cutoff)

#obtain  MACS2    Combo1 mean_Tx_bin50_P2 for HEK293
#run bedtools intensity
#HEK293
#prepare the stranded.BAM.Depth
Samples <- unique(dt.Method.BestCombo.FPR20.m6A.HEK293$Sample) %>% grep(pattern="HEK",value=T)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/Figure3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID
selected.option <- dt.TPR20.MACS2.BestOption.cutoff$optionID
#load the old method m6A
dt.Method.m6a <- LoadPeakMethod(Method=dt.TPR20.MACS2.BestOption.cutoff$Method,
                                Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB/",dt.TPR20.MACS2.BestOption.cutoff$Method,"_",dt.TPR20.MACS2.BestOption.cutoff$ComboName)) %>%
  mutate(Method=dt.TPR20.MACS2.BestOption.cutoff$Method, ComboName=dt.TPR20.MACS2.BestOption.cutoff$ComboName)

dt.Method.m6a <- dt.Method.m6a %>% filter(str_detect(seqnames,"chr"))  %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf"
org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite"
org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed"
org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/hg38_gencode.gene.bed"
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo=dt.Method.m6a[Method==M,unique(ComboName)], .combine='rbind')%do%{
    foreach(s=Samples, .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      # selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",
                 intron.bed = "/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
  }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
# intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)

#obtain m6a at FPR cutoff
dt.method.combo.selected <- dt.TPR20.MACS2.BestOption.cutoff
dt.MACS2.BestOption.m6A.FPR20.HEK293 <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
  foreach(combo = dt.method.combo.selected[Method==M,unique(ComboName)], .combine='rbind')%do%{
    # message(paste0("Obtain peak at FPR cutoff at ",ifelse(is.null(Model.FPR.cutoff), "default", Model.FPR.cutoff)," for ", M, " ", combo, " ", dt.method.combo.selected[Method==M & ComboName==combo,optionID]))
    intensity.files <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = T) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T)
    names(intensity.files) <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = F) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T) %>%
      strsplit(split="__",fixed=T) %>% sapply("[",1)
    intensity.files <- intensity.files[Samples]
    dt.intensity <- foreach(s=names(intensity.files),.combine='rbind')%do%{data.table::fread(file = intensity.files[s]) %>% mutate(Sample=s)}#load intensity for all sample
    dt.intensity <- dt.intensity %>% dplyr::select(name,Sample,pos,neg)
    #pull out m6a with intensity higher than cutoff at FPR20
    #filtered and un overlapped  m6a at pos strand
    dt.m6a.filtered.pos <- dt.intensity %>% dplyr::filter(pos >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.pos)>0){
      dt.m6a.filtered.pos <- dt.m6a.filtered.pos %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="+") %>%
        dplyr::select(seqnames,start,end,name,score=pos,strand,Sample)
      dt.m6a.filtered.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.pos,keep.extra.columns = T)
      dt.m6a.filtered.pos.merged <- foreach(sam = unique(dt.m6a.filtered.pos$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.pos.gr[dt.m6a.filtered.pos.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.pos.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    #filtered  and un overlapped m6a at neg strand
    dt.m6a.filtered.neg <- dt.intensity %>% dplyr::filter(neg >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.neg)>0){
      dt.m6a.filtered.neg <- dt.m6a.filtered.neg %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="-") %>%
        dplyr::select(seqnames,start,end,name,score=neg,strand,Sample)
      dt.m6a.filtered.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.neg,keep.extra.columns = T)
      dt.m6a.filtered.neg.merged <- foreach(sam = unique(dt.m6a.filtered.neg$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.neg.gr[dt.m6a.filtered.neg.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.neg.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    dt.m6a.filtered <- rbind(dt.m6a.filtered.pos.merged, dt.m6a.filtered.neg.merged)
    if(nrow(dt.m6a.filtered)>0){
      #calculate exonic sum width to filter extreme long peak (<=6kb)
      dt.m6a.filtered.sumwidth <- LongPeakRemoveIntron(dt.peak = dt.m6a.filtered %>% distinct(seqnames,start,end,name,score,strand),#unique m6a peak
                                                       dt.intron = data.table::fread(org.intron.bed),
                                                       genome_file = org.genome_file)
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.sumwidth  %>% dplyr::distinct(name,SumWidth), by=c("name")) %>%
        dplyr::mutate(SumWidth=case_when(is.na(SumWidth) ~ 0, .default = SumWidth))
      #the overlap gene
      dt.m6a.filtered.overlapped.gene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=data.table::fread(org.genebed),
                                                                 s=T, f=0.5, F=0.8, c=T,e = T) %>% as.data.table()
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.overlapped.gene %>% dplyr::select(name=V4,Sample=V7,OverlappedNoGene=V9) %>%
                                                         dplyr::distinct(Sample,name,OverlappedNoGene), by=c("Sample","name"))
      dt.m6a.filtered %>% mutate(Method=M, ComboName=combo, optionID=dt.method.combo.selected[Method==M & ComboName==combo,optionID])
    }
  }
}
dt.MACS2.BestOption.m6A.FPR20.HEK293  <- dt.MACS2.BestOption.m6A.FPR20.HEK293 %>% mutate(name=paste(Method,ComboName,optionID,Sample,name,sep="_"))


#mESC
#prepare the stranded.BAM.Depth
Samples <- unique(dt.Method.BestCombo.FPR20.m6A.mESC$Sample)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/Figure3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID
selected.option <- dt.TPR20.MACS2.BestOption.cutoff$optionID
#load the old method m6A
dt.Method.m6a <- LoadPeakMethod(Method=dt.TPR20.MACS2.BestOption.cutoff$Method,
                                Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB_mESC/",dt.TPR20.MACS2.BestOption.cutoff$Method,"_",dt.TPR20.MACS2.BestOption.cutoff$ComboName)) %>%
  mutate(Method=dt.TPR20.MACS2.BestOption.cutoff$Method, ComboName=dt.TPR20.MACS2.BestOption.cutoff$ComboName)

dt.Method.m6a <- dt.Method.m6a %>% filter(str_detect(seqnames,"chr"))  %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf"
org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.mm39.knownGene.sqlite"
org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_mm39.txt"
org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.mm39.bed"
org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/mm39_gencode.gene.bed"
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo=dt.Method.m6a[Method==M,unique(ComboName)], .combine='rbind')%do%{
    foreach(s=Samples, .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      # selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = org.genome_file,
                 intron.bed = org.intron.bed,
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
  }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
# intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)

#obtain m6a at FPR cutoff
dt.method.combo.selected <- dt.TPR20.MACS2.BestOption.cutoff
dt.MACS2.BestOption.m6A.FPR20.mESC <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
  foreach(combo = dt.method.combo.selected[Method==M,unique(ComboName)], .combine='rbind')%do%{
    # message(paste0("Obtain peak at FPR cutoff at ",ifelse(is.null(Model.FPR.cutoff), "default", Model.FPR.cutoff)," for ", M, " ", combo, " ", dt.method.combo.selected[Method==M & ComboName==combo,optionID]))
    intensity.files <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = T) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T)
    names(intensity.files) <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = F) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T) %>%
      strsplit(split="__",fixed=T) %>% sapply("[",1)
    intensity.files <- intensity.files[Samples]
    dt.intensity <- foreach(s=names(intensity.files),.combine='rbind')%do%{data.table::fread(file = intensity.files[s]) %>% mutate(Sample=s)}#load intensity for all sample
    dt.intensity <- dt.intensity %>% dplyr::select(name,Sample,pos,neg)
    #pull out m6a with intensity higher than cutoff at FPR20
    #filtered and un overlapped  m6a at pos strand
    dt.m6a.filtered.pos <- dt.intensity %>% dplyr::filter(pos >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.pos)>0){
      dt.m6a.filtered.pos <- dt.m6a.filtered.pos %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="+") %>%
        dplyr::select(seqnames,start,end,name,score=pos,strand,Sample)
      dt.m6a.filtered.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.pos,keep.extra.columns = T)
      dt.m6a.filtered.pos.merged <- foreach(sam = unique(dt.m6a.filtered.pos$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.pos.gr[dt.m6a.filtered.pos.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.pos.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    #filtered  and un overlapped m6a at neg strand
    dt.m6a.filtered.neg <- dt.intensity %>% dplyr::filter(neg >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.neg)>0){
      dt.m6a.filtered.neg <- dt.m6a.filtered.neg %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="-") %>%
        dplyr::select(seqnames,start,end,name,score=neg,strand,Sample)
      dt.m6a.filtered.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.neg,keep.extra.columns = T)
      dt.m6a.filtered.neg.merged <- foreach(sam = unique(dt.m6a.filtered.neg$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.neg.gr[dt.m6a.filtered.neg.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.neg.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    dt.m6a.filtered <- rbind(dt.m6a.filtered.pos.merged, dt.m6a.filtered.neg.merged)
    if(nrow(dt.m6a.filtered)>0){
      #calculate exonic sum width to filter extreme long peak (<=6kb)
      dt.m6a.filtered.sumwidth <- LongPeakRemoveIntron(dt.peak = dt.m6a.filtered %>% distinct(seqnames,start,end,name,score,strand),#unique m6a peak
                                                       dt.intron = data.table::fread(org.intron.bed),
                                                       genome_file = org.genome_file)
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.sumwidth  %>% dplyr::distinct(name,SumWidth), by=c("name")) %>%
        dplyr::mutate(SumWidth=case_when(is.na(SumWidth) ~ 0, .default = SumWidth))
      #the overlap gene
      dt.m6a.filtered.overlapped.gene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=data.table::fread(org.genebed),
                                                                 s=T, f=0.5, F=0.8, c=T,e = T) %>% as.data.table()
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.overlapped.gene %>% dplyr::select(name=V4,Sample=V7,OverlappedNoGene=V9) %>%
                                                         dplyr::distinct(Sample,name,OverlappedNoGene), by=c("Sample","name"))
      dt.m6a.filtered %>% mutate(Method=M, ComboName=combo, optionID=dt.method.combo.selected[Method==M & ComboName==combo,optionID])
    }
  }
}
dt.MACS2.BestOption.m6A.FPR20.mESC  <- dt.MACS2.BestOption.m6A.FPR20.mESC %>% mutate(name=paste(Method,ComboName,optionID,Sample,name,sep="_"))

##pull out the m6A peaks from exomePeak2 BestCombo at FPR20 cutoff
dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks <- readRDS("dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks.RDS")

#check the overlap between MACS2 best option peaks with exomePeak2 best combo peaks
dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap <- foreach(c = unique(dt.MACS2.BestOption.m6A.FPR20.HEK293$Sample), .combine='rbind')%do%{
  bed.peak <- dt.MACS2.BestOption.m6A.FPR20.HEK293 %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.exomePeak2 <- dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.exomePeak2, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
  dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(exomePeak2overlap=ifelse(sum.overlap/width>=0.0001,"exomePeak2+","exomePeak2-"))
  bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,exomePeak2overlap), by="name") %>% mutate(Sample=c)
}
dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap[,.N,by=.(Sample,exomePeak2overlap)] %>% dplyr::arrange(Sample,exomePeak2overlap)
#           Sample exomePeak2overlap     N
# 1: HEK_NEB_mRNA1       exomePeak2+ 10526
# 2: HEK_NEB_mRNA1       exomePeak2-  2733
# 3: HEK_NEB_mRNA2       exomePeak2+ 14313
# 4: HEK_NEB_mRNA2       exomePeak2-  2746

dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap <- foreach(c = unique(dt.MACS2.BestOption.m6A.FPR20.mESC$Sample), .combine='rbind')%do%{
  bed.peak <- dt.MACS2.BestOption.m6A.FPR20.mESC %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.exomePeak2 <- dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.exomePeak2, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
  dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(exomePeak2overlap=ifelse(sum.overlap/width>=0.0001,"exomePeak2+","exomePeak2-"))
  bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,exomePeak2overlap), by="name") %>% mutate(Sample=c)
}
dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap[,.N,by=.(Sample,exomePeak2overlap)] %>% dplyr::arrange(Sample,exomePeak2overlap)
#     Sample exomePeak2overlap     N
# 1: mESC_WT1       exomePeak2+ 10375
# 2: mESC_WT1       exomePeak2-  3505
# 3: mESC_WT2       exomePeak2+ 11091
# 4: mESC_WT2       exomePeak2-  3618

#check the exomePeak2- and exomePeak2+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
#HEK293
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
dt.MACS2.BestOption.HEK293.peaks.Benchmarkoverlap <- foreach(c = unique(dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap$Sample), .combine='rbind')%do%{
  bed.peak <- dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#obtain the MACS2 Best option only m6A+ gene
dt.MACS2.BestOption.HEK293.peaks.annotgene <- annot_peak(peak.bed=dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                             strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                             annot_type = "gene")
dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap <- dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap %>% left_join(x=.,y=dt.MACS2.BestOption.HEK293.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
dt.MACS2.BestOption.specific.m6Agene.HEK293 <-   dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestOption.specific=all(exomePeak2overlap=="exomePeak2-")) %>% as.data.table() %>%#if all peak of one gene is exomePeak2-, then this gene is BestOption specific m6A+ gene
  dplyr::filter(BestOption.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
#obtain the BestOption specific peaks overlapped with (high confident) benchmark sites
dt.MACS2.BestOption.specific.peak.benchmark.HEK293 <-  dt.MACS2.BestOption.specific.m6Agene.HEK293 %>% filter(exomePeak2overlap=="exomePeak2-" & BestOption.specific==TRUE) %>%
  left_join(x=.,dt.MACS2.BestOption.HEK293.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
  dplyr::filter(nBenchSite>0)
#rank the BestCombo specific m6A+ genes with benchmark m6A sites support
dt.MACS2.BestOption.specific.m6Agene.benchmark.HEK293 <- dt.MACS2.BestOption.specific.peak.benchmark.HEK293 %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
  group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
  as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))

#mESC
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
dt.MACS2.BestOption.mESC.peaks.Benchmarkoverlap <- foreach(c = unique(dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap$Sample), .combine='rbind')%do%{
  bed.peak <- dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"mESC")){
    bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=F,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#obtain the MACS2 Best option only m6A+ gene
dt.MACS2.BestOption.mESC.peaks.annotgene <- annot_peak(peak.bed=dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                         strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                                         annot_type = "gene")
dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap <- dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap %>% left_join(x=.,y=dt.MACS2.BestOption.mESC.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
dt.MACS2.BestOption.specific.m6Agene.mESC <-   dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestOption.specific=all(exomePeak2overlap=="exomePeak2-")) %>% as.data.table() %>%#if all peak of one gene is exomePeak2-, then this gene is BestOption specific m6A+ gene
  dplyr::filter(BestOption.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
#obtain the BestOption specific peaks overlapped with (high confident) benchmark sites
dt.MACS2.BestOption.specific.peak.benchmark.mESC <-  dt.MACS2.BestOption.specific.m6Agene.mESC %>% filter(exomePeak2overlap=="exomePeak2-" & BestOption.specific==TRUE) %>%
  left_join(x=.,dt.MACS2.BestOption.mESC.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
  dplyr::filter(nBenchSite>0)
#rank the BestCombo specific m6A+ genes with benchmark m6A sites support
dt.MACS2.BestOption.specific.m6Agene.benchmark.mESC <- dt.MACS2.BestOption.specific.peak.benchmark.mESC %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
  group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
  as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))


#visualization of ratio of exomePeak2 ratio of MACS2 BestOption peaks in HEK293 and mESC
dt.MACS2.BestOption.specific.ratio.m6A.peak <- rbind(dt.MACS2.BestOption.m6A.FPR20.mESC.exomePeak2Bestoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestOptionOnly = n_distinct(name[exomePeak2overlap=="exomePeak2-"])/n_distinct(name),
                                                                                                                                             nPeak.BestOptionOnly= n_distinct(name[exomePeak2overlap=="exomePeak2-"]), nPeak.BestOption=n_distinct(name)) %>%
                                                           as.data.table() %>% distinct(Sample,Ratio.Peak.BestOptionOnly,nPeak.BestOptionOnly,nPeak.BestOption),
                                                     dt.MACS2.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestOptionOnly = n_distinct(name[exomePeak2overlap=="exomePeak2-"])/n_distinct(name),
                                                                                                                                              nPeak.BestOptionOnly= n_distinct(name[exomePeak2overlap=="exomePeak2-"]), nPeak.BestOption=n_distinct(name)) %>%
                                                       as.data.table() %>% distinct(Sample,Ratio.Peak.BestOptionOnly,nPeak.BestOptionOnly,nPeak.BestOption))

dt.MACS2.BestOption.specific.ratio.m6A.peak <- dt.MACS2.BestOption.specific.ratio.m6A.peak %>% mutate(label1=paste0(nPeak.BestOptionOnly,"/",nPeak.BestOption), label2=paste0(round(Ratio.Peak.BestOptionOnly*100,1),"%")) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_NEB_mRNA1","HEK_NEB_mRNA2","mESC_WT1","mESC_WT2"), labels=c("NEB_HEK293_mRNA1","NEB_HEK293_mRNA2","NEB_mESC_mRNA1","NEB_mESC_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293",.default = "mESC"))
p.MACS2.BestOption.specific.ratio.m6A.peak <-
  ggplot(dt.MACS2.BestOption.specific.ratio.m6A.peak, aes(y = Sample, x = Ratio.Peak.BestOptionOnly, fill = Cell)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  # percent label to the right of the bar
  geom_text(aes(label = label2, x =  Ratio.Peak.BestOptionOnly + 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  #N peak label to the left of the bar
  geom_text(aes(label = label1, x = Ratio.Peak.BestOptionOnly/2),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  scale_fill_manual(values = c("#7b93c6","#acd9ee"),breaks = c("HEK293","mESC")) +
  labs(x = "% of novel m6A peaks in MACS2S v.s exomePeak2 best combo", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

#visualization of count of BestCombo only m6A+ gene with benchmark m6As and highly confident m6As
dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark <- rbind(dt.MACS2.BestOption.specific.peak.benchmark.HEK293 %>% group_by(Sample) %>%
                                                                         mutate(No.BestOptionOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                         distinct(Sample,No.BestOptionOnly.m6AGene.BenchmarkSite,No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite),
                                                                       dt.MACS2.BestOption.specific.peak.benchmark.mESC %>% group_by(Sample) %>%
                                                                         mutate(No.BestOptionOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                         distinct(Sample,No.BestOptionOnly.m6AGene.BenchmarkSite,No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite)
)
dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark <- dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark %>%
  pivot_longer(cols = c("No.BestOptionOnly.m6AGene.BenchmarkSite", "No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite"),values_to = "No.BestOptionOnly.m6AGene", names_to = "GeneGroup",names_prefix = "No.BestOptionOnly.m6AGene.") %>%
  as.data.table() %>% mutate(GeneGroup=factor(GeneGroup, levels=c("BenchmarkSite","ConfidentBenchmarkSite"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_NEB_mRNA1","HEK_NEB_mRNA2","mESC_WT1","mESC_WT2"), labels=c("NEB_HEK293_mRNA1","NEB_HEK293_mRNA2","NEB_mESC_mRNA1","NEB_mESC_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293",.default = "mESC"))

dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark <- dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark %>% mutate(label.count=paste0(paste0("No.Gene=",No.BestOptionOnly.m6AGene)))
p.count.MACS2.BestOption.specific.m6Agene.with.benchmark  <- ggplot(data=dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark, aes(x=Sample,y=No.BestOptionOnly.m6AGene, fill=GeneGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark %>% dplyr::filter(GeneGroup=="Benchmark m6As"),
            aes(label = label.count, y =  No.BestOptionOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 1.8, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.MACS2.BestOption.specific.m6Agene.with.benchmark %>% dplyr::filter(GeneGroup=="Confident benchmark m6As"),
            aes(label = label.count, y =  No.BestOptionOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 1.8, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y=str_wrap("N of novel m6A+ genes in MACS2S v.s exomePeak2 best combo",width=40))+
  guides(fill=guide_legend(title = "Peaks contain"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.MACS2.BestOption.specific.m6Agene.with.benchmark
#visualization of IGV plot of MACS2 Best Option specific m6A+ gene in HEK293 and mESC

dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

selected.Novel.m6Agene.NEB.HEK <-  dt.MACS2.BestOption.specific.m6Agene.benchmark.HEK293  %>% dplyr::filter(nRep==2) %>%
  distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=5)
dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==c("DHCR7","EEF1D","ZER1")[2])
StudySamples <- c("HEK_NEB_mRNA1", "HEK_NEB_mRNA2")
## Create a page (7.5*7.5cm)
window.size=10#100kb
pseudocount=1
track.height=0.5
pdf("Fig3_NEB_HEK293_EEF1D_MACS2_BestOption_specific_peak.pdf")

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
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

# #add peak in this region
dt.peak.to.see <- rbind(dt.MACS2.BestOption.m6A.FPR20.HEK293 %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="MACS2_BestOption"),
                        dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="exomePeak2_BestCombo")) %>%
                mutate(Group=factor(Group, levels=c("MACS2_BestOption","exomePeak2_BestCombo"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)

for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  exomePeak2.bestcombo <- plotRanges(data = dt.peak.region, params = region, order = "random",
                                     fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),strandSplit = F,
                                     x = 0.5, y = 1.5+(s-1+length(StudySamples))*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm"
  )
  sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5+length(StudySamples))*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
  sample_rect <- plotRect(x=0.5,y=1.5+(s-1+length(StudySamples))*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#track for benchmark m6A sites (perturbation effect)
#FTO GLORI HEK293
signal_confident.benchmark.FTO <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.FTO), params = region,
                                             x = 0.5, y = 1.5+2*length(StudySamples)*track.height, width = 5, height = track.height,
                                             just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="FTO treatment", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+2*length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
#METTL3i GLORI HEK293
signal_confident.benchmark.METTL3i <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.METTL3i), params = region,
                                                 x = 0.5, y = 1.5+(2*length(StudySamples)+1)*track.height, width = 5, height = track.height,
                                                 just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+1.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="METTL3 inhibition", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+1.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fontface="plain")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in GLORI_HEK293"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fontface="plain")
legendPlot.exomePeak2.peak <- plotLegend(legend = c("MACS2_BestOption","exomePeak2_BestCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 2.5, height = 1,
                                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                         title = expression("PeakGroup"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")
dev.off()

#mESC exomePeak2 best combo specific m6A+ gene with high confident benchmark m6As
#select well-studied genes
selected.Novel.m6Agene.NEB.mESC <-  dt.MACS2.BestOption.specific.m6Agene.benchmark.mESC  %>% dplyr::filter(nRep==2) %>%
  distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"Gm")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=10)
dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==c("Cdh1","Ubr7","Mtdh","Zfp207","Golt1b","Rdx")[5])
StudySamples <- c("mESC_WT1","mESC_WT2")
library(BSgenome.Mmusculus.UCSC.mm39)

## Create a page (7.5*7.5cm)
window.size=50#100kb
pseudocount=1
track.height=0.5
pdf("Fig3_NEB_mESC_Golt1b_exomePeak2_bestcombo_specific_peak.pdf")
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
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
# #add peak in this region
dt.peak.to.see <- rbind(dt.MACS2.BestOption.m6A.FPR20.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="MACS2_BestOption"),
                        dt.exomePeak2.BestCombo.HEK.mESC.FPR20.peaks %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="exomePeak2_BestCombo")) %>%
  mutate(Group=factor(Group, levels=c("MACS2_BestOption","exomePeak2_BestCombo"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)

for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  exomePeak2.bestcombo <- plotRanges(data = dt.peak.region, params = region, order = "random",
                                     fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),strandSplit = F,
                                     x = 0.5, y = 1.5+(s-1+length(StudySamples))*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm"
  )
  sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5+length(StudySamples))*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
  sample_rect <- plotRect(x=0.5,y=1.5+(s-1+length(StudySamples))*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#track for benchmark m6A sites (perturbation effect)
#METTL3cko
signal_confident.benchmark.Mettl3cko <- plotSignal(data = dt.eTAMseq.mESC.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect), params = region,
                                                   x = 0.5, y = 1.5+2*length(StudySamples)*track.height, width = 5, height = track.height,
                                                   just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="Mettl3_CKO", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+2*length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(2*length(StudySamples)+1)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in eTAMseq_mESC"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+1)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fondface="bold")
legendPlot.exomePeak2.peak <- plotLegend(legend = c("MACS2_BestOption","exomePeak2_BestCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                         title = expression("PeakType"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#save intermediate variables for re-use
save.image("sm6APeak_Figure3.intermediate.results.RDS")
# q(save="no")
####################### combine all figure 3 together #############################
fig3.list <- list(p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count=p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count,
                  p.bothstrand.Method.BestCombo.FPR20.FP=p.bothstrand.Method.BestCombo.FPR20.FP,
                  p.bothstrand.Method.BestCombo.FPR20.FN=p.bothstrand.Method.BestCombo.FPR20.FN,
                  p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN=p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN,
                  p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio=p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio,
                  p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark= p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark,#this plot need to be output in raster pdf
                  p.AUC.NEB.all.option=p.AUC.NEB.all.option,
                  p.pruned.TPR.FPR.NEB.Method.BestOption=p.pruned.TPR.FPR.NEB.Method.BestOption,
                  p.MACS2.BestOption.specific.ratio.m6A.peak=p.MACS2.BestOption.specific.ratio.m6A.peak,
                  p.count.MACS2.BestOption.specific.m6Agene.with.benchmark=p.count.MACS2.BestOption.specific.m6Agene.with.benchmark
                  )
saveRDS(fig3.list, file="Figure3.plot.list.RDS")


ggsave(plot=p.bothstrand.Method.BestCombo.FPR20.FP, dpi = 300,units = "cm",device = cairo_pdf,
        width=5.8,height = 4.5,filename = "Figure3_p.bothstrand.Method.BestCombo.FPR20.FP.pdf")
ggsave(plot=p.bothstrand.Method.BestCombo.FPR20.FN, dpi = 300,units = "cm",device = cairo_pdf,
        width=5.8,height = 4.5,filename = "Figure3_p.bothstrand.Method.BestCombo.FPR20.FN.pdf")

ggsave(plot=p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 4.5, filename = "Figure3_Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.pdf")


pdf("Figure3.Stranded_intensity_facilitate more_accurate_and_sensitive_peak_calling.pdf",width = 8.2677, height = 11.693)

pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count, x = 0.05, y=0.2, default.units = "cm",width = 6, height = 4.5)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 6.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.bothstrand.Method.BestCombo.FPR20.FP, x = 7, y=0.2, default.units = "cm",width = 6.3, height = 4.5)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 12.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.bothstrand.Method.BestCombo.FPR20.FN, x = 12.5, y=0.2, default.units = "cm",width = 6.3, height = 4.5)
#row2
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN, x = 0.05, y=4.7, default.units = "cm",width = 8, height = 4.5)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 8.2, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio, x = 8.2, y=4.7, default.units = "cm",width = 4, height = 4.5)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 12.2, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark, x = 14, y=5, default.units = "cm",width = 5.8, height = 4.5)
#row3
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05, y = 9.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.AUC.NEB.all.option, x = 0.05, y=9.2, default.units = "cm",width = 10, height = 3.8)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 10.5, y = 9.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.MACS2.BestOption.specific.ratio.m6A.peak, x = 10.5, y=9.2, default.units = "cm",width = 5.8, height = 3.8)
#row4
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05, y = 13.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.pruned.TPR.FPR.NEB.Method.BestOption, x = 0.05, y=13.2, default.units = "cm",width = 10, height = 4.5)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 10.2, y = 13.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.MACS2.BestOption.specific.m6Agene.with.benchmark, x = 10.2, y=13.2, default.units = "cm",width = 6.2, height = 4.5)
#row5
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05, y = 17.7,just = c("top","left"), default.units = "cm",draw=T)
plotText(label = "L", fontsize = 8, fontface = "bold",x = 9, y = 17.7,just = c("top","left"), default.units = "cm",draw=T)



dev.off()


### Part3S Incorporation of stranded peak intensity lead to higher performance of peak calling for stranded MeRIPseq data regardless of anti-m6A###############
#S3a The stranded RIP/Input enrichment of false positive peaks of best combo (Abcam_HEK293 and SYSY_HEK293)
#load the exomePeak best combo m6A peak in HEK293
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE)
dt.pruned.FPR.TPR.AUC.Abcam.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.all.combo.RDS") %>% mutate(Cell="Abcam_HEK293")
dt.pruned.FPR.TPR.AUC.SYSY.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.all.combo.RDS") %>% mutate(Cell="SYSY_HEK293")
dt.GLORI.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.AUC.Abcam.HEK %>% distinct(Cell,Method,ComboName,Benchmark,AUC.scaled,TPR20) ,
                             dt.pruned.FPR.TPR.AUC.SYSY.HEK %>% distinct(Cell,Method,ComboName,Benchmark,AUC.scaled,TPR20)) %>%
  left_join(x=.,y=dt.method.combo %>% dplyr::select(Method,ComboName,IsDefaultCombo))

dt.BestCombo.nonNEB <- dt.GLORI.AUC.nonNEB %>% filter(str_detect(Benchmark,"HEK293")) %>% group_by(Cell,Benchmark,Method) %>% dplyr::arrange(desc(AUC.scaled), desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()

dt.method.best.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.AUC.Abcam.HEK %>% inner_join(x=.,y=dt.BestCombo.nonNEB %>% dplyr::distinct(Cell,Method,ComboName),by=c("Cell","Method","ComboName")),
                                   dt.pruned.FPR.TPR.AUC.SYSY.HEK %>% inner_join(x=.,y=dt.BestCombo.nonNEB %>% dplyr::distinct(Cell,Method,ComboName),by=c("Cell","Method","ComboName")))
dt.method.best.AUC.nonNEB %>% distinct(Cell,Method,ComboName,Benchmark)
#highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
dt.method.best.AUC.nonNEB <- dt.method.best.AUC.nonNEB %>% mutate(Cell=factor(Cell,levels=c("Abcam_HEK293","SYSY_HEK293"))) %>%
  dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))
dt.method.best.AUC.nonNEB %>% distinct(Cell,Method,ComboName,Benchmark) %>% filter(Method=="exomePeak") %>% distinct(Method,ComboName)

#obtain the corresponding cutoff and method combo for exomePeak best combo
dt.method.best.AUC.nonNEB %>% filter(Method=="exomePeak") %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
# FPR  TPR cutoff    Method ComboName                          Benchmark  AUC AUC.scaled maxTPRmFPR TPR20   Cell                    label
# 1: 0.2 0.39   9.42 exomePeak   Combo26 dt.benchmark.GLORI_HEK293.m6As.bed 0.18       0.36       0.21  0.39 HEK293 exomePeak:ScaledAUC=0.36
dt.selected.method.combo.cutoff.FPR20 <- dt.method.best.AUC.nonNEB %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>% filter(!is.na(cutoff) & TPR>0) %>%
  dplyr::select(Cell,Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Best",Combo=paste0(Method,"_",ComboName))
#            Cell     Method ComboName cutoff FPR  TPR Group              Combo
# 1:  SYSY_HEK293 exomePeak2   Combo15   1.93 0.2 0.53  Best exomePeak2_Combo15
# 2:  SYSY_HEK293  exomePeak   Combo26   9.02 0.2 0.36  Best  exomePeak_Combo26
# 3:  SYSY_HEK293    MeTPeak   Combo18   6.33 0.2 0.37  Best    MeTPeak_Combo18
# 4:  SYSY_HEK293      TRESS   Combo27   1.44 0.2 0.17  Best      TRESS_Combo27
# 5:  SYSY_HEK293 MeRIPtools   Combo13   3.36 0.2 0.18  Best MeRIPtools_Combo13
# 6: Abcam_HEK293 exomePeak2   Combo15   1.59 0.2 0.38  Best exomePeak2_Combo15
# 7: Abcam_HEK293    MeTPeak    Combo2   4.17 0.2 0.21  Best     MeTPeak_Combo2
# 8: Abcam_HEK293  exomePeak   Combo18   3.79 0.2 0.20  Best  exomePeak_Combo18
# 9: Abcam_HEK293 MeRIPtools   Combo13   2.74 0.2 0.06  Best MeRIPtools_Combo13
#load selected method best combo FPR20 peak in Abcam_HEK and SYSY_HEK
#Abcam_HEK293
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB <- foreach(C = 1:nrow(dt.selected.method.combo.cutoff.FPR20),.combine = rbind)%do%{
  if(str_detect(dt.selected.method.combo.cutoff.FPR20$Cell[C],"Abcam")){
    LoadPeakMethod(Method = dt.selected.method.combo.cutoff.FPR20$Method[C], Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_Abcam/", dt.selected.method.combo.cutoff.FPR20$Combo[C])) %>%
      mutate(Group=paste0(dt.selected.method.combo.cutoff.FPR20$Method[C],"Best")) %>% mutate(PeakOverCutoff=score>= dt.selected.method.combo.cutoff.FPR20$cutoff[C])
  }else{
    if(str_detect(dt.selected.method.combo.cutoff.FPR20$Cell[C],"SYSY")){
      LoadPeakMethod(Method = dt.selected.method.combo.cutoff.FPR20$Method[C], Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_SYSY/", dt.selected.method.combo.cutoff.FPR20$Combo[C])) %>%
        mutate(Group=paste0(dt.selected.method.combo.cutoff.FPR20$Method[C],"Best")) %>% mutate(PeakOverCutoff=score>= dt.selected.method.combo.cutoff.FPR20$cutoff[C])
    }
  }
}
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>% group_by(Group) %>% slice_head(n=1)
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB[,.N,by=.(Group,Sample,PeakOverCutoff)] %>% dplyr::arrange(Sample,Group,PeakOverCutoff)

#divided into false positive peak(FP peak), false negative peak (FP peak), and other ture positive peak
#overlap with benchmark sites to define FP, TP, and FN
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
dt.Method.BestCombo.FPR20.HEK293.nonNEB.peaks.Benchmarkoverlap <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand) %>% distinct()
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#combine benchmark overlap into peak
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.HEK293.nonNEB.peaks.Benchmarkoverlap %>% filter(benchmark_name=="GLORI_HEK293") %>% dplyr::distinct(name,Sample,nBenchSite,nConfidentBenchSite), by=c("name","Sample")) %>%
  mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite))
#FP and TP, and FN
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap %>%
  mutate(PeakGroup=case_when(PeakOverCutoff==TRUE & nBenchSite>0 ~ "TP", PeakOverCutoff==FALSE & nBenchSite>0 ~ "FN", PeakOverCutoff==TRUE & nBenchSite==0 ~ "FP", .default = "TN"))
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap[,.N,by=.(Sample,Group,PeakGroup)] %>% dplyr::arrange(Sample,Group,PeakGroup)

#plot the count of FP and FN peaks of all method best combo at FPR20
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB  <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN")%>%
  group_by(Sample,Method,Group,PeakGroup) %>% mutate(nPeak=n_distinct(name)) %>% as.data.table() %>% distinct(Sample,Method,Group,PeakGroup,nPeak)
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB %>% mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam_HEK293",str_detect(Sample,"SYSY")~"SYSY_HEK293"))
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB %>% mutate(Method=factor(Method,levels=rev(c("MeRIPtools","TRESS","MeTPeak","exomePeak","exomePeak2"))),
                                                                                                                      PeakGroup=factor(PeakGroup,levels=c("TP","FP","FN")))
#add grey rect
# indices to highlight
pos <- seq_along(levels(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB$Method))
idx <- pos[pos %% 2!=0]
rects <- data.frame(xmin = pos[idx] - 0.5,xmax = pos[idx] + 0.5,ymin = -Inf,ymax = Inf)
p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB <- ggplot(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB, aes(x = Method, y = nPeak)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(aes(color=PeakGroup),fill="transparent",position = position_dodge2(padding = 0.05,width=0.6),width = 0.55,outlier.shape = NA,size = 0.4,alpha=0.7) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(aes(x = Method, y = nPeak, color=PeakGroup, shape=Cell), position = position_dodge2(padding = 0.05,width = 0.6),size=1.5, stroke = 0.6,alpha=0.7) +
  geom_rect(data = rects,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),inherit.aes = FALSE,fill = "grey92", alpha = 0.2)+
  # scales for shapes and fills
  scale_shape_manual(values = c(18,22),breaks=c("Abcam_HEK293","SYSY_HEK293"), name = "nonNEB sample") +
  scale_color_manual(values = c("#e9abac","#b5d2e8","#cee4b4"), breaks = c("TP","FP","FN"), name = "PeakType") +
  labs(x = NULL, y = "N of m6A peak from method's best combo") +
  guides(shape=guide_legend(nrow=2),color=guide_legend(nrow=2))+
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )
p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB

#check the positive and negative RIP/Input enrichment of FP peak on positive and negative strand
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB  <-dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB.Benchmarkoverlap %>% dplyr::filter(PeakGroup!="TN")
tmpdir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmpdir)){dir.create(tmpdir)}
#prepare method m6A
dt.Method.m6a.HEK <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB %>% dplyr::filter(str_detect(Sample,"HEK"))
#prepare the stranded.BAM.Depth
#calculate stranded bam depth
Samples <- unique(dt.Method.m6a.HEK$Sample)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/Figure3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID:mean_Ex_bin100_P1
selected.option <- "mean_Ex_bin100_P1"
dt.Method.m6a <- dt.Method.m6a.HEK %>% mutate(ComboName="Best") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
# org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
# org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite"
# org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
# org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed"
# org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/hg38_gencode.gene.bed"
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo="Best", .combine='rbind')%do%{
    foreach(s=Samples, .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      # selected.option <- dt.method.combo.selected[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",
                 intron.bed = "/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
  }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
# intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)


#load calculate bedtools-based peak intensity in both pos strand and neg strand
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB <- foreach(S = unique(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB$Sample), .combine='rbind')%do%{
  foreach(M = as.character(unique(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB[Sample==S,Method])), .combine = 'rbind')%do%{
    dt.intensity <- fread(paste0(stranded.bam.dir,"/",S,"__",M,"__Best_",selected.option,"_bedtools_intensity.tsv"),header=T) %>%
      left_join(x=.,y=dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB %>% dplyr::select(name,strand,Sample,Method,PeakGroup), by="name")
    dt.intensity
  }
}
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB  %>% group_by(Sample,Method,PeakGroup) %>% mutate(nPeak=n_distinct(name), nPeak.samestrand=n_distinct(name[strand.new==strand])) %>%
  as.data.table() %>% distinct(Sample,Method,PeakGroup,nPeak,nPeak.samestrand) %>% dplyr::arrange(Sample,Method,PeakGroup)
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>%
  mutate(score.samestrand=case_when(strand=="+" ~ pos, strand=="-" ~ neg), score.reversestrand = case_when(strand=="-" ~ pos, strand=="+" ~ neg))

dt.bothstrand.Method.BestCombo.FPR20.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>% dplyr::filter(PeakGroup=="FP")
p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB <-
  ggplot(dt.bothstrand.Method.BestCombo.FPR20.nonNEB, aes(x = score.samestrand, y = score.reversestrand)) +
  geom_hline(yintercept = 2,color = "grey60", linetype = "dashed")+
  geom_vline(xintercept = 2,color = "grey60", linetype = "dashed")+
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d'),dpi=300) +
  viridis::scale_color_viridis(name = "Density", option = "magma", trans = "log",  labels = function(x) scales::scientific(x,  digits = 2)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  # facet_wrap(~PeakGroup,ncol=3)+
  labs(x = "Intensity on same strand of peak", y = "Intensity on reverse strand of peak", title="FP peaks") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = c(0.8,0.5),legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
        plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1)),
        plot.subtitle = element_text(hjust = 0.5))

set.seed(888)  # Set the random seed for reproducibility
dt.bothstrand.Method.BestCombo.FPR20.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>% dplyr::filter(PeakGroup=="FN") %>%
  group_by(Sample,Method) %>% slice_sample(prop = 0.5) %>% as.data.table() %>%
  mutate(score.reversestrand=case_when(is.na(score.reversestrand) ~ 0, .default = score.reversestrand),score.samestrand=case_when(is.na(score.samestrand) ~ 0, .default = score.samestrand))
p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB <-
  ggplot(dt.bothstrand.Method.BestCombo.FPR20.nonNEB, aes(x = score.samestrand, y = score.reversestrand)) +
  geom_hline(yintercept = 2,color = "grey60", linetype = "dashed")+
  geom_vline(xintercept = 2,color = "grey60", linetype = "dashed")+
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.8,method='kde2d'),dpi = 300) +
  viridis::scale_color_viridis(name = "Density", option = "magma", trans = "log",  labels = function(x) scales::scientific(x,  digits = 2)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  # facet_wrap(~PeakGroup,ncol=3)+
  labs(x = "Intensity on same strand of peak", y = "Intensity on reverse strand of peak", title="FN peaks") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(legend.position = c(0.8,0.5),legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
        plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.15)), plot.subtitle = element_text(hjust = 0.5))

cowplot::plot_grid( p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB, p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB,nrow = 1)

# dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% group_by(PeakGroup,Method,Sample) %>%
#   mutate(nPeak=n_distinct(name), nPeak.samestrand.enrichment = n_distinct(name[score.samestrand>=2]), nPeak.reversestrand.enrichment = n_distinct(name[score.reversestrand>=2])) %>% as.data.table() %>%
#   dplyr::distinct(PeakGroup,Method,Sample,nPeak,nPeak.samestrand.enrichment,nPeak.reversestrand.enrichment) %>%
#   mutate(ratio.enrich.samestrand=nPeak.samestrand.enrichment/nPeak, ratio.enrich.reversestrand=nPeak.reversestrand.enrichment/nPeak) %>% dplyr::arrange(Sample,Method,PeakGroup)
# dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity <- dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity %>% group_by(PeakGroup,Method) %>%
#   mutate(avg.ratio.enrich.samestrand=mean(ratio.enrich.samestrand), avg.ratio.enrich.reversestrand=mean(ratio.enrich.reversestrand)) %>% as.data.table() %>%
#   arrange(Method,PeakGroup)
dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>% group_by(PeakGroup,Method,Sample) %>%
  mutate(EnrichGroup=case_when(score.reversestrand>=2 & (score.samestrand<2 | is.na(score.samestrand)) ~ "ReverseStrandOnly",
                               score.samestrand>=2 & (score.reversestrand<2 | is.na(score.reversestrand)) ~ "SameStrandOnly",
                               score.samestrand>=2 & score.reversestrand>=2 ~ "BothStrand", .default = "NeitherStrand")) %>%
  mutate(nPeak=n_distinct(name)) %>% mutate(Ratio.sames.only = n_distinct(name[EnrichGroup=="SameStrandOnly"])/nPeak, Ratio.reverse.only = n_distinct(name[EnrichGroup=="ReverseStrandOnly"])/nPeak,
                                            Ratio.both = n_distinct(name[EnrichGroup=="BothStrand"])/nPeak, Ratio.neither=n_distinct(name[EnrichGroup=="NeitherStrand"])/nPeak) %>% as.data.table() %>%
  dplyr::distinct(PeakGroup,Method,Sample,nPeak,Ratio.sames.only,Ratio.reverse.only,Ratio.both,Ratio.neither) %>% dplyr::arrange(Sample,Method,PeakGroup)


dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB  <- dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>%
  mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam_HEK293", str_detect(Sample,"SYSY")~"SYSY_HEK293")) %>% mutate(PeakGroup=factor(PeakGroup, levels=c("TP","FP","FN"))) %>%
  pivot_longer(cols = c("Ratio.sames.only", "Ratio.reverse.only",   "Ratio.both", "Ratio.neither"), names_to = "EnrichedStrandnessGroup",names_prefix = "Ratio.", values_to = "Ratio") %>% as.data.table() %>%
  mutate(Ratio=Ratio*100) %>% mutate(EnrichedStrandnessGroup=factor(EnrichedStrandnessGroup,levels=c("sames.only", "both", "reverse.only", "neither")))
#add grey rect
# indices to highlight
pos <- seq_along(levels(dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB$Method))
idx <- pos[pos %% 2!=0]
rects <- data.frame(xmin = pos[idx] - 0.5,xmax = pos[idx] + 0.5,ymin = -Inf,ymax = Inf)
p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB <-
  ggplot(data=dt.ratio.enrichment.bothstrand.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB, aes(x=Method, y=Ratio))+
  geom_boxplot(aes(color=EnrichedStrandnessGroup),fill="transparent",position = position_dodge2(padding = 0.01,width=0.6),width = 0.55,outlier.shape = NA,size = 0.4,alpha=0.7) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(aes(color=EnrichedStrandnessGroup, shape=Cell), position = position_dodge2(padding = 0.01,width = 0.6),size=1.5, stroke = 0.6,alpha=0.7) +
  geom_rect(data = rects,aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),inherit.aes = FALSE,fill = "grey92", alpha = 0.35)+
  # scales for shapes and fills
  scale_shape_manual(values = c(18,22),breaks=c("Abcam_HEK293","SYSY_HEK293"), name = "nonNEB sample") +
  # scale_color_manual(values = c("#d14e51","#5f93c4","#7da14e"), breaks = c("TP","FP","FN"), name = "PeakType") +
  scale_color_manual(values = c("#eaabac","#caafd4","#7bb0d5","#cfe4b6"), breaks = c("sames.only", "both", "reverse.only", "neither"), name = "EnrichedStrandnessGroup") +
  guides(shape=guide_legend(nrow = 2),color=guide_legend(nrow=2,title = "EnrichedType"))+
  labs(x = NULL, y = str_wrap("% of m6A peak display enrichment (intensity>2)",width=40)) +
  facet_wrap(~PeakGroup)+
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85)),
    axis.text.x = element_text(face = "plain", angle = 30,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )
p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB
cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB[PeakGroup=="TP",score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB[PeakGroup=="TP",score.samestrand],method = c("pearson","spearman")[1])

#whether benchmark sites m6A level correlate better with the new intensity?
dt.Method.BestCombo.FPR20.HEK293.peaks.BenchmarkScore.nonNEB <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB$Sample), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,site.score=V11,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".") %>% mutate(site.score=as.numeric(site.score))
  dt.overlap <- dt.overlap %>%  group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),max.BenchScore=max(site.score,na.rm=T), mean.BenchScore=mean(site.score,na.rm=T),median.BenchScore=median(site.score,na.rm=T)) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,max.BenchScore,mean.BenchScore,median.BenchScore)
  dt.overlap %>% mutate(Sample=c)
}

#incorporate the benchmark score to the peak
dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>% filter(str_detect(Sample,"HEK")) %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.HEK293.peaks.BenchmarkScore.nonNEB %>% dplyr::select(name,contains("BenchScore")),by="name")#calculate the cor (spearman, pearson) of score.old and score.samestrand with benchmarkscore(max,mean,median) in TP peaks

dt.cor.with.benchmarkscore.nonNEB <- foreach(P=c("TP","FN"), .combine='rbind')%do%{
  foreach(M = c("pearson","spearman"), .combine='rbind')%do%{
    cor.maxscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,max.BenchScore],method = M)$estimate
    cor.maxscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,max.BenchScore],method = M)$estimate
    cor.meanscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,mean.BenchScore],method = M)$estimate
    cor.meanscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,mean.BenchScore],method = M)$estimate
    cor.medianscore.old <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.old], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,median.BenchScore],method = M)$estimate
    cor.medianscore.new <- cor.test(dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,score.samestrand], dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB[PeakGroup==P,median.BenchScore],method = M)$estimate
    data.table(PeakGroup=P, CorMethod=M, cor.maxscore.old, cor.maxscore.new, cor.meanscore.old, cor.meanscore.new, cor.medianscore.old, cor.medianscore.new)
  }
}
dt.cor.with.benchmarkscore.nonNEB
# PeakGroup CorMethod cor.maxscore.old cor.maxscore.new cor.meanscore.old cor.meanscore.new cor.medianscore.old cor.medianscore.new
# 1:        TP   pearson        0.1405485        0.2611852        0.08081702      -0.014576266          0.07217005         -0.03663198
# 2:        TP  spearman        0.1719300        0.3087800        0.12489991      -0.004111702          0.10672274         -0.02619653
# 3:        FN   pearson        0.1813279        0.3297462        0.09739662       0.091006005          0.08122568          0.05837863
# 4:        FN  spearman        0.1567752        0.3415342        0.12638400       0.144531585          0.11208434          0.09969597

dt.cor <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.benchmarkscore.nonNEB  %>% dplyr::select(name,Sample,Method,PeakGroup,score.old,score.samestrand,max.BenchScore) %>%
  pivot_longer(cols = c("score.samestrand","score.old"),values_to = "score.tocompare", names_to = "MetricType", names_prefix = "score.") %>% as.data.table()
plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB <- vector(mode="list",length=4)
names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB) <- paste(rep(c("TP","FN"),each=2), rep(c("old","samestrand"),2), sep="_")
for(i in 1:length(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB)){
  P <- strsplit(names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB)[i],split="_") %>% sapply("[",1)
  M <- strsplit(names(plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB)[i],split="_") %>% sapply("[",2)
  dt.cor.input <- dt.cor %>% dplyr::filter(MetricType== M & PeakGroup==P)
  p <- ggplot(dt.cor.input, aes(x = score.tocompare, y = max.BenchScore)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(adjust = 1.0, alpha = 0.9, size = 0.1,method='kde2d'),dpi=300) +
    # geom_bin2d(inherit.aes = T,bins=200) +
    viridis::scale_color_viridis(name = "Point density", option = "magma", trans = "log") +
    # geom_density2d(aes(color=MetricType))+
    stat_cor(method = "spearman",label.x=quantile(dt.cor.input$score.tocompare,0.99,na.rm=T)*0.1, size=1.8)+
    # scale_color_manual(values =c("#f0766d","#839cd1"),breaks=c("samestrand","old"))+
    coord_cartesian(xlim = c(1,quantile(dt.cor.input$score.tocompare,0.99,na.rm=T)))+
    # facet_wrap(~MetricType,nrow=1)+
    labs(x = ifelse(M=="old","Original peak score", "New peak intensity"), y = "m6A level of bechmark sites in peak", title=paste0(P," peaks")) +
    theme_minimal(base_size = 6,base_family = "Helvetica") +
    theme(legend.position = "none", plot.title = element_text(face = "plain",hjust = 0.5,size = rel(1.)), plot.subtitle = element_text(hjust = 0.5))
  if(i %in% c(2,4) ){p <- p+theme(axis.title.y = element_blank())}
  if(i %in% c(1,2)){p <- p+theme(axis.title.x= element_blank())}
  plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB[[i]] <- p
}
p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB <- cowplot::plot_grid(plotlist = plist.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB, nrow = 2)

#the count of FP peaks (ratio of this type FP peaks with benchmark sites) due to enrichment in the reverse strand
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.nonNEB <- dt.Method.BestCombo.FPR20.m6A.TP.FP.FN.intensity.nonNEB %>%
  filter(PeakGroup=="FP" & score.reversestrand >=2 & score.samestrand<=2 & Method %in% c("MeTPeak","exomePeak","TRESS"))
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.nonNEB[,.N,by=.(Sample,Method)] %>% dplyr::arrange(Sample,Method)
#overlap with benchmark
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.HEK.nonNEB <- foreach(c = unique(dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.nonNEB$Sample) %>% grep(pattern="HEK",value=T), .combine='rbind')%do%{
  bed.peak <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.nonNEB %>% dplyr::filter(Sample==c) %>% mutate(seqnames=strsplit(name,split="_") %>% sapply("[",5),
                                                                                                              start=strsplit(name,split="_") %>% sapply("[",6),
                                                                                                              end=strsplit(name,split="_") %>% sapply("[",7),score=score.new,strand=strand.new) %>%
    dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}

dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.nonNEB <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.nonNEB %>% dplyr::filter(str_detect(Sample,"HEK")) %>%
  left_join(x=.,y=dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.HEK.nonNEB %>% dplyr::select(name,nBenchSite,nConfidentBenchSite),by="name")%>%
  mutate(nBenchSite=case_when(is.na(nBenchSite) ~ 0, .default = nBenchSite), nConfidentBenchSite=case_when(is.na(nConfidentBenchSite) ~ 0, .default = nConfidentBenchSite))

#calculate the ratio of those FP reverse
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.nonNEB %>% group_by(Sample,Method) %>% mutate(nPeak=n_distinct(name)) %>%
  mutate(ratio.Benchmark=n_distinct(name[nBenchSite>0])/nPeak, ratio.ConfidentBenchmark = n_distinct(name[nConfidentBenchSite>0])/nPeak) %>%
  as.data.table() %>% distinct(Sample,Method,nPeak,ratio.Benchmark,ratio.ConfidentBenchmark) %>% dplyr::arrange(Sample,Method)
dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB <- dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB %>% group_by(Method) %>%
  mutate(avg.ratio.Benchmark=mean(ratio.Benchmark), avg.ratio.ConfidentBenchmark=mean(ratio.ConfidentBenchmark)) %>% as.data.table() %>% distinct(Method,avg.ratio.Benchmark,avg.ratio.ConfidentBenchmark) %>%
  pivot_longer(cols = c("avg.ratio.Benchmark", "avg.ratio.ConfidentBenchmark"),names_to = "ContainSite",values_to = "Ratio",names_prefix = "avg.ratio.") %>% as.data.table() %>%
  mutate(ContainSite=factor(ContainSite,levels=c("Benchmark","ConfidentBenchmark"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>% mutate(label.pt=paste0(round(Ratio*100,1),"%")) %>%
  mutate(Method=factor(Method,levels=unique(Method)))
p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB <-
  ggplot(data=dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB,
         aes(x=Method,y=Ratio, fill=ContainSite))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB %>% dplyr::filter(ContainSite=="Benchmark m6As"),
            aes(label = label.pt, y =  Ratio+0.02, x=Method),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB %>% dplyr::filter(ContainSite=="Confident benchmark m6As"),
            aes(label = label.pt, y =  Ratio+0.02, x=Method),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y=str_wrap("% of FP peaks with intensity(>2) at reverse strand",width=40))+
  guides(fill=guide_legend(title = expression("Peaks\ncontain"),nrow = 2))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  # scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB

#S3c The distribution of AUC of six method's  all combo with all alternative options (nonNEB_HEK293 and nonNEB_mESC)
#load all combo option AUC
dt.pruned.FPR.TPR.option.AUC.Abcam.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.RDS") %>% mutate(Cell="Abcam_HEK293")
dt.pruned.FPR.TPR.option.AUC.Abcam.HEK %>% distinct(Cell,Method,ComboName,optionID,AUC.scaled)
dt.pruned.FPR.TPR.option.AUC.SYSY.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.RDS") %>% mutate(Cell="SYSY_HEK293")
dt.pruned.FPR.TPR.option.AUC.SYSY.HEK %>% distinct(Cell,Method,ComboName,optionID,AUC.scaled)
#load all combo AUC
dt.pruned.FPR.TPR.combo.AUC.Abcam.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.all.combo.RDS") %>% mutate(Cell="Abcam_HEK293")
dt.pruned.FPR.TPR.combo.AUC.SYSY.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.all.combo.RDS") %>% mutate(Cell="SYSY_HEK293")
dt.pruned.FPR.TPR.combo.AUC.nonNEB.HEK <- rbind(dt.pruned.FPR.TPR.combo.AUC.Abcam.HEK,dt.pruned.FPR.TPR.combo.AUC.SYSY.HEK )
dt.BestCombo.nonNEB <- dt.pruned.FPR.TPR.combo.AUC.nonNEB.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Cell,Method,ComboName,AUC.scaled,TPR20) %>%
  group_by(Cell,Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.DefaultCombo.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.combo.AUC.Abcam.HEK %>% distinct(Cell,Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")),
                                    dt.pruned.FPR.TPR.combo.AUC.SYSY.HEK %>% distinct(Cell,Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")))
dt.BestCombo.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.combo.AUC.Abcam.HEK %>% distinct(Cell,Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.nonNEB %>% dplyr::distinct(Cell,Method,ComboName),by=c("Cell","Method","ComboName")),
                                 dt.pruned.FPR.TPR.combo.AUC.SYSY.HEK %>% distinct(Cell,Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.nonNEB %>% dplyr::distinct(Cell,Method,ComboName),by=c("Cell","Method","ComboName")))

dt.AUC.nonNEB.all.option <- rbind(dt.pruned.FPR.TPR.option.AUC.Abcam.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Cell,Method,ComboName,optionID,AUC.scaled,TPR20),
                                  dt.pruned.FPR.TPR.option.AUC.SYSY.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Cell,Method,ComboName,optionID,AUC.scaled,TPR20))

dt.BestComboOption.nonNEB <- dt.AUC.nonNEB.all.option  %>% group_by(Cell,Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20))

dt.AUC.nonNEB.all.option <- dt.AUC.nonNEB.all.option %>% left_join(x=.,y=dt.BestComboOption.nonNEB %>% mutate(IsBestComboOption=TRUE) %>% dplyr::select(Cell,Method,ComboName,optionID,IsBestComboOption),by=c("Cell","Method","ComboName","optionID")) %>%
  mutate(ComboGroup=case_when(IsBestComboOption==TRUE ~ "BestComboOption", .default = "OthersComboOption")) %>% dplyr::select(-IsBestComboOption)
dt.AUC.nonNEB.all.option <- rbind(dt.AUC.nonNEB.all.option %>% mutate(Group=ComboGroup) %>% dplyr::select(Cell,Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group),
                                  dt.DefaultCombo.AUC.nonNEB %>% mutate(optionID=NA, Group="DefaultCombo") %>% dplyr::select(Cell,Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group),
                                  dt.BestCombo.AUC.nonNEB  %>% mutate(optionID=NA, Group="BestCombo") %>% dplyr::select(Cell,Method,ComboName,optionID,AUC.scaled,TPR20,Cell,Group))
dt.AUC.nonNEB.all.option <- dt.AUC.nonNEB.all.option %>% mutate(Method=factor(Method,levels=unique(dt.BestComboOption.nonNEB$Method)),Group=factor(Group,levels=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo")))
p.AUC.nonNEB.all.option <- ggplot(dt.AUC.nonNEB.all.option, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.AUC.nonNEB.all.option %>% dplyr::filter(Group =="OthersComboOption"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=0.1,alpha=0.15, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.AUC.nonNEB.all.option %>% dplyr::filter(Group =="DefaultCombo"), aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=1.5)+
  geom_jitter(data = dt.AUC.nonNEB.all.option %>% dplyr::filter(Group =="BestCombo") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=1.5)+
  geom_point(data = dt.AUC.nonNEB.all.option %>% dplyr::filter(Group =="BestComboOption") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#eb4601",size=3)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled)", title = paste0("9120 (285 Combo x 32 Option) ComboOption combination")) +
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9)),
    axis.title = element_text(face = "plain", size = rel(0.95)),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 30,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    legend.position = c(0.92, 0.85),  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )
p.AUC.nonNEB.all.option
# legend-only plot for color: hide linetype legend
p_col_legend <- ggplot(data=dt.AUC.nonNEB.all.option,aes(x = Method, y = AUC.scaled)) +
  geom_jitter(aes(shape=Group,color=Group,size=Group), width = 0.18, height = 0,  stroke = 0.6,alpha=1) +
  scale_shape_manual(values = c(18,22,18,18),breaks=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  scale_size_manual(values = c(3,0.5,2,2), breaks = c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  scale_color_manual(values = c("#eb4601","grey70","#e9abac","#7bb0d5"), breaks = c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo"), name = "Group") +
  theme_minimal() +
  guides(linetype = "none") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "right",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"), legend.title = element_blank(),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
# extract legend grobs
leg_col <- cowplot::get_legend(p_col_legend)
# assemble: main plot + legends at arbitrary positions
p.AUC.nonNEB.all.option  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.AUC.nonNEB.all.option , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1.15, y = 0.8, width = 0.22, height = 0.25, hjust = 1.2, vjust = 0.5)
p.AUC.nonNEB.all.option

#S3d The AUC of six method's best combo and option (nonNEB_HEK293 and nonNEB_mESC)
dt.pruned.TPR.FPR.nonNEB.Method.BestOption <- rbind(dt.pruned.FPR.TPR.option.AUC.Abcam.HEK %>% inner_join(x=.,y=dt.BestComboOption.nonNEB %>% dplyr::select(Cell,Method,ComboName,optionID),by=c("Cell","Method","ComboName","optionID")) %>%
                                                      filter(str_detect(Benchmark,"GLORI")),
                                                    dt.pruned.FPR.TPR.option.AUC.SYSY.HEK %>% inner_join(x=.,y=dt.BestComboOption.nonNEB %>% dplyr::select(Cell,Method,ComboName,optionID),by=c("Cell","Method","ComboName","optionID")) %>%
                                                      filter(str_detect(Benchmark,"GLORI")))

dt.pruned.TPR.FPR.nonNEB.Method.BestOption.label <- dt.pruned.TPR.FPR.nonNEB.Method.BestOption %>% distinct(Cell,Method,ComboName,optionID,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>%
  mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.pruned.TPR.FPR.nonNEB.Method.BestOption$TPR),2)) %>% mutate(xpos=0.6)
p.pruned.TPR.FPR.nonNEB.Method.BestOption  <-  ggplot(data=dt.pruned.TPR.FPR.nonNEB.Method.BestOption ,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of m6A peak from best combo and option combination",y=str_wrap("Surrogate TPR of m6A peak from best combo and option combination",width = 40))+
  geom_text(data=dt.pruned.TPR.FPR.nonNEB.Method.BestOption.label, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  facet_wrap(~Cell)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.minor.ticks.x.bottom = element_blank(),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.pruned.TPR.FPR.nonNEB.Method.BestOption

########################## S3e the novel m6A+ gene in TRESSBestOption vs exomePeak2Best at (FPR20 cutoff) in Abcam_HEK and SYSY_HEK

#calling  m6A peaks from TRESS BestOption at FPR20 cutoff
dt.pruned.TPR.FPR.nonNEB.Method.BestOption %>% filter(Method=="TRESS" & FPR==0.2)
# FPR  TPR cutoff Method ComboName            optionID                          Benchmark  AUC AUC.scaled maxTPRmFPR TPR20         Cell
# 1: 0.2 0.56   1.78  TRESS   Combo40 mean_Tx_bin100_P0.5 dt.benchmark.GLORI_HEK293.m6As.bed 0.27       0.54       0.38  0.56 Abcam_HEK293
# 2: 0.2 0.63   1.89  TRESS   Combo32 mean_Tx_bin100_P0.5 dt.benchmark.GLORI_HEK293.m6As.bed 0.29       0.58       0.44  0.63  SYSY_HEK293
dt.TPR20.TRESS.BestOption.cutoff <- dt.pruned.TPR.FPR.nonNEB.Method.BestOption %>% filter(Method=="TRESS" & FPR==0.2 & str_detect(Benchmark,"GLORI_HEK293")) %>% dplyr::select(Cell,Method,ComboName,optionID,cutoff)
# Cell Method ComboName            optionID cutoff
# 1: Abcam_HEK293  TRESS   Combo40 mean_Tx_bin100_P0.5   1.78
# 2:  SYSY_HEK293  TRESS   Combo32 mean_Tx_bin100_P0.5   1.89
#obtain  TRESS    Combo1 mean_Tx_bin50_P2 for HEK293
#run bedtools intensity
#HEK293
#prepare the stranded.BAM.Depth
Samples <- unique(dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB$Sample) %>% grep(pattern="HEK",value=T)
n.cores <- 80
bin.dir <- "~/anaconda3/envs/MyPackage/bin/"
stranded.bam.dir = "/data/m6A_calling_strategy/Analysis/FigureS3_stranded_bam/"
if(!dir.exists(stranded.bam.dir)){dir.create(stranded.bam.dir)}
registerDoParallel(cl=min(floor(n.cores/4), length(Samples)))
dt.stranded.BAM.Depth <- foreach(s=Samples,.combine='rbind')%dopar%{
  sample.input.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_Input",value=T)
  sample.ip.s.bam <- list.files("/data/m6A_calling_strategy/6_alignment_star/",pattern="(.R2.bam)$",full.names = T) %>% grep(pattern=s,value=T) %>% grep(pattern="_RIP",value=T)
  system(command = paste0("cp ",sample.input.s.bam, " ", stranded.bam.dir))
  system(command = paste0("cp ",sample.ip.s.bam, " ", stranded.bam.dir))
  Depth <- c(as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.input.s.bam), wait = T, intern = T)),
             as.integer(system(command = paste0(bin.dir,"samtools view --threads 4 -F 256 -c ", sample.ip.s.bam), wait = T, intern = T)))
  data.table(Sample=paste0(s,c("_Input","_RIP")), Depth=Depth) %>% mutate(Replicate=s)
}
stopImplicitCluster()
dt.stranded.BAM.Depth  <- dt.stranded.BAM.Depth %>% mutate(Mapped=Depth/1000000)
dt.stranded.BAM.Depth <- dt.stranded.BAM.Depth %>% mutate(Library=ifelse(str_detect(Sample,"_Input"),"Input","RIP")) %>%
  tidyr::pivot_wider(id_cols = c("Replicate"),names_from = c("Library"),values_from = c("Mapped")) %>% as.data.table() %>%
  mutate(RIPLibraryScaleFactor=paste0(round(RIP,2),"/",round(Input,2)))
print(dt.stranded.BAM.Depth)
tmp.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
if(!dir.exists(tmp.dir)){dir.create(tmp.dir)}
script.dir <- "~/software/mysterypackage/mysterypackage/inst/scripts/"
#prepare the optionID
selected.option <- dt.TPR20.TRESS.BestOption.cutoff$optionID
#load the old method m6A
dt.Method.m6a <- foreach(i = 1:nrow(dt.TPR20.TRESS.BestOption.cutoff),.combine = 'rbind')%do%{
  if(str_detect(dt.TPR20.TRESS.BestOption.cutoff$Cell[i],"Abcam")){
    LoadPeakMethod(Method=dt.TPR20.TRESS.BestOption.cutoff$Method[i],
                   Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_Abcam/",dt.TPR20.TRESS.BestOption.cutoff$Method[i],"_",dt.TPR20.TRESS.BestOption.cutoff$ComboName[i])) %>%
      mutate(Method=dt.TPR20.TRESS.BestOption.cutoff$Method[i], ComboName=dt.TPR20.TRESS.BestOption.cutoff$ComboName[i])
  }else{
    if(str_detect(dt.TPR20.TRESS.BestOption.cutoff$Cell[i],"SYSY")){
      LoadPeakMethod(Method=dt.TPR20.TRESS.BestOption.cutoff$Method[i],
                     Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_SYSY/",dt.TPR20.TRESS.BestOption.cutoff$Method[i],"_",dt.TPR20.TRESS.BestOption.cutoff$ComboName[i])) %>%
        mutate(Method=dt.TPR20.TRESS.BestOption.cutoff$Method[i], ComboName=dt.TPR20.TRESS.BestOption.cutoff$ComboName[i])
    }
  }
}

dt.Method.m6a <- dt.Method.m6a %>% filter(str_detect(seqnames,"chr"))  %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
org.gtf_file="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf"
org.sqlite_file="/data/m6A_calling_strategy/RIPPeakS_resource/UCSC.hg38.knownGene.sqlite"
org.genome_file="/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt"
org.intron.bed="/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed"
org.genebed="/data/m6A_calling_strategy/RIPPeakS_resource/hg38_gencode.gene.bed"
log.dir <- "/data/m6A_calling_strategy/Analysis/tmp"
dt.parameter.bedtools.parallel <- foreach(M=unique(dt.Method.m6a$Method), .combine='rbind')%do%{
  foreach(combo=dt.Method.m6a[Method==M,unique(ComboName)], .combine='rbind')%do%{
    foreach(s=dt.Method.m6a[Method==M & ComboName==combo,unique(Sample)], .combine = 'rbind')%do%{
      fwrite(dt.Method.m6a %>% dplyr::filter(Method==M & ComboName==combo & Sample==s), file=paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),row.names=F,col.names=T,sep="\t")
      selected.option <- dt.TPR20.TRESS.BestOption.cutoff[Method==M & ComboName==combo,optionID]
      data.table(input.m6a = paste0(tmp.dir, "/", M,"_",combo, "_bedtools_intensity_input_m6A_",s,".tsv"),
                 prefix=s,
                 RIPLibraryScaleFactor=dt.stranded.BAM.Depth[match(s,Replicate),RIPLibraryScaleFactor],
                 stranded.bam.dir = stranded.bam.dir,
                 bedtools_path = bin.dir,
                 dir.tmp = paste0(tmp.dir,"/", M, "_", combo, "_bedtools_intensity_tmp_",s),
                 genome_file = "/data/m6A_calling_strategy/RIPPeakS_resource/chrNameLength_hg38.txt",
                 intron.bed = "/data/m6A_calling_strategy/RIPPeakS_resource/dt.intron.hg38.bed",
                 intronremoved = str_detect(selected.option, "Ex"),
                 bedtools_mode = ifelse(str_detect(selected.option, "mean"),"mean","counts"),
                 binsize = as.integer(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",3) %>% str_replace(pattern="bin",replacement = "")),
                 Pseudocount = as.numeric(selected.option %>% strsplit(split="_",fixed = T) %>% sapply("[",4) %>% str_replace(pattern="P",replacement = "")),
                 SummaryMethod = "max",
                 ComboName = paste0(M, "__", combo, "_",selected.option),
                 Strandness = "s",
                 Logfile=paste0(log.dir,"/",M, "_", combo,"_bedtools_intensity_",s,".log")
      )
    }
  }
}
fwrite(dt.parameter.bedtools.parallel, file=paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt"),sep="\t",row.names = F,col.names=F)
intensity.nthread <- max(floor(n.cores/(quantile(c(dt.stranded.BAM.Depth$Input,dt.stranded.BAM.Depth$RIP),0.75)/30)/8),1)
# intensity.nthread  <- 10
parallel.Method.bedtools.cmd <- paste0(bin.dir,"/parallel --workdir ",  tmp.dir, " -j ", min(intensity.nthread,10)," --will-cite -a ", paste0(tmp.dir,"/parameter.Method.bedtools.parallel.txt") ," --colsep '\t' '", paste0(bin.dir,"/Rscript "),
                                       paste0(script.dir,"/Rscript_bedtools_intensity.R")," {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} 1>{16}  2>&1 '")
system(command = parallel.Method.bedtools.cmd, wait = T)

#obtain m6a at FPR cutoff
dt.method.combo.selected <- dt.TPR20.TRESS.BestOption.cutoff
dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
  foreach(combo = dt.method.combo.selected[Method==M,unique(ComboName)], .combine='rbind')%do%{
    # message(paste0("Obtain peak at FPR cutoff at ",ifelse(is.null(Model.FPR.cutoff), "default", Model.FPR.cutoff)," for ", M, " ", combo, " ", dt.method.combo.selected[Method==M & ComboName==combo,optionID]))
    intensity.files <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = T) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T)
    names(intensity.files) <- list.files(stranded.bam.dir,pattern = paste0(M,"__",combo),full.names = F) %>% grep(pattern="(_bedtools_intensity.tsv)$",value=T) %>%
      strsplit(split="__",fixed=T) %>% sapply("[",1)
    dt.intensity <- foreach(s=names(intensity.files),.combine='rbind')%do%{data.table::fread(file = intensity.files[s]) %>% mutate(Sample=s)}#load intensity for all sample
    dt.intensity <- dt.intensity %>% dplyr::select(name,Sample,pos,neg)
    #pull out m6a with intensity higher than cutoff at FPR20
    #filtered and un overlapped  m6a at pos strand
    dt.m6a.filtered.pos <- dt.intensity %>% dplyr::filter(pos >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.pos)>0){
      dt.m6a.filtered.pos <- dt.m6a.filtered.pos %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="+") %>%
        dplyr::select(seqnames,start,end,name,score=pos,strand,Sample)
      dt.m6a.filtered.pos.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.pos,keep.extra.columns = T)
      dt.m6a.filtered.pos.merged <- foreach(sam = unique(dt.m6a.filtered.pos$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.pos.gr[dt.m6a.filtered.pos.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.pos.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    #filtered  and un overlapped m6a at neg strand
    dt.m6a.filtered.neg <- dt.intensity %>% dplyr::filter(neg >= dt.method.combo.selected[Method==M & ComboName==combo,cutoff])
    if(nrow(dt.m6a.filtered.neg)>0){
      dt.m6a.filtered.neg <- dt.m6a.filtered.neg %>% dplyr::mutate(name = name %>% strsplit(split=paste0(M,"_"),fixed=T) %>% sapply(tail,1)) %>%
        dplyr::mutate(seqnames=name %>% strsplit(split="_",fixed=T) %>% sapply("[",1),
                      start=name %>% strsplit(split="_",fixed=T) %>% sapply("[",2),
                      end=name %>% strsplit(split="_",fixed=T) %>% sapply("[",3),
                      strand="-") %>%
        dplyr::select(seqnames,start,end,name,score=neg,strand,Sample)
      dt.m6a.filtered.neg.gr <- GenomicRanges::makeGRangesFromDataFrame(dt.m6a.filtered.neg,keep.extra.columns = T)
      dt.m6a.filtered.neg.merged <- foreach(sam = unique(dt.m6a.filtered.neg$Sample), .combine='rbind')%do%{
        GenomicRanges::reduce(dt.m6a.filtered.neg.gr[dt.m6a.filtered.neg.gr$Sample==sam,]) %>% as.data.table() %>%
          mutate(score=1000,name=paste(seqnames,start,end,strand,sep="_"),Sample=sam) %>%
          dplyr::select(seqnames,start,end,name,score,strand,Sample)
      }
    }else{
      dt.m6a.filtered.neg.merged <- data.table(seqnames=NULL,start=NULL,end=NULL,name=NULL,score=NULL,strand=NULL,Sample=NULL)
    }
    dt.m6a.filtered <- rbind(dt.m6a.filtered.pos.merged, dt.m6a.filtered.neg.merged)
    if(nrow(dt.m6a.filtered)>0){
      #calculate exonic sum width to filter extreme long peak (<=6kb)
      dt.m6a.filtered.sumwidth <- LongPeakRemoveIntron(dt.peak = dt.m6a.filtered %>% distinct(seqnames,start,end,name,score,strand),#unique m6a peak
                                                       dt.intron = data.table::fread(org.intron.bed),
                                                       genome_file = org.genome_file)
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.sumwidth  %>% dplyr::distinct(name,SumWidth), by=c("name")) %>%
        dplyr::mutate(SumWidth=case_when(is.na(SumWidth) ~ 0, .default = SumWidth))
      #the overlap gene
      dt.m6a.filtered.overlapped.gene <- bedtoolsr::bt.intersect(a=dt.m6a.filtered, b=data.table::fread(org.genebed),
                                                                 s=T, f=0.5, F=0.8, c=T,e = T) %>% as.data.table()
      dt.m6a.filtered <- dt.m6a.filtered %>% left_join(x=.,y=dt.m6a.filtered.overlapped.gene %>% dplyr::select(name=V4,Sample=V7,OverlappedNoGene=V9) %>%
                                                         dplyr::distinct(Sample,name,OverlappedNoGene), by=c("Sample","name"))
      dt.m6a.filtered %>% mutate(Method=M, ComboName=combo, optionID=dt.method.combo.selected[Method==M & ComboName==combo,optionID])
    }
  }
}
dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB  <- dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB %>% mutate(name=paste(Method,ComboName,optionID,Sample,name,sep="_"))

##pull out the m6A peaks from exomePeak2 BestCombo at FPR20 cutoff
dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>% filter(PeakOverCutoff==TRUE & Method=="exomePeak2")
dt.exomePeak2.BestCombo.HEK.FPR20.peaks.nonNEB <- dt.Method.BestCombo.FPR20.m6A.HEK293.nonNEB %>% filter(PeakOverCutoff==TRUE & Method=="exomePeak2")

#check the overlap between TRESS best option peaks with exomePeak2 best combo peaks
dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB <- foreach(c = unique(dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB$Sample), .combine='rbind')%do%{
  bed.peak <- dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.exomePeak2 <- dt.exomePeak2.BestCombo.HEK.FPR20.peaks.nonNEB %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.exomePeak2, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
  dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(exomePeak2overlap=ifelse(sum.overlap/width>=0.0001,"exomePeak2+","exomePeak2-"))
  bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,exomePeak2overlap), by="name") %>% mutate(Sample=c)
}
dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB[,.N,by=.(Sample,exomePeak2overlap)] %>% dplyr::arrange(Sample,exomePeak2overlap)
#           Sample exomePeak2overlap     N
# 1: HEK_nonNEB_mRNA1       exomePeak2+ 10526
# 2: HEK_nonNEB_mRNA1       exomePeak2-  2733
# 3: HEK_nonNEB_mRNA2       exomePeak2+ 14313
# 4: HEK_nonNEB_mRNA2       exomePeak2-  2746

#check the exomePeak2- and exomePeak2+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
#HEK293
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
dt.TRESS.BestOption.HEK293.peaks.Benchmarkoverlap.nonNEB <- foreach(c = unique(dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB$Sample), .combine='rbind')%do%{
  bed.peak <- dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  if(str_detect(c,"HEK")){
    bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  }
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  if(str_detect(c,"HEK")){dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.GLORI.HEK293.confident[PerturbEffect.FTO< -0.1 & PerturbEffect.METTL3i < -0.1,name])}else{
    dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  }
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
#obtain the TRESS Best option only m6A+ gene
dt.TRESS.BestOption.HEK293.peaks.annotgene.nonNEB <- annot_peak(peak.bed=dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                                strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                annot_type = "gene")
dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB <- dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB %>% left_join(x=.,y=dt.TRESS.BestOption.HEK293.peaks.annotgene.nonNEB %>% dplyr::select(name,OverlappedGenes),by="name")
dt.TRESS.BestOption.specific.m6Agene.HEK293.nonNEB <-   dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestOption.specific=all(exomePeak2overlap=="exomePeak2-")) %>% as.data.table() %>%#if all peak of one gene is exomePeak2-, then this gene is BestOption specific m6A+ gene
  dplyr::filter(BestOption.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
#obtain the BestOption specific peaks overlapped with (high confident) benchmark sites
dt.TRESS.BestOption.specific.peak.benchmark.HEK293.nonNEB <-  dt.TRESS.BestOption.specific.m6Agene.HEK293.nonNEB %>% filter(exomePeak2overlap=="exomePeak2-" & BestOption.specific==TRUE) %>%
  left_join(x=.,dt.TRESS.BestOption.HEK293.peaks.Benchmarkoverlap.nonNEB,by=c("name","Sample")) %>%
  dplyr::filter(nBenchSite>0)
#rank the BestCombo specific m6A+ genes with benchmark m6A sites support
dt.TRESS.BestOption.specific.m6Agene.benchmark.HEK293.nonNEB <- dt.TRESS.BestOption.specific.peak.benchmark.HEK293.nonNEB %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
  group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
  as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))

#visualization of ratio of exomePeak2 ratio of TRESS BestOption peaks in Abcam_HEK293 and SYSY_HEK293
dt.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB <- dt.TRESS.BestOption.m6A.FPR20.HEK293.exomePeak2Bestoverlap.nonNEB %>% group_by(Sample) %>% mutate(Ratio.Peak.BestOptionOnly = n_distinct(name[exomePeak2overlap=="exomePeak2-"])/n_distinct(name),
                                                                                                                                                        nPeak.BestOptionOnly= n_distinct(name[exomePeak2overlap=="exomePeak2-"]), nPeak.BestOption=n_distinct(name)) %>%
  as.data.table() %>% distinct(Sample,Ratio.Peak.BestOptionOnly,nPeak.BestOptionOnly,nPeak.BestOption)

dt.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB <- dt.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB %>% mutate(label1=paste0(nPeak.BestOptionOnly,"/",nPeak.BestOption), label2=paste0(round(Ratio.Peak.BestOptionOnly*100,1),"%")) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam_HEK293",str_detect(Sample,"SYSY")~"SYSY_HEK293"))
p.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB <-
  ggplot(dt.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB, aes(y = Sample, x = Ratio.Peak.BestOptionOnly, fill = Cell)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  # percent label to the right of the bar
  geom_text(aes(label = label2, x =  Ratio.Peak.BestOptionOnly + 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  #N peak label to the left of the bar
  geom_text(aes(label = label1, x = Ratio.Peak.BestOptionOnly/2),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  scale_fill_manual(values = c("#7b93c6","#acd9ee"),breaks = c("SYSY_HEK293","Abcam_HEK293")) +
  labs(x = str_wrap("% of novel m6A peaks in TRESS's best option v.s exomePeak2's best combo",width = 40), y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB

#visualization of count of BestCombo only m6A+ gene with benchmark m6As and highly confident m6As
dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB <- dt.TRESS.BestOption.specific.peak.benchmark.HEK293.nonNEB %>% group_by(Sample) %>%
  mutate(No.BestOptionOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
         No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
  distinct(Sample,No.BestOptionOnly.m6AGene.BenchmarkSite,No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite)

dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB <- dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB %>%
  pivot_longer(cols = c("No.BestOptionOnly.m6AGene.BenchmarkSite", "No.BestOptionOnly.m6AGene.ConfidentBenchmarkSite"),values_to = "No.BestOptionOnly.m6AGene", names_to = "GeneGroup",names_prefix = "No.BestOptionOnly.m6AGene.") %>%
  as.data.table() %>% mutate(GeneGroup=factor(GeneGroup, levels=c("BenchmarkSite","ConfidentBenchmarkSite"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam_HEK293",str_detect(Sample,"SYSY")~"SYSY_HEK293"))

dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB <- dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB %>% mutate(label.count=paste0(paste0("No.Gene=",No.BestOptionOnly.m6AGene)))
p.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB  <- ggplot(data=dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB, aes(x=Sample,y=No.BestOptionOnly.m6AGene, fill=GeneGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB %>% dplyr::filter(GeneGroup=="Benchmark m6As"),
            aes(label = label.count, y =  No.BestOptionOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 1.8, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB %>% dplyr::filter(GeneGroup=="Confident benchmark m6As"),
            aes(label = label.count, y =  No.BestOptionOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 1.8, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y=str_wrap("N of novel m6A+ genes in TRESSS v.s exomePeak2 best combo",width=40))+
  guides(fill=guide_legend(title = "Peaks contain"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB

#visualization of IGV plot of TRESS Best Option specific m6A+ gene in HEK293 and mESC
dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
# dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

#Abcam
selected.Novel.m6Agene.Abcam.HEK <-  dt.TRESS.BestOption.specific.m6Agene.benchmark.HEK293.nonNEB %>% filter(str_detect(Sample,"Abcam"))  %>% dplyr::filter(nRep==2) %>%
  distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=10)
dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==unique(selected.Novel.m6Agene.Abcam.HEK$OverlappedGenes)[8])
StudySamples <- c("HEK_Abcam_mRNA1", "HEK_Abcam_mRNA2")
## Create a page (7.5*7.5cm)
window.size=20#100kb
pseudocount=1
track.height=0.5
pdf(paste0("Fig3S_Abcam_HEK293_", dt.selected.gene$gene_name,"_TRESS_BestOption_specific_peak.pdf"))
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
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

# #add peak in this region
dt.peak.to.see <- rbind(dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="TRESS_BestOption"),
                        dt.exomePeak2.BestCombo.HEK.FPR20.peaks.nonNEB %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="exomePeak2_BestCombo")) %>%
  mutate(Group=factor(Group, levels=c("TRESS_BestOption","exomePeak2_BestCombo"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)

for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  exomePeak2.bestcombo <- plotRanges(data = dt.peak.region, params = region, order = "random",
                                     fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),strandSplit = F,
                                     x = 0.5, y = 1.5+(s-1+length(StudySamples))*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm"
  )
  sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5+length(StudySamples))*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
  sample_rect <- plotRect(x=0.5,y=1.5+(s-1+length(StudySamples))*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#track for benchmark m6A sites (perturbation effect)
#FTO GLORI HEK293
signal_confident.benchmark.FTO <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.FTO,name,strand) %>% filter(strand==dt.selected.gene$strand & !is.na(score)), params = region,
                                             x = 0.5, y = 1.5+2*length(StudySamples)*track.height, width = 5, height = track.height,
                                             just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="FTO treatment", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+2*length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
#METTL3i GLORI HEK293
signal_confident.benchmark.METTL3i <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.METTL3i), params = region,
                                                 x = 0.5, y = 1.5+(2*length(StudySamples)+1)*track.height, width = 5, height = track.height,
                                                 just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+1.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="METTL3 inhibition", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+1.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fontface="plain")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in GLORI_HEK293"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fontface="plain")
legendPlot.exomePeak2.peak <- plotLegend(legend = c("TRESS_BestOption","exomePeak2_BestCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 1.5, height = 1,
                                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                         title = expression("PeakGroup"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


#SYSY
selected.Novel.m6Agene.SYSY.HEK <-  dt.TRESS.BestOption.specific.m6Agene.benchmark.HEK293.nonNEB %>% filter(str_detect(Sample,"SYSY"))  %>% dplyr::filter(nRep==2) %>%
  distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=20)
dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==unique(selected.Novel.m6Agene.SYSY.HEK$OverlappedGenes)[19])#ACAD8
StudySamples <- c("HEK_SYSY_mRNA1", "HEK_SYSY_mRNA2")
## Create a page (7.5*7.5cm)
window.size=30#100kb
pseudocount=1
track.height=0.5
# pdf(paste0("Fig3S_SYSY_HEK293_", dt.selected.gene$gene_name,"_TRESS_BestOption_specific_peak.pdf"))
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
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

# #add peak in this region
dt.peak.to.see <- rbind(dt.TRESS.BestOption.m6A.FPR20.HEK293.nonNEB %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="TRESS_BestOption"),
                        dt.exomePeak2.BestCombo.HEK.FPR20.peaks.nonNEB %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="exomePeak2_BestCombo")) %>%
  mutate(Group=factor(Group, levels=c("TRESS_BestOption","exomePeak2_BestCombo"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)

for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  exomePeak2.bestcombo <- plotRanges(data = dt.peak.region, params = region, order = "random",
                                     fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),strandSplit = F,
                                     x = 0.5, y = 1.5+(s-1+length(StudySamples))*track.height, width = 5, height = track.height,
                                     just = c("left", "top"), default.units = "cm"
  )
  sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5+length(StudySamples))*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
  sample_rect <- plotRect(x=0.5,y=1.5+(s-1+length(StudySamples))*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#track for benchmark m6A sites (perturbation effect)
#FTO GLORI HEK293
signal_confident.benchmark.FTO <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.FTO,name,strand) %>% filter(strand==dt.selected.gene$strand & !is.na(score)), params = region,
                                             x = 0.5, y = 1.5+2*length(StudySamples)*track.height, width = 5, height = track.height,
                                             just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+0.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="FTO treatment", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+0.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+2*length(StudySamples)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
#METTL3i GLORI HEK293
signal_confident.benchmark.METTL3i <- plotSignal(data = dt.GLORI.HEK293.confident %>% dplyr::select(seqnames,start,end,score=PerturbEffect.METTL3i), params = region,
                                                 x = 0.5, y = 1.5+(2*length(StudySamples)+1)*track.height, width = 5, height = track.height,
                                                 just = c("left", "top"), default.units = "cm",linecolor = "#68bd48",fill="#68bd48",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
plotText(label = paste0("[0-1]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+1.8)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
sample_label <- plotText(label="METTL3 inhibition", fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+1.5)*track.height,x=5.6,draw=T,just = "left", fontcolor = "black",default.units = "cm")
sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+1)*track.height,width=5, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fontface="plain")
legendPlot.confident.benchmark <- plotLegend(legend = c("Pertubation delta m6A level in GLORI_HEK293"), fill = c("#68bd48"),border = F, x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height,
                                             width = 2.5, height = 1,
                                             just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = NULL,fontface="plain")
legendPlot.exomePeak2.peak <- plotLegend(legend = c("TRESS_BestOption","exomePeak2_BestCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 1.5, height = 1,
                                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                         title = expression("PeakGroup"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


#save intermediate variables for re-use
save.image("sm6APeak_FigureS3.intermediate.results.RDS")
# q(save="no")
####################### combine all figure S2 together
figS3.list <- list(p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB=p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB,
                   p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB=p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB,
                   p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB=p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB,
                   p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB=p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB,
                   p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB=p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB,
                   p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB= p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB,#this plot need to be output in raster pdf
                   p.AUC.nonNEB.all.option=p.AUC.nonNEB.all.option,
                   p.pruned.TPR.FPR.nonNEB.Method.BestOption=p.pruned.TPR.FPR.nonNEB.Method.BestOption,
                   p.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB=p.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB,
                   p.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB=p.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB
)
saveRDS(figS3.list, file="FigureS3.plot.list.RDS")


ggsave(plot=p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 4.5,filename = "FigureS3_p.bothstrand.Method.BestCombo.FPR20.FP.nonNEB.pdf")
ggsave(plot=p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 4.5,filename = "FigureS3_p.bothstrand.Method.BestCombo.FPR20.FN.nonNEB.pdf")

ggsave(plot=p.Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB, dpi = 300,units = "cm",device = cairo_pdf,
       width=5.8,height = 4.5, filename = "FigureS3_Method.BestCombo.FPR20.m6A.TP.FN.cor.intensity.benchmark.nonNEB.pdf")


pdf("FigureS3.Stranded_intensity_facilitate more_accurate_and_sensitive_peak_calling_regardless_antim6A.pdf",width = 8.2677, height = 11.693)

pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Method.BestCombo.FPR20.m6A.TP.FP.FN.count.nonNEB, x = 0.05, y=0.2, default.units = "cm",width = 6, height = 4.5)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 6.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x = 6.2, y=0.2, width = 5.8, height = 4.5,just = c("top","left"), default.units = "cm",draw=T)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 12.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x = 12.2, y=0.2, width = 5.8, height = 4.5,just = c("top","left"), default.units = "cm",draw=T)
#row2
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.bothstrand.enrichment.ratio.Method.BestCombo.FPR20.m6A.TP.FP.FN.nonNEB, x = 0.05, y=4.7, default.units = "cm",width = 8, height = 4.5)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 8.2, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Method.BestCombo.FPR20.m6A.FP.reverseenriched.Benchmarkoverlap.ratio.nonNEB, x = 8.2, y=4.7, default.units = "cm",width = 4, height = 4.5)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 12.2, y = 4.7,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x = 12.2, y=4.9, width = 5.8, height = 4.5,just = c("top","left"), default.units = "cm",draw=T)
#row3
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05, y = 9.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.AUC.nonNEB.all.option, x = 0.05, y=9.2, default.units = "cm",width = 10, height = 3.8)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 11.5, y = 9.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TRESS.BestOption.specific.ratio.m6A.peak.nonNEB, x = 11.5, y=9.2, default.units = "cm",width = 6, height = 3.8)
#row4
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05, y = 13.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.pruned.TPR.FPR.nonNEB.Method.BestOption, x = 0.05, y=13.2, default.units = "cm",width = 10, height = 4.5)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 10.2, y = 13.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.TRESS.BestOption.specific.m6Agene.with.benchmark.nonNEB, x = 10.2, y=13.2, default.units = "cm",width = 6.2, height = 4.5)
#row5
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05, y = 17.7,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x = 0.05, y = 17.9, width = 7.9, height = 5.9,just = c("top","left"), default.units = "cm",draw=T)
plotText(label = "L", fontsize = 8, fontface = "bold",x = 9.5, y = 17.7,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x = 9.5, y = 17.9, width = 7.9, height = 5.9,just = c("top","left"), default.units = "cm",draw=T)
dev.off()


