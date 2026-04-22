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
source("/data/m6A_calling_strategy/Script/Wrapper_function_for_m6A_calling_using_MACS2SMePeak.R")


#### Figure 5 sm6APeak identified extensive aberrant m6A modification  mediated transcription dysregulation in CRC #########
############5.1 Compared with exomePeak2/MeTPeak/MeRIPtools, sm6APeak identify extensive novel m6A+ gene in (GSE179042, another dataset was used in validation) ##############
#the count of peak and m6AGene of sm6APeak compared with  those from selected method with same level FPR
#load peaks
dt.m6a.all <- readRDS("/data2/CRC_m6A/Analysis/dt.CRC.m6a.combined.RDS")
dt.m6a.combined <- dt.m6a.all %>% filter(StudyID=="GSE179042") %>%
  mutate(Origin=factor(Origin,levels=c("M6APeakS","exomePeak2","MeTPeak","MeRIPtools"),labels=c("sm6APeak","exomePeak2","MeTPeak","MeRIPtools")))
#load study and sample info
dt.CRC.sample <- fread("/data2/CRC_m6A/rawdata/CRC_m6A_Sample_info.csv") %>% mutate(DonorName=SampleName %>% str_replace(pattern="_Input",replacement = "") %>% str_replace(pattern = "_RIP",replacement = ""))
dt.m6a.combined <- dt.m6a.combined %>% left_join(x=.,y=dt.CRC.sample %>% dplyr::distinct(DonorName,StudyID), by=c("Sample"="DonorName","StudyID"))
dt.m6a.combined %>% distinct(Origin,StudyID)
dt.m6a.combined <- dt.m6a.combined %>% dplyr::filter(str_count(OverlappedGenes,",")<=2)
#count of peak and m6A gene
dt.m6a.peak.gene.count <- dt.m6a.combined %>% dplyr::filter(str_detect(Sample, "CRC") & Origin != "GSE247632") %>% separate_longer_delim(cols = c("OverlappedGenes"), delim = ",") %>%
  group_by(Sample,Origin) %>% mutate(nPeak=n_distinct(name), nGene=n_distinct(OverlappedGenes[!is.na(OverlappedGenes)])) %>% as.data.table() %>% distinct(StudyID,Sample,Origin,nPeak,nGene) %>%
  mutate(Condition=ifelse(str_detect(Sample,"NT"),"Adjacent","Tumor")) %>% dplyr::arrange(StudyID,Condition,Sample,Origin)
dt.m6a.peak.gene.count <- dt.m6a.peak.gene.count %>% mutate(Method=Origin) %>%
  pivot_longer(cols = c("nPeak","nGene"), names_to = "Metric", values_to = "Count") %>% as.data.table() %>%
  mutate(Metric=factor(Metric,levels=c("nPeak","nGene"), labels=c("m6A peak", "m6A+ gene")), Sample=factor(Sample,levels=unique(Sample)))
dt.m6a.peak.gene.count <- dt.m6a.peak.gene.count %>% group_by(StudyID,Method,Metric) %>% mutate(avgCount=as.integer(mean(Count))) %>% as.data.table() %>% mutate(Label=paste(paste0("AvgN = ",avgCount)))
library(ggplot2)
library(tidyr)
library(dplyr)

#peak.count.test.res <- wilcox.test(Count ~ Method, data = dt.m6a.peak.gene.count %>% dplyr::filter(Metric=="m6A peak"), paired = FALSE)
p.m6a.peak.count.sm6APeakvsOrigin <- dt.m6a.peak.gene.count %>% dplyr::filter(Metric=="m6A peak") %>%
  ggplot(data=., aes(x = Method, y = Count/1000)) +
  geom_col(data=dt.m6a.peak.gene.count %>% dplyr::filter(Metric=="m6A peak") %>% distinct(Method,avgCount,Label),
           aes(x=Method,y=avgCount/1000, fill=Method),alpha=0.7, width = 0.6,show.legend = FALSE)+
  geom_text(data=dt.m6a.peak.gene.count %>% dplyr::filter(Metric=="m6A peak") %>% distinct(Method,avgCount,Label),
            aes(x=Method,y=3.5,label=Label),size=2,angle=90)+
  geom_line(aes(group = Sample), color = "gray40", alpha = 0.3, linewidth = 0.6) +
  geom_jitter(aes(fill = Method), width = 0.1, height = 0,
              shape = 21, size = 1., stroke = 0.3, show.legend = FALSE) +
  labs(y = "Count (x1000) of m6A peak", x=NULL)+
  scale_fill_manual(values=c("#a03d48", "#efb091", "#c0c1e3", "#b5cfe2"),breaks = c("sm6APeak","exomePeak2","MeTPeak","MeRIPtools")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),n.breaks = 6) +
  # facet_wrap(~StudyID,scales = "free_y")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.m6a.gene.count.sm6APeakvsOrigin <- dt.m6a.peak.gene.count %>% dplyr::filter(Metric!="m6A peak") %>%
  ggplot(data=., aes(x = Method, y = Count/1000)) +
  geom_col(data=dt.m6a.peak.gene.count %>% dplyr::filter(Metric!="m6A peak") %>% distinct(Method,avgCount,Label),
           aes(x=Method,y=avgCount/1000, fill=Method),alpha=0.7, width = 0.6,show.legend = FALSE)+
  geom_text(data=dt.m6a.peak.gene.count %>% dplyr::filter(Metric!="m6A peak") %>% distinct(Method,avgCount,Label),
            aes(x=Method,y=3,label=Label),size=2,angle=90)+
  geom_line(aes(group = Sample), color = "gray40", alpha = 0.3, linewidth = 0.6) +
  geom_jitter(aes(fill = Method), width = 0.1, height = 0,
              shape = 21, size = 1., stroke = 0.3, show.legend = FALSE) +
  labs(y = "Count (x1000) of m6A+ gene", x=NULL)+
  scale_fill_manual(values=c("#a03d48", "#efb091", "#c0c1e3", "#b5cfe2"),breaks = c("sm6APeak","exomePeak2","MeTPeak","MeRIPtools")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),n.breaks = 6) +
  # facet_wrap(~StudyID,scales = "free_y")+
  ggpubr::theme_pubr(base_size = 7,base_family = "Times") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.m6a.peak.count.sm6APeakvsOrigin.combined <- cowplot::plot_grid(a=p.m6a.peak.count.sm6APeakvsOrigin,b=p.m6a.gene.count.sm6APeakvsOrigin,nrow = 1,align = "tblr")

#the count of strictly novel peaks and m6A+ gene by sm6APeaks
#identification of strictly novel peak
dt.novel.m6A.M6APeakS <- foreach(study =c("GSE247632","GSE179042")[2], .combine='rbind')%do%{
  foreach(sam1 = dt.m6a.combined[StudyID ==study,unique(Sample)], .combine = 'rbind')%do%{
    dt.m6a.new <- dt.m6a.combined %>% dplyr::filter(Sample==sam1 & Origin=="M6APeakS")
    foreach(method = c("exomePeak2","MeTPeak","MeRIPtools"), .combine = 'rbind')%do%{
      message(paste0("Detect novel m6a of M6APeakS vs ", method, " in ", sam1))
      foreach(sam2 = dt.m6a.combined[StudyID ==study,unique(Sample)], .combine = 'rbind')%do%{
        dt.m6a.old <- dt.m6a.combined %>% dplyr::filter(Sample==sam2 & Origin==method)
        # dt.m6a.novel <- bt.intersect(a=dt.m6a.new, b=dt.m6a.old[,1:6], s=T, f=0.5, v=T) %>% as.data.table()
        # colnames(dt.m6a.novel) <- colnames(dt.m6a.new)
        dt.m6a.novel <- bt.intersect(a=dt.m6a.new, b=dt.m6a.old[,1:6], s=T, f=0.05, v=T) %>% as.data.table()
        colnames(dt.m6a.novel) <- colnames(dt.m6a.new)
        dt.m6a.novel %>% mutate(Sample=sam1, OldSample=sam2,ComparedMethod=method)
      }
    }
  }
}
# saveRDS(dt.novel.m6A.M6APeakS, file="dt.novel.m6A.M6APeakS.RDS")
#strict novel m6a peak (peak novel in at least 90% of total sample with old peak)
# dt.strict.novel.m6A.M6APeakS <- dt.novel.m6A.M6APeakS %>% group_by(name,Sample,StudyID,ComparedMethod) %>% mutate(nOldSample=n_distinct(OldSample)) %>%
#   as.data.table() %>% distinct(name,Sample,OverlappedGenes,StudyID,ComparedMethod,nOldSample) %>% dplyr::arrange(StudyID,ComparedMethod,desc(nOldSample)) %>%
#   left_join(x=.,y=dt.m6a.combined %>% distinct(StudyID,Sample) %>% group_by(StudyID) %>% mutate(TotalSample=n_distinct(Sample)) %>% as.data.table() %>% distinct(StudyID,TotalSample),by=c("StudyID")) %>%
#   dplyr::filter(nOldSample>=0.9*TotalSample) #novel in at least 90%  sample under corresponding compared method
# dt.strict.novel.m6A.M6APeakS <- dt.strict.novel.m6A.M6APeakS  %>% group_by(name,StudyID) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>% dplyr::arrange(StudyID,desc(nSample))
# dt.strict.novel.m6A.M6APeakS[,list(nNovelPeak=n_distinct(name)),by=.(StudyID,ComparedMethod)]

dt.strict.novel.m6A.M6APeakS.deredundant <- readRDS("/data2/CRC_m6A/Analysis/dt.strict.novel.m6A.M6APeakS.deredundant.RDS") %>% dplyr::filter(StudyID=="GSE179042")
#assess the detection of strict novel m6a in two cohort
dt.unique.m6a <- dt.strict.novel.m6A.M6APeakS.deredundant %>% distinct(seqnames,start,end,name,score,strand)
dt.strict.novel.m6A.M6APeakS.deredundant.detailed <- foreach(study=unique(dt.strict.novel.m6A.M6APeakS.deredundant$StudyID),.combine='rbind')%do%{
  AllSamples <- unique(dt.m6a.combined[StudyID==study, Sample])
  dt.harbor.absent <- foreach(s = AllSamples, .combine ='rbind')%do%{
    dt.overlap <- bt.intersect(a=dt.unique.m6a, b=dt.m6a.combined %>% dplyr::filter(StudyID==study & Sample==s & Origin=="M6APeakS"), s=T, f=0.9, wo=T) %>% as.data.table()
    dt.absent <- bt.intersect(a=dt.unique.m6a, b=dt.m6a.combined %>% dplyr::filter(StudyID==study & Sample==s & Origin=="M6APeakS"), s=T, f=0.1, wo=T, v=T) %>% as.data.table()
    dt.unique.m6a %>% mutate(Sample=s, Harbor=name %in% dt.overlap$V4, Absent=name %in% dt.absent$V4)
  }
  dt.harbor.absent %>% group_by(name) %>% mutate(StudyID=study,TotalSample=length(AllSamples), HarborSample=n_distinct(Sample[Harbor==TRUE]), AbsentSample=n_distinct(Sample[Absent==TRUE])) %>%
    as.data.table() %>% distinct(seqnames,start,end,strand,name,StudyID,TotalSample,HarborSample,AbsentSample) %>% as.data.table() %>% dplyr::arrange(desc(HarborSample),AbsentSample)
}
#annotate the strict.novel.m6A.M6APeakS.deredundant.detailed to gene
dt.strict.novel.m6A.M6APeakS.deredundant.detailed.annotate <- annot_peak(peak.bed=dt.strict.novel.m6A.M6APeakS.deredundant.detailed %>% mutate(score=1000) %>%  dplyr::distinct(seqnames,start,end,name,score,strand),
                                                                strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                annot_type = "gene")
dt.strict.novel.m6A.M6APeakS.deredundant <- dt.strict.novel.m6A.M6APeakS.deredundant %>% left_join(x=.,y=dt.strict.novel.m6A.M6APeakS.deredundant.detailed %>%
                                                                                                     mutate(HarborRatio=round(HarborSample/TotalSample,2),AbsentRatio=round(AbsentSample/TotalSample,2)) %>%
                                                                                                     distinct(name,StudyID,HarborSample,HarborRatio,AbsentRatio),by=c("name","StudyID"))
dt.strict.novel.m6A.M6APeakS.deredundant <- dt.strict.novel.m6A.M6APeakS.deredundant %>% left_join(x=.,y=dt.strict.novel.m6A.M6APeakS.deredundant.detailed.annotate %>% dplyr::select(name,OverlappedGenes),by="name")

dt.strict.novel.m6A.M6APeakS.deredundant %>% dplyr::filter(str_detect(OverlappedGenes,"PTRH1")) %>% dplyr::arrange(ComparedMethod,desc(HarborRatio)) %>% filter(HarborRatio>=1/3)
# seqnames     start       end                       name score strand   StudyID ComparedMethod HarborSample HarborRatio AbsentRatio OverlappedGenes
# 1:     chr9 127695052 127713553 chr9_127695052_127713553_-  1000      - GSE179042     MeRIPtools           10         0.5        0.50           PTRH1
# 2:     chr9 127711888 127714837 chr9_127711888_127714837_-  1000      - GSE179042     MeRIPtools           10         0.5        0.45           PTRH1
# 3:     chr9 127695052 127713553 chr9_127695052_127713553_-  1000      - GSE179042        MeTPeak           10         0.5        0.50           PTRH1
# 4:     chr9 127711888 127714837 chr9_127711888_127714837_-  1000      - GSE179042        MeTPeak           10         0.5        0.45           PTRH1
# 5:     chr9 127695052 127713553 chr9_127695052_127713553_-  1000      - GSE179042     exomePeak2           10         0.5        0.50           PTRH1
# 6:     chr9 127711888 127714837 chr9_127711888_127714837_-  1000      - GSE179042     exomePeak2           10         0.5        0.45           PTRH1
dt.strict.novel.m6A.M6APeakS.deredundant %>% dplyr::filter(str_detect(OverlappedGenes,"C5AR1")) %>% dplyr::arrange(ComparedMethod,desc(HarborRatio)) %>% filter(HarborRatio>=1/3)
#     seqnames    start      end                      name score strand   StudyID ComparedMethod HarborSample HarborRatio AbsentRatio OverlappedGenes
# 1:    chr19 47307465 47322066 chr19_47307465_47322066_+  1000      + GSE179042     MeRIPtools            8        0.40        0.05           C5AR1
# 2:    chr19 47290022 47307465 chr19_47290022_47307465_+  1000      + GSE179042     MeRIPtools            7        0.35        0.60           C5AR1
# 3:    chr19 47290022 47307465 chr19_47290022_47307465_+  1000      + GSE179042        MeTPeak            7        0.35        0.60           C5AR1
# 4:    chr19 47290022 47307465 chr19_47290022_47307465_+  1000      + GSE179042     exomePeak2            7        0.35        0.60           C5AR1

#replicate in another cohort
dt.m6a.all <- readRDS("/data2/CRC_m6A/Analysis/dt.CRC.m6a.combined.RDS")
dt.m6a.validation <- dt.m6a.all %>% dplyr::filter(StudyID=="GSE271081" & Origin=="M6APeakS")

# dt.overlap <- bt.intersect(a=dt.unique.m6a, b=dt.m6a.combined %>% dplyr::filter(StudyID==study & Sample==s & Origin=="M6APeakS"), s=T, f=0.9, wo=T) %>% as.data.table()
# dt.absent <- bt.intersect(a=dt.unique.m6a, b=dt.m6a.combined %>% dplyr::filter(StudyID==study & Sample==s & Origin=="M6APeakS"), s=T, f=0.1, wo=T, v=T) %>% as.data.table()
dt.strict.novel.m6A.M6APeakS.deredundant.validation <- bt.intersect(a=dt.strict.novel.m6A.M6APeakS.deredundant %>%
                                                                      # filter(str_detect(OverlappedGenes,"PTRH1") | str_detect(OverlappedGenes,"C5AR1")) %>%
                                                                      distinct(seqnames,start,end,name,score,strand),
                                                                    b=dt.m6a.validation %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                                                    s=T, f=0.25,wo=T) %>% as.data.table() %>% dplyr::arrange(V4,V10,V13) %>%
    group_by(name=V4) %>% mutate(nValidation=n_distinct(V13)) %>% as.data.table() %>% mutate(ValidationRatio=nValidation/n_distinct(dt.m6a.validation$Sample)) %>% dplyr::distinct(name,ValidationRatio,nValidation)
dt.strict.novel.m6A.M6APeakS.deredundant.validation[,.N,by=.(nValidation)]
dt.strict.novel.m6A.M6APeakS.deredundant.validation <- dt.strict.novel.m6A.M6APeakS.deredundant %>% left_join(x=.,y=dt.strict.novel.m6A.M6APeakS.deredundant.validation,by=c("name")) %>%
                                                       mutate(ValidationRatio=case_when(is.na(ValidationRatio) ~ 0, .default = ValidationRatio), nValidation=case_when(is.na(nValidation) ~ 0, .default = nValidation))
dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% filter(str_detect(OverlappedGenes,"PTRH1") | str_detect(OverlappedGenes,"C5AR1")) %>% filter(ValidationRatio>=0.5)
dt.strict.novel.m6A.M6APeakS.deredundant.validation[HarborRatio>=1/3,.(nPeak=n_distinct(name)),by=.(ComparedMethod,ValidationRatio>=0.5)]
dt.strict.novel.m6A.M6APeakS.deredundant.validation <- dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>% as.data.table()

#determine the count of these strict novel m6A located gene that were not currently (m6A+ in less than 10% samples) m6A modification in old method
dt.m6AGene.oldmethod <- dt.m6a.combined %>% filter(Origin %in% c("MeRIPtools","exomePeak2","MeTPeak") & StudyID == "GSE179042" & !is.na(OverlappedGenes)) %>%
  dplyr::distinct(Origin,Sample,OverlappedGenes) %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>% as.data.table() %>% distinct(Origin,Sample,OverlappedGenes)
dt.m6AGene.oldmethod <- dt.m6AGene.oldmethod %>% group_by(OverlappedGenes,Origin) %>% mutate(nSample=n_distinct(Sample)) %>% as.data.table() %>% distinct(Origin,OverlappedGenes,nSample)
dt.m6AGene.oldmethod <- dt.m6AGene.oldmethod %>% mutate(ratio.m6AGene=nSample/n_distinct(dt.m6a.combined$Sample))
dt.m6AGene.oldmethod %>% filter(ratio.m6AGene>0.1)#29701
dt.m6AGene.oldmethod %>% filter(OverlappedGenes %in% c("PTRH1","C5AR1"))

dt.strict.novel.m6A.M6APeakS.deredundant.validation <- dt.strict.novel.m6A.M6APeakS.deredundant.validation %>%
  left_join(x=.,y=dt.m6AGene.oldmethod %>% dplyr::select(Origin,OverlappedGenes,ratio.m6AGene),by=c("ComparedMethod"="Origin", "OverlappedGenes")) %>%
  mutate(ratio.m6AGene=case_when(is.na(ratio.m6AGene) ~ 0,.default = ratio.m6AGene))
dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% filter(OverlappedGenes %in% c("PTRH1","C5AR1") & HarborRatio>=1/3 & ValidationRatio>=0.5 )

#the count of strict novel peaks (<5% overlap with old peaks from at least 90% samples) and deredundant for highly overlapped peaks
dt.count.strict.novel.concurrent.m6a.peak <- dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% filter(HarborRatio>=1/3) %>% group_by(ComparedMethod,Validated=ValidationRatio>=0.5) %>%
  mutate(nPeak=n_distinct(name)) %>% as.data.table() %>% distinct(ComparedMethod,Validated,nPeak)
dt.count.strict.novel.concurrent.m6a.peak <- dt.count.strict.novel.concurrent.m6a.peak %>% group_by(ComparedMethod) %>% mutate(TotalPeak=sum(nPeak)) %>% filter(Validated==TRUE) %>%
                                              mutate(ValidationPeakRatio=nPeak/TotalPeak*100, Label=paste0(nPeak,"/",TotalPeak)) %>% as.data.table() %>%
                                             dplyr::mutate(ComparedMethod=factor(ComparedMethod,levels=c("exomePeak2","MeTPeak","MeRIPtools")))
p.count.strict.novel.concurrent.m6a.peak <- dt.count.strict.novel.concurrent.m6a.peak %>%
  ggplot(data=., aes(x = ComparedMethod, y = ValidationPeakRatio)) +
  geom_col(aes(fill=ComparedMethod),alpha=0.7, width = 0.6,show.legend = FALSE)+
  geom_text(aes(x = ComparedMethod, y = 30, label=Label),size=2,angle=90)+
  geom_text(aes(x = ComparedMethod, y = ValidationPeakRatio, label=paste0(round(ValidationPeakRatio,1),"%")),size=2,angle=0)+
  labs(y = str_wrap("% of strict novel and concurrent peak replicated in GSE271081",35), x="Method compared with sm6APeak")+
  scale_fill_manual(values=c("#a03d48", "#efb091", "#c0c1e3", "#b5cfe2")[-1],breaks = c("sm6APeak","exomePeak2","MeTPeak","MeRIPtools")[-1]) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)),n.breaks = 6) +
  # facet_wrap(~StudyID,scales = "free_y")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#the count of located genes of validated strict novel peaks(<5% overlap with old peaks from at least 90% samples) and deredundant for highly overlapped peaks
dt.count.strict.novel.concurrent.m6a.peak.locatedgene <- dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% filter(HarborRatio>=1/3 & ValidationRatio>=0.5) %>%
  group_by(ComparedMethod) %>% mutate(nGene=n_distinct(OverlappedGenes[!is.na(OverlappedGenes)])) %>% as.data.table() %>% distinct(ComparedMethod,nGene) %>%
  mutate(ComparedMethod=factor(ComparedMethod,levels=c("exomePeak2","MeTPeak","MeRIPtools")))

p.count.strict.novel.concurrent.m6a.peak.locatedgene <- dt.count.strict.novel.concurrent.m6a.peak.locatedgene %>%
  ggplot(data=., aes(x = ComparedMethod, y = nGene)) +
  geom_col(aes(fill=ComparedMethod),alpha=0.7, width = 0.6,show.legend = FALSE)+
  geom_text(aes(label=nGene),size=2,angle=0)+
  labs(y = str_wrap("N of replicated strict novel and concurrent peak located gene",35), x="Method compared with sm6APeak")+
  scale_fill_manual(values=c("#a03d48", "#efb091", "#c0c1e3", "#b5cfe2")[-1],breaks = c("sm6APeak","exomePeak2","MeTPeak","MeRIPtools")[-1]) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)),n.breaks = 6) +
  # facet_wrap(~StudyID,scales = "free_y")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#the venn plot display the count of gene harbor those strict novel and concurrent peaks, especially the COSMIC driver gene
dt.strict.novel.concurrent.m6a.gene <- dt.strict.novel.m6A.M6APeakS.deredundant.validation %>% filter(HarborRatio>=1/3 & ValidationRatio>=0.5 & ratio.m6AGene<=0.1) %>%
  filter(!is.na(OverlappedGenes)) %>% distinct(ComparedMethod,OverlappedGenes,HarborRatio, ValidationRatio, ratio.m6AGene)
dt.strict.novel.concurrent.m6a.gene[,.N,by=.(ComparedMethod)]
d_fname = file.path(system.file('extdata', package='TargetOsteoAnalysis'), "CosmicCensus.csv.gz")
dt.cosmic.cancer.gene.census <- read_csv(d_fname) %>% as.data.table()
colnames(dt.cosmic.cancer.gene.census) <- colnames(dt.cosmic.cancer.gene.census) %>% str_replace_all(pattern=" ",replacement = "")
dt.cosmic.cancer.gene.census %>% dplyr::filter(GeneSymbol=="CD74")
dt.strict.novel.concurrent.m6a.gene[OverlappedGenes %in% dt.cosmic.cancer.gene.census$GeneSymbol,] %>% View()
dt.strict.novel.concurrent.m6a.gene %>% filter(OverlappedGenes %in% c("PRKG2","C5AR1","PCNA","PTRH1","LMNB2","LGALS4","SYNE3","PRIM1","MCM10","PAPSS2")) %>% dplyr::arrange(OverlappedGenes)
strict.novel.concurrent.m6a.gene.list <- NULL
strict.novel.concurrent.m6a.gene.list$exomePeak2 <- dt.strict.novel.concurrent.m6a.gene[ComparedMethod=="exomePeak2",unique(OverlappedGenes)]
strict.novel.concurrent.m6a.gene.list$MeTPeak <- dt.strict.novel.concurrent.m6a.gene[ComparedMethod=="MeTPeak",unique(OverlappedGenes)]
strict.novel.concurrent.m6a.gene.list$MeRIPtools <- dt.strict.novel.concurrent.m6a.gene[ComparedMethod=="MeRIPtools",unique(OverlappedGenes)]
# strict.novel.concurrent.m6a.gene.list$COSMIC.CRC <- dt.cosmic.cancer.gene.census[Tier==1 & str_detect(`TumourTypes(Germline)`,"color") | str_detect(`TumourTypes(Somatic)`,"color"),unique(GeneSymbol)]
# strict.novel.concurrent.m6a.gene.list$COSMIC.UnkownCRC <- dt.cosmic.cancer.gene.census[Tier==1 & !GeneSymbol %in% strict.novel.concurrent.m6a.gene.list$COSMIC.CRC,unique(GeneSymbol)]
strict.novel.concurrent.m6a.gene.list$COSMIC <- dt.cosmic.cancer.gene.census[Tier==1,unique(GeneSymbol)]

names(strict.novel.concurrent.m6a.gene.list) <- paste0(names(strict.novel.concurrent.m6a.gene.list),"\n(N=",sapply(strict.novel.concurrent.m6a.gene.list,length),")")
Reduce(x=strict.novel.concurrent.m6a.gene.list[1:3],intersect) %>% str()#N=636
# Reduce(x=strict.novel.concurrent.m6a.gene.list[1:4],intersect) %>% str()#N=2
# Reduce(x=strict.novel.concurrent.m6a.gene.list[c(1:3,5)],intersect) %>% str()#N=24

require(ggVennDiagram)
p.Strict.Novel.m6AGene <- ggVennDiagram(strict.novel.concurrent.m6a.gene.list[1:3],label = c("count","percent","both","none")[1],
                                                              category.names=names(strict.novel.concurrent.m6a.gene.list)[-4],
                                                              set_color = c("#efb091", "#c0c1e3", "#b5cfe2","#b9dfe0","#71b7ab")[1:3],
                                                              edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[c(1)],
                                                              set_size = 2,label_size = 2,label_alpha = 0.6,
                                                              edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  labs(title = "Strict novel m6A+ gene via sm6APeak")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#the function enrichment of those novel m6A+ gene
dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

#the IGV plots of novel m6A in C5AR1
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name %in% c("PTRH1","C5AR1")[2])

pdf("Figure5_subplot_strict_novel_m6A_in_C5AR1_in_CRC.pdf",width = 8.2677, height = 11.693)

## Create a page (7.5*7.5cm)
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3.2
pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
# region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = 47318*1000, chromend = 47323*1000,assembly = "hg38")
#plot gene
# paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38", width=3.0)
## Plot genes small
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width=track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

#track for bigwig files of study samples
StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")), "Primary_m6A","LiverMeta_m6A")[3:6]
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos and neg strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height, width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
peak_height=0.05
#add m6a peak
# Selected_m6A <- dt.mergepeak %>% filter(StudyID=="GSE179042") %>%  dplyr::filter(seqnames==dt.selected.gene$seqnames & strand == dt.selected.gene$strand) %>%
#   dplyr::filter(start>=(floor(dt.selected.gene$start/1000)-window.size)*1000 & end<= (ceiling(dt.selected.gene$end/1000)+window.size)*1000) %>%
#   dplyr::arrange(Method,Condition,seqnames,start,end) %>% dplyr::select(seqnames,start,end,name,score,strand,everything())## Define region
# Selected_m6A <- Selected_m6A %>% distinct(seqnames,start,end,name,score,strand,Method,Condition)
# dt.peak.track <- data.table(Condition=rep(c("NT","T"),each=4),
#                             Method=rep(c("MeRIPtools","MeTPeak","exomePeak2","M6APeakS"), 2)) %>%
#   mutate(PeakName=paste(Condition,Method,sep=".")) %>%
#   mutate(trackcolor=case_match(Method,"M6APeakS" ~ "#a03d48",  "exomePeak2" ~ "#efb091","MeTPeak" ~ "#c0c1e3", "MeRIPtools" ~ "#b5cfe2"))
# dt.peak.track <- dt.peak.track %>% distinct(Method,Condition,trackcolor) %>% mutate(PeakName=paste(Method,Condition,sep="."))

Selected_m6A <- dt.m6a.combined
# Selected_m6A <- dt.m6a.all
Selected_m6A <- Selected_m6A %>% filter(strand==as.character(dt.selected.gene$strand)) %>% mutate(Method=Origin) %>% distinct(seqnames,start,end,name,score,strand,Method,Sample) %>%
  mutate(Sample=factor(Sample,levels=c(paste0("CRC",1:10,"T"),paste0("CRC",1:10,"NT")))) %>% dplyr::arrange(Method,Sample,name)
# Selected_m6A <- Selected_m6A %>% filter(seqnames == region$chrom & start>=region$chromstart & end<= region$chromend & strand==as.character(dt.selected.gene$strand))
dt.peak.track <- data.table(Method=rep(unique(Selected_m6A$Method),each=length(unique(Selected_m6A$Sample))), Sample=rep(unique(Selected_m6A$Sample), n=n_distinct(Selected_m6A$Method))) %>%
  mutate(trackcolor=case_match(Method,"M6APeakS" ~ "#a03d48",  "exomePeak2" ~ "#efb091","MeTPeak" ~ "#c0c1e3", "MeRIPtools" ~ "#b5cfe2"),
         PeakName=paste0(Method,"_","_peak"))

for(i in (1:nrow(dt.peak.track))){
  #print(Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & StudyID==dt.peak.track$StudyID[i] & Condition==dt.peak.track$Condition[i]))
  plotRanges(data=Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & Sample==dt.peak.track$Sample[i]) %>%dplyr::select(seqnames:strand),
             params = region, just = c("left","top"), x = 0.5, y = 1.5+length(StudySamples)*track.height+(i-1)*peak_height,
             default.units = "cm", width=track.width, height=peak_height, fill=dt.peak.track$trackcolor[i],
             boxHeight = unit(0.3,"mm"),linecolor = "grey5",linewidth=0.1,
             collapse = F, order = "random",draw = T)
  if(i %in% c(10,30,50,70)){
    sample_label <- plotText(label=dt.peak.track$PeakName[i], fontsize = 6, fontface="plain",y=1.5+length(StudySamples)*track.height+(i-0.5)*peak_height,x=0.6+track.width,draw=T,just = "left", fontcolor = dt.peak.track$trackcolor[i],default.units = "cm")
  }
  sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height+(i-1)*peak_height,width=track.width, height=peak_height,linecolor = dt.peak.track$trackcolor[i], lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#add legend for MeRIP
legendPlot <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = FALSE,x =1+track.width, y = 0.1, width = 1, height = 1,
                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "v",
                         title = expression("log"["2"]~"(CPM)"),fondface="bold")
#add rect for whole panel
plotRect(x=0.05,y=0.01,width=5.7, height = 7.6,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


#the IGV plot of PTRH1 novel m6A

#load peaks
dt.m6a.combined <- readRDS("/data2/CRC_m6A/Analysis/dt.CRC.m6a.combined.RDS") %>% filter(StudyID=="GSE179042")
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name %in% c("PTRH1","C5AR1")[1])

pdf("Figure5_subplot_strict_novel_m6A_in_PTRH1_in_CRC.pdf",width = 8.2677, height = 11.693)
#PTRH1
## Create a page (7.5*7.5cm)
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3.2
pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
# region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = 47318*1000, chromend = 47323*1000,assembly = "hg38")
#plot gene
# paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38", width=3.0)
## Plot genes small
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width=track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

#track for bigwig files of study samples
StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")), "Primary_m6A","LiverMeta_m6A")[3:6]
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos and neg strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height, width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
peak_height=0.05
#add m6a peak
# Selected_m6A <- dt.mergepeak %>% filter(StudyID=="GSE179042") %>%  dplyr::filter(seqnames==dt.selected.gene$seqnames & strand == dt.selected.gene$strand) %>%
#   dplyr::filter(start>=(floor(dt.selected.gene$start/1000)-window.size)*1000 & end<= (ceiling(dt.selected.gene$end/1000)+window.size)*1000) %>%
#   dplyr::arrange(Method,Condition,seqnames,start,end) %>% dplyr::select(seqnames,start,end,name,score,strand,everything())## Define region
# Selected_m6A <- Selected_m6A %>% distinct(seqnames,start,end,name,score,strand,Method,Condition)
# dt.peak.track <- data.table(Condition=rep(c("NT","T"),each=4),
#                             Method=rep(c("MeRIPtools","MeTPeak","exomePeak2","M6APeakS"), 2)) %>%
#   mutate(PeakName=paste(Condition,Method,sep=".")) %>%
#   mutate(trackcolor=case_match(Method,"M6APeakS" ~ "#a03d48",  "exomePeak2" ~ "#efb091","MeTPeak" ~ "#c0c1e3", "MeRIPtools" ~ "#b5cfe2"))
# dt.peak.track <- dt.peak.track %>% distinct(Method,Condition,trackcolor) %>% mutate(PeakName=paste(Method,Condition,sep="."))

Selected_m6A <- dt.m6a.combined
# Selected_m6A <- dt.m6a.all
Selected_m6A <- Selected_m6A %>% filter(strand==as.character(dt.selected.gene$strand)) %>% mutate(Method=Origin) %>% distinct(seqnames,start,end,name,score,strand,Method,Sample) %>%
  mutate(Sample=factor(Sample,levels=c(paste0("CRC",1:10,"T"),paste0("CRC",1:10,"NT")))) %>% dplyr::arrange(Method,Sample,name)
# Selected_m6A <- Selected_m6A %>% filter(seqnames == region$chrom & start>=region$chromstart & end<= region$chromend & strand==as.character(dt.selected.gene$strand))
dt.peak.track <- data.table(Method=rep(unique(Selected_m6A$Method),each=length(unique(Selected_m6A$Sample))), Sample=rep(unique(Selected_m6A$Sample), n=n_distinct(Selected_m6A$Method))) %>%
  mutate(trackcolor=case_match(Method,"M6APeakS" ~ "#a03d48",  "exomePeak2" ~ "#efb091","MeTPeak" ~ "#c0c1e3", "MeRIPtools" ~ "#b5cfe2"),
         PeakName=paste0(Method,"_","_peak"))

for(i in (1:nrow(dt.peak.track))){
  #print(Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & StudyID==dt.peak.track$StudyID[i] & Condition==dt.peak.track$Condition[i]))
  plotRanges(data=Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & Sample==dt.peak.track$Sample[i]) %>%dplyr::select(seqnames:strand),
             params = region, just = c("left","top"), x = 0.5, y = 1.5+length(StudySamples)*track.height+(i-1)*peak_height,
             default.units = "cm", width=track.width, height=peak_height, fill=dt.peak.track$trackcolor[i],
             boxHeight = unit(0.3,"mm"),linecolor = "grey5",linewidth=0.1,
             collapse = F, order = "random",draw = T)
  if(i %in% c(10,30,50,70)){
    sample_label <- plotText(label=dt.peak.track$PeakName[i], fontsize = 6, fontface="plain",y=1.5+length(StudySamples)*track.height+(i-0.5)*peak_height,x=0.6+track.width,draw=T,just = "left", fontcolor = dt.peak.track$trackcolor[i],default.units = "cm")
  }
  sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height+(i-1)*peak_height,width=track.width, height=peak_height,linecolor = dt.peak.track$trackcolor[i], lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#add legend for MeRIP
legendPlot <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = FALSE,x =1+track.width, y = 0.1, width = 1, height = 1,
                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "v",
                         title = expression("log"["2"]~"(CPM)"),fondface="bold")
#add rect for whole panel
plotRect(x=0.05,y=0.01,width=5.7, height = 7.6,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name == Reduce(x=strict.novel.concurrent.m6a.gene.list[1:4],intersect)[1])
pdf("Figure5_subplot_strict_novel_m6A_in_TPM3_in_CRC.pdf",width = 8.2677, height = 11.693)
#TPM3
## Create a page (7.5*7.5cm)
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3.2
pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
# region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = 47318*1000, chromend = 47323*1000,assembly = "hg38")
#plot gene
# paramsgene <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38", width=3.0)
## Plot genes small
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width=track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

#track for bigwig files of study samples
StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")), "Primary_m6A","LiverMeta_m6A")[3:6]
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.neg %>% mutate(score=-(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos and neg strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both,-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height, width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    Input <- Input.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    RIP <- RIP.pos %>% mutate(score=(log2(score+pseudocount)-log2(pseudocount))) %>% dplyr::select(seqnames,start,end,score)
    rbind(Input %>% mutate(library="Input"), RIP %>% mutate(library="RIP")) %>% mutate(sample=StudySamples[s])
  }
  #determine the signal range for pos strand
  range_both <- ceiling(quantile(abs(dt.density$score),1))+1
  #plot
  for(s in 1:length(StudySamples)){
    signal_RIP <- plotSignal(data = dt.density %>% dplyr::filter(library=="RIP" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                             params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both), x = 0.5, y = 1.5+(s-1)*track.height, width=track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.25)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
peak_height=0.05
#add m6a peak
Selected_m6A <- dt.m6a.combined
# Selected_m6A <- dt.m6a.all
Selected_m6A <- Selected_m6A %>% filter(strand==as.character(dt.selected.gene$strand)) %>% mutate(Method=Origin) %>% distinct(seqnames,start,end,name,score,strand,Method,Sample) %>%
  mutate(Sample=factor(Sample,levels=c(paste0("CRC",1:10,"T"),paste0("CRC",1:10,"NT")))) %>% dplyr::arrange(Method,Sample,name)
# Selected_m6A <- Selected_m6A %>% filter(seqnames == region$chrom & start>=region$chromstart & end<= region$chromend & strand==as.character(dt.selected.gene$strand))
dt.peak.track <- data.table(Method=rep(unique(Selected_m6A$Method),each=length(unique(Selected_m6A$Sample))), Sample=rep(unique(Selected_m6A$Sample), n=n_distinct(Selected_m6A$Method))) %>%
  mutate(trackcolor=case_match(Method,"M6APeakS" ~ "#a03d48",  "exomePeak2" ~ "#efb091","MeTPeak" ~ "#c0c1e3", "MeRIPtools" ~ "#b5cfe2"),
         PeakName=paste0(Method,"_","_peak"))

for(i in (1:nrow(dt.peak.track))){
  #print(Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & StudyID==dt.peak.track$StudyID[i] & Condition==dt.peak.track$Condition[i]))
  plotRanges(data=Selected_m6A %>% dplyr::filter(Method==dt.peak.track$Method[i] & Sample==dt.peak.track$Sample[i]) %>%dplyr::select(seqnames:strand),
             params = region, just = c("left","top"), x = 0.5, y = 1.5+length(StudySamples)*track.height+(i-1)*peak_height,
             default.units = "cm", width=track.width, height=peak_height, fill=dt.peak.track$trackcolor[i],
             boxHeight = unit(0.3,"mm"),linecolor = "grey5",linewidth=0.1,
             collapse = F, order = "random",draw = T)
  if(i %in% c(10,30,50,70)){
    sample_label <- plotText(label=dt.peak.track$PeakName[i], fontsize = 6, fontface="plain",y=1.5+length(StudySamples)*track.height+(i-0.5)*peak_height,x=0.6+track.width,draw=T,just = "left", fontcolor = dt.peak.track$trackcolor[i],default.units = "cm")
  }
  sample_rect <- plotRect(x=0.5,y=1.5+length(StudySamples)*track.height+(i-1)*peak_height,width=track.width, height=peak_height,linecolor = dt.peak.track$trackcolor[i], lwd = 0.5,just = c("left","top"),default.units = "cm")
}

#add legend for MeRIP
legendPlot <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = FALSE,x =1+track.width, y = 0.1, width = 1, height = 1,
                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "v",
                         title = expression("log"["2"]~"(CPM)"),fondface="bold")
#add rect for whole panel
plotRect(x=0.05,y=0.01,width=5.7, height = 7.6,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

############# 5.2. sm6APeak identified aberrant m6A mediated transcription dysregulation in Tumor #################
#the count of DMEG and ratio of novel DMG
DMG.DEG.list <- readRDS("/data2/CRC_m6A/Analysis/CRC_m6A_DMG_DEG.list.RDS")
DMG.DEG.list <- DMG.DEG.list[grep(names(DMG.DEG.list),pattern="GSE179042")]
require(ggVennDiagram)
p.venn.DMG.DEG <- ggVennDiagram(DMG.DEG.list[1:4],label = c("count","percent","both","none")[1],
                                          category.names=names(DMG.DEG.list)[1:4] %>% strsplit(split=".",fixed=T) %>% sapply("[",1),
                                          set_color = c("#fcbb44","#839cd1", "#f0766d","#71b7ab"),
                                          edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[1],
                                          set_size = 2.,label_size = 2.,label_alpha = 0.6,
                                          edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#GSEA enrichment of DEG for DMG
#build TERM2NAME and TERM2GENE data frame
dt.term2name.gr <- list(data.table(gs_id="HyperDMG",gs_name="Hypermethylated genes"),
                        data.table(gs_id="HypoDMG",gs_name="Hypomethylated genes")) %>% Reduce(x=.,rbind) %>% distinct()
dt.term2gene.gr <- list(data.table(gs_id="HyperDMG",gene_name=DMG.DEG.list$HyperDMG.GSE179042),
                        data.table(gs_id="HypoDMG",gene_name=DMG.DEG.list$HypoDMG.GSE179042)) %>%
  Reduce(x=.,rbind) %>% distinct()
dt.term2gene.gr <- genekitr::transId(id = unique(dt.term2gene.gr$gene_name),transTo = c("entrez"),unique = T,org = "hs") %>% as.data.table() %>%
  inner_join(x=.,dt.term2gene.gr,by=c("input_id"="gene_name")) %>% dplyr::select(gs_id,entrezid)
dt.term2gene.gr[,.N,by=gs_id]
#load DEG
dt.DEG.two.study <- readRDS("/data2/CRC_m6A/Analysis/dt.CRC.DEG.two.study.RDS") %>% filter(StudyID=="GSE179042")
dt.DEG <- dt.DEG.two.study %>% dplyr::filter(StudyID==study) %>% dplyr::arrange(desc(LFC))
dt.DEG <- genekitr::transId(id = unique(dt.DEG$Symbol),transTo = c("entrez"),unique = T,org = "hs") %>% as.data.table() %>%
  inner_join(x=.,dt.DEG,by=c("input_id"="Symbol")) %>%
  dplyr::filter(!is.na(entrezid) & entrezid != "" & entrezid != " ") %>% dplyr::mutate(entrezid=as.character(entrezid)) %>% dplyr::arrange(desc(LFC))
input.deg.list <- dt.DEG$LFC
names(input.deg.list) <- dt.DEG$entrezid
# GSEA
library(clusterProfiler)
library(GseaVis)
library(org.Hs.eg.db)
gsea_result <- GSEA(geneList = input.deg.list,TERM2GENE = dt.term2gene.gr,  TERM2NAME = dt.term2name.gr,
                    pvalueCutoff = 1,pAdjustMethod = "BH", verbose = FALSE,maxGSSize = 2000 )
gsea_result <- setReadable(gsea_result,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
library(BiocParallel)
library(GseaVis)
register(SerialParam())#turn off parallel
p.gsea.DEGenrichDMG <- gseaNb(object = gsea_result, lineSize=0.6,geneSetID = unique(dt.term2name.gr$gs_id),
                                        subPlot = 2,termWidth = 45,legend.position = c(0.8,0.5),markTopgene = 5,
                                        addPval = F,pvalX = 0.8,pvalY=0.8,pvalSize = 0.5,pCol = "grey5",pHjust = 1,
                                        curveCol = scico::scico(palette = "vik",n = 2,direction = -1),
                                        htHeight=0.3,nesDigit=1,base_size = 6)+
  labs(x=paste0("Ranked of DEGs by descending LFC"),title = paste0("Enrichment of DEG for DMG"))+
  theme(plot.margin = margin(),legend.text = element_text(size=6), legend.key.height = unit(6,"point"), legend.key.width = unit(6,"point"),
        axis.title = element_text(size=6), legend.background = element_rect(fill="transparent"),legend.title = element_blank(),
        plot.title = element_blank())

# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = T)
# plotGG(plot = p.gsea.DEGenrichDMG, x = 0.2, y=4.9, default.units = "cm",width = 5.5, height = 4.7)

## functional enrichment of DMEG
enrichment_res.DMEG  <- readRDS("/data2/CRC_m6A/Analysis/Functiona_enrichment_res.DMEG.RDS") %>% dplyr::filter(StudyID=="GSE179042")
enrichment_res.DMEG %>% dplyr::filter(p.adjust<=0.05 & Count>=3) %>% dplyr::arrange(StudyID,DiffMethy,DiffExpr,p.adjust) %>%  group_by(StudyID,DiffMethy,DiffExpr) %>% slice_head(n=10) %>%
  as.data.table() %>% View()
enrichment_res.DMEG %>% dplyr::filter(p.adjust<=0.05 & Count>=3) %>%  group_by(DiffMethy,DiffExpr,Description) %>% filter(n_distinct(StudyID)>1) %>% mutate(avg.lgFDR=mean(-log10(p.adjust))) %>%
  as.data.table() %>% dplyr::arrange(DiffMethy,DiffExpr,desc(avg.lgFDR),Description) %>% as.data.table() %>% View()#consistently enriched terms
#visualization of functional enrichment of DMEG
mixedToFloat <- function(x){
  foreach(i=1:length(x), .combine = 'c')%do%{
    x1 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",1) %>% as.integer()
    x2 <- strsplit(x[i],split = "/",fixed=T) %>% sapply("[",2) %>% as.integer()
    x1/x2*100
  }
}
dt.top.enriched.term.DMEG <- enrichment_res.DMEG %>% dplyr::filter(p.adjust<=0.05 & Count>=3) %>%
  #dplyr::arrange(StudyID,DiffMethy,DiffExpr,p.adjust) %>% group_by(StudyID,DiffMethy,DiffExpr) %>% slice_head(n=10) %>% as.data.table() %>%
  group_by(DiffMethy,DiffExpr) %>% dplyr::arrange(p.adjust) %>% slice_head(n=8) %>% as.data.table() %>%
  mutate(GeneRatio=mixedToFloat(GeneRatio), BgRatio=mixedToFloat(BgRatio), Count=as.integer(Count), Description=paste0(Ontology,":",Description)) %>%
  mutate(EnrichFold=GeneRatio/BgRatio) %>%
  #group_by(DiffMethy,DiffExpr,Description) %>% mutate(SharedTerm=factor(n_distinct(StudyID)>1,levels=c(FALSE,TRUE))) %>%
  mutate(lgFDR= -log10(p.adjust), Description=stringr::str_wrap(Description,width = 40)) %>% as.data.table() %>%
  dplyr::arrange(DiffMethy,DiffExpr,desc(lgFDR),Description) %>% mutate(Group=paste0(DiffMethy,DiffExpr,"_DMEG"),Description=factor(Description,levels=unique(Description)))

library(ggplot2)
library(ggrepel)
require(ggsci)
require(ggpubr)
p.top.enriched.term.DMEG  <- ggplot(dt.top.enriched.term.DMEG , aes(x = EnrichFold, y = Description, color = lgFDR, size = Count)) +
  geom_point(alpha = 1) +
  # scale_color_gradient2(low = "#f0a043",mid = "#b1396d", high = "#282a60", name = "−log10(FDR)") +
  # scale_color_gradient(low = "blue", high = "red", name = "−log10(FDR)") +
  scale_color_gradient(low = "#626fb8", high = "#e6493d", name = "−log10(FDR)") +
  scale_size(range = c(1, 3), name = "Gene count") +
  scale_x_continuous(limits = c(1,10),n.breaks = 6,na.value = 10,expand = expansion(add=c(0,0.4)))+
  guides(size=guide_legend(nrow = 2))+
  labs(x="Enrichment fold",y=NULL)+
  facet_wrap( ~Group, scales = "free_y",ncol = 1) +
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

# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
# plotGG(plot = p.top.enriched.term.DMEG , x = 0.2, y=4.9, default.units = "cm",width = 6, height = 8)

#IGV of DMEG of PTRH1
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name=="PTRH1")
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3
pdf("Figure5_subplot_DMEG_of_PTRH1_in_CRC.pdf")
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
#track for bigwig files of study samples
StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")))[3:4]
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIPvsInputRatio.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIPvsInputRatio.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
#add track for gene expression
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(2*length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_expr[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(2*length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_expr[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)*3)*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)*3)*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
legendPlot.mRNA <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)*3)*track.height, width = 2.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)*3)*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#IGV of DMEG of C5AR1
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name=="C5AR1")
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3
pdf("Figure5_subplot_DMEG_of_C5AR1_in_CRC.pdf")
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
#track for bigwig files of study samples
StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")))[3:4]
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIPvsInputRatio.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIPvsInputRatio.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}
#add track for gene expression
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(2*length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_expr[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(StudySamples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=StudySamples[s],Study=strsplit(StudySamples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(StudySamples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(2*length(StudySamples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(StudySamples[s],"NT"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_expr[sample==StudySamples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(2*length(StudySamples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(2*length(StudySamples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(2*length(StudySamples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)*3)*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)*3)*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
legendPlot.mRNA <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)*3)*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)*3)*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

########################## combine plots together ############
save.image("sm6APeak_Figure5_intermediate.results.RDS")
fig5.list <- list(p.m6a.peak.count.sm6APeakvsOrigin.combined=p.m6a.peak.count.sm6APeakvsOrigin.combined,
                  p.count.strict.novel.concurrent.m6a.peak=p.count.strict.novel.concurrent.m6a.peak,
                  p.count.strict.novel.concurrent.m6a.peak.locatedgene=p.count.strict.novel.concurrent.m6a.peak.locatedgene,
                  p.Strict.Novel.m6AGene=p.Strict.Novel.m6AGene,
                  p.venn.DMG.DEG=p.venn.DMG.DEG,
                  p.gsea.DEGenrichDMG=p.gsea.DEGenrichDMG,
                  p.top.enriched.term.DMEG=p.top.enriched.term.DMEG
                  )

saveRDS(fig5.list, file="Figure5.plot.list.RDS")

pdf("Figure5.sm6APeak_identified_extensive_aberrant_m6A_modification_mediated_transcription_dysregulation_in_CRC.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.m6a.peak.count.sm6APeakvsOrigin.combined, x = 0.2, y=0.05, default.units = "cm",width = 8.4, height = 3.8)
#row2
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05, y = 4,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.strict.novel.concurrent.m6a.peak, x = 0.2, y=4, default.units = "cm",width = 4.0, height = 3.8)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 4.4, y = 4,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.strict.novel.concurrent.m6a.peak.locatedgene, x = 4.6, y=4, default.units = "cm",width = 4.0, height = 3.8)
#row3
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 8,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Strict.Novel.m6AGene, x = 0.2, y=8, default.units = "cm",width = 4.0, height = 4)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 4.2, y = 8,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.gsea.DEGenrichDMG, x = 4.4, y=8, default.units = "cm",width = 5.5, height = 5.5)
#row4
plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05, y = 12,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.venn.DMG.DEG, x = 0.2, y=12, default.units = "cm",width = 4.0, height = 4)
#row5
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05, y = 15.2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.top.enriched.term.DMEG, x = 0.2, y=15.2, default.units = "cm",width = 5, height = 8.6)


dev.off()


############# 5.3 in silico validation of DMEG  ########################################

##########  5.4 Experiment validation of DMEG of PTRH1 and C5AR1 #################################


############Figure S5  In silico validation and functional prioritation of sm6APeak identified DMEG in CRC#############
#load dt.DMEG.CRC.m6APerturbation.validation
dt.DMEG.CRC.m6APerturbation.validation <- readRDS("/data2/CRC_m6A/Analysis/dt.DMEG.CRC.m6APerturbation.validation.RDS") %>% dplyr::filter(Study=="GSE179042")
attr(dt.DMEG.CRC.m6APerturbation.validation,"index") <- NULL
PerturbRNA.DMEG.gene.list <- NULL
PerturbRNA.DMEG.gene.list$DMEG <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042",Symbol])
# PerturbRNA.DMEG.gene.list$Caco2.Perturb.m6A <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbMeRIPseq.study,"Caco2.M3KD"),Symbol])
# PerturbRNA.DMEG.gene.list$CT26.Perturb.m6A <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbMeRIPseq.study,"CT26"),Symbol])
PerturbRNA.DMEG.gene.list$HCT116.Perturb.mRNA <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbRNAseq.study,"HCT116"),Symbol])
PerturbRNA.DMEG.gene.list$Caco2.Perturb.mRNA <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbRNAseq.study,"Caco2"),Symbol])
PerturbRNA.DMEG.gene.list$CT26.Perturb.mRNA <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbRNAseq.study,"CT26"),Symbol])
PerturbRNA.DMEG.gene.list$DKO1.Perturb.mRNA <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbRNAseq.study,"DKO1"),Symbol])
names(PerturbRNA.DMEG.gene.list) <- paste0(names(PerturbRNA.DMEG.gene.list), "\n(N=", sapply(PerturbRNA.DMEG.gene.list,length),")")
p.PerturbRNA.DMEG <- ggVennDiagram(PerturbRNA.DMEG.gene.list,label = c("count","percent","both","none")[1],
                                              category.names=names(PerturbRNA.DMEG.gene.list),
                                              set_color =c("#f29f3b", "#9089c2","#8fc750","#a6539f","#337435"),
                                              edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[c(1)],
                                              set_size = 2.5,label_size = 2.5,label_alpha = 0.6,
                                              edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  labs(subtitle=str_wrap("GSE179042",width = 40),title = "DMEG validated m6A writer/eraser perturbation")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.PerturbRNA.DMEG

#b DMEG validation by m6A write/eraser pertubation in m6A level
Perturbm6A.DMEG.gene.list <- NULL
Perturbm6A.DMEG.gene.list$DMEG <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042",Symbol])
Perturbm6A.DMEG.gene.list$Caco2.Perturb.m6A <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbMeRIPseq.study,"Caco2.M3KD"),Symbol])
Perturbm6A.DMEG.gene.list$CT26.Perturb.m6A <- unique(dt.DMEG.CRC.m6APerturbation.validation[Study=="GSE179042" & str_detect(PerturbMeRIPseq.study,"CT26"),Symbol])
names(Perturbm6A.DMEG.gene.list) <- paste0(names(Perturbm6A.DMEG.gene.list), "\n(N=", sapply(Perturbm6A.DMEG.gene.list,length),")")

p.Perturbm6A.DMEG <- ggVennDiagram(Perturbm6A.DMEG.gene.list,label = c("count","percent","both","none")[1],
                                   category.names=names(Perturbm6A.DMEG.gene.list),
                                   set_color =c("#f29f3b","#a6daef", "#7b93c6"),
                                   edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[c(1)],
                                   set_size = 2.5,label_size = 2.5,label_alpha = 0.6,
                                   edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  labs(subtitle=str_wrap("GSE179042",width = 40),title = "DMEG validated m6A writer/eraser perturbation")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.Perturbm6A.DMEG

#c High ratio of DMEG validated by mRNA or m6A level under pertubation of writer or eraser
dt.count.DMEG.perturb.validated <-  dt.DMEG.CRC.m6APerturbation.validation %>% mutate(nPerturbScore=nPerturbRNAseq+nPerturbMeRIPseq) %>%
  group_by(Study,PerturbRNA=nPerturbRNAseq>0,PerturbMeRIP=nPerturbMeRIPseq>0) %>% mutate(nDMEG=n_distinct(Symbol)) %>% as.data.table() %>%
  dplyr::select(Study,DiffMethy,PerturbRNA,PerturbMeRIP,nDMEG,Symbol,contains("nPerturb")) %>% distinct()
dt.count.DMEG.perturb.validated[,.N,by=.(Study,nPerturbScore)]
#        Study nPerturbScore   N
# 9: GSE179042             7   1
# 10: GSE179042             6  29
# 11: GSE179042             5 104
# 12: GSE179042             4 109
# 13: GSE179042             3 138
# 14: GSE179042             2 145
# 15: GSE179042             1 329
# 16: GSE179042             0 138
dt.count.DMEG.perturb.validated <- dt.count.DMEG.perturb.validated %>% mutate(DiffMethy=factor(DiffMethy,levels=c("Hyper","Hypo"), labels=c("HyperUpDMEG","HypoDownDMEG"))) %>%
  mutate(Group=case_when(PerturbRNA & PerturbMeRIP ~ "mRNA & m6A", PerturbRNA ~ "mRNA", PerturbMeRIP ~ "m6A",.default = "Neither") %>%
                                                                                factor(levels=c("Neither","mRNA","m6A","mRNA & m6A")))
dt.count.DMEG.perturb.validated[,.N,by=.(Study,DiffMethy,Group)]
#label for count of both validated
dt.label.bothvalidated.count.ratio <- dt.count.DMEG.perturb.validated[,list(Count=n_distinct(Symbol)),by=.(Study,DiffMethy,Group)] %>%
  group_by(Study,DiffMethy) %>% mutate(Ratio=round(Count/sum(Count),2)) %>% as.data.table() %>% dplyr::filter(Group %in% c("mRNA & m6A")) %>%
  mutate(label=paste0(Count,"\n(",Ratio*100,"%)"), ylabel=Ratio)

#label for count of either validated
dt.label.cumulative.validated.count.ratio <- dt.count.DMEG.perturb.validated[,list(Count=n_distinct(Symbol)),by=.(Study,DiffMethy,Group)] %>%
  group_by(Study,DiffMethy) %>% mutate(Ratio=round(Count/sum(Count),2)) %>% as.data.table() %>% dplyr::filter(!Group %in% c("Neither")) %>%
  as.data.table() %>% group_by(DiffMethy) %>% mutate(CumulativeRatio=sum(Ratio),CumulativeCount=sum(Count)) %>%
  mutate(label=paste0(CumulativeCount,"\n(",CumulativeRatio*100,"%)"), ylabel=CumulativeRatio) %>% as.data.table() %>%
  distinct(Study,DiffMethy,CumulativeRatio, CumulativeCount, label, ylabel)


p.DMEG.perturb.validated <- dt.count.DMEG.perturb.validated[,list(Count=n_distinct(Symbol)),by=.(Study,DiffMethy,Group)] %>%
  ggplot(data=.)+
  geom_bar(mapping = aes(x=DiffMethy,fill=Group,y=Count), position = "fill",stat = "identity",width=0.7,alpha=0.8)+
  geom_text(data = dt.label.bothvalidated.count.ratio, aes(x=DiffMethy, y=ylabel, label=label),size=2)+
  geom_text(data = dt.label.cumulative.validated.count.ratio, aes(x=DiffMethy, y=ylabel, label=label),size=2)+
  # facet_wrap(~Study)+
  scale_fill_manual(values = c("#C1CDCD","#ef756c", "#839dcf", "#fcbb46"), breaks = c("Neither","mRNA","m6A","mRNA & m6A"))+
  guides(fill=guide_legend(title = "ValidationGroup",nrow = 2))+
  labs(x=NULL,y="Proportion of DMEGs", title="DMEGs validted by m6A perturbation")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.grid.major.x = element_blank()
  )
p.DMEG.perturb.validated

#d venn plot display overlap of perturb validated DMEG and novel m6A+ gene identified by sm6APeak
#load Novel m6A+ gene (contain strict novel peak)
StrictNovel.DMEG.list <- readRDS("dt.StrictNovel.DMEG.list.RDS")
names(StrictNovel.DMEG.list)  <- names(StrictNovel.DMEG.list) %>% str_replace(pattern="GSE179042.",replacement="")
p.venn.strict.novel.DMEG <- ggVennDiagram(StrictNovel.DMEG.list[1:3],label = c("count","percent","both","none")[1],
                                          category.names=names(StrictNovel.DMEG.list[1:3]),
                                          set_color = c("#b5cfe2","#c0c1e3","#efb091","#a03d48")[1:3],
                                          edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[c(1)],
                                          set_size = 2.5,label_size = 2.5,label_alpha = 0.6,
                                          edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  # Add rectangle to enclose all circles
  annotate("rect",xmin = -6, xmax = 11, ymin = -10, ymax = 6,  # Adjust rectangle limits as needed
    color = "#a03d48", fill = NA, size = 0.5) +
  labs(subtitle=str_wrap("GSE179042",width = 40),title  = names(StrictNovel.DMEG.list)[4])+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(colour = "#a03d48"),
    plot.subtitle = element_blank(),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p.venn.strict.novel.DMEG
#pie plot of Novel DMEG validation type
dt.selected.strict.novel.DMEG <- dt.count.DMEG.perturb.validated %>% filter(Symbol %in% Reduce(x=StrictNovel.DMEG.list[1:3],intersect)) %>%
  group_by(Group) %>% mutate(Count=n_distinct(Symbol)) %>% as.data.table() %>% distinct(Group,Count)
# Compute percentages for the pie chart
dt.selected.strict.novel.DMEG$Percent <- round(dt.selected.strict.novel.DMEG$Count / sum(dt.selected.strict.novel.DMEG$Count) * 100, 1)
dt.selected.strict.novel.DMEG$LabelPos <- dt.selected.strict.novel.DMEG$Count / 2 + cumsum(c(0, head(dt.selected.strict.novel.DMEG$Count, -1)))

# Create the pie chart
p.pie.selected.strict.novel.DMEG <- ggplot(dt.selected.strict.novel.DMEG, aes(x = "", y = Count, fill = Group)) +
  # Create the pie slices
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  # Use a custom fill color palette
  scale_fill_manual(values = c("#C1CDCD", "#ef756c", "#839dcf", "#fcbb46"), breaks = c("Neither", "mRNA", "m6A", "mRNA & m6A"),guide="none") +
  # Add percentage labels inside the pie slices
  geom_text(aes(label = paste0(Count," (",Percent, "%)")),position = position_stack(vjust = 0.5),size = 6 / .pt) +
  # Add group labels just outside the pie slices
  geom_text(aes(label = Group, x = 1.2, y = LabelPos), hjust = 0.5,size = 6 / .pt) +
  # Clean up the chart layout
  theme_minimal(base_size = 6, base_family = "Helvetica") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
p.pie.selected.strict.novel.DMEG

#e scatter plot of sm6APeak identified novel and pertubation validated DMEGs
dt.strict.novel.DMEG.perturb.validated <- dt.DMEG.CRC.m6APerturbation.validation %>% mutate(nPerturbScore=nPerturbRNAseq+nPerturbMeRIPseq) %>% dplyr::filter(nPerturbRNAseq>0 | nPerturbMeRIPseq>0) %>%
  group_by(Study,Symbol,DiffMethy) %>% mutate(NovelbyMeRIPtools=sum(StrictNovel.MeRIPtools == FALSE)==0, NovelbyMeTPeak=sum(StrictNovel.MeTPeak == FALSE)==0, NovelbyexomePeak2=sum(StrictNovel.exomePeak2 == FALSE)==0) %>%
  as.data.table() %>% dplyr::select(Study,Symbol,DiffMethy,contains("Novelby"), nPerturbScore,nPerturbRNAseq, nPerturbMeRIPseq,contains("study"),PerturbRNAseq.logP, PerturbRNAseq.LFC.abs,DiffMethyPerturb) %>% distinct() %>% dplyr::filter(NovelbyMeRIPtools | NovelbyMeTPeak | NovelbyMeTPeak)
dt.strict.novel.DMEG.perturb.validated %>% dplyr::filter(NovelbyMeRIPtools & NovelbyMeTPeak & NovelbyexomePeak2) %>% filter(Symbol %in% dt.count.DMEG.perturb.validated[Group=="mRNA & m6A",Symbol])

dt.selected.strict.novel.DMEG.perturb.validated <- dt.strict.novel.DMEG.perturb.validated %>% dplyr::filter(NovelbyMeRIPtools & NovelbyMeTPeak & NovelbyexomePeak2) %>% dplyr::arrange(desc(nPerturbScore)) %>%#45 DMEGs
                                                    filter(Symbol %in% dt.count.DMEG.perturb.validated[Group=="mRNA & m6A",Symbol])
#the scatter plot display the LFC and DiffMethy of selected strict novel DMEG validated by perturb
dt.scatterplot.strict.novel.DMEG.perturb.validated <- dt.selected.strict.novel.DMEG.perturb.validated %>% mutate(DiffMethyPerturb=case_when(is.na(DiffMethyPerturb) ~ 0, .default = DiffMethyPerturb),
                                                                                                        PerturbRNAseq.LFC.abs=case_when(is.na(PerturbRNAseq.LFC.abs) ~ 0, .default = PerturbRNAseq.LFC.abs)) %>%
  distinct(Symbol,DiffMethy,NovelbyMeRIPtools,NovelbyMeTPeak,NovelbyexomePeak2,nPerturbScore,DiffMethyPerturb,PerturbRNAseq.LFC.abs) %>%
  mutate(DiffMethy=factor(DiffMethy,levels=c("Hyper","Hypo"), labels=c("HyperUpDMEG","HypoDownDMEG")))
p.scatterplot.strict.novel.DMEG.perturb.validated <- ggplot(data=dt.scatterplot.strict.novel.DMEG.perturb.validated)+
  geom_point(aes(x=DiffMethyPerturb,y=PerturbRNAseq.LFC.abs,size=nPerturbScore,color=DiffMethy),shape=20)+
  ggrepel::geom_text_repel(data=dt.scatterplot.strict.novel.DMEG.perturb.validated %>% dplyr::filter(NovelbyMeRIPtools & NovelbyMeTPeak & NovelbyexomePeak2),
                           aes(x=DiffMethyPerturb,y=PerturbRNAseq.LFC.abs,label=Symbol,color=DiffMethy),
                           size=6 / .pt,max.overlaps = 20, fontface="italic",show.legend = F)+
  scale_color_manual(values = c("#ee6f6f", "#0d8b41"),breaks = c("HyperUpDMEG","HypoDownDMEG"))+
  scale_size_continuous(range = c(1,4),breaks = range(dt.scatterplot.strict.novel.DMEG.perturb.validated$nPerturbScore))+
  scale_y_continuous(expand = expansion(add = c(0.05,0.05)))+
  guides(size=guide_legend(nrow = 2),color=guide_legend(nrow=2))+
  labs(x=str_wrap("Absolute delta RIP/Input ratio under perturbation",width=40), y=str_wrap("Absolute LFC of mRNA under perturbation",width = 30), title="Novel perturbation validated DMEGs by M6APeakS")+
  # annotate("text",x=max(dt.scatterplot.strict.novel.DMEG.perturb.validated$DiffMethyPerturb)*0.5, y=ceiling(max(dt.scatterplot.strict.novel.DMEG.perturb.validated$PerturbRNAseq.LFC.abs))-0.2,
  #          label="HyperUpDMEG",size=2,color=c("#ee6f6f", "#0d8b41")[1])+
  # annotate("text",x=max(dt.scatterplot.strict.novel.DMEG.perturb.validated$DiffMethyPerturb)*0.8, y=ceiling(max(dt.scatterplot.strict.novel.DMEG.perturb.validated$PerturbRNAseq.LFC.abs))-0.2,
  #          label="HypoDownDMEG",size=2,color=c("#ee6f6f", "#0d8b41")[2])+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
  legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
  axis.text.y = element_text(face = "plain"),
  axis.title.x = element_text(face = "plain"),
  # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
  plot.title = element_blank(),
  plot.subtitle = element_text(size = rel(1),hjust = 0.5),
  # panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank()
  )
p.scatterplot.strict.novel.DMEG.perturb.validated

## f IGV plot of perturb validated DMEG of PTRH1 and C5AR1
dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

#IGV of DMEG of C5AR1
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name==c("C5AR1","PTRH1")[1])
dt.selected.strict.novel.DMEG.perturb.validated %>% dplyr::filter(Symbol == dt.selected.gene$gene_name) %>% dplyr::select(Symbol,PerturbRNAseq.study,PerturbMeRIPseq.study)
# Symbol     PerturbRNAseq.study           PerturbMeRIPseq.study
# 1:  PTRH1 M3KO.GSE142589.CRC.CT26 Caco2.M3KD|CT26.M3KO|CT26.M14KO

#human Caco2 MeRIPseq DMG
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name==c("C5AR1","PTRH1")[1])
pdf("Figure5S_Human_Caco2_METTL3KD_perturbation_C5AR1.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
#track for bigwig files of study samples
# StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")))[3:4]
StudySamples <- c("GSE167075_Caco2shCTL","GSE167075_Caco2shM3")
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio for Caco2.M3KD RIPvsInputRatio
m6A_samples <- c("GSE167075_Caco2shCTL","GSE167075_Caco2shM3")
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"CTL"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"CTL"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)*2)*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)*2)*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
# legendPlot.mRNA <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)*3)*track.height, width = 2.5, height = 1,
#                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)*2)*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#Add mouse CT26 m6A and mRNA under peturbation
dt.selected.gene <- dt.gene.mm39 %>% filter(gene_name=="C5ar1")
window.size=5#100kb
pseudocount=1
track.height=0.5
pdf("FigS5_Mouse_CT26_Mettl3_Mettl14_pertubation_C5ar1.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
mm39.assembly <- assembly(BSgenome = BSgenome.Mmusculus.UCSC.mm39,Genome = "mm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",OrgDb = "org.Mm.eg.db")
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly)
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

#add CT26 M3KO M14KO MeRIP (Input and RIP is strand unspecific)
StudySamples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio for Caco2.M3KD RIPvsInputRatio (strand unspecific)
m6A_samples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"control"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"control"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add track for gene expression (not strand specific)
mRNA_samples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")[-3]
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(mRNA_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",mRNA_samples[s],"_Input.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=mRNA_samples[s],Study=strsplit(mRNA_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(mRNA_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==mRNA_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==mRNA_samples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(mRNA_samples[s],"control"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_expr[sample==mRNA_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=mRNA_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+ length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(mRNA_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",mRNA_samples[s],"_Input.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=mRNA_samples[s],Study=strsplit(mRNA_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(mRNA_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==mRNA_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==mRNA_samples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(mRNA_samples[s],"control"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_expr[sample==mRNA_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=mRNA_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+ length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("Control", "Peturbation"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
legendPlot.mRNA <- plotLegend(legend = c("Control", "Peturbation"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


#IGV of DMEG of C5AR1
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name==c("C5AR1","PTRH1")[2])
dt.selected.strict.novel.DMEG.perturb.validated %>% dplyr::filter(Symbol == dt.selected.gene$gene_name) %>% dplyr::select(Symbol,PerturbRNAseq.study,PerturbMeRIPseq.study)
# Symbol     PerturbRNAseq.study           PerturbMeRIPseq.study
# 1:  PTRH1 M3KO.GSE142589.CRC.CT26 Caco2.M3KD|CT26.M3KO|CT26.M14KO

#human Caco2 MeRIPseq DMG
window.size=2#100kb
pseudocount=1
track.height=0.5
track.width=3
dt.selected.gene <- dt.gene.hg38 %>% filter(gene_name==c("C5AR1","PTRH1")[2])
pdf("Figure5S_Human_Caco2_METTL3KD_perturbation_PTRH1.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = "hg38")
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")
#track for bigwig files of study samples
# StudySamples <- c(paste0("GSE247632_", c("T","NT")),paste0("GSE179042_", c("T","NT")))[3:4]
StudySamples <- c("GSE167075_Caco2shCTL","GSE167075_Caco2shM3")
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio for Caco2.M3KD RIPvsInputRatio
m6A_samples <- c("GSE167075_Caco2shCTL","GSE167075_Caco2shM3")
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.neg.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"CTL"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.pos.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"CTL"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)*2)*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)*2)*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
# legendPlot.mRNA <- plotLegend(legend = c("NonTumor", "Tumor"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)*3)*track.height, width = 2.5, height = 1,
#                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)*2)*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#Add mouse CT26 m6A and mRNA under peturbation
dt.selected.gene <- dt.gene.mm39 %>% filter(gene_name=="Ptrh1")
window.size=2#100kb
pseudocount=1
track.height=0.5
pdf("FigS5_Mouse_CT26_Mettl3_Mettl14_pertubation_Ptrh1.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 8.1, height = 7.6, default.units = "cm",showGuides = F)
mm39.assembly <- assembly(BSgenome = BSgenome.Mmusculus.UCSC.mm39,Genome = "mm39", TxDb = "TxDb.Mmusculus.UCSC.mm39.knownGene",OrgDb = "org.Mm.eg.db")
region <- pgParams(chrom = as.character(dt.selected.gene$seqnames),chromstart = (floor(dt.selected.gene$start/1000)-window.size)*1000, chromend = (ceiling(dt.selected.gene$end/1000)+window.size)*1000,assembly = mm39.assembly)
#plot gene
genesPlot <- plotGenes(params = region,geneHighlights = data.frame( "gene" = dt.selected.gene$gene_name, "color" = c("#eb4601")),geneBackground = "grey",
                       x = 0.5, y = 0, height = 1.0,width = track.width,just = c("left", "top"), default.units = "cm",fontsize = 7)
annoGenomeLabel(plot = genesPlot, x = 0.5, y = 1,scale = "Kb", just = c("left", "top"), fontsize = 6,default.units = "cm")

#add CT26 M3KO M14KO MeRIP (Input and RIP is strand unspecific)
StudySamples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")
if(dt.selected.gene$strand=="-"){
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.bw"),params=region)
    # Input <- Input.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.neg %>% mutate(score=-score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-1*range_both[sample==StudySamples[s],range],-log2(pseudocount)), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.3)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}else{
  dt.density <- foreach(s = 1:length(StudySamples), .combine = 'rbind')%do%{
    Input.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_Input.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",StudySamples[s],"_RIP.bw"),params=region)
    # Input <- Input.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
    # RIP <- RIP.pos %>% mutate(score=score) %>% dplyr::select(seqnames,start,end,score)
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
                             params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                             just = c("left", "top"), default.units = "cm",linecolor = "#FF8C00",fill="#FF8C00",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    signal_Input <- plotSignal(data = dt.density %>% dplyr::filter(library=="Input" & sample==StudySamples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(-log2(pseudocount),range_both[sample==StudySamples[s],range]), x = 0.5, y = 1.5+(s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = "grey80",fill="grey80",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_both[sample==StudySamples[s],range]+log2(pseudocount)),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(s-0.6)*track.height,just = c("bottom","left"), default.units = "cm",draw=T)
    sample_label <- plotText(label=StudySamples[s], fontsize = 6, fontface="plain",y=1.5+(s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#track for bigwig files of RIPvsInput ratio for Caco2.M3KD RIPvsInputRatio (strand unspecific)
m6A_samples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"control"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(m6A_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",m6A_samples[s],"_RIPvsInputRatio.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=m6A_samples[s],Study=strsplit(m6A_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_ratio <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(m6A_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==m6A_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(1,range_ratio[sample==m6A_samples[s],range]), x = 0.5, y = 1.5+(length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(m6A_samples[s],"control"),"#add9ee" ,"#7b72b7"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [1 - ",(range_ratio[sample==m6A_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=m6A_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add track for gene expression (not strand specific)
mRNA_samples <- c("GSE142589_control","GSE142589_M3","GSE142589_M14")[-3]
if(dt.selected.gene$strand == "-"){
  dt.ratio <- foreach(s=1:length(mRNA_samples),.combine='rbind')%do%{
    ratio.neg <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",mRNA_samples[s],"_Input.bw"),params=region)
    ratio.neg %>% mutate(score=score,strand="-") %>% mutate(sample=mRNA_samples[s],Study=strsplit(mRNA_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(mRNA_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==mRNA_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==mRNA_samples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(mRNA_samples[s],"control"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("-  [0 - ",(range_expr[sample==mRNA_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-0.35)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=mRNA_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+ length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }

}else{
  dt.ratio <- foreach(s=1:length(mRNA_samples),.combine='rbind')%do%{
    ratio.pos <- plotgardener::readBigwig(file = paste0("/data2/CRC_m6A/BigWig/",mRNA_samples[s],"_Input.bw"),params=region)
    ratio.pos %>% mutate(score=score,strand="+") %>% mutate(sample=mRNA_samples[s],Study=strsplit(mRNA_samples[s],split="_") %>% sapply("[",1))
  }
  #plot
  range_expr <- dt.ratio %>% group_by(Study) %>% mutate(range=ceiling(max(score))) %>% as.data.table() %>%  distinct(sample,Study,range)
  for(s in 1:length(mRNA_samples)){
    signal_Tumor <- plotSignal(data = dt.ratio %>% dplyr::filter(sample==mRNA_samples[s]) %>% dplyr::select(seqnames,start,end,score),
                               params = region,range = c(0,range_expr[sample==mRNA_samples[s],range]), x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height, width = track.width, height = track.height,
                               just = c("left", "top"), default.units = "cm",linecolor = ifelse(str_detect(mRNA_samples[s],"control"),"#f1aea7" ,"#ee6f6f"),fill="transparent",baseline = FALSE, baseline.color = NULL, baseline.lwd = 0, draw = T)
    plotText(label = paste0("+  [0 - ",(range_expr[sample==mRNA_samples[s],range]),"]"), fontsize = 6, fontface = "plain",x = 0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+s-0.65)*track.height,just = c("bottom","left"),
             default.units = "cm",draw=T)
    sample_label <- plotText(label=mRNA_samples[s], fontsize = 6, fontface="plain",y=1.5+(length(StudySamples)+ length(m6A_samples)+s-0.5)*track.height,x=0.6+track.width,draw=T,just = "left", fontcolor = "black",default.units = "cm")
    sample_rect <- plotRect(x=0.5,y=1.5+(length(StudySamples)+ length(m6A_samples)+s-1)*track.height,width=track.width, height = track.height,linecolor = "black", lwd = 0.5,just = c("left","top"),default.units = "cm")
  }
}

#add legend for MeRIP
legendPlot.MeRIP <- plotLegend(legend = c("Input", "RIP"), fill = c("grey80", "#FF8C00"),border = F,x =0.5, y = 1.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 1.5, height = 1,
                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                               title = expression("log"["2"]~"(CPM)"),fondface="bold")
legendPlot.RIPvsInput <- plotLegend(legend = c("Control", "Peturbation"), fill = c("#add9ee" ,"#7b72b7"),border = F, x =0.5 + 3, y = 1.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 1.5, height = 1,
                                    just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "RIP/Input enrichment",fondface="bold")
legendPlot.mRNA <- plotLegend(legend = c("Control", "Peturbation"), fill = c("#f1aea7" ,"#ee6f6f"),border = F, x =0.5 , y = 2.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",title = "mRNA expression",fondface="bold")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 2.5+(length(StudySamples)+ length(m6A_samples)+ length(mRNA_samples))*track.height+1,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#g DEG profiling of above novel and pertubation validated DMEGs
#load all TumorDEG, StageDEG, SurvivalGene result
TumorDEG.StageDEG.SurvivalGene.all.list <- readRDS("dt.TumorDEG.StageDEG.SurvivalGene.all.list.RDS")
names(TumorDEG.StageDEG.SurvivalGene.all.list)
dt.TumorDEG.all <- TumorDEG.StageDEG.SurvivalGene.all.list$TumorDEG
#visualization of validation of DMEG as consistent TumorDEG
dt.TumorDEG.selected.strict.novel.DMEG.perturb.validated <- dt.selected.strict.novel.DMEG.perturb.validated %>% dplyr::filter(nPerturbScore>=3) %>% dplyr::select(Study,Symbol,DiffMethy,contains("nPerturb")) %>%
  left_join(x=.,y=dt.TumorDEG.all %>% dplyr::filter(Symbol %in% dt.selected.strict.novel.DMEG.perturb.validated[nPerturbScore>=3,Symbol]) %>%
              dplyr::filter(FDR<=0.05) %>% group_by(Symbol,Direction) %>% mutate(nCohort=n_distinct(StudyID[FDR<=0.05]), LFC.avg=mean(LFC[FDR<=0.05]), lgFDR.avg=mean(-log10(FDR[FDR<=0.05]))) %>%
              as.data.table() %>% distinct(Symbol,Direction,nCohort,LFC.avg,lgFDR.avg,StudyID,LFC,FDR), by=c("Symbol")) %>% dplyr::filter( (DiffMethy=="Hyper" & Direction=="Up") | (DiffMethy=="Hypo" & Direction=="Down") ) %>%
  dplyr::arrange(DiffMethy,desc(nCohort),desc(abs(LFC.avg))) %>% as.data.table() %>% mutate(lgFDR = -log10(FDR), Symbol=factor(Symbol,levels=unique(Symbol))) %>%
  mutate(DiffMethy=factor(DiffMethy,levels=c("Hyper","Hypo"), labels=c("HyperUpDMEG","HypoDownDMEG")))

dt.label <- dt.TumorDEG.selected.strict.novel.DMEG.perturb.validated %>% distinct(DiffMethy,Symbol,nCohort) %>% mutate(label=paste0(nCohort,"/",15)) %>%
  mutate(ypos=case_when(DiffMethy=="Hyper" ~ 2.5, DiffMethy=="Hypo" ~ 0))

p.TumorDEG.selected.strict.novel.DMEG.perturb.validated <- dt.TumorDEG.selected.strict.novel.DMEG.perturb.validated  %>%
  ggplot(data=.,aes(x=Symbol, y=LFC))+
  geom_boxplot(fill="transparent", color="grey5",width=0.6, alpha=0.7, outlier.shape = NA,linewidth=0.3)+
  geom_jitter(aes(color=lgFDR),width=0.3,height = 0,alpha=1,size=1.)+
  geom_text(aes(x=Symbol,y=ypos,label=label),data=dt.label,size=2)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_color_gradient(low = c("#E9967A", "#FF3030")[1], high = c("#E9967A", "#FF3030")[2], limits=quantile(dt.TumorDEG.selected.strict.novel.DMEG.perturb.validated$lgFDR,c(0,0.95)),
                       na.value = "#EE2C2C",name = "-log10 FDR")+
  labs(x=str_wrap("Validated Strict novel DMEGs (nPerturbScore>=3)",width=45), y=str_wrap("LFC of DEG in CRC tumor vs adjacent",width=30))+
  facet_wrap(~DiffMethy,nrow = 1,scales = "free")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    strip.background = element_rect(fill="white"),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.text.x  = element_text(face = "italic"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.TumorDEG.selected.strict.novel.DMEG.perturb.validated

##h Survival Association significance of above novel and validated DMEGs
dt.SurvivalGene.all <- TumorDEG.StageDEG.SurvivalGene.all.list$SurvivalGene
#visualization of validation of DMEG as consistent survival markers
dt.SurvivalGene.selected.strict.novel.DMEG.perturb.validated <- dt.selected.strict.novel.DMEG.perturb.validated %>% dplyr::filter(nPerturbScore>=3) %>% dplyr::select(Study,Symbol,DiffMethy,contains("nPerturb")) %>%
  left_join(x=.,y=dt.SurvivalGene.all %>% dplyr::filter(Symbol %in% dt.selected.strict.novel.DMEG.perturb.validated[nPerturbScore>=3,Symbol]) %>%
              dplyr::filter(Pvalue<=0.05) %>% mutate(SurvivalMetricClass=factor(SurvivalMetric, levels=c("OS","PFS","PFI","DFS","DSS"), labels=c("OS", rep("PFS",4)))) %>%
              group_by(Symbol,SurvivalMetricClass,Direction) %>% mutate(nCohort=n_distinct(StudyID), lgP.avg=mean(-log10(Pvalue))) %>%
              group_by(Symbol,SurvivalMetricClass,Direction,StudyID) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>%
              as.data.table() %>% distinct(Symbol,SurvivalMetricClass,Direction,StudyID,SurvivalMetric,Pvalue,nCohort,lgP.avg), by=c("Symbol")) %>% dplyr::filter( (DiffMethy=="Hyper" & Direction=="Risk") | (DiffMethy=="Hypo" & Direction=="Protect") ) %>%
  dplyr::arrange(DiffMethy,desc(nCohort)) %>% mutate(Symbol=factor(Symbol, levels=unique(Symbol)))
p.SurvivalGene.selected.strict.novel.DMEG.perturb.validated <- dt.SurvivalGene.selected.strict.novel.DMEG.perturb.validated %>% mutate(logP=-log10(Pvalue)) %>%
  ggplot(data=.,aes(x=Symbol,y=logP))+
  geom_bar(aes(color=StudyID, linetype=Direction),fill="transparent",width=0.5,stat = "identity",position = position_dodge2())+
  geom_hline(yintercept = -log10(0.05),linetype=3, linewidth=0.5,color="grey")+
  scale_color_nejm()+
  scale_y_continuous(expand = expansion(mult = c(0,0.03)))+
  scale_linetype_manual(values = c(1,2),breaks = c("Risk","Protect"))+
  labs(x="Validated strict novel DMEGs (nPerturbScore>=3)", y="-log10 log-rank test p-value")+
  facet_wrap(~SurvivalMetricClass,nrow = 1,scales = "free")+
  guides(color=guide_legend(nrow = 2),linetype=guide_legend(nrow = 2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    strip.background = element_rect(fill="white"),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.text.x  = element_text(face = "italic"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.SurvivalGene.selected.strict.novel.DMEG.perturb.validated

#i Tumor progression association significance of above novel and validated DMEGs
dt.StageDEG.all <- TumorDEG.StageDEG.SurvivalGene.all.list$StageDEG
#visualization of stageDMEGs: PTRH1, PAPSS2, and LGALS4
dt.StageDEG.selected.strict.novel.DMEG.perturb.validated <- dt.selected.strict.novel.DMEG.perturb.validated %>% dplyr::filter(nPerturbScore>=3) %>% dplyr::select(Study,Symbol,DiffMethy,contains("nPerturb")) %>%
  left_join(x=.,y=dt.StageDEG.all %>% dplyr::filter(Symbol %in% dt.selected.strict.novel.DMEG.perturb.validated[nPerturbScore>=3,Symbol]) %>%
              dplyr::filter(Pvalue<=0.05) %>% group_by(Symbol,Direction) %>% mutate(nCohort=n_distinct(StudyID), LFC.avg=mean(LFC), lgFDR.avg=mean(-log10(FDR))) %>%
              as.data.table() %>% distinct(Symbol,Direction,StudyID,LFC,Pvalue,nCohort,LFC.avg,lgFDR.avg), by=c("Symbol")) %>% dplyr::filter( (DiffMethy=="Hyper" & Direction=="Up") | (DiffMethy=="Hypo" & Direction=="Down") ) %>%
  dplyr::arrange(desc(nCohort)) %>% mutate(Symbol=factor(Symbol,levels=unique(Symbol))) %>% mutate(logP=-log10(Pvalue))
p.StageDEG.selected.strict.novel.DMEG.perturb.validated <- dt.StageDEG.selected.strict.novel.DMEG.perturb.validated %>%
  ggplot(data=.,aes(x=Symbol,y=LFC))+
  geom_segment(aes(x=Symbol,xend=Symbol,y=0,yend=LFC),size=0.7,linetype="dotdash",color="grey5")+
  geom_point(size=3, color="transparent", aes(fill=Pvalue), shape=21, stroke=2) +
  scale_fill_gradient(low = c("#E9967A", "#FF3030")[2], high = c("#E9967A", "#FF3030")[1], limits=quantile(dt.StageDEG.selected.strict.novel.DMEG.perturb.validated$Pvalue,c(0,0.95)),
                      na.value = "#EE2C2C",name = "StageDEG p-value")+
  labs(x=str_wrap("Validated Strict novel DMEGs (nPerturbScore>=3)",width = 40), y="LFC of StageDEG")+
  facet_wrap(~StudyID,scales = "free",nrow = 1)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    strip.background = element_rect(fill="white"),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.text.x  = element_text(face = "italic"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.StageDEG.selected.strict.novel.DMEG.perturb.validated

##j Venn plot highlighted PTRH1 and C5AR1 for further validation
#venn diagram of StageDEG, TumorDEG, and SurvivalGenes
TumorDEG.StageDEG.SurvivalGene.DMEG.list <- NULL
TumorDEG.StageDEG.SurvivalGene.DMEG.list$TumorDEG <- unique(dt.TumorDEG.selected.strict.novel.DMEG.perturb.validated[nCohort>=2,Symbol])
TumorDEG.StageDEG.SurvivalGene.DMEG.list$StageDEG <- unique(dt.StageDEG.selected.strict.novel.DMEG.perturb.validated[nCohort>=1,Symbol])
TumorDEG.StageDEG.SurvivalGene.DMEG.list$SurvivalGene <- unique(dt.SurvivalGene.selected.strict.novel.DMEG.perturb.validated[nCohort>=2,Symbol])
names(TumorDEG.StageDEG.SurvivalGene.DMEG.list) <- c(paste0("TumorDEG","\n","(nCohort>=2)"),"StageDEG", paste0("SurvivalGene", "\n", "(nCohort>=2)"))
p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene <- ggVennDiagram(TumorDEG.StageDEG.SurvivalGene.DMEG.list,label = c("count","percent","both","none")[1],
                                                              category.names=names(TumorDEG.StageDEG.SurvivalGene.DMEG.list),
                                                              set_color =c("#D15FEE", "#EE30A7", "#FF6A6A"),
                                                              edge_lty = c("solid","dashed","dotted","dotdash","longdash","twodash")[c(1)],
                                                              set_size = 2.5,label_size = 2.,label_alpha = 0.6,
                                                              edge_size = 0.5)+
  scale_fill_gradient(high = "white",low="white",guide = "none")+
  scale_x_continuous(expand = expansion(mult = c(0.1,0.1)))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  labs(title = "Functional relevant of novel and m6A perturbation validated DMEGs")+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    axis.line=element_blank(),axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    # axis.text.y = element_text(face = "plain"),
    # axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    # panel.grid.major.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
center_x <- mean(ggplot_build(p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene)$layout$panel_params[[1]]$x.range)
center_y <- quantile(ggplot_build(p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene)$layout$panel_params[[1]]$y.range,0.5)
p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene <- p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene+
  annotate("text", x=center_x, y=center_y,
           label = paste(Reduce(x=TumorDEG.StageDEG.SurvivalGene.DMEG.list,intersect),collapse = ", "),
           size = 2., vjust = 0.5,color="grey5",fontface="bold.italic",family="Helvetica")
p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene


##combine all subplot except IGV plots together
fig5s.list <- list(p.Perturbm6A.DMEG=p.Perturbm6A.DMEG,
                   p.PerturbRNA.DMEG=p.PerturbRNA.DMEG,
                   p.DMEG.perturb.validated=p.DMEG.perturb.validated,
                   p.venn.strict.novel.DMEG=p.venn.strict.novel.DMEG,
                   p.pie.selected.strict.novel.DMEG=p.pie.selected.strict.novel.DMEG,
                   p.scatterplot.strict.novel.DMEG.perturb.validated=p.scatterplot.strict.novel.DMEG.perturb.validated,
                   p.TumorDEG.selected.strict.novel.DMEG.perturb.validated=p.TumorDEG.selected.strict.novel.DMEG.perturb.validated,
                   p.StageDEG.selected.strict.novel.DMEG.perturb.validated=p.StageDEG.selected.strict.novel.DMEG.perturb.validated,
                   p.SurvivalGene.selected.strict.novel.DMEG.perturb.validated=p.SurvivalGene.selected.strict.novel.DMEG.perturb.validated,
                   p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene=p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene
                   )
saveRDS(fig5s.list, file="FigureS5.plot.list.RDS")
save.image("sm6APeak_FigureS5_intermediate.results.RDS")

pdf("Figure5S_In_silico_validation_of_DMEG_highlighted_PTRH1_and_C5AR1_in_CRC.pdf",width = 8.2677, height = 11.693)
pageCreate(width =18.1, height =24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.Perturbm6A.DMEG, x = 0.05, y=0.2, default.units = "cm",width = 6.5, height = 4.6)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05+5.8, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.PerturbRNA.DMEG, x = 0.05+5.5, y=0.2, default.units = "cm",width = 7.5, height = 4.6)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 12, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DMEG.perturb.validated, x = 12, y=0.2, default.units = "cm",width = 4.5, height = 4.6)
#row2
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05+4.8*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.venn.strict.novel.DMEG, x = 0.05, y=0.2+4.8*1, default.units = "cm",width = 6.5, height = 4.6)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 5.2, y = 0.05+4.8*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.pie.selected.strict.novel.DMEG, x = 5.2, y=0.2+4.8*1, default.units = "cm",width = 6.5, height = 4.6)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 10, y = 0.05+4.8*1,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x=10,y=0.2+4.8*1, width = 7.8, height = 9.2,just = c("top","left"), default.units = "cm",draw=T)
#row3
plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05+4.8*2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.scatterplot.strict.novel.DMEG.perturb.validated, x = 0.05, y=0.2+4.8*2, default.units = "cm",width = 4.7, height = 4.6)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 4.8, y = 0.05+4.8*2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.StageDEG.selected.strict.novel.DMEG.perturb.validated, x = 4.8, y=0.2+4.8*2, default.units = "cm",width = 5.2, height = 4.6)
plotText(label = "H", fontsize = 8, fontface = "bold",x = 10, y = 14.4,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x=10,y=14.6, width = 7.8, height = 9.2,just = c("top","left"), default.units = "cm",draw=T)
#row4
plotText(label = "I", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05+4.8*3,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TumorDEG.selected.strict.novel.DMEG.perturb.validated, x = 0.05, y=0.2+4.8*3, default.units = "cm",width = 5.8, height = 4.6)
plotText(label = "L", fontsize = 8, fontface = "bold",x = 6, y = 0.05+4.8*3,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.venn.PerturbDMEG.TumorDEG.StageDEG.SurvivalGene, x = 6, y=0.2+4.8*3, default.units = "cm",width = 4.5, height = 4.6)
#row5
plotText(label = "J", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05+4.8*4,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.SurvivalGene.selected.strict.novel.DMEG.perturb.validated, x = 0.05, y=0.2+4.8*4, default.units = "cm",width = 9.5, height = 4.6)


dev.off()
