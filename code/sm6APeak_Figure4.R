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

# nature_seq <- c("#440154","#482173","#433E85","#31688E","#21908C","#35B779","#8FD744")
# nejm_div  <- c("#8C1515","#D04A4A","#E58585","#F0F0F0","#85C3E5","#4A96D0","#2B5F8C")

### Part4 Integration via machine learning enable comprehensive m6A calling ###############
#a the jaccard based similarity between each method's best option
##################### load the six method's best option #####################
dt.pruned.FPR.TPR.option.AUC.NEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HEK.RDS")
dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% distinct(Method,ComboName,optionID,AUC.scaled)
dt.pruned.FPR.TPR.option.AUC.NEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.mESC.RDS")
dt.pruned.FPR.TPR.option.AUC.NEB.mESC %>% distinct(Method,ComboName,optionID,AUC.scaled)

dt.AUC.NEB.all.option <- rbind(dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>% mutate(Cell="HEK293"),
                               dt.pruned.FPR.TPR.option.AUC.NEB.mESC %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>% mutate(Cell="mESC"))

dt.BestComboOption.NEB <- dt.AUC.NEB.all.option %>% filter(Cell=="HEK293") %>% group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>%
  dplyr::arrange(desc(AUC.scaled),desc(TPR20))
#obtain the cutoff for each best option
dt.BestComboOption.NEB <- dt.BestComboOption.NEB %>% left_join(x=.,y=dt.pruned.FPR.TPR.option.AUC.NEB.HEK %>% dplyr::select(FPR,TPR,cutoff,Method,ComboName,optionID,Benchmark) %>% dplyr::filter(FPR==0.2),by=c("Method","ComboName","optionID"))
#obtain the m6A peak according to the cutoff for each best method's option
#HEK293
#run bedtools intensity
#prepare the stranded.BAM.Depth
Samples <- c("HEK_NEB_mRNA1","HEK_NEB_mRNA2")
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


dt.Method.m6a <- foreach(o = 1:nrow(dt.BestComboOption.NEB), .combine='rbind')%do%{
  #load the old method m6A
  dt.m6a <- LoadPeakMethod(Method=dt.BestComboOption.NEB$Method[o],
                                  Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB/",dt.BestComboOption.NEB$Method[o],"_",dt.BestComboOption.NEB$ComboName[o])) %>%
    mutate(Method=dt.BestComboOption.NEB$Method[o], ComboName=dt.BestComboOption.NEB$ComboName[o])
  dt.m6a %>% filter(str_detect(seqnames,"chr"))  %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
}
dt.Method.m6a %>% distinct(Method,ComboName,Sample)
  #prepare the optionID
  selected.option <- dt.BestComboOption.NEB$optionID[o]
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
        selected.option <- dt.BestComboOption.NEB[Method==M & ComboName==combo,optionID]
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
dt.method.combo.selected <- dt.BestComboOption.NEB
dt.BestOption.m6A.FPR20.HEK293 <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
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
dt.BestOption.m6A.FPR20.HEK293  <- dt.BestOption.m6A.FPR20.HEK293 %>% mutate(name=paste(Method,ComboName,optionID,Sample,name,sep="_"))


#mESC
#prepare the stranded.BAM.Depth
Samples <- c("mESC_WT1", "mESC_WT2")
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
dt.Method.m6a <- foreach(o = 1:nrow(dt.BestComboOption.NEB), .combine='rbind')%do%{
  #load the old method m6A
  dt.m6a <- LoadPeakMethod(Method=dt.BestComboOption.NEB$Method[o],
                           Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB_mESC/",dt.BestComboOption.NEB$Method[o],"_",dt.BestComboOption.NEB$ComboName[o])) %>%
    mutate(Method=dt.BestComboOption.NEB$Method[o], ComboName=dt.BestComboOption.NEB$ComboName[o])
  dt.m6a %>% filter(str_detect(seqnames,"chr"))  %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,Method,ComboName)
}
dt.Method.m6a %>% distinct(Method,ComboName,Sample)

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
        selected.option <- dt.BestComboOption.NEB[Method==M & ComboName==combo,optionID]
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
dt.method.combo.selected <- dt.BestComboOption.NEB
dt.BestOption.m6A.FPR20.mESC <- foreach(M=unique(dt.method.combo.selected$Method),.combine='rbind')%do%{
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
dt.BestOption.m6A.FPR20.mESC  <- dt.BestOption.m6A.FPR20.mESC %>% mutate(name=paste(Method,ComboName,optionID,Sample,name,sep="_"))
dt.BestOption.m6A.FPR20.mESC %>% distinct(Method,ComboName,optionID)

#obtain the overlapped benchmark sites for these BestOption.m6A.FPR20
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
dt.BestOption.HEK293.peaks.overlappedSites <- foreach(M = unique(dt.BestOption.m6A.FPR20.HEK293$Method), .combine='rbind')%do%{
  foreach(S = unique(dt.BestOption.m6A.FPR20.HEK293$Sample), .combine='rbind')%do%{
    bed.peak <- dt.BestOption.m6A.FPR20.HEK293 %>% dplyr::filter(Method==M & Sample==S) %>% dplyr::select(seqnames,start,end,name,score,strand)
    if(str_detect(c,"HEK")){
      bed.bench <- dt.benchmark.m6A.HEK293 %>% dplyr::filter(benchmark_name!="SACseq_HEK293") %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
    }
    dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
    dt.overlap <- dt.overlap %>% dplyr::select(site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
    dt.overlap %>% mutate(Method=M, Sample=S)
  }
}
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
dt.BestOption.mESC.peaks.overlappedSites <- foreach(M = unique(dt.BestOption.m6A.FPR20.mESC$Method), .combine='rbind')%do%{
  foreach(S = unique(dt.BestOption.m6A.FPR20.mESC$Sample), .combine='rbind')%do%{
    bed.peak <- dt.BestOption.m6A.FPR20.mESC %>% dplyr::filter(Method==M & Sample==S) %>% dplyr::select(seqnames,start,end,name,score,strand)
    if(str_detect(c,"HEK")){
      bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
    }
    dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
    dt.overlap <- dt.overlap %>% dplyr::select(site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
    dt.overlap %>% mutate(Method=M, Sample=S)
  }
}

dt.BestOption.mESC.peaks.overlappedSites[,.(Nsite=n_distinct(site)),by=.(Method,benchmark_name,Sample)] %>% dplyr::arrange(Method,benchmark_name,Sample,Nsite)

#combine overlapped Sites together
dt.BestOption.peaks.overlappedSites  <- rbind(dt.BestOption.HEK293.peaks.overlappedSites,dt.BestOption.mESC.peaks.overlappedSites)
#calculate the jaccard similarity between Sites of each method's best option
dt.jaccard.BestOption.overlappedBenchmarkSite <- foreach(S=unique(dt.BestOption.peaks.overlappedSites$Sample), .combine='rbind')%do%{
  dt.jaccard.benchmark <- foreach(b=unique(dt.BestOption.peaks.overlappedSites[Sample==S,benchmark_name]), .combine='rbind')%do%{
    foreach(M1 = unique(dt.BestOption.peaks.overlappedSites$Method), .combine='rbind')%do%{
      site.list1 <- dt.BestOption.peaks.overlappedSites[benchmark_name==b & Sample==S & Method==M1,site] %>% unique()
      foreach(M2=unique(dt.BestOption.peaks.overlappedSites[Method!=M1,Method]), .combine='rbind')%do%{
        site.list2  <- dt.BestOption.peaks.overlappedSites[benchmark_name==b & Sample==S & Method==M2,site] %>% unique()
        nShare <- intersect(site.list1,site.list2) %>% unique() %>% length()
        nUnion <- c(site.list1,site.list2) %>% unique() %>% length()
        data.table(Sample=S,Benchmark=b, Method1=M1, Method2=M2, nShare=nShare, Jaccard=nShare/nUnion, RatioImproveMethod1=nUnion/length(site.list1), RatioImproveMethod2=nUnion/length(site.list2))
      }
    }
  }
  dt.jaccard.benchmark
}
dt.jaccard.BestOption.overlappedBenchmarkSite %>%
  filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
  group_by(Benchmark,Method1,Method2) %>% mutate(avg.Jaccard=mean(Jaccard),avg.RatioImproveMethod1=mean(RatioImproveMethod1), avg.RatioImproveMethod2=mean(RatioImproveMethod2)) %>%
  as.data.table() %>% distinct(Benchmark,Method1,Method2,avg.Jaccard,avg.RatioImproveMethod1,avg.RatioImproveMethod2) %>%
  dplyr::arrange(Benchmark,Method1,Method2)

#assess the overlap similarity between FPR20 peaks
dt.jaccard.BestOption.FPR20.peaks.HEK293 <- foreach(S=unique(dt.BestOption.m6A.FPR20.HEK293$Sample), .combine='rbind')%do%{
  foreach(M1 = unique(dt.BestOption.m6A.FPR20.HEK293$Method), .combine='rbind')%do%{
    bed.peak1 <- dt.BestOption.m6A.FPR20.HEK293 %>% dplyr::filter(Sample==S & Method==M1) %>% dplyr::select(seqnames,start,end,name,score,strand)
    foreach(M2 = unique(dt.BestOption.peaks.overlappedSites[Method!=M1,Method]), .combine='rbind')%do%{
      bed.peak2 <- dt.BestOption.m6A.FPR20.HEK293 %>% dplyr::filter(Sample==S & Method==M2) %>% dplyr::select(seqnames,start,end,name,score,strand)
      dt.overlap <- bt.intersect(a=bed.peak2, b=bed.peak1, s=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
      dt.overlap <- dt.overlap %>% dplyr::select(peak2=V4,peak1=V10) %>% distinct()
      nNovel1=bed.peak1 %>% filter(!name %in% dt.overlap$peak1) %>% nrow()# novel peak in peak1
      nNovel2=bed.peak2 %>% filter(!name %in% dt.overlap[peak1 != ".",peak2]) %>% nrow()# novel peak in peak2
      data.table(Sample=S,Method1=M1, Method2=M2, nPeak1=nrow(bed.peak1), nPeak2=nrow(bed.peak2),nNovelPeak1=nNovel1, nNovelPeak2=nNovel2) %>%
        mutate(NoveltyRatio1=nNovelPeak1/nPeak1, NoveltyRatio2=nNovelPeak2/nPeak2) %>%
        mutate(NoveltyPeakRatio=(NoveltyRatio1+NoveltyRatio2)/2)
    }
  }
}
dt.jaccard.BestOption.FPR20.peaks.HEK293 %>%
  filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
  group_by(Method1,Method2) %>% mutate(avg.NoveltyPeakRatio=mean(NoveltyRatio1)) %>%
  as.data.table() %>% distinct(Method1,Method2,avg.NoveltyPeakRatio) %>%
  dplyr::arrange(Method1,Method2,avg.NoveltyPeakRatio)

dt.jaccard.BestOption.FPR20.peaks.mESC <- foreach(S=unique(dt.BestOption.m6A.FPR20.mESC$Sample), .combine='rbind')%do%{
  foreach(M1 = unique(dt.BestOption.m6A.FPR20.mESC$Method), .combine='rbind')%do%{
    bed.peak1 <- dt.BestOption.m6A.FPR20.mESC %>% dplyr::filter(Sample==S & Method==M1) %>% dplyr::select(seqnames,start,end,name,score,strand)
    foreach(M2 = unique(dt.BestOption.peaks.overlappedSites[Method!=M1,Method]), .combine='rbind')%do%{
      bed.peak2 <- dt.BestOption.m6A.FPR20.mESC %>% dplyr::filter(Sample==S & Method==M2) %>% dplyr::select(seqnames,start,end,name,score,strand)
      dt.overlap <- bt.intersect(a=bed.peak2, b=bed.peak1, s=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
      dt.overlap <- dt.overlap %>% dplyr::select(peak2=V4,peak1=V10) %>% distinct()
      nNovel1=bed.peak1 %>% filter(!name %in% dt.overlap$peak1) %>% nrow()# novel peak in peak1
      nNovel2=bed.peak2 %>% filter(!name %in% dt.overlap[peak1 != ".",peak2]) %>% nrow()# novel peak in peak2
      data.table(Sample=S,Method1=M1, Method2=M2, nPeak1=nrow(bed.peak1), nPeak2=nrow(bed.peak2),nNovelPeak1=nNovel1, nNovelPeak2=nNovel2) %>%
        mutate(NoveltyRatio1=nNovelPeak1/nPeak1, NoveltyRatio2=nNovelPeak2/nPeak2) %>%
        mutate(NoveltyPeakRatio=(NoveltyRatio1+NoveltyRatio2)/2)
    }
  }
}
dt.jaccard.BestOption.FPR20.peaks.mESC %>%
  filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
  group_by(Method1,Method2) %>% mutate(avg.NoveltyPeakRatio=mean(NoveltyRatio1)) %>%
  as.data.table() %>% distinct(Method1,Method2,avg.NoveltyPeakRatio) %>%
  dplyr::arrange(Method1,Method2,avg.NoveltyPeakRatio)

#heatmap to visualize the novelty ratio of top 4 BestOption FPR20 peaks
dt.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC <- rbind(dt.jaccard.BestOption.FPR20.peaks.mESC %>%
                                                                   filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
                                                                   group_by(Method1,Method2) %>% mutate(avg.NoveltyPeakRatio=mean(NoveltyRatio1)) %>%
                                                                   as.data.table() %>% distinct(Method1,Method2,avg.NoveltyPeakRatio) %>%
                                                                   dplyr::arrange(Method1,Method2,avg.NoveltyPeakRatio) %>% mutate(Cell="mESC"),
                                                                 dt.jaccard.BestOption.FPR20.peaks.HEK293 %>%
                                                                   filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
                                                                   group_by(Method1,Method2) %>% mutate(avg.NoveltyPeakRatio=mean(NoveltyRatio1)) %>%
                                                                   as.data.table() %>% distinct(Method1,Method2,avg.NoveltyPeakRatio) %>%
                                                                   dplyr::arrange(Method1,Method2,avg.NoveltyPeakRatio) %>% mutate(Cell="HEK293")) %>%
  mutate(Cell=factor(Cell, levels=c("HEK293","mESC"), labels=c("NEB_HEK293","NEB_mESC")), Method1=factor(Method1,levels=c("MACS2","TRESS","exomePeak2","exomePeak")),
         Method2=factor(Method2,levels=c("MACS2","TRESS","exomePeak2","exomePeak"))) %>% mutate(avg.NoveltyPeakRatio=avg.NoveltyPeakRatio*100)
# compute limits and midpoint if needed
zlim <- range(dt.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC$avg.NoveltyPeakRatio, na.rm = TRUE)
zmid <- median(dt.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC$avg.NoveltyPeakRatio, na.rm = TRUE)
p.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC <- ggplot(data=dt.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC,aes(x=Method1,y=Method2,fill=avg.NoveltyPeakRatio))+
  geom_tile(color = "white", size = 0.4) +
  # numeric labels inside tiles with 2 decimal places
  geom_text(aes(label = sprintf("%.2f", avg.NoveltyPeakRatio)), color = "black", size = 1.8) +
  coord_fixed() +# keep tiles square
  scale_fill_gradient2(
      low = "#fbeddf",  mid = "#edb5c1", high = "#ce3341",
      midpoint = zmid,
      name = "% of NovelPeak",
      limits = zlim,
      labels = scales::label_number(accuracy = 0.01),
      na.value = "grey90"
    )+
  facet_wrap(~Cell, nrow = 1)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.title = element_blank(),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    panel.grid = element_blank(),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(15,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )

#visualization of jaccard similarity of those overlapped Benchmark Sites
dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC <- dt.jaccard.BestOption.overlappedBenchmarkSite %>%
  filter(Method1 %in% c("MACS2","exomePeak2","exomePeak","TRESS") & Method2 %in% c("MACS2","exomePeak2","exomePeak","TRESS")) %>%
  group_by(Benchmark,Method1,Method2) %>% mutate(avg.Jaccard=mean(Jaccard),avg.RatioImproveMethod=mean(RatioImproveMethod1)) %>%
  as.data.table() %>% distinct(Benchmark,Method1,Method2,avg.Jaccard,avg.RatioImproveMethod) %>%
  dplyr::arrange(Benchmark,Method1,Method2)
dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC <- dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC %>%
  mutate(Benchmark=factor(Benchmark, levels=c("GLORI_HEK293","GLORI_mESC","eTAMseq_mESC"), labels=c("HEK293_GLORI","mESC_GLORI","mESC_eTAMseq")),
         Method1=factor(Method1,levels=c("MACS2","TRESS","exomePeak2","exomePeak")),
         Method2=factor(Method2,levels=c("MACS2","TRESS","exomePeak2","exomePeak"))) %>% mutate(avg.RatioImproveMethod=(avg.RatioImproveMethod-1)*100)
# compute limits and midpoint if needed
zlim <- range(dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC$avg.RatioImproveMethod, na.rm = TRUE)
zmid <- median(dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC$avg.RatioImproveMethod, na.rm = TRUE)
p.heatmap.TPRImprovement.BestOption.FPR20.peak.HEK293.mESC <- ggplot(data=dt.ImprovedBenchmarkSite.BestOption.FPR20.peak.HEK293.mESC,aes(x=Method1,y=Method2,fill=avg.RatioImproveMethod))+
  geom_tile(color = "white", size = 0.4) +
  # numeric labels inside tiles with 2 decimal places
  geom_text(aes(label = sprintf("%.2f", avg.RatioImproveMethod)), color = "black", size = 1.8) +
  coord_fixed() +# keep tiles square
  scale_fill_gradient2(
    low = "#fbeddf",  mid = "#edb5c1", high = "#ce3341",
    midpoint = zmid,
    name = "% of TPR Improvement",
    limits = zlim,
    labels = scales::label_number(accuracy = 0.01),
    na.value = "grey90"
  )+
  facet_wrap(~Benchmark, nrow = 1)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.title = element_blank(),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    panel.grid = element_blank(),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(15,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )

#b the MCA -based dimension reduction and integration
#load MCA result of HEK_NEB
MCA.files <- list.files(pattern = "^res\\.mca.*\\.RDS$") %>% grep(pattern="NEB",value=T)
names(MCA.files) <- MCA.files %>% strsplit(split=".",fixed=T) %>% sapply("[",3)
dt.MCA.top.dimension <- foreach(i = 1:length(MCA.files), .combine='rbind')%do%{
  mca.res <- readRDS(MCA.files[i])
  dt.percent.of.var <- mca.res$eig %>% as.data.table(keep.rownames = "Dim") %>% dplyr::select(Dim,`percentage of variance`, `cumulative percentage of variance`) %>%
    dplyr::filter(Dim %in% paste0("dim ",1:20)) %>% mutate(Sample=names(MCA.files)[i])
  dt.percent.of.var
}
dt.MCA.top.dim.BestOption.FPR20.BenchmarkSites <- dt.MCA.top.dimension %>% group_by(Dim) %>% mutate(POV = mean(`percentage of variance`), CumulativePOV=mean(`cumulative percentage of variance`)) %>%
  as.data.table() %>% dplyr::distinct(Dim,POV,CumulativePOV) %>% mutate(Dimension=str_replace(Dim,pattern="dim ",replacement="") %>% as.integer())

p.MCA.top.dimension.NEB <- ggplot(data=dt.MCA.top.dim.BestOption.FPR20.BenchmarkSites, aes(x=Dimension, y =POV))+
  geom_bar(width = 0.7,fill="#4682b4",alpha=1,stat = "identity")+
  geom_line(aes(y=CumulativePOV),linecolor="grey5",linetype=1)+
  geom_point(aes(y=CumulativePOV),size=1,shape=22)+
  coord_cartesian(xlim = c(0,15))+
  scale_y_continuous(breaks = seq(0,80,20))+
  labs(x="Dimension of multiple correspondence analysis", y="% of explained variance", title=str_wrap("MCA dimension reduction of benchmark sites similarity",width=60))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    plot.title = element_blank(),
    text = element_text(family = "Helvetica", colour = "black"),
    # axis.title = element_blank(),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    # axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    panel.grid = element_blank(),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(15,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )

dt.avg.TPR.FPR.MCA.topn.combo <- readRDS(file="dt.average.TPR.FPR.MCA.topn.combo.HEK.NEB.RDS")
dt.avg.TPR.FPR.MCA.topn.combo <- dt.avg.TPR.FPR.MCA.topn.combo %>% mutate(TopNDim=strsplit(Strategy,split=".",fixed=T) %>% sapply("[",2)) %>% group_by(TopNDim) %>%
  mutate(FPR=mean(avg.FPR), TPR=mean(avg.TPR)) %>% as.data.table() %>% distinct(TopNDim,FPR,TPR) %>% mutate(TopNDim=factor(TopNDim,levels=paste0("top",seq(2,10,by=2),"Dim"))) %>%
  pivot_longer(cols = c("FPR","TPR"),names_to = "Metric",values_to = "Value") %>% as.data.table()

p.FPR.TPR.top.dimension.NEB <- ggplot(data=dt.avg.TPR.FPR.MCA.topn.combo,aes(x=TopNDim,y=Value,fill=Metric))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.avg.TPR.FPR.MCA.topn.combo %>% dplyr::filter(Metric=="FPR"),
            aes(label = round(Value,2), y = Value+0.02, x=TopNDim),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.avg.TPR.FPR.MCA.topn.combo %>% dplyr::filter(Metric=="TPR"),
           aes(label = round(Value,2), y = Value+0.02, x=TopNDim),
           color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
  scale_fill_manual(values=c("#86a959","#a03d48"),breaks=c("FPR","TPR"))+
  labs(x=str_wrap("Integration of ComboOption most correlated with topN dimension",width=50),y=str_wrap("Value of surrogate FPR/TPR",width=40))+
  guides(fill=guide_legend(title = "Metric"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.FPR.TPR.top.dimension.NEB

#c the better performance AUC of sm6APeak over each method's best option
#load TPR FPR of NEB_HEK and NEB_mESC of sm6APeak
dt.FPR.TPR.FourMode.sm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
dt.FPP.TPR.NEB.sm6APeak.HEK293.mESC <- dt.FPR.TPR.FourMode.sm6APeak %>% dplyr::filter(SelectedMode=="NEB" & str_detect(SelectedSample,"NEB")) %>%
  mutate(Cell=factor(SelectedSample,levels=c("NEB","NEB_mESC"), labels=c("HEK293","mESC")))
#load the Method's  best option of NEB
dt.pruned.TPR.FPR.NEB.Method.BestOption <- readRDS("dt.pruned.TPR.FPR.NEB.Method.BestOption.RDS")
#combine sm6APeak with those top4 best combo
dt.TPR.FPR.sm6APeak.Method.BestOption <- rbind(dt.FPP.TPR.NEB.sm6APeak.HEK293.mESC %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell) %>% mutate(Method="sm6APeak"),
                                               dt.pruned.TPR.FPR.NEB.Method.BestOption %>% dplyr::filter(!Method %in% c("MeTPeak","MeRIPtools")) %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell,Method) ) %>%
  mutate(Cell=factor(Cell,levels=c("HEK293","mESC"), labels=c("NEB_HEK293","NEB_mESC")))

dt.TPR.FPR.sm6APeak.Method.BestOption.label <- dt.TPR.FPR.sm6APeak.Method.BestOption %>% distinct(Method,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>%
  mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/5),5)*max(dt.TPR.FPR.sm6APeak.Method.BestOption$TPR),2)) %>% mutate(xpos=0.6)
p.TPR.FPR.sm6APeak.Method.BestOption  <-  ggplot(data=dt.TPR.FPR.sm6APeak.Method.BestOption ,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of sm6APeak or method's BestCombo",y=str_wrap("Surrogate TPR of sm6APeak or method's BestCombo",width = 40))+
  geom_text(data=dt.TPR.FPR.sm6APeak.Method.BestOption.label, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#eb4601",c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0")[c(1:3,5)]),breaks = c("sm6APeak",c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools")[c(1:3,5)]),guide="none")+
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
p.TPR.FPR.sm6APeak.Method.BestOption

#d the better performance of sm6APeak is highly consistent on different anti-m6A
dt.FPR.TPR.FourMode.ISm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
#heatmap display the scaled AUC when apply ISm6APeak from different model display highly consistent performance
dt.AUC.ISm6APeak.FourModel <- dt.FPR.TPR.FourMode.ISm6APeak %>% dplyr::select(Model=SelectedMode, Sample=SelectedSample,AUC.scaled) %>% distinct() %>%
  mutate(Model=factor(Model,levels=c("NEB","SYSY","Abcam","Millipore","NEB_HeLa"), labels=paste0("ISm6APeak_",c("NEB","SYSY","Abcam","Millipore","NEB_HeLa"))),
         Sample=factor(Sample, levels=c("NEB","NEB_mESC","SYSY","Abcam","Millipore","NEB_HeLa"), labels=c("NEB_HEK293","NEB_mESC","SYSY_HEK293","Abcam_HEK293","Millipore_mESC","NEB_HeLa")))

# compute limits and midpoint if needed
zlim <- range(dt.AUC.ISm6APeak.FourModel$AUC.scaled, na.rm = TRUE)
# zmid <- median(dt.AUC.ISm6APeak.FourModel$AUC.scaled, na.rm = TRUE)
zmid <- (zlim[1]+zlim[2])/2
p.heatmap.AUC.ISm6APeak.FourModel <- ggplot(data=dt.AUC.ISm6APeak.FourModel,aes(y=Model,x=Sample,fill=AUC.scaled))+
  geom_tile(color = "white", size = 0.4) +
  # numeric labels inside tiles with 2 decimal places
  geom_text(aes(label = sprintf("%.2f", AUC.scaled)), color = "black", size = 1.8) +
  coord_fixed() +# keep tiles square
  scale_fill_gradient2(
    low = "#fbeddf",  mid = "#edb5c1", high = "#ce3341",
    midpoint = zmid,
    name = "Scaled AUC",
    limits = zlim,
    labels = scales::label_number(accuracy = 0.01),
    na.value = "grey90"
  )+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    axis.title = element_blank(),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    axis.text.x = element_text(face = "plain", angle = 30,vjust = 1,hjust = 1),
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey92", size = 0.35),
    panel.grid = element_blank(),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(15,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )
p.heatmap.AUC.ISm6APeak.FourModel

pdf("new_Fig4F_ISm6APeak_of_all_model_in_all_samples.PDF",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot =p.heatmap.AUC.ISm6APeak.FourModel, x =12.1, y=5, default.units = "cm",width = 5.6, height = 4.7)
dev.off()

#determine the optimal cutoff for sm6APeak_NEB model
dt.cutoff.parameter.sm6APeak  <- fread("/data/m6A_calling_strategy/RIPPeakS_resource/dt.cutoff.at.all.FPR.M6APeakS.tsv")
dt.TPR.FPR.sm6APeak.NEB <- dt.cutoff.parameter.sm6APeak %>% distinct(M6APeakModel, MergedFPR,FPR,  TPR,TPRmFPR, maxTPRmFPR, Mode) %>% filter(M6APeakModel=="NEB") %>%
  mutate(cutoff=MergedFPR) %>% mutate(near.maxTPRmFPR=maxTPRmFPR-TPRmFPR<0.05)
dt.TPR.FPR.tolabel <- dt.TPR.FPR.sm6APeak.NEB %>% filter(Mode %in% c("Strict","Moderate","Sensitive")) %>% distinct(cutoff,FPR,TPR,Mode)

p.dumbbell.optimal.cutoff.sm6APeak_NEB <- dt.TPR.FPR.sm6APeak.NEB %>%
  ggplot(data=.)+
  geom_segment(aes(x = cutoff, y=FPR, yend=TPR, color=near.maxTPRmFPR))+
  geom_point(aes(x=cutoff,y=FPR),shape=16,size=1,color="grey70")+
  geom_point(aes(x=cutoff,y=TPR),shape=21,size=1,color="grey70")+
  geom_step(aes(x=cutoff,y=TPRmFPR),linewidth=0.5)+
  geom_label(data=dt.TPR.FPR.tolabel[Mode=="Sensitive",], aes(x=cutoff, y=TPR, label=Mode),color=c("#2f9842","#df7930","#266b9f")[1],size=1.8,nudge_y = 0.05,nudge_x = 0,angle=90)+
  geom_label(data=dt.TPR.FPR.tolabel[Mode=="Moderate",], aes(x=cutoff, y=TPR, label=Mode),color=c("#2f9842","#df7930","#266b9f")[2],size=1.8,nudge_y = 0.05,nudge_x = 0,angle=90)+
  geom_label(data=dt.TPR.FPR.tolabel[Mode=="Strict",], aes(x=cutoff, y=TPR, label=Mode),color=c("#2f9842","#df7930","#266b9f")[3],size=1.8,nudge_y = 0.05,nudge_x = 0,angle=90)+
  scale_color_manual(breaks=c(FALSE,TRUE),values=c("grey70", "#e31a1c"),guide="none")+
  coord_cartesian(xlim = c(0,0.6))+
  labs(x="Cutoff of sm6APeak", y="Surrogate TPR (up) and FPR (down)")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    # axis.title = element_blank(),
    axis.text = element_text(size = rel(0.85), colour = "black"),
    # axis.text.x = element_text(face = "plain", angle = 15,vjust = 1,hjust = 1),
    # panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92", size = 0.35),
    # panel.grid = element_blank(),
    legend.position = "top",  # small inset legend (adjust as desired)
    legend.background = element_blank(),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.height = unit(7,"pt"),legend.key.width = unit(15,"pt"),
    legend.title = element_text(face = "plain"),
    plot.caption = element_text(size = rel(0.75), hjust = 0)
  )

#e sm6APeak identify novel m6A in mESC and HEK293 compare with MACS2_BestComboOption (keep this in supplementary figures)

#f sm6APeak reveal novel m6A modification in blastocyst than reported in published study (GSE192440), could be validated using GLORI and eTAMseq
#load published m6A peaks
GSE192440.bed.files <- list.files("/data2/mESC_m6A/Analysis/",pattern = "(.bed)$",full.names = T)
names(GSE192440.bed.files) <- list.files("/data2/mESC_m6A/Analysis/",pattern = "(.bed)$",full.names = F) %>% strsplit(split="_MS_",fixed=T) %>%
  sapply("[",2) %>% str_replace(pattern=".bed",replacement = "") %>% str_replace(pattern="_IP",replacement = "")
GSE192440.bed.files <- GSE192440.bed.files[grep(names(GSE192440.bed.files),pattern="Blastocyst")]
dt.GSE192440.m6a <- foreach(i=1:length(GSE192440.bed.files), .combine = 'rbind')%do%{
  dt.m6a <- fread(GSE192440.bed.files[i]) %>% dplyr::select(seqnames=V1,start=V2,end=V3,score=V7) %>% distinct() %>% group_by(seqnames,start,end) %>% mutate(score=mean(score)) %>% as.data.table() %>% distinct()
  dt.m6a <- dt.m6a %>% mutate(strand="*",Sample=names(GSE192440.bed.files)[i]) %>% dplyr::filter(str_detect(seqnames,"chr"))
  dt.m6a
}
#lifover mm10 peaks to mm39 peaks
m6A_liftover <- function(
    chainfile="~/genome_db/LiftOverChain/mm10ToMm39.over.chain",
    input.bed=dt.GSE192440.m6a[Sample=="2cell_rep1",]%>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(seqnames,start,end,name,score,strand)
){
  chain <- rtracklayer::import.chain(chainfile)
  input.GR <- GenomicRanges::makeGRangesFromDataFrame(input.bed,keep.extra.columns = T)
  as.data.table(unlist(rtracklayer::liftOver(input.GR, chain))) %>% dplyr::select(seqnames,start,end,name,score,strand)
}
dt.GSE192440.m6a.liftover <- m6A_liftover(chainfile="~/genome_db/LiftOverChain/mm10ToMm39.over.chain",
                                          input.bed=dt.GSE192440.m6a %>%
                                            mutate(name=paste(seqnames,start,end,sep="_")) %>% distinct(seqnames,start,end,name,score,strand))
dt.GSE192440.m6a <- dt.GSE192440.m6a %>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(name,score,Sample) %>%
  inner_join(x=.,y=dt.GSE192440.m6a.liftover %>% dplyr::select(seqnames,start,end,name,strand) %>% distinct(), by=c("name")) %>%
  mutate(mm10.name=name) %>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,everything())

#load M6APeakS for GSE192440
# sm6APeak.peak.files <- list.files("/data2/mESC_m6A/M6APeakS/",pattern="(GSE192440).+(Millipore_default.tsv)$",full.names = T)
# names(sm6APeak.peak.files) <- list.files("/data2/mESC_m6A/M6APeakS/",pattern="(GSE192440).+(Millipore_default.tsv)$",full.names = F) %>%
#   strsplit(split="GSE192440_") %>% sapply("[",2) %>% strsplit(split="_M6APeakS_peak") %>% sapply("[",1)
sm6APeak.peak.files <- list.files("/data3/mESC_m6A/sm6APeak_GSE192440_Bc/",pattern="(GSE192440).+(NEB_0.03.tsv)$",full.names = T)
names(sm6APeak.peak.files) <- list.files("/data3/mESC_m6A/sm6APeak_GSE192440_Bc/",pattern="(GSE192440).+(NEB_0.03.tsv)$",full.names = F) %>%
  strsplit(split="GSE192440_") %>% sapply("[",2) %>% strsplit(split="_M6APeakS_peak") %>% sapply("[",1)
sm6APeak.peak.files <- sm6APeak.peak.files[grep(names(sm6APeak.peak.files),pattern="Bc")]
dt.sm6APeak.m6a <- foreach(i = 1:length(sm6APeak.peak.files), .combine='rbind')%do%{
  fread(sm6APeak.peak.files[i]) %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,OverlappedGenes) %>% distinct() %>% mutate(Sample=names(sm6APeak.peak.files)[i])
}
dt.sm6APeak.m6a <- dt.sm6APeak.m6a %>%  mutate(Sample=factor(Sample,levels=c("Bc_rep1","Bc_rep2"),labels=c("Blastocyst_rep1","Blastocyst_rep2")))
dt.sm6APeak.m6a[,.N,by=.(Sample)]

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

dt.sm6APeak.m6a.Benchmarkoverlap <- foreach(c = unique(dt.sm6APeak.m6a$Sample), .combine='rbind')%do%{
  bed.peak <- dt.sm6APeak.m6a %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
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
dt.sm6APeak.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident)

dt.GSE192440.m6a.Benchmarkoverlap <- foreach(c = unique(dt.GSE192440.m6a$Sample), .combine='rbind')%do%{
  bed.peak <- dt.GSE192440.m6a %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
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
dt.GSE192440.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident)
dt.TPR.FPR.GSE192440.sm6APeak <- rbind(dt.GSE192440.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                        nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                         as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                       left_join(x=.,y=dt.GSE192440.m6a[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                         left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                         mutate(Group="GSE192440"),
                                       dt.sm6APeak.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                       nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                         as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                         left_join(x=.,y=dt.sm6APeak.m6a[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                         left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                         mutate(Group="sm6APeak")
                                         )
dt.TPR.FPR.GSE192440.sm6APeak %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)
# Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite     Group       TPR       FPR
# 1: Blastocyst_rep1     GLORI_mESC   27462 6076          2096      9135              96094 GSE192440 0.2857827 0.3348659
# 2: Blastocyst_rep1     GLORI_mESC   54096 6825          2655     10110              96094  sm6APeak 0.5629488 0.3249258
# 3: Blastocyst_rep1   eTAMseq_mESC   18973 4845          2162      9135              61365 GSE192440 0.3091828 0.4696223
# 4: Blastocyst_rep1   eTAMseq_mESC   36329 5938          2725     10110              61365  sm6APeak 0.5920150 0.4126607
# 5: Blastocyst_rep2     GLORI_mESC   17176 4313          1430      7948              96094 GSE192440 0.1787416 0.4573478
# 6: Blastocyst_rep2     GLORI_mESC   47784 8025          2658     20601              96094  sm6APeak 0.4972631 0.6104558
# 7: Blastocyst_rep2   eTAMseq_mESC   11931 3375          1482      7948              61365 GSE192440 0.1944268 0.5753649
# 8: Blastocyst_rep2   eTAMseq_mESC   33112 6639          2735     20601              61365  sm6APeak 0.5395910 0.6777341

#identify the replicate shared peaks for sm6APeak and GSE192440
dt.GSE192440.shared.m6A <- bt.intersect(a=dt.GSE192440.m6a %>% filter(Sample=="Blastocyst_rep1") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        b=dt.GSE192440.m6a %>% filter(Sample=="Blastocyst_rep2") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        s=F, wo=T) %>% filter(V15>=50) %>%
                           dplyr::mutate(peak1=V4,peak2=V11,seqnames=V1,start1=as.integer(V2),start2=as.integer(V9),end1=as.integer(V3),end2=as.integer(V10)) %>%
                          mutate(start=case_when(start1>=start2 ~ start1, .default = start2),end=case_when(end1<=end2 ~ end1, .default = end2)) %>%
                          dplyr::select(seqnames,start,end) %>% as.data.table() %>% mutate(strand="*",score=1000,name=paste(seqnames,start,end,sep="_")) %>%
                         dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample="shared")
dt.GSE192440.shared.m6a.Benchmarkoverlap <- foreach(c = "shared", .combine='rbind')%do%{
  bed.peak <- dt.GSE192440.shared.m6A  %>% dplyr::select(seqnames,start,end,name,score,strand)
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
dt.GSE192440.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                 nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
  left_join(x=.,y=dt.GSE192440.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
  left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
  mutate(Group="GSE192440.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)
#    Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite            Group       TPR       FPR
# 1: shared     GLORI_mESC   14471 3591          1263      5173              96094 GSE192440.shared 0.1505921 0.3058187
# 2: shared   eTAMseq_mESC   10168 2856          1312      5173              61365 GSE192440.shared 0.1656971 0.4479026
#sm6APeak
dt.sm6APeak.shared.m6A <- bt.intersect(a=dt.sm6APeak.m6a %>% filter(Sample=="Blastocyst_rep1") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        b=dt.sm6APeak.m6a %>% filter(Sample=="Blastocyst_rep2") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        s=T, wo=T) %>% filter(V15>=50) %>%
  dplyr::mutate(peak1=V4,peak2=V11,seqnames=V1,strand=V6,start1=as.integer(V2),start2=as.integer(V9),end1=as.integer(V3),end2=as.integer(V10)) %>%
  mutate(start=case_when(start1>=start2 ~ start1, .default = start2),end=case_when(end1<=end2 ~ end1, .default = end2)) %>%
  dplyr::select(seqnames,start,end,strand) %>% as.data.table() %>% mutate(score=1000,name=paste(seqnames,start,end,sep="_")) %>%
  dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample="shared")
dt.sm6APeak.shared.m6a.Benchmarkoverlap <- foreach(c = "shared", .combine='rbind')%do%{
  bed.peak <- dt.sm6APeak.shared.m6A  %>% dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(name=paste(name,strand,sep="_"))
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
dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                        nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
  left_join(x=.,y=dt.sm6APeak.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
  left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
  mutate(Group="sm6APeak.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)
# Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite           Group       TPR       FPR
# 1: shared     GLORI_mESC   41640 6016          2375      8613              96094 sm6APeak.shared 0.4333257 0.3015210
# 2: shared   eTAMseq_mESC   28339 5135          2442      8613              61365 sm6APeak.shared 0.4618105 0.4038082

#annotate of GSE and sm6APeak gene to gene
dt.GSE192440.peaks.annotgene <- annot_peak(peak.bed= rbind(dt.GSE192440.m6a %>% dplyr::distinct(seqnames,start,end,name), dt.GSE192440.shared.m6A %>% dplyr::distinct(seqnames,start,end,name)) %>%
                                             mutate(score=1000,strand="*") %>% dplyr::select(seqnames,start,end,name,score,strand),
                                                       strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                                       annot_type = "gene")
dt.sm6APeak.peaks.annotgene <- annot_peak(peak.bed= rbind(dt.sm6APeak.m6a %>% dplyr::distinct(seqnames,start,end,name,strand), dt.sm6APeak.shared.m6A %>% dplyr::distinct(seqnames,start,end,name,strand)) %>%
                                            mutate(score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand),
                                           strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                           annot_type = "gene")

dt.sm6APeak.shared.peaks.annotgene <- dt.sm6APeak.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.peaks.annotgene %>% dplyr::select(seqnames,start,end,strand,OverlappedGenes)) %>%
  tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()
dt.GSE192440.shared.peaks.annotgene <- dt.GSE192440.shared.m6A %>% left_join(x=.,y=dt.GSE192440.peaks.annotgene %>% dplyr::select(seqnames,start,end,strand,OverlappedGenes)) %>%
  tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()

#the comparison of count of m6A peaks (rep1, rep2 and shared peak) of sm6APeak with GSE192440
dt.count.shared.m6A.m6AGene.Blastocyst <- rbind(data.table(Group=c("sm6APeak"),nPeak=nrow(dt.sm6APeak.shared.m6A), nGene=n_distinct(dt.sm6APeak.shared.peaks.annotgene[OverlappedGenes!="." & !is.na(OverlappedGenes),OverlappedGenes])),
                                     data.table(Group=c("published"),nPeak=nrow(dt.GSE192440.shared.m6A), nGene=n_distinct(dt.GSE192440.shared.peaks.annotgene[OverlappedGenes!="." & !is.na(OverlappedGenes),OverlappedGenes]))) %>%
  pivot_longer(cols = c("nPeak","nGene"), names_to = "Metric", values_to = "Count") %>% as.data.table() %>% mutate(Label=paste0("N=",Count))

p.count.shared.m6A.m6AGene.Blastocyst <- ggplot(data=dt.count.shared.m6A.m6AGene.Blastocyst , aes(x=Group,y=Count,fill=Metric))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.shared.m6A.m6AGene.Blastocyst  %>% dplyr::filter(Metric=="nGene"),
            aes(label = Label, y = Count*0.5,x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.shared.m6A.m6AGene.Blastocyst  %>% dplyr::filter(Metric=="nPeak"),
            aes(label = Label, y = Count*0.5,x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#256da0","#a4c8d5"),breaks=c("nGene","nPeak"))+
  scale_y_continuous(expand = expansion(add = c(0,0.05)))+
  labs(x=NULL,y=str_wrap("Count of consistent m6A peak or m6A+ gene",width=45))+
  guides(fill=guide_legend(title = NULL))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.shared.m6A.m6AGene.Blastocyst
#the comparison of TPR and FPR of sm6APeak with GSE192440
dt.TPR.FPR.shared.m6A.Blastocyst <- rbind(dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                                 nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                            as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                            left_join(x=.,y=dt.sm6APeak.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                            left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                            mutate(Group="sm6APeak.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group) %>% dplyr::select(benchmark_name,Group,TPR,FPR),
                                          dt.GSE192440.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                                  nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                            as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                            left_join(x=.,y=dt.GSE192440.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                            left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                            mutate(Group="GSE192440.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)%>% dplyr::select(benchmark_name,Group,TPR,FPR))
dt.TPR.FPR.shared.m6A.Blastocyst <- dt.TPR.FPR.shared.m6A.Blastocyst %>% mutate(Group=factor(Group, levels=c("GSE192440.shared","sm6APeak.shared"), labels=c("published","sm6APeak"))) %>%
  pivot_longer(cols = c("FPR","TPR"),names_to = "Metric",values_to = "Value") %>% as.data.table() %>% mutate(benchmark_name=factor(benchmark_name,levels=c("GLORI_mESC","eTAMseq_mESC")))


p.TPR.FPR.shared.m6A.Blastocyst <- ggplot(data=dt.TPR.FPR.shared.m6A.Blastocyst ,aes(x=Metric,y=Value,fill=Group))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.TPR.FPR.shared.m6A.Blastocyst %>% dplyr::filter(Group=="published"),
            aes(label = round(Value,2), y = Value+0.02, x=Metric),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data =  dt.TPR.FPR.shared.m6A.Blastocyst %>% dplyr::filter(Group=="sm6APeak"),
            aes(label = round(Value,2), y = Value+0.02, x=Metric),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
  scale_fill_manual(values=c("#efb091","#eb4601"),breaks=c("published","sm6APeak"))+
  labs(x=NULL,y=str_wrap("Value of surrogate FPR/TPR ",width=40))+
  guides(fill=guide_legend(title = "Peak"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  facet_wrap(~benchmark_name)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.TPR.FPR.shared.m6A.Blastocyst
#the ratio of shared novel m6A peaks and m6A genes identified by sm6APeak
dt.sm6APeak.novel.shared.m6A <- bt.intersect(a=dt.sm6APeak.shared.m6A %>% dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(name=paste(name,strand,sep="_")),
                                             b=dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand) %>% distinct(),
                                             s=F, wao=T) %>% as.data.table() %>% dplyr::filter(V10 == ".") %>% dplyr::select(seqnames=V1,start=V2,end=V3,name=V4,score=V5,strand=V6)
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation <- dt.sm6APeak.novel.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(name) %>%
                                                                                                 mutate(nBenchSite=max(nBenchSite,na.rm=T),nConfidentBenchSite=max(nConfidentBenchSite,na.rm=T)) %>%
                                                                                                 as.data.table() %>% distinct(name,nBenchSite,nConfidentBenchSite), by="name")
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[,.N,by=.(nBenchSite>0)]#1133/2498 = 45%
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[,.N,by=.(nConfidentBenchSite>0)]#201/2498 = 8%

dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.shared.peaks.annotgene %>% mutate(name=paste(name,strand,sep="_")) %>% dplyr::distinct(name,OverlappedGenes),by="name" )
dt.GSE192440.peaks.annotgene.expanded <- dt.GSE192440.peaks.annotgene %>% tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()
dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6AGene %>% mutate(IsNovelGene=!OverlappedGenes %in% dt.GSE192440.peaks.annotgene.expanded$OverlappedGenes)
dt.sm6APeak.novel.shared.m6AGene[,.(NGene=n_distinct(OverlappedGenes)),by=.(IsNovelGene)]
#    IsNovelGene NGene
# 1:       FALSE  1505
# 2:        TRUE   866
dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6AGene %>% left_join(x=.,y=dt.sm6APeak.novel.shared.m6A.BenchmarkValidation %>% dplyr::select(name,nBenchSite, nConfidentBenchSite))
dt.sm6APeak.novel.shared.m6AGene[,.(NGene=n_distinct(OverlappedGenes)),by=.(BenchValidated=nBenchSite>0, ConfidentBenchValidated=nConfidentBenchSite>0,IsNovelGene)] %>% dplyr::arrange(IsNovelGene,BenchValidated,ConfidentBenchValidated)
#    BenchValidated ConfidentBenchValidated IsNovelGene NGene
# 1:           TRUE                   FALSE       FALSE   602
# 2:           TRUE                    TRUE       FALSE   126
# 3:             NA                      NA       FALSE   938
# 4:           TRUE                   FALSE        TRUE   406
# 5:           TRUE                    TRUE        TRUE   101
# 6:             NA                      NA        TRUE   415

dt.sm6APeak.novel.m6A.peak.m6AGene <- rbind(data.table(Type=c("BenchmarkSite","ConfidentBenchmarkSite"),Count=c(nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[nBenchSite>0,]),
                                                                                                                nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[nConfidentBenchSite>0,]))) %>%
                                              mutate(Group="NovelPeak", Ratio=Count/nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation), Label=paste0(Count,"/",nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation))),
                                            data.table(Type=c("BenchmarkSite","ConfidentBenchmarkSite"), Count=c(n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE & nBenchSite>0,OverlappedGenes]),
                                                                                                                 n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE & nConfidentBenchSite>0,OverlappedGenes])))%>%
                                              mutate(Group="Novelm6A+Gene",Ratio=Count/n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE,OverlappedGenes]), Label=paste0(Count,"/",n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE,OverlappedGenes])))
                                            ) %>%
  mutate(Type=factor(Type,levels=c("BenchmarkSite","ConfidentBenchmarkSite")),Group=factor(Group,levels=c("NovelPeak","Novelm6A+Gene")))

p.count.sm6APeak.novel.peak.m6AGene  <-
  ggplot(data=dt.sm6APeak.novel.m6A.peak.m6AGene, aes(x=Group,y=Ratio*100, fill=Type))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.sm6APeak.novel.m6A.peak.m6AGene %>% dplyr::filter(Type=="BenchmarkSite"),
            aes(label = Label, y =  Ratio*100/2+0.02, x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.sm6APeak.novel.m6A.peak.m6AGene %>% dplyr::filter(Type=="ConfidentBenchmarkSite"),
            aes(label = Label, y =  Ratio*100/2+0.02, x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("BenchmarkSite","ConfidentBenchmarkSite"))+
  labs(x=NULL,y=str_wrap("% of novel peak or m6A+ genes",width=40))+
  guides(fill=guide_legend(title = "ValidatedBy",direction = "vertical"))+
  scale_y_continuous(expand = expansion(add=c(0,0)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.sm6APeak.novel.peak.m6AGene

#the IGV plots of selected novel m6A gene identified by sm6APeak
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

dt.sm6APeak.novel.shared.m6AGene


selected.Novel.m6Agene.Blastocyst <-  dt.sm6APeak.novel.shared.m6AGene  %>% dplyr::filter(nConfidentBenchSite>0 & IsNovelGene==TRUE) %>%
  distinct(OverlappedGenes,IsNovelGene, nBenchSite, nConfidentBenchSite) %>%  dplyr::arrange(desc(nConfidentBenchSite)) %>%  slice_head(n=10)
dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==c("Ing5","Dda1","Llgl1","Srr","Srsf5")[3])
StudySamples <- c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2")
library(BSgenome.Mmusculus.UCSC.mm39)

## Create a page (7.5*7.5cm)
window.size=20#100kb
pseudocount=1
track.height=0.5
pdf("Fig4_sm6APeak_novel_m6A_Llgl1_in_GSE192440_Blastocyst.pdf")
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
    Input.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
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
    Input.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
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
dt.peak.to.see <- rbind(dt.sm6APeak.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="sm6APeak"),
                        dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="published") %>% mutate(strand=dt.selected.gene$strand)) %>%
  mutate(Group=factor(Group, levels=c("sm6APeak","published"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand) %>%
  mutate(Sample=factor(Sample,levels=c("Blastocyst_rep1","Blastocyst_rep2"),labels=c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2"))) %>%
  mutate(strand=as.character(strand))
dt.peak.to.see %>% filter(Group=="published")
for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  plotRanges(data = dt.peak.region, params = region, order = "random",
                                     fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),
                                     strandSplit = F,
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
legendPlot.peak <- plotLegend(legend = c("sm6APeak","published"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                                         just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                         title = expression("PeakType"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==c("Ing5","Dda1","Llgl1","Srr","Srsf5")[2])
StudySamples <- c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2")
library(BSgenome.Mmusculus.UCSC.mm39)

## Create a page (7.5*7.5cm)
window.size=20#100kb
pseudocount=1
track.height=0.5
pdf("Fig4_sm6APeak_novel_m6A_Dda1_in_GSE192440_Blastocyst.pdf")
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
    Input.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
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
    Input.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
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
dt.peak.to.see <- rbind(dt.sm6APeak.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="sm6APeak"),
                        dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="published") %>% mutate(strand=dt.selected.gene$strand)) %>%
  mutate(Group=factor(Group, levels=c("sm6APeak","published"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand) %>%
  mutate(Sample=factor(Sample,levels=c("Blastocyst_rep1","Blastocyst_rep2"),labels=c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2"))) %>%
  mutate(strand=as.character(strand))
dt.peak.to.see %>% filter(Group=="published")
for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  plotRanges(data = dt.peak.region, params = region, order = "random",
             fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),
             strandSplit = F,
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
legendPlot.peak <- plotLegend(legend = c("sm6APeak","published"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                              title = expression("PeakType"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()


#save intermediate variables for re-use
save.image("sm6Peak_Figure4.intermediate.results.RDS")
# q(save="no")
####################### combine all figure 4 together #############################
fig4.list <- list(p.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC=p.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC,
                  p.heatmap.TPRImprovement.BestOption.FPR20.peak.HEK293.mESC=p.heatmap.TPRImprovement.BestOption.FPR20.peak.HEK293.mESC,
                  p.MCA.top.dimension.NEB=p.MCA.top.dimension.NEB,
                  p.FPR.TPR.top.dimension.NEB=p.FPR.TPR.top.dimension.NEB,
                  p.TPR.FPR.sm6APeak.Method.BestOption=p.TPR.FPR.sm6APeak.Method.BestOption,
                  p.heatmap.AUC.sm6APeak.FourModel=p.heatmap.AUC.sm6APeak.FourModel,
                  p.dumbbell.optimal.cutoff.sm6APeak_NEB=p.dumbbell.optimal.cutoff.sm6APeak_NEB,
                  p.count.shared.m6A.m6AGene.Blastocyst=p.count.shared.m6A.m6AGene.Blastocyst,
                  p.TPR.FPR.shared.m6A.Blastocyst=p.TPR.FPR.shared.m6A.Blastocyst,
                  p.count.sm6APeak.novel.peak.m6AGene=p.count.sm6APeak.novel.peak.m6AGene
                  )
saveRDS(fig4.list, file="Figure4.plot.list.RDS")

pdf("Figure4.sm6APeak_enable_comprehensive_m6A_calling_with_high_accuracy_for_MeRIPseq_data.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.heatmap.similarity.BestOption.FPR20.peak.HEK293.mESC, x = 0.2, y=0.05, default.units = "cm",width = 6.8, height = 4.7)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 7.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.heatmap.TPRImprovement.BestOption.FPR20.peak.HEK293.mESC, x = 7.4, y=0.05, default.units = "cm",width = 10, height = 4.7)
#row2
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 4.9,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.MCA.top.dimension.NEB, x = 0.2, y=4.9, default.units = "cm",width = 5.5, height = 4.7)
plotText(label = "D", fontsize = 8, fontface = "bold",x =5.9, y = 4.9,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.FPR.TPR.top.dimension.NEB, x = 6.1, y=4.9, default.units = "cm",width = 5.5, height = 4.7)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 11.9, y = 4.9,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot =p.heatmap.AUC.sm6APeak.FourModel, x =12.1, y=5, default.units = "cm",width = 5.2, height = 4.7)
#row3
plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05, y = 9.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.sm6APeak.Method.BestOption, x = 0.2, y=9.7, default.units = "cm",width = 10, height = 4.7)
plotText(label = "G", fontsize = 8, fontface = "bold",x = 10.2, y = 9.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.dumbbell.optimal.cutoff.sm6APeak_NEB, x = 10.4, y=9.7, default.units = "cm",width = 5.8, height = 4.7)
#row4
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05, y = 14.5,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.shared.m6A.m6AGene.Blastocyst, x = 0.2, y=14.5, default.units = "cm",width = 4.5, height = 4.7)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 4.9, y = 14.5,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.shared.m6A.Blastocyst, x = 5.1, y=14.5, default.units = "cm",width = 6.5, height = 4.7)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 11.6, y = 14.5,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.sm6APeak.novel.peak.m6AGene, x = 11.8, y=14.5, default.units = "cm",width = 4.5, height = 4.7)
#row5
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05, y = 19.2,just = c("top","left"), default.units = "cm",draw=T)

dev.off()


### PartS4 Integration via machine learning enable comprehensive m6A calling ###############

#S4A the better performance AUC of sm6APeak over each method's best option
#load TPR FPR of NEB_HEK and NEB_mESC of sm6APeak
dt.FPR.TPR.FourMode.sm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
dt.FPP.TPR.sm6APeak.nonNEB <- dt.FPR.TPR.FourMode.sm6APeak %>% dplyr::filter(SelectedMode!="NEB" & !str_detect(SelectedSample,"NEB") & SelectedMode==SelectedSample) %>%
  mutate(Cell=factor(SelectedSample,levels=c("Abcam","SYSY","Millipore"), labels=c("Abcam_HEK293","SYSY_HEK293","Millipore_mESC")))
#load the Method's  best option of nonNEB
#Millipore
dt.pruned.TPR.FPR.Method.AllOption.Millipore <- readRDS("dt.pruned.FPR.TPR.AUC.Millipore.mESC.RDS")
dt.pruned.TPR.FPR.Method.BestOption.Millipore <- dt.pruned.TPR.FPR.Method.AllOption.Millipore %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=4) %>%
  inner_join(x=.,y=dt.pruned.TPR.FPR.Method.AllOption.Millipore)
dt.pruned.TPR.FPR.Method.BestOption.Millipore %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20)
#combine sm6APeak with those top4 best combo
dt.TPR.FPR.sm6APeak.Method.BestOption.Millipore <- rbind(dt.FPP.TPR.sm6APeak.nonNEB %>% filter(Cell=="Millipore_mESC") %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell) %>% mutate(Method="sm6APeak"),
                                                         dt.pruned.TPR.FPR.Method.BestOption.Millipore %>% mutate(Cell="Millipore_mESC") %>%  dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell,Method) )
#Abcam
dt.pruned.TPR.FPR.Method.AllOption.Abcam <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.RDS")
dt.pruned.TPR.FPR.Method.BestOption.Abcam <- dt.pruned.TPR.FPR.Method.AllOption.Abcam %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=4) %>%
  inner_join(x=.,y=dt.pruned.TPR.FPR.Method.AllOption.Abcam)
dt.pruned.TPR.FPR.Method.BestOption.Abcam %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20)
#combine sm6APeak with those top4 best combo
dt.TPR.FPR.sm6APeak.Method.BestOption.Abcam <- rbind(dt.FPP.TPR.sm6APeak.nonNEB %>% filter(Cell=="Abcam_HEK293") %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell) %>% mutate(Method="sm6APeak"),
                                                     dt.pruned.TPR.FPR.Method.BestOption.Abcam %>% mutate(Cell="Abcam_HEK293") %>%  dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell,Method) )
#SYSY
dt.pruned.TPR.FPR.Method.AllOption.SYSY <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.RDS")
dt.pruned.TPR.FPR.Method.BestOption.SYSY <- dt.pruned.TPR.FPR.Method.AllOption.SYSY %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=4) %>%
  inner_join(x=.,y=dt.pruned.TPR.FPR.Method.AllOption.SYSY)
dt.pruned.TPR.FPR.Method.BestOption.SYSY %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20)
#combine sm6APeak with those top4 best combo
dt.TPR.FPR.sm6APeak.Method.BestOption.SYSY <- rbind(dt.FPP.TPR.sm6APeak.nonNEB %>% filter(Cell=="SYSY_HEK293") %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell) %>% mutate(Method="sm6APeak"),
                                                    dt.pruned.TPR.FPR.Method.BestOption.SYSY %>% mutate(Cell="SYSY_HEK293") %>%  dplyr::select(FPR,TPR,AUC.scaled,TPR20,Cell,Method) )

#combine three non NEB together
dt.TPR.FPR.sm6APeak.Method.BestOption.nonNEB <- rbind(dt.TPR.FPR.sm6APeak.Method.BestOption.Abcam,dt.TPR.FPR.sm6APeak.Method.BestOption.SYSY,dt.TPR.FPR.sm6APeak.Method.BestOption.Millipore)

dt.TPR.FPR.sm6APeak.Method.BestOption.label.nonNEB <- dt.TPR.FPR.sm6APeak.Method.BestOption.nonNEB %>% distinct(Method,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>% group_by(Cell) %>%
  mutate(ypos=tail(seq(0.7,0,by= -0.6/5),5)*max(dt.TPR.FPR.sm6APeak.Method.BestOption.nonNEB$TPR)) %>% mutate(xpos=0.6) %>% as.data.table()
p.TPR.FPR.sm6APeak.Method.BestOption.nonNEB  <-  ggplot(data=dt.TPR.FPR.sm6APeak.Method.BestOption.nonNEB ,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of sm6APeak or method's BestCombo",y=str_wrap("Surrogate TPR of sm6APeak or method's BestCombo",width = 40))+
  geom_text(data=dt.TPR.FPR.sm6APeak.Method.BestOption.label.nonNEB, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#eb4601",c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0")),breaks = c("sm6APeak",c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools")),guide="none")+
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
p.TPR.FPR.sm6APeak.Method.BestOption.nonNEB

#sm6APeak_NEB display comparable AUC in nonNEB samples
dt.FPR.TPR.FourMode.sm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
dt.FPP.TPR.sm6APeak.NEBvsOriginal <- dt.FPR.TPR.FourMode.sm6APeak %>% dplyr::filter(!str_detect(SelectedSample,"NEB")) %>% filter((SelectedMode==SelectedSample) | SelectedMode=="NEB") %>%
  mutate(Cell=factor(SelectedSample,levels=c("Abcam","SYSY","Millipore"), labels=c("Abcam_HEK293","SYSY_HEK293","Millipore_mESC"))) %>%
  mutate(ModelType=case_when(SelectedMode==SelectedSample ~ "Original", .default = SelectedMode))
dt.TPR.FPR.sm6APeak.label.NEBvsOriginal <- dt.FPP.TPR.sm6APeak.NEBvsOriginal %>% distinct(Cell,ModelType,AUC.scaled,TPR20) %>%  mutate(label=paste0(ModelType,"Model,AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>% group_by(Cell) %>%
  mutate(ypos=head(seq(0.7,0,by= -0.6/5),2)*max(dt.FPP.TPR.sm6APeak.NEBvsOriginal$TPR)) %>% mutate(xpos=0.5) %>% as.data.table()
p.TPR.FPR.sm6APeak.NEBvsOriginal  <-  ggplot(data=dt.FPP.TPR.sm6APeak.NEBvsOriginal ,aes(x=FPR,y=TPR,color=ModelType))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of sm6APeak from different anti-m6A model",y=str_wrap("Surrogate TPR of sm6APeak from different anti-m6A model",width = 40))+
  geom_text(data=dt.TPR.FPR.sm6APeak.label.NEBvsOriginal, aes(x=xpos,y=ypos,color=ModelType,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#fcbb44","#f0766d"),breaks = c("NEB","Original"),guide="none")+
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
p.TPR.FPR.sm6APeak.NEBvsOriginal

#e sm6APeak identify novel m6A in mESC and HEK293 compare with MACS2_BestComboOption (keep this in supplementary figures)

#f sm6APeak reveal novel m6A modification in blastocyst than reported in published study (GSE192440), could be validated using GLORI and eTAMseq
#load published m6A peaks
GSE192440.bed.files <- list.files("/data2/mESC_m6A/Analysis/",pattern = "(.bed)$",full.names = T)
names(GSE192440.bed.files) <- list.files("/data2/mESC_m6A/Analysis/",pattern = "(.bed)$",full.names = F) %>% strsplit(split="_MS_",fixed=T) %>%
  sapply("[",2) %>% str_replace(pattern=".bed",replacement = "") %>% str_replace(pattern="_IP",replacement = "")
GSE192440.bed.files <- GSE192440.bed.files[grep(names(GSE192440.bed.files),pattern="Blastocyst")]
dt.GSE192440.m6a <- foreach(i=1:length(GSE192440.bed.files), .combine = 'rbind')%do%{
  dt.m6a <- fread(GSE192440.bed.files[i]) %>% dplyr::select(seqnames=V1,start=V2,end=V3,score=V7) %>% distinct() %>% group_by(seqnames,start,end) %>% mutate(score=mean(score)) %>% as.data.table() %>% distinct()
  dt.m6a <- dt.m6a %>% mutate(strand="*",Sample=names(GSE192440.bed.files)[i]) %>% dplyr::filter(str_detect(seqnames,"chr"))
  dt.m6a
}
#lifover mm10 peaks to mm39 peaks
m6A_liftover <- function(
    chainfile="~/genome_db/LiftOverChain/mm10ToMm39.over.chain",
    input.bed=dt.GSE192440.m6a[Sample=="2cell_rep1",]%>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(seqnames,start,end,name,score,strand)
){
  chain <- rtracklayer::import.chain(chainfile)
  input.GR <- GenomicRanges::makeGRangesFromDataFrame(input.bed,keep.extra.columns = T)
  as.data.table(unlist(rtracklayer::liftOver(input.GR, chain))) %>% dplyr::select(seqnames,start,end,name,score,strand)
}
dt.GSE192440.m6a.liftover <- m6A_liftover(chainfile="~/genome_db/LiftOverChain/mm10ToMm39.over.chain",
                                          input.bed=dt.GSE192440.m6a %>%
                                            mutate(name=paste(seqnames,start,end,sep="_")) %>% distinct(seqnames,start,end,name,score,strand))
dt.GSE192440.m6a <- dt.GSE192440.m6a %>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(name,score,Sample) %>%
  inner_join(x=.,y=dt.GSE192440.m6a.liftover %>% dplyr::select(seqnames,start,end,name,strand) %>% distinct(), by=c("name")) %>%
  mutate(mm10.name=name) %>% mutate(name=paste(seqnames,start,end,sep="_")) %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,everything())

#load M6APeakS for GSE192440
# sm6APeak.peak.files <- list.files("/data2/mESC_m6A/M6APeakS/",pattern="(GSE192440).+(Millipore_default.tsv)$",full.names = T) %>% grep(pattern="Bc",value=T)
# names(sm6APeak.peak.files) <- list.files("/data2/mESC_m6A/M6APeakS/",pattern="(GSE192440).+(Millipore_default.tsv)$",full.names = F) %>% grep(pattern="Bc",value=T) %>%
#   strsplit(split="GSE192440_") %>% sapply("[",2) %>% strsplit(split="_M6APeakS_peak") %>% sapply("[",1)
# dt.sm6APeak.m6a <- foreach(i = 1:length(sm6APeak.peak.files), .combine='rbind')%do%{
#   fread(sm6APeak.peak.files[i]) %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,OverlappedGenes) %>% distinct() %>% mutate(Sample=names(sm6APeak.peak.files)[i])
# }
# dt.sm6APeak.m6a <- dt.sm6APeak.m6a %>%  mutate(Sample=factor(Sample,levels=c("Bc_rep1","Bc_rep2"),labels=c("Blastocyst_rep1","Blastocyst_rep2")))
# dt.sm6APeak.m6a[,.N,by=.(Sample)]
sm6APeak.peak.files <- list.files("/data3/mESC_m6A/sm6APeak_GSE192440_Bc/",pattern="(GSE192440).+(NEB_default.tsv)$",full.names = T)
names(sm6APeak.peak.files) <- list.files("/data3/mESC_m6A/sm6APeak_GSE192440_Bc/",pattern="(GSE192440).+(NEB_default.tsv)$",full.names = F) %>%
  strsplit(split="GSE192440_") %>% sapply("[",2) %>% strsplit(split="_M6APeakS_peak") %>% sapply("[",1)
dt.sm6APeak.m6a <- foreach(i = 1:length(sm6APeak.peak.files), .combine='rbind')%do%{
  fread(sm6APeak.peak.files[i]) %>% dplyr::select(seqnames,start,end,name,score,strand,Sample,OverlappedGenes) %>% distinct() %>% mutate(Sample=names(sm6APeak.peak.files)[i])
}
dt.sm6APeak.m6a <- dt.sm6APeak.m6a %>%  mutate(Sample=factor(Sample,levels=c("Bc_rep1","Bc_rep2"),labels=c("Blastocyst_rep1","Blastocyst_rep2")))
dt.sm6APeak.m6a[,.N,by=.(Sample)]


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

dt.sm6APeak.m6a.Benchmarkoverlap <- foreach(c = unique(dt.sm6APeak.m6a$Sample), .combine='rbind')%do%{
  bed.peak <- dt.sm6APeak.m6a %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.bench <- dt.benchmark.m6A.mESC %>% dplyr::select(seqnames,start,end,name,score,strand,benchmark_name)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.bench, s=T,F=1, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% dplyr::select(name=V4,site=V10,benchmark_name=V13) %>% distinct() %>% dplyr::filter(benchmark_name!=".")
  #benchmark sites with high confidence
  dt.overlap <- dt.overlap %>% mutate(confident_site=site %in% dt.eTAMseq.mESC.confident[PerturbEffect < -0.1,name])
  dt.overlap <- dt.overlap %>% group_by(name,benchmark_name) %>% mutate(nBenchSite=n_distinct(site),nConfidentBenchSite=n_distinct(site[confident_site])) %>%
    as.data.table() %>% distinct(name,benchmark_name,nBenchSite,nConfidentBenchSite)
  dt.overlap %>% mutate(Sample=c)
}
dt.sm6APeak.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident)

dt.GSE192440.m6a.Benchmarkoverlap <- foreach(c = unique(dt.GSE192440.m6a$Sample), .combine='rbind')%do%{
  bed.peak <- dt.GSE192440.m6a %>% dplyr::filter(Sample==c) %>% dplyr::select(seqnames,start,end,name,score,strand)
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
dt.GSE192440.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                 nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident)
dt.TPR.FPR.GSE192440.sm6APeak <- rbind(dt.GSE192440.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                        nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                         as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                         left_join(x=.,y=dt.GSE192440.m6a[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                         left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                         mutate(Group="GSE192440"),
                                       dt.sm6APeak.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                       nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                         as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                         left_join(x=.,y=dt.sm6APeak.m6a[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                         left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                         mutate(Group="sm6APeak")
)
dt.TPR.FPR.GSE192440.sm6APeak %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)

#identify the replicate shared peaks for sm6APeak and GSE192440
dt.GSE192440.shared.m6A <- bt.intersect(a=dt.GSE192440.m6a %>% filter(Sample=="Blastocyst_rep1") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        b=dt.GSE192440.m6a %>% filter(Sample=="Blastocyst_rep2") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                        s=F, wo=T) %>% filter(V15>=50) %>%
  dplyr::mutate(peak1=V4,peak2=V11,seqnames=V1,start1=as.integer(V2),start2=as.integer(V9),end1=as.integer(V3),end2=as.integer(V10)) %>%
  mutate(start=case_when(start1>=start2 ~ start1, .default = start2),end=case_when(end1<=end2 ~ end1, .default = end2)) %>%
  dplyr::select(seqnames,start,end) %>% as.data.table() %>% mutate(strand="*",score=1000,name=paste(seqnames,start,end,sep="_")) %>%
  dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample="shared")
dt.GSE192440.shared.m6a.Benchmarkoverlap <- foreach(c = "shared", .combine='rbind')%do%{
  bed.peak <- dt.GSE192440.shared.m6A  %>% dplyr::select(seqnames,start,end,name,score,strand)
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
dt.GSE192440.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                        nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
  left_join(x=.,y=dt.GSE192440.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
  left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
  mutate(Group="GSE192440.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)
#    Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite            Group       TPR       FPR
# 1: shared     GLORI_mESC   14471 3591          1263      5173              96094 GSE192440.shared 0.1505921 0.3058187
# 2: shared   eTAMseq_mESC   10168 2856          1312      5173              61365 GSE192440.shared 0.1656971 0.4479026
#sm6APeak
dt.sm6APeak.shared.m6A <- bt.intersect(a=dt.sm6APeak.m6a %>% filter(Sample=="Blastocyst_rep1") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                       b=dt.sm6APeak.m6a %>% filter(Sample=="Blastocyst_rep2") %>% dplyr::select(seqnames,start,end,name,score,strand,Sample),
                                       s=T, wo=T) %>% filter(V15>=50) %>%
  dplyr::mutate(peak1=V4,peak2=V11,seqnames=V1,strand=V6,start1=as.integer(V2),start2=as.integer(V9),end1=as.integer(V3),end2=as.integer(V10)) %>%
  mutate(start=case_when(start1>=start2 ~ start1, .default = start2),end=case_when(end1<=end2 ~ end1, .default = end2)) %>%
  dplyr::select(seqnames,start,end,strand) %>% as.data.table() %>% mutate(score=1000,name=paste(seqnames,start,end,sep="_")) %>%
  dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(Sample="shared")
dt.sm6APeak.shared.m6a.Benchmarkoverlap <- foreach(c = "shared", .combine='rbind')%do%{
  bed.peak <- dt.sm6APeak.shared.m6A  %>% dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(name=paste(name,strand,sep="_"))
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
dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                       nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
  as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
  left_join(x=.,y=dt.sm6APeak.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
  left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
  mutate(Group="sm6APeak.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)
#    Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite           Group       TPR       FPR
# 1: shared     GLORI_mESC   41640 6016          2375      8613              96094 sm6APeak.shared 0.4333257 0.3015210
# 2: shared   eTAMseq_mESC   28339 5135          2442      8613              61365 sm6APeak.shared 0.4618105 0.4038082
#    Sample benchmark_name SumSite  nTP nTP.confident TotalPeak TotalBenchmarkSite           Group       TPR       FPR
# 1: shared     GLORI_mESC   50539 8707          2830     17018              96094 sm6APeak.shared 0.5259329 0.4883653
# 2: shared   eTAMseq_mESC   35051 7196          2907     17018              61365 sm6APeak.shared 0.5711888 0.5771536

#annotate of GSE and sm6APeak gene to gene
dt.GSE192440.peaks.annotgene <- annot_peak(peak.bed= rbind(dt.GSE192440.m6a %>% dplyr::distinct(seqnames,start,end,name), dt.GSE192440.shared.m6A %>% dplyr::distinct(seqnames,start,end,name)) %>%
                                             mutate(score=1000,strand="*") %>% dplyr::select(seqnames,start,end,name,score,strand),
                                           strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                           annot_type = "gene")
dt.sm6APeak.peaks.annotgene <- annot_peak(peak.bed= rbind(dt.sm6APeak.m6a %>% dplyr::distinct(seqnames,start,end,name,strand), dt.sm6APeak.shared.m6A %>% dplyr::distinct(seqnames,start,end,name,strand)) %>%
                                            mutate(score=1000) %>% dplyr::select(seqnames,start,end,name,score,strand),
                                          strand=F, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                          annot_type = "gene")

dt.sm6APeak.shared.peaks.annotgene <- dt.sm6APeak.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.peaks.annotgene %>% dplyr::select(seqnames,start,end,strand,OverlappedGenes)) %>%
  tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()
dt.GSE192440.shared.peaks.annotgene <- dt.GSE192440.shared.m6A %>% left_join(x=.,y=dt.GSE192440.peaks.annotgene %>% dplyr::select(seqnames,start,end,strand,OverlappedGenes)) %>%
  tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()

#the comparison of count of m6A peaks (rep1, rep2 and shared peak) of sm6APeak with GSE192440
dt.count.shared.m6A.m6AGene.Blastocyst <- rbind(data.table(Group=c("sm6APeak"),nPeak=nrow(dt.sm6APeak.shared.m6A), nGene=n_distinct(dt.sm6APeak.shared.peaks.annotgene[OverlappedGenes!="." & !is.na(OverlappedGenes),OverlappedGenes])),
                                                data.table(Group=c("published"),nPeak=nrow(dt.GSE192440.shared.m6A), nGene=n_distinct(dt.GSE192440.shared.peaks.annotgene[OverlappedGenes!="." & !is.na(OverlappedGenes),OverlappedGenes]))) %>%
  pivot_longer(cols = c("nPeak","nGene"), names_to = "Metric", values_to = "Count") %>% as.data.table() %>% mutate(Label=paste0("N=",Count))

p.count.shared.m6A.m6AGene.Blastocyst <- ggplot(data=dt.count.shared.m6A.m6AGene.Blastocyst , aes(x=Group,y=Count,fill=Metric))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.shared.m6A.m6AGene.Blastocyst  %>% dplyr::filter(Metric=="nGene"),
            aes(label = Label, y = Count*0.5,x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.shared.m6A.m6AGene.Blastocyst  %>% dplyr::filter(Metric=="nPeak"),
            aes(label = Label, y = Count*0.5,x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#256da0","#a4c8d5"),breaks=c("nGene","nPeak"))+
  scale_y_continuous(expand = expansion(add = c(0,0.05)))+
  labs(x=NULL,y=str_wrap("Count of consistent m6A peak or m6A+ gene",width=45))+
  guides(fill=guide_legend(title = NULL))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.shared.m6A.m6AGene.Blastocyst
#the comparison of TPR and FPR of sm6APeak with GSE192440
dt.TPR.FPR.shared.m6A.Blastocyst <- rbind(dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                                 nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                            as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                            left_join(x=.,y=dt.sm6APeak.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                            left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                            mutate(Group="sm6APeak.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group) %>% dplyr::select(benchmark_name,Group,TPR,FPR),
                                          dt.GSE192440.shared.m6a.Benchmarkoverlap %>% group_by(Sample,benchmark_name) %>% mutate(SumSite=sum(nBenchSite),nTP=n_distinct(name[nBenchSite>0]),
                                                                                                                                  nTP.confident=n_distinct(name[nConfidentBenchSite>0])) %>%
                                            as.data.table() %>% distinct(Sample,benchmark_name,SumSite,nTP,nTP.confident) %>%
                                            left_join(x=.,y=dt.GSE192440.shared.m6A[,.(TotalPeak=n_distinct(name)),by=.(Sample)], by="Sample") %>%
                                            left_join(x=.,y=dt.benchmark.m6A.mESC[,.(TotalBenchmarkSite=n_distinct(name)),by=.(benchmark_name)],by="benchmark_name") %>%
                                            mutate(Group="GSE192440.shared") %>% mutate(TPR=SumSite/TotalBenchmarkSite,FPR=1-nTP/TotalPeak) %>% dplyr::arrange(Sample,benchmark_name,Group)%>% dplyr::select(benchmark_name,Group,TPR,FPR))
dt.TPR.FPR.shared.m6A.Blastocyst <- dt.TPR.FPR.shared.m6A.Blastocyst %>% mutate(Group=factor(Group, levels=c("GSE192440.shared","sm6APeak.shared"), labels=c("published","sm6APeak"))) %>%
  pivot_longer(cols = c("FPR","TPR"),names_to = "Metric",values_to = "Value") %>% as.data.table() %>% mutate(benchmark_name=factor(benchmark_name,levels=c("GLORI_mESC","eTAMseq_mESC")))

p.TPR.FPR.shared.m6A.Blastocyst <- ggplot(data=dt.TPR.FPR.shared.m6A.Blastocyst ,aes(x=Metric,y=Value,fill=Group))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.TPR.FPR.shared.m6A.Blastocyst %>% dplyr::filter(Group=="published"),
            aes(label = round(Value,2), y = Value+0.02, x=Metric),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data =  dt.TPR.FPR.shared.m6A.Blastocyst %>% dplyr::filter(Group=="sm6APeak"),
            aes(label = round(Value,2), y = Value+0.02, x=Metric),
            color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
  scale_fill_manual(values=c("#efb091","#eb4601"),breaks=c("published","sm6APeak"))+
  labs(x=NULL,y=str_wrap("Value of surrogate FPR/TPR ",width=40))+
  guides(fill=guide_legend(title = "Peak"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  facet_wrap(~benchmark_name)+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.TPR.FPR.shared.m6A.Blastocyst
#the ratio of shared novel m6A peaks and m6A genes identified by sm6APeak
dt.sm6APeak.novel.shared.m6A <- bt.intersect(a=dt.sm6APeak.shared.m6A %>% dplyr::select(seqnames,start,end,name,score,strand) %>% mutate(name=paste(name,strand,sep="_")),
                                             b=dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand) %>% distinct(),
                                             s=F, wao=T) %>% as.data.table() %>% dplyr::filter(V10 == ".") %>% dplyr::select(seqnames=V1,start=V2,end=V3,name=V4,score=V5,strand=V6)
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation <- dt.sm6APeak.novel.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.shared.m6a.Benchmarkoverlap %>% group_by(name) %>%
                                                                                                 mutate(nBenchSite=max(nBenchSite,na.rm=T),nConfidentBenchSite=max(nConfidentBenchSite,na.rm=T)) %>%
                                                                                                 as.data.table() %>% distinct(name,nBenchSite,nConfidentBenchSite), by="name")
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[,.N,by=.(nBenchSite>0)]#1133/2498 = 45%
dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[,.N,by=.(nConfidentBenchSite>0)]#201/2498 = 8%

dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6A %>% left_join(x=.,y=dt.sm6APeak.shared.peaks.annotgene %>% mutate(name=paste(name,strand,sep="_")) %>% dplyr::distinct(name,OverlappedGenes),by="name" )
dt.GSE192440.peaks.annotgene.expanded <- dt.GSE192440.peaks.annotgene %>% tidyr::separate_longer_delim(cols = "OverlappedGenes", delim = ",") %>% as.data.table()
dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6AGene %>% mutate(IsNovelGene=!OverlappedGenes %in% dt.GSE192440.peaks.annotgene.expanded$OverlappedGenes)
dt.sm6APeak.novel.shared.m6AGene[,.(NGene=n_distinct(OverlappedGenes)),by=.(IsNovelGene)]
#    IsNovelGene NGene
# 1:       FALSE  1505
# 2:        TRUE   866
dt.sm6APeak.novel.shared.m6AGene <- dt.sm6APeak.novel.shared.m6AGene %>% left_join(x=.,y=dt.sm6APeak.novel.shared.m6A.BenchmarkValidation %>% dplyr::select(name,nBenchSite, nConfidentBenchSite))
dt.sm6APeak.novel.shared.m6AGene[,.(NGene=n_distinct(OverlappedGenes)),by=.(BenchValidated=nBenchSite>0, ConfidentBenchValidated=nConfidentBenchSite>0,IsNovelGene)] %>% dplyr::arrange(IsNovelGene,BenchValidated,ConfidentBenchValidated)
#    BenchValidated ConfidentBenchValidated IsNovelGene NGene
# 1:           TRUE                   FALSE       FALSE   602
# 2:           TRUE                    TRUE       FALSE   126
# 3:             NA                      NA       FALSE   938
# 4:           TRUE                   FALSE        TRUE   406
# 5:           TRUE                    TRUE        TRUE   101
# 6:             NA                      NA        TRUE   415

dt.sm6APeak.novel.m6A.peak.m6AGene <- rbind(data.table(Type=c("BenchmarkSite","ConfidentBenchmarkSite"),Count=c(nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[nBenchSite>0,]),
                                                                                                                nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation[nConfidentBenchSite>0,]))) %>%
                                              mutate(Group="NovelPeak", Ratio=Count/nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation), Label=paste0(Count," / ",nrow(dt.sm6APeak.novel.shared.m6A.BenchmarkValidation))),
                                            data.table(Type=c("BenchmarkSite","ConfidentBenchmarkSite"), Count=c(n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE & nBenchSite>0,OverlappedGenes]),
                                                                                                                 n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE & nConfidentBenchSite>0,OverlappedGenes])))%>%
                                              mutate(Group="Novelm6A+Gene",Ratio=Count/n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE,OverlappedGenes]), Label=paste0(Count," / ",n_distinct(dt.sm6APeak.novel.shared.m6AGene[IsNovelGene==TRUE,OverlappedGenes])))
) %>%
  mutate(Type=factor(Type,levels=c("BenchmarkSite","ConfidentBenchmarkSite")),Group=factor(Group,levels=c("NovelPeak","Novelm6A+Gene")))

p.count.sm6APeak.novel.peak.m6AGene  <-
  ggplot(data=dt.sm6APeak.novel.m6A.peak.m6AGene, aes(x=Group,y=Ratio*100, fill=Type))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.sm6APeak.novel.m6A.peak.m6AGene %>% dplyr::filter(Type=="BenchmarkSite"),
            aes(label = Label, y =  Ratio*100/2+0.02, x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.sm6APeak.novel.m6A.peak.m6AGene %>% dplyr::filter(Type=="ConfidentBenchmarkSite"),
            aes(label = Label, y =  max(Ratio*100/2+0.02,10), x=Group),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("BenchmarkSite","ConfidentBenchmarkSite"))+
  labs(x=NULL,y=str_wrap("% of novel peak or m6A+ genes",width=40))+
  guides(fill=guide_legend(title = "ValidatedBy",direction = "vertical"))+
  scale_y_continuous(expand = expansion(add=c(0,0)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.sm6APeak.novel.peak.m6AGene

#the IGV plots of selected novel m6A gene identified by sm6APeak
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

dt.sm6APeak.novel.shared.m6AGene


selected.Novel.m6Agene.Blastocyst <-  dt.sm6APeak.novel.shared.m6AGene  %>% dplyr::filter(nConfidentBenchSite>0 & IsNovelGene==TRUE) %>%
  distinct(OverlappedGenes,IsNovelGene, nBenchSite, nConfidentBenchSite) %>%  dplyr::arrange(desc(nConfidentBenchSite)) %>%  slice_head(n=10)
dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==unique(selected.Novel.m6Agene.Blastocyst$OverlappedGenes)[2])
StudySamples <- c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2")
library(BSgenome.Mmusculus.UCSC.mm39)

## Create a page (7.5*7.5cm)
window.size=20#100kb
pseudocount=1
track.height=0.5
pdf(paste0("FigS4_sm6APeak_novel_m6A_", dt.selected.gene$gene_name, "_in_GSE192440_Blastocyst.pdf"))
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
    Input.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
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
    Input.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
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
dt.peak.to.see <- rbind(dt.sm6APeak.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="sm6APeak"),
                        dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="published") %>% mutate(strand=dt.selected.gene$strand)) %>%
  mutate(Group=factor(Group, levels=c("sm6APeak","published"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand) %>%
  mutate(Sample=factor(Sample,levels=c("Blastocyst_rep1","Blastocyst_rep2"),labels=c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2"))) %>%
  mutate(strand=as.character(strand))
dt.peak.to.see %>% filter(Group=="published")
for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  plotRanges(data = dt.peak.region, params = region, order = "random",
             fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),
             strandSplit = F,
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
legendPlot.peak <- plotLegend(legend = c("sm6APeak","published"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                              title = expression("PeakType"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")


dev.off()

selected.Novel.m6Agene.Blastocyst <-  dt.sm6APeak.novel.shared.m6AGene  %>% dplyr::filter(nConfidentBenchSite>0 & IsNovelGene==TRUE) %>%
  distinct(OverlappedGenes,IsNovelGene, nBenchSite, nConfidentBenchSite) %>%  dplyr::arrange(desc(nConfidentBenchSite)) %>%  slice_head(n=10)
dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==unique(selected.Novel.m6Agene.Blastocyst$OverlappedGenes)[3])
StudySamples <- c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2")
library(BSgenome.Mmusculus.UCSC.mm39)

## Create a page (7.5*7.5cm)
window.size=20#100kb
pseudocount=1
track.height=0.5
pdf(paste0("FigS4_sm6APeak_novel_m6A_", dt.selected.gene$gene_name, "_in_GSE192440_Blastocyst.pdf"))
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
    Input.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.neg.bw"),params=region)
    RIP.neg <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.neg.bw"),params=region)
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
    Input.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_Input.pos.bw"),params=region)
    RIP.pos <- plotgardener::readBigwig(file = paste0("/data3/mESC_m6A/BigWig/",StudySamples[s],"_RIP.pos.bw"),params=region)
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
dt.peak.to.see <- rbind(dt.sm6APeak.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="sm6APeak"),
                        dt.GSE192440.m6a %>% dplyr::select(seqnames,start,end,name,score,strand,Sample) %>% mutate(Group="published") %>% mutate(strand=dt.selected.gene$strand)) %>%
  mutate(Group=factor(Group, levels=c("sm6APeak","published"))) %>% dplyr::filter(seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand) %>%
  mutate(Sample=factor(Sample,levels=c("Blastocyst_rep1","Blastocyst_rep2"),labels=c("Millipore_GSE192440_Bc_rep1","Millipore_GSE192440_Bc_rep2"))) %>%
  mutate(strand=as.character(strand))
dt.peak.to.see %>% filter(Group=="published")
for(s in 1:length(StudySamples)){
  dt.peak.region <- dt.peak.to.see %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend) & strand==dt.selected.gene$strand)
  plotRanges(data = dt.peak.region, params = region, order = "random",
             fill = colorby("Group", palette =colorRampPalette(c("#e9abac","#7bb0d5"))),
             strandSplit = F,
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
legendPlot.peak <- plotLegend(legend = c("sm6APeak","published"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                              just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                              title = expression("PeakType"),fontface="plain")
#add rect for whole panel
plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

dev.off()

#identify the replicate-consistent m6A for published m6A and sm6APeak


#save intermediate variables for re-use
save.image("sm6APeak_FigureS4.intermediate.results.RDS")
# q(save="no")

####################### combine all figure S4 together
figS4.list <- list(p.TPR.FPR.sm6APeak.Method.BestOption.nonNEB=p.TPR.FPR.sm6APeak.Method.BestOption.nonNEB,
                   p.TPR.FPR.sm6APeak.NEBvsOriginal=p.TPR.FPR.sm6APeak.NEBvsOriginal,
                   p.count.shared.m6A.m6AGene.Blastocyst=p.count.shared.m6A.m6AGene.Blastocyst,
                   p.count.sm6APeak.novel.peak.m6AGene=p.count.sm6APeak.novel.peak.m6AGene
)
saveRDS(figS4.list, file="FigureS4.plot.list.RDS")

pdf("FigureS4.sm6APeak_enable_comprehensive_m6A_calling_with_high_accuracy_for_MeRIPseq_data_regardless_of_anti_m6A.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.sm6APeak.Method.BestOption.nonNEB, x = 0.05, y=0.2, default.units = "cm",width = 17, height = 5.5)
#row2
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05, y = 5.9,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.sm6APeak.NEBvsOriginal, x = 0.05, y=6.1, default.units = "cm",width = 17, height = 5.5)
#row3
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 11.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.shared.m6A.m6AGene.Blastocyst, x = 0.05, y=11.9, default.units = "cm",width = 4.5, height = 5.5)
# plotText(label = "D", fontsize = 8, fontface = "bold",x = 4.9, y = 14.5,just = c("top","left"), default.units = "cm",draw=T)
# plotGG(plot = p.TPR.FPR.shared.m6A.Blastocyst, x = 5.1, y=14.5, default.units = "cm",width = 6.5, height = 4.7)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 4.4, y = 11.7,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.count.sm6APeak.novel.peak.m6AGene, x = 4.4, y=11.9, default.units = "cm",width = 4.5, height = 5.5)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 9, y = 11.7,just = c("top","left"), default.units = "cm",draw=T)
plotRect( x = 9, y=11.9, just = c("top","left"), default.units = "cm",width = 8, height = 5.5)
# #row5
# plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05, y = 19.2,just = c("top","left"), default.units = "cm",draw=T)

dev.off()
