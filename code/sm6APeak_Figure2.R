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

### Part2 The parameter refinement  significantly improved the performance of m6A calling algorithm for MeRIP-seq data ###############
#2a The surrogate FPR and TPR when use different fold change cutoff
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE)
dt.FPR.TPR.NEB.HEK <- readRDS("dt.TPR.FPR.NEB.HEK.combo.RDS")
dt.FPR.TPR.NEB.mESC <- readRDS("dt.TPR.FPR.NEB.mESC.combo.RDS")
dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 <- rbind(dt.FPR.TPR.NEB.HEK %>% dplyr::filter(Method=="exomePeak2" & ComboName==dt.method.combo[Method=="exomePeak2" & IsDefaultCombo==TRUE,ComboName]),
                                                     dt.FPR.TPR.NEB.mESC %>% dplyr::filter(Method=="exomePeak2" & ComboName==dt.method.combo[Method=="exomePeak2" & IsDefaultCombo==TRUE,ComboName]))
dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 %>% distinct(Sample,Benchmark)
#           Sample                          Benchmark
# 1: HEK_NEB_mRNA1 dt.benchmark.GLORI_HEK293.m6As.bed
# 2: HEK_NEB_mRNA2 dt.benchmark.GLORI_HEK293.m6As.bed
# 3:      mESC_WT1   dt.benchmark.GLORI_mESC.m6As.bed
# 4:      mESC_WT2   dt.benchmark.GLORI_mESC.m6As.bed
dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 <- dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 %>% mutate(Cell=case_when(str_detect(Sample,"HEK")~ "HEK293", .default = "mESC"),
                                               Benchmark=case_when(str_detect(Benchmark,"GLORI_HEK293") ~ "GLORI_HEK293", str_detect(Benchmark,"GLORI_mESC") ~ "GLORI_mESC")) %>%
  group_by(Cell,Benchmark,cutoff) %>% mutate(avg.FPRsurrogate=mean(FPR_site), avg.TPRsurrogate=mean(TPR_site)) %>% as.data.table() %>% distinct(Cell,Benchmark,cutoff,avg.FPRsurrogate,avg.TPRsurrogate)

dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 <- dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 %>% mutate(cutoff=as.character(cutoff)) %>%
  pivot_longer(cols = c("avg.FPRsurrogate", "avg.TPRsurrogate"), names_to = "PerformanceMetric", values_to = "Value") %>% as.data.table() %>% mutate(cutoff=as.numeric(cutoff)) %>%
  mutate(PerformanceMetric=factor(PerformanceMetric, levels=c("avg.FPRsurrogate","avg.TPRsurrogate"), labels=c("surrogateFPR = % of peak without benchmark sites", "surrogateTPR = % of total benchmark sites in peaks"))) %>%
  mutate(Cell=paste0("NEB_",Cell))

p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 <-  ggplot(data=dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2,aes(x=cutoff,y=Value,color=Cell))+
  # geom_smooth(aes(linetype=PerformanceMetric),linewidth=0.5,se = F)+
  geom_step(aes(linetype=PerformanceMetric),linewidth=0.5)+
  # geom_smooth(aes(linetype=PerformanceMetric),linewidth=0.5,se = F)+
  coord_cartesian(xlim = c(1,12))+
  labs(x="Fold change cutoff for exomePeak2 m6A peak",y="Performance metric value")+
  scale_color_manual(values = c("#A7D2E6","#E9909B"),breaks = c("NEB_mESC","NEB_HEK293"),guide="none")+
  scale_linetype_manual(values=c(2,1),breaks = c("surrogateFPR = % of peak without benchmark sites","surrogateTPR = % of total benchmark sites in peaks"))+
  guides(linetype=guide_legend(nrow = 2,override.aes = list(color = "black", size = 0.6)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
# legend-only plot for color: hide linetype legend
p_col_legend <- ggplot(data=dt.FPR.TPR.changeby.foldchange.cutoff.exomePeak2,aes(x=cutoff,y=Value,color=Cell)) +
  geom_line() +
  scale_color_manual(values = c("#A7D2E6","#E9909B"),breaks = c("NEB_mESC","NEB_HEK293"))+
  theme_minimal() +
  guides(linetype = "none") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "right",
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
p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2 <- cowplot::ggdraw() +
  cowplot::draw_plot(p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2, 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 0.98, y = 0.55, width = 0.22, height = 0.35, hjust = 1, vjust = 0.5)

#2b The performance of six method default parameter combo in NEB HEK293 and mESC using GLORI and eTAMseq benchmark sites
dt.pruned.FPR.TPR.AUC.NEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HEK.all.combo.RDS")
dt.pruned.FPR.TPR.AUC.NEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.mESC.all.combo.RDS")
# dt.pruned.FPR.TPR.AUC.SYSY.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.all.combo.RDS")
# dt.pruned.FPR.TPR.AUC.Abcam.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.all.combo.RDS")
# dt.pruned.FPR.TPR.AUC.Milipore.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.Millipore.mESC.all.combo.RDS")

dt.method.defaultcombo.AUC.NEB <- rbind(dt.pruned.FPR.TPR.AUC.NEB.HEK %>% inner_join(x=.,y=dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE) %>% dplyr::select(Method,ComboName,Combo,ComboID,IsDefaultCombo),by=c("Method","ComboName")),
                                        dt.pruned.FPR.TPR.AUC.NEB.mESC %>% inner_join(x=.,y=dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE) %>% dplyr::select(Method,ComboName,Combo,ComboID,IsDefaultCombo),by=c("Method","ComboName")))
dt.method.defaultcombo.AUC.NEB %>% distinct(Method,ComboName,IsDefaultCombo,Benchmark)
#        Method ComboName IsDefaultCombo                          Benchmark
# 1:      MACS2   Combo15           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 2: MeRIPtools   Combo22           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 3:    MeTPeak    Combo6           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 4:      TRESS   Combo19           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 5:  exomePeak   Combo11           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 6: exomePeak2    Combo2           TRUE dt.benchmark.GLORI_HEK293.m6As.bed
# 7:      MACS2   Combo15           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
# 8: MeRIPtools   Combo22           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
# 9:    MeTPeak    Combo6           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
# 10:      TRESS   Combo19           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
# 11:  exomePeak   Combo11           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
# 12: exomePeak2    Combo2           TRUE   dt.benchmark.GLORI_mESC.m6As.bed
#highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
dt.method.defaultcombo.AUC.NEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::filter(abs(cutoff-2)<=1) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
dt.method.defaultcombo.AUC.NEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
dt.method.defaultcombo.AUC.NEB <- dt.method.defaultcombo.AUC.NEB %>% mutate(Cell = case_when(str_detect(Benchmark,"HEK")~ "HEK293", .default = "mESC")) %>% mutate(Cell=factor(Cell,levels=c("HEK293","mESC"))) %>%
  dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))

#add the ScaledAUC and TPR20 for each method default combo
dt.ScaledAUC.TPR20.Method.DefaultCombo.NEB <- dt.method.defaultcombo.AUC.NEB %>% distinct(Method,ComboName,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
                          mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.method.defaultcombo.AUC.NEB$TPR),2)) %>% mutate(xpos=0.7)
p.GLORI.AUC.Method.DefaultCombo.NEB.HEK.mESC <-  ggplot(data=dt.method.defaultcombo.AUC.NEB,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of m6A peak from default combo",y="Surrogate TPR of m6A peak from default combo")+
  geom_text(data=dt.ScaledAUC.TPR20.Method.DefaultCombo.NEB, aes(x=xpos,y=ypos,color=Method,label=label),size=1.6,fontface="bold")+
  scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                     labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                 as.character(x), ""))+
  coord_cartesian(ylim = c(0,0.6))+
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

#2c the improved performance by parameter refinement
dt.GLORI.AUC.NEB <- rbind(dt.pruned.FPR.TPR.AUC.NEB.HEK %>% distinct(Method,ComboName,Benchmark,AUC.scaled, TPR20),
                           dt.pruned.FPR.TPR.AUC.NEB.mESC %>% distinct(Method,ComboName,Benchmark,AUC.scaled, TPR20)) %>%
                    left_join(x=.,y=dt.method.combo %>% dplyr::select(Method,ComboName,IsDefaultCombo))

dt.BestCombo.NEB <- dt.GLORI.AUC.NEB %>% filter(str_detect(Benchmark,"HEK293")) %>% dplyr::arrange(Benchmark,Method) %>% group_by(Benchmark,Method) %>% dplyr::arrange(desc(AUC.scaled), desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.GLORI.AUC.NEB <- dt.GLORI.AUC.NEB %>% left_join(x=.,y=dt.BestCombo.NEB %>% dplyr::select(Method,ComboName) %>% mutate(IsBestCombo=TRUE),by=c("Method","ComboName"))
dt.GLORI.AUC.NEB <- dt.GLORI.AUC.NEB %>% mutate(IsBestCombo=case_when(is.na(IsBestCombo) ~ FALSE, .default = IsBestCombo))
dt.GLORI.AUC.NEB %>% filter(IsDefaultCombo==TRUE | IsBestCombo==TRUE) %>% dplyr::arrange(Benchmark,Method,IsDefaultCombo,IsBestCombo)
dt.GLORI.AUC.NEB <- dt.GLORI.AUC.NEB %>% mutate(ComboGroup=case_when(IsBestCombo==TRUE ~ "Best", IsDefaultCombo==TRUE ~ "Default", .default = "Others"))
dt.GLORI.AUC.NEB <- dt.GLORI.AUC.NEB %>% mutate(Cell=case_when(str_detect(Benchmark,"HEK") ~ "HEK293", .default = "mESC"))
dt.GLORI.AUC.NEB <- dt.GLORI.AUC.NEB %>% mutate(Method=factor(Method,levels=rev(c("MeRIPtools","MACS2","TRESS","MeTPeak","exomePeak","exomePeak2"))))
p.allcombo.GLORI.AUC.NEB <- ggplot(dt.GLORI.AUC.NEB, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.GLORI.AUC.NEB %>% dplyr::filter(ComboGroup =="Others"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=1, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.GLORI.AUC.NEB %>% dplyr::filter(ComboGroup =="Default"),aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=2)+
  geom_jitter(data = dt.GLORI.AUC.NEB %>% dplyr::filter(ComboGroup =="Best"),aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=4)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled)") +
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
p_col_legend <- ggplot(data=dt.GLORI.AUC.NEB,aes(x = Method, y = AUC.scaled)) +
  geom_jitter(aes(shape=ComboGroup,color=ComboGroup,size=ComboGroup), width = 0.18, height = 0,  stroke = 0.6) +
  scale_shape_manual(values = c(18,18,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  scale_size_manual(values = c(4,2,1), breaks = c("Best","Default","Others"), name = "Combo group") +
  scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  theme_minimal() +
  guides(linetype = "none") +
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "right",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
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
p.allcombo.GLORI.AUC.NEB  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.allcombo.GLORI.AUC.NEB , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1, y = 0.75, width = 0.22, height = 0.25, hjust = 1, vjust = 0.5)

#2d the AUC curve of Default and Best Combo of Method in NEB HEK and mESC
dt.method.best.AUC.NEB <- rbind(dt.pruned.FPR.TPR.AUC.NEB.HEK %>% inner_join(x=.,y=dt.GLORI.AUC.NEB %>% dplyr::filter(IsBestCombo==TRUE) %>% dplyr::distinct(Method,ComboName,IsBestCombo),by=c("Method","ComboName")),
                                dt.pruned.FPR.TPR.AUC.NEB.mESC %>% inner_join(x=.,y=dt.GLORI.AUC.NEB %>% dplyr::filter(IsBestCombo==TRUE) %>% dplyr::distinct(Method,ComboName,IsBestCombo),by=c("Method","ComboName")))
dt.method.best.AUC.NEB %>% distinct(Method,ComboName,IsBestCombo,Benchmark)
#highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
dt.method.best.AUC.NEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::filter(abs(cutoff-2)<=1) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
dt.method.best.AUC.NEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
dt.method.best.AUC.NEB <- dt.method.best.AUC.NEB %>% mutate(Cell = case_when(str_detect(Benchmark,"HEK")~ "HEK293", .default = "mESC")) %>% mutate(Cell=factor(Cell,levels=c("HEK293","mESC"))) %>%
  dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))

#add the ScaledAUC and TPR20 for each method default combo
dt.ScaledAUC.TPR20.Method.best.NEB <- dt.method.best.AUC.NEB %>% distinct(Method,ComboName,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.method.best.AUC.NEB$TPR),2)) %>% mutate(xpos=0.7)
p.GLORI.AUC.Method.best.NEB.HEK.mESC <-  ggplot(data=dt.method.best.AUC.NEB,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of m6A peak from best combo",y="Surrogate TPR of m6A peak from best combo")+
  geom_text(data=dt.ScaledAUC.TPR20.Method.best.NEB , aes(x=xpos,y=ypos,color=Method,label=label),size=1.6,fontface="bold")+
  scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
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

#2e the added TP and removed FP by used the optimal parameter
#compare the exomePeak2_default and exomePeak2_best using same FPR<=20% cutoff
#obtain the corresponding cutoff and method combo for exomePeak2 best
dt.method.best.AUC.NEB %>% distinct(Method,ComboName,IsBestCombo,Benchmark) %>% filter(Method=="exomePeak2") %>% distinct(Method,ComboName)
dt.method.best.AUC.NEB %>% filter(Method=="exomePeak2" & IsBestCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
#obtain the corresponding cutoff and method combo for exomePeak2 default
dt.method.defaultcombo.AUC.NEB %>% dplyr::filter(Method=="exomePeak2" & IsDefaultCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
dt.selected.method.combo.cutoff.FPR20 <- rbind(dt.method.best.AUC.NEB %>% filter(Method=="exomePeak2" & IsBestCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>%
                                                 dplyr::select(Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Best"),
                                               dt.method.defaultcombo.AUC.NEB %>% dplyr::filter(Method=="exomePeak2" & IsDefaultCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>%
                                                 dplyr::select(Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Default")) %>% mutate(Combo = paste0(Method,"_",ComboName))
#load exomePeak2 (Default and Best Combo) for HEK_NEB and mESC_NEB
#HEK293
dt.exomePeak2.BestDefault.m6A.HEK293 <- rbind(LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB/", dt.selected.method.combo.cutoff.FPR20[Group=="Best",Combo])) %>% mutate(Group="Best"),
                                            LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB/", dt.selected.method.combo.cutoff.FPR20[Group=="Default",Combo])) %>% mutate(Group="Default")
                                            )

#filter according to cutoff at FPR20
dt.exomePeak2.BestDefault.m6A.HEK293 <- dt.exomePeak2.BestDefault.m6A.HEK293 %>% dplyr::filter( (Group=="Best" & score >= dt.selected.method.combo.cutoff.FPR20[Group=="Best",cutoff]) | (Group=="Default" & score >= dt.selected.method.combo.cutoff.FPR20[Group=="Default",cutoff]))
#divide BestCombo m6A as Default- and Default+ (************************ strictly, use 0% overlap as cutoff for Default- peak definition **************************)
dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap <- foreach(c = unique(dt.exomePeak2.BestDefault.m6A.HEK293$Sample), .combine='rbind')%do%{
    bed.peak <- dt.exomePeak2.BestDefault.m6A.HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
    bed.default <- dt.exomePeak2.BestDefault.m6A.HEK293 %>% dplyr::filter(Sample==c & Group=="Default") %>% dplyr::select(seqnames,start,end,name,score,strand)
    dt.overlap <- bt.intersect(a=bed.peak, b=bed.default, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
    dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
    dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(Defaultoverlap=ifelse(sum.overlap/width>=0.0001,"Default+","Default-"))
    bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,Defaultoverlap), by="name") %>% mutate(Sample=c)
}
dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap[,.N,by=.(Sample,Defaultoverlap)] %>% dplyr::arrange(Sample,Defaultoverlap)
#          Sample Defaultoverlap    N
# 1: HEK_NEB_mRNA1       Default+  7595
# 2: HEK_NEB_mRNA1       Default-  3443
# 3: HEK_NEB_mRNA2       Default+ 11378
# 4: HEK_NEB_mRNA2       Default-  3146

#check the Default- and Default+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
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
dt.exomePeak2.BestCombo.HEK293.peaks.Benchmarkoverlap <- foreach(c = unique(dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap$Sample), .combine='rbind')%do%{
  bed.peak <- dt.exomePeak2.BestDefault.m6A.HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
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
#obtain the BestCombo only m6A+ gene
dt.exomePeak2.BestCombo.HEK293.peaks.annotgene <- annot_peak(peak.bed=dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                                 strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                 annot_type = "gene")
dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap <- dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap %>% left_join(x=.,y=dt.exomePeak2.BestCombo.HEK293.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
dt.exomePeak2.BestCombo.specific.m6Agene.HEK293 <-   dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestCombo.specific=all(Defaultoverlap=="Default-")) %>% as.data.table() %>%#if all peak of one gene is Default-, then this gene is BestCombo specific m6A+ gene
   dplyr::filter(BestCombo.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
#obtain the BestCombo specific peaks overlapped with (high confident) benchmark sites
dt.exomePeak2.BestCombo.specific.peak.benchmark.HEK293 <-  dt.exomePeak2.BestCombo.specific.m6Agene.HEK293 %>% filter(Defaultoverlap=="Default-" & BestCombo.specific==TRUE) %>%
                                                  left_join(x=.,dt.exomePeak2.BestCombo.HEK293.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
                                                  dplyr::filter(nBenchSite>0)
#rank the BestCombo specific m6A+ genes with benchmark m6A sites support
dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.HEK293 <- dt.exomePeak2.BestCombo.specific.peak.benchmark.HEK293 %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
                                                             group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
  as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))



#mESC
dt.exomePeak2.BestDefault.m6A.mESC <- rbind(LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB_mESC//", dt.selected.method.combo.cutoff.FPR20[Group=="Best",Combo])) %>% mutate(Group="Best"),
                                            LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_NEB_mESC//", dt.selected.method.combo.cutoff.FPR20[Group=="Default",Combo])) %>% mutate(Group="Default")
                                            )
#filter according to cutoff at FPR20
dt.exomePeak2.BestDefault.m6A.mESC <- dt.exomePeak2.BestDefault.m6A.mESC %>% dplyr::filter( (Group=="Best" & score >= dt.selected.method.combo.cutoff.FPR20[Group=="Best",cutoff]) | (Group=="Default" & score >= dt.selected.method.combo.cutoff.FPR20[Group=="Default",cutoff]))
#divide BestCombo m6A as Default- and Default+
dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap <- foreach(c = unique(dt.exomePeak2.BestDefault.m6A.mESC$Sample), .combine='rbind')%do%{
  bed.peak <- dt.exomePeak2.BestDefault.m6A.mESC %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
  bed.default <- dt.exomePeak2.BestDefault.m6A.mESC %>% dplyr::filter(Sample==c & Group=="Default") %>% dplyr::select(seqnames,start,end,name,score,strand)
  dt.overlap <- bt.intersect(a=bed.peak, b=bed.default, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
  dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
  dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(Defaultoverlap=ifelse(sum.overlap/width>=0.0001,"Default+","Default-"))
  bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,Defaultoverlap), by="name") %>% mutate(Sample=c)
}
dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap[,.N,by=.(Sample,Defaultoverlap)] %>% dplyr::arrange(Sample,Defaultoverlap)
#       Sample Defaultoverlap    N
# 1: mESC_WT1       Default+ 7372
# 2: mESC_WT1       Default- 2990
# 3: mESC_WT2       Default+ 6897
# 4: mESC_WT2       Default- 5758

#check the Default- and Default+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
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

dt.exomePeak2.BestCombo.mESC.peaks.Benchmarkoverlap <- foreach(c = unique(dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap$Sample), .combine='rbind')%do%{
  bed.peak <- dt.exomePeak2.BestDefault.m6A.mESC %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
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
#obtain the BestCombo only m6A+ gene
dt.exomePeak2.BestCombo.mESC.peaks.annotgene <- annot_peak(peak.bed=dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                             strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.vM33.annotation.gtf",
                                                             annot_type = "gene")
dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap <- dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap %>% left_join(x=.,y=dt.exomePeak2.BestCombo.mESC.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
dt.exomePeak2.BestCombo.specific.m6Agene.mESC <-   dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
  group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestCombo.specific=all(Defaultoverlap=="Default-")) %>% as.data.table() %>%#if all peak of one gene is Default-, then this gene is BestCombo specific m6A+ gene
  dplyr::filter(BestCombo.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
#obtain the BestCombo specific peaks overlapped with (high confident) benchmark sites
dt.exomePeak2.BestCombo.specific.peak.benchmark.mESC <-  dt.exomePeak2.BestCombo.specific.m6Agene.mESC %>% filter(Defaultoverlap=="Default-" & BestCombo.specific==TRUE) %>%
  left_join(x=.,dt.exomePeak2.BestCombo.mESC.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
  dplyr::filter(nBenchSite>0)
#rank the BestCombo specific m6A+ genes with benchmark m6A sites support
dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.mESC <- dt.exomePeak2.BestCombo.specific.peak.benchmark.mESC %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
  group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
  as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))

#visualization of ratio of Default- ratio of exomePeak2 peaks in HEK293 and mESC
dt.exomePeak2.BestCombo.specific.ratio.m6A.peak <- rbind(dt.exomePeak2.BestCombo.HEK293.peaks.Defaultoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestComboOnly = n_distinct(name[Defaultoverlap=="Default-"])/n_distinct(name),
                                                                                                                                             nPeak.BestComboOnly= n_distinct(name[Defaultoverlap=="Default-"]), nPeak.BestCombo=n_distinct(name)) %>%
                                                           as.data.table() %>% distinct(Sample,Ratio.Peak.BestComboOnly,nPeak.BestComboOnly,nPeak.BestCombo),
                                                         dt.exomePeak2.BestCombo.mESC.peaks.Defaultoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestComboOnly = n_distinct(name[Defaultoverlap=="Default-"])/n_distinct(name),
                                                                                                                                           nPeak.BestComboOnly= n_distinct(name[Defaultoverlap=="Default-"]), nPeak.BestCombo=n_distinct(name)) %>%
                                                           as.data.table() %>% distinct(Sample,Ratio.Peak.BestComboOnly,nPeak.BestComboOnly,nPeak.BestCombo)
                                                        )
dt.exomePeak2.BestCombo.specific.ratio.m6A.peak <- dt.exomePeak2.BestCombo.specific.ratio.m6A.peak %>% mutate(label1=paste0(nPeak.BestComboOnly,"/",nPeak.BestCombo), label2=paste0(round(Ratio.Peak.BestComboOnly*100,1),"%")) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_NEB_mRNA1","HEK_NEB_mRNA2","mESC_WT1","mESC_WT2"), labels=c("NEB_HEK293_mRNA1","NEB_HEK293_mRNA2","NEB_mESC_mRNA1","NEB_mESC_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293",.default = "mESC"))
p.exomePeak2.BestCombo.specific.ratio.m6A.peak <-
  ggplot(dt.exomePeak2.BestCombo.specific.ratio.m6A.peak, aes(y = Sample, x = Ratio.Peak.BestComboOnly, fill = Cell)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  # percent label to the right of the bar
  geom_text(aes(label = label2, x =  Ratio.Peak.BestComboOnly + 0.02),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  #N peak label to the left of the bar
  geom_text(aes(label = label1, x = Ratio.Peak.BestComboOnly/2),
            color = "black", size = 2, fontface = "plain", hjust = 0) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
  scale_fill_manual(values = c("#7b93c6","#acd9ee"),breaks = c("HEK293","mESC")) +
  labs(x = "% of exomePeak2 best combo specific peaks", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
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
dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark <- rbind(dt.exomePeak2.BestCombo.specific.peak.benchmark.HEK293 %>% group_by(Sample) %>%
                                                                         mutate(No.BestComboOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                No.BestComboOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                         distinct(Sample,No.BestComboOnly.m6AGene.BenchmarkSite,No.BestComboOnly.m6AGene.ConfidentBenchmarkSite),
                                                                       dt.exomePeak2.BestCombo.specific.peak.benchmark.mESC %>% group_by(Sample) %>%
                                                                         mutate(No.BestComboOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                No.BestComboOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                         distinct(Sample,No.BestComboOnly.m6AGene.BenchmarkSite,No.BestComboOnly.m6AGene.ConfidentBenchmarkSite)
                                                                       )
dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark <- dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark %>%
  pivot_longer(cols = c("No.BestComboOnly.m6AGene.BenchmarkSite", "No.BestComboOnly.m6AGene.ConfidentBenchmarkSite"),values_to = "No.BestComboOnly.m6AGene", names_to = "GeneGroup",names_prefix = "No.BestComboOnly.m6AGene.") %>%
  as.data.table() %>% mutate(GeneGroup=factor(GeneGroup, levels=c("BenchmarkSite","ConfidentBenchmarkSite"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>%
  mutate(Sample=factor(Sample, levels=c("HEK_NEB_mRNA1","HEK_NEB_mRNA2","mESC_WT1","mESC_WT2"), labels=c("NEB_HEK293_mRNA1","NEB_HEK293_mRNA2","NEB_mESC_mRNA1","NEB_mESC_mRNA2"))) %>%
  mutate(Cell=case_when(str_detect(Sample,"HEK")~"HEK293",.default = "mESC"))

dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark <- dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark %>% mutate(label.count=paste0(paste0("No.Gene=",No.BestComboOnly.m6AGene)))
p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark  <- ggplot(data=dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark, aes(x=Sample,y=No.BestComboOnly.m6AGene, fill=GeneGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  # percent label to the left  bar
  geom_text(data = dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark %>% dplyr::filter(GeneGroup=="Benchmark m6As"),
            aes(label = label.count, y =  No.BestComboOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
  # percent label to the right  bar
  geom_text(data = dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark %>% dplyr::filter(GeneGroup=="Confident benchmark m6As"),
            aes(label = label.count, y =  No.BestComboOnly.m6AGene/2+0.02, x=Sample),
            color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
  scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
  labs(x=NULL,y=str_wrap("N of genes harbor only exomePeak2 best combo specific peaks",width=40))+
  guides(fill=guide_legend(title = "m6A peaks contain"))+
  scale_y_continuous(expand = expansion(add=c(0,0.05)))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "inside",legend.position.inside = c(0.5,0.9),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1.15)),
    plot.subtitle = element_text(size = rel(0.95)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark

#visualization of IGV plot of exomePeak2 BestCombo specific m6A+ gene in HEK293 and mESC

dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

selected.exomePeak2.BestComboOnly.m6Agene.NEB.HEK <-  dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.HEK293  %>% dplyr::filter(nRep==2) %>%
  distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=5)
dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==c("GSK3A","SOCS5","PPM1B")[1])
StudySamples <- c("HEK_NEB_mRNA1", "HEK_NEB_mRNA2")
      ## Create a page (7.5*7.5cm)
      window.size=5#100kb
      pseudocount=1
      track.height=0.5
      pdf("Fig2_NEB_HEK293_GSK3A_exomePeak2_bestcombo_specific_peak.pdf")

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

      # #add exomePeak2 default combo peak in this region
      for(s in 1:length(StudySamples)){
        dt.peak.region <- dt.exomePeak2.BestDefault.m6A.HEK293 %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
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
      legendPlot.exomePeak2.peak <- plotLegend(legend = c("BestCombo","DefaultCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 2.5, height = 1,
                                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                               title = expression("exomePeak2 peak"),fontface="plain")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")
      dev.off()

      #mESC exomePeak2 best combo specific m6A+ gene with high confident benchmark m6As
      #select well-studied genes
      selected.exomePeak2.BestComboOnly.m6Agene.NEB.mESC <-  dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.mESC  %>% dplyr::filter(nRep==2) %>%
        distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"Gm")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=5)
      dt.selected.gene <- dt.gene.mm39 %>% dplyr::filter(gene_name==c("Lrp2","Foxm1","Acot1")[2])
      StudySamples <- c("mESC_WT1","mESC_WT2")
      library(BSgenome.Mmusculus.UCSC.mm39)

      ## Create a page (7.5*7.5cm)
      window.size=10#100kb
      pseudocount=1
      track.height=0.5
      pdf("Fig2_NEB_mESC_Foxm1_exomePeak2_bestcombo_specific_peak.pdf")
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
      # #add exomePeak2 default combo peak in this region
      for(s in 1:length(StudySamples)){
        dt.peak.region <- dt.exomePeak2.BestDefault.m6A.mESC %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
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
      legendPlot.exomePeak2.peak <- plotLegend(legend = c("BestCombo","DefaultCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+2)*track.height, width = 2.5, height = 1,
                                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                               title = expression("exomePeak2 peak"),fontface="plain")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+3)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

      dev.off()

#save intermediate variables for re-use
save.image("sm6APeak_Figure2.intermediate.results.RDS")

 ####################### combine all figure 2 together #############################
      fig2.list <- list(p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2=p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2,
                        p.GLORI.AUC.Method.DefaultCombo.NEB.HEK.mESC=p.GLORI.AUC.Method.DefaultCombo.NEB.HEK.mESC,
                        p.allcombo.GLORI.AUC.NEB=p.allcombo.GLORI.AUC.NEB,
                        p.GLORI.AUC.Method.best.NEB.HEK.mESC=p.GLORI.AUC.Method.best.NEB.HEK.mESC,
                        p.exomePeak2.BestCombo.specific.ratio.m6A.peak=p.exomePeak2.BestCombo.specific.ratio.m6A.peak,
                        p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark=p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark)
      saveRDS(fig2.list, file="Figure2.plot.list.RDS")

      pdf("Figure2.Parameter_refinement_could_improve_performance_of_m6A_calling.pdf",width = 8.2677, height = 11.693)
      pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
      plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.FPR.TPR.changeby.foldchange.cutoff.exomePeak2, x = 0.05, y=0.2, default.units = "cm",width = 5.8, height = 5.8)
      plotText(label = "B", fontsize = 8, fontface = "bold",x = 6+2.3, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.GLORI.AUC.Method.DefaultCombo.NEB.HEK.mESC, x = 6+2.3, y=0.2, default.units = "cm",width = 9.8, height = 5.8)
      plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.allcombo.GLORI.AUC.NEB, x = 0.05, y=6.2, default.units = "cm",width = 8.0, height = 5.8)
      plotText(label = "D", fontsize = 8, fontface = "bold",x = 6+2.3, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.GLORI.AUC.Method.best.NEB.HEK.mESC, x = 6+2.3, y=6.2, default.units = "cm",width = 9.8, height = 5.8)
      plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.exomePeak2.BestCombo.specific.ratio.m6A.peak, x = 0.05, y=12.2, default.units = "cm",width = 7.5, height = 5.8)
      plotText(label = "G", fontsize = 8, fontface = "bold",x = 8.5, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)

      plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05, y = 18.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = fig2.list$p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark, x = 0.05, y=18.2, default.units = "cm",width = 7.5, height = 5.8)
      plotText(label = "H", fontsize = 8, fontface = "bold",x = 8.5, y = 18.2,just = c("top","left"), default.units = "cm",draw=T)


      dev.off()

### Part2S The parameter refinement significantly improved the performance of m6A calling algorithm for MeRIP-seq data (non NEB anti-m6A, SYSY, Abcam, and Milipore) ###############
      dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
      dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE)
      #S2A The performance of six method default parameter combo in SYSY_HEK293, Abcam_HEK293
      # dt.pruned.FPR.TPR.AUC.nonNEB.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.nonNEB.HEK.all.combo.RDS")
      # dt.pruned.FPR.TPR.AUC.nonNEB.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.nonNEB.mESC.all.combo.RDS")
      dt.pruned.FPR.TPR.AUC.SYSY.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.SYSY.HEK.all.combo.RDS") %>% mutate(Cell="SYSY_HEK293")
      dt.pruned.FPR.TPR.AUC.Abcam.HEK <- readRDS("dt.pruned.FPR.TPR.AUC.Abcam.HEK.all.combo.RDS") %>% mutate(Cell="Abcam_HEK293")
      # dt.pruned.FPR.TPR.AUC.Milipore.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.Millipore.mESC.all.combo.RDS")

      dt.method.defaultcombo.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.AUC.SYSY.HEK %>% mutate(Cell="SYSY_HEK293") %>% inner_join(x=.,y=dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE) %>% dplyr::select(Method,ComboName,Combo,ComboID,IsDefaultCombo),by=c("Method","ComboName")),
                                                 dt.pruned.FPR.TPR.AUC.Abcam.HEK %>% mutate(Cell="Abcam_HEK293") %>% inner_join(x=.,y=dt.method.combo %>% dplyr::filter(IsDefaultCombo==TRUE) %>% dplyr::select(Method,ComboName,Combo,ComboID,IsDefaultCombo),by=c("Method","ComboName")))
      dt.method.defaultcombo.AUC.nonNEB %>% distinct(Cell,Method,ComboName,IsDefaultCombo,Benchmark)
      #highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
      dt.method.defaultcombo.AUC.nonNEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::filter(abs(cutoff-2)<=1) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
      dt.method.defaultcombo.AUC.nonNEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
      dt.method.defaultcombo.AUC.nonNEB <- dt.method.defaultcombo.AUC.nonNEB %>% mutate(Cell=factor(Cell,levels=c("SYSY_HEK293","Abcam_HEK293"))) %>%
        dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))

      #add the ScaledAUC and TPR20 for each method default combo
      dt.ScaledAUC.TPR20.Method.DefaultCombo.NEB <- dt.method.defaultcombo.AUC.nonNEB %>% distinct(Method,ComboName,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
        mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.method.defaultcombo.AUC.nonNEB$TPR),2)) %>% mutate(xpos=0.7)

      p.GLORI.AUC.Method.DefaultCombo.nonNEB <-  ggplot(data=dt.method.defaultcombo.AUC.nonNEB,aes(x=FPR,y=TPR,color=Method))+
        geom_step(linewidth=0.5)+
        labs(x="Surrogate FPR of m6A peak from default combo",y="Surrogate TPR of m6A peak from default combo")+
        geom_text(data=dt.ScaledAUC.TPR20.Method.DefaultCombo.NEB, aes(x=xpos,y=ypos,color=Method,label=label),size=1.6,fontface="bold")+
        scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
        scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
                           labels = function(x) ifelse(x %in% seq(0.2, 1, by = 0.2),    # show labels only for these
                                                       as.character(x), ""))+
        coord_cartesian(ylim = c(0,0.6))+
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
      p.GLORI.AUC.Method.DefaultCombo.nonNEB

      #2c the improved performance by parameter refinement
      dt.GLORI.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.AUC.SYSY.HEK %>% mutate(Cell="SYSY_HEK293") %>% distinct(Cell,Method,ComboName,Benchmark,AUC.scaled,TPR20),
                                   dt.pruned.FPR.TPR.AUC.Abcam.HEK%>% mutate(Cell="Abcam_HEK293") %>% distinct(Cell,Method,ComboName,Benchmark,AUC.scaled, TPR20)) %>%
        left_join(x=.,y=dt.method.combo %>% dplyr::select(Method,ComboName,IsDefaultCombo))
      dt.BestCombo.nonNEB <- dt.GLORI.AUC.nonNEB %>% filter(str_detect(Benchmark,"GLORI")) %>% dplyr::arrange(Cell,Benchmark,Method) %>% group_by(Cell,Benchmark,Method) %>% dplyr::arrange(desc(AUC.scaled), desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
      dt.GLORI.AUC.nonNEB <- dt.GLORI.AUC.nonNEB %>% left_join(x=.,y=dt.BestCombo.nonNEB %>% dplyr::select(Cell,Method,ComboName,Benchmark) %>% mutate(IsBestCombo=TRUE),by=c("Cell","Method","ComboName","Benchmark"))
      dt.GLORI.AUC.nonNEB <- dt.GLORI.AUC.nonNEB %>% mutate(IsBestCombo=case_when(is.na(IsBestCombo) ~ FALSE, .default = IsBestCombo))
      dt.GLORI.AUC.nonNEB %>% filter(IsDefaultCombo==TRUE | IsBestCombo==TRUE) %>% dplyr::arrange(Benchmark,Method,IsDefaultCombo,IsBestCombo)
      dt.GLORI.AUC.nonNEB <- dt.GLORI.AUC.nonNEB %>% mutate(ComboGroup=case_when(IsBestCombo==TRUE ~ "Best", IsDefaultCombo==TRUE ~ "Default", .default = "Others"))
      dt.GLORI.AUC.nonNEB <- dt.GLORI.AUC.nonNEB %>% mutate(Method=factor(Method,levels=rev(c("MeRIPtools","MACS2","TRESS","MeTPeak","exomePeak","exomePeak2"))))

      p.allcombo.GLORI.AUC.nonNEB <- ggplot(dt.GLORI.AUC.nonNEB, aes(x = Method, y = AUC.scaled)) +
        # boxplots without showing outliers (we plot points separately)
        geom_boxplot(width = 0.55,
                     outlier.shape = NA,
                     fill = "transparent",
                     colour = "grey5",
                     size = 0.4) +
        # individual points: shape by ComboGroup, fill by ComboGroup, black border
        geom_jitter(data = dt.GLORI.AUC.nonNEB %>% dplyr::filter(ComboGroup =="Others"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=1, width = 0.18, height = 0,  stroke = 0.6) +
        geom_jitter(data = dt.GLORI.AUC.nonNEB %>% dplyr::filter(ComboGroup =="Default"),aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=2)+
        geom_jitter(data = dt.GLORI.AUC.nonNEB %>% dplyr::filter(ComboGroup =="Best"),aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=4)+
        # show mean as a black diamond
        # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
        # scales for shapes and fills
        # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
        # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
        # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
        # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
        facet_wrap(~Cell)+
        # labels
        labs(x = NULL, y = "AUC (scaled)") +
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
      p_col_legend <- ggplot(data=dt.GLORI.AUC.nonNEB,aes(x = Method, y = AUC.scaled)) +
        geom_jitter(aes(shape=ComboGroup,color=ComboGroup,size=ComboGroup), width = 0.18, height = 0,  stroke = 0.6) +
        scale_shape_manual(values = c(18,18,22),breaks=c("Best","Default","Others"), name = "Combo group") +
        scale_size_manual(values = c(4,2,1), breaks = c("Best","Default","Others"), name = "Combo group") +
        scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
        theme_minimal() +
        guides(linetype = "none") +
        theme_minimal(base_size = 6,base_family = "Helvetica") +
        theme(
          legend.position = "right",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
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
      p.allcombo.GLORI.AUC.nonNEB  <- cowplot::ggdraw() +
        cowplot::draw_plot(p.allcombo.GLORI.AUC.nonNEB , 0, 0, 1, 1) +
        # place color legend (right)
        cowplot::draw_grob(leg_col, x = 1, y = 0.75, width = 0.22, height = 0.25, hjust = 1, vjust = 0.5)
      p.allcombo.GLORI.AUC.nonNEB

      #S2C the AUC curve of Default and Best Combo of Method in SYSY_HEK and Millipore_mESC
      dt.method.best.AUC.nonNEB <- rbind(dt.pruned.FPR.TPR.AUC.SYSY.HEK %>% inner_join(x=.,y=dt.GLORI.AUC.nonNEB %>% dplyr::filter(IsBestCombo==TRUE) %>% dplyr::distinct(Cell,Benchmark,Method,ComboName,IsBestCombo),by=c("Cell","Method","ComboName","Benchmark")),
                                         dt.pruned.FPR.TPR.AUC.Abcam.HEK %>% inner_join(x=.,y=dt.GLORI.AUC.nonNEB %>% dplyr::filter(IsBestCombo==TRUE) %>% dplyr::distinct(Cell,Benchmark,Method,ComboName,IsBestCombo),by=c("Cell","Method","ComboName","Benchmark")))
      dt.method.best.AUC.nonNEB %>% distinct(Cell,Method,ComboName,IsBestCombo,Benchmark)
      #highlight the FPR/TPR at cutoff~2, and the FPR/TPR use the optimal maxTPRmFPR for six default combo
      dt.method.best.AUC.nonNEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::filter(abs(cutoff-2)<=1) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
      dt.method.best.AUC.nonNEB %>% group_by(Method,ComboName,Benchmark) %>% dplyr::arrange(desc(TPR-FPR)) %>% slice_head(n=1) %>% dplyr::arrange(Benchmark,Method,ComboName)
      dt.method.best.AUC.nonNEB <- dt.method.best.AUC.nonNEB %>%
        dplyr::arrange(desc(Cell),desc(AUC.scaled)) %>% mutate(label=paste0(Method,":ScaledAUC=",AUC.scaled)) %>%  mutate(Method=factor(Method,levels=unique(Method)))
      #add the ScaledAUC and TPR20 for each method default combo
      dt.ScaledAUC.TPR20.Method.best.nonNEB <- dt.method.best.AUC.nonNEB %>% distinct(Method,ComboName,Cell,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
        mutate(ypos=rep(tail(seq(0.6,0,by= -0.6/6),6)*max(dt.method.best.AUC.nonNEB$TPR),2)) %>% mutate(xpos=0.7)
      p.GLORI.AUC.Method.best.nonNEB <-  ggplot(data=dt.method.best.AUC.nonNEB,aes(x=FPR,y=TPR,color=Method))+
        geom_step(linewidth=0.5)+
        labs(x="Surrogate FPR of m6A peak from best combo",y="Surrogate TPR of m6A peak from best combo")+
        geom_text(data=dt.ScaledAUC.TPR20.Method.best.nonNEB , aes(x=xpos,y=ypos,color=Method,label=label),size=1.6,fontface="bold")+
        scale_color_manual(values = c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0"),breaks = c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools"),guide="none")+
        scale_x_continuous(breaks = seq(0, 1, by = 0.1),                      # grid & ticks at every 0.2
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
      p.GLORI.AUC.Method.best.nonNEB


      #S2D the optimal intensity cutoff is not equal to arbitrary cutoff like 2
      #three antibody (SYSY_HEK293, Abcam_HEK293) x two method (exomePeak2, MeTPeak,exomePeak,TRESS)
      dt.select.bestcombo <- dt.method.best.AUC.nonNEB %>% distinct(Cell,Method,ComboName,IsBestCombo) %>% filter(Method %in% c("exomePeak2","exomePeak","MeTPeak","TRESS"))
      #            Cell     Method ComboName IsBestCombo
      # 1:  SYSY_HEK293 exomePeak2   Combo15        TRUE
      # 2:  SYSY_HEK293  exomePeak   Combo26        TRUE
      # 3:  SYSY_HEK293    MeTPeak   Combo18        TRUE
      # 4:  SYSY_HEK293      TRESS   Combo27        TRUE
      # 5: Abcam_HEK293 exomePeak2   Combo15        TRUE
      # 6: Abcam_HEK293    MeTPeak    Combo2        TRUE
      # 7: Abcam_HEK293      TRESS   Combo27        TRUE
      # 8: Abcam_HEK293  exomePeak   Combo18        TRUE
      #obtain the arbitrary cutoff at 2x
      #load exomePeak2 and MeTPeak best combo peak
      dt.method.bestcombo.m6A <- foreach(i = 1:nrow(dt.select.bestcombo),.combine = 'rbind')%do%{
        if(str_detect(dt.select.bestcombo$Cell[i],"SYSY")){
          LoadPeakMethod(Method = dt.select.bestcombo$Method[i],Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_SYSY/",dt.select.bestcombo$Method[i],"_",dt.select.bestcombo$ComboName[i],"/"))
        }else{
          if(str_detect(dt.select.bestcombo$Cell[i],"Abcam")){
            LoadPeakMethod(Method = dt.select.bestcombo$Method[i],Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_Abcam/",dt.select.bestcombo$Method[i],"_",dt.select.bestcombo$ComboName[i],"/"))
          }
        }
      }
      #filter at cutoff 2x
      dt.method.bestcombo.m6A.artitrary <- dt.method.bestcombo.m6A %>% filter(!str_detect(Sample,"IVT")) %>% filter(score>=2)
      dt.method.bestcombo.m6A.artitrary[,.N,by=.(Sample,Method)]
      #            Sample     Method     N
      # 1:  HEK_SYSY_mRNA1 exomePeak2 11736
      # 2:  HEK_SYSY_mRNA2 exomePeak2 11381
      # 3:  HEK_SYSY_mRNA1  exomePeak 37595
      # 4:  HEK_SYSY_mRNA2  exomePeak 38500
      # 5:  HEK_SYSY_mRNA1    MeTPeak 24307
      # 6:  HEK_SYSY_mRNA2    MeTPeak 24140
      # 7:  HEK_SYSY_mRNA1      TRESS    83
      # 8:  HEK_SYSY_mRNA2      TRESS    47
      # 9: HEK_Abcam_mRNA1 exomePeak2  4786
      # 10: HEK_Abcam_mRNA2 exomePeak2  5667
      # 11: HEK_Abcam_mRNA1    MeTPeak 15413
      # 12: HEK_Abcam_mRNA2    MeTPeak 15552
      # 13: HEK_Abcam_mRNA1  exomePeak 13346
      # 14: HEK_Abcam_mRNA2  exomePeak 13898
      #calculate TPR and FPR for 2x m6a
      dt.Benchmark.GLORI.HEK293 <- fread("/data/m6A_calling_strategy/RIPPeakS_resource/dt.benchmark.GLORI_HEK293.m6As.bed")
      dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary <- foreach(S = unique(dt.method.bestcombo.m6A.artitrary$Sample), .combine='rbind')%do%{
        foreach(M = unique(dt.method.bestcombo.m6A.artitrary[Sample==S,Method]), .combine='rbind')%do%{
          dt.m6a.input <- dt.method.bestcombo.m6A %>% filter(Sample==S & Method==M & score>=2)
          if(nrow(dt.m6a.input)>0){
            dt.TPR.FPR <- EstimatePerformance(InputBed = dt.m6a.input %>% dplyr::select(seqnames,start,end,name,score,strand),
                                              PositiveBed = dt.Benchmark.GLORI.HEK293[,1:6], strandness = "s")
            dt.TPR.FPR %>% head(n=1) %>% mutate(Sample=S, Method=M) %>% mutate(Cutoff="2X",FPR=FalsePositiveRate, TPR=TruePositiveRate) %>%
              dplyr::select(Sample,Method,Cutoff,FPR,TPR)
          }else{
            data.table(Sample=S,Method=M,Cutoff="2X",FPR=0,TPR=0)
          }
        }
      }
      dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary <- rbind(dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary,data.table(Sample=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2"),Method="TRESS",Cutoff="2X",FPR=0,TPR=0))
      dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method)

      #obtain the optimal cutoff
      dt.method.best.AUC.nonNEB.optimalcutoff <- dt.method.best.AUC.nonNEB %>% group_by(Cell,label) %>% filter(Method %in% c("exomePeak2","MeTPeak","exomePeak","TRESS")) %>%
        filter(abs(TPR-FPR-maxTPRmFPR)<=0.02 & !is.na(cutoff)) %>% dplyr::arrange(desc(TPR)) %>% slice_head(n=1) %>% as.data.table()
      dt.method.best.AUC.nonNEB.optimalcutoff
      #     FPR  TPR cutoff     Method ComboName                          Benchmark  AUC AUC.scaled maxTPRmFPR TPR20         Cell IsBestCombo                     label
      # 1: 0.20 0.21   4.17    MeTPeak    Combo2 dt.benchmark.GLORI_HEK293.m6As.bed 0.09       0.18       0.01  0.21 Abcam_HEK293        TRUE    MeTPeak:ScaledAUC=0.18
      # 2: 0.21 0.40   1.57 exomePeak2   Combo15 dt.benchmark.GLORI_HEK293.m6As.bed 0.18       0.36       0.19  0.38 Abcam_HEK293        TRUE exomePeak2:ScaledAUC=0.36
      # 3: 0.25 0.42   5.14    MeTPeak   Combo18 dt.benchmark.GLORI_HEK293.m6As.bed 0.16       0.32       0.17  0.37  SYSY_HEK293        TRUE    MeTPeak:ScaledAUC=0.32
      # 4: 0.22 0.56   1.81 exomePeak2   Combo15 dt.benchmark.GLORI_HEK293.m6As.bed 0.25       0.50       0.34  0.53  SYSY_HEK293        TRUE  exomePeak2:ScaledAUC=0.5
      #calculate TPR and FPR for optimal cutted m6A
      dt.Benchmark.GLORI.HEK293 <- fread("/data/m6A_calling_strategy/RIPPeakS_resource/dt.benchmark.GLORI_HEK293.m6As.bed")
      dt.pruned.TPR.FPR.method.bestcombo.m6A.optimal <- foreach(S = unique(dt.method.bestcombo.m6A$Sample), .combine='rbind')%do%{
        foreach(M = unique(dt.method.bestcombo.m6A[Sample==S,Method]), .combine='rbind')%do%{
          #determine the cutoff
          if(str_detect(S,"SYSY")){optimal.cutoff <- dt.method.best.AUC.nonNEB.optimalcutoff[Method==M & Cell=="SYSY_HEK293",cutoff]}else{
            if(str_detect(S,"Abcam")){optimal.cutoff <- dt.method.best.AUC.nonNEB.optimalcutoff[Method==M & Cell=="Abcam_HEK293",cutoff]}
          }
          dt.m6a.input <- dt.method.bestcombo.m6A %>% filter(Sample==S & Method==M & score>=optimal.cutoff)
          dt.TPR.FPR <- EstimatePerformance(InputBed = dt.m6a.input %>% dplyr::select(seqnames,start,end,name,score,strand),
                                            PositiveBed = dt.Benchmark.GLORI.HEK293[,1:6], strandness = "s")
          dt.TPR.FPR %>% head(n=1) %>% mutate(Sample=S, Method=M) %>% mutate(Cutoff=paste0(optimal.cutoff,"X"),FPR=FalsePositiveRate, TPR=TruePositiveRate) %>%
            dplyr::select(Sample,Method,Cutoff,FPR,TPR)
        }
      }
      dt.pruned.TPR.FPR.method.bestcombo.m6A.optimal %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method)
      rbind(dt.pruned.TPR.FPR.method.bestcombo.m6A.optimal %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method),
            dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method)) %>% dplyr::arrange(Sample,Method,Cutoff=="2X") %>%
        mutate(TPRmFPR=round(TPR-FPR,2),FPR=round(FPR,2), TPR=round(TPR,2),Cell=case_when(str_detect(Sample,"SYSY") ~ "SYSY_HEK293",str_detect(Sample,"Abcam") ~ "Abcam_HEK293")) %>%
        group_by(Cell,Method,Cutoff) %>% mutate(FPR=mean(FPR),TPR=mean(TPR),TPRmFPR=mean(TPRmFPR)) %>% as.data.table() %>% dplyr::select(-Sample) %>% distinct() %>%
        dplyr::arrange(Cell,Method,Cutoff=="2X")

      #Abcam subplot
      dt.TPR.FPR.bestcombo.arbitraryvsoptimal <- rbind(dt.pruned.TPR.FPR.method.bestcombo.m6A.optimal %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method),
                                                       dt.pruned.TPR.FPR.method.bestcombo.m6A.artitrary %>% dplyr::arrange(str_detect(Sample,"Abcam"),Method)) %>% dplyr::arrange(Sample,Method,Cutoff=="2X") %>%
        mutate(TPRmFPR=round(TPR-FPR,2),FPR=round(FPR,2), TPR=round(TPR,2),Cell=case_when(str_detect(Sample,"SYSY") ~ "SYSY_HEK293",str_detect(Sample,"Abcam") ~ "Abcam_HEK293")) %>%
        group_by(Cell,Method,Cutoff) %>% mutate(FPR=mean(FPR),TPR=mean(TPR),TPRmFPR=mean(TPRmFPR)) %>% as.data.table() %>% dplyr::select(-Sample) %>% distinct() %>%
        dplyr::arrange(Cell,Method,Cutoff=="2X")  %>%  mutate(Group=factor(Cutoff=="2X", levels=c(TRUE,FALSE), labels=c("Arbitrary2X","OptimalCutoff"))) %>%
        pivot_longer(cols = c("FPR","TPR"),names_to = "Metric",values_to = "Value") %>% as.data.table()

      p.TPR.FPR.bestcombo.arbitraryvsoptimal <- ggplot(data=dt.TPR.FPR.bestcombo.arbitraryvsoptimal ,aes(x=Metric,y=Value,fill=Group))+
        geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
        # # percent label to the left  bar
        # geom_text(data = dt.TPR.FPR.bestcombo.arbitraryvsoptimal %>% dplyr::filter(Group=="Arbitrary2X"),
        #           aes(label = Cutoff, y = Value+0.02, x=Metric),
        #           color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = -0.15) +
        # percent label to the right  bar
        geom_text(data =  dt.TPR.FPR.bestcombo.arbitraryvsoptimal %>% dplyr::filter(Group=="OptimalCutoff"),
                  aes(label = Cutoff, y = Value+0.02, x=Metric),
                  color = "black", size = 1.8, fontface = "plain",angle=0,nudge_x = 0.15) +
        scale_fill_manual(values=c("#efb091","#eb4601"),breaks=c("Arbitrary2X","OptimalCutoff"))+
        labs(x=NULL,y=str_wrap("Value of surrogate FPR/TPR ",width=40))+
        guides(fill=guide_legend(title = "Peak"))+
        scale_y_continuous(expand = expansion(add=c(0,0.05)))+
        # facet_wrap(~Cell+Method,nrow = 2)+
        facet_grid(Cell~Method)+
        theme_minimal(base_size = 6,base_family = "Helvetica") +
        theme(
          strip.background = element_rect(fill="transparent"),
          legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
          axis.text.y = element_text(face = "plain"),
          axis.title.x = element_text(face = "plain"),
          # axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
          plot.title = element_text(face = "plain", size = rel(1.15)),
          plot.subtitle = element_text(size = rel(0.95)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
        )
      p.TPR.FPR.bestcombo.arbitraryvsoptimal

      pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
      plotGG(plot = p.TPR.FPR.bestcombo.arbitraryvsoptimal, x = 0.05, y = 0.2, default.units = "cm",width = 8, height = 5.8)
      # plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

      #S2E the added TP by used the optimal parameter
      #compare the exomePeak2_default and exomePeak2_best using same FPR<=20% cutoff
      #obtain the corresponding cutoff and method combo for exomePeak2 best
      dt.method.best.AUC.nonNEB %>% distinct(Cell,Method,ComboName,IsBestCombo,Benchmark) %>% filter(Method=="exomePeak2") %>% distinct(Cell,Method,ComboName)
      dt.method.best.AUC.nonNEB %>% filter(Method=="exomePeak2" & IsBestCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
      #obtain the corresponding cutoff and method combo for exomePeak2 default
      dt.method.defaultcombo.AUC.nonNEB %>% dplyr::filter(Method=="exomePeak2" & IsDefaultCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293"))
      dt.selected.method.combo.cutoff.FPR20 <- rbind(dt.method.best.AUC.nonNEB %>% filter(Method=="exomePeak2" & IsBestCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>%
                                                       dplyr::select(Cell,Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Best"),
                                                     dt.method.defaultcombo.AUC.nonNEB %>% dplyr::filter(Method=="exomePeak2" & IsDefaultCombo==TRUE) %>% dplyr::filter(FPR==0.2 & str_detect(Benchmark,"HEK293")) %>%
                                                       dplyr::select(Cell,Method,ComboName,cutoff,FPR,TPR) %>% mutate(Group="Default")) %>% mutate(Combo = paste0(Method,"_",ComboName))
      #load exomePeak2 (Default and Best Combo) for HEK_SYSY and HEK_Abcam

      #Abcam
      dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 <- rbind(LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_Abcam/", dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"Abcam") & Group=="Best",Combo])) %>% mutate(Group="Best"),
                                                          LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_Abcam/", dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"Abcam") & Group=="Default",Combo])) %>% mutate(Group="Default")
      )
      #filter according to cutoff at FPR20
      dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 <- dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 %>% dplyr::filter( (Group=="Best" & score >= dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"Abcam") & Group=="Best",cutoff]) | (Group=="Default" & score >= dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"Abcam") & Group=="Default",cutoff]))
      #divide BestCombo m6A as Default- and Default+ (************************ strictly, use 0% overlap as cutoff for Default- peak definition **************************)
      dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap <- foreach(c = unique(dt.exomePeak2.BestDefault.m6A.Abcam_HEK293$Sample), .combine='rbind')%do%{
        bed.peak <- dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
        bed.default <- dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 %>% dplyr::filter(Sample==c & Group=="Default") %>% dplyr::select(seqnames,start,end,name,score,strand)
        dt.overlap <- bt.intersect(a=bed.peak, b=bed.default, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
        dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
        dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(Defaultoverlap=ifelse(sum.overlap/width>=0.0001,"Default+","Default-"))
        bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,Defaultoverlap), by="name") %>% mutate(Sample=c)
      }
      dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap[,.N,by=.(Sample,Defaultoverlap)] %>% dplyr::arrange(Sample,Defaultoverlap)
      #          Sample Defaultoverlap    N
      # 1: HEK_Abcam_mRNA1       Default+ 5178
      # 2: HEK_Abcam_mRNA1       Default- 5944
      # 3: HEK_Abcam_mRNA2       Default+ 5826
      # 4: HEK_Abcam_mRNA2       Default- 6657

      #check the Default- and Default+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
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
      dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Benchmarkoverlap <- foreach(c = unique(dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap$Sample), .combine='rbind')%do%{
        bed.peak <- dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
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
      #obtain the BestCombo only m6A+ gene
      dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.annotgene <- annot_peak(peak.bed=dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                                         strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                         annot_type = "gene")
      dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap <- dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap %>% left_join(x=.,y=dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
      dt.exomePeak2.BestCombo.specific.m6Agene.Abcam_HEK293 <-   dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
        group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestCombo.specific=all(Defaultoverlap=="Default-")) %>% as.data.table() %>%#if all peak of one gene is Default-, then this gene is BestCombo specific m6A+ gene
        dplyr::filter(BestCombo.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
      #obtain the BestCombo specific peaks overlapped with (high confident) benchmark sites
      dt.exomePeak2.BestCombo.specific.peak.benchmark.Abcam_HEK293 <-  dt.exomePeak2.BestCombo.specific.m6Agene.Abcam_HEK293 %>% filter(Defaultoverlap=="Default-" & BestCombo.specific==TRUE) %>%
        left_join(x=.,dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
        dplyr::filter(nBenchSite>0)
      #rank the BestCombo specific m6A+ genes with benchmark m6A sites support
      dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.Abcam_HEK293 <- dt.exomePeak2.BestCombo.specific.peak.benchmark.Abcam_HEK293 %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
        group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
        as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))

      #SYSY
      dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 <- rbind(LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_SYSY/", dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"SYSY") & Group=="Best",Combo])) %>% mutate(Group="Best"),
                                                         LoadPeakMethod(Method = "exomePeak2", Method.peak.dir = paste0("/data/m6A_calling_strategy/RIPPeakS_SYSY/", dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"SYSY") & Group=="Default",Combo])) %>% mutate(Group="Default")
      )
      #filter according to cutoff at FPR20
      dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 <- dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 %>% dplyr::filter( (Group=="Best" & score >= dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"SYSY") & Group=="Best",cutoff]) | (Group=="Default" & score >= dt.selected.method.combo.cutoff.FPR20[str_detect(Cell,"SYSY") & Group=="Default",cutoff]))
      #divide BestCombo m6A as Default- and Default+ (************************ strictly, use 0% overlap as cutoff for Default- peak definition **************************)
      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap <- foreach(c = unique(dt.exomePeak2.BestDefault.m6A.SYSY_HEK293$Sample), .combine='rbind')%do%{
        bed.peak <- dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
        bed.default <- dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 %>% dplyr::filter(Sample==c & Group=="Default") %>% dplyr::select(seqnames,start,end,name,score,strand)
        dt.overlap <- bt.intersect(a=bed.peak, b=bed.default, s=T, f=0.0001,F=0.0001, e=T, wao=T) %>% as.data.table() %>% dplyr::arrange(V4,V10)
        dt.overlap <- dt.overlap %>% mutate(width=V3-V2) %>% dplyr::select(name=V4,width,overlap=V13)
        dt.peak <- dt.overlap %>% group_by(name) %>% mutate(sum.overlap=sum(overlap)) %>% as.data.table() %>% mutate(Defaultoverlap=ifelse(sum.overlap/width>=0.0001,"Default+","Default-"))
        bed.peak %>% left_join(x=.,y=dt.peak %>% distinct(name,Defaultoverlap), by="name") %>% mutate(Sample=c)
      }
      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap[,.N,by=.(Sample,Defaultoverlap)] %>% dplyr::arrange(Sample,Defaultoverlap)
      #          Sample Defaultoverlap    N
      # 1: HEK_SYSY_mRNA1       Default+ 8166
      # 2: HEK_SYSY_mRNA1       Default- 4213
      # 3: HEK_SYSY_mRNA2       Default+ 7987
      # 4: HEK_SYSY_mRNA2       Default- 4071

      #check the Default- and Default+ peak overlapped with BenchmarkSites (Confident Benchmark Sites)
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
      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Benchmarkoverlap <- foreach(c = unique(dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap$Sample), .combine='rbind')%do%{
        bed.peak <- dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 %>% dplyr::filter(Sample==c & Group=="Best") %>% dplyr::select(seqnames,start,end,name,score,strand)
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
      #obtain the BestCombo only m6A+ gene
      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.annotgene <- annot_peak(peak.bed=dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap %>% dplyr::distinct(seqnames,start,end,name,score,strand),
                                                                        strand=T, fract = 0.5, gtf="/data/m6A_calling_strategy/RIPPeakS_resource/gencode.v44.annotation.gtf",
                                                                        annot_type = "gene")
      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap <- dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap %>% left_join(x=.,y=dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.annotgene %>% dplyr::select(name,OverlappedGenes),by="name")
      dt.exomePeak2.BestCombo.specific.m6Agene.SYSY_HEK293 <-   dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap %>% tidyr::separate_longer_delim(OverlappedGenes,delim = ",") %>%
        group_by(SampleName=Sample,OverlappedGenes) %>% mutate(BestCombo.specific=all(Defaultoverlap=="Default-")) %>% as.data.table() %>%#if all peak of one gene is Default-, then this gene is BestCombo specific m6A+ gene
        dplyr::filter(BestCombo.specific==TRUE) %>% dplyr::arrange(SampleName,OverlappedGenes)# keep all BestCombo specific m6A+ gene and  corresponding m6A peaks
      #obtain the BestCombo specific peaks overlapped with (high confident) benchmark sites
      dt.exomePeak2.BestCombo.specific.peak.benchmark.SYSY_HEK293 <-  dt.exomePeak2.BestCombo.specific.m6Agene.SYSY_HEK293 %>% filter(Defaultoverlap=="Default-" & BestCombo.specific==TRUE) %>%
        left_join(x=.,dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Benchmarkoverlap,by=c("name","Sample")) %>%
        dplyr::filter(nBenchSite>0)
      #rank the BestCombo specific m6A+ genes with benchmark m6A sites support
      dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.SYSY_HEK293 <- dt.exomePeak2.BestCombo.specific.peak.benchmark.SYSY_HEK293 %>% dplyr::distinct(Sample,OverlappedGenes,nBenchSite,nConfidentBenchSite) %>%
        group_by(OverlappedGenes) %>% mutate(nRep=n_distinct(Sample),avg.BenchSite=mean(nBenchSite), avg.ConfidentBenchSite=mean(nConfidentBenchSite)) %>%
        as.data.table() %>%  dplyr::arrange(desc(nRep),desc(avg.ConfidentBenchSite),desc(avg.BenchSite))

      #visualization of ratio of Default- ratio of exomePeak2 peaks in Abcam_HEK and SYSY_HEK
      dt.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB <- rbind(dt.exomePeak2.BestCombo.Abcam_HEK293.peaks.Defaultoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestComboOnly = n_distinct(name[Defaultoverlap=="Default-"])/n_distinct(name),
                                                                                                                                                                nPeak.BestComboOnly= n_distinct(name[Defaultoverlap=="Default-"]), nPeak.BestCombo=n_distinct(name)) %>%
                                                                        as.data.table() %>% distinct(Sample,Ratio.Peak.BestComboOnly,nPeak.BestComboOnly,nPeak.BestCombo),
                                                                      dt.exomePeak2.BestCombo.SYSY_HEK293.peaks.Defaultoverlap %>% group_by(Sample) %>% mutate(Ratio.Peak.BestComboOnly = n_distinct(name[Defaultoverlap=="Default-"])/n_distinct(name),
                                                                                                                                                               nPeak.BestComboOnly= n_distinct(name[Defaultoverlap=="Default-"]), nPeak.BestCombo=n_distinct(name)) %>%
                                                                        as.data.table() %>% distinct(Sample,Ratio.Peak.BestComboOnly,nPeak.BestComboOnly,nPeak.BestCombo)
      )
      dt.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB <- dt.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB %>% mutate(label1=paste0(nPeak.BestComboOnly,"/",nPeak.BestCombo), label2=paste0(round(Ratio.Peak.BestComboOnly*100,1),"%")) %>%
        mutate(Sample=factor(Sample, levels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"), labels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"))) %>%
        mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam",str_detect(Sample,"SYSY")~"SYSY"))
      p.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB <-
        ggplot(dt.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB, aes(y = Sample, x = Ratio.Peak.BestComboOnly, fill = Cell)) +
        geom_col(width = 0.5, show.legend = FALSE) +
        # percent label to the right of the bar
        geom_text(aes(label = label2, x =  Ratio.Peak.BestComboOnly + 0.02),
                  color = "black", size = 2, fontface = "plain", hjust = 0) +
        #N peak label to the left of the bar
        geom_text(aes(label = label1, x = Ratio.Peak.BestComboOnly/2),
                  color = "black", size = 2, fontface = "plain", hjust = 0) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1),expand=expansion(add=c(0,0.05))) +
        scale_fill_manual(values = c("#7b93c6","#acd9ee"),breaks = c("SYSY","Abcam")) +
        labs(x = "% of exomePeak2 best combo specific peaks", y = NULL, title = NULL, subtitle = NULL,caption = NULL) +
        theme_minimal(base_size = 6,base_family = "Helvetica") +
        theme(
          axis.text.y = element_text(face = "plain"),
          axis.title.x = element_text(face = "plain"),
          plot.title = element_text(face = "plain", size = rel(1.15)),
          plot.subtitle = element_text(size = rel(0.95)),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()
        )
      p.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB
      #visualization of count of BestCombo only m6A+ gene with benchmark m6As and highly confident m6As
      dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB <- rbind(dt.exomePeak2.BestCombo.specific.peak.benchmark.Abcam_HEK293 %>% group_by(Sample) %>%
                                                                                      mutate(No.BestComboOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                             No.BestComboOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                                      distinct(Sample,No.BestComboOnly.m6AGene.BenchmarkSite,No.BestComboOnly.m6AGene.ConfidentBenchmarkSite),
                                                                                    dt.exomePeak2.BestCombo.specific.peak.benchmark.SYSY_HEK293 %>% group_by(Sample) %>%
                                                                                      mutate(No.BestComboOnly.m6AGene.BenchmarkSite=n_distinct(OverlappedGenes[nBenchSite>0]),
                                                                                             No.BestComboOnly.m6AGene.ConfidentBenchmarkSite=n_distinct(OverlappedGenes[nConfidentBenchSite>0])) %>% as.data.table() %>%
                                                                                      distinct(Sample,No.BestComboOnly.m6AGene.BenchmarkSite,No.BestComboOnly.m6AGene.ConfidentBenchmarkSite)
      )
      dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB <- dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB %>%
        pivot_longer(cols = c("No.BestComboOnly.m6AGene.BenchmarkSite", "No.BestComboOnly.m6AGene.ConfidentBenchmarkSite"),values_to = "No.BestComboOnly.m6AGene", names_to = "GeneGroup",names_prefix = "No.BestComboOnly.m6AGene.") %>%
        as.data.table() %>% mutate(GeneGroup=factor(GeneGroup, levels=c("BenchmarkSite","ConfidentBenchmarkSite"), labels=c("Benchmark m6As","Confident benchmark m6As"))) %>%
        mutate(Sample=factor(Sample, levels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"), labels=c("HEK_Abcam_mRNA1","HEK_Abcam_mRNA2","HEK_SYSY_mRNA1","HEK_SYSY_mRNA2"))) %>%
        mutate(Cell=case_when(str_detect(Sample,"Abcam")~"Abcam",str_detect(Sample,"SYSY")~"SYSY"))
      dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB <- dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB %>% mutate(label.count=paste0(paste0("No.Gene=",No.BestComboOnly.m6AGene)))
      p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB  <- ggplot(data=dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB, aes(x=Sample,y=No.BestComboOnly.m6AGene, fill=GeneGroup))+
        geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
        # percent label to the left  bar
        geom_text(data = dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB %>% dplyr::filter(GeneGroup=="Benchmark m6As"),
                  aes(label = label.count, y =  No.BestComboOnly.m6AGene/2+0.02, x=Sample),
                  color = "black", size = 2, fontface = "plain",angle=90,nudge_x = -0.15) +
        # percent label to the right  bar
        geom_text(data = dt.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB %>% dplyr::filter(GeneGroup=="Confident benchmark m6As"),
                  aes(label = label.count, y =  No.BestComboOnly.m6AGene/2+0.02, x=Sample),
                  color = "black", size = 2, fontface = "plain",angle=90,nudge_x = 0.15) +
        scale_fill_manual(values=c("#cac0e3","#6766b1"),breaks=c("Benchmark m6As","Confident benchmark m6As"))+
        labs(x=NULL,y=str_wrap("N of genes harbor only exomePeak2 best combo specific peaks",width=40))+
        guides(fill=guide_legend(title = "m6A peaks contain"))+
        scale_y_continuous(expand = expansion(add=c(0,0.05)))+
        theme_minimal(base_size = 6,base_family = "Helvetica") +
        theme(
          legend.position = "inside",legend.position.inside = c(0.8,0.9),legend.key.height = unit(8,"pt"),legend.key.width = unit(8,"pt"),
          axis.text.y = element_text(face = "plain"),
          axis.title.x = element_text(face = "plain"),
          axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
          plot.title = element_text(face = "plain", size = rel(1.15)),
          plot.subtitle = element_text(size = rel(0.95)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
        )
      p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB

      #visualization of IGV plot of exomePeak2 BestCombo specific m6A+ gene in Abcam_HEK293 and SYSY_HEK293
      dt.gene.hg38 <- genomation::gffToGRanges("~/genome_db/gencode.v44.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)
      # dt.gene.mm39 <- genomation::gffToGRanges("~/genome_db/gencode.vM33.annotation.gtf") %>% as.data.table() %>% filter(type=="gene") %>% dplyr::select(seqnames,start,end,gene_id,score,strand,gene_name)

      #Abcam_HEK
      selected.exomePeak2.BestComboOnly.m6Agene.Abcam.HEK <-  dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.Abcam_HEK293 %>% dplyr::filter(nRep==2) %>%
        distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=5)
      selected.exomePeak2.BestComboOnly.m6Agene.Abcam.HEK
      dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==unique(selected.exomePeak2.BestComboOnly.m6Agene.Abcam.HEK$OverlappedGenes)[1])
      StudySamples <- c("HEK_Abcam_mRNA1", "HEK_Abcam_mRNA2")
      ## Create a page (7.5*7.5cm)
      window.size=20#100kb
      pseudocount=1
      track.height=0.5
      pdf(paste0("Fig2S_Abcam_HEK293_",dt.selected.gene$gene_name,"_exomePeak2_bestcombo_specific_peak.pdf"))
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

      # #add exomePeak2 default combo peak in this region
      for(s in 1:length(StudySamples)){
        dt.peak.region <- dt.exomePeak2.BestDefault.m6A.Abcam_HEK293 %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
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
      legendPlot.exomePeak2.peak <- plotLegend(legend = c("BestCombo","DefaultCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 2.5, height = 1,
                                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                               title = expression("exomePeak2 peak"),fontface="plain")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

      dev.off()

      #SYSY_HEK
      selected.exomePeak2.BestComboOnly.m6Agene.SYSY.HEK <-  dt.exomePeak2.BestCombo.specific.m6Agene.benchmark.SYSY_HEK293 %>% dplyr::filter(nRep==2) %>%
        distinct(OverlappedGenes,nRep,avg.BenchSite, avg.ConfidentBenchSite) %>% filter(!str_detect(OverlappedGenes,"ENSG")) %>% dplyr::arrange(desc(avg.ConfidentBenchSite)) %>%  slice_head(n=5)
      selected.exomePeak2.BestComboOnly.m6Agene.SYSY.HEK
      dt.selected.gene <- dt.gene.hg38 %>% dplyr::filter(gene_name==unique(selected.exomePeak2.BestComboOnly.m6Agene.SYSY.HEK$OverlappedGenes)[1])
      StudySamples <- c("HEK_SYSY_mRNA1", "HEK_SYSY_mRNA2")
      ## Create a page (7.5*7.5cm)
      window.size=25#100kb
      pseudocount=1
      track.height=0.5
      pdf(paste0("Fig2S_SYSY_HEK293_",dt.selected.gene$gene_name,"_exomePeak2_bestcombo_specific_peak.pdf"))
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

      # #add exomePeak2 default combo peak in this region
      for(s in 1:length(StudySamples)){
        dt.peak.region <- dt.exomePeak2.BestDefault.m6A.SYSY_HEK293 %>% dplyr::filter(Sample==StudySamples[s] & seqnames==region$chrom & (start >=region$chromstart &  end <= region$chromend))
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
      legendPlot.exomePeak2.peak <- plotLegend(legend = c("BestCombo","DefaultCombo"), fill = c("#e9abac","#7bb0d5"),border = F,x =0.5 + 2.5, y = 1.5+(2*length(StudySamples)+3)*track.height, width = 2.5, height = 1,
                                               just = c("left", "top"),default.units = "cm",fontsize = 6,fontface="plain",orientation = "h",draw = T,
                                               title = expression("exomePeak2 peak"),fontface="plain")
      #add rect for whole panel
      plotRect(x=0.5-0.45,y=0.01,width=7.5, height = 1.5+(2*length(StudySamples)+4)*track.height+0.5,linecolor = "black", lwd = 0.5,lty = 2,just = c("left","top"),default.units = "cm")

      dev.off()

      #save intermediate variables for re-use
      save.image("sm6APeak_FigureS2.intermediate.results.RDS")


      ####################### combine all figure S2 together
      figS2.list <- list(p.GLORI.AUC.Method.DefaultCombo.nonNEB=p.GLORI.AUC.Method.DefaultCombo.nonNEB,
                         p.allcombo.GLORI.AUC.nonNEB=p.allcombo.GLORI.AUC.nonNEB,
                         p.GLORI.AUC.Method.best.nonNEB=p.GLORI.AUC.Method.best.nonNEB,
                         p.TPR.FPR.bestcombo.arbitraryvsoptimal=p.TPR.FPR.bestcombo.arbitraryvsoptimal,
                         p.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB=p.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB,
                         p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB=p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB)
      saveRDS(figS2.list, file="FigureS2.plot.list.RDS")

      pdf("FigureS2.Parameter_refinement_could_improve_performance_of_m6A_calling_regardless_of_m6A_antibody.pdf",width = 8.2677, height = 11.693)
      pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
      plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.GLORI.AUC.Method.DefaultCombo.nonNEB, x = 0.05, y=0.2, default.units = "cm",width = 9.8, height = 5.8)
      plotText(label = "B", fontsize = 8, fontface = "bold",x = 10, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.allcombo.GLORI.AUC.nonNEB, x = 10, y=0.2, default.units = "cm",width = 8.0, height = 5.8)
      plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.GLORI.AUC.Method.best.nonNEB, x = 0.05, y=6.2, default.units = "cm",width = 9.8, height = 5.8)
      plotText(label = "D", fontsize = 8, fontface = "bold",x = 10, y = 6.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.TPR.FPR.bestcombo.arbitraryvsoptimal, x = 10, y=6.2, default.units = "cm",width = 8.0, height = 5.8)
      plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.exomePeak2.BestCombo.specific.ratio.m6A.peak.nonNEB, x = 0.05, y=12.2, default.units = "cm",width = 7.5, height = 5.8)
      plotText(label = "G", fontsize = 8, fontface = "bold",x = 8.5, y = 12.2,just = c("top","left"), default.units = "cm",draw=T)
      plotRect(x = 8.5, y = 12.4,just = c("top","left"),width = 8.5,height = 5.6,default.units = "cm",draw=T)
      plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05, y = 18.2,just = c("top","left"), default.units = "cm",draw=T)
      plotGG(plot = p.count.exomePeak2.BestCombo.specific.m6Agene.with.benchmark.nonNEB, x = 0.05, y=18.2, default.units = "cm",width = 7.5, height = 5.8)
      plotText(label = "H", fontsize = 8, fontface = "bold",x = 8.5, y = 18.2,just = c("top","left"), default.units = "cm",draw=T)
      plotRect(x = 8.5, y = 18.4,just = c("top","left"),width = 8.5,height = 5.6,default.units = "cm",draw=T)

      dev.off()
