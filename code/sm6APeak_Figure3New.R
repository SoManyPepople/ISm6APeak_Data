
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

################# FigSNew_Benchmark_and_ISm6APeak_training_in_independent_Millipore_mESC_dataset_and_NEB_HeLa_datatset##################
#load all AUC of 285 method combo using original score of Millipore_mESC

#load all AUC of 285x32 method comboxoption  of Miilipore_mESC
#load all combo option AUC
dt.pruned.FPR.TPR.option.AUC.Millipore.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.Millipore.mESC.RDS")
dt.pruned.FPR.TPR.option.AUC.Millipore.mESC %>% distinct(Method,ComboName,optionID,AUC.scaled)
#load all combo AUC
dt.pruned.FPR.TPR.combo.AUC.Millipore.mESC <- readRDS("dt.pruned.FPR.TPR.AUC.Millipore.mESC.all.combo.RDS")
dt.BestCombo.Millipore <- dt.pruned.FPR.TPR.combo.AUC.Millipore.mESC %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.DefaultCombo.AUC.Millipore <- dt.pruned.FPR.TPR.combo.AUC.Millipore.mESC %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName"))
                                 
dt.BestCombo.AUC.Millipore <- dt.pruned.FPR.TPR.combo.AUC.Millipore.mESC %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.Millipore %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) 

dt.AUC.Millipore.all.option <- dt.pruned.FPR.TPR.option.AUC.Millipore.mESC %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) 

dt.BestComboOption.Millipore <- dt.AUC.Millipore.all.option %>%  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20))

dt.AUC.Millipore.all.option <- dt.AUC.Millipore.all.option %>% left_join(x=.,y=dt.BestComboOption.Millipore %>% mutate(IsBestComboOption=TRUE) %>% dplyr::select(Method,ComboName,optionID,IsBestComboOption),by=c("Method","ComboName","optionID")) %>%
  mutate(ComboGroup=case_when(IsBestComboOption==TRUE ~ "BestComboOption", .default = "OthersComboOption")) %>% dplyr::select(-IsBestComboOption)
dt.AUC.Millipore.all.option <- rbind(dt.AUC.Millipore.all.option %>% mutate(Group=ComboGroup) %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                               dt.DefaultCombo.AUC.Millipore %>% mutate(optionID=NA, Group="DefaultCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                               dt.BestCombo.AUC.Millipore %>% mutate(optionID=NA, Group="BestCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group))
dt.AUC.Millipore.all.option <- dt.AUC.Millipore.all.option %>% mutate(Method=factor(Method,levels=unique(dt.BestComboOption.Millipore$Method)),Group=factor(Group,levels=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo")))
p.AUC.Millipore.all.option <- ggplot(dt.AUC.Millipore.all.option, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.AUC.Millipore.all.option %>% dplyr::filter(Group =="OthersComboOption"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=0.1,alpha=0.15, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.AUC.Millipore.all.option %>% dplyr::filter(Group =="DefaultCombo"), aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=1.5)+
  geom_jitter(data = dt.AUC.Millipore.all.option %>% dplyr::filter(Group =="BestCombo") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=1.5)+
  geom_point(data = dt.AUC.Millipore.all.option %>% dplyr::filter(Group =="BestComboOption") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#eb4601",size=3)+
  scale_y_continuous(n.breaks = 6)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  # facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled)", title = paste0("9120 ComboOption combination"),subtitle = "Millipore_mESC") +
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9),hjust = 0.5),
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
p_col_legend <- ggplot(data=dt.AUC.Millipore.all.option,aes(x = Method, y = AUC.scaled)) +
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
    plot.subtitle = element_text(size = rel(0.95),hjust = 0.5),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
# extract legend grobs
leg_col <- cowplot::get_legend(p_col_legend)
# assemble: main plot + legends at arbitrary positions
p.AUC.Millipore.all.option  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.AUC.Millipore.all.option , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1.15, y = 0.8, width = 0.22, height = 0.25, hjust = 1.2, vjust = 0.5)

p.AUC.Millipore.all.option



#load all AUC of 285 method combo using original score of NEB_HeLa

#load all AUC of 285x32 method comboxoption  of NEB_HeLa
#load all combo option AUC
dt.pruned.FPR.TPR.option.AUC.NEB.HeLa <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HeLa.RDS")
dt.pruned.FPR.TPR.option.AUC.NEB.HeLa %>% distinct(Method,ComboName,optionID,AUC.scaled)
#load all combo AUC
dt.pruned.FPR.TPR.combo.AUC.NEB.HeLa <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HeLa.all.combo.RDS")
dt.BestCombo.NEBHeLa <- dt.pruned.FPR.TPR.combo.AUC.NEB.HeLa %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.DefaultCombo.AUC.NEBHeLa <- dt.pruned.FPR.TPR.combo.AUC.NEB.HeLa %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName"))

dt.BestCombo.AUC.NEBHeLa <- dt.pruned.FPR.TPR.combo.AUC.NEB.HeLa %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.NEBHeLa %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) 

dt.AUC.NEBHeLa.all.option <- dt.pruned.FPR.TPR.option.AUC.NEB.HeLa %>% filter(str_detect(Benchmark,"GLORI")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) 

dt.BestComboOption.NEBHeLa <- dt.AUC.NEBHeLa.all.option %>%  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20))

dt.AUC.NEBHeLa.all.option <- dt.AUC.NEBHeLa.all.option %>% left_join(x=.,y=dt.BestComboOption.NEBHeLa %>% mutate(IsBestComboOption=TRUE) %>% dplyr::select(Method,ComboName,optionID,IsBestComboOption),by=c("Method","ComboName","optionID")) %>%
  mutate(ComboGroup=case_when(IsBestComboOption==TRUE ~ "BestComboOption", .default = "OthersComboOption")) %>% dplyr::select(-IsBestComboOption)
dt.AUC.NEBHeLa.all.option <- rbind(dt.AUC.NEBHeLa.all.option %>% mutate(Group=ComboGroup) %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                                     dt.DefaultCombo.AUC.NEBHeLa %>% mutate(optionID=NA, Group="DefaultCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                                     dt.BestCombo.AUC.NEBHeLa %>% mutate(optionID=NA, Group="BestCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group))
dt.AUC.NEBHeLa.all.option <- dt.AUC.NEBHeLa.all.option %>% mutate(Method=factor(Method,levels=unique(dt.BestComboOption.NEBHeLa$Method)),Group=factor(Group,levels=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo")))
p.AUC.NEBHeLa.all.option <- ggplot(dt.AUC.NEBHeLa.all.option, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="OthersComboOption"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=0.1,alpha=0.15, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="DefaultCombo"), aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=1.5)+
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="BestCombo") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=1.5)+
  geom_point(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="BestComboOption") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#eb4601",size=3)+
  scale_y_continuous(n.breaks = 5)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  # facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled)", title = paste0("9120 ComboOption combination"),subtitle = "NEB_HeLa") +
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9),hjust = 0.5),
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
p_col_legend <- ggplot(data=dt.AUC.NEBHeLa.all.option,aes(x = Method, y = AUC.scaled)) +
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
p.AUC.NEBHeLa.all.option  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.AUC.NEBHeLa.all.option , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1.15, y = 0.8, width = 0.22, height = 0.25, hjust = 1.2, vjust = 0.5)

p.AUC.NEBHeLa.all.option


#AUC display supeiror ISm6APeak training using NEB_HeLa 
dt.FPR.TPR.FourMode.ISm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
dt.FPR.TPR.ISm6APeak.NEBHeLa <- dt.FPR.TPR.FourMode.ISm6APeak %>% dplyr::filter(SelectedMode=="NEB_HeLa" & SelectedSample=="NEB_HeLa")
dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa <- readRDS("dt.pruned.FPR.TPR.AUC.NEB.HeLa.RDS")
dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa <- dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=4) %>%
  inner_join(x=.,y=dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa)
dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20)
#combine ISm6APeak with those top4 best combo
dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa <- rbind(dt.FPR.TPR.ISm6APeak.NEBHeLa  %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20) %>% mutate(Method="ISm6APeak"),
                                                    dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa %>%  dplyr::select(FPR,TPR,AUC.scaled,TPR20,Method) )

dt.TPR.FPR.ISm6APeak.Method.BestOption.label.NEBHeLa <- dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa %>% distinct(Method,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% 
  # group_by(Cell) %>%
  mutate(ypos=tail(seq(0.7,0,by= -0.6/5),5)*max(dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa$TPR)) %>% mutate(xpos=0.6) %>% as.data.table()
p.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa  <-  ggplot(data=dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa %>% mutate(Cell="NEB_HeLa"),aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of ISm6APeak or method's BestCombo",y=str_wrap("Surrogate TPR of ISm6APeak or method's BestCombo",width = 40))+
  geom_text(data=dt.TPR.FPR.ISm6APeak.Method.BestOption.label.NEBHeLa, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#eb4601",c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0")),breaks = c("ISm6APeak",c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools")),guide="none")+
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
    plot.title = element_text(face = "plain", size = rel(1.15),hjust = 0.5),
    plot.subtitle = element_text(size = rel(0.95)),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
p.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa

#AUC display ISm6APeak of NEB model and NEB_HeLa model in NEB_HeLa
dt.FPR.TPR.FourMode.ISm6APeak <- readRDS("dt.pruned.TPR.FPR.fivesamples.using.four.different.M6APeakS.mode.RDS")
dt.FPP.TPR.ISm6APeak.NEBvsOriginal <- dt.FPR.TPR.FourMode.ISm6APeak %>% dplyr::filter(SelectedSample == "NEB_HeLa") %>% filter((SelectedMode==SelectedSample) | SelectedMode=="NEB") %>%
  mutate(ModelType=case_when(SelectedMode==SelectedSample ~ "Original", .default = SelectedMode)) %>% mutate(Cell="NEB_HeLa")
dt.TPR.FPR.ISm6APeak.label.NEBvsOriginal <- dt.FPP.TPR.ISm6APeak.NEBvsOriginal %>% distinct(Cell,ModelType,AUC.scaled,TPR20) %>%  mutate(label=paste0(ModelType,"Model,AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>% group_by(Cell) %>%
  mutate(ypos=head(seq(0.7,0,by= -0.6/5),2)*max(dt.FPP.TPR.ISm6APeak.NEBvsOriginal$TPR)) %>% mutate(xpos=0.5) %>% as.data.table()
p.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa  <-  ggplot(data=dt.FPP.TPR.ISm6APeak.NEBvsOriginal ,aes(x=FPR,y=TPR,color=ModelType))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of ISm6APeak from different anti-m6A model",y=str_wrap("Surrogate TPR of ISm6APeak from different anti-m6A model",width = 40))+
  geom_text(data=dt.TPR.FPR.ISm6APeak.label.NEBvsOriginal, aes(x=xpos,y=ypos,color=ModelType,label=label),size=1.8,fontface="bold")+
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
p.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa


####  AUC of NEB_HeLa based on alternative eTAM_HeLa ####################

dt.pruned.alternative.FPR.TPR.option.AUC.NEB.HeLa <- readRDS("dt.pruned.alternative.FPR.TPR.AUC.NEB.HeLa.RDS")
dt.pruned.alternative.FPR.TPR.option.AUC.NEB.HeLa %>% distinct(Method,ComboName,optionID,AUC.scaled)
#load all combo AUC
dt.pruned.alternative.FPR.TPR.combo.AUC.NEB.HeLa <- readRDS("dt.pruned.alternative.FPR.TPR.AUC.NEB.HeLa.all.combo.RDS")
dt.BestCombo.NEBHeLa <- dt.pruned.alternative.FPR.TPR.combo.AUC.NEB.HeLa %>% filter(str_detect(Benchmark,"eTAM")) %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table()
dt.method.combo <- readRDS("dt.method.combo.filtered.RDS")
dt.DefaultCombo.AUC.NEBHeLa <- dt.pruned.alternative.FPR.TPR.combo.AUC.NEB.HeLa %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.method.combo[IsDefaultCombo==TRUE,] %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName"))

dt.BestCombo.AUC.NEBHeLa <- dt.pruned.alternative.FPR.TPR.combo.AUC.NEB.HeLa %>% distinct(Method,ComboName,AUC.scaled,TPR20) %>% inner_join(x=.,y=dt.BestCombo.NEBHeLa %>% dplyr::distinct(Method,ComboName),by=c("Method","ComboName")) 

dt.AUC.NEBHeLa.all.option <- dt.pruned.alternative.FPR.TPR.option.AUC.NEB.HeLa %>% filter(str_detect(Benchmark,"eTAM")) %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) 

dt.BestComboOption.NEBHeLa <- dt.AUC.NEBHeLa.all.option %>%  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20))

dt.AUC.NEBHeLa.all.option <- dt.AUC.NEBHeLa.all.option %>% left_join(x=.,y=dt.BestComboOption.NEBHeLa %>% mutate(IsBestComboOption=TRUE) %>% dplyr::select(Method,ComboName,optionID,IsBestComboOption),by=c("Method","ComboName","optionID")) %>%
  mutate(ComboGroup=case_when(IsBestComboOption==TRUE ~ "BestComboOption", .default = "OthersComboOption")) %>% dplyr::select(-IsBestComboOption)
dt.AUC.NEBHeLa.all.option <- rbind(dt.AUC.NEBHeLa.all.option %>% mutate(Group=ComboGroup) %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                                   dt.DefaultCombo.AUC.NEBHeLa %>% mutate(optionID=NA, Group="DefaultCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group),
                                   dt.BestCombo.AUC.NEBHeLa %>% mutate(optionID=NA, Group="BestCombo") %>% dplyr::select(Method,ComboName,optionID,AUC.scaled,TPR20,Group))
dt.AUC.NEBHeLa.all.option <- dt.AUC.NEBHeLa.all.option %>% mutate(Method=factor(Method,levels=unique(dt.BestComboOption.NEBHeLa$Method)),Group=factor(Group,levels=c("BestComboOption","OthersComboOption","BestCombo","DefaultCombo")))
p.alternative.AUC.NEBHeLa.all.option <- ggplot(dt.AUC.NEBHeLa.all.option, aes(x = Method, y = AUC.scaled)) +
  # boxplots without showing outliers (we plot points separately)
  geom_boxplot(width = 0.55,
               outlier.shape = NA,
               fill = "transparent",
               colour = "grey5",
               size = 0.4) +
  # individual points: shape by ComboGroup, fill by ComboGroup, black border
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="OthersComboOption"),aes(x = Method, y = AUC.scaled),shape=22, color="grey70", size=0.1,alpha=0.15, width = 0.18, height = 0,  stroke = 0.6) +
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="DefaultCombo"), aes(x = Method, y = AUC.scaled), shape=18, color= "#7bb0d5",width = 0.18, size=1.5)+
  geom_jitter(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="BestCombo") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#e9abac",width = 0.18, size=1.5)+
  geom_point(data = dt.AUC.NEBHeLa.all.option %>% dplyr::filter(Group =="BestComboOption") ,aes(x = Method, y = AUC.scaled), shape=18, color= "#eb4601",size=3)+
  scale_y_continuous(n.breaks = 5)+
  # show mean as a black diamond
  # stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6, colour = "black", stroke = 0.8) +
  # scales for shapes and fills
  # scale_shape_manual(values = c(18,21,22),breaks=c("Best","Default","Others"), name = "Combo group") +
  # # scale_size_manual(values = c(2,2,0.5), breaks = c("Best","Default","Others"), name = "Combo group") +
  # scale_alpha_manual(values = c(1,1,0.5), breaks = c("Best","Default","Others"), name = "Combo group")+
  # scale_color_manual(values = c("#e9abac","#7bb0d5","grey88"), breaks = c("Best","Default","Others"), name = "Combo group") +
  # facet_wrap(~Cell)+
  # labels
  labs(x = NULL, y = "AUC (scaled) based on eTAM_HeLa m6As", title = paste0("9120 ComboOption combination"),subtitle = "NEB_HeLa") +
  # Nature-like theme tweaks
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    text = element_text(family = "Helvetica", colour = "black"),
    plot.title = element_text(face = "plain", size = rel(1.05)),
    plot.subtitle = element_text(size = rel(0.9),hjust = 0.5),
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
p_col_legend <- ggplot(data=dt.AUC.NEBHeLa.all.option,aes(x = Method, y = AUC.scaled)) +
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
p.alternative.AUC.NEBHeLa.all.option  <- cowplot::ggdraw() +
  cowplot::draw_plot(p.alternative.AUC.NEBHeLa.all.option , 0, 0, 1, 1) +
  # place color legend (right)
  cowplot::draw_grob(leg_col, x = 1.15, y = 0.8, width = 0.22, height = 0.25, hjust = 1.2, vjust = 0.5)

p.alternative.AUC.NEBHeLa.all.option

#AUC display supeiror ISm6APeak training using NEB_HeLa 
#prune the FPR/TPR and calculate AUC for alternative FPR/FPR in NEB_HeLa using NEB_HeLa model and NEB model based on eTAM_HeLa m6As
#prune TPR/FPR
dt.pruned.alternative.TPR.FPR.M6APeakS.NEBHeLa <- foreach(MODE=c("NEB","SYSY","Abcam","Millipore","NEB_HeLa")[c(1,5)],.combine='rbind')%do%{
  foreach(SAMPLE=c("NEB","SYSY","Abcam","Millipore","NEB_mESC","NEB_HeLa")[6],.combine='rbind')%do%{
    if(!file.exists(paste0("dt.",SAMPLE,".TPR.FPR.",MODE,".M6APeakS.RDS"))){
      message(paste0("Error! Can't find TPR FPR result for ", SAMPLE, " using ", MODE))
    }else{
      dt.FPR.TPR <- readRDS(paste0("dt.",SAMPLE,".TPR.FPR.",MODE,".M6APeakS.using.eTAM.RDS"))
      dt.avg.FPR.TPR <- dt.FPR.TPR %>% dplyr::select(cutoff=MergedFPR,Sample,FPR,TPR) %>%
        group_by(cutoff) %>% mutate(FPR.avg=mean(FPR),TPR.avg=mean(TPR)) %>%
        as.data.table() %>% distinct(cutoff,FPR.avg,TPR.avg) %>% dplyr::rename(FPR=FPR.avg, TPR=TPR.avg) %>%
        dplyr::arrange(cutoff)
      attr(dt.avg.FPR.TPR,"groups") <- NULL
      #prune TPR and FPR and calculate AUC
      message(paste0("TPR FPR prune for ", SAMPLE, " using ",MODE, " M6APeakS mode"))
      dt.pruned.FPR.TPR <-  TPR.FPR.pruned(dt.TPR.FPR.fc = dt.avg.FPR.TPR[,.(TPR,FPR,cutoff)]) %>%
        mutate(SelectedMode=MODE,SelectedSample=SAMPLE) %>% distinct()
      dt.pruned.FPR.TPR.AUC <- dt.pruned.FPR.TPR %>% mutate(AUC=TargetAUC(FPR = dt.pruned.FPR.TPR$FPR, TPR=dt.pruned.FPR.TPR$TPR,TargetRegion = 0.5)) %>%
        mutate(AUC=round(AUC,2),AUC.scaled=round(AUC/0.5,2), maxTPRmFPR=max(TPR-FPR), TPR20=max(TPR[FPR<=0.2]))
      dt.pruned.FPR.TPR.AUC
    }
  }
}
dt.pruned.alternative.TPR.FPR.M6APeakS.NEBHeLa %>% distinct(SelectedMode,SelectedSample,AUC.scaled,maxTPRmFPR,TPR20) %>% dplyr::arrange(SelectedSample,desc(AUC.scaled))
# SelectedMode SelectedSample AUC.scaled maxTPRmFPR TPR20
# 1:     NEB_HeLa       NEB_HeLa       0.36       0.27  0.34
# 2:          NEB       NEB_HeLa       0.34       0.26  0.00

dt.alternative.FPR.TPR.ISm6APeak.NEBHeLa <- dt.pruned.alternative.TPR.FPR.M6APeakS.NEBHeLa %>% dplyr::filter(SelectedMode=="NEB_HeLa")

dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa <- readRDS("dt.pruned.alternative.FPR.TPR.AUC.NEB.HeLa.RDS")
dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa <- dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa %>% distinct(Method,ComboName,optionID,AUC.scaled,TPR20) %>%
  group_by(Method) %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=1) %>% as.data.table() %>% dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% slice_head(n=4) %>%
  inner_join(x=.,y=dt.pruned.TPR.FPR.Method.AllOption.NEBHeLa)
dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa %>% distinct(Method,ComboName,optionID,Benchmark,AUC.scaled,TPR20)
#combine ISm6APeak with those top4 best combo
dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa <- rbind(dt.alternative.FPR.TPR.ISm6APeak.NEBHeLa  %>% dplyr::select(FPR,TPR,AUC.scaled,TPR20) %>% mutate(Method="ISm6APeak"),
                                                        dt.pruned.TPR.FPR.Method.BestOption.NEBHeLa %>%  dplyr::select(FPR,TPR,AUC.scaled,TPR20,Method) )

dt.TPR.FPR.ISm6APeak.Method.BestOption.label.NEBHeLa <- dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa %>% distinct(Method,AUC.scaled,TPR20) %>%  mutate(label=paste0(Method,",AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(desc(AUC.scaled),desc(TPR20)) %>% 
  # group_by(Cell) %>%
  mutate(ypos=tail(seq(0.7,0,by= -0.6/5),5)*max(dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa$TPR)) %>% mutate(xpos=0.6) %>% as.data.table()
p.alternative.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa  <-  ggplot(data=dt.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa %>% mutate(Cell="NEB_HeLa") ,aes(x=FPR,y=TPR,color=Method))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR based on eTAM_HeLa m6As",y=str_wrap("Surrogate TPR of based on eTAM_HeLa m6As",width = 40))+
  geom_text(data=dt.TPR.FPR.ISm6APeak.Method.BestOption.label.NEBHeLa, aes(x=xpos,y=ypos,color=Method,label=label),size=1.8,fontface="bold")+
  scale_color_manual(values = c("#eb4601",c("#e9abac","#7bb0d5","#cfe4b6","#cab0d2","#f6c780","#84cdc0")),breaks = c("ISm6APeak",c("MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools")),guide="none")+
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
p.alternative.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa

#AUC display ISm6APeak of NEB model and NEB_HeLa model in NEB_HeLa based on alternative eTAM_HeLa

dt.FPP.TPR.ISm6APeak.NEBvsOriginal <- dt.pruned.alternative.TPR.FPR.M6APeakS.NEBHeLa %>% dplyr::filter(SelectedSample == "NEB_HeLa") %>% filter((SelectedMode==SelectedSample) | SelectedMode=="NEB") %>%
  mutate(ModelType=case_when(SelectedMode==SelectedSample ~ "Original", .default = SelectedMode)) %>% mutate(Cell="NEB_HeLa")
dt.TPR.FPR.ISm6APeak.label.NEBvsOriginal <- dt.FPP.TPR.ISm6APeak.NEBvsOriginal %>% distinct(Cell,ModelType,AUC.scaled,TPR20) %>%  mutate(label=paste0(ModelType,"Model,AUC=",AUC.scaled, ",TPR20%=",TPR20*100,"%")) %>%
  dplyr::arrange(Cell,desc(AUC.scaled),desc(TPR20)) %>% group_by(Cell) %>%
  mutate(ypos=head(seq(0.7,0,by= -0.6/5),2)*max(dt.FPP.TPR.ISm6APeak.NEBvsOriginal$TPR)) %>% mutate(xpos=0.5) %>% as.data.table()
p.alternative.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa  <-  ggplot(data=dt.FPP.TPR.ISm6APeak.NEBvsOriginal ,aes(x=FPR,y=TPR,color=ModelType))+
  geom_step(linewidth=0.5)+
  labs(x="Surrogate FPR of ISm6APeak based on eTAM_HeLa m6As",y=str_wrap("Surrogate TPR of ISm6APeak based on eTAM_HeLa m6As",width = 40))+
  geom_text(data=dt.TPR.FPR.ISm6APeak.label.NEBvsOriginal, aes(x=xpos,y=ypos,color=ModelType,label=label),size=1.8,fontface="bold")+
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
p.alternative.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa



save.image("sm6APeak_FigureS_HeLa.intermediate.results.RDS")


############### combine plots together #############
pdf("NewFigureS3.ISm6APeak_traning_from_independent_datasets_display_robustness.pdf",width = 8.2677, height = 11.693)

pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
#row1
# ggsave(plot=p.HeLa.shared.benchmark.m6A.sites, dpi = 300,units = "cm",device = cairo_pdf,
#        width=5.8,height = 5.8,filename = "FigureS1_p.HeLa.shared.benchmark.m6A.sites.pdf")
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotRect(x=0.05,y=0.2,width = 5.8,height = 5.8,just = c("top","left"), default.units = "cm",draw=T)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 6, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.AUC.Millipore.all.option, x = 6, y=0.2, default.units = "cm",width = 5.9, height = 5.8)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 12, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.AUC.NEBHeLa.all.option, x = 12, y=0.2, default.units = "cm",width = 5.9, height = 5.8)

#row2
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 6,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa, x = 0.05, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 6, y = 6,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa, x = 6, y=6.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 12, y = 6,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.alternative.AUC.NEBHeLa.all.option, x = 12, y=6.2, default.units = "cm",width = 5.9, height = 5.8)
#row3
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05, y = 12,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.alternative.TPR.FPR.ISm6APeak.Method.BestOption.NEBHeLa, x = 0.05, y=12.2, default.units = "cm",width = 5.8, height = 5.8)
plotText(label = "H", fontsize = 8, fontface = "bold",x = 6, y = 12,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.alternative.TPR.FPR.ISm6APeak.NEBvsOriginal.in.NEBHeLa, x = 6, y=12.2, default.units = "cm",width = 5.8, height = 5.8)

dev.off()
