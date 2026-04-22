########### Figure of computational resource of six methods and ISm6APeak ###################
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

#load the peak time, peak CPU, and runtime of six methods

#MACS2
dt.result.MACS2 <- fread("/data2/CRC_m6A/BenchMACS2/MACS2_Ncore20.mem.log",fill = T)
dt.result.MACS2 <- dt.result.MACS2[,1:3]
colnames(dt.result.MACS2) <- c("TimeSecond","CPU","MemoryMb")
#fixed BAM index and depth calculation
fixtime <- max(dt.result.MACS2[CPU>701,TimeSecond])
dt.result.MACS2 <- dt.result.MACS2 %>% filter(TimeSecond>fixtime)
dt.resource.MACS2 <- data.table(RunMinute= max(dt.result.MACS2$TimeSecond)/60 - min(dt.result.MACS2$TimeSecond)/60, PeakCPU=max(dt.result.MACS2$CPU), PeakRAMGB=max(dt.result.MACS2$MemoryMb)/1024)
#TRESS
dt.result.TRESS <- fread("/data2/CRC_m6A/BenchTRESS/TRESS_Ncore20.mem.log",fill = T)
dt.result.TRESS <- dt.result.TRESS[,1:3]
colnames(dt.result.TRESS) <- c("TimeSecond","CPU","MemoryMb")
dt.result.TRESS <- dt.result.TRESS %>% filter(TimeSecond>fixtime)
dt.resource.TRESS <- data.table(RunMinute= max(dt.result.TRESS$TimeSecond)/60 - min(dt.result.TRESS$TimeSecond)/60, PeakCPU=max(dt.result.TRESS$CPU), PeakRAMGB=max(dt.result.TRESS$MemoryMb)/1024)
#exomePeak
dt.result.exomePeak <- fread("/data2/CRC_m6A/BenchexomePeak/exomePeak_Ncore20.mem.log",fill = T)
dt.result.exomePeak <- dt.result.exomePeak[,1:3]
colnames(dt.result.exomePeak) <- c("TimeSecond","CPU","MemoryMb")
dt.result.exomePeak <- dt.result.exomePeak %>% filter(TimeSecond>fixtime)
dt.resource.exomePeak <- data.table(RunMinute= max(dt.result.exomePeak$TimeSecond)/60 - min(dt.result.exomePeak$TimeSecond)/60, PeakCPU=max(dt.result.exomePeak$CPU), PeakRAMGB=max(dt.result.exomePeak$MemoryMb)/1024)
#MeTPeak
dt.result.MeTPeak <- fread("/data2/CRC_m6A/BenchMeTPeak/MeTPeak_Ncore20.mem.log",fill = T)
dt.result.MeTPeak <- dt.result.MeTPeak[,1:3]
colnames(dt.result.MeTPeak) <- c("TimeSecond","CPU","MemoryMb")
dt.result.MeTPeak <- dt.result.MeTPeak %>% filter(TimeSecond>fixtime)
dt.resource.MeTPeak <- data.table(RunMinute= max(dt.result.MeTPeak$TimeSecond)/60 - min(dt.result.MeTPeak$TimeSecond)/60, PeakCPU=max(dt.result.MeTPeak$CPU), PeakRAMGB=max(dt.result.MeTPeak$MemoryMb)/1024)
#exomePeak2
#RAM256GB
dt.result.exomePeak2.RAM256 <- fread("/data2/CRC_m6A/BenchExomePeak2_RAM256/ExomePeak2_RAM256.mem.log",fill = T)
dt.result.exomePeak2.RAM256 <- dt.result.exomePeak2.RAM256[,1:3]
colnames(dt.result.exomePeak2.RAM256) <- c("TimeSecond","CPU","MemoryMb")
dt.result.exomePeak2.RAM256 <- dt.result.exomePeak2.RAM256 %>% filter(TimeSecond>fixtime)
#RAM512
dt.result.exomePeak2.RAM512 <- fread("/data2/CRC_m6A/BenchExomePeak2_RAM512/ExomePeak2_RAM512.mem.log",fill = T)
dt.result.exomePeak2.RAM512 <- dt.result.exomePeak2.RAM512[,1:3]
colnames(dt.result.exomePeak2.RAM512) <- c("TimeSecond","CPU","MemoryMb")
dt.result.exomePeak2.RAM512 <- dt.result.exomePeak2.RAM512 %>% filter(TimeSecond>fixtime)
dt.resource.exomePeak2 <- rbind(data.table(RunMinute= max(dt.result.exomePeak2.RAM256$TimeSecond)/60 - min(dt.result.exomePeak2.RAM256$TimeSecond)/60,
                                           PeakCPU=max(dt.result.exomePeak2.RAM256$CPU), PeakRAMGB=max(dt.result.exomePeak2.RAM256$MemoryMb)/1024) %>% mutate(RAMGBLimit=256),
                                data.table(RunMinute= max(dt.result.exomePeak2.RAM512$TimeSecond)/60 - min(dt.result.exomePeak2.RAM512$TimeSecond)/60,
                                           PeakCPU=max(dt.result.exomePeak2.RAM512$CPU), PeakRAMGB=max(dt.result.exomePeak2.RAM512$MemoryMb)/1024) %>% mutate(RAMGBLimit=512))
#MeRIPtools
#Ncore20
dt.result.MeRIPtools.Ncore20 <- fread("/data2/CRC_m6A/BenchMeRIPtools_Ncore20/MeRIPtools_Ncore20.mem.log",fill = T)
dt.result.MeRIPtools.Ncore20 <- dt.result.MeRIPtools.Ncore20[,1:3]
colnames(dt.result.MeRIPtools.Ncore20) <- c("TimeSecond","CPU","MemoryMb")
dt.result.MeRIPtools.Ncore20 <- dt.result.MeRIPtools.Ncore20 %>% filter(TimeSecond>fixtime)
#Ncore40
dt.result.MeRIPtools.Ncore40 <- fread("/data2/CRC_m6A/BenchMeRIPtools_Ncore40/MeRIPtools_Ncore40.mem.log",fill = T)
dt.result.MeRIPtools.Ncore40 <- dt.result.MeRIPtools.Ncore40[,1:3]
colnames(dt.result.MeRIPtools.Ncore40) <- c("TimeSecond","CPU","MemoryMb")
dt.result.MeRIPtools.Ncore40 <- dt.result.MeRIPtools.Ncore40 %>% filter(TimeSecond>fixtime)
dt.resource.MeRIPtools <- rbind(data.table(RunMinute= max(dt.result.MeRIPtools.Ncore20$TimeSecond)/60 - min(dt.result.MeRIPtools.Ncore20$TimeSecond)/60,
                                           PeakCPU=max(dt.result.MeRIPtools.Ncore20$CPU), PeakRAMGB=max(dt.result.MeRIPtools.Ncore20$MemoryMb)/1024) %>% mutate(CPULimit=20),
                                data.table(RunMinute= max(dt.result.MeRIPtools.Ncore40$TimeSecond)/60 - min(dt.result.MeRIPtools.Ncore40$TimeSecond)/60,
                                           PeakCPU=max(dt.result.MeRIPtools.Ncore40$CPU), PeakRAMGB=max(dt.result.MeRIPtools.Ncore40$MemoryMb)/1024) %>% mutate(CPULimit=40)
                                )
#ISm6APeak
dt.resource.ISm6APeak <- foreach(Mode=c("Default","Sensitive","Strict"), .combine='rbind')%do%{
  foreach(CPULimit=c(20,40),.combine='rbind')%do%{
    foreach(RAMGBLimit=c(256,512),.combine='rbind')%do%{
      res_file <- paste0("/data2/CRC_m6A/BenchNcore",CPULimit,"MaxRAM",RAMGBLimit,"_ISm6APeak_",Mode,"/ISm6APeak_Ncore",CPULimit,"MaxRAM",RAMGBLimit,".mem.log")
      if(file.exists(res_file)){
        dt.result.ISm6APeak <- fread(res_file,fill = T)
        dt.result.ISm6APeak <- dt.result.ISm6APeak[,1:3]
        colnames(dt.result.ISm6APeak) <- c("TimeSecond","CPU","MemoryMb")
        data.table(RunMinute= max(dt.result.ISm6APeak$TimeSecond)/60 - min(dt.result.ISm6APeak$TimeSecond)/60,
                   PeakCPU=max(dt.result.ISm6APeak$CPU), PeakRAMGB=max(dt.result.ISm6APeak$MemoryMb)/1024) %>% mutate(RAMGBLimit=RAMGBLimit, CPULimit=CPULimit,Mode=Mode)
      }else{
        message(paste0("Error!, Can't find file for ",res_file))
      }
    }
  }
}
dt.resource.ISm6APeak <- dt.resource.ISm6APeak %>% mutate(RunMinute = case_when( !(RAMGBLimit==512 & CPULimit== 40) & Mode=="Default" ~ RunMinute+4.130401*60,
                                                                                 !(RAMGBLimit==512 & CPULimit== 40) & Mode=="Sensitive" ~ RunMinute+4.020663*60,
                                                                                 !(RAMGBLimit==512 & CPULimit== 40) & Mode=="Strict" ~ RunMinute+4.212626*60, .default = RunMinute))

#### visualization of six tools (CPU, RAM, Runtime) among four condition (CPU40_RAM512GB, CPU40_RAM256GB, CPU20_RAM512GB, CPU20_RAM256GB)
dt.resource.sixtools <- dplyr::bind_rows(
                             x1=dplyr::bind_rows(dt.resource.MACS2 %>% mutate(Method="MACS2",Mode="Default",CPULimit=20,RAMGBLimit=256),
                                      dt.resource.TRESS %>% mutate(Method="TRESS",Mode="Default",CPULimit=20,RAMGBLimit=256),
                                      dt.resource.exomePeak %>% mutate(Method="exomePeak",Mode="Default",CPULimit=20,RAMGBLimit=256),
                                      dt.resource.MeTPeak %>% mutate(Method="MeTPeak",Mode="Default",CPULimit=20,RAMGBLimit=256),
                                      dt.resource.exomePeak2 %>% filter(RAMGBLimit==256) %>%  mutate(Method="exomePeak2",Mode="Default",CPULimit=20,RAMGBLimit=256),
                                      dt.resource.MeRIPtools %>% filter(CPULimit==20) %>% mutate(Method="MeRIPtools",Mode="Default",CPULimit=20,RAMGBLimit=256)),
                             x2=dplyr::bind_rows(dt.resource.MACS2 %>% mutate(Method="MACS2",Mode="Default",CPULimit=20,RAMGBLimit=512),
                                      dt.resource.TRESS %>% mutate(Method="TRESS",Mode="Default",CPULimit=20,RAMGBLimit=512),
                                      dt.resource.exomePeak %>% mutate(Method="exomePeak",Mode="Default",CPULimit=20,RAMGBLimit=512),
                                      dt.resource.MeTPeak %>% mutate(Method="MeTPeak",Mode="Default",CPULimit=20,RAMGBLimit=512),
                                      dt.resource.exomePeak2 %>% filter(RAMGBLimit==512) %>%  mutate(Method="exomePeak2",Mode="Default",CPULimit=20,RAMGBLimit=512),
                                      dt.resource.MeRIPtools %>% filter(CPULimit==20) %>% mutate(Method="MeRIPtools",Mode="Default",CPULimit=20,RAMGBLimit=512)),
                             x3=dplyr::bind_rows(dt.resource.MACS2 %>% mutate(Method="MACS2",Mode="Default",CPULimit=40,RAMGBLimit=256),
                                      dt.resource.TRESS %>% mutate(Method="TRESS",Mode="Default",CPULimit=40,RAMGBLimit=256),
                                      dt.resource.exomePeak %>% mutate(Method="exomePeak",Mode="Default",CPULimit=40,RAMGBLimit=256),
                                      dt.resource.MeTPeak %>% mutate(Method="MeTPeak",Mode="Default",CPULimit=40,RAMGBLimit=256),
                                      dt.resource.exomePeak2 %>% filter(RAMGBLimit==256) %>%  mutate(Method="exomePeak2",Mode="Default",CPULimit=40,RAMGBLimit=256),
                                      dt.resource.MeRIPtools %>% filter(CPULimit==40) %>% mutate(Method="MeRIPtools",Mode="Default",CPULimit=40,RAMGBLimit=256)),
                             x4=dplyr::bind_rows(dt.resource.MACS2 %>% mutate(Method="MACS2",Mode="Default",CPULimit=40,RAMGBLimit=512),
                                      dt.resource.TRESS %>% mutate(Method="TRESS",Mode="Default",CPULimit=40,RAMGBLimit=512),
                                      dt.resource.exomePeak %>% mutate(Method="exomePeak",Mode="Default",CPULimit=40,RAMGBLimit=512),
                                      dt.resource.MeTPeak %>% mutate(Method="MeTPeak",Mode="Default",CPULimit=40,RAMGBLimit=512),
                                      dt.resource.exomePeak2 %>% filter(RAMGBLimit==512) %>%  mutate(Method="exomePeak2",Mode="Default",CPULimit=40,RAMGBLimit=512),
                                      dt.resource.MeRIPtools %>% filter(CPULimit==40) %>% mutate(Method="MeRIPtools",Mode="Default",CPULimit=40,RAMGBLimit=512)),
                              )
dt.resource.sixtools <- dt.resource.sixtools %>% mutate(CPURAMGroup = paste0("CPU",CPULimit,"_RAM",RAMGBLimit,"Gb") %>% factor(levels=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))) %>%
  group_by(Method) %>% mutate(avg.RunMinute=mean(RunMinute)) %>% as.data.table() %>% dplyr::arrange(avg.RunMinute) %>% mutate(Method=factor(Method,levels=unique(Method)))
p.resource.sixtools.RunMinute  <- ggplot(data=dt.resource.sixtools,
                               aes(x=Method,y=RunMinute, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Running time (Minute)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,max(dt.resource.sixtools$RunMinute),by=60))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
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
p.resource.sixtools.RunMinute
p.resource.sixtools.PeakCPU  <- ggplot(data=dt.resource.sixtools,
                                         aes(x=Method,y=PeakCPU/100, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Peak CPU usage (cores)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,40,by=5))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
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
p.resource.sixtools.PeakCPU
p.resource.sixtools.PeakRAM  <- ggplot(data=dt.resource.sixtools,
                                         aes(x=Method,y=PeakRAMGB, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Peak memory usage (Gb)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,512,by=50))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
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
p.resource.sixtools.PeakRAM
keep_space_drop_legend <- function(p) {
  p +
    theme(
      # legend.position = "right",           
      legend.title = element_blank(),
      legend.text  = element_blank(),
      legend.key   = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0)),
           colour = guide_legend(override.aes = list(alpha = 0)),
           color  = guide_legend(override.aes = list(alpha = 0)),
           shape  = guide_legend(override.aes = list(alpha = 0)),
           linetype = guide_legend(override.aes = list(alpha = 0)))
}
p.resource.sixtools <- cowplot::plot_grid(plotlist = list(keep_space_drop_legend(p.resource.sixtools.RunMinute),
                                                          p.resource.sixtools.PeakCPU,
                                                          keep_space_drop_legend(p.resource.sixtools.PeakRAM)), 
                                          align = "h",axis = "tblr",nrow = 1)

#visualization of ISm6APeak
dt.resource.ISm6APeak <- dt.resource.ISm6APeak %>% mutate(CPURAMGroup = paste0("CPU",CPULimit,"_RAM",RAMGBLimit,"Gb") %>% factor(levels=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))) %>%
  group_by(Method=paste0("ISm6APeak_",Mode)) %>% mutate(avg.RunMinute=mean(RunMinute)) %>% as.data.table() %>% dplyr::arrange(avg.RunMinute) %>% mutate(Method=factor(Method,levels=unique(Method)))
p.resource.ISm6APeak.RunMinute  <- ggplot(data=dt.resource.ISm6APeak,
                                         aes(x=Method,y=RunMinute, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Running time (Minute)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,max(dt.resource.ISm6APeak$RunMinute),by=60))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
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
p.resource.ISm6APeak.RunMinute
p.resource.ISm6APeak.PeakCPU  <- ggplot(data=dt.resource.ISm6APeak,
                                       aes(x=Method,y=PeakCPU/100, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Peak CPU usage (cores)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,40,by=5))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
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
p.resource.ISm6APeak.PeakCPU
p.resource.ISm6APeak.PeakRAM  <- ggplot(data=dt.resource.ISm6APeak,
                                       aes(x=Method,y=PeakRAMGB, fill=CPURAMGroup))+
  geom_bar(position = position_dodge(width = 0.55),stat = "identity",width=0.5)+
  scale_fill_manual(values=c("#BCD6AD","#91CEEC","#F4D8A2","#F09B9A"),breaks=c("CPU20_RAM256Gb","CPU20_RAM512Gb","CPU40_RAM256Gb","CPU40_RAM512Gb"))+
  labs(x=NULL,y="Peak memory usage (Gb)")+
  guides(fill=guide_legend(title = "Resource limit"))+
  scale_y_continuous(expand = expansion(add=c(0,0.1)),breaks = seq(0,512,by=50))+
  # facet_wrap(~IVToverlapped,scales = "free_y")+
  ggpubr::theme_pubr(base_size = 6,base_family = "Helvetica") +
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
p.resource.ISm6APeak.PeakRAM
keep_space_drop_legend <- function(p) {
  p +
    theme(
      # legend.position = "right",           
      legend.title = element_blank(),
      legend.text  = element_blank(),
      legend.key   = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = 0)),
           colour = guide_legend(override.aes = list(alpha = 0)),
           color  = guide_legend(override.aes = list(alpha = 0)),
           shape  = guide_legend(override.aes = list(alpha = 0)),
           linetype = guide_legend(override.aes = list(alpha = 0)))
}
p.resource.ISm6APeak <- cowplot::plot_grid(plotlist = list(keep_space_drop_legend(p.resource.ISm6APeak.RunMinute),
                                                          p.resource.ISm6APeak.PeakCPU,
                                                          keep_space_drop_legend(p.resource.ISm6APeak.PeakRAM)), 
                                          align = "h",axis = "tblr",nrow = 1)

#3D plot to visualization of ISm6APeak vs six tools  computational resource
breaks <- c("ISm6APeak","MACS2","exomePeak","exomePeak2","MeTPeak","TRESS","MeRIPtools")
values <- c("#eb4601", "#e9abac", "#7bb0d5", "#cfe4b6", "#cab0d2", "#f6c780", "#84cdc0")
col_map <- setNames(values, breaks)


dt.resource.comparison <- bind_rows(dt.resource.sixtools %>% mutate(PeakRAM=PeakRAMGB,RunTime=RunMinute) %>% dplyr::select(Method,RunTime,PeakCPU,PeakRAM,CPURAMGroup),
                                    dt.resource.ISm6APeak %>% dplyr::filter(Method=="ISm6APeak_Default") %>% mutate(Method="ISm6APeak",PeakRAM=PeakRAMGB,RunTime=RunMinute) %>%  dplyr::select(Method,RunTime,PeakCPU,PeakRAM,CPURAMGroup)) %>%
                          mutate(Method = factor(Method, levels = breaks),Color  = unname(col_map[as.character(Method)]),Label  = as.character(Method))
#scaled to max of sixtools
dt.resource.comparison <- dt.resource.comparison %>% group_by(CPURAMGroup) %>% mutate(RunTime.scaled=RunTime/max(RunTime[Method!="ISm6APeak"]),
                                                                                      PeakCPU.scaled=PeakCPU/max(PeakCPU[Method!="ISm6APeak"]),
                                                                                      PeakRAM.scaled=PeakRAM/max(PeakRAM[Method!="ISm6APeak"])) %>% as.data.table()
#function to make dropline to three trace
make_dropline_trace <- function(df, xcol="RunTime", ycol="PeakCPU", zcol="PeakRAM",
                                x0=NULL, y0=NULL, z0=NULL) {
  x <- df[[xcol]]; y <- df[[ycol]]; z <- df[[zcol]]
  if (is.null(x0)) x0 <- min(x, na.rm = TRUE)
  if (is.null(y0)) y0 <- min(y, na.rm = TRUE)
  if (is.null(z0)) z0 <- min(z, na.rm = TRUE)
  
  X <- Y <- Z <- c()
  for (i in seq_along(x)) {
    # 线1：投影到 z=z0 (xy平面)
    X <- c(X, x[i], x[i], NA);  Y <- c(Y, y[i], y[i], NA);  Z <- c(Z, z[i], z0,   NA)
    # 线2：投影到 y=y0 (xz平面)
    X <- c(X, x[i], x[i], NA);  Y <- c(Y, y[i], y0,   NA);  Z <- c(Z, z[i], z[i], NA)
    # 线3：投影到 x=x0 (yz平面)
    X <- c(X, x[i], x0,   NA);  Y <- c(Y, y[i], y[i], NA);  Z <- c(Z, z[i], z[i], NA)
  }
  list(X=X, Y=Y, Z=Z, x0=x0, y0=y0, z0=z0)
}
#3D plot
font_family <- "Arial"
axis_title_size <- 10
axis_tick_size  <- 10
label_size      <- 10   # 点标签大小

for(g in unique(dt.resource.comparison$CPURAMGroup)[c(1,4)]){
  df <- dt.resource.comparison %>% dplyr::filter(CPURAMGroup==g) %>% dplyr::select(-RunTime,-PeakCPU, -PeakRAM) %>%
    dplyr::rename(RunTime=RunTime.scaled, PeakCPU=PeakCPU.scaled,PeakRAM=PeakRAM.scaled)
  # 投影线落到“各轴最小值平面”
  dl <- make_dropline_trace(df, "RunTime","PeakCPU","PeakRAM")
  
  p3d <- plot_ly() %>%
    # (A) 投影线：统一灰色
    add_trace(
      data = NULL, inherit = FALSE,
      type = "scatter3d", mode = "lines",
      x = dl$X, y = dl$Y, z = dl$Z,
      showlegend = FALSE,
      line = list(color = "grey8", width = 2)
    ) %>%
    # (B) 点
    add_trace(
      data = df,
      type = "scatter3d", mode = "markers+text",
      x = ~RunTime, y = ~PeakCPU, z = ~PeakRAM,
      text = ~Label,
      textposition = "top center",
      textfont = list(family = font_family, size = label_size),
      marker = list(size = 6, color = ~Color),
      showlegend = FALSE
    ) %>%
    layout(
      showlegend = FALSE,
      scene = list(
        xaxis = list(
          title = "Normalized RunTime",
          titlefont = list(family = font_family, size = axis_title_size),
          tickfont  = list(family = font_family, size = axis_tick_size),
          zeroline = FALSE
        ),
        yaxis = list(
          title = "Normalized PeakCPU",
          titlefont = list(family = font_family, size = axis_title_size),
          tickfont  = list(family = font_family, size = axis_tick_size),
          zeroline = FALSE
        ),
        zaxis = list(
          title = "Normalized PeakRAM",
          titlefont = list(family = font_family, size = axis_title_size),
          tickfont  = list(family = font_family, size = axis_tick_size),
          zeroline = FALSE
        )
      ),
      margin = list(l = 0, r = 0, b = 0, t = 0)
    )
  p3d
  
  # views <- list(
  #   front = list(x = 1.7, y = 1.7, z = 1.2),
  #   side  = list(x = 2.2, y = 0.2, z = 1.0),
  #   top   = list(x = 0.1, y = 0.1, z = 3.0)
  # )
  # #save pdf
  # for (nm in names(views)[1]) {
  #   p_tmp <- p3d %>% layout(scene = list(camera = list(eye = views[[nm]])))
  #   
  #   plotly::save_image(
  #     p_tmp,
  #     file  = paste0(g,"_3d_scatter_", nm, ".pdf"),
  #     # width = w_px,
  #     # height = h_px,
  #     width=400*1.3,
  #     height=300*1.3,
  #     scale = 0.8
  #   )
  # }
  cm_to_px <- function(cm, dpi) round(cm / 2.54 * dpi)
  
  # 目标版面尺寸（cm）
  w_cm <- 10
  h_cm <- 8
  
  # 投稿建议：600dpi + scale=2（等效~1200dpi）
  dpi   <- 600
  scale <- 2
  
  w_px <- cm_to_px(w_cm, dpi)   # ~2362
  h_px <- cm_to_px(h_cm, dpi)   # ~1890
  
  views <- list(
    front = list(x = 1.7, y = 1.7, z = 1.2),
    side  = list(x = 2.2, y = 0.2, z = 1.0),
    top   = list(x = 0.1, y = 0.1, z = 3.0)
  )
  
  for (nm in names(views)[1]) {
    p_tmp <- p3d %>% layout(scene = list(camera = list(eye = views[[nm]])))
    
    plotly::save_image(
      p_tmp,
      file   = paste0(g,"_3d_scatter_", nm, "_10x8cm.pdf"),
      width  = w_px,
      height = h_px,
      scale  = scale
    )
  }
  
}



#solution use scatterplot3d
pdf("Figure_S_3D_usability.pdf", width = 8.5, height = 7, bg = "white")  # ← Vector PDF
for(g in unique(dt.resource.comparison$CPURAMGroup)[c(1,4)]){
  df <- dt.resource.comparison %>% dplyr::filter(CPURAMGroup == g)
  s3d <- scatterplot3d(
    x        = df$PeakRAM.scaled,      # x-axis
    y        = df$PeakCPU.scaled,      # y-axis
    z        = df$RunTime.scaled,  # z-axis (vertical)
    color    = df$Color,
    pch      = 19,
    cex.symbols = 2,
    type     = "h",                  # blue drop lines (exactly like your figure)
    lty.hplot = 1,
    angle    = 45,                   # ← change to 30–60 to match your view angle
    scale.y  = 1,
    xlab     = "Normalized peak memory",
    ylab     = "Normalized peak CPU usage",
    zlab     = "Normalized running time",
    grid     = TRUE,
    box      = TRUE,
    ticktype = "detailed"
  )
  
  # Add tool name labels
  text(s3d$xyz.convert(df$PeakRAM.scaled, df$PeakCPU.scaled, df$RunTime.scaled),
       labels = df$Method,
       cex = 0.95,
       pos = 3,                     # 1=below, 2=left, 3=above, 4=right
       col = "black")
}
dev.off()


pdf("FigureS_computational_resource_of_six_tools_vs_ISm6APeak.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18, height = 24, default.units = "cm",showGuides = F)
plotGG(plot = p.resource.sixtools, x = 0.2, y=0.05, default.units = "cm",width = 17.8, height = 7)
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.resource.ISm6APeak, x = 0.2, y=7.2, default.units = "cm",width = 17.8, height = 7)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05, y = 7.2,just = c("top","left"), default.units = "cm",draw=T)


dev.off()

