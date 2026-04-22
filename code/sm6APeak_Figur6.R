

#### Figure 6 Bioinformatic and experimental validation revealed Hyper-methylated driven activation of PTRH1 and C5AR1 were associated with tumor promotive roles in CRC #########

#a. Dysregulated and prognostic association of novel DMEGs identified by sm6APeak highlight potential importance of  C5AR1 and PTRH1
#the DEG of PTRH1 and C5AR1 among all CRC cohort
#load the DEG result of PRTH1 and C5AR1
dt.TumorDEG.C5AR1.PRTH1.CRC <- readRDS("/data2/CRC_m6A/Analysis/dt.TumorDEG.C5AR1.PRTH1.RDS")
#for each gene, select the top3 most DEG dataset
dt.TumorDEG.C5AR1.PRTH1.CRC <- dt.TumorDEG.C5AR1.PRTH1.CRC %>% filter(LFC>0 & FDR<0.05)
unique(dt.TumorDEG.C5AR1.PRTH1.CRC$StudyID)
# "GSE50760"  "GSE50421"  "GSE134525" "TCGA_COAD" "GSE113513" "TCGA_READ" "GSE50117"  "GSE103512" "GSE110224" "GSE41657"  "GSE156355"
dt.top2.TumorDGE.C5AR1.PTRH1 <- dt.TumorDEG.C5AR1.PRTH1.CRC %>% group_by(Symbol) %>% dplyr::arrange(desc(LFC)) %>% slice_head(n=2) %>% as.data.table()
dt.top2.TumorDGE.C5AR1.PTRH1
#    Direction Symbol                Pvalue                  FDR      LFC   StudyID
# 1:        Up  C5AR1 0.0000043715411802511 0.000024165980576109 1.239100  GSE50760
# 2:        Up  C5AR1 0.0001371930524394302 0.002344334161578928 1.173712  GSE50421
# 4:        Up  PTRH1 0.0000553894177880381 0.000266966161724311 2.435079 TCGA_READ
# 5:        Up  PTRH1 0.0000000000002359187 0.000000000001327596 2.400279 TCGA_COAD

#load C5AR1 and PTRH1 expr for DEG visualization
C5AR1.PTRH1.expr.files <- list.files(path = "/data2/CRC_m6A/Analysis/",pattern = "^(dt).+(C5AR1).+(PTRH1).+(RDS)$",full.names = T)
names(C5AR1.PTRH1.expr.files) <- list.files(path = "/data2/CRC_m6A/Analysis/",pattern = "^(dt).+(C5AR1).+(PTRH1).+(RDS)$",full.names = F) %>% strsplit(split=".",fixed=T) %>% sapply("[",2)

#C5AR1 in GSE50760 and GSE50421
dt.DEG.expr.C5AR1.GSE50760 <- readRDS(C5AR1.PTRH1.expr.files["GSE50760"]) %>% filter(Symbol=="C5AR1") %>% mutate(logExpr=log2(FPM)) %>%
  mutate(Label2=sprintf("FDR = %.2g",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="C5AR1" & StudyID=="GSE50760",FDR]),
         Label1=sprintf("LFC = %.2f",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="C5AR1" & StudyID=="GSE50760",LFC]),
         Condition=condition2)
dt.DEG.expr.C5AR1.GSE50421 <- readRDS(C5AR1.PTRH1.expr.files["GSE50421"]) %>% filter(Gene.symbol=="C5AR1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="C5AR1" & StudyID=="GSE50421",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="C5AR1" & StudyID=="GSE50421",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.C5AR1.GSE50421 <- dt.DEG.expr.C5AR1.GSE50421 %>% inner_join(x=.,y=dt.DEG.expr.C5AR1.GSE50421 %>% filter(logFC>0) %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.C5AR1.GSE50760 <- dt.DEG.expr.C5AR1.GSE50760 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.C5AR1.GSE50760$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of C5AR1",title="GSE50760", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.C5AR1.GSE50760$Label1)),caption  = unique(dt.DEG.expr.C5AR1.GSE50760$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.C5AR1.GSE50421<- dt.DEG.expr.C5AR1.GSE50421 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.C5AR1.GSE50421$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("Wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of C5AR1",title="GSE50421", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.C5AR1.GSE50421$Label1)),caption = unique(dt.DEG.expr.C5AR1.GSE50421$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.DEG.expr.C5AR1.GSE50760, x = 0.2, y=4.9, default.units = "cm",width = 2.5, height = 4.2)
plotGG(plot = p.DEG.expr.C5AR1.GSE50421, x = 2.8, y=4.9, default.units = "cm",width = 2.5, height = 4.2)


#PTRH1 in TCGA_COAD and TCGA READ
dt.DEG.expr.PTRH1.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="PTRH1") %>% mutate(logExpr=log2(tpm_unstranded),Condition=SampleCondition) %>%
  dplyr::select(Symbol=gene_name,logExpr,Condition,SampleID) %>% distinct() %>% mutate(Condition=factor(Condition, levels=c("PrimaryTumor", "SolidTissueNormal"),labels=c("Tumor","Adjacent"))) %>%
  mutate(Label2=sprintf("FDR = %.2g",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="PTRH1" & StudyID=="TCGA_COAD",FDR]),
         Label1=sprintf("LFC = %.2f",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="PTRH1" & StudyID=="TCGA_COAD",LFC])) %>% filter(!is.na(Condition))
dt.DEG.expr.PTRH1.TCGA_READ <- readRDS(C5AR1.PTRH1.expr.files["TCGA_READ"]) %>% filter(gene_name=="PTRH1") %>% mutate(logExpr=log2(tpm_unstranded),Condition=SampleCondition) %>%
  dplyr::select(Symbol=gene_name,logExpr,Condition,SampleID) %>% distinct() %>% mutate(Condition=factor(Condition, levels=c("PrimaryTumor", "SolidTissueNormal"),labels=c("Tumor","Adjacent"))) %>%
  mutate(Label2=sprintf("FDR = %.2g",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="PTRH1" & StudyID=="TCGA_READ",FDR]),
         Label1=sprintf("LFC = %.2f",dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="PTRH1" & StudyID=="TCGA_READ",LFC])) %>% filter(!is.na(Condition))

p.DEG.expr.PTRH1.TCGA_READ <- dt.DEG.expr.PTRH1.TCGA_READ %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.TCGA_READ$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="TCGA_READ", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.TCGA_READ$Label1)), caption = unique(dt.DEG.expr.PTRH1.TCGA_READ$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.TCGA_COAD <- dt.DEG.expr.PTRH1.TCGA_COAD %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.TCGA_COAD$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("Wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="TCGA_COAD", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.TCGA_COAD$Label1)), caption = unique(dt.DEG.expr.PTRH1.TCGA_COAD$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.DEG.expr.PTRH1.TCGA_READ, x = 0.2, y=9.2, default.units = "cm",width = 2.5, height = 4.2)
plotGG(plot = p.DEG.expr.PTRH1.TCGA_COAD, x = 2.8, y=9.2, default.units = "cm",width = 2.5, height = 4.2)


#the stage DGE of PTRH1 and C5AR1 in TCGA_COAD and SidraLUMC2022
#PTRH1 in TCGA_COAD
dt.stageDEG.TCGA_COAD <- readRDS("/data/COAD/Analysis//dt.DEG.TCGA_COAD.StageDEG.RDS")
dt.stageDEG.TCGA_COAD[Symbol=="PTRH1" & StudyID=="TCGA_COAD_StageDEG",]
# GeneID Symbol       LFC      pvalue       padj            StudyID
# 1: ENSG00000187024.15  PTRH1 0.7680169 0.002510133 0.03720088 TCGA_COAD_StageDEG
dt.StageDEG.expr.PTRH1.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="PTRH1" & SampleCondition=="PrimaryTumor") %>% mutate(logExpr=log2(tpm_unstranded+0.01)) %>%
  dplyr::select(Symbol=gene_name,logExpr,ajcc_stage,SampleID) %>% distinct() %>% mutate(Condition=dplyr::case_match(ajcc_stage,
                                                                                                                   c("I","IA") ~ "I",
                                                                                                                   c("II","IIA","IIB","IIC") ~ "II",
                                                                                                                   c("III", "IIIA","IIIB","IIIC") ~ "III",
                                                                                                                   c("IV","IVA","IVB") ~ "IV")) %>%
  mutate(Condition=factor(Condition,levels=c("I","II","III","IV"))) %>%
  mutate(Label1=sprintf("LFC = %.2f",dt.stageDEG.TCGA_COAD[Symbol=="PTRH1" & StudyID=="TCGA_COAD_StageDEG",LFC]),
         Label2=sprintf("Pvalue = %.2g",dt.stageDEG.TCGA_COAD[Symbol=="PTRH1" & StudyID=="TCGA_COAD_StageDEG",pvalue])) %>% filter(!is.na(Condition) & !is.na(logExpr))
table(dt.StageDEG.expr.PTRH1.TCGA_COAD$Condition)
p.StageDEG.expr.PTRH1.TCGA_COAD <- dt.StageDEG.expr.PTRH1.TCGA_COAD %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent")+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.StageDEG.expr.PTRH1.TCGA_COAD$Condition)[c(1,4)],2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))+
  # coord_cartesian(ylim = c(-6.5,-1))+
  scale_color_manual(values=c("#a7c9e7","#feee84","#f9b56d","#f48485"),breaks = c("I","II","III","IV"),guide="none")+
  labs(x="Tumor stage", y="log2 expression of PTRH1",title="TCGA_COAD", subtitle = paste0("ProgressionDEG ",unique(dt.StageDEG.expr.PTRH1.TCGA_COAD$Label1)), caption = unique(dt.StageDEG.expr.PTRH1.TCGA_COAD$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor = element_line(linewidth = 0.6,linetype = 4)
  )

#C5AR1 progression DEG in SidraLUMC
dt.stageDEG.SidraLUMC <- readRDS("/data/COAD/Analysis//dt.DEG.SidraLUMC2022.StageDEG.RDS")
dt.stageDEG.SidraLUMC[Symbol=="C5AR1" & StudyID=="SidraLUMC2022_StageDEG",]
#   Symbol      padj      pvalue       LFC                StudyID
# 1:  C5AR1 0.1116859 0.009490311 0.6273614 SidraLUMC2022_StageDEG
dt.StageDEG.expr.C5AR1.SidraLUMC <- readRDS(C5AR1.PTRH1.expr.files["SidraLUMC2022"]) %>% filter(Symbol=="C5AR1") %>% mutate(logExpr=fpkm) %>%
  dplyr::select(Symbol,logExpr,ajcc_stage=ajcc_path_stage,clinic_stage,SampleID) %>% distinct() %>%
  mutate(Condition=dplyr::case_match(ajcc_stage %>% as.character(),
                                     c("I","IA","1") ~ "I",
                                     c("II","IIA","IIB","IIC","2") ~ "II",
                                     c("III", "IIIA","IIIB","IIIC","3") ~ "III",
                                     c("IV","IVA","IVB","4") ~ "IV") %>% factor(levels = c("I","II","III","IV"))) %>%
  # mutate(Condition=dplyr::case_match(clinic_stage %>% as.character(),
  #                                    c("I","IA","1") ~ "I",
  #                                    c("II","IIA","IIB","IIC","2") ~ "II",
  #                                    c("III", "IIIA","IIIB","IIIC","3") ~ "III",
  #                                    c("IV","IVA","IVB","4") ~ "IV") %>% factor(levels = c("I","II","III","IV"))) %>%
  # mutate(Condition=factor(Condition,levels=c("I","II","III","IV"))) %>%
  mutate(Label1=sprintf("LFC = %.2f",dt.stageDEG.SidraLUMC[Symbol=="C5AR1" & StudyID=="SidraLUMC2022_StageDEG",LFC]),
         Label2=sprintf("Pvalue = %.2f",dt.stageDEG.SidraLUMC[Symbol=="C5AR1" & StudyID=="SidraLUMC2022_StageDEG",pvalue])) %>% filter(!is.na(Condition) & !is.na(logExpr))
table(dt.StageDEG.expr.C5AR1.SidraLUMC$Condition)
# I  II III  IV
# 55 122 110  61
p.StageDEG.expr.C5AR1.SidraLUMC<- dt.StageDEG.expr.C5AR1.SidraLUMC %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent")+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.StageDEG.expr.C5AR1.SidraLUMC$Condition)[c(1,4)],2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("Wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))+
  scale_color_manual(values=c("#a7c9e7","#feee84","#f9b56d","#f48485"),breaks = c("I","II","III","IV"),guide="none")+
  # coord_cartesian(ylim = c(5,12))+
  labs(x="Tumor stage", y="log2 expression of C5AR1",title="SidraLUMC", subtitle = paste0("ProgressionDEG ",unique(dt.StageDEG.expr.C5AR1.SidraLUMC$Label1)), caption = unique(dt.StageDEG.expr.C5AR1.SidraLUMC$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.StageDEG.expr.C5AR1.SidraLUMC

# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.StageDEG.expr.PTRH1.TCGA_COAD, x = 5.4, y=4.9, default.units = "cm",width = 4.2, height = 4.2)
plotGG(plot = p.StageDEG.expr.C5AR1.SidraLUMC, x = 5.4, y=9.2, default.units = "cm",width = 4.2, height = 4.2)

#the survival association of PTRH1 and C5AR1 among all CRC cohort
dt.SurvivalGene.C5AR1.PTRH1.CRC <- readRDS("/data2/CRC_m6A/Analysis/dt.SurvivalGene.res.C5AR1.PRTH1.RDS") %>%
  dplyr::filter(Pvalue<=0.05 & Direction=="Risk") %>% mutate(SurvivalMetricClass=factor(SurvivalMetric, levels=c("OS","PFS","PFI","DFS","DSS"), labels=c("OS", rep("PFS",4)))) %>%
  dplyr::arrange(Symbol,SurvivalMetricClass,Pvalue)
unique(dt.SurvivalGene.C5AR1.PTRH1.CRC$StudyID)
# "CPTAC2_Colon"   "TCGA_COAD"      "TCGA_READ"      "SidraLUMC2022"  "Rectal_MSK2022"
dt.SurvivalGene.C5AR1.PTRH1.CRC %>% group_by(Symbol,SurvivalMetricClass) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table()
#    SurvivalMetric Direction Symbol    Expr   Pvalue      StudyID SurvivalMetricClass
# 1:             OS      Risk  C5AR1 fpkm_uq 0.000024 CPTAC2_Colon                  OS
# 2:            PFI      Risk  C5AR1 fpkm_uq 0.000840    TCGA_COAD                 PFS
# 3:             OS      Risk  PTRH1    fpkm 0.004900    TCGA_COAD                  OS
# 4:            PFI      Risk  PTRH1     tpm 0.000620    TCGA_READ                 PFS

#PTRH1 OS in TCGA_COAD
dt.SurvivalGene.PTRH1.OS.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="PTRH1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_unstranded,SampleID,OS,OS.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="TCGA_COAD" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=OS.time/30)
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.OS.TCGA_COAD, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   0.0379  2.461373
dt.SurvivalGene.PTRH1.OS.TCGA_COAD <- dt.SurvivalGene.PTRH1.OS.TCGA_COAD %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.TCGA_COAD)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.OS.TCGA_COAD.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.TCGA_COAD)
p.PTRH1.OS.TCGA_COAD <- ggsurvplot(PTRH1.OS.TCGA_COAD.fit, data = dt.SurvivalGene.PTRH1.OS.TCGA_COAD,
                   risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                   pval = TRUE, pval.method = T,pval.size=2,
                   censor.size=2,size=0.5,
                   legend.title="PTRH1",legend=c(0.8,0.8),legend.labs=c("low","high"),
                   palette = c("#80bcb0",  "#e9786a"),
                   break.time.by =12,
                   ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                   theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                           legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                           legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                           strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="TCGA_COAD")
# p.PTRH1.OS.TCGA_COAD$plot <- p.PTRH1.OS.TCGA_COAD$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.OS.TCGA_COAD$table <- p.PTRH1.OS.TCGA_COAD$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.OS.TCGA_COAD <- cowplot::plot_grid(p.PTRH1.OS.TCGA_COAD$plot, p.PTRH1.OS.TCGA_COAD$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.OS.TCGA_COAD <- p.PTRH1.OS.TCGA_COAD$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.OS.TCGA_COAD
# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.PTRH1.OS.TCGA_COAD, x = 10, y=4.9, default.units = "cm",width = 4.2, height = 4.2)

#PTRH1 PFS in TCGA_READ
dt.SurvivalGene.PTRH1.PFS.TCGA_READ <- readRDS(C5AR1.PTRH1.expr.files["TCGA_READ"]) %>% filter(gene_name=="PTRH1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_unstranded,SampleID,PFS=PFI,PFS.time=PFI.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="TCGA_READ" & SurvivalMetric=="PFI" & Direction=="Risk",Pvalue])) %>%
  mutate(PFS.time=PFS.time/30)
#determine the optimal cutoff for PFS
require(survminer)
require(survival)
PFS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.PFS.TCGA_READ, time="PFS.time",event = "PFS",variables = c("Expr"))
PFS.dt.cutpoint <- summary(PFS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
PFS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   0.0094  3.239751
dt.SurvivalGene.PTRH1.PFS.TCGA_READ <- dt.SurvivalGene.PTRH1.PFS.TCGA_READ %>% mutate(ExprGroup=if_else(Expr > PFS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
PFS.logrank.res <- survdiff(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.PFS.TCGA_READ)
PFS.logrank.p = 1 - pchisq(PFS.logrank.res$chisq, length(PFS.logrank.res$n) - 1)
PFS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.PFS.TCGA_READ.fit <- survfit(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.PFS.TCGA_READ)
p.PTRH1.PFS.TCGA_READ <- ggsurvplot(PTRH1.PFS.TCGA_READ.fit, data = dt.SurvivalGene.PTRH1.PFS.TCGA_READ,
                                   risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                   pval = TRUE, pval.method = T,pval.size=2,
                                   censor.size=2,size=0.5,
                                   legend.title="PTRH1",legend=c(0.8,0.7),legend.labs=c("low","high"),
                                   palette = c("#80bcb0",  "#e9786a"),
                                   break.time.by =12,
                                   ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                     theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                           legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                           legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                           strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="PFS probability",x="Survival time (months)", title="TCGA_READ")
# p.PTRH1.PFS.TCGA_READ$plot <- p.PTRH1.PFS.TCGA_READ$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.PFS.TCGA_READ$table <- p.PTRH1.PFS.TCGA_READ$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.PFS.TCGA_READ <- cowplot::plot_grid(p.PTRH1.PFS.TCGA_READ$plot, p.PTRH1.PFS.TCGA_READ$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.PFS.TCGA_READ <- p.PTRH1.PFS.TCGA_READ$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.PFS.TCGA_READ
# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.PTRH1.PFS.TCGA_READ , x = 14.3, y=4.9, default.units = "cm",width = 4.2, height = 4.2)

#C5AR1 OS in CPTAC2_Colon
dt.SurvivalGene.C5AR1.OS.CPTAC2 <- readRDS(C5AR1.PTRH1.expr.files["CPTAC2_Colon"]) %>% filter(Symbol=="C5AR1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol,Expr=fpkm_uq,SampleID=Case_ID,OS=os_status,OS.time=os_months) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="CPTAC2_Colon" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=as.numeric(OS.time),OS=case_match(OS,"0:LIVING" ~ 0,"1:DECEASED" ~ 1)) %>% filter(!is.na(OS.time) & !is.na(OS))
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.OS.CPTAC2, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   8.7017  4.146039
dt.SurvivalGene.C5AR1.OS.CPTAC2 <- dt.SurvivalGene.C5AR1.OS.CPTAC2 %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.CPTAC2)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.OS.CPTAC2.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.CPTAC2)
p.C5AR1.OS.CPTAC2 <- ggsurvplot(C5AR1.OS.CPTAC2.fit, data = dt.SurvivalGene.C5AR1.OS.CPTAC2,
                                   risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                   pval = TRUE, pval.method = T,pval.size=2,
                                   censor.size=2,size=0.5,
                                   legend.title="C5AR1",legend=c(0.8,0.4),legend.labs=c("low","high"),
                                   palette = c("#80bcb0",  "#e9786a"),
                                   break.time.by =12,
                                   ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                     theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                           legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                           legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                           strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="CPTAC2_Colon")
# p.C5AR1.OS.CPTAC2$plot <- p.C5AR1.OS.CPTAC2$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.OS.CPTAC2$table <- p.C5AR1.OS.CPTAC2$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.OS.CPTAC2 <- cowplot::plot_grid(p.C5AR1.OS.CPTAC2$plot, p.C5AR1.OS.CPTAC2$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.OS.CPTAC2 <- p.C5AR1.OS.CPTAC2$plot
p.C5AR1.OS.CPTAC2
# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.C5AR1.OS.CPTAC2$plot, x = 10, y=4.9, default.units = "cm",width = 4.2, height = 4.2)



#C5AR1 PFS in TCGA_COAD
dt.SurvivalGene.C5AR1.PFS.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="C5AR1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_unstranded,SampleID,PFS=PFI,PFS.time=PFI.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="TCGA_COAD" & SurvivalMetric=="PFI" & Direction=="Risk",Pvalue])) %>%
  mutate(PFS.time=PFS.time/30)
#determine the optimal cutoff for PFS
require(survminer)
require(survival)
PFS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.PFS.TCGA_COAD, time="PFS.time",event = "PFS",variables = c("Expr"))
PFS.dt.cutpoint <- summary(PFS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
PFS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   0.0094  3.239751
dt.SurvivalGene.C5AR1.PFS.TCGA_COAD <- dt.SurvivalGene.C5AR1.PFS.TCGA_COAD %>% mutate(ExprGroup=if_else(Expr > PFS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
PFS.logrank.res <- survdiff(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.PFS.TCGA_COAD)
PFS.logrank.p = 1 - pchisq(PFS.logrank.res$chisq, length(PFS.logrank.res$n) - 1)
PFS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.PFS.TCGA_COAD.fit <- survfit(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.PFS.TCGA_COAD)
p.C5AR1.PFS.TCGA_COAD <- ggsurvplot(C5AR1.PFS.TCGA_COAD.fit, data = dt.SurvivalGene.C5AR1.PFS.TCGA_COAD,
                                    risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                    pval = TRUE, pval.method = T,pval.size=2,
                                    censor.size=2,size=0.5,
                                    legend.title="C5AR1",legend=c(0.8,0.8),legend.labs=c("low","high"),
                                    palette = c("#80bcb0",  "#e9786a"),
                                    break.time.by =12,
                                    ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                      theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                            legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                            legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                            strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="PFS probability",x="Survival time (months)", title="TCGA_COAD")
# p.C5AR1.PFS.TCGA_COAD$plot <- p.C5AR1.PFS.TCGA_COAD$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.PFS.TCGA_COAD$table <- p.C5AR1.PFS.TCGA_COAD$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.PFS.TCGA_COAD <- cowplot::plot_grid(p.C5AR1.PFS.TCGA_COAD$plot, p.C5AR1.PFS.TCGA_COAD$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.PFS.TCGA_COAD <- p.C5AR1.PFS.TCGA_COAD$plot+coord_cartesian(xlim = c(0,96))
p.C5AR1.PFS.TCGA_COAD
# pageCreate(width = 8.1, height =20, default.units = "cm",showGuides = F)
plotGG(plot = p.C5AR1.PFS.TCGA_COAD, x = 14.3, y=4.9, default.units = "cm",width = 4.2, height = 4.2)


save.image("sm6APeak_Figure6_intermediate.results.RDS")
fig6.list <- list(p.DEG.expr.PTRH1.TCGA_READ=p.DEG.expr.PTRH1.TCGA_READ,
                  p.DEG.expr.PTRH1.TCGA_COAD=p.DEG.expr.PTRH1.TCGA_COAD,
                  p.StageDEG.expr.PTRH1.TCGA_COAD=p.StageDEG.expr.PTRH1.TCGA_COAD,
                  p.PTRH1.PFS.TCGA_READ=p.PTRH1.PFS.TCGA_READ,
                  p.PTRH1.OS.TCGA_COAD=p.PTRH1.OS.TCGA_COAD,
                  p.DEG.expr.C5AR1.GSE50760=p.DEG.expr.C5AR1.GSE50760,
                  p.DEG.expr.C5AR1.GSE50421=p.DEG.expr.C5AR1.GSE50421,
                  p.StageDEG.expr.C5AR1.SidraLUMC=p.StageDEG.expr.C5AR1.SidraLUMC,
                  p.C5AR1.PFS.TCGA_COAD=p.C5AR1.PFS.TCGA_COAD,
                  p.C5AR1.OS.CPTAC2
)

saveRDS(fig6.list, file="Figure6.plot.list.RDS")


#combine four plots together
pdf("Figure6_Activation_of_PTRH1_and_C5AR1_associatied_with_tumor_promoting_roles_in_CRC.pdf",width = 8.2677, height = 11.693)
pageCreate(width = 18.2, height =24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.TCGA_READ, x = 0.05, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotGG(plot = p.DEG.expr.PTRH1.TCGA_COAD, x = 2.6, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 5.2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.StageDEG.expr.PTRH1.TCGA_COAD, x = 5.2, y=0.2, default.units = "cm",width = 4.1, height = 4.4)
plotGG(plot = p.PTRH1.PFS.TCGA_READ, x = 9.5, y=0.2, default.units = "cm",width = 4.2, height = 4.2)
plotGG(plot = p.PTRH1.OS.TCGA_COAD, x = 9.5+4.2+0.1, y=0.2, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 9.5, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)

#row2
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05, y = 4.5,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.C5AR1.GSE50760, x = 0.05, y=4.7, default.units = "cm",width = 2.5, height = 4.4)
plotGG(plot = p.DEG.expr.C5AR1.GSE50421, x = 2.6, y=4.7, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 5.2, y = 4.5,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.StageDEG.expr.C5AR1.SidraLUMC, x = 5.2, y=4.7, default.units = "cm",width = 4.1, height = 4.4)
plotGG(plot = p.C5AR1.PFS.TCGA_COAD, x = 9.5, y=4.7, default.units = "cm",width = 4.2, height = 4.2)
plotGG(plot = p.C5AR1.OS.CPTAC2, x = 9.5+4.2+0.1, y=4.7, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 9.5, y = 4.5,just = c("top","left"), default.units = "cm",draw=T)

dev.off()

#b. Tumor migration, invasion, proliferation, and CDX tumor growth association of PTRH1 and C5AR1


########################## S6_1 In silico validation of sm6APeak revealed DMEG ########################
######################### S6_2 Full results of Tumor DEG results and prognostic association of PTRH1 and C5AR1 #####################
dt.TumorDEG.C5AR1.PRTH1.CRC <- readRDS("/data2/CRC_m6A/Analysis/dt.TumorDEG.C5AR1.PRTH1.RDS")
#for each gene, select the top3 most DEG dataset
dt.TumorDEG.C5AR1.PRTH1.CRC <- dt.TumorDEG.C5AR1.PRTH1.CRC %>% filter(LFC>0 & FDR<0.05)
unique(dt.TumorDEG.C5AR1.PRTH1.CRC$StudyID)
# "GSE50760"  "GSE50421"  "GSE134525" "TCGA_COAD" "GSE113513" "TCGA_READ" "GSE50117"  "GSE103512" "GSE110224" "GSE41657"  "GSE156355"
dt.top2.TumorDGE.C5AR1.PTRH1 <- dt.TumorDEG.C5AR1.PRTH1.CRC %>% group_by(Symbol) %>% dplyr::arrange(desc(LFC)) %>% slice_head(n=2) %>% as.data.table()
dt.top2.TumorDGE.C5AR1.PTRH1
#    Direction Symbol                Pvalue                  FDR      LFC   StudyID
# 1:        Up  C5AR1 0.0000043715411802511 0.000024165980576109 1.239100  GSE50760
# 2:        Up  C5AR1 0.0001371930524394302 0.002344334161578928 1.173712  GSE50421
# 4:        Up  PTRH1 0.0000553894177880381 0.000266966161724311 2.435079 TCGA_READ
# 5:        Up  PTRH1 0.0000000000002359187 0.000000000001327596 2.400279 TCGA_COAD

#load C5AR1 and PTRH1 expr for DEG visualization
C5AR1.PTRH1.expr.files <- list.files(path = "/data2/CRC_m6A/Analysis/",pattern = "^(dt).+(C5AR1).+(PTRH1).+(RDS)$",full.names = T)
names(C5AR1.PTRH1.expr.files) <- list.files(path = "/data2/CRC_m6A/Analysis/",pattern = "^(dt).+(C5AR1).+(PTRH1).+(RDS)$",full.names = F) %>% strsplit(split=".",fixed=T) %>% sapply("[",2)

#C5AR1 in remain cohort (GSE134525)
dt.DEG.expr.C5AR1.GSE134525 <- readRDS(C5AR1.PTRH1.expr.files["GSE134525"]) %>% filter(Gene.symbol=="C5AR1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="C5AR1" & StudyID=="GSE134525",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="C5AR1" & StudyID=="GSE134525",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.C5AR1.GSE134525 <- dt.DEG.expr.C5AR1.GSE134525 %>% inner_join(x=.,y=dt.DEG.expr.C5AR1.GSE134525  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.C5AR1.GSE134525 <- dt.DEG.expr.C5AR1.GSE134525 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.C5AR1.GSE134525$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of C5AR1",title="GSE134525", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.C5AR1.GSE134525$Label1)),caption  = unique(dt.DEG.expr.C5AR1.GSE134525$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#PTRH1 in remain cohort ()
dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & !(StudyID %in% dt.top2.TumorDGE.C5AR1.PTRH1[Symbol=="PTRH1",StudyID]),]
#    Direction Symbol             Pvalue              FDR       LFC   StudyID
# 1:        Up  PTRH1 0.0000000000993405 0.00000001628186 0.7836069  GSE50421
# 2:        Up  PTRH1 0.0000018262559915 0.00001138278121 0.9089827  GSE50760
# 3:        Up  PTRH1 0.0000058101946466 0.00007343434023 0.9319651 GSE113513
# 4:        Up  PTRH1 0.0000398824511915 0.00071753731752 1.4735995  GSE50117
# 5:        Up  PTRH1 0.0000879773466835 0.00073886117019 0.2875475 GSE103512
# 6:        Up  PTRH1 0.0000022869929981 0.00157906976631 0.9259380 GSE110224
# 7:        Up  PTRH1 0.0001331035452092 0.00592989286661 1.1782485  GSE41657
# 8:        Up  PTRH1 0.0047919823865727 0.03203495409821 0.8125803 GSE156355

#PTRH1 in GSE50421
dt.DEG.expr.PTRH1.GSE50421 <- readRDS(C5AR1.PTRH1.expr.files["GSE50421"]) %>% filter(Gene.symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50421",FDR]),
                                                                                                                         Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50421",LFC]),
                                                                                                                         Condition=condition)
dt.DEG.expr.PTRH1.GSE50421 <- dt.DEG.expr.PTRH1.GSE50421 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE50421  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                          by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE50421 <- dt.DEG.expr.PTRH1.GSE50421 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE50421$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE50421", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE50421$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE50421$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE50421
#PTRH1 in GSE50760
dt.DEG.expr.PTRH1.GSE50760 <- readRDS(C5AR1.PTRH1.expr.files["GSE50760"]) %>% filter(Symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50760",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50760",LFC]),
                                                                                                                       Condition=condition2, logExpr=log2(FPM))
p.DEG.expr.PTRH1.GSE50760 <- dt.DEG.expr.PTRH1.GSE50760 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE50760$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE50760", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE50760$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE50760$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE50760

#PTRH1 in GSE113513
dt.DEG.expr.PTRH1.GSE113513 <- readRDS(C5AR1.PTRH1.expr.files["GSE113513"]) %>% filter(Gene.Symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE113513",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE113513",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.PTRH1.GSE113513 <- dt.DEG.expr.PTRH1.GSE113513 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE113513  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE113513 <- dt.DEG.expr.PTRH1.GSE113513 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE113513$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE113513", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE113513$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE113513$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE113513

#PTRH1 in GSE50117
dt.DEG.expr.PTRH1.GSE50117 <- readRDS(C5AR1.PTRH1.expr.files["GSE50117"]) %>% filter(Gene.symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50117",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE50117",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.PTRH1.GSE50117 <- dt.DEG.expr.PTRH1.GSE50117 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE50117  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE50117 <- dt.DEG.expr.PTRH1.GSE50117 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE50117$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE50117", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE50117$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE50117$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE50117

#PTRH1 in GSE103512
dt.DEG.expr.PTRH1.GSE103512 <- readRDS(C5AR1.PTRH1.expr.files["GSE103512"]) %>% filter(Gene.symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE103512",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE103512",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.PTRH1.GSE103512 <- dt.DEG.expr.PTRH1.GSE103512 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE103512  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE103512 <- dt.DEG.expr.PTRH1.GSE103512 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE103512$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE103512", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE103512$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE103512$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE103512

#PTRH1 in GSE110224
dt.DEG.expr.PTRH1.GSE110224 <- readRDS(C5AR1.PTRH1.expr.files["GSE110224"]) %>% filter(Gene.symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE110224",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE110224",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.PTRH1.GSE110224 <- dt.DEG.expr.PTRH1.GSE110224 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE110224  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE110224 <- dt.DEG.expr.PTRH1.GSE110224 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE110224$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE110224", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE110224$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE110224$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE110224

#PTRH1 in GSE41657
dt.DEG.expr.PTRH1.GSE41657 <- readRDS(C5AR1.PTRH1.expr.files["GSE41657"]) %>% filter(Gene.symbol=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE41657",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE41657",LFC]),
                                                                                                                       Condition=condition2) %>% filter(condition1 %in% c("Normal","Tumor"))
dt.DEG.expr.PTRH1.GSE41657 <- dt.DEG.expr.PTRH1.GSE41657 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE41657  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE41657 <- dt.DEG.expr.PTRH1.GSE41657 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE41657$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","NonTumor"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE41657", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE41657$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE41657$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE41657

#PTRH1 in GSE156355
dt.DEG.expr.PTRH1.GSE156355 <- readRDS(C5AR1.PTRH1.expr.files["GSE156355"]) %>% filter(GENE_SYMBOL=="PTRH1") %>%  mutate(Label2=sprintf("FDR = %.2g",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE156355",FDR]),
                                                                                                                       Label1=sprintf("LFC = %.2f",dt.TumorDEG.C5AR1.PRTH1.CRC[Symbol=="PTRH1" & StudyID=="GSE156355",LFC]),
                                                                                                                       Condition=condition)
dt.DEG.expr.PTRH1.GSE156355 <- dt.DEG.expr.PTRH1.GSE156355 %>% inner_join(x=.,y=dt.DEG.expr.PTRH1.GSE156355  %>% distinct(ProbeID,adj.P.Val) %>% dplyr::arrange(adj.P.Val) %>% slice_head(n=1),
                                                                        by=c("ProbeID","adj.P.Val"))

p.DEG.expr.PTRH1.GSE156355 <- dt.DEG.expr.PTRH1.GSE156355 %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent",scale = c("area","count","width")[3])+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.DEG.expr.PTRH1.GSE156355$Condition),2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.2)))+
  scale_color_manual(values=c("#ef756c","#839dcf"),breaks = c("Tumor","Adjacent"),guide="none")+
  # facet_wrap(~GeneID)+
  labs(x=NULL, y="log2 expression of PTRH1",title="GSE156355", subtitle = paste0("TumorDEG ",unique(dt.DEG.expr.PTRH1.GSE156355$Label1)),caption  = unique(dt.DEG.expr.PTRH1.GSE156355$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
p.DEG.expr.PTRH1.GSE156355

pageCreate(width =18.1, height =24, default.units = "cm",showGuides = T)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.C5AR1.GSE134525, x = 0.05, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05+2.7*1, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50421, x = 0.05+2.7*1, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05+2.7*2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50760, x = 0.05+2.7*2, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05+2.7*3, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE113513, x = 0.05+2.7*3, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50117, x = 0.05+2.7*4, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05+2.7*5, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE103512, x = 0.05+2.7*5, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
#row2
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05+2.7*0, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE110224, x = 0.05+2.7*0, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05+2.7*1, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE41657, x = 0.05+2.7*1, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 0.05+2.7*2, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE156355, x = 0.05+2.7*2, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)

#remianing stage plot for C5AR1 (TCGA_COAD)
#C5AR1 in TCGA_COAD
dt.stageDEG.TCGA_COAD <- readRDS("/data/COAD/Analysis//dt.DEG.TCGA_COAD.StageDEG.RDS")
dt.stageDEG.TCGA_COAD[Symbol=="C5AR1" & StudyID=="TCGA_COAD_StageDEG",]
# GeneID Symbol       LFC      pvalue       padj            StudyID
# 1: ENSG00000197405.8  C5AR1 0.3927598 0.03586189 0.1918104 TCGA_COAD_StageDEG
dt.StageDEG.expr.C5AR1.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="C5AR1" & SampleCondition=="PrimaryTumor") %>% mutate(logExpr=log2(tpm_unstranded+0.01)) %>%
  dplyr::select(Symbol=gene_name,logExpr,ajcc_stage,SampleID) %>% distinct() %>% mutate(Condition=dplyr::case_match(ajcc_stage,
                                                                                                                    c("I","IA") ~ "I",
                                                                                                                    c("II","IIA","IIB","IIC") ~ "II",
                                                                                                                    c("III", "IIIA","IIIB","IIIC") ~ "III",
                                                                                                                    c("IV","IVA","IVB") ~ "IV")) %>%
  mutate(Condition=factor(Condition,levels=c("I","II","III","IV"))) %>%
  mutate(Label1=sprintf("LFC = %.2f",dt.stageDEG.TCGA_COAD[Symbol=="C5AR1" & StudyID=="TCGA_COAD_StageDEG",LFC]),
         Label2=sprintf("Pvalue = %.2g",dt.stageDEG.TCGA_COAD[Symbol=="C5AR1" & StudyID=="TCGA_COAD_StageDEG",pvalue])) %>% filter(!is.na(Condition) & !is.na(logExpr))
table(dt.StageDEG.expr.C5AR1.TCGA_COAD$Condition)
p.StageDEG.expr.C5AR1.TCGA_COAD <- dt.StageDEG.expr.C5AR1.TCGA_COAD %>%
  ggplot(data=.,aes(x=Condition,y=logExpr))+
  geom_violin(aes(color=Condition),fill="transparent")+
  # ggbeeswarm::geom_beeswarm(aes(color=Condition),size=0.8,dodge.width = 0.8,alpha=1)+
  #geom_jitter(size=0.2,width = 0.3,alpha=0.5,aes(color=Tissue))+
  geom_boxplot(aes(color=Condition),fill="transparent",width=0.1,linewidth=0.2,outlier.shape = NA)+
  ggsignif::geom_signif(test="wilcox.test",comparisons = combn(levels(dt.StageDEG.expr.C5AR1.TCGA_COAD$Condition)[c(1,4)],2,simplify = F), test.args = list(paired=FALSE),
                        map_signif_level = function(p) sprintf("wilcox p = %.2g", p), textsize = 2, size = 0.2,step_increase = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0.05,0.05)))+
  # coord_cartesian(ylim = c(-6.5,-1))+
  scale_color_manual(values=c("#a7c9e7","#feee84","#f9b56d","#f48485"),breaks = c("I","II","III","IV"),guide="none")+
  labs(x="Tumor stage", y="log2 expression of C5AR1",title="TCGA_COAD", subtitle = paste0("ProgressionDEG ",unique(dt.StageDEG.expr.C5AR1.TCGA_COAD$Label1)), caption = unique(dt.StageDEG.expr.C5AR1.TCGA_COAD$Label2))+
  theme_minimal(base_size = 6,base_family = "Helvetica") +
  theme(
    legend.position = "top",legend.key.height = unit(7,"pt"),legend.key.width = unit(7,"pt"),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(face = "plain"),
    # axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
    plot.title = element_text(face = "plain", size = rel(1),hjust = 0.5),
    plot.subtitle = element_text(size = rel(1),hjust = 0.5),
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_blank(),
    # panel.grid.minor = element_line(linewidth = 0.6,linetype = 4)
  )
p.StageDEG.expr.C5AR1.TCGA_COAD

#remaining survival plot for PTRH1 and C5AR1
dt.SurvivalGene.C5AR1.PTRH1.CRC <- readRDS("/data2/CRC_m6A/Analysis/dt.SurvivalGene.res.C5AR1.PRTH1.RDS") %>%
  dplyr::filter(Pvalue<=0.05 & Direction=="Risk") %>% mutate(SurvivalMetricClass=factor(SurvivalMetric, levels=c("OS","PFS","PFI","DFS","DSS"), labels=c("OS", rep("PFS",4)))) %>%
  dplyr::arrange(Symbol,SurvivalMetricClass,Pvalue)
unique(dt.SurvivalGene.C5AR1.PTRH1.CRC$StudyID)
# "CPTAC2_Colon"   "TCGA_COAD"      "TCGA_READ"      "SidraLUMC2022"  "Rectal_MSK2022"
dt.top.SurvivalGene.C5AR1.PTRH1.CRC <- dt.SurvivalGene.C5AR1.PTRH1.CRC %>% group_by(Symbol,SurvivalMetricClass) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table()
#remaining plot
dt.remaining.SurvivalGene.C5AR1.PTRH1 <- rbind(dt.SurvivalGene.C5AR1.PTRH1.CRC %>% filter(Symbol=="C5AR1" & SurvivalMetricClass=="PFS" & !StudyID %in% dt.top.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & SurvivalMetricClass=="PFS" ,StudyID]) %>%
                                           dplyr::arrange(Pvalue) %>% group_by(StudyID) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table(),
                                           dt.SurvivalGene.C5AR1.PTRH1.CRC %>% filter(Symbol=="C5AR1" & SurvivalMetricClass=="OS" & !StudyID %in% dt.top.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & SurvivalMetricClass=="OS" ,StudyID]) %>%
                                             dplyr::arrange(Pvalue) %>% group_by(StudyID) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table(),
                                           dt.SurvivalGene.C5AR1.PTRH1.CRC %>% filter(Symbol=="PTRH1" & SurvivalMetricClass=="PFS" & !StudyID %in% dt.top.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & SurvivalMetricClass=="PFS" ,StudyID]) %>%
                                             dplyr::arrange(Pvalue) %>% group_by(StudyID) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table(),
                                           dt.SurvivalGene.C5AR1.PTRH1.CRC %>% filter(Symbol=="PTRH1" & SurvivalMetricClass=="OS" & !StudyID %in% dt.top.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & SurvivalMetricClass=="OS" ,StudyID]) %>%
                                             dplyr::arrange(Pvalue) %>% group_by(StudyID) %>% dplyr::arrange(Pvalue) %>% slice_head(n=1) %>% as.data.table())
#PTRH1 for PFS of TCGA_COAD
dt.SurvivalGene.PTRH1.PFS.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="PTRH1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_unstranded,SampleID,PFS=DSS,PFS.time=DSS.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="TCGA_COAD" & SurvivalMetric=="PFI" & Direction=="Risk",Pvalue])) %>%
  mutate(PFS.time=PFS.time/30)
#determine the optimal cutoff for PFS
require(survminer)
require(survival)
PFS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.PFS.TCGA_COAD, time="PFS.time",event = "PFS",variables = c("Expr"))
PFS.dt.cutpoint <- summary(PFS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
PFS.dt.cutpoint
dt.SurvivalGene.PTRH1.PFS.TCGA_COAD <- dt.SurvivalGene.PTRH1.PFS.TCGA_COAD %>% mutate(ExprGroup=if_else(Expr > PFS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
PFS.logrank.res <- survdiff(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.PFS.TCGA_COAD)
PFS.logrank.p = 1 - pchisq(PFS.logrank.res$chisq, length(PFS.logrank.res$n) - 1)
PFS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.PFS.TCGA_COAD.fit <- survfit(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.PFS.TCGA_COAD)
p.PTRH1.PFS.TCGA_COAD <- ggsurvplot(PTRH1.PFS.TCGA_COAD.fit, data = dt.SurvivalGene.PTRH1.PFS.TCGA_COAD,
                                    risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                    pval = TRUE, pval.method = T,pval.size=2,
                                    censor.size=2,size=0.5,
                                    legend.title="PTRH1",legend=c(0.8,0.4),legend.labs=c("low","high"),
                                    palette = c("#80bcb0",  "#e9786a"),
                                    break.time.by =12,
                                    ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                      theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                            legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                            legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                            strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="PFS probability",x="Survival time (months)", title="TCGA_COAD")
# p.PTRH1.PFS.TCGA_COAD$plot <- p.PTRH1.PFS.TCGA_COAD$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.PFS.TCGA_COAD$table <- p.PTRH1.PFS.TCGA_COAD$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.PFS.TCGA_COAD <- cowplot::plot_grid(p.PTRH1.PFS.TCGA_COAD$plot, p.PTRH1.PFS.TCGA_COAD$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.PFS.TCGA_COAD <- p.PTRH1.PFS.TCGA_COAD$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.PFS.TCGA_COAD
# pageCreate(width = 18, height =24, default.units = "cm",showGuides = F)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 0.05+2.7*3, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.PTRH1.PFS.TCGA_COAD, x = 0.05+2.7*3, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)

#PTRH1 for OS of CPTAC2_Colon
dt.SurvivalGene.PTRH1.OS.CPTAC2 <- readRDS(C5AR1.PTRH1.expr.files["CPTAC2_Colon"]) %>% filter(Symbol=="PTRH1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol,Expr=fpkm_uq,SampleID=Case_ID,OS=os_status,OS.time=os_months) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="CPTAC2_Colon" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=as.numeric(OS.time),OS=case_match(OS,"0:LIVING" ~ 0,"1:DECEASED" ~ 1)) %>% filter(!is.na(OS.time) & !is.na(OS))
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.OS.CPTAC2, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
dt.SurvivalGene.PTRH1.OS.CPTAC2 <- dt.SurvivalGene.PTRH1.OS.CPTAC2 %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.CPTAC2)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.OS.CPTAC2.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.CPTAC2)
p.PTRH1.OS.CPTAC2 <- ggsurvplot(PTRH1.OS.CPTAC2.fit, data = dt.SurvivalGene.PTRH1.OS.CPTAC2,
                                risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                pval = TRUE, pval.method = T,pval.size=2,
                                censor.size=2,size=0.5,
                                legend.title="PTRH1",legend=c(0.8,0.3),legend.labs=c("low","high"),
                                palette = c("#80bcb0",  "#e9786a"),
                                break.time.by =12,
                                ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                  theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                        legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                        legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                        strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="CPTAC2_Colon")
# p.PTRH1.OS.CPTAC2$plot <- p.PTRH1.OS.CPTAC2$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.OS.CPTAC2$table <- p.PTRH1.OS.CPTAC2$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.OS.CPTAC2 <- cowplot::plot_grid(p.PTRH1.OS.CPTAC2$plot, p.PTRH1.OS.CPTAC2$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.OS.CPTAC2 <- p.PTRH1.OS.CPTAC2$plot
p.PTRH1.OS.CPTAC2
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.PTRH1.OS.CPTAC2, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#PTRH1 for OS of Rectal_MSK2022
dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022 <- readRDS(C5AR1.PTRH1.expr.files["Rectal_MSK2022"]) %>% filter(Symbol=="PTRH1") %>%
  dplyr::select(Symbol,Expr=fpkm,SampleID=CaseID,OS,OS.time=OS_month) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="Rectal_MSK2022" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=as.numeric(OS.time),OS=case_match(OS,"0:LIVING" ~ 0,"1:DECEASED" ~ 1)) %>% filter(!is.na(OS.time) & !is.na(OS))
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022 <- dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022 %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.OS.Rectal_MSK2022.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022)
p.PTRH1.OS.Rectal_MSK2022 <- ggsurvplot(PTRH1.OS.Rectal_MSK2022.fit, data = dt.SurvivalGene.PTRH1.OS.Rectal_MSK2022,
                                risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                pval = TRUE, pval.method = T,pval.size=2,
                                censor.size=2,size=0.5,
                                legend.title="PTRH1",legend=c(0.8,0.3),legend.labs=c("low","high"),
                                palette = c("#80bcb0",  "#e9786a"),
                                break.time.by =12,
                                ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                  theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                        legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                        legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                        strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="Rectal_MSK2022_Colon")
# p.PTRH1.OS.Rectal_MSK2022$plot <- p.PTRH1.OS.Rectal_MSK2022$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.OS.Rectal_MSK2022$table <- p.PTRH1.OS.Rectal_MSK2022$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.OS.Rectal_MSK2022 <- cowplot::plot_grid(p.PTRH1.OS.Rectal_MSK2022$plot, p.PTRH1.OS.Rectal_MSK2022$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.OS.Rectal_MSK2022 <- p.PTRH1.OS.Rectal_MSK2022$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.OS.Rectal_MSK2022
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.PTRH1.OS.Rectal_MSK2022, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#PTRH1 for OS of SidraLUMC2022
dt.SurvivalGene.PTRH1.OS.SidraLUMC2022 <- readRDS(C5AR1.PTRH1.expr.files["SidraLUMC2022"]) %>% filter(Symbol=="PTRH1") %>%
  dplyr::select(Symbol,Expr=fpkm,SampleID=CaseID,OS,OS.time=OS.month) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="SidraLUMC2022" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=as.numeric(OS.time),OS=case_match(OS,"0:LIVING" ~ 0,"1:DECEASED" ~ 1)) %>% filter(!is.na(OS.time) & !is.na(OS))
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.OS.SidraLUMC2022, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
dt.SurvivalGene.PTRH1.OS.SidraLUMC2022 <- dt.SurvivalGene.PTRH1.OS.SidraLUMC2022 %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.SidraLUMC2022)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.OS.SidraLUMC2022.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.SidraLUMC2022)
p.PTRH1.OS.SidraLUMC2022 <- ggsurvplot(PTRH1.OS.SidraLUMC2022.fit, data = dt.SurvivalGene.PTRH1.OS.SidraLUMC2022,
                                        risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                        pval = TRUE, pval.method = T,pval.size=2,
                                        censor.size=2,size=0.5,
                                        legend.title="PTRH1",legend=c(0.8,0.8),legend.labs=c("low","high"),
                                        palette = c("#80bcb0",  "#e9786a"),
                                        break.time.by =12,
                                        ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                          theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                                legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                                legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                                strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="SidraLUMC2022")
# p.PTRH1.OS.SidraLUMC2022$plot <- p.PTRH1.OS.SidraLUMC2022$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.OS.SidraLUMC2022$table <- p.PTRH1.OS.SidraLUMC2022$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.OS.SidraLUMC2022 <- cowplot::plot_grid(p.PTRH1.OS.SidraLUMC2022$plot, p.PTRH1.OS.SidraLUMC2022$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.OS.SidraLUMC2022 <- p.PTRH1.OS.SidraLUMC2022$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.OS.SidraLUMC2022
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.PTRH1.OS.SidraLUMC2022, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#PTRH1 OS in TCGA_READ
dt.SurvivalGene.PTRH1.OS.TCGA_READ <- readRDS(C5AR1.PTRH1.expr.files["TCGA_READ"]) %>% filter(gene_name=="PTRH1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_uq_unstranded,SampleID,OS,OS.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="PTRH1" & StudyID=="TCGA_READ" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=OS.time/30)
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.PTRH1.OS.TCGA_READ, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   0.0094  3.239751
dt.SurvivalGene.PTRH1.OS.TCGA_READ <- dt.SurvivalGene.PTRH1.OS.TCGA_READ %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.TCGA_READ)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
PTRH1.OS.TCGA_READ.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.PTRH1.OS.TCGA_READ)
p.PTRH1.OS.TCGA_READ <- ggsurvplot(PTRH1.OS.TCGA_READ.fit, data = dt.SurvivalGene.PTRH1.OS.TCGA_READ,
                                    risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                    pval = TRUE, pval.method = T,pval.size=2,
                                    censor.size=2,size=0.5,
                                    legend.title="PTRH1",legend=c(0.8,0.7),legend.labs=c("low","high"),
                                    palette = c("#80bcb0",  "#e9786a"),
                                    break.time.by =12,
                                    ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                      theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                            legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                            legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                            strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="TCGA_READ")
# p.PTRH1.OS.TCGA_READ$plot <- p.PTRH1.OS.TCGA_READ$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.PTRH1.OS.TCGA_READ$table <- p.PTRH1.OS.TCGA_READ$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.PTRH1.OS.TCGA_READ <- cowplot::plot_grid(p.PTRH1.OS.TCGA_READ$plot, p.PTRH1.OS.TCGA_READ$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.PTRH1.OS.TCGA_READ <- p.PTRH1.OS.TCGA_READ$plot+coord_cartesian(xlim = c(0,96))
p.PTRH1.OS.TCGA_READ
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.PTRH1.OS.TCGA_READ, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#C5AR1 PFS in SidraLUMC2022
dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022 <- readRDS(C5AR1.PTRH1.expr.files["SidraLUMC2022"]) %>% filter(Symbol=="C5AR1") %>%
  dplyr::select(Symbol,Expr=fpkm,SampleID=CaseID,PFS,PFS.time=PFS.month) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="SidraLUMC2022" & SurvivalMetric=="PFS" & Direction=="Risk",Pvalue])) %>%
  mutate(PFS.time=as.numeric(PFS.time),PFS=case_match(PFS,"0:DiseaseFree" ~ 0,"1:Recurred" ~ 1)) %>% filter(!is.na(PFS.time) & !is.na(PFS))
#determine the optimal cutoff for PFS
require(survminer)
require(survival)
PFS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022, time="PFS.time",event = "PFS",variables = c("Expr"))
PFS.dt.cutpoint <- summary(PFS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
PFS.dt.cutpoint
dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022 <- dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022 %>% mutate(ExprGroup=if_else(Expr > PFS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
PFS.logrank.res <- survdiff(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022)
PFS.logrank.p = 1 - pchisq(PFS.logrank.res$chisq, length(PFS.logrank.res$n) - 1)
PFS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.PFS.SidraLUMC2022.fit <- survfit(Surv(PFS.time, PFS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022)
p.C5AR1.PFS.SidraLUMC2022 <- ggsurvplot(C5AR1.PFS.SidraLUMC2022.fit, data = dt.SurvivalGene.C5AR1.PFS.SidraLUMC2022,
                                       risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                       pval = TRUE, pval.method = T,pval.size=2,
                                       censor.size=2,size=0.5,
                                       legend.title="C5AR1",legend=c(0.8,0.3),legend.labs=c("low","high"),
                                       palette = c("#80bcb0",  "#e9786a"),
                                       break.time.by =12,
                                       ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                         theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                               legend.pPFSition=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                               legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                               strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="PFS probability",x="Survival time (months)", title="SidraLUMC2022")
# p.C5AR1.PFS.SidraLUMC2022$plot <- p.C5AR1.PFS.SidraLUMC2022$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.PFS.SidraLUMC2022$table <- p.C5AR1.PFS.SidraLUMC2022$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.PFS.SidraLUMC2022 <- cowplot::plot_grid(p.C5AR1.PFS.SidraLUMC2022$plot, p.C5AR1.PFS.SidraLUMC2022$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.PFS.SidraLUMC2022 <- p.C5AR1.PFS.SidraLUMC2022$plot+coord_cartesian(xlim = c(0,96))
p.C5AR1.PFS.SidraLUMC2022
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.C5AR1.PFS.SidraLUMC2022, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#C5AR1 OS in SidraLUMC2022
dt.SurvivalGene.C5AR1.OS.SidraLUMC2022 <- readRDS(C5AR1.PTRH1.expr.files["SidraLUMC2022"]) %>% filter(Symbol=="C5AR1") %>%
  dplyr::select(Symbol,Expr=fpkm,SampleID=CaseID,OS,OS.time=OS.month) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="SidraLUMC2022" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=as.numeric(OS.time),OS=case_match(OS,"0:LIVING" ~ 0,"1:DECEASED" ~ 1)) %>% filter(!is.na(OS.time) & !is.na(OS))
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.OS.SidraLUMC2022, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
dt.SurvivalGene.C5AR1.OS.SidraLUMC2022 <- dt.SurvivalGene.C5AR1.OS.SidraLUMC2022 %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.SidraLUMC2022)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.OS.SidraLUMC2022.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.SidraLUMC2022)
p.C5AR1.OS.SidraLUMC2022 <- ggsurvplot(C5AR1.OS.SidraLUMC2022.fit, data = dt.SurvivalGene.C5AR1.OS.SidraLUMC2022,
                                       risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                       pval = TRUE, pval.method = T,pval.size=2,
                                       censor.size=2,size=0.5,
                                       legend.title="C5AR1",legend=c(0.8,0.3),legend.labs=c("low","high"),
                                       palette = c("#80bcb0",  "#e9786a"),
                                       break.time.by =12,
                                       ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                         theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                               legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                               legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                               strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="SidraLUMC2022")
# p.C5AR1.OS.SidraLUMC2022$plot <- p.C5AR1.OS.SidraLUMC2022$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.OS.SidraLUMC2022$table <- p.C5AR1.OS.SidraLUMC2022$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.OS.SidraLUMC2022 <- cowplot::plot_grid(p.C5AR1.OS.SidraLUMC2022$plot, p.C5AR1.OS.SidraLUMC2022$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.OS.SidraLUMC2022 <- p.C5AR1.OS.SidraLUMC2022$plot+coord_cartesian(xlim = c(0,96))
p.C5AR1.OS.SidraLUMC2022
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.C5AR1.OS.SidraLUMC2022, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#C5AR1 OS in TCGA_COAD
dt.SurvivalGene.C5AR1.OS.TCGA_COAD <- readRDS(C5AR1.PTRH1.expr.files["TCGA_COAD"]) %>% filter(gene_name=="C5AR1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=fpkm_uq_unstranded,SampleID,OS,OS.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="TCGA_COAD" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=OS.time/30)
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.OS.TCGA_COAD, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint

dt.SurvivalGene.C5AR1.OS.TCGA_COAD <- dt.SurvivalGene.C5AR1.OS.TCGA_COAD %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.TCGA_COAD)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.OS.TCGA_COAD.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.TCGA_COAD)
p.C5AR1.OS.TCGA_COAD <- ggsurvplot(C5AR1.OS.TCGA_COAD.fit, data = dt.SurvivalGene.C5AR1.OS.TCGA_COAD,
                                   risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                   pval = TRUE, pval.method = T,pval.size=2,
                                   censor.size=2,size=0.5,
                                   legend.title="C5AR1",legend=c(0.8,0.2),legend.labs=c("low","high"),
                                   palette = c("#80bcb0",  "#e9786a"),
                                   break.time.by =12,
                                   ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                     theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                           legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                           legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                           strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="TCGA_COAD")
# p.C5AR1.OS.TCGA_COAD$plot <- p.C5AR1.OS.TCGA_COAD$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.OS.TCGA_COAD$table <- p.C5AR1.OS.TCGA_COAD$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.OS.TCGA_COAD <- cowplot::plot_grid(p.C5AR1.OS.TCGA_COAD$plot, p.C5AR1.OS.TCGA_COAD$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.OS.TCGA_COAD <- p.C5AR1.OS.TCGA_COAD$plot+coord_cartesian(xlim = c(0,96))
p.C5AR1.OS.TCGA_COAD
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.C5AR1.OS.TCGA_COAD, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)

#C5AR1 OS in TCGA_READ
dt.SurvivalGene.C5AR1.OS.TCGA_READ <- readRDS(C5AR1.PTRH1.expr.files["TCGA_READ"]) %>% filter(gene_name=="C5AR1" & SampleCondition=="PrimaryTumor") %>%
  dplyr::select(Symbol=gene_name,Expr=tpm_unstranded,SampleID,OS,OS.time) %>% distinct()  %>%
  mutate(Label=sprintf("Pvalue = %.2g",dt.SurvivalGene.C5AR1.PTRH1.CRC[Symbol=="C5AR1" & StudyID=="TCGA_READ" & SurvivalMetric=="OS" & Direction=="Risk",Pvalue])) %>%
  mutate(OS.time=OS.time/30)
#determine the optimal cutoff for OS
require(survminer)
require(survival)
OS.res.cut <- surv_cutpoint(dt.SurvivalGene.C5AR1.OS.TCGA_READ, time="OS.time",event = "OS",variables = c("Expr"))
OS.dt.cutpoint <- summary(OS.res.cut) %>% as.data.table(keep.rownames="var") %>% arrange(desc(statistic))
OS.dt.cutpoint
# var cutpoint statistic
# 1: Expr   0.0094  3.239751
dt.SurvivalGene.C5AR1.OS.TCGA_READ <- dt.SurvivalGene.C5AR1.OS.TCGA_READ %>% mutate(ExprGroup=if_else(Expr > OS.dt.cutpoint[var=="Expr",cutpoint], "high", "low") %>% factor(levels=c("low","high")))
OS.logrank.res <- survdiff(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.TCGA_READ)
OS.logrank.p = 1 - pchisq(OS.logrank.res$chisq, length(OS.logrank.res$n) - 1)
OS.logrank.p
library(survminer)
library(ggpubr)
library(ggplot2)
C5AR1.OS.TCGA_READ.fit <- survfit(Surv(OS.time, OS) ~ ExprGroup, data = dt.SurvivalGene.C5AR1.OS.TCGA_READ)
p.C5AR1.OS.TCGA_READ <- ggsurvplot(C5AR1.OS.TCGA_READ.fit, data = dt.SurvivalGene.C5AR1.OS.TCGA_READ,
                                   risk.table = F, fontsize=2, tables.height = 0.2,table.y.text=TRUE,
                                   pval = TRUE, pval.method = T,pval.size=2,
                                   censor.size=2,size=0.5,
                                   legend.title="C5AR1",legend=c(0.8,0.3),legend.labs=c("low","high"),
                                   palette = c("#80bcb0",  "#e9786a"),
                                   break.time.by =12,
                                   ggtheme = theme_pubr(base_size = 6,base_family = "Helvetica")+
                                     theme(plot.title = element_text(size=6,hjust = 0.5), plot.subtitle = element_text(size=6,hjust = 0.5),
                                           legend.position=c(0.1,0.8),legend.key.height = unit(7,"point"),legend.key.width=unit(7,"point"),
                                           legend.text = element_text(size = 6),legend.title = element_text(size=6),legend.spacing = unit(1,"point"),
                                           strip.background = element_rect(fill="white")))+
  # scale_color_manual(values = c("#8fc750",  "#f2a036"),breaks = c("low","high"))+
  labs(y="OS probability",x="Survival time (months)", title="TCGA_READ")
# p.C5AR1.OS.TCGA_READ$plot <- p.C5AR1.OS.TCGA_READ$plot+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
# p.C5AR1.OS.TCGA_READ$table <- p.C5AR1.OS.TCGA_READ$table+labs(y="")+theme(plot.title = element_blank(),axis.line=element_blank())
# p.C5AR1.OS.TCGA_READ <- cowplot::plot_grid(p.C5AR1.OS.TCGA_READ$plot, p.C5AR1.OS.TCGA_READ$table, axis = "l", align="v", ncol=1, byrow = T,rel_heights = c(0.75,0.25))
p.C5AR1.OS.TCGA_READ <- p.C5AR1.OS.TCGA_READ$plot+coord_cartesian(xlim = c(0,36))
p.C5AR1.OS.TCGA_READ
pageCreate(width = 16, height =20, default.units = "cm",showGuides = T)
plotGG(plot = p.C5AR1.OS.TCGA_READ, x = 0.05+2.7*4, y = 0.05+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)


figs6.list <- list(p.DEG.expr.PTRH1.GSE50421=p.DEG.expr.PTRH1.GSE50421,p.DEG.expr.PTRH1.GSE50760=p.DEG.expr.PTRH1.GSE50760,
                   p.DEG.expr.PTRH1.GSE113513=p.DEG.expr.PTRH1.GSE113513,p.DEG.expr.PTRH1.GSE50117=p.DEG.expr.PTRH1.GSE50117,
                   p.DEG.expr.PTRH1.GSE103512=p.DEG.expr.PTRH1.GSE103512,p.DEG.expr.PTRH1.GSE110224=p.DEG.expr.PTRH1.GSE110224,
                   p.DEG.expr.PTRH1.GSE41657=p.DEG.expr.PTRH1.GSE41657,p.DEG.expr.PTRH1.GSE156355=p.DEG.expr.PTRH1.GSE156355,
                   p.DEG.expr.C5AR1.GSE134525=p.DEG.expr.C5AR1.GSE134525,p.StageDEG.expr.C5AR1.TCGA_COAD=p.StageDEG.expr.C5AR1.TCGA_COAD,
                   p.PTRH1.PFS.TCGA_COAD=p.PTRH1.PFS.TCGA_COAD,
                   p.PTRH1.OS.CPTAC2=p.PTRH1.OS.CPTAC2,p.PTRH1.OS.Rectal_MSK2022=p.PTRH1.OS.Rectal_MSK2022,
                   p.PTRH1.OS.SidraLUMC2022=p.PTRH1.OS.SidraLUMC2022,p.PTRH1.OS.TCGA_READ=p.PTRH1.OS.TCGA_READ,
                   p.C5AR1.PFS.SidraLUMC2022=p.C5AR1.PFS.SidraLUMC2022,p.C5AR1.OS.SidraLUMC2022=p.C5AR1.OS.SidraLUMC2022,
                   p.C5AR1.OS.TCGA_COAD=p.C5AR1.OS.TCGA_COAD,p.C5AR1.OS.TCGA_READ=p.C5AR1.OS.TCGA_READ
                   )
saveRDS(figs6.list,file="FigureS6.plot.list.RDS")
save.image("sm6APeak_FigureS6_intermediate.results.RDS")

#combine all plots together
pdf("Figure6S_Activation_and_prognostic_association_of_PTRH1_and_C5AR1_in_CRC.pdf",width = 8.2677, height = 11.693)
pageCreate(width =18.1, height =24, default.units = "cm",showGuides = F)
#row1
plotText(label = "A", fontsize = 8, fontface = "bold",x = 0.05, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50421, x = 0.05, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "B", fontsize = 8, fontface = "bold",x = 0.05+2.7*1, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50760, x = 0.05+2.7*1, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "C", fontsize = 8, fontface = "bold",x = 0.05+2.7*2, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE113513, x = 0.05+2.7*2, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "D", fontsize = 8, fontface = "bold",x = 0.05+2.7*3, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE50117, x = 0.05+2.7*3, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "E", fontsize = 8, fontface = "bold",x = 0.05+2.7*4, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE103512, x = 0.05+2.7*4, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "F", fontsize = 8, fontface = "bold",x = 0.05+2.7*5, y = 0.05,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE110224, x = 0.05+2.7*5, y=0.2, default.units = "cm",width = 2.5, height = 4.4)
#row2
plotText(label = "G", fontsize = 8, fontface = "bold",x = 0.05+2.7*0, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE41657, x = 0.05+2.7*0, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "H", fontsize = 8, fontface = "bold",x = 0.05+2.7*1, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.PTRH1.GSE156355, x = 0.05+2.7*1, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "I", fontsize = 8, fontface = "bold",x = 0.05+2.7*2, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.DEG.expr.C5AR1.GSE134525, x = 0.05+2.7*2, y=0.2+4.4*1, default.units = "cm",width = 2.5, height = 4.4)
plotText(label = "J", fontsize = 8, fontface = "bold",x = 0.05+2.7*3, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.StageDEG.expr.C5AR1.TCGA_COAD, x = 0.05+2.7*3, y=0.2+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotGG(plot = p.PTRH1.PFS.TCGA_COAD, x = 0.05+2.7*3+4.4*1, y=0.2+4.4*1, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "K", fontsize = 8, fontface = "bold",x = 0.05+2.7*3+4.4*1, y = 0.05+4.4*1,just = c("top","left"), default.units = "cm",draw=T)
#row3
plotGG(plot = p.PTRH1.OS.CPTAC2, x = 0.05+4.4*0, y=0.2+4.4*2, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "L", fontsize = 8, fontface = "bold",x = 0.05+4.4*0, y = 0.05+4.4*2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.PTRH1.OS.Rectal_MSK2022, x = 0.05+4.4*1, y=0.2+4.4*2, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "M", fontsize = 8, fontface = "bold",x = 0.05+4.4*1, y = 0.05+4.4*2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.PTRH1.OS.SidraLUMC2022, x = 0.05+4.4*2, y=0.2+4.4*2, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "N", fontsize = 8, fontface = "bold",x = 0.05+4.4*2, y = 0.05+4.4*2,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.PTRH1.OS.TCGA_READ, x = 0.05+4.4*3, y=0.2+4.4*2, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "O", fontsize = 8, fontface = "bold",x = 0.05+4.4*3, y = 0.05+4.4*2,just = c("top","left"), default.units = "cm",draw=T)
#row4
plotGG(plot = p.C5AR1.PFS.SidraLUMC2022, x = 0.05+4.4*0, y=0.2+4.4*3, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "P", fontsize = 8, fontface = "bold",x = 0.05+4.4*0, y = 0.05+4.4*3,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.C5AR1.OS.SidraLUMC2022, x = 0.05+4.4*1, y=0.2+4.4*3, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "Q", fontsize = 8, fontface = "bold",x = 0.05+4.4*1, y = 0.05+4.4*3,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.C5AR1.OS.TCGA_COAD, x = 0.05+4.4*2, y=0.2+4.4*3, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "R", fontsize = 8, fontface = "bold",x = 0.05+4.4*2, y = 0.05+4.4*3,just = c("top","left"), default.units = "cm",draw=T)
plotGG(plot = p.C5AR1.OS.TCGA_READ, x = 0.05+4.4*3, y=0.2+4.4*3, default.units = "cm",width = 4.2, height = 4.2)
plotText(label = "S", fontsize = 8, fontface = "bold",x = 0.05+4.4*3, y = 0.05+4.4*3,just = c("top","left"), default.units = "cm",draw=T)

dev.off()
