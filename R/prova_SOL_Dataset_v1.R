library(readr)
library(readxl)
library(dplyr)
library(tibble)
library(DESeq2)
library(GEOquery)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(ggpubr)  
library(patchwork)
library(igraph)
library(network)



info_pat = read_excel("GSE156902.xlsx",skip = 2)
gene_annotation = read.delim("Human.GRCh38.p13.annot.tsv", row.names = 1)
counts = read.table("GSE156902_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1)


dim(counts)# 39376   802 

# 802 samples
# 366 cancer: 126 Brain Metastases (15 breast + 1 colorectal + 1 esophageal + 15 melanoma + 85 NSCLC + 1 pancreatic + 7 renal + 1 unkown)
#  + 240 Glioblastoma
# 151 follow up in Glioblastoma
# 439 non cancer: 355 Asym Controls + 86 multiple scleroris

# stats 
##### 
table(info_pat$`Classification group`)
table(info_pat$`Classification group`,info_pat$`Follow-up`)

table(info_pat$`Patient group`)
table(info_pat$`Patient group`,info_pat$`Classification group`)
table(info_pat$`Patient group`,info_pat$`Follow-up`)


freq_table <- table(info_pat[info_pat$`Patient group`=="GBM",]$"Baseline sample ID")
freq_of_freq <- table(freq_table)
sum(freq_of_freq)
info_pat$Age <- as.numeric(info_pat$Age)
hist(info_pat$Age)

age <- as.numeric(as.character(info_pat$Age))
######
GBM=info_pat[info_pat$`Patient group`%in%c('GBM',"AC"),] 
nrow(GBM) # 593
353+89+151 # 353 AC + 89 Reesection + 151 follow up

colnames(GBM)
GBM=GBM[c("Patient group","Age","Gender","Baseline sample ID","Sample name")]

dim(GBM) # 593 x 6
length(unique(GBM$`Baseline sample ID`)) # 55

GBM[duplicated(GBM$`Sample name`),]
# GBM[duplicated(GBM$`Sample name`, fromLast = TRUE),]
# GBM=GBM[!duplicated(GBM$`Sample name`) & ! duplicated(GBM$`Sample name`, fromLast = TRUE),]
GBM=GBM[!duplicated(GBM$`Sample name`),]
dim(GBM) # 592 , removed 1 or 2 duplicates

mentioned_in_baseline <- unique(GBM$`Sample name`[GBM$`Sample name` %in% GBM$`Baseline sample ID`[!is.na(GBM$`Baseline sample ID`)]])

GBM$stage=ifelse(GBM$`Baseline sample ID`=="NA" & GBM$`Sample name` %in% mentioned_in_baseline,0,-1)
GBM$stage <- as.integer(sub(".*[Tt](\\d+)$", "\\1", GBM$`Sample name`))
table(GBM$stage)
table(GBM$stage>0)

table(is.na(GBM$stage))

GBM$stage=ifelse(is.na(GBM$stage),-1,GBM$stage)

table(GBM$stage,GBM$`Patient group`)

dim(GBM[GBM$stage==-1 & GBM$`Patient group`=="GBM",]) # 47 6
GBM$stage=ifelse(GBM$stage==-1 & GBM$`Patient group`=="GBM",0,GBM$stage)
table(GBM$stage)
table(GBM$stage>0) # 150 true= 150 followups, nel paper 151!

# GBM$stage=as.character(GBM$stage+1) # control =0 t0=1, t1=2,...
table(GBM$stage)
# Ac     0     1  2   3   4   5   6   7   8   9  10 
# 353  89 48  34  20  21   8   7   5   4   1   2

############
# DATA PREPROCESSING
gse <- getGEO("GSE156902", GSEMatrix = T)

meta1 <- pData(gse[[1]])
meta2 <- pData(gse[[2]])

clean_titles <- function(meta) {
  meta$title_clean <- sub("^Blood_Platelets_(.*?) \\[.*$", "\\1", meta$title)
  meta$title_clean <- sub("-TR\\d+$", "", meta$title_clean)
  # Prefix “VU” if title starts with a digit and has “VUMC”
  idx_vumc <- grepl("^\\d", meta$title_clean) & grepl("VUMC", meta$title, ignore.case = TRUE)
  meta$title_clean[idx_vumc] <- paste0("VU", meta$title_clean[idx_vumc])
  # Helper for replacing title_clean based on keyword and bracket contents
  replace_from_bracket <- function(keyword, pattern = ".*\\[([^]]+)\\].*", subpattern = NULL) {
    idx <- grepl(keyword, meta$title, ignore.case = TRUE)
    val <- sub(pattern, "\\1", meta$title[idx])
    if (!is.null(subpattern)) val <- sub(subpattern, "\\1", val, ignore.case = TRUE)
    val <- sub("-TR\\d+$", "", val)
    meta$title_clean[idx] <<- val }
  # Apply to different site patterns
  replace_from_bracket("Vienna")
  replace_from_bracket("VUMC-HC")
  replace_from_bracket("HD-NKI")
  replace_from_bracket("DOT-HD", subpattern = "^VUMC-(.*?)-DOT-HD.*$")
  
  meta
}

meta1 <- clean_titles(meta1)
meta2 <- clean_titles(meta2)

meta1 <- meta1[, intersect(colnames(meta1), colnames(meta2))]
meta2 <- meta2[, intersect(colnames(meta1), colnames(meta2))]

m <- rbind(meta1, meta2)
colnames(m)
m=m[c("title","geo_accession","title_clean")]
m["GSM4746640",]$title_clean="VU609 T0"



###########################
meta <- GBM %>% left_join(m %>% select(title_clean, geo_accession),by = c("Sample name" = "title_clean"))
table(is.na(meta$geo_accession))
meta  =  meta %>%column_to_rownames(var = "geo_accession")
table(meta$stage)
table(meta$stage,meta$`Patient group`)
# meta=meta[as.numeric(meta$stage)<6,]
# meta=meta[meta$stage<5,]
nrow(meta) # 592

#### ANALISI runnare questo, spiegazione dopo 
meta["GSM4746493",]$`Sample name`="Vumc-GBM-477t0" # prima era "Vumc-GBM-477t0.2" ,corretto per "Vumc-GBM-follow-up-VU477t1"
meta["GSM4746495",]$`Sample name`="Vumc-GBM-491t0"  # prima era Vumc-GBM-491t0-2 , corretto per Vumc-GBM-follow-up-VU491t1
meta["GSM4746484",]$`Sample name`="Vumc-GBM-493t0" # prima era "Vumc-GBM-493t0-2". corretto per Vumc-GBM-follow-up-VU493t1, Vumc-GBM-follow-up-VU493t1
meta["GSM4746637",]$`Baseline sample ID`="Vumc-GBM-493t0" # aggiustato anche per VU493 T4, era l'unico che corrispondeva con l'originale t0
meta["GSM4746506",]$`Sample name`="Vumc-GBM-follow-up-VU501t0" # prima era Vumc-GBM-follow-up-VU501t0-2 , corretto per t1,t3,t4
meta["GSM4746638",]$`Baseline sample ID`="VU635 T0" # prima era "VU653 T0" hanno invertito 635 con 653 

################
pazienti_T=unique(meta$`Baseline sample ID`) # pazienti con followup

stage_4=meta[meta$stage==4,]
stage_3=meta[meta$stage==3,]
stage_2=meta[meta$stage==2,]
stage_1=meta[meta$stage==1,]
stage_0=meta[meta$stage==0,]

sum(stage_0$`Sample name`%in% pazienti_T) # 48 ! IMPOSSIBILE
pazienti_og_3=stage_3$`Baseline sample ID`

sum(stage_3$`Baseline sample ID` %in% pazienti_og) # 20
sum(stage_2$`Baseline sample ID` %in% pazienti_og) # 14
sum(stage_1$`Baseline sample ID` %in% pazienti_og) # 19
sum(stage_0$`Sample name` %in% pazienti_og) # 18


stage_3$sta_in_t0=ifelse(stage_3$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_3$sta_in_t0)
stage_3[stage_3$sta_in_t0=="NO",]

stage_4$sta_in_t0=ifelse(stage_4$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_4$sta_in_t0)


stage_1$sta_in_t0=ifelse(stage_1$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_1$sta_in_t0)
# 5 NO.... INDAGHIAMO
## 	 Vumc-GBM-429t1 questo non ha nessun sample a tempo 0
meta["GSM4746493",]$`Sample name`="Vumc-GBM-477t0" # prima era "Vumc-GBM-477t0.2" ,corretto per "Vumc-GBM-follow-up-VU477t1"
meta["GSM4746495",]$`Sample name`="Vumc-GBM-491t0"  # prima era Vumc-GBM-491t0-2 , corretto per Vumc-GBM-follow-up-VU491t1
meta["GSM4746484",]$`Sample name`="Vumc-GBM-493t0" # prima era "Vumc-GBM-493t0-2". corretto per Vumc-GBM-follow-up-VU493t1, Vumc-GBM-follow-up-VU493t1
meta["GSM4746637",]$`Baseline sample ID`="Vumc-GBM-493t0" # aggiustato anche per VU493 T4, era l'unico che corrispondeva con l'originale t0
meta["GSM4746506",]$`Sample name`="Vumc-GBM-follow-up-VU501t0" # prima era Vumc-GBM-follow-up-VU501t0-2 , corretto per t1,t3,t4


stage_2$sta_in_t0=ifelse(stage_2$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_2$sta_in_t0)
# 1 NO ... INDAGHIAMO!
meta["GSM4746638",]$`Baseline sample ID`="VU635 T0" # prima er "VU653 T0" hanno invertito 635 con 653 


stage_3$sta_in_t0=ifelse(stage_3$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_3$sta_in_t0)
# 2 NO, ma dovrebbero essere corretti ora!

stage_4$sta_in_t0=ifelse(stage_4$`Baseline sample ID` %in% stage_0$`Sample name`,"YES","NO")
table(stage_4$sta_in_t0)
# 2 NO.... INDAGHIAMO!
# Vumc-GBM-follow-up-VU429t4 non ha baseline a tempo 0









#################

valid_samples <- intersect(rownames(meta), colnames(counts)) # 564

counts <- counts[, valid_samples]
meta = meta[valid_samples,]
colnames(meta)[colnames(meta) == "Patient group"] <- "Group"
colnames(meta)[colnames(meta) == "stage"] <- "Stage"


##########
meta$Gender=ifelse(meta$Gender=='NA',NA,meta$Gender)
meta$Age=as.numeric(meta$Age)

table(meta$Stage)
table(meta$Gender)
table(meta$Gender,meta$Stage)
table(meta$Group,meta$Gender)



summary(meta$Age)
mean_age <- mean(meta$Age, na.rm = TRUE)
sd_age   <- sd(meta$Age, na.rm = TRUE)
median_age=median(meta$Age,na.rm=TRUE)
p1=ggplot(meta, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.3, fill = "black",color='black') +
  geom_density(color = "red", size = 0.5) +
  geom_vline(aes(xintercept = mean_age),color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Age Distribution") +
  annotate("text",x = mean_age,y = 0.021,
           label = paste0("Mean = ", round(mean_age, 2),"\n SD = ", round(sd_age, 2), "\nMedian = ",round(median_age,2)),
           vjust = +0.5,hjust = +2.8)



p2 = ggplot(meta, aes(x = Stage, y = Age, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age by Patient Group")


by(meta$Age, meta$Gender, shapiro.test) # Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta[!is.na(meta$Gender), ], aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Gender)


p3=ggplot(meta[!is.na(meta$Gender), ], aes(x = Gender, y = Age, fill = Gender)) +
  geom_boxplot(na.rm = TRUE) +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age Distribution by Gender") +
  stat_compare_means(method = "wilcox.test",label = "p.format",
                     label.y = 75,size = 5) +
  theme(legend.position = "none")


library(ggridges)
ggplot(meta, aes(x = Age, y = Stage, fill = Stage)) +
  geom_density_ridges(alpha = 1.7) +
  theme_ridges() +
  labs(title = "Age Density by Stage")


ggplot(meta[!is.na(meta$Gender), ], aes(x = interaction(Group, Gender), y = Age, fill = Gender)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Age Distribution by Group and Gender", x = "Group_Gender")




(p1 | p2) / (p3 | p4) + plot_annotation(title = "Patient Cohort Overview")




by(meta$Age, meta$Gender, shapiro.test)
# Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta, aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Gender)




comparisons <- list(c("GBM.F", "GBM.M"),c("AC.F", "AC.M"),c("AC.M","GBM.M"),c("AC.F","GBM.F"))

p4=ggplot(meta[!is.na(meta$Gender), ],
          aes(x = interaction(Group, Gender), y = Age, fill = Gender)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink", "lightblue"))+
  coord_flip() +
  theme_minimal() +
  labs(title = "Age Distribution by Group and Gender",x = "GROUP x GENDER") +
  stat_compare_means(comparisons = comparisons,method = "wilcox.test",
                     label = "p.format",size = 4,label.y = c(85, 90,72,72)) +
  theme(legend.position = "none")


(p1 / p2 / (p3|p4))  +
  plot_annotation(title = "Patient's Characteristics",
                  theme = theme(plot.title = element_text(size = 16, face = "bold",family="helvetica")))
######################

colnames(meta)

meta_AC_vs_t0=meta[meta$Stage%in%c(-1,0),]
table(meta_AC_vs_t0$Group,meta_AC_vs_t0$Stage) # 353 vs 89

valid_samples_2 <- intersect(rownames(meta_AC_vs_t0), colnames(counts))
counts_AC_vs_t0=counts[, valid_samples_2]
meta_AC_vs_t0 = meta_AC_vs_t0[valid_samples_2,]

#dds <- DESeqDataSetFromMatrix(countData = counts,colData=meta[, c("Group"), drop = FALSE],design= ~ Group)
dds <- DESeqDataSetFromMatrix(countData = counts_AC_vs_t0,colData=meta_AC_vs_t0[, c("Group"), drop = FALSE],design= ~ Group)
summary(rowMeans(counts(dds)))

keep <- rowSums(counts(dds) >= 20) >= 0.8*ncol(counts_AC_vs_t0)
sum(keep) #5348, 6501

sum(keep)/length(keep) # 0.14
dds <- dds[keep,]

dds <- DESeq(dds)

res_contr_vs_gbm=results(dds,contrast = c("Group","AC","GBM"))
res_contr_vs_gbm=res_contr_vs_gbm[order(res_contr_vs_gbm$padj),]
deg_contr_vs_gbm=subset(res_contr_vs_gbm,padj<0.05 & abs(log2FoldChange)>0.33)
dim(deg_contr_vs_gbm) #2602, 2381


norm_counts <- counts(dds, normalized = TRUE)
expr_log2 <- log2(norm_counts + 1)
expr_log2_filtr = expr_log2[rownames(deg_contr_vs_gbm),]

# PRENDIAMO SOLO QUELLI CHE SI TROVANO DAL T0 AL T3
#########
pazienti_T_2=unique(meta$`Baseline sample ID`) # pazienti con followup
length(pazienti_T_2) # 53,incluso 1 na, quindi 52 pazienti con followup, nel paper 55

stage_3=meta[meta$Stage==3,]
stage_2=meta[meta$Stage==2,]
stage_1=meta[meta$Stage==1,]
stage_0=meta[meta$Stage==0,]

meta$t0_t3=ifelse((meta$`Sample name` %in% stage_1$`Baseline sample ID`) &
                    (meta$`Sample name` %in% stage_2$`Baseline sample ID`) &
                    (meta$`Sample name` %in% stage_3$`Baseline sample ID`),1,0)
table(meta$t0_t3) # 13


meta$t1_t3=ifelse((meta$Stage==1) & ((meta$`Baseline sample ID` %in% stage_2$`Baseline sample ID`)
                                     & (meta$`Baseline sample ID` %in% stage_3$`Baseline sample ID`)
                                     & (meta$`Baseline sample ID` %in% stage_0$`Sample name`)),1,0)
table(meta$t1_t3)

meta$t2_t3=as.numeric(meta$Stage==2)
table(meta$t2_t3)
meta$t2_t3=ifelse((meta$Stage==2) & ((meta$`Baseline sample ID` %in% stage_3$`Baseline sample ID`)
                                     &(meta$`Baseline sample ID` %in% stage_1$`Baseline sample ID`)
                                     & (meta$`Baseline sample ID` %in% stage_0$`Sample name`)),1,0)
table(meta$t2_t3)

meta$t3_t3=ifelse(meta$Stage==3 & ((meta$`Baseline sample ID` %in% stage_2$`Baseline sample ID`)
                                   &(meta$`Baseline sample ID` %in% stage_1$`Baseline sample ID`) 
                                   & (meta$`Baseline sample ID` %in% stage_0$`Sample name`)),1,0)
table(meta$t3_t3)

meta$in_all_stages=ifelse(meta$t0_t3 + meta$t1_t3 + meta$t2_t3 + meta$t3_t3>0,1,0)
table(meta$in_all_stages) # 13*4=52
meta_T=meta[meta$in_all_stages==1,]
table(meta_T$Stage) 

#######

table(meta_T$Gender)
meta_T$Age=as.numeric(meta_T$Age)

meta_T$Age
table(meta_T$Stage)
table(meta_T$Gender)
table(meta_T$Gender,meta_T$Stage)
table(meta_T$Group,meta_T$Gender)

meta_T[meta_T$Stage=='0' & meta_T$Gender=='F',]$Age
meta_T[meta_T$Stage=='0' & meta_T$Gender=='M',]$Age



summary(meta_T[meta_T$Stage==0,]$Age)
mean_age <- mean(meta_T$Age, na.rm = TRUE)
sd_age   <- sd(meta_T$Age, na.rm = TRUE)
median_age=median(meta_T$Age,na.rm=TRUE)
p1=ggplot(meta_T[meta_T$Stage==0,], aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 10, alpha = 0.3, fill = "black",color='black') +
  geom_density(color = "red", size = 0.5) +
  geom_vline(aes(xintercept = mean_age),color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Age Distribution") +
  annotate("text",x = mean_age,y = 0.021,
           label = paste0("Mean = ", round(mean_age, 2),"\n SD = ", round(sd_age, 2), "\nMedian = ",round(median_age,2)),
           vjust = +0.5,hjust = +2.8)



p2 = ggplot(meta_T[meta_T$Stage==0,], aes(x = Stage, y = Age, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age by Patient Group")


by(meta_T$Age, meta_T$Gender, shapiro.test) # Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta_T[!is.na(meta$Gender), ], aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Gender)


p3=ggplot(meta_T[meta_T$Stage==0,], aes(x = Gender, y = Age, fill = Gender)) +
  geom_boxplot(na.rm = TRUE) +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age Distribution by Gender") +
  stat_compare_means(method = "wilcox.test",label = "p.format",
                     label.y = 75,size = 5) +
  theme(legend.position = "none")


library(ggridges)
ggplot(meta, aes(x = Age, y = Stage, fill = Stage)) +
  geom_density_ridges(alpha = 1.7) +
  theme_ridges() +
  labs(title = "Age Density by Stage")


ggplot(meta[!is.na(meta$Gender), ], aes(x = interaction(Group, Gender), y = Age, fill = Gender)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Age Distribution by Group and Gender", x = "Group_Gender")




(p1 | p2) / (p3 | p4) + plot_annotation(title = "Patient Cohort Overview")




by(meta$Age, meta$Gender, shapiro.test)
# Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta, aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Gender)


table(meta_T[meta_T$Stage==0,]$Age)

comparisons <- list(c("GBM.F", "GBM.M"))

p4=ggplot(meta_T[meta_T$Stage==0,],
          aes(x = interaction(Group, Gender), y = Age, fill = Gender)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink", "lightblue"))+
  coord_flip() +
  theme_minimal() +
  labs(title = "Age Distribution by Group and Gender",x = "GROUP x GENDER") +
  stat_compare_means(comparisons = comparisons,method = "wilcox.test",
                     label = "p.format",size = 4,label.y = c(85, 90,72,72)) +
  theme(legend.position = "none")


(p1 / p2 / (p3|p4))  +
  plot_annotation(title = "Patient's Characteristics",
                  theme = theme(plot.title = element_text(size = 16, face = "bold",family="helvetica")))




############
deg_genes <- rownames(deg_contr_vs_gbm)
samples_t0_t3 <- rownames(meta_T)
counts_t0_t3 <- counts[, samples_t0_t3]
dds_all <- DESeqDataSetFromMatrix(countData = counts_t0_t3,colData = meta[samples_t0_t3, , drop = FALSE],design = ~ 1)
dds_all <- estimateSizeFactors(dds_all)
norm_counts_all <- counts(dds_all, normalized = TRUE)
expr_log2_all <- log2(norm_counts_all + 1)
expr_log2_filtr <- expr_log2_all[deg_genes, , drop = FALSE]


#########

stage_samples <- split(rownames(meta_T), meta_T$Stage)

X_list <- list()
for (stage_name in names(stage_samples)) {
  stage_cells <- stage_samples[[stage_name]]
  mat <- expr_log2_filtr[, stage_cells, drop = FALSE]
  mat_centered <- scale(t(mat), center = TRUE, scale = FALSE)
  X_list[[stage_name]] <- mat_centered
  S_hat <- cov(mat_centered)
  X_path <- file.path("dataset_sol", paste0("X_", stage_name, ".csv"))
  write.csv(mat_centered, file = X_path, row.names = TRUE)
}
dim(X_list$'1')
# guardare fino al 3
###
gene_annotation$Symbol[match(rownames(expr_log2_filtr), rownames(gene_annotation))]
gene_names <- gene_annotation$Symbol[match(rownames(expr_log2_filtr), rownames(gene_annotation))]
gene_names <- ifelse(is.na(gene_names), rownames(expr_log2_filtr), gene_names)
write.csv(data.frame(genes =gene_names ),file=file.path("dataset_sol", "genes_keep.csv") , row.names = FALSE)
