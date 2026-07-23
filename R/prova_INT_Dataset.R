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

## Dataset
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE183635

# https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00370-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1535610822003701%3Fshowall%3Dtrue

info_pat = read_excel("mmc3.xlsx", skip = 2)
counts=read.table("GSE183635_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, row.names = 1, check.names = FALSE)
gene_annotation = read.delim("Human.GRCh38.p13.annot.tsv", row.names = 1)

#
# 39376 genes
# 390 Asyptomatic controls
# 522 NSCLC (16,17,64,400 + 25 na)
# 144 Ovarian Cancer  (36,18,49,35 + 6 na)
# 132 Glioma (0 stages)
# 126 Pancreatic (3,55,31,25 + 1 na)
# 101 HeadNeck (4,7,30,60)
# 93 Breast (12,18,4,57 + 2 na)
# 85 Cholangiocarcinoma (6,22,9,45 + 3 na)
# 85 Colorectal (1,3,5,58 + 18 na)

##

# info_pat  =  info_pat %>%column_to_rownames(var = "Sample ID") NON RUNNARE QUESTO!
rownames(info_pat)
colnames(counts)

################ 
# stats
table(info_pat$Sex)
sort(table(info_pat$Group),decreasing = T)
table(info_pat$Stage,info_pat$Group)
table(info_pat$Stage)
table(info_pat$Stage,info_pat$Sex)

table(info_pat$Age)
age <- as.numeric(as.character(info_pat$Age))

hist(age, breaks=20, main='AGE DISTRIBUTION', col='pink', xlab='', ylim = c(0,400),
     cex.lab = 1.2,cex.main = 1.4, font.main = 2,family = "helvetica",xaxt='n')

axis(1, at=seq(10, 90, by=5))
abline(v=mean(age,na.rm=T), col='red', lwd=3,lty=2)

mean_age <- mean(age,na.rm=T)
sd_age <- sd(age,na.rm=T)

legend("topleft",legend = paste("Mean =", round(mean_age, 2), "\nSD =", round(sd_age, 2)), 
       col = NA,lwd = 3,bty = "n")



################ 
# get metadata
gse <- getGEO("GSE183635", GSEMatrix = T)

meta1 <- pData(gse[[1]])
head(meta1,1)
meta2= pData(gse[[2]])

colnames(meta1)
common_cols <- intersect(colnames(meta1), colnames(meta2))
meta1 <- meta1[, common_cols]
meta2 <- meta2[, common_cols]
m <- rbind(meta1, meta2) # sarebbe il nuovo info_pat

common_rows=intersect(rownames(m),colnames(counts))
nrow(m) # 2351
ncol(counts) # 2334
length(common_rows) # 2334

m=m[common_rows,]
head(m)



########################

table(info_pat$Group)
144+390

info_pat <- info_pat %>%left_join(m %>% select(title, geo_accession),by = c("Sample ID" = "title"))

#### Ovarian
meta = info_pat[info_pat$Group%in%c('Ovarian cancer',"Asymptomatic controls"),]
dim(meta) # 534 12
table(meta$Group) # 534=390 AC+144
table(meta$Stage,meta$Group)
meta= meta[!(meta$Stage=="n.a." & meta$Group=="Ovarian cancer"),]
dim(meta) # 528 12
meta = meta[!(meta$Stage=="n.a." & meta$Sex!="F"),]
dim(meta) # 355
table(meta$Group)

########
meta$Stage=ifelse(meta$Stage=="n.a.","0"
                  ,ifelse(meta$Stage=="I","1"
                              ,ifelse(meta$Stage=="II","2"
                                      ,ifelse(meta$Stage=="III","3","4"))))


colnames(meta)
meta = meta[c("Group","Stage","Age","Sex","geo_accession")]
dim(meta)
length(meta$geo_accession)
length(unique(meta$geo_accession))

dup_accessions <- meta$geo_accession[duplicated(meta$geo_accession) |
                                       duplicated(meta$geo_accession, fromLast = TRUE)]

cat("Number of duplicated geo_accession entries:",length(unique(dup_accessions)), "\n")
dup_accessions
meta<- meta[!is.na(meta$geo_accession), ]
dim(meta)
meta  =  meta %>%column_to_rownames(var = "geo_accession")
table(meta$Group)
### PLOT OVARIAN
#######

dim(meta)
table(meta$Stage)
table(meta$Group)

summary(meta[meta$Stage=="0",]$Age)
summary(meta[meta$Stage=="1",]$Age)
summary(meta[meta$Stage=="2",]$Age)
summary(meta[meta$Stage=="3",]$Age)
summary(meta[meta$Stage=="4",]$Age)

meta$Age=ifelse(meta$Age=="n.a.",NA,as.integer(meta$Age))
summary(meta$Age)
mean_age <- mean(meta$Age, na.rm = TRUE)
sd_age   <- sd(meta$Age, na.rm = TRUE)

sum(is.na(meta$Age))


p1=ggplot(meta, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.1, fill = "black",color='black') +
  geom_density(color = "red", size = 0.5) +
  geom_vline(aes(xintercept = mean_age),color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Age Distribution") +
  annotate("text",x = mean_age,y = 0.03,
           label = paste0("Mean = ", round(mean_age, 2),"\nSD = ", round(sd_age, 2)),
           vjust = +0.5,hjust = +2.8)

p2 = ggplot(meta, aes(x = Stage, y = Age, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age by Patient Group")


(p1 / p2)  +
  plot_annotation(title = "Clinical and Demographic Characteristics of the Cohort",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))


by(meta$Age, meta$Group, shapiro.test)
# Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta, aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Group)




p3= ggplot(meta[!is.na(meta$Group), ], aes(x = Group, y = Age, fill = Group)) +
  geom_boxplot(na.rm = TRUE) +
  #scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age Distribution by Gender") +
  stat_compare_means(method = "wilcox.test",label = "p.format",
                     label.y = 20,size = 5) +
  theme(legend.position = "none")

(p1 / (p3|p2))  +
  plot_annotation(title = "Patient's Characteristics",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

######

########
valid_samples <- intersect(rownames(meta), colnames(counts))
counts <- counts[, valid_samples]
meta <- meta[valid_samples, ]
table(meta$Group)

dds <- DESeqDataSetFromMatrix(countData = counts,colData   = meta[, c("Group"), drop = FALSE],design= ~ Group)
summary(rowMeans(counts(dds)))
keep <- rowSums(counts(dds) >= 20) >= 0.8*ncol(counts) # 
sum(keep) #4797
sum(keep)/length(keep)
dds <- dds[keep,]

########
dds <- DESeq(dds)

res_contr_vs_ov=results(dds,contrast = c("Group","Asymptomatic controls","Ovarian cancer"))
res_contr_vs_ov=res_contr_vs_ov[order(res_contr_vs_ov$padj),]
deg_contr_vs_ov=subset(res_contr_vs_ov,padj<0.05 & abs(log2FoldChange)>0.33)
dim(deg_contr_vs_ov) # 1228

########
# dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
expr_log2 <- log2(norm_counts + 1)
expr_log2_filtr = expr_log2[rownames(deg_contr_vs_ov),]

############
table(meta$Stage)
stage_samples <- split(rownames(meta), meta$Stage)

sigma_hat_list <- list()
X_list <- list()
for (stage_name in names(stage_samples)) {
  stage_cells <- stage_samples[[stage_name]]
  mat <- expr_log2_filtr[, stage_cells, drop = FALSE]
  mat_centered <- scale(t(mat), center = TRUE, scale = FALSE)
  X_list[[stage_name]] <- mat_centered
  S_hat <- cov(mat_centered)
  sigma_hat_list[[stage_name]] <- S_hat
  # Build file paths inside the folder
  sigma_path <- file.path("dataset_int_veld", paste0("sigma_hat_", stage_name, ".csv"))
  X_path <- file.path("dataset_int_veld", paste0("X_", stage_name, ".csv"))
  
  # Write files
  write.csv(S_hat, file = sigma_path, row.names = TRUE)
  write.csv(mat_centered, file = X_path, row.names = TRUE)
}


gene_names <- gene_annotation$Symbol[match(rownames(expr_log2_filtr), rownames(gene_annotation))]
gene_names <- ifelse(is.na(gene_names), rownames(expr_log2_filtr), gene_names)
write.csv(data.frame(genes =gene_names ),file=file.path("dataset_int_veld", "genes_keep.csv") , row.names = FALSE)

####################################################################################








