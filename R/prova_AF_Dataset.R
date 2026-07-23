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




supp_tab_AF  =  read_excel("41598_2023_35818_MOESM2_ESM.xlsx",skip = 2)
gene_annotation = read.delim("Human.GRCh38.p13.annot.tsv", row.names = 1)
counts_AF = read.table("GSE207586_NSCLC_CountMatrix_2ndVersion.txt", header = TRUE, row.names = 1)


# data preprocessing
###########
new_colnames  =  gsub("\\.", "-", colnames(counts_AF)) # Replace all points with hyphens
new_colnames  =  sub("-$", "", new_colnames) # Remove trailing hyphen (from final point)
colnames(counts_AF)= new_colnames

supp_tab_AF  =  supp_tab_AF %>%column_to_rownames(var = "Sample ID")
length(intersect(rownames(supp_tab_AF),colnames(counts_AF))) # dovrebbe essere 766 (ncol di counts)

### stats
##############
colnames(supp_tab_AF)
meta_AF=supp_tab_AF[,c("Patient Group","Gender","Age","Smoking","Cancer Stage")]
head(meta_AF,5)
colnames(meta_AF)=c("Group","Gender","Age","Smoking","Stage")
meta_AF$Group=as.factor(meta_AF$Group)
meta_AF$Gender=as.factor(meta_AF$Gender)
meta_AF$Smoking=as.factor(meta_AF$Smoking)
meta_AF$Stage=as.factor(meta_AF$Stage)

table(meta_AF$Group) 
table(meta_AF$Gender)
table(meta_AF$Stage,meta_AF$Gender)
table(meta_AF$Stage,meta_AF$Group)
table(meta_AF$Smoking,meta_AF$Gender,meta_AF$Stage)
table(meta_AF$Stage,meta_AF$Smoking)


#############

hist(meta_AF$Age,col='white',border = 'black',breaks = 10)
abline(v=mean(meta_AF$Age),col='red')

summary(meta_AF$Age)
mean_age <- mean(meta_AF$Age, na.rm = TRUE)
sd_age   <- sd(meta_AF$Age, na.rm = TRUE)


p1=ggplot(meta_AF, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), bins = 15, alpha = 0.1, fill = "black",color='black') +
  geom_density(color = "red", size = 0.5) +
  geom_vline(aes(xintercept = mean_age),color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Age Distribution") +
  annotate("text",x = mean_age,y = 0.03,
           label = paste0("Mean = ", round(mean_age, 2),"\nSD = ", round(sd_age, 2)),
            vjust = +0.5,hjust = +2.8)

p2 = ggplot(meta_AF, aes(x = Stage, y = Age, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age by Patient Group")


by(meta_AF$Age, meta_AF$Gender, shapiro.test)
# Age is not normally distributed in either gender group (p-value<0.05)
ggplot(meta_AF, aes(sample = Age)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Gender)


p3=ggplot(meta_AF[!is.na(meta_AF$Gender), ], aes(x = Gender, y = Age, fill = Gender)) +
  geom_boxplot(na.rm = TRUE) +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_minimal() +
  coord_flip() +
  scale_y_continuous(breaks = seq(10, 90, 10)) +
  labs(title = "Age Distribution by Gender") +
  stat_compare_means(method = "wilcox.test",label = "p.format",
                     label.y = 85,size = 5) +
  theme(legend.position = "none")

comparisons <- list(c("NSCLC.F", "NSCLC.M"),c("Control.F", "Control.M"),c("Control.M","NSCLC.M"),c("Control.F","NSCLC.F"))

p4=ggplot(meta_AF[!is.na(meta_AF$Gender), ],
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


### filtering
#####################
valid_samples_AF <- intersect(rownames(meta_AF), colnames(counts_AF))
counts_AF <- counts_AF[, valid_samples_AF]
meta_AF <- meta_AF[valid_samples_AF, ]

dds_AF <- DESeqDataSetFromMatrix(countData = counts_AF,colData   = meta_AF[, c("Group"), drop = FALSE],design= ~ Group)
summary(rowMeans(counts(dds_AF)))
keep_AF <- rowSums(counts(dds_AF) >= 20) >= 0.80*ncol(counts_AF) #

sum(keep_AF) # 1445
sum(keep_AF)/length(keep_AF) # 0.34
dds_AF <- dds_AF[keep_AF,]

dds_AF <- DESeq(dds_AF)

res_contr_vs_nsclc=results(dds_AF,contrast = c("Group","Control","NSCLC"))
res_contr_vs_nsclc=res_contr_vs_nsclc[order(res_contr_vs_nsclc$padj),]
deg_contr_vs_nsclc=subset(res_contr_vs_nsclc,padj<0.05 & abs(log2FoldChange) > 0.33)
dim(deg_contr_vs_nsclc) #283

norm_counts <- counts(dds_AF, normalized = TRUE)

expr_log2 <- log2(norm_counts + 1) 

expr_log2_filtr = expr_log2[rownames(deg_contr_vs_nsclc),] # keep only DE genes
dim(expr_log2_filtr)




#######

### covariance matrices 
############################
stage_samples <- split(rownames(meta_AF), meta_AF$Stage)
table(meta_AF$Stage)
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
  sigma_path <- file.path("dataset_antunes_ferreira", paste0("sigma_hat_", stage_name, ".csv"))
  X_path <- file.path("dataset_antunes_ferreira", paste0("X_", stage_name, ".csv"))
  
  # Write files
  write.csv(S_hat, file = sigma_path, row.names = TRUE)
  write.csv(mat_centered, file = X_path, row.names = TRUE)
}

gene_names <- gene_annotation$Symbol[match(rownames(expr_log2_filtr), gene_annotation$EnsemblGeneID)]
gene_names <- ifelse(is.na(gene_names), rownames(expr_log2_filtr), gene_names)
write.csv(data.frame(genes =gene_names ),file=file.path("dataset_antunes_ferreira", "genes_keep.csv") , row.names = FALSE)


