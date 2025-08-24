###IVPL_Analysis

###protein ranking plot
rm(list=ls())
library(ggplot2)
#final input file for MG：
##df<-read.csv("MGinputlog10mean.csv",header=T,check.names = F)
df$Rank <- rank(-(df$log10mean))
ggplot(df, aes(x = Rank, y = log10mean, label = Accession)) +
  geom_point(size = 1,color="#4472c4") +  # Add points
  # geom_text(vjust = 1.5, color = "blue") +  # Add gene symbols as text labels above the points
  labs(x = "Rank", y = "Log10 Average Abundance",
       title = "Protein Ranking Plot",size=8) +
  theme(axis.text=element_text(size=8,color = "black"),
        axis.title=element_text(size=8,color = "black"))+
  theme_minimal()+  # Use a minimal theme
  theme(axis.line = element_line(colour = "black",size=0.5))

###enrichment analysis visualization
library(ggplot2)
data<-read.csv("mgenrich.csv",header=T,check.names=F)

library(dplyr)
data <- data %>%
  mutate(ShortTerm = substr(Term, 1, 30))
data <- data %>%
  arrange(-log10P) %>%
  mutate(ShortTerm = factor(ShortTerm, levels = unique(ShortTerm))) 

ggplot(data, aes(x=enrichratio, y=ShortTerm))+
  geom_point(aes(size=number, fill=-log10P),shape=21,color="black")+
  theme_minimal() +theme(panel.background = element_blank(),  # Remove background
                         panel.border = element_rect(colour = "black", fill=NA, size=1),  # Add border around the plot
                         panel.grid.major = element_blank(),  # Remove major grid lines
                         panel.grid.minor = element_blank(),legend.position = "right",
                         # Customize plot title
                         axis.title = element_text(color = "black", size = 8),  # Customize axis titles
                         axis.text = element_text(color = "black", size = 8))+ 
  scale_fill_gradient(name="-log10P", low="yellow", high="red")


###GSEA visualization for Treg signature
data<-read.csv("TregGSEAwebinput.csv",header=T,check.names = F,row.names = 1)
geneList<-data$log2FC
names(geneList)<-rownames(data)
geneList<-sort(geneList,decreasing=T)
geneList<-geneList[geneList!=0]
head(geneList)

set<-read.csv("tregupcoregeneset.csv",header=T,check.names=F)
colnames(set)<-c("term","gene")

library(clusterProfiler)
set.seed(123456)
egmt<-GSEA(geneList, TERM2GENE=set,verbose=F,scoreType="pos",pvalueCutoff = 0.05)
gsea_res<-egmt@result

library(enrichplot)
gseaplot2(egmt, geneSetID = 1,title = "Treg core signature")

###GOplot
library(GOplot)
enrich<-read.csv("uniq991GOBP10terms.csv",header = T)
deg<-read.csv("uniq991deps2.csv",header = T,stringsAsFactors = F)
circ <- circle_dat(enrich, deg)
GOCircle(circ, nsub = 10, label.size = 3, rad1 = 3, rad2 = 4, table.legend = T)
process <- head(enrich$term,10)
chord <- chord_dat(data = circ, genes = deg, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 2, limit = c(0,2))
GOCluster(circ, process,clust.by = 'term', term.width = 2, lfc.max=3)
GOCluster(circ, process,clust.by = 'logFC', term.width = 2, lfc.max=9)

###z_score plot
oveovl<-read.csv("z score for seven glycolysis proteins new.csv",header=T,check.names = F)
d1 <- oveovl$protein
oveovl$protein <- seq(from=dim(oveovl)[1],to=1,by=-1)
###oveovl

#1
opar <- par(no.readonly=TRUE)
par(pin=c(4,12),mar=c(5, 15, 4, 8) + 0.1)

plot(x=oveovl[,2],y=oveovl[,1],yaxt="n",xlab="z score", ylab="",pch=20,cex=1.5,col="#D95F02",xlim = c(-3,3))
for (i in 5:7) {
  k <- i
  points(x=oveovl[,k],y=oveovl[,1],yaxt="n",xlab="", ylab="",pch=20,cex=1.5,col="#ff0000")
}
for (i in 2:4) {
  j <- i
  points(x=oveovl[,j],y=oveovl[,1],yaxt="n",xlab="", ylab="",pch=20,cex=1.5,col="#0000ff")
}
axis(2, at=oveovl$protein, labels=d1, col.axis="black",las=2,cex.axis=0.7)
for (k in 1:dim(oveovl)[1]) {
  abline(h = k,col = rgb(0, 0, 0, 50, maxColorValue=255),lty=3)
}
abline(v = 0,col = rgb(0, 0, 0, 50, maxColorValue=255),lty=1)
abline(v = 2,col = rgb(0, 0, 0, 50, maxColorValue=255),lty=1)
axis(side = 1, at = c(-5,-2,0,2,10,20,30,40,50))


###ligand<-read.csv("predicted_ligands.csv",header=T,check.names = F)
DB<-read.csv("PanglaoDB.csv",header=T,check.names = F)
library(dplyr)
result <- DB %>%
  semi_join(ligand, by = "gene_symbol") %>%  # 保留在 df_B 中出现的基因
  distinct(gene_symbol, cell_type)         # 去重（可选：避免重复行）

receptor<-read.csv("edgesreceptortop20input.csv",header=T,check.names = F)
library(dplyr)
merged_df <- receptor %>%
  left_join(result, by = c("ligand" = "gene_symbol")) %>%
  rename(ligand_cell_type = cell_type)
merged_df <- merged_df %>%
  mutate(receptor_cell_type = "Treg")

merged_df_clean <- merged_df[!is.na(merged_df$ligand_cell_type), ]

library(ggplot2)
library(ggforce)
library(dplyr)

# ligand_gene_symbol, receptor_gene_symbol, ligand_cell_type, receptor_cell_type
df<-read.csv("ligand receptor cell type by PanglaoDB omit impossible cell type.csv",header=T,check.names = F)
# 简化为 cell-cell 连接 + annotation
df_clean <- df %>%
  filter(!is.na(ligand_cell_type)) %>%
  mutate(pair = paste0(ligand_gene_symbol, "-", receptor_gene_symbol)) %>%
  group_by(ligand_cell_type, receptor_cell_type, pair) %>%
  summarise(n = n(), .groups = "drop")

# 所有细胞类型（左边为ligand source，右边为Treg）
cell_types <- unique(c(df_clean$ligand_cell_type, df_clean$receptor_cell_type))

# 分配坐标（左边细胞类型x=1，右边x=3）
node_df <- data.frame(
  cell_type = unique(c(df_clean$ligand_cell_type, df_clean$receptor_cell_type)),
  x = ifelse(unique(c(df_clean$ligand_cell_type, df_clean$receptor_cell_type)) == "Treg", 3, 1)
) %>%
  arrange(x, cell_type) %>%
  mutate(y = seq(from = 1, to = length(cell_type)))

# 合并坐标信息
df_plot <- df_clean %>%
  left_join(node_df, by = c("ligand_cell_type" = "cell_type")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(node_df, by = c("receptor_cell_type" = "cell_type")) %>%
  rename(x2 = x, y2 = y)


ggplot() +
  # 曲线连线
  geom_curve(data = df_plot,
             aes(x = x1, y = y1, xend = x2, yend = y2, color = pair),
             curvature = 0.2, alpha = 0.8, size = 1) +
  
  # 节点标签
  geom_text(data = node_df,
            aes(x = x, y = y, label = cell_type),
            hjust = ifelse(node_df$x == 1, 1.1, -0.1),
            size = 4) +
  
  xlim(0.5, 3.5) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Cell-cell communication via ligand-receptor pairs")

node_df$y[node_df$cell_type=="Treg"]<-8
df_plot$y2<-16

# 添加弯曲方向列
df_plot <- df_plot %>%
  mutate(curve_direction = ifelse(y1 > y2, "down", "up"))

# 上弯（curvature = 0.3）
ggplot() +
  geom_curve(
    data = df_plot %>% filter(curve_direction == "up"),
    aes(x = x1, y = y1, xend = x2, yend = y2, color = pair),
    curvature = 0.2, alpha = 0.8, size = 1
  ) +
  
  # 下弯（curvature = -0.3）
  geom_curve(
    data = df_plot %>% filter(curve_direction == "down"),
    aes(x = x1, y = y1, xend = x2, yend = y2, color = pair),
    curvature = -0.2, alpha = 0.8, size = 1
  ) +
  
  geom_text(data = node_df,
            aes(x = x, y = y, label = cell_type),
            hjust = ifelse(node_df$x == 1, 1.1, -0.1),
            size = 4) +
  
  xlim(0.5, 3.5) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Cell–cell communication via ligand–receptor pairs")

#ligand enrichment
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(readr)

go_results_list <- list()
cell_types <- unique(df$ligand_cell_type)

for (cell in cell_types) {
  gene_list <- unique(df$ligand_gene_symbol[df$ligand_cell_type == cell])
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  if (!is.null(entrez_ids) && nrow(entrez_ids) > 0) {
    ego <- enrichGO(gene = entrez_ids$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    ont = "BP",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE)
    if (!is.null(ego) && nrow(ego) > 0) {
      ego_df <- as.data.frame(ego)
      ego_df$CellType <- cell
      go_results_list[[cell]] <- ego_df
    }
  }
}

go_all <- bind_rows(go_results_list)

####exosome analysis
library(tidyverse)
df<-read.csv("EXO-spectronaut.csv",header=T,check.names = F)
updf<-read.csv("uptakeexo.csv",header=T,check.names = F)
upTNBC <- updf[!(updf$TNBC_1 == 0 & updf$TNBC_2 == 0 & updf$TNBC_3 == 0), ]
upTNBC_proteins<-upTNBC[[2]]
df[ , 4:7] <- lapply(df[ , 4:7], function(x) {
  x[x == "NaN"] <- 0       
  as.numeric(x)            
})
df <- df[which(rowSums(df[,4:7]) > 0),] #删除行全为0的行

TNBC2<- df[df$`TNBCEXO-1` != 0 & df$`TNBCEXO-2` != 0, ]
TNBC_proteins<-TNBC2[[2]]
intersect_genes <- intersect(upTNBC_proteins, TNBC_proteins)

diffinput_only <- setdiff(upTNBC_proteins,TNBC_proteins)
diffuptake_only <- setdiff(TNBC_proteins,upTNBC_proteins)
cat("Input proteins: ", length(TNBC_proteins), "\n") 
cat("Uptake proteins: ", length(upTNBC_proteins), "\n")
cat("Overlap: ", length(intersect_genes), "\n")
cat("Uptake-specific: ", length(diffuptake_only), "\n")

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list("TNBC-exo Input" = TNBC_proteins, "Uptake Proteome" = upTNBC_proteins),
  filename = NULL,
  fill = c("skyblue", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)
grid::grid.newpage()
grid::grid.draw(venn.plot)

diffupdf<-read.csv("diffuptakeexo.csv",header=T,check.names = F)
diffupTNBC<-diffupdf[[2]]
intersect_genes2 <- intersect(diffinput_only, diffupTNBC) 
IVPLTNBCexouniq151<-as.data.frame(intersect_genes2)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list("Sigup TNBCexo" = diffupTNBC, "Unique Uptake TNBCexo" = diffinput_only ),
  filename = NULL,
  fill = c("skyblue", "tomato"),
  alpha = 0.6,
  cex = 1.5,
  cat.cex = 1.5,
  cat.pos = 0
)
grid::grid.newpage()
grid::grid.draw(venn.plot)


##MA plot
df$Avalue <- log2(df$Mean)  
library(ggplot2)
library(ggrepel)
top10 <- df[order(-df$log2FC), ][1:10, ]

ggplot(df, aes(x = Avalue, y = log2FC)) +
  geom_point(aes(color = log2FC), alpha = 0.7) +
  scale_color_gradient(low = "skyblue", high = "firebrick") +
  geom_text_repel(data = top10, aes(label = intersect_genes),
                  size = 3.5, color = "black", max.overlaps = Inf) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  labs(
    x = "log2(Mean Abundance)", 
    y = "log2FC", 
    title = ""
  ) +
  theme_minimal()+
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )


mat_diff <- as.matrix(sigdf[ , 4:7])

rownames(mat_diff) <- sigdf[[2]]
mat_diff_scaled <- t(scale(t(mat_diff)))

colnames(mat_diff_scaled)
# [1] "NTexo1" "NTexo2" "TNBCexo1" "TNBCexo2"
group <- data.frame(Group = c("NT", "NT", "TNBC", "TNBC"))
rownames(group) <- colnames(mat_diff_scaled)
ann_colors <- list(
  Group = c(NT = "blue", TNBC = "red")
)

library(pheatmap)
pheatmap(mat_diff_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,annotation_col = group,       # ✅ 添加分组条
         annotation_colors = ann_colors,  # ✅ 设置颜色（可选）
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "TNBCexo vs NTexo Heatmap")



#enrichment
selected_names <- df[df[[9]] >= 1.5, 2] 
library(clusterProfiler)
library(org.Hs.eg.db)  
library(AnnotationDbi)

gene_symbols <- unique(selected_names)

# 转换为 ENTREZ ID
gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)
ego_bp <- enrichGO(gene          = gene_entrez$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = "ENTREZID",
                   ont           = "BP",           # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2,
                   readable      = TRUE)  # 输出结果中带 SYMBOL 名


gene_ratio <- as.numeric(sapply(strsplit(ego_bp$GeneRatio, "/"), `[`, 1)) /
  as.numeric(sapply(strsplit(ego_bp$GeneRatio, "/"), `[`, 2))

bg_ratio <- as.numeric(sapply(strsplit(ego_bp$BgRatio, "/"), `[`, 1)) /
  as.numeric(sapply(strsplit(ego_bp$BgRatio, "/"), `[`, 2))
ego_bp<- as.data.frame(ego_bp)
ego_bp$EnrichmentRatio <- gene_ratio / bg_ratio
ego_bp_simplified <- simplify(ego_bp)
gene_ratio <- as.numeric(sapply(strsplit(ego_bp_simplified$GeneRatio, "/"), `[`, 1)) /
  as.numeric(sapply(strsplit(ego_bp_simplified$GeneRatio, "/"), `[`, 2))

bg_ratio <- as.numeric(sapply(strsplit(ego_bp_simplified$BgRatio, "/"), `[`, 1)) /
  as.numeric(sapply(strsplit(ego_bp_simplified$BgRatio, "/"), `[`, 2))
ego_bp_simplified<- as.data.frame(ego_bp_simplified)

ego_bp_simplified$EnrichmentRatio <- gene_ratio / bg_ratio

library(rrvgo)
library(org.Hs.eg.db)
go_ids <- ego_bp$ID
scores <- setNames(-log10(ego_bp$p.adjust), go_ids)  # 用 p.adjust 作为评分

str(scores)

simMatrix <- calculateSimMatrix(go_ids,
                                orgdb = "org.Hs.eg.db",
                                ont = "BP",
                                method = "Rel")


reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold = 0.7,
                                orgdb = "org.Hs.eg.db")

treemapPlot(reducedTerms)

###sorted vs IVPL
sorted<-read.csv("sorted.csv",header=T,check.names = F)
IVPL<-read.csv("IVPL.csv",header=T,check.names = F)
length(intersect(IVPL$`Gene Symbol`, sorted$`Gene Symbol`)) / length(sorted$`Gene Symbol`)
overlap<-intersect(IVPL$`Gene Symbol`, sorted$`Gene Symbol`)

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    IVPLInt = IVPL$`Gene Symbol`,
    sortedInt = sorted$`Gene Symbol`
  ),
  filename = NULL,
  fill = c("lightpink","skyblue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.pos = 0
)

grid.newpage()
grid.draw(venn.plot)

intersect_IVPL3828 <- IVPL[IVPL$`Gene Symbol` %in% overlap, ]
intersect_sorted3828 <- sorted[sorted$`Gene Symbol` %in% overlap, ]

IVPL3828 <- intersect_IVPL3828[, c("Gene Symbol","Int1","Int2")]
sorted3828 <-intersect_sorted3828[, c("Gene Symbol","sorted1","sorted2")]
all3828 <- merge(IVPL3828, sorted3828, by = "Gene Symbol")
#all1875[is.na(all1875)] <- 0
##归一化处理
cols <- c("Int1","Int2", "sorted1","sorted2")
all3828[cols] <- sweep(all3828[cols], 2, colSums(all3828[cols], na.rm = TRUE), FUN = "/")

df<-all3828
expr_mat <- df[, (ncol(df)-3):ncol(df)]

rownames(expr_mat) <- df$Gene_Symbol

expr_mat <- expr_mat[complete.cases(expr_mat), ]

sample_cor <- cor(expr_mat, method = "pearson")

library(pheatmap)

pheatmap(sample_cor,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         main = "Sample-wise Pearson Correlation of Protein Expression")


expr_mat[, (ncol(expr_mat)-3):ncol(expr_mat)] <- expr_mat[, (ncol(expr_mat)-3):ncol(expr_mat)] * 1e5
expr_mat[expr_mat == 0] <- 0.01
expr_matlog2<- log2(expr_mat)

sample_cor <- cor(expr_matlog2, method = "pearson")

pheatmap(sample_cor,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         main = "Sample-wise Pearson Correlation of Protein Expression")

expr_mat$Int_avg <- rowMeans(expr_mat[, c("Int1", "Int2")], na.rm = TRUE)
expr_mat$sorted_avg <- rowMeans(expr_mat[, c("sorted1", "sorted2")], na.rm = TRUE)
expr_mat$log2FC <- log2(expr_mat$Int_avg ) - log2(expr_mat$sorted_avg)  

cor(log2(expr_mat$sorted_avg), log2(expr_mat$Int_avg), method = "pearson")
cor_result <- cor.test(log2(expr_mat$sorted_avg), log2(expr_mat$Int_avg), method = "pearson")
r_val <- round(cor_result$estimate, 3)  
p_val <- signif(cor_result$p.value, 3)  
plot(log2(expr_mat$sorted_avg), log2(expr_mat$Int_avg),
     xlab = "log2(Sorted Intestine Avg)", ylab = "log2(IVPL Intestine Avg)",
     main = "Correlation of Expression Trends")
abline(lm(log2(Int_avg) ~ log2(sorted_avg), data = expr_mat), col = "red")
text(x = max(log2(expr_mat$sorted_avg)) * 0.5,
     y = max(log2(expr_mat$Int_avg)) * 0.95,
     labels = paste0("r = ", r_val, ", p ", format.pval(p_val, eps = 0.0001, digits = 2)),
     cex = 1.2)
dev.new()

library(ggplot2)
library(ggpmisc)
library(viridis)  
library(ggpubr)
# 计算log2值
expr_mat$log_sorted <- log2(expr_mat$sorted_avg)
expr_mat$log_Int <- log2(expr_mat$Int_avg)

ggplot(expr_mat, aes(x = log_sorted, y = log_Int)) +
  geom_point(aes(color = log_sorted), alpha = 0.6, size = 1.5) +
  scale_color_viridis(option = "plasma", direction = -1) +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE, size = 1) +
  stat_cor(
    method = "pearson",
    label.x = min(expr_mat$log_sorted) + 0.5,
    label.y = max(expr_mat$log_Int) - 0.5,
    size = 5,
    color = "black",
    parse = FALSE
  ) +
  labs(
    title = "Correlation of IVPL and Sorted Intestine",
    x = "log2(Sorted Intestine Average Abundance)",
    y = "log2(IVPL Intestine Average Abundance)",
    color = "log2(Abundance)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )
