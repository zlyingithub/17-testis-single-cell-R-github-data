library(vegan)
library(ComplexHeatmap)
Xgene<-read.csv("chr_X_gene.csv",stringsAsFactors = FALSE)
Xgene<-Xgene[,1]
Ygene<-read.csv("chr_Y_gene.csv",stringsAsFactors = FALSE)
Ygene<-Ygene[,1]
#计算相似度
OAKSMATRIX<-read.csv("iNOA KS 6群细胞基因平均表达.csv",stringsAsFactors = FALSE, row.names = 1)
OAKSMATRIX<-read.csv("OA KS 9群细胞基因平均表达.csv",stringsAsFactors = FALSE, row.names = 1)
OAKSMATRIX<-OAKSMATRIX[Xgene,]
OAKSMATRIX<-t(OAKSMATRIX)
OAKSMATRIX[is.na(OAKSMATRIX)] <- 0
OAKS_jaccardDIS<- vegdist(OAKSMATRIX, method = 'jaccard', binary = TRUE)
OAKS_jaccardDIS<-as.matrix(OAKS_jaccardDIS)
OAKS_jaccardDIS=1-OAKS_jaccardDIS

#计算相关性
OAKSMATRIX<-read.csv("OA KS 9群细胞基因平均表达.csv",stringsAsFactors = FALSE, row.names = 1)
OAKSMATRIX<-OAKSMATRIX[Xgene,]
OAKSMATRIX[is.na(OAKSMATRIX)] <- 0
OAKSMATRIX_cor <- cor(OAKSMATRIX,method = c("spearman"))
OAKSMATRIX_cor_p <- matrix(0, nrow = ncol(OAKSMATRIX), ncol = ncol(OAKSMATRIX))
rownames(OAKSMATRIX_cor_p) <- colnames(OAKSMATRIX)
colnames(OAKSMATRIX_cor_p) <- colnames(OAKSMATRIX)

for (i in 1:ncol(OAKSMATRIX)){
  for (j in 1:ncol(OAKSMATRIX)){
    p <- cor.test(OAKSMATRIX[,i],OAKSMATRIX[,j],method = c("spearman"))
    OAKSMATRIX_cor_p[i,j] <- p$p.value
  }
}

corrplot(corr=OAKS_jaccardDIS)
# 定义右上部分图形的颜色
colCorRight <-  circlize::colorRamp2(c(0, 1), c("gray", "#ef3b2c"))
OAKS_jaccardDIS<-read.csv("各年龄SC细胞jaccard差异度.csv",stringsAsFactors = FALSE, row.names = 1)

## 绘制基本图形
Heatmap(OAKS_jaccardDIS, rect_gp = gpar(type = "none"), 
              show_heatmap_legend = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(col = "grey", fill = NA))
                if(i == j) {
                  grid.circle(x = x, y = y, r = 0.01 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
                } else {
                  grid.circle(x = x, y = y, r = (abs(OAKS_jaccardDIS[i, j])/2) * min(unit.c(width, height)), 
                              gp = gpar(fill = colCorRight(OAKSMATRIX_cor[i, j]), col = NA))
                }
              },
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = T, show_column_names = T, 
              row_names_side = "right", 
              row_names_rot = 0,
              row_names_gp = gpar(fontsize = 16),
              column_names_gp = gpar(fontsize = 16)
)
p2
lgdRight <- Legend(col_fun = colCorRight, title = "KICH", 
                   direction = "horizontal")
pd = list(lgdRight)
draw(p2, annotation_legend_list = pd,
     annotation_legend_side = "top")


inputdata<-read.csv("OA KS 相关性 相似性.csv")
library(ggplot2)
ggplot(inputdata,aes(x=Type,y=Cell))+geom_point(aes(size=Similarity,color=Correlation))+
  scale_size_continuous(range=c(10,15))+
  scale_colour_gradient(low="gray",high="red")+theme_bw()+theme(text = element_text(size=20))
p
