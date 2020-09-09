setwd("F:/学习/数据分析/single cell sequence/测序数据/human/发育过程分析/result/3病理综合")
save.image(file = "6age加3病理综合allcells.Rdata")
save.image(file = "adult加3病理综合allcells_mastcells.Rdata")
library(cowplot)
library(Seurat)
library(dplyr)
library(monocle)
library(stringr)
library(pheatmap)
library(vegan)
library(biomaRt)
#建立seurat文件
#规范命名规则，project一律小写，细胞名前一律小写，xxx$sample中sample小写，sample内的标号大写。
#新建的seurat文件没有进行质量筛选
lz013_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/BD/SEURAT/SampleTag_04_SC089_DBEC_MolsPerCell.csv",row.names=1,header=T)
lz013_expression<-t(lz013_expression)
lz013 <- CreateSeuratObject(counts = lz013_expression, project = "lz013")
lz013$sample <- "LZ013"
lz013$tech<-"10X"
lz013$age<-"adult"
lz013[["percent.mt"]] <- PercentageFeatureSet(object = lz013, pattern = "^MT.")
remove(lz013_expression)

lz014_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/BD/SEURAT/SampleTag_05_SC089_DBEC_MolsPerCell.csv",row.names=1,header=T)
lz014_expression<-t(lz014_expression)
lz014 <- CreateSeuratObject(counts = lz014_expression, project = "lz014")
lz014$sample <- "LZ014"
lz014$tech<-"10X"
lz014$age<-"adult"
lz014[["percent.mt"]] <- PercentageFeatureSet(object = lz014, pattern = "^MT.")
remove(lz014_expression)


lz015_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/BD/SEURAT/SampleTag_06_SC089_DBEC_MolsPerCell.csv",row.names=1,header=T)
lz015_expression<-t(lz015_expression)
lz015 <- CreateSeuratObject(counts = lz015_expression, project = "lz015")
lz015$sample <- "LZ015"
lz015$tech<-"10X"
lz015$age<-"adult"
lz015[["percent.mt"]] <- PercentageFeatureSet(object = lz015, pattern = "^MT.")
remove(lz015_expression)

#lz005
lz005_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ_005_L/filtered_feature_bc_matrix")
lz005 <- CreateSeuratObject(counts = lz005_expr, project = "LZ005")
lz005$sample <- "LZ005"
lz005$tech<-"10X"
remove(lz005_expr)
lz005$age<-"5 years"
lz005[["percent.mt"]] <- PercentageFeatureSet(object = lz005, pattern = "^MT-")

#lz003
lz003_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ003/P101SC18111759-01-F007-B4-37_10X_results/2.analysis/LZ003/outs/filtered_feature_bc_matrix")
lz003 <- CreateSeuratObject(counts = lz003_expr, project = "LZ003")
lz003$sample <- "LZ003"
lz003$age<-"adult"
lz003$tech<-"10X"
remove(lz003_expr)
lz003[["percent.mt"]] <- PercentageFeatureSet(object = lz003, pattern = "^MT-")

#lz007
lz007_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ007/filtered_feature_bc_matrix")
lz007 <- CreateSeuratObject(counts = lz007_expr, project = "LZ007")
lz007$sample <- "LZ007"
lz007$age<-"adult"
lz007$tech<-"10X"
remove(lz007_expr)
lz007[["percent.mt"]] <- PercentageFeatureSet(object = lz007, pattern = "^MT-")

#lz008
lz008_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ008/P101SC18111759-01-F004-B9-37_10X_results/2.analysis/LZ008/outs/filtered_feature_bc_matrix")
lz008 <- CreateSeuratObject(counts = lz008_expr, project = "LZ008")
lz008$sample <- "LZ008"
lz008$age<-"11 years"
lz008$tech<-"10X"
remove(lz008_expr)
lz008[["percent.mt"]] <- PercentageFeatureSet(object = lz008, pattern = "^MT-")

#lz009
lz009_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ009/outs/filtered_feature_bc_matrix")
lz009 <- CreateSeuratObject(counts = lz009_expr, project = "LZ009")
lz009$sample <- "LZ009"
lz009$age<-"8 years"
lz009$tech<-"10X"
remove(lz009_expr)
lz009[["percent.mt"]] <- PercentageFeatureSet(object = lz009, pattern = "^MT-")

#lz011
lz011_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ011/filtered_feature_bc_matrix")
lz011 <- CreateSeuratObject(counts = lz011_expr, project = "LZ011")
lz011$sample <- "LZ011"
lz011$age<-"2 years"
lz011$tech<-"10X"
remove(lz011_expr)
lz011[["percent.mt"]] <- PercentageFeatureSet(object = lz011, pattern = "^MT-")
remove(lz011_expr)

#lz016
lz016_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ016/filtered_feature_bc_matrix")
lz016 <- CreateSeuratObject(counts = lz016_expr, project = "LZ016")
lz016$sample <- "LZ016"
lz016$age<-"17 years"
lz016$tech<-"10X"
remove(lz016_expr)
lz016[["percent.mt"]] <- PercentageFeatureSet(object = lz016, pattern = "^MT-")
remove(lz016_expr)

#lz002
lz002_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ002/P101SC18111759-01-F006-B3-37_10X_results/2.analysis/LZ002/outs/filtered_feature_bc_matrix")
lz002 <- CreateSeuratObject(counts = lz002_expr, project = "LZ002")
lz002$sample <- "LZ002"
lz002$age<-"AZFa_Del"
lz002$tech<-"10X"
lz002[["percent.mt"]] <- PercentageFeatureSet(object = lz002, pattern = "^MT-")
remove(lz002_expr)

lz017_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/LZ017LZ018LZ019/07 样本标记分析/C 样本独立文件夹/SampleTag_05_SC162_DBEC_MolsPerCell.csv",row.names=1,header=T)
row.names(lz017_expression)<-paste0("LZ017",row.names(lz017_expression))
lz017_expression<-t(lz017_expression)
lz017 <- CreateSeuratObject(counts = lz017_expression, project = "lz017")
lz017$sample <- "LZ017"
lz017$age<-"iNOA"
lz017$tech<-"BD"
lz017[["percent.mt"]] <- PercentageFeatureSet(object = lz017, pattern = "^MT.")
remove(lz017_expression)

lz018_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/lz017LZ018LZ019/07 样本标记分析/C 样本独立文件夹/SampleTag_06_SC162_DBEC_MolsPerCell.csv",row.names=1,header=T)
row.names(lz018_expression)<-paste0("LZ018",row.names(lz018_expression))
lz018_expression<-t(lz018_expression)
lz018 <- CreateSeuratObject(counts = lz018_expression, project = "lz018")
lz018$sample <- "LZ018"
lz018$age<-"iNOA"
lz018$tech<-"BD"
lz018[["percent.mt"]] <- PercentageFeatureSet(object = lz018, pattern = "^MT.")
remove(lz018_expression)

lz019_expression<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/lz017lz018LZ019/07 样本标记分析/C 样本独立文件夹/SampleTag_07_SC162_DBEC_MolsPerCell.csv",row.names=1,header=T)
row.names(lz019_expression)<-paste0("LZ019",row.names(lz019_expression))
lz019_expression<-t(lz019_expression)
lz019 <- CreateSeuratObject(counts = lz019_expression, project = "lz019")
lz019$sample <- "LZ019"
lz019$age<-"iNOA"
lz019$tech<-"BD"
lz019[["percent.mt"]] <- PercentageFeatureSet(object = lz019, pattern = "^MT.")
remove(lz019_expression)

#lz004
lz004_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ004/P101SC18111759-01-F003-B8-18_10X_results/2.analysis/LZ004/filtered_feature_bc_matrix")
lz004 <- CreateSeuratObject(counts = lz004_expr, project = "lZ004")
lz004$sample <- "LZ004"
lz004$age<- "KS"
lz004$tech<-"10X"
remove(lz004_expr)
lz004[["percent.mt"]] <- PercentageFeatureSet(object = lz004, pattern = "^MT-")

#lz010
lz010_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ010/outs/filtered_feature_bc_matrix")
lz010 <- CreateSeuratObject(counts = lz010_expr, project = "lZ010")
lz010$sample <- "LZ010"
lz010$age<-"KS"
lz010$tech<-"10X"
remove(lz010_expr)
lz010[["percent.mt"]] <- PercentageFeatureSet(object = lz010, pattern = "^MT-")

#lz012
lz012_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/LZ012/P101SC18111759-01-F013-B13-37_10X_results/2.analysis/LZ012/outs/filtered_feature_bc_matrix")
lz012 <- CreateSeuratObject(counts = lz012_expr, project = "lz012")
lz012$sample <- "LZ012"
lz012$age<-"KS"
lz012$tech<-"10X"
remove(lz012_expr)
lz012[["percent.mt"]] <- PercentageFeatureSet(object = lz012, pattern = "^MT-")

allsamples_anchors <- FindIntegrationAnchors(object.list = list(lz002,lz003,lz005,lz007,lz008,lz009,lz011,lz013,lz014,lz015,lz016,lz017,lz018,lz019,lz004,lz010,lz012), dims = 1:20)
allsample_combined_cca <- IntegrateData(anchorset = allsamples_anchors, dims = 1:20)

mypalette<-c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")
mypalette<-c("#e78ac3","#66c2a5","#fc8d62","#8da0cb","#e78ac3")
#质量控制
allsample_combined_cca <- subset(x = allsample_combined_cca, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 40 & nCount_RNA <80000)
#批次效应处理
allsample_combined_cca_list <- SplitObject(allsample_combined_cca, split.by = "tech")
for (i in 1:length(allsample_combined_cca_list)) {
  allsample_combined_cca_list[[i]] <- NormalizeData(allsample_combined_cca_list[[i]], verbose = FALSE)
  allsample_combined_cca_list[[i]] <- FindVariableFeatures(allsample_combined_cca_list[[i]], selection.method = "vst", 
                                                           nfeatures = 1000, verbose = FALSE)
}
allsample_combined_cca_anchors <- FindIntegrationAnchors(object.list = allsample_combined_cca_list, dims = 1:30,anchor.features = 1000)
allsample_combined_cca <- IntegrateData(anchorset = allsample_combined_cca_anchors, dims = 1:30)
#标准化
#DefaultAssay(allsample_combined_cca) <- "RNA"
allsample_combined_cca <- ScaleData(allsample_combined_cca)
allsample_combined_cca <- RunPCA(allsFindClustersFindClustersample_combined_cca, npcs = 30)
#聚类
allsample_combined_cca <- FindNeighbors(object = allsample_combined_cca, dims = 1:30)
allsample_combined_cca <- FindClusters(object = allsample_combined_cca, resolution = 0.2)
#降维
allsample_combined_cca <- RunUMAP(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)
allsample_combined_cca <- RunTSNE(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)

# Visualization
DimPlot(object = allsample_combined_cca, reduction = "tsne", group.by = "age",label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "umap", group.by = "agetype",label.size = 8,pt.size = 1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "tsne",  label = F,label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

FeaturePlot(allsample_combined_cca, features = c( "FGFR3","SYCP3", "TNP1", "SOX9", "DLK1","MYH11","VWF","NOTCH2", "CD163","TPSAB1","CD3D"), cols = c("blue", "yellow"),min.cutoff = "q9",pt.size = 0.1, reduction = "tsne") #SSC
FeaturePlot(allsample_combined_cca, features = c( "DDX4","FGFR3","SYCP3", "TNP1","VIM", "SOX9", "DLK1","MYH11","VWF","NOTCH3","PTPRC", "CD163","TPSAB1","CD3D"), cols = c("blue", "yellow"),min.cutoff = "q9",pt.size = 0.1, reduction = "tsne") #SSC
VlnPlot(allsample_combined_cca_10X_KSandOA_germcell,features = c("VWF"),pt.size = 0.1,group.by = "seurat_GERMclusters",split.by = "age")+scale_fill_manual(values =mypalette)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 30),plot.title =element_text(size = 45,face="plain") )
VlnPlot(allsample_combined_cca_10X_OA,features = c("TIMP1"),pt.size = 0,group.by = "seurat_11clusters",cols = mypalette)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )
#带趋势线
VlnPlot(allsample_combined_cca_germcells,features = c("percent.glycolysis_gene"),pt.size = 0.1)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )+stat_summary(aes(group=1),fun.y=mean, geom="smooth", shape=1,size=2,color="gray")+ scale_y_continuous(limits = c(0,2))

allsample_combined_cca<-subset(x = allsample_combined_cca, idents = c("0","1","2","3","4","5","6","7","8","9","11","12","13","14","15","16"))
allsample_combined_cca <- RenameIdents(object = allsample_combined_cca, `0` = "MIX_cells", `1` = "MIX_cells",`2` = "spermatid",
                                       `3` = "spermatogonia",`4` = "Sertoli_cells",`5` = "Sertoli_cells",
                                       `6` = "spermatocyte",`7` = "spermatocyte",`8` = "spermatocyte",`9` = "macrophages",`11` = "VSM_cells",
                                       `12` = "endotheliocyte",`13` = "spermatid",`14` = "mast_cells",`15` = "T_cells",`16` = "MIX_cells")
allsample_combined_cca$seurat_noa9clusters<-Idents(allsample_combined_cca)
allsample_combined_cca$age<-factor(allsample_combined_cca$age,levels = c("OA","iNOA","KS","AZFa_Del"),ordered = F)
allsample_combined_cca$seurat_11clusters<-factor(allsample_combined_cca$seurat_11clusters,levels = c("SPGs","SPCs","SPTs","SCs","LCs","PTMs","ECs","VSMs","Mo&Mφ","MCs","TCs"),ordered = F)
MARKERS <- FindAllMarkers(object = allsample_combined_cca_10X_KS, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
MARKERS <- FindMarkers(object = allsample_combined_cca_10X_OA, only.pos = F,ident.1 = "SCs", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(MARKERS,"KS 内部11群DEGs.csv")
#计算基因平均值
MARKERS <- AverageExpression(allsample_combined_cca, return.seurat = T)
#计算metadata里的平均值（表格，按照xx分组list属性，平均值或其他）
MARKERS<-aggregate(allsample_combined_cca_10X_KSandOA@meta.data,list(allsample_combined_cca_10X_KSandOA@meta.data$seurat_11clusters,allsample_combined_cca_10X_KSandOA@meta.data$age),mean)

#代谢分析
#读取基因
X_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_X_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)
Y_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_Y_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)
eXi_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/X染色体失活逃逸基因.csv",header = F,sep = "\t",,stringsAsFactors=F)
X_gene<-X_gene[,1]
Y_gene<-Y_gene[,1]
eXi_gene<-eXi_gene[,1]
OXPHOS_genes<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/氧化磷酸化-无线粒体基因 go 基因列表.csv")
glycolysis_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/经典糖酵解 go 基因列表.csv")
triglyceride_metabolic_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/甘油三酯代谢.csv")
OXPHOS_genes<-as.character(OXPHOS_genes[,1])
glycolysis_gene<-as.character(glycolysis_gene[,1])
triglyceride_metabolic_gene<-as.character(triglyceride_metabolic_gene[,1])
#计算X，Y染色体基因表达，注意这里使用的genelist应该都要在object的RNA@meta.features里面存在，所以应当筛选
allgene<-row.names(allsample_combined_cca@assays$RNA@meta.features)
OXPHOS_genes<-intersect(x=OXPHOS_genes, y = allgene)
glycolysis_gene<-intersect(x=glycolysis_gene, y = allgene)
triglyceride_metabolic_gene<-intersect(x=triglyceride_metabolic_gene, y = allgene)
X_gene<-subset(X_gene,X_gene %in% allgene)
Y_gene<-subset(Y_gene,Y_gene %in% allgene)
eXi_gene<-subset(eXi_gene,eXi_gene %in% allgene)
Xi_gene<-setdiff(X_gene, eXi_gene)
#计算比例
allsample_combined_cca[["percent.glycolysis_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= glycolysis_gene)
allsample_combined_cca[["percent.oxidative_phosphorylation_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= OXPHOS_genes)
allsample_combined_cca[["percent.triglyceride_metabolic_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)
allsample_combined_cca[["percent.X_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)
allsample_combined_cca[["percent.Y_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)
allsample_combined_cca[["percent.eXi_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)
allsample_combined_cca[["percent.Xi_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)


#基因表达相关性回归分析
MARKER<-allsample_combined_cca_10X_KS_SC@meta.data
MARKERS<-allsample_combined_cca_10X_KS_SC@assays$RNA@data
MARKERS<-data.frame(MARKERS)
MARKERS<-t(MARKERS)
MARKER<-cbind(MARKER,MARKERS)

ggplot(data = MARKER, mapping = aes(x = XIST, y = percent.chr_neXi_gene))+ 
  geom_point(aes(color = percent.chr_X_gene,size= nCount_RNA))+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+labs(caption ="y = 4.092-0.698x  p-value = 8.185e-14")
MARKERS<-lm(percent.chr_eXi_gene~XIST,MARKER)
summary(MARKERS) 

#绘制火山图
MARKER<-read.csv("KS SSC 所有基因与all X回归分析.csv",header = T,stringsAsFactors=F)

ggplot(data = MARKER, aes(x = r, y = -log10(p.value),size=mean,color=UPDOWN)) +
  geom_point(alpha=0.6)+scale_size(limits=c(0,4)) + scale_color_manual(values=c("#90bff9", "grey","#f2b77c"))+
  xlim(c(-1.1, 0.7)) +
  ylim(c(0, 4.7)) +
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8)


#各组平均表达热图
MARKERS <- AverageExpression(allsample_combined_cca_10X_OA, return.seurat = T)
DoHeatmap(MARKERS, features = c("USP9Y","DDX3Y","UTY","HSFY1","HSFY2","RBMY1A1","RBMY1B","RBMY1C","RBMY1D","RBMY1E","RBMY1J","DAZ1","DAZ2","DAZ3","DAZ4","BPY2","CDY1B","PRY","CSPG4P1Y","DAZL","BOLL","SRY","ZFX","ZFY"), size = 3,slot = "data", disp.min = 0, disp.max = 5)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")
DoHeatmap(MARKERS, features = c("AMH","INHBB","INHA"), size = 6)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")

#计算多个基因打分
IMMUNE_MARKER<-list(c("PTPRC","TPSAB1","CD163","CD3D","CD4","CD8A","MS4A1", "CD79A","GNLY","ITGAX"))
lz024<-AddModuleScore(object = lz024,features = IMMUNE_MARKER,name = 'IMMUNE_Features')