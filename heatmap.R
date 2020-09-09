heatmap_gene<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/发育过程分析/result/SERTOLI CELL/steroid hormone_gene.csv",header = F)
heatmap_gene<-as.character(heatmap_gene[,1])


#建立细胞分组信息
heatmap_idents<-pData(allsample_combined_cca_sertolicells_monocle)
heatmap_idents<-heatmap_idents[,c("age","State","Pseudotime")]
#提取表达矩阵
allsample_combined_cca_sertolicells_expr<-as.matrix(allsample_combined_cca_sertolicells[["RNA"]]@data)
heatmap_gene<-intersect(x=heatmap_gene,y=sertolicells_5percentexpr)#剔除矩阵中不存在的基因(为了避免scale后相对表达影响绝对表达，删除部分和线图不一样的基因)
#DELGENE<-c("LAMA2","CDH1","OCLN","DSC2","DSG2","IL1A","ICAM1")
#index<-heatmap_gene %in% DELGENE
#heatmap_gene<-heatmap_gene[!index]
#heatmap_gene<-c(heatmap_gene,"NECTIN3")
#length(heatmap_gene)
#提取所需基因，构建表达矩阵
heatmap_exprmatrix<-allsample_combined_cca_sertolicells_expr[heatmap_gene,]
#heatmap_exprmatrix<-scale(heatmap_exprmatrix)
heatmap_exprmatrix_preprocess<-preProcess(heatmap_exprmatrix)
heatmap_exprmatrix<-predict(heatmap_exprmatrix_preprocess,heatmap_exprmatrix)

dim(heatmap_exprmatrix)
heatmap_exprmatrix<-rbind(heatmap_exprmatrix,t(heatmap_idents$Pseudotime))#合并age属性
heatmap_exprmatrix<-rbind(heatmap_exprmatrix,t(heatmap_idents$State))#合并state属性
heatmap_exprmatrix<-t(heatmap_exprmatrix)
heatmap_exprmatrix<-heatmap_exprmatrix[order(heatmap_exprmatrix[,ncol(heatmap_exprmatrix)-1]),]#按照age重排矩阵
heatmap_exprmatrix<-heatmap_exprmatrix[order(heatmap_exprmatrix[,ncol(heatmap_exprmatrix)]),]#按照state重排矩阵
heatmap_exprmatrix<-heatmap_exprmatrix[,-ncol(heatmap_exprmatrix)]
heatmap_exprmatrix<-heatmap_exprmatrix[,-ncol(heatmap_exprmatrix)]#连续删除age和state列
heatmap_exprmatrix<-t(heatmap_exprmatrix)
write.csv(heatmap_exprmatrix,"heatmap_exprmatrix.csv")#为了去除矩阵的list属性
heatmap_exprmatrix<-read.csv("heatmap_exprmatrix.csv",header = T,row.names = 1)
heatmap_exprmatrix<-log2(heatmap_exprmatrix+2)

index<-heatmap_exprmatrix>3 #限定表达值范围
heatmap_exprmatrix[index]=3
index<-heatmap_exprmatrix<'-3'
heatmap_exprmatrix[index]=-3

breaks<-c(seq(-2,2,by=0.1))

#细胞分组信息重排
heatmap_idents<-heatmap_idents[order(heatmap_idents[,"Pseudotime"]),]
heatmap_idents<-heatmap_idents[order(heatmap_idents[,"State"]),]
#heatmap_idents<-heatmap_idents[,1:2]
row.names(heatmap_idents)<-row.names(t(heatmap_exprmatrix))
heatmap_idents$age<-factor(heatmap_idents$age,levels = c("2 years","5 years","8 years","11 years","17 years","adult"),ordered = T)

pheatmap(heatmap_exprmatrix, cellwidth = 0.02, cellheight = 1.5, fontsize = 8,
         method="spearman", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T, #为基因做聚类
         cluster_cols=F, #为sample做聚类
         color = colorRampPalette(c("steelblue", "white", "darkorange2"))(41), #自定义颜色
         breaks=breaks,#控制颜色范围
         show_colnames=F, #显示样品名称就改为T
         show_rownames =F, #显示基因名就改为T
         annotation_col = heatmap_idents,
         #如果想显示树状结构，就删掉下面这行前面的#
         treeheight_row = "0",treeheight_col = "0",#不画树
         #filename = "类固醇反应.pdf", #直接保存到pdf文件
         border_color = "NA") #不给每个小格子画边框，如果画边框，可以把NA改为你想要的颜色
