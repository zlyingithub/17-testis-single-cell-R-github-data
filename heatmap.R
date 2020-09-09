heatmap_gene<-read.csv("F:/ѧϰ/���ݷ���/single cell sequence/��������/human/�������̷���/result/SERTOLI CELL/steroid hormone_gene.csv",header = F)
heatmap_gene<-as.character(heatmap_gene[,1])


#����ϸ��������Ϣ
heatmap_idents<-pData(allsample_combined_cca_sertolicells_monocle)
heatmap_idents<-heatmap_idents[,c("age","State","Pseudotime")]
#��ȡ�������
allsample_combined_cca_sertolicells_expr<-as.matrix(allsample_combined_cca_sertolicells[["RNA"]]@data)
heatmap_gene<-intersect(x=heatmap_gene,y=sertolicells_5percentexpr)#�޳������в����ڵĻ���(Ϊ�˱���scale����Ա���Ӱ����Ա��ɾ�����ֺ���ͼ��һ���Ļ���)
#DELGENE<-c("LAMA2","CDH1","OCLN","DSC2","DSG2","IL1A","ICAM1")
#index<-heatmap_gene %in% DELGENE
#heatmap_gene<-heatmap_gene[!index]
#heatmap_gene<-c(heatmap_gene,"NECTIN3")
#length(heatmap_gene)
#��ȡ������򣬹����������
heatmap_exprmatrix<-allsample_combined_cca_sertolicells_expr[heatmap_gene,]
#heatmap_exprmatrix<-scale(heatmap_exprmatrix)
heatmap_exprmatrix_preprocess<-preProcess(heatmap_exprmatrix)
heatmap_exprmatrix<-predict(heatmap_exprmatrix_preprocess,heatmap_exprmatrix)

dim(heatmap_exprmatrix)
heatmap_exprmatrix<-rbind(heatmap_exprmatrix,t(heatmap_idents$Pseudotime))#�ϲ�age����
heatmap_exprmatrix<-rbind(heatmap_exprmatrix,t(heatmap_idents$State))#�ϲ�state����
heatmap_exprmatrix<-t(heatmap_exprmatrix)
heatmap_exprmatrix<-heatmap_exprmatrix[order(heatmap_exprmatrix[,ncol(heatmap_exprmatrix)-1]),]#����age���ž���
heatmap_exprmatrix<-heatmap_exprmatrix[order(heatmap_exprmatrix[,ncol(heatmap_exprmatrix)]),]#����state���ž���
heatmap_exprmatrix<-heatmap_exprmatrix[,-ncol(heatmap_exprmatrix)]
heatmap_exprmatrix<-heatmap_exprmatrix[,-ncol(heatmap_exprmatrix)]#����ɾ��age��state��
heatmap_exprmatrix<-t(heatmap_exprmatrix)
write.csv(heatmap_exprmatrix,"heatmap_exprmatrix.csv")#Ϊ��ȥ�������list����
heatmap_exprmatrix<-read.csv("heatmap_exprmatrix.csv",header = T,row.names = 1)
heatmap_exprmatrix<-log2(heatmap_exprmatrix+2)

index<-heatmap_exprmatrix>3 #�޶�����ֵ��Χ
heatmap_exprmatrix[index]=3
index<-heatmap_exprmatrix<'-3'
heatmap_exprmatrix[index]=-3

breaks<-c(seq(-2,2,by=0.1))

#ϸ��������Ϣ����
heatmap_idents<-heatmap_idents[order(heatmap_idents[,"Pseudotime"]),]
heatmap_idents<-heatmap_idents[order(heatmap_idents[,"State"]),]
#heatmap_idents<-heatmap_idents[,1:2]
row.names(heatmap_idents)<-row.names(t(heatmap_exprmatrix))
heatmap_idents$age<-factor(heatmap_idents$age,levels = c("2 years","5 years","8 years","11 years","17 years","adult"),ordered = T)

pheatmap(heatmap_exprmatrix, cellwidth = 0.02, cellheight = 1.5, fontsize = 8,
         method="spearman", #����gene��sample֮�������Եķ�������ѡ"pearson" (default), "kendall", or "spearman"
         scale="row", #Ϊ������scale
         cluster_rows=T, #Ϊ����������
         cluster_cols=F, #Ϊsample������
         color = colorRampPalette(c("steelblue", "white", "darkorange2"))(41), #�Զ�����ɫ
         breaks=breaks,#������ɫ��Χ
         show_colnames=F, #��ʾ��Ʒ���ƾ͸�ΪT
         show_rownames =F, #��ʾ�������͸�ΪT
         annotation_col = heatmap_idents,
         #�������ʾ��״�ṹ����ɾ����������ǰ���#
         treeheight_row = "0",treeheight_col = "0",#������
         #filename = "��̴���Ӧ.pdf", #ֱ�ӱ��浽pdf�ļ�
         border_color = "NA") #����ÿ��С���ӻ��߿�������߿򣬿��԰�NA��Ϊ����Ҫ����ɫ