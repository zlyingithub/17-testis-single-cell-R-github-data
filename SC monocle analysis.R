#提取文件做monocle分析（sertolicells）
#基因表达矩阵
allsample_combined_cca_sertolicells_expr<-as.matrix(allsample_combined_cca_sertolicells[["RNA"]]@counts)
#分组信息
allsample_combined_cca_sertolicells_idents<-allsample_combined_cca_sertolicells@meta.data
#基因信息
allsample_combined_cca_sertolicells_features <-as.matrix(allsample_combined_cca_sertolicells_expr[,1])
allsample_combined_cca_sertolicells_features[,1] <-row.names(allsample_combined_cca_sertolicells_features)
colnames(allsample_combined_cca_sertolicells_features)[1] <-c("gene_short_name")
allsample_combined_cca_sertolicells_features<-data.frame(allsample_combined_cca_sertolicells_features)
#创建monocle文件
allsample_combined_cca_sertolicells_idents <- new("AnnotatedDataFrame", data = allsample_combined_cca_sertolicells_idents)
allsample_combined_cca_sertolicells_features <- new("AnnotatedDataFrame", data = allsample_combined_cca_sertolicells_features)
allsample_combined_cca_sertolicells_monocle <- newCellDataSet(as.matrix(allsample_combined_cca_sertolicells_expr),phenoData = allsample_combined_cca_sertolicells_idents, featureData = allsample_combined_cca_sertolicells_features, expressionFamily=negbinomial.size())
#Estimate size factors and dispersions
allsample_combined_cca_sertolicells_monocle <- estimateSizeFactors(allsample_combined_cca_sertolicells_monocle)
allsample_combined_cca_sertolicells_monocle <- estimateDispersions(allsample_combined_cca_sertolicells_monocle)
#质量控制
allsample_combined_cca_sertolicells_monocle <- detectGenes(allsample_combined_cca_sertolicells_monocle, min_expr = 0.1)
print(head(fData(allsample_combined_cca_sertolicells_monocle)))
expressed_genes <- row.names(subset(fData(allsample_combined_cca_sertolicells_monocle), num_cells_expressed >= 10))
pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs <- Matrix::colSums(exprs(allsample_combined_cca_sertolicells_monocle))
allsample_combined_cca_sertolicells_monocle <- allsample_combined_cca_sertolicells_monocle[,pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs)) +
                     2*sd(log10(pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs)) -
                     2*sd(log10(pData(allsample_combined_cca_sertolicells_monocle)$Total_mRNAs)))
qplot(Total_mRNAs, data = pData(allsample_combined_cca_sertolicells_monocle), color = cluster, geom = "density") +geom_vline(xintercept = lower_bound) +geom_vline(xintercept = upper_bound)
allsample_combined_cca_sertolicells_monocle <- detectGenes(allsample_combined_cca_sertolicells_monocle, min_expr = 0.1)
#选差异基因
disp_table <- dispersionTable(allsample_combined_cca_sertolicells_monocle)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.001)
allsample_combined_cca_sertolicells_monocle <- setOrderingFilter(allsample_combined_cca_sertolicells_monocle, unsup_clustering_genes$gene_id)
plot_ordering_genes(allsample_combined_cca_sertolicells_monocle)
plot_pc_variance_explained(allsample_combined_cca_sertolicells_monocle, return_all = F)
diff_test_res <- differentialGeneTest(allsample_combined_cca_sertolicells_monocle,fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.001))
allsample_combined_cca_sertolicells_monocle <- setOrderingFilter(allsample_combined_cca_sertolicells_monocle, ordering_genes)
plot_ordering_genes(allsample_combined_cca_sertolicells_monocle)
#降维,排序，作图
allsample_combined_cca_sertolicells_monocle <- reduceDimension(allsample_combined_cca_sertolicells_monocle, max_components = 2, reduction_method = 'DDRTree',residualModelFormulaStr="~tech")
allsample_combined_cca_sertolicells_monocle <- orderCells(allsample_combined_cca_sertolicells_monocle)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "state",cell_size = 1,cell_link_size = 1.5)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "sample",cell_size = 1,cell_link_size = 1.5)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "age",cell_size = 1,cell_link_size = 1.5)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, cell_size = 1,cell_link_size = 1.5)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "cluster",cell_size = 1,cell_link_size = 1.5)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),plot.background  = element_blank())

#基因表达时序图
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, markers = "EGR3",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
#只做各个群的点图
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "age") + facet_wrap(~age, nrow = 2)
plot_cell_trajectory(allsample_combined_cca_sertolicells_monocle, color_by = "sample") + facet_wrap(~sample, nrow = 2)
#更改顺序名称后将state信息和pseudotime信息加入seurat文件
allsample_combined_cca_sertolicells_monocle_State<-allsample_combined_cca_sertolicells_monocle$State
allsample_combined_cca_sertolicells_monocle_State<-str_replace_all(allsample_combined_cca_sertolicells_monocle_State,c("1"="Stage_a","2"="Stage_b","3"="Stage_c"))
allsample_combined_cca_sertolicells_monocle$State<-allsample_combined_cca_sertolicells_monocle_State
allsample_combined_cca_sertolicells$seurat_sertolicells_7clusters<-Idents(allsample_combined_cca_sertolicells)
allsample_combined_cca_sertolicells_monocle$cluster<-allsample_combined_cca_sertolicells$seurat_sertolicells_5clusters
allsample_combined_cca_sertolicells$state<-allsample_combined_cca_sertolicells_monocle_State
allsample_combined_cca_sertolicells$state<-factor(allsample_combined_cca_sertolicells$state,levels = c("Stage_a","Stage_b","Stage_c","Stage_d"),ordered = T)

#查看随着pseudotime变化的基因
allsample_combined_cca_sertolicells_monocle_pseudotimegenes <- differentialGeneTest(allsample_combined_cca_sertolicells_monocle,fullModelFormulaStr = "~sm.ns(Pseudotime)")
allsample_combined_cca_sertolicells_monocle_sigpseudotimegenes <- subset(allsample_combined_cca_sertolicells_monocle_pseudotimegenes, qval < 0.1)
write.csv(allsample_combined_cca_sertolicells_monocle_sigpseudotimegenes,"随着pseudotime变化的基因.csv")
plot_pseudotime_heatmap(allsample_combined_cca_sertolicells_monocle[c("ISG15","NOC2L","ZFX","SOX9","AMH"),],num_clusters = 3, cores = 1,show_rownames = T)
#按照各群各时期分组选取某个基因做点图
allsample_combined_cca_sertolicells_monocle_sub <- allsample_combined_cca_sertolicells_monocle[c("ISG15","NOC2L","ZFX","SOX9","AMH"),]
plot_genes_jitter(allsample_combined_cca_sertolicells_monocle_sub,grouping = "State", color_by = "age", plot_trend = T) + facet_wrap( ~ feature_label, scales= "free_y")+ stat_summary(fun.y=mean, geom="point", shape=18)

#查看在分叉处的差异基因
allsample_combined_cca_sertolicells_monocle_BEAM_res <- BEAM(allsample_combined_cca_sertolicells_monocle, branch_point = 1, cores = 1)
allsample_combined_cca_sertolicells_monocle_BEAM_res <- allsample_combined_cca_sertolicells_monocle_BEAM_res[order(allsample_combined_cca_sertolicells_monocle_BEAM_res$qval),]
allsample_combined_cca_sertolicells_monocle_BEAM_res <- allsample_combined_cca_sertolicells_monocle_BEAM_res[,c("gene_short_name", "pval", "qval")]
write.csv(allsample_combined_cca_sertolicells_monocle_BEAM_res,file = "sertolicells_BEAM_res分叉处基因.csv")
plot_genes_branched_heatmap(allsample_combined_cca_sertolicells_monocle[row.names(subset(allsample_combined_cca_sertolicells_monocle_BEAM_res, qval < 1e-120)),],branch_point = 1,num_clusters = 6,cores = 1,use_gene_short_name = T,show_rownames = T)
#使用分叉处基因作图
A<-c("HOPX","EGR1","YBX3","FOS","JUNB","FOSB","ID3","MEF2C","JUND","JUN","NR4A1","ZFP36L2","TSC22D1","RORA","TMF1","NFIX")
B<-c("ZFP36L1","PIAS2","ZNF677","ATF3","ARID2","SMARCA1","ZNF91","PMS1","LRRFIP2","HMGB1","KLF2","BBX","NR0B1","CDC5L","TBX22","CREB3L4")
C<-c("ISG15","TCEA3", "IFI6","S100A13","TMEM176B", "BEX1")
plot_genes_branched_pseudotime(allsample_combined_cca_sertolicells_monocle[C,],branch_point = 1,color_by = "State",ncol = 3,cell_size=1.5,branch_labels = c("to Stage_b","to Stage_c"))
#左右基因分化热图
plot_genes_branched_heatmap(allsample_combined_cca_sertolicells_monocle[row.names(subset(allsample_combined_cca_sertolicells_monocle_BEAM_res, qval < 1e-5)),],
                            branch_point = 1,
                            branch_labels = c("to Stage_b", "to Stage_c"),
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = F,
                            hmcols = colorRampPalette(c("steelblue", "white","darkorange2", "orangered4"))(100))
