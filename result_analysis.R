##### Patient stratification
options(stringsAsFactors = F)

library(NbClust)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(maftools)
#print(cancer_type)
cancer_type<-"SKCM"
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
print(cancer_type)
#samples <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))


#length(samples_f2)
cmtScores <- read.csv(paste0("./data/python_related/result/community/combine_",cancer_type,"_score_profile_test.csv"), check.names = F, header = F)
#dim(cmtScores)
#View(cmtScores)
#exp_intgr <- readRDS("./data/tcga_data_processed/SKCM_exp_intgr_all01.RData")
#mty_intgr <- readRDS("./data/tcga_data_processed/SKCM_mty_intgr_all01.RData")
#snv_intgr <- readRDS("./data/tcga_data_processed/SKCM_snv_intgr_test06.RData")
#cnv_intgr <- readRDS("./data/tcga_data_processed/SKCM_cnv_intgr.RData")

#clinicalInfo <- readRDS("./data/tcga_data_processed/SKCM_clinical_info_test06.RData")
#write.csv(as.data.frame(clinicalInfo),"./data/tcga_data_processed/SKCM_clinical_info_test06.csv")
# 
# exp_intgr<-
# mty_intgr<-
# snv_intgr<-
# cnv_intgr<-

#therapy <- readRDS(paste0("./data/tcga_data/",cancer_type,"_therapy_test06.RData"))
#radiation <- readRDS(paste0("./data/tcga_data/",cancer_type,"_radiation_test06.RData"))
#melanet_cmt <- readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

#View(therapy)

#print(melanet_cmt)
#dim(samples_f2)
#print(samples_f2)
#View(clinicalInfo)

new_clinicalInfo<-readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test01.RData"))
clinicalInfo_tmp<-new_clinicalInfo

# View(clinicalInfo_tmp)
# #View(new_clinicalInfo)
# row.names(cmtScores) <- samples_f2
# colnames(cmtScores) <- paste0("cmt", 1:ncol(cmtScores))
# print(cmtScores)
# View(cmtScores)
#saveRDS(cmtScores, paste0("./data/",cancer_type,"_community_scores.RData"))

#print(cmtScores)
### Determine the best number of clusters
nc <- NbClust(scale(cmtScores), distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")

best_clusters <- nc$Best.nc[1, ]  # 获取所有指数的推荐集群数量

# 计算每个集群数量的出现频率
cluster_counts <- table(best_clusters)
View(cluster_counts)
freq_cluster<-names(which.max(cluster_counts))
samples_name<-readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))
####################################################################################

if (!require(kernlab)) install.packages("kernlab")
library(kernlab)  # 加载kernlab包
# 可视化（如果数据是二维或三维的话）
library(ggplot2)

if (!require("Rtsne")) install.packages("Rtsne")
library(Rtsne)


n_rows<-nrow(cmtScores)
normalized_data<-scale(cmtScores)


rownames(cmtScores)<-samples_name
#View(cmtScores)
rownames(normalized_data)<-samples_name
#View(normalized_data)
# set.seed(42)  # 为了可重复性设置随机种子
# 
# tsne_results <- Rtsne(normalized_data, dims = 2, perplexity = ceiling(sqrt(n_rows)), check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
# 
# 
# tsne_data <- as.data.frame(tsne_results$Y)
# 
# # 使用 ggplot2 进行可视化
# 
# ggplot(tsne_data, aes(x = V1, y = V2)) +
#   geom_point(alpha = 0.6) +  # 点的透明度
#   theme_minimal() +
#   labs(title = "t-SNE Visualization of cmtScores")


#-------------------------------------------------------------------------------------------
tsne_results01 <- Rtsne(normalized_data, dims = 3, perplexity =29 , check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
tsne_data01 <- as.data.frame(tsne_results01$Y)
tsne_data_knn<-as.data.frame(tsne_results01$Y)

#View(tsne_data_knn)
#tsne_dbscan<-as.data.frame(tsne_results01$Y)
# ggplot(tsne_data01, aes(x = V1, y = V2,z=v3)) +
#   geom_point(alpha = 0.6) +  # 点的透明度
#   theme_minimal() +
#   labs(title = "t-SNE Visualization of cmtScores")

#nrow(cmtScores)

#row_names <- rownames(cmtScores)
#View(cmtScores)
#rownames(tsne_data01) <- clinicalInfo_tmp$sample

rownames(tsne_data_knn)<-samples_name
# 
# plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers',
#         marker = list(size = 5, opacity = 0.6)) %>%
#   layout(title = "3D t-SNE Visualization of cmtScores",
#          scene = list(xaxis = list(title = 'Component 1'),
#                       yaxis = list(title = 'Component 2'),
#                       zaxis = list(title = 'Component 3')))
# 
# #print(tsne_data01)
# #sc_result <- specc(scale(cmtScores), centers = 3)
# 
# 
# # 计算欧式距离矩阵
# distance_matrix <- as.matrix(dist(tsne_data01))
# 
# # 转换为相似性矩阵，sigma 是高斯核的宽度参数
# sigma <- mean(distance_matrix)  # 示例，实际可能需要调整
# similarity_matrix <- exp(-distance_matrix^2 / (2 * sigma^2))
# k <- 4
# spectral_clust <- specc(similarity_matrix, centers = k)
# clusters <- spectral_clust@.Data
# 
# #rownames(final_data) <- row_names
# tsne_data01$SampleID <- row_names
# print(row_names)
# tsne_data01$cluster <- as.factor(clusters)  # 添加聚类结果到数据框
# 
# # plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3, color = ~cluster, type = 'scatter3d', mode = 'markers',
# #         marker = list(size = 5, opacity = 0.6)) %>%
# #   layout(title = "3D t-SNE with Spectral Clustering",
# #          scene = list(xaxis = list(title = 'Component 1'),
# #                       yaxis = list(title = 'Component 2'),
# #                       zaxis = list(title = 'Component 3')))
# # 
# # print(tsne_data01)
# #---------------------------------------------
# 
# plot_3d <- plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3,
#                    text = ~SampleID,  # 显示行名称
#                    type = 'scatter3d', mode = 'markers',
#                    marker = list(size = 5, opacity = 0.6),
#                    color = ~cluster,  # 根据聚类结果上色
#                    hoverinfo = 'text+name')  # 显示文本和簇名称
# 
# plot_3d <- plot_3d %>%
#   layout(title = "3D Visualization of Spectral Clustering",
#          scene = list(xaxis = list(title = 'Component 1'),
#                       yaxis = list(title = 'Component 2'),
#                       zaxis = list(title = 'Component 3')))
# 
# # 显示图形
# plot_3d
# 
# 
# #View(tsne_data01)
##############################k-means#####################################

#--------------------------------------------------------
# k1 <- 4  # 设定聚类数量
# kmeans_result <- kmeans(tsne_data_knn, centers = k1)
# 
# # 
# if (!require("scatterplot3d")) install.packages("scatterplot3d")
#  library(scatterplot3d)
# # 
# # # 查看聚类结果
# print(kmeans_result$cluster)
# 
# # 可视化
# plot(tsne_data_knn, col = kmeans_result$cluster, pch = 19, cex = 1.5)
# points(kmeans_result$centers, col = 1:k, pch = 8, cex = 2)
# #---------------------------------------------------------------------------------
# s3d <- scatterplot3d(tsne_data$V1, tsne_data$V2, tsne_data$V3, pch=19, cex.symbols=1.5, color=kmeans_result$cluster)
# 
# # 标记聚类中心
# s3d$points3d(kmeans_result$centers[,1], kmeans_result$centers[,2], kmeans_result$centers[,3], col=1:k, pch=8, cex=2)
# 
# 
# 
# if (!require("rgl")) install.packages("rgl")
# library(rgl)
# 
# plot3d(tsne_data$V1, tsne_data$V2, tsne_data$V3, col = rainbow(3)[kmeans_result$cluster], size = 5)
# 
# # 添加聚类中心点
# centers <- kmeans_result$centers
# spheres3d(centers[,1], centers[,2], centers[,3], radius = 0.1, color = "black", alpha = 0.8)
# # 将三维视图转换为交互式网页组件
# if (!require("htmlwidgets")) install.packages("htmlwidgets")
# library(htmlwidgets)
# rglwidget()


#-------------------------------------------------------------
# 创建一个简单的三维散点图
# plot3d(x = tsne_data_knn$V1, y = tsne_data_knn$V2, z = tsne_data_knn$V3,
#        col = kmeans_result$cluster, size = 3)
# 
# # 添加 Kmeans 聚类中心
# points3d(kmeans_result$centers[,1], kmeans_result$centers[,2], kmeans_result$centers[,3],
#          col = "red", size = 10)
# 
# # 添加坐标轴标签
# #xlabel3d("Component 1")
# #ylabel3d("Component 2")
# #zlabel3d("Component 3")
# 
# # 添加标题
# #title3d("3D Visualization of Kmeans Clustering")
# 
# # 添加图例
# legend3d("topright", legend = unique(kmeans_result$cluster), col = 1:length(unique(kmeans_result$cluster)), pch = 16, cex = 1.5)
# 
# # 旋转和交互
# rglwidget()

#-----------------------------------------------------




kmeans_result <- kmeans(tsne_data_knn, centers = freq_cluster)  # 假定我们有 4 个聚类中心

# 将聚类结果添加到数据框中
tsne_data_knn$cluster <- as.factor(kmeans_result$cluster)

centroids <- kmeans_result$centers
centroids_df <- as.data.frame(centroids)
centroids_df$cluster <- as.factor(1:nrow(centroids))  # 为每个中心点创建一个聚类标签

View(tsne_data_knn)
# 使用 plot_ly 创建 3D 散点图来展示聚类结果
plot_kmeans_3d <- plot_ly(data = tsne_data_knn, x = ~V1, y = ~V2, z = ~V3,
                          type = 'scatter3d', mode = 'markers',
                          marker = list(size = 5, opacity = 0.6),
                          color = ~cluster,  # 根据聚类结果上色
                          hoverinfo = 'text+name')  # 显示行名称和聚类

# 添加聚类中心到图形
plot_kmeans_3d <- plot_kmeans_3d %>%
  add_trace(data = centroids_df, x = ~V1, y = ~V2, z = ~V3,
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 6, color = 'black', symbol = 'star'),  # 更大更显眼的中心点
            name = 'Centroids') %>%
  layout(title = "3D Visualization of k-Means Clustering with Centroids",
         scene = list(xaxis = list(title = 'Component 1'),
                      yaxis = list(title = 'Component 2'),
                      zaxis = list(title = 'Component 3')))

# 显示图形
plot_kmeans_3d


View(tsne_data_knn)

#-----------------------------------------------
# install.packages("dbscan")
# library(dbscan)
# 
# library(rgl)
# 
# # 进行 DBSCAN 聚类
# dbscan_result <- dbscan(tsne_dbscan, eps = 5, MinPts = 5)
# 
# # 获取聚类结果
# clusters <- dbscan_result$cluster
# 
# print(clusters)
# 
# # 创建3D散点图
# plot3d(tsne_dbscan, col = clusters, size = 5)
# 
# # 添加标签
# text3d(tsne_dbscan, texts = rownames(tsne_dbscan), adj = c(-0.5,0), cex = 0.7)
# 
# # 设置3D场景
# rglwidget()  # 打开3D视图窗口，可以旋转和缩放
###############################################################
print(cancer)
# 
# best_clusters <- nc$Best.nc[1, ]  # 获取所有指数的推荐集群数量
# 
# # 计算每个集群数量的出现频率
# cluster_counts <- table(best_clusters)

# View(cluster_counts)
# # 找到频率最高的集群数量
# most_frequent_cluster <- as.numeric(names(which.max(cluster_counts)))
# 
# # 打印最频繁推荐的集群数量
# print(most_frequent_cluster)
# 
# pdf(paste0("./figure/",cancer_type,"_best_number_of_clusters01.pdf", width = 7, height = 7))
# ggplot(data.frame(cluster = factor(nc$Best.nc[1,])), aes(x = cluster)) + 
#   geom_bar(stat = "count", fill = "#C1BFBF") + 
#   labs(x = "Number of clusters", y = "Number of criteria", title = "Number of clusters chosen by 26 criteria") + 
#   theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) + 
#   scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14))
# dev.off()

#View(clinicalInfo)
### Community scores of clustered patients
tumorType <- clinicalInfo_tmp[,"shortLetterCode"]

View(clinicalInfo_tmp)

#tumorStage <- clinicalInfo_tmp[,"tumor_stage"]
clinicalInfo_col<-as.data.frame(colnames(clinicalInfo))
#write.csv(clinicalInfo_col,'./data/clinicalInfo_col.csv')


# tumorStage[-grep("^stage", tumorStage)] <- NA
# tumorStage <- gsub("^stage ", "", tumorStage)
# tumorStage <- gsub("[a-c]$", "", tumorStage)
# 
# tumorStage <- clinicalInfo_tmp[,"tumor_stage"]
# tumorStage[-grep("^stage", tumorStage)] <- NA
# tumorStage <- gsub("^stage ", "", tumorStage)
# tumorStage <- gsub("[a-c]$", "", tumorStage)

# therapy <- therapy[3:nrow(therapy),]
# ifTherapy <- substr(samples, 1, 12) %in% therapy$bcr_patient_barcode
# ifTherapy <- ifelse(ifTherapy, "Yes", "No")
# 
# radiation <- radiation[3:nrow(radiation),]
# ifRadiation <- substr(samples, 1, 12) %in% radiation$bcr_patient_barcode
# ifRadiation <- ifelse(ifRadiation, "Yes", "No")
# 
# therapy_type <- sapply(substr(samples, 1, 12), function(x){paste(sort(unique(therapy$pharmaceutical_therapy_type[therapy$bcr_patient_barcode == x])),collapse = ";")})
# 
# tr_df <- data.frame(sample = samples, theray = ifTherapy, radiation = ifRadiation, therapy_type = therapy_type)
# saveRDS(tr_df, paste0("./data/",cancer_type,"_therapy_radiation_df.RData"))
# 
# tumorType_col_fun <- c("TM" = "#CC79A7", "TP" = "#0072B2")
# tumorStage_col_fun <- c("0" = "#FAFCC2", "i" = "#FFEFA0", "ii" = "#FFD57E", "iii" = "#FCA652", "iv" = "#AC4B1C")
# ifTherapy_col_fun <- c("Yes" = "red", "No" = "gray")
# ifRadiation_col_fun <- c("Yes" = "red", "No" = "gray")
# #cancer_type<-"COAD"
# #print(most_frequent_cluster)
# print(cancer_type)
# #topAnno <- HeatmapAnnotation(Therapy = ifTherapy, Radiation = ifRadiation, `Tumor type` = tumorType,  col = list(Therapy = ifTherapy_col_fun, Radiation = ifRadiation_col_fun, `Tumor type` = tumorType_col_fun, `Tumor stage` = tumorStage_col_fun), border = T, show_annotation_name = T)
# ht = Heatmap(t(scale(cmtScores)), 
#              name = "Community score", 
#              show_column_names = F,
#              # top_annotation = topAnno,
#              clustering_distance_columns = "euclidean",
#              #clustering_distance_columns = "spearman",
#              clustering_method_columns = "complete",
#              column_split = most_frequent_cluster,
#              #column_split = 2,
#              column_title = "%s",
# )#`Tumor stage` = tumorStage, removed
# 
# pdf(paste0("./figure/",cancer_type,"_heatmap_cmtScores01_05.pdf", width = 7, height = 7))
# draw(ht, merge_legends = TRUE)
# dev.off()
# 
# 
# ht = draw(ht)
# rowOrder <- row_order(ht)
# colOrder <- column_order(ht)
# print(rowOrder)
# print(colOrder)
# View(as.data.frame(rowOrder))
# View(as.data.frame(colOrder))

#View(ht)
#View(colOrder)
#View(rowOrder)

#View(ht)
##############################################################################################
# samplePartition <- data.frame(cluster = rep(1:length(colOrder), lengths(colOrder)), sampleID = unlist(colOrder))
# samplePartition <- samplePartition[order(samplePartition$sampleID),]
# saveRDS(samplePartition, paste0("./data/",cancer_type,"_sample_partition01.RData"))
#View(data.frame(cluster = rep(1:length(colOrder), lengths(colOrder)), sampleID = unlist(colOrder)))
#View(samplePartition)
#nrow(samplePartition)
#View(samplePartition)
################################################################################
##
##
tsne_data_knn$SampleID<-rownames(tsne_data_knn)
#View(tsne_data_knn)
#View(survivalInfo)
survivalInfo_tmp_test<-survivalInfo
survivalInfo_tmp_test$PatientID <- rownames(survivalInfo)
#View(as.data.frame(substr((survivalInfo_tmp_test$PatientID),1,15)))
survivalInfo_tmp_test$PatientID<-substr((survivalInfo_tmp_test$PatientID),1,16)
#View(survivalInfo_tmp_test)
#View(tsne_data01)
combined_data <- merge(survivalInfo_tmp_test, tsne_data_knn, by.x = "PatientID", by.y = "SampleID")
combined_data$V1<-NULL
combined_data$V2<-NULL
combined_data$V3<-NULL

#colnames(combined_data)[colnames(data) == "cluster"] <- "hc"
#View(combined_data)
#View(survivalInfo_tmp_test)
###################################Section 5################################################################
##### Survival analysis
options(stringsAsFactors = F)
library(TCGAbiolinks)
#cancer_type<-"BLCA"
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
print(cancer_type)
#setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19")
#clinicalInfo <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_new_clinicalInfo_test06.RData"))
#View(new_clinicalInfo)
#clinicalInfo_tmp<-new_clinicalInfo
#samplePartition <- readRDS(paste0("./data/",cancer_type,"_sample_partition01.RData"))
#mty_colData<-readRDS(paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData"))
#mty_colData<-paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData")
#clinicalInfo <- mty_colData
#View(as.data.frame(clinicalInfo))
#row.names(clinicalInfo) <- clinicalInfo$sample
#clinicalInfo <- clinicalInfo[samples,]


#saveRDS(clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
#write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
#dim(samplePartition)
#View(clinicalInfo)
#dim(clinicalInfo)
#View(samplePartition)
#View(samplePartition)
#dim(samplePartition)
#View(clinicalInfo_tmp)
#nrow(clinicalInfo_tmp)
survivalInfo <- clinicalInfo_tmp[,c("shortLetterCode", "vital_status","days_to_death","days_to_last_follow_up")]
# nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308
# > nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308

View(survivalInfo)

survivalInfo_tmp <- survivalInfo[1:nrow(samplePartition),]

View(survivalInfo_tmp)

survivalInfo_tmp$hc <- samplePartition$cluster # ???˲??ξ???????




#View(combined_data)


#View(samplePartition)
#View(survivalInfo_tmp)

# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage 0")] <- 0
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage i")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ia")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ib")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ic")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ii")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iia")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iib")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iic")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iii")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiia")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiib")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiic")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iv")] <- 4
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "not reported")] <- "NA"
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "i/ii nos")] <- "NA"
# survivalInfo$tumor_stage[which(is.na(survivalInfo$tumor_stage))] <- "NA"
#print(cancer_type)
# TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "hc", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update.pdf"), conf.int = F, width = 7, height = 7)
# 
# TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update.pdf"), conf.int = F, width = 7, height = 7)
# 
# #TCGAanalyze_survival(survivalInfo, clusterCol = "tumor_stage", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage.pdf"), conf.int = F, width = 7, height = 7)
# #print(cancer_type)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TCGAanalyze_survival(combined_data, clusterCol = "cluster", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565","#FB9A99", 
                                                                      "#CAB2D6", "#FF7F00", "#C2C2F0", "#B0E2FF"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update.pdf"), conf.int = F, width = 7, height = 7)

TCGAanalyze_survival(combined_data, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update.pdf"), conf.int = F, width = 7, height = 7)

