##### Patient stratification
options(stringsAsFactors = F)

library(NbClust)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(maftools)
#install.packages("umap")
library(umap)
library(cluster)
#if (!require(kernlab)) install.packages("kernlab")
library(kernlab)  # 加载kernlab包
# 可视化（如果数据是二维或三维的话）
library(ggplot2)

#if (!require("Rtsne")) install.packages("Rtsne")
library(Rtsne)
#print(cancer_type)

dataset_TCGA <- c("BLCA", "BRCA", "COAD", "ESCA", "KICH", "KIRC", 
                  "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", 
                  "READ", "SKCM", "STAD", "THCA", "THYM", "UCEC")

cancer_type<-"COAD"
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
# nc <- NbClust(scale(cmtScores), distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")
# 
# #class(scale(cmtScores))
# best_clusters <- nc$Best.nc[1, ]  # 获取所有指数的推荐集群数量
# 
# # 计算每个集群数量的出现频率
# cluster_counts <- table(best_clusters)
# View(cluster_counts)
# freq_cluster<-names(which.max(cluster_counts))
samples_name<-readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))
####################################################################################




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
set.seed(123)

#-------------------------------------------------------------------------------------------
tsne_results01 <- Rtsne(normalized_data, dims = 3, perplexity =30 , check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
tsne_data01 <- as.data.frame(tsne_results01$Y)
tsne_data_knn<-as.data.frame(tsne_results01$Y)





#View(tsne_data_numeric)

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



tsne_data_knn_test<-tsne_data_knn

tsne_data_numeric <- tsne_data_knn[, c("V1", "V2", "V3")]

# 确保所有列都是数值型
tsne_data_numeric[] <- lapply(tsne_data_numeric, as.numeric)

tsne_matrix<-as.matrix(tsne_data_numeric)

#class(tsne_matrix)
#View(tsne_matrix)

#if (!require("NbClust")) install.packages("NbClust", dependencies = TRUE)


# 使用NbClust进行聚类数的确定
# nc_result <- NbClust(tsne_matrix, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")
# 
# 
# #View(nc_result)
# best_clusters_kmeans <- nc_result$Best.nc[1, ]  # 获取所有指数的推荐集群数量
# #View(best_clusters_kmeans)
# # 计算每个集群数量的出现频率
# cluster_counts_kmeans <- table(best_clusters_kmeans)
# #View(cluster_counts_kmeans)
# freq_cluster_means<-names(which.max(cluster_counts_kmeans))
# print(freq_cluster_means)

################################################################################
num_iterations <- 20  # 例如，迭代100次
all_clusters <- integer(num_iterations)  # 预分配存储空间





# 循环执行NbClust
# for (i in 1:num_iterations) {
#   #nc_result <- NbClust(tsne_matrix, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")
#   con_result <- file("nul", open = "wt")
#   sink(con_result, type = "message")
#   result <- capture.output(
#     print(paste0("iteration :",i)) ,
#     #cat("iteration :",i)
#     nc_result <- NbClust(tsne_matrix, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all"),
#     type = "message"
#   )
#   
#   best_clusters_kmeans <- nc_result$Best.nc[1, -1]  # 获取所有指数的推荐集群数量
#   cluster_counts_kmeans <- table(best_clusters_kmeans)  # 计算每个集群数量的出现频率
#   freq_cluster_means <- as.integer(names(which.max(cluster_counts_kmeans)))  # 找到频率最高的集群数量
#   all_clusters[i] <- freq_cluster_means  # 存储每次迭代的结果
#   sink(type = "message")
#   close(con_result)
#   print(all_clusters[i])
# }

nc_result <- NbClust(tsne_matrix, distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")
best_clusters_kmeans <- nc_result$Best.nc[1, -1]  # 获取所有指数的推荐集群数量
cluster_counts_kmeans <- table(best_clusters_kmeans)  # 计算每个集群数量的出现频率
freq_cluster_means <- as.integer(names(which.max(cluster_counts_kmeans)))  # 找到频率最高的集群数量
most_freq_cluster<-freq_cluster_means
#print(cancer_type)
View(cluster_counts_kmeans)
#nc_result<-NULL
# 计算所有迭代中出现最频繁的聚类数
# final_cluster_counts <- table(freq_cluster_means)
# View(final_cluster_counts)
# most_freq_cluster <- as.integer(names(which.max(final_cluster_counts)))
# print(most_freq_cluster)

# sink(type = "message")
# close(con_result)

#all_clusters<-NULL
#final_cluster_counts<-NULL
################################################################################



# wcss <- sapply(1:10, function(k) {
#   kmeans_result_test <- kmeans(tsne_data_numeric, centers = k, nstart = 30)
#   kmeans_result_test$tot.withinss
# })
# 
# plot(1:10, wcss, type = 'b', xlab = 'Number of Clusters', ylab = 'WCSS', main = 'Elbow Method for Optimal K')
# 
# 
# 
# wcss_diff <- diff(wcss)
# 
# # 找到差分最大的点，这通常对应于WCSS下降最快的位置，即肘部
# elbow_point <- which.min(wcss_diff) + 1  # 加1是因为差分向量比原始WCSS向量短1
# 
# # 输出肘部点
# print(paste("The Elbow point is at k =", elbow_point))
# 
# # 绘制WCSS曲线
# plot(1:length(wcss), wcss, type = 'b', xlab = 'Number of Clusters', ylab = 'WCSS', main = 'Elbow Method for Optimal K')
# # 添加一个点以突出显示肘部位置
# points(elbow_point, wcss[elbow_point], col='red', pch=19, cex=1.5)
# 
# print(elbow_point)



#summary(tsne_data_knn_test)
#any(is.na(tsne_data_knn_test)) # 检查是否有NA
#any(is.nan(tsne_data_knn_test)) # 检查是否有NaN
#any(is.infinite(tsne_data_knn_test)) # 检查是否有Inf
#######################################################

# sil_width <- sapply(2:10, function(k) {
#   model <- pam(tsne_data_numeric, k)
#   model$silinfo$avg.width
# })
# 
# # 找到最大平均轮廓系数的索引
# optimal_clusters <- which.max(sil_width)
# 
# # 由于簇数的范围是从2到10，计算得到的索引对应的簇数是索引+1（因为sapply的起始索引是2）
# optimal_clusters <- optimal_clusters + 1
# 
# # 输出最优簇数
# print(paste("The optimal number of clusters is:", optimal_clusters))
# 
# # 可视化平均轮廓宽度与簇数的关系
# plot(2:10, sil_width, type = 'b', xlab = 'Number of Clusters', ylab = 'Average Silhouette Width', main = 'Optimal number of clusters')
# points(optimal_clusters, max(sil_width), col='red', pch=19, cex=1.5)  # 突出显

###########################################################
kmeans_result <- kmeans(tsne_data_knn, centers = most_freq_cluster)  # 假定我们有 4 个聚类中心
#View(as.data.frame(kmeans_result))
# 将聚类结果添加到数据框中

tsne_data_knn$cluster <- as.factor(kmeans_result$cluster)

centroids <- kmeans_result$centers
centroids_df <- as.data.frame(centroids)
centroids_df$cluster <- as.factor(1:nrow(centroids))  # 为每个中心点创建一个聚类标签

#View(tsne_data_knn)
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


#View(tsne_data_knn)

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
#print(cancer)
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

#View(clinicalInfo_tmp)

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
#print(cancer_type)
survivalInfo <- clinicalInfo_tmp[,c("shortLetterCode", "vital_status","days_to_death","days_to_last_follow_up","ajcc_pathologic_stage")]
#View(clinicalInfo_tmp)
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
#survivalInfo <- clinicalInfo_tmp[,c("shortLetterCode", "vital_status","days_to_death","days_to_last_follow_up")]
# nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308
# > nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308

# View(survivalInfo)
# 
# survivalInfo_tmp <- survivalInfo[1:nrow(samplePartition),]
# 
# View(survivalInfo_tmp)
# 
# survivalInfo_tmp$hc <- samplePartition$cluster # ???˲??ξ???????


survivalInfo_tmp_test_record<-survivalInfo_tmp_test

#View(combined_data)

View(survivalInfo_tmp_test_record)
#View(samplePartition)
#View(survivalInfo_tmp)

survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage 0")] <- 0
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage i")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ia")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ib")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ic")] <- 1
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage ii")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iia")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iib")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iic")] <- 2
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iii")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiia")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiib")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iiic")] <- 3
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "stage iv")] <- 4
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "not reported")] <- "NA"
survivalInfo_tmp_test$ajcc_pathologic_stage[which(tolower(survivalInfo_tmp_test$ajcc_pathologic_stage) == "i/ii nos")] <- "NA"
survivalInfo_tmp_test$ajcc_pathologic_stage[which(is.na(survivalInfo_tmp_test$ajcc_pathologic_stage))] <- "NA"





# #0
# indices <- grepl("stage 0", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 0
# 
# #1
# indices <- grepl("stage i", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ia", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ib", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ic", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 1
# 
# #2
# indices <- grepl("stage ii", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iia", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iib", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iic", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 2
# #3
# indices <- grepl("stage iii", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiia", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiib", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiic", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 3
# #4
# indices <- grepl("stage iv", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- 4
# #NA
# indices <- grepl("not reported", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- NA
# 
# indices <- grepl("i/ii nos", survivalInfo_tmp_test$ajcc_pathologic_stage, ignore.case = TRUE)
# survivalInfo_tmp_test$ajcc_pathologic_stage[indices] <- NA
# 
# survivalInfo_tmp_test$ajcc_pathologic_stage[which(is.na(survivalInfo_tmp_test$ajcc_pathologic_stage))] <- "NA"
# 
# #View(survivalInfo_tmp_test)
# 
# 
# #View(combined_data)
# 
# #0
# indices <- grepl("stage 0", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 0
# 
# #1
# indices <- grepl("stage i", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ia", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ib", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 1
# 
# indices <- grepl("stage ic", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 1
# 
# #2
# indices <- grepl("stage ii", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iia", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iib", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 2
# indices <- grepl("stage iic", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 2
# #3
# indices <- grepl("stage iii", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiia", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiib", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 3
# indices <- grepl("stage iiic", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 3
# #4
# indices <- grepl("stage iv", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- 4
# #NA
# indices <- grepl("not reported", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- NA
# 
# indices <- grepl("i/ii nos", combined_data$ajcc_pathologic_stage, ignore.case = TRUE)
# combined_data$ajcc_pathologic_stage[indices] <- NA
# 
# combined_data$ajcc_pathologic_stage[which(is.na(combined_data$ajcc_pathologic_stage))] <- "NA"
# 
# 
# 
# View(survivalInfo_tmp_test)
# TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "hc", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update.pdf"), conf.int = F, width = 7, height = 7)
# 
# TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update.pdf"), conf.int = F, width = 7, height = 7)
# 
# #TCGAanalyze_survival(survivalInfo, clusterCol = "tumor_stage", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage.pdf"), conf.int = F, width = 7, height = 7)
# #print(cancer_type)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TCGAanalyze_survival(combined_data, clusterCol = "cluster", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565","#FB9A99", 
                                                                      "#CAB2D6", "#FF7F00", "#C2C2F0", "#B0E2FF"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update01.pdf"), conf.int = F, width = 7, height = 7)



TCGAanalyze_survival(survivalInfo_tmp_test, clusterCol = "ajcc_pathologic_stage", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage_update03.pdf"), conf.int = F, width = 7, height = 7)

TCGAanalyze_survival(combined_data, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update01.pdf"), conf.int = F, width = 7, height = 7)
