library(Seurat)
library(ggplot2)

rank_by_chi_square <- function(
  adata,
  cluster_field,
  control_guide = 'sgNon-targeting',
  guide_list = NULL,
  return_fig = FALSE
) {
  if (!(cluster_field %in% colnames(adata@meta.data))) {
    print("Error, cluster field not in Seurat object metadata, please perform clustering first!")
    return(NULL)
  }
  
  if (is.null(guide_list)) {
    guide_list <- c(rownames(adata), control_guide)
  }
  
  c_df <- data.frame(t(adata@assays$RNA@data[guide_list, ]))
  c_df$cluster <- adata@meta.data[[cluster_field]]
  c_df <- aggregate(. ~ cluster, data = c_df, sum)
  rownames(c_df) <- c_df$cluster
  c_df$cluster <- NULL
  c_df <- sweep(c_df, 1, rowSums(c_df), FUN = "/") * sum(c_df[control_guide, ])
  
  chi_dict <- list()
  
  for (guide in guide_list) {
    if (sum(c_df[guide, ]) == 0) {
      chi_dict[[guide]] <- 1
      next
    }
    if (guide == 'sgNon-targeting') next
    contingency_table <- rbind(c_df[control_guide, , drop = FALSE], c_df[guide, , drop = FALSE])
    chi_dict[[guide]] <- chisq.test(contingency_table)$p.value
  }
  
  pdf <- data.frame(Guide = names(chi_dict), `Chi2 p-value` = unlist(chi_dict))
  pdf <- pdf[order(pdf$`Chi2 p-value`), ]
  
  plot <- ggplot(pdf, aes(x = reorder(Guide, `Chi2 p-value`), y = `Chi2 p-value`)) +
    geom_bar(stat = 'identity', fill = scales::brewer_pal(palette = "RdBu")(length(guide_list))) +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Cluster Chi-square Test") +
    xlab("Guide") +
    ylab("Chi2 p-value")
  
  if (return_fig) {
    return(plot)
  } else {
    print(plot)
  }
}

rank_by_cluster_permanova <- function(
  adata, 
  cluster_field, 
  control_guide = 'sgNon-targeting', 
  n_permutation = 999, 
  guide_list = NULL, 
  return_fig = FALSE
) {
  
  if (!(cluster_field %in% colnames(adata@meta.data))) {
    stop("Error, cluster field not in adata obs, please perform clustering first!")
  }
  
  if (is.null(guide_list)) {
    guide_list <- rownames(adata)
  }
  
  # 提取 count 数据和 cluster 信息
  c_df <- data.frame(adata@assays$RNA@counts[guide_list, ], 
                     cluster = adata@meta.data[[cluster_field]])
  
  n_clusters <- length(unique(c_df$cluster))
  p_values <- list()
  
  # 遍历 guide_list 中的每个 guide，排除 control_guide
  for (guide in guide_list) {
    
    if (guide == control_guide) next
    
    # 生成 guide 和 control 的数据框
    g_df <- c_df %>%
      select(all_of(c(guide, control_guide, 'cluster')))
    
    # 计算每个 cluster 中 guide 和 control 的数量
    guide_counts <- table(g_df$cluster, g_df[[guide]])
    control_counts <- table(g_df$cluster, g_df[[control_guide]])
    
    # 合并 guide 和 control 的数据用于距离计算
    data_matrix <- rbind(guide_counts, control_counts)
    
    # 计算欧几里得距离矩阵
    dist_matrix <- dist(data_matrix, method = "euclidean")
    
    # 创建元数据，指定 PERMANOVA 分组
    metadata <- data.frame(group = rep(c('guide', 'control'), each = n_clusters))
    
    # 执行 PERMANOVA 测试
    results <- adonis2(as.dist(dist_matrix) ~ group, data = metadata, permutations = n_permutation)
    
    # 存储 p-value
    p_values[[guide]] <- results$aov.tab$`Pr(>F)`[1]
  }
  
  # 创建 p 值数据框并排序
  p_df <- data.frame(guide = names(p_values), p_value = unlist(p_values))
  p_df <- p_df %>% arrange(desc(p_value))
  
  # 可视化 p 值
  p <- ggplot(p_df, aes(x = reorder(guide, -p_value), y = p_value, fill = guide)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = min(p_df$p_value), linetype = "dashed", color = "black", alpha = 0.4) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Cluster PERMANOVA", x = "Guide", y = "p-value")
  
  print(p)
  
  if (return_fig) {
    return(p)
  } else {
    return(NULL)
  }
}
