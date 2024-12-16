library(parallel)
library(ggplot2)
library(dplyr)
library(transport)
library(entropy)

rank_by_relative_entropy <- function(
  adata, 
  control_guide = 'sgNon-targeting', 
  guide_list = NULL, 
  return_fig = FALSE
) {
  # 初始化熵存储列表
  entro <- list()
  
  # 计算 control_guide 的 qk 分布
  qk_guide <- colSums(adata[, control_guide]) / sum(adata[, control_guide])
  
  # 如果 guide_list 未提供，使用所有变量名称
  if (is.null(guide_list)) {
    guide_list <- colnames(adata)
  }
  
  # 计算每个 guide 的 pk 分布并计算 KL divergence
  for (guide in guide_list) {
    pk_guide <- colSums(adata[, guide]) / sum(adata[, guide])
    entro[[guide]] <- KL.plugin(pk_guide, qk_guide)
  }
  
  # 将熵数据转换为数据框，并排序
  df <- data.frame(Guide = names(entro), Shannon = unlist(entro)) %>%
    arrange(Shannon)
  
  # 绘制条形图
  p <- ggplot(df, aes(x = reorder(Guide, Shannon), y = Shannon, fill = Shannon)) +
    geom_bar(stat = "identity") +
    scale_fill_distiller(palette = "RdBu", direction = -1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "KL distance", x = "Guide", y = "Shannon Index")
  
  print(p)
  
  # 返回图形或数据框
  if (return_fig) {
    return(p)
  } else {
    return(df)
  }
}

rank_by_kernel_estimated_distance <- function(
  adata, 
  control_guide = 'sgNon-targeting', 
  guide_list = NULL, 
  n_permutation = 50, 
  n_process = 8, 
  return_fig = FALSE, 
  return_dataframe = TRUE, 
  sort_by_replicate = '_'
) {
  
  # 定义 KDE 计算函数
  compute_kde_statsmodels <- function(data) {
    kde <- density(data)
    return(kde$y)  # 返回密度值
  }
  
  # 定义工作者函数来计算实际和置换的 Wasserstein 距离
  kde_worker <- function(args) {
    guide <- args$guide
    control <- args$control
    spatial_coords <- args$spatial_coords
    data <- args$data
    n_permutation <- args$n_permutation
    
    if (guide == control) return(NULL)
    
    gene_a_expression <- as.numeric(data[, guide])
    gene_b_expression <- as.numeric(data[, control])
    
    # 将空间坐标和表达数据结合
    data_a <- cbind(spatial_coords[,1], spatial_coords[,2], gene_a_expression)
    data_b <- cbind(spatial_coords[,1], spatial_coords[,2], gene_b_expression)
    
    # 计算 KDE
    Za <- compute_kde_statsmodels(data_a)
    Zb <- compute_kde_statsmodels(data_b)
    
    # 计算 Wasserstein 距离
    actual_distance <- wasserstein1d(Za, Zb)
    
    combined_expression <- c(gene_a_expression, gene_b_expression)
    permuted_distances <- numeric(n_permutation)
    
    # 随机打乱数据并进行置换
    set.seed(42)
    for (i in 1:n_permutation) {
      permuted_expression <- sample(combined_expression)
      permuted_a <- permuted_expression[1:length(gene_a_expression)]
      permuted_b <- permuted_expression[(length(gene_a_expression)+1):length(permuted_expression)]
      
      permuted_Za <- compute_kde_statsmodels(cbind(spatial_coords[,1], spatial_coords[,2], permuted_a))
      permuted_Zb <- compute_kde_statsmodels(cbind(spatial_coords[,1], spatial_coords[,2], permuted_b))
      
      permuted_distances[i] <- wasserstein1d(permuted_Za, permuted_Zb)
    }
    
    # 计算 p-value
    p_value <- sum(permuted_distances >= actual_distance) / length(permuted_distances)
    
    return(c(guide, actual_distance, p_value))
  }
  
  # 如果 guide_list 未提供，使用所有变量
  if (is.null(guide_list)) {
    guide_list <- rownames(adata)
  }
  
  # 检查是否同时返回图像和数据框
  if (return_fig && return_dataframe) {
    stop('Error, can only return one of fig or dataframe!')
  }
  
  # 筛选数据
  filtered_data <- adata[rowSums(adata[, guide_list] > 0) > 0, ]
  spatial_coords <- as.matrix(filtered_data$spatial)
  
  # 使用 parallel 包的 mclapply 进行并行处理
  args_list <- lapply(guide_list, function(guide) list(guide = guide, control = control_guide, spatial_coords = spatial_coords, data = filtered_data, n_permutation = n_permutation))
  results <- mclapply(args_list, kde_worker, mc.cores = n_process)
  
  # 过滤结果并转换为数据框
  d_df <- do.call(rbind, results)
  d_df <- as.data.frame(d_df, stringsAsFactors = FALSE)
  colnames(d_df) <- c('guide', 'wd', 'p_value')
  d_df$wd <- as.numeric(d_df$wd)
  d_df$p_value <- as.numeric(d_df$p_value)
  
  d_df <- d_df %>% arrange(wd)
  
  # 可视化
  if (!is.null(sort_by_replicate)) {
    d_df$replicate <- sapply(strsplit(d_df$guide, sort_by_replicate), `[`, 1)
    positions <- seq(1, nrow(d_df)) * 2
    width <- 0.5
    
    p <- ggplot(d_df, aes(x = reorder(guide, -wd), y = wd)) +
      geom_bar(stat = "identity", width = width, fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Guide", y = "WD", title = "Kernel Estimated Distance")
    
    print(p)
    
  } else {
    p <- ggplot(d_df, aes(x = reorder(guide, wd), y = wd, fill = guide)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Guide", y = "WD", title = "Kernel Estimated Distance")
    
    print(p)
  }
  
  if (return_dataframe) {
    return(d_df)
  } else if (return_fig) {
    return(p)
  } else {
    return(NULL)
  }
}

rank_by_proportion_within_threshold <- function(
  adata, 
  thresholds = seq(1, 501, length.out = 21),
  guide_list = NULL,
  control_guide = 'sgNon-targeting',
  method = 'bin',
  kernel_func = 'log',
  kernel_scoring = NULL,
  return_fig = FALSE
) {
  
  proportion <- list()
  
  # 如果 guide_list 没有提供，使用所有变量名称
  if (is.null(guide_list)) {
    guide_list <- rownames(adata)
  }
  
  # 二分法（bin）或者计数法（count）
  if (method == 'bin') {
    for (threshold in thresholds) {
      proportion[[as.character(threshold)]] <- c()
      for (guide in guide_list) {
        if (guide == control_guide) next
        
        guide_data <- adata[, adata[guide, ] > 0]
        ntc_data <- adata[, adata[control_guide, ] > 0]
        
        # 计算距离矩阵，并检查 guide 数据是否在 threshold 距离内
        in_distance <- sum(rowSums(as.matrix(dist(guide_data@obsm$spatial, ntc_data@obsm$spatial)) < threshold) > 0)
        proportion[[as.character(threshold)]] <- c(proportion[[as.character(threshold)]], in_distance / nrow(guide_data))
      }
    }
    
  } else if (method == 'count') {
    for (threshold in thresholds) {
      proportion[[as.character(threshold)]] <- c()
      for (guide in guide_list) {
        if (guide == control_guide) next
        
        guide_data <- adata[, adata[guide, ] > 0]
        ntc_data <- adata[, adata[control_guide, ] > 0]
        
        # 计算距离并对 guide 数据进行计数
        in_distance <- sum(guide_data[(rowSums(as.matrix(dist(guide_data@obsm$spatial, ntc_data@obsm$spatial)) < threshold)) > 0, guide])
        proportion[[as.character(threshold)]] <- c(proportion[[as.character(threshold)]], in_distance / nrow(guide_data))
      }
    }
    
  } else {
    stop('Error, method must be one of "bin" or "count"!')
  }
  
  # 将 proportion 转换为 data frame
  d_df <- melt(as.data.frame(proportion))
  d_df$variable <- rep(guide_list, each = length(thresholds))
  
  # 根据不同的 kernel 计算评分
  if (!is.null(kernel_scoring)) {
    sc <- apply(d_df, 1, function(x) sum(kernel_scoring * (1 - as.numeric(x))))
  } else {
    if (kernel_func == 'log') {
      ln_scoring <- log2(seq(1, exp2(10), length.out = length(thresholds)))
      sc <- apply(d_df, 1, function(x) sum(ln_scoring * (1 - as.numeric(x))))
    } else if (kernel_func == 'linear') {
      linear_scoring <- seq(0, 10, length.out = length(thresholds))
      sc <- apply(d_df, 1, function(x) sum(linear_scoring * (1 - as.numeric(x))))
    } else if (kernel_func == 'exp') {
      exp_scoring <- 2 ^ seq(1e-32, log2(10), length.out = length(thresholds))
      sc <- apply(d_df, 1, function(x) sum(exp_scoring * (1 - as.numeric(x))))
    }
  }
  
  # 创建 dataframe
  df <- data.frame(dist = sc, guide = guide_list)
  
  # 可视化
  p <- ggplot(df, aes(x = reorder(guide, dist), y = dist)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_log10() +
    labs(x = "Gene", y = "Distance score", title = "Distance Score by Guide")
  
  print(p)
  
  if (return_fig) {
    return(p)
  } else {
    return(NULL)
  }
}
