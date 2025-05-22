# 定义绘制散点图并添加相关性信息的函数
plot_scatter_correlation <- function(
  data,           # 输入的数据框
  x_var,          # X轴变量名（字符串）
  y_var,          # Y轴变量名（字符串）
  method = "spearman", # 计算相关性的方法 ("pearson", "kendall", "spearman")
  title = "",       # 图表标题（字符串）
  verbose = TRUE,   # 是否显示进度信息（逻辑值）
  opre = NULL,      # 输出文件路径前缀（字符串，可选）
  pw = 3.5,           # 保存图形的宽度 (英寸)
  ph = 3            # 保存图形的高度 (英寸)
) {
  # 加载所需包 ----
  pacman::p_load(tidyverse, ggpubr, cli, writexl)
  
  if (verbose) cli::cli_alert_info("开始绘制散点图并计算相关性...")
  
  # 检查输入参数 ----
  if (!is.data.frame(data)) {
    cli::cli_abort("输入 'data' 必须是数据框。")
  }
  if (!all(c(x_var, y_var) %in% names(data))) {
    cli::cli_abort("'{x_var}' 和 '{y_var}' 必须是 'data' 中的列名。")
  }
  if (!is.numeric(data[[x_var]]) || !is.numeric(data[[y_var]])) {
    cli::cli_abort("'{x_var}' 和 '{y_var}' 列必须是数值类型。")
  }
  if (!method %in% c("pearson", "kendall", "spearman")) {
    cli::cli_warn("无效的相关性方法 '{method}'。将使用默认的 'pearson'。")
    method <- "pearson"
  }
  
  # 构建基本图表----
  if (verbose) cli::cli_alert_info("构建 ggplot 图表...")
  p <- ggplot2::ggplot(data, ggplot2::aes(x = !!sym(x_var),
                                          y = !!sym(y_var))) +
    ggplot2::geom_point(alpha = 0.8) # 添加散点
  
  # 添加拟合线----
  p <- p + ggplot2::geom_smooth(
    method = "lm", se = FALSE,
    color = "blue", formula = y ~ x) # 添加线性拟合线，不显示置信区间
  
  # 添加相关性系数和p值----
  # 使用 ggpubr::stat_cor 添加相关性信息
  p <- p + ggpubr::stat_cor(
    method = method,
    ggplot2::aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), # 同时显示 R 值和 p 值
    label.x.npc = "left",  # 将标签放在左侧
    label.y.npc = "top",   # 将标签放在顶部
    hjust = 0              # 左对齐标签
  )
  
  # 添加标题和标签 ----
  p <- p + ggplot2::ggtitle(title) +
    ggplot2::xlab(x_var) +
    ggplot2::ylab(y_var) +
    ggplot2::theme_classic() # 使用简洁主题
  
  # 打印图表----
  #print(p)
  if (verbose) cli::cli_alert_success("散点图绘制完成。")
  
  # 关键结果和数据 ----
  # 计算相关性以供返回
  cor_test_result <- stats::cor.test(data[[x_var]], data[[y_var]], method = method)
  result_data <- tibble::tibble(
    x_variable = x_var,
    y_variable = y_var,
    correlation_method = method,
    correlation_coefficient = cor_test_result$estimate,
    p_value = cor_test_result$p.value
  )
  
  ol <- list(
    dat = result_data, # 包含相关性结果的数据框
    plot = p           # ggplot 对象
  )
  
  # 保存结果 ----
  if (!is.null(opre)) {
    opre <- as.character(opre)
    dir.create(dirname(opre), showWarnings = F, recursive = T)
    if (verbose) cli::cli_alert_info("正在保存结果到文件，前缀为: {opre}")
    
    # 保存图表
    ggplot2::ggsave(filename = paste0(opre, ".pdf"), plot = p, width = pw, height = ph)
    cli::cli_alert_success("图表已保存为: {opre}.pdf")
    
    # 保存数据
    writexl::write_xlsx(list(correlation_results = ol$dat), path = paste0(opre, ".xlsx"))
    cli::cli_alert_success("相关性结果已保存为: {opre}.xlsx")
    
    if (verbose) cli::cli_alert_success("结果文件保存成功。")
  }
  
  return(ol)
}
