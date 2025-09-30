#'
#'
plot_scatter_correlation <- function(
  data,  
  x_var, 
  y_var, 
  method = "spearman", #  ("pearson", "kendall", "spearman")
  title = "",      
  verbose = TRUE,  
  opre = NULL,   
  pw = 3.5,     
  ph = 3    
) {
  pacman::p_load(tidyverse, ggpubr, cli, writexl)

  if (!method %in% c("pearson", "kendall", "spearman")) {
    cli::cli_warn(" '{method}' 'pearson'。")
    method <- "pearson"
  }
  
  # ----
  if (verbose) cli::cli_alert_info("构建 ggplot 图表...")
  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = !!sym(x_var), y = !!sym(y_var))) +
    ggplot2::geom_point(alpha = 0.8) # 
  
  # lm ----
  p <- p + ggplot2::geom_smooth(
    method = "lm", se = FALSE,
    color = "blue", formula = y ~ x) #
  
  # ----
  #  ggpubr::stat_cor 
  p <- p + ggpubr::stat_cor(
    method = method,
    ggplot2::aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
    label.x.npc = "left", 
    label.y.npc = "top", 
    hjust = 0             
  )
  
  #  ----
  p <- p + ggplot2::ggtitle(title) +
    ggplot2::xlab(x_var) +
    ggplot2::ylab(y_var) +
    ggplot2::theme_classic() 
  
  # ----
  cor_test_result <- stats::cor.test(data[[x_var]], data[[y_var]], method = method)
  result_data <- tibble::tibble(
    x_variable = x_var,
    y_variable = y_var,
    correlation_method = method,
    correlation_coefficient = cor_test_result$estimate,
    p_value = cor_test_result$p.value
  )
  
  ol <- list(
    dat = result_data,
    plot = p
  )
  
  #  ----
  if (!is.null(opre)) {
    opre <- as.character(opre)
    dir.create(dirname(opre), showWarnings = F, recursive = T)
    # 
    ggplot2::ggsave(filename = paste0(opre, ".pdf"), plot = p, width = pw, height = ph)
    cli::cli_alert_success("to: {opre}.pdf")
    # 
    writexl::write_xlsx(list(correlation_results = ol$dat), path = paste0(opre, ".xlsx"))
    cli::cli_alert_success("{opre}.xlsx")
  }
  
  return(ol)
}
