gg_varImpPlot <- function(mrf, main = 'RF') {
  imp_df <- randomForest::importance(mrf, scale = T) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable")
  # randomForest::varImpPlot(mrf)
  # plot
  dplyr::arrange(imp_df, MeanDecreaseAccuracy) %>%
    ggplot(aes(x = MeanDecreaseAccuracy, y = fct_inorder(variable))) +
    geom_point(size = 2.6, shape = 21) +
    theme_bw(base_size = 13) +
    labs(y = NULL) +
    theme(
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, l=0.5)
    ) -> pmda

  dplyr::arrange(imp_df, MeanDecreaseGini) %>%
    ggplot(aes(x = MeanDecreaseGini, y = fct_inorder(variable))) +
    geom_point(size = 2.6, shape = 21) +
    theme_bw(base_size = 13) +
    labs(y = NULL) +
    theme(
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, l=0.5)
    ) -> pgini

  list(
    pmda = pmda,
    pgini = pgini,
    pmerge = pmda + pgini
  )
}
