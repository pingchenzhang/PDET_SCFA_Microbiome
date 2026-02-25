#'
#'
#'
ancom_BC <- function(
  df, sgc, opre=NULL
){
  # start ----
  # https://bioconductor.riken.jp/packages/3.13/bioc/html/ANCOMBC.html
  library(ANCOMBC122)
  library(phyloseq)
  # to phyloseq 
  tibble::enframe(sgc$sg, value = 'Group') -> df_group
  phyloseq(
    otu_table(df[,df_group$name], taxa_are_rows = T),
    sample_data(tibble::column_to_rownames(df_group,'name'))
  ) -> datphy
  # run anocm BC
  ANCOMBC122::ancombc(datphy, formula = "Group", group = "Group",
                   neg_lb = T, global = TRUE) -> outa
  if(length(unique(df_group$Group))==2){
    purrr::imap_dfr(outa$res,function(x, i){
      tibble::rownames_to_column(x) %>%
        dplyr::mutate(.i = i)
    }) %>% tidyr::spread(.i, colnames(.)[2]) %>%
      list(Result = .) -> ot
  }else{
    list(
      res_global =   tibble::rownames_to_column(outa$res_global)
    ) %>% append(
      purrr::map(outa$res, ~{tibble::rownames_to_column(.x)})
    ) -> ot
  }
  ot$feature_table = tibble::rownames_to_column(as.data.frame(outa$feature_table))
  ot$meta_tab <- df_group
  ## 
  if(!is.null(opre)){
    dir.create(dirname(opre),showWarnings = F, recursive = T)
    saveRDS(ot, paste0(opre, '.rds'))
    writexl::write_xlsx(ot, paste0(opre, '.xlsx'))
    cli::cli_alert_success('{opre}.xlsx')
  }else{ return(ot) }
  
}

