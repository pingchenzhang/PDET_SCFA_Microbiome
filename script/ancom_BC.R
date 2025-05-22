#'
#'
#'
ancom_BC <- function(
  df, sgc, opre=NULL
){
  #' https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
  #' 带有偏差校正的微生物组成分分析 (ANCOM-BC) （Lin 和 Peddada 2020）
  #' 是一种对微生物绝对丰度进行差异丰度 (DA) 分析的方法。
  #'  R_LIBS_USER='-' singularity shell  ~/disk10/sif/sandbox/tidyverse4.3.2/
  if(F){
    # test ----
    pacman::p_load(tidyverse, ANCOMBC)
    library(phyloseq)
    library(microbiome)
    library(tidyverse)
    data(GlobalPatterns)
    
    # Aggregate to phylum level
    phylum_data = aggregate_taxa(GlobalPatterns, "Phylum")
    # The taxonomy table
    tax_mat = as(tax_table(phylum_data), "matrix")
    
    # Run ancombc function
    out = ancombc(phyloseq = phylum_data, formula = "SampleType",
                  p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                  group = "SampleType", struc_zero = TRUE, neg_lb = FALSE,
                  tol = 1e-5, max_iter = 100, conserve = TRUE,
                  alpha = 0.05, global = TRUE)
    
    res = out$res
    res_global = out$res_global
    pkgload::load_all('~/script/rbin/ampv1/',quiet = T)
    sample_data()
    # test ampv1 eg data ----
    dplyr::left_join(enframe(eg.sgc$sg,value='gp3'),
                     enframe(eg.sgc2$sg,value = 'gp2'), by='name') %>%
      tibble::column_to_rownames('name') %>%
      phyloseq::sample_data() -> sap_dat
    phyloseq(
      otu_table(round(eg.phylum*30000, 0), taxa_are_rows = T),
      sap_dat
    ) -> datp
    # Run ancombc function
    out = ancombc(phyloseq = datp, formula = "gp3", group = "gp3",
                  p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                  struc_zero = TRUE, neg_lb = FALSE,
                  tol = 1e-5, max_iter = 100, conserve = TRUE,
                  alpha = 0.05, global = TRUE)
    out$res_global
    
    pkgload::load_all('/mnt/script/rbin/ampv1/',quiet = T)
    df=round(eg.phylum*30000, 0); sgc=eg.sgc2
  }
  # start ----
  # https://bioconductor.riken.jp/packages/3.13/bioc/html/ANCOMBC.html
  library(ANCOMBC122)
  library(phyloseq)
  # 转化为 phyloseq 对象
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

