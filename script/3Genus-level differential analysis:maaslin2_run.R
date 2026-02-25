#'
#' 
maaslin2_run <- function(data, sgc=NULL, meta=NULL, odir){
  library(Maaslin2)
  # ?Maaslin2
  # https://github.com/biobakery/Maaslin2?tab=readme-ov-file#options
  # run MaAsLin2
  dir.create(odir, recursive = T)
  if(!is.null(sgc)){
    tibble::enframe(sgc$sg, 'name', 'Group') %>%
      tibble::column_to_rownames('name') -> meta
  }
  tibble::rownames_to_column(as.data.frame(t(data[,rownames(meta)]))) %>%
    write_tsv(file.path(odir,'input_data.tsv'))
  tibble::rownames_to_column(as.data.frame(meta)) %>%
    write_tsv(file.path(odir,'input_metadata.tsv'))
  fit_data = Maaslin2(
    input_data     = file.path(odir,'input_data.tsv'), 
    input_metadata = file.path(odir,'input_metadata.tsv'), 
    output = odir, fixed_effects = colnames(meta))
  # file.remove(paste0(finput,'data.tsv'), paste0(finput,'meta.tsv'))
  
  list(
    result        =read_tsv(glue('{odir}/all_results.tsv'), show_col_types = F),
    input_data    =read_tsv(glue('{odir}/input_data.tsv'), show_col_types = F),
    input_metadata=read_tsv(glue('{odir}/input_metadata.tsv'), show_col_types = F)
  ) %>%
    writexl::write_xlsx(glue('{odir}.MaAsLin2.xlsx'))
  cli::cli_alert_success('{odir}.MaAsLin2.xlsx')
}

