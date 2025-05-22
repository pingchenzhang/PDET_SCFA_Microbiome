
maaslin2_run <- function(data, sgc=NULL, meta=NULL, odir){
  library(Maaslin2) # ?Maaslin2
  #' https://github.com/biobakery/Maaslin2?tab=readme-ov-file#options
  if(F){
    # 示例数据来自 HMP2 可从 https://ibdmdb.org/ 下载。
    (finput_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2"))
    (finput_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2"))
    # •HMP2_taxonomy.tsv 是以制表符分隔的物种丰度矩阵，列为物种，行为样本。该文件只包含 species 水平的丰度信息。
    # •HMP2_metadata.tsv 是以制表符分隔的元数据矩阵。
    # 1. 读入输入文件
    input_data = read.table(file = finput_data, header = TRUE, sep = "\t",
                            row.names = 1, stringsAsFactors = FALSE)
    input_data[1:5, 1:5]
    input_metadata = read.table(file = finput_metadata, header = TRUE, sep = "\t",
                                row.names = 1, stringsAsFactors = FALSE)
    input_metadata[1:5, ]    
    data=eg.phylum; meta=eg.env; odir='demo_output'
    sgc = eg.sgc
  }
  # 2. 运行 MaAsLin2
  # 下面将基于示例数据运行 MaAsLin2，使用多元回归模型来分析微生物丰度
  # 与 IBD 诊断和菌群失调评分（fixed_effects = c("diagnosis", "dysbiosis")）之间的关系。
  # 结果将保存在工作目录下的 demo_output 文件夹中（output = "demo_output"）
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

