
#' sample group color list
#' 
sgc <- function(mapdf, color = NULL, shapes = NULL, orderbyGroup = F){
  # mapdf <- eg.mapw; color = NULL; nn=3; orderbyGroup = F
  gpns   <- colnames(mapdf)[-1]
  gpuniq <- na.omit(unique(unlist(mapdf[,-1])))
  if(is.null(names(color))){
    cc    <- c(color, ggsci::pal_igv()(51), rainbow(100))
    color <- setNames(cc[1:length(gpuniq)], gpuniq)
  }
  ## group
  data.frame(
    gpn   = names(color),
    gpnod = seq_along(color),
    stringsAsFactors = F) -> godf
  ## map group
  purrr::map(gpns,function(nn){
    dplyr::select(mapdf, .id=1, .gp=nn) %>%
      dplyr::filter(!is.na(.gp), .gp!='', .gp!=' ') %>%
      dplyr::left_join(godf, by=c('.gp'='gpn'))  -> gpdf
    if(orderbyGroup) dplyr::arrange(gpdf, gpnod) -> gpdf  # order
    list(sg=setNames(gpdf$.gp,gpdf$.id), gc=color[unique(gpdf$.gp)]) -> ol
    if(!is.null(shapes)) ol$gs <- shapes[unique(gpdf$.gp)]
    ol
  }) %>% purrr::set_names(gpns)
}

