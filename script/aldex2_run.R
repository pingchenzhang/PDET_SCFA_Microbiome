#'
#'
aldex2_run <- function(dat, sgc, opre=NULL){
  library(ALDEx2)
  ALDEx2::aldex(
    reads = dat, conditions = sgc$sg[colnames(dat)],
    denom = 'all',effect = TRUE ) -> res
  list(
    tibble::rownames_to_column(res),
    data = tibble::rownames_to_column(dat),
    Group= tibble::enframe(sgc$sg)
  )  -> odat
  
  ##
  list(
    dat = dat, sgc=sgc, res= res, odat = odat
  ) -> ol
  if(is.null(opre)){ return(ol) }else{
    dir.create(dirname(opre), showWarnings = F, recursive = T)
    writexl::write_xlsx(ol$odat, paste0(opre, '.xlsx'))
    cli::cli_alert_success('{opre}.xlsx')  
  }
}



