#'
#' randomforest
#' 
randomforest <- function(
  dat, sgc, ntree=999, seed=NULL,
  main = "Importance",
  opre = NULL, n.var = 20, pw=8, ph=7,
  ...
){
  # dat=eg.phylum; sgc=eg.sgc;ntree=2000;seed=NULL;...=NULL
  require(randomForest)
  tibble::rownames_to_column(dat,'.id') %>%
    tidyr::gather(SampleID, val, -.id) %>%
    tidyr::spread(.id, val) %>%
    dplyr::mutate(Group = as.factor(sgc$sg[SampleID]),.after=1) -> dat2
  ## 
  if(!is.null(seed)) set.seed(seed)
  randomForest::randomForest(
    Group~., data = dplyr::select(dat2, -SampleID),
    ntree=ntree, importance=T, ...) -> rf.res

  ## result
  if(!is.null(opre)){
    dir.create(dirname(opre), showWarnings = F, recursive = T)
    list(
      Importance = tibble::rownames_to_column(as.data.frame(
        randomForest::importance(rf.res))),
      Data = dat2
    ) -> odat
    writexl::write_xlsx(odat, paste0(opre,'Data.xlsx'))
    capture.output(print(rf.res), file = paste0(opre,'Result.txt'))
    n.var <- min(nrow(odat$Importance), n.var)
    if(n.var>1){
      cairo_pdf(paste0(opre,'Importance.pdf'), width = pw, height = ph)
      # mean decrease accuracy (MDA) and mean decrease Gini (MDG)
      randomForest::varImpPlot(rf.res, n.var = n.var, main = main)
      dev.off()     
    }
    message('to:', opre)
  }
  return(list(
    res = rf.res,
    importance = randomForest::importance(rf.res)
  ))
}

