#' Alpha多样性计算
#' vegan/ picante
alpha_div <- function(dat, tree = NULL, base = exp(1),
                        include.root = F, gini_simpson=F) {
  if(F){
    # x=t(res$otu.count);tree=res$otufa.tree;base = exp(1);include.root = F
  }
  ###
  if(!suppressMessages(require(vegan)))   install.packages('vegan')
  if(!suppressMessages(require(picante))) install.packages('picante')
  x <- t(dat)
  est  <-  vegan::estimateR(x)
  Observed <- est[1, ]
  Shannon  <- vegan::diversity(x, index = 'shannon', base = base)
  if(gini_simpson){
    Simpson  <- vegan::diversity(x, index = 'simpson')    # Gini-Simpson 指数
  }else{
    Simpson  <- 1 - vegan::diversity(x, index = 'simpson') # 
  }
  Invsimpson<- vegan::diversity(x, index = 'invsimpson')
  Pielou <-   Shannon / log(Observed, base)
  Goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  #Fisher_alpha <- vegan::fisher.alpha(x)
  result <- data.frame(Reads=rowSums(x),Observed, Shannon, 
                       Simpson,Invsimpson,Pielou,
                       Chao1=est[2,], ACE=est[4,],
                       Goods_coverage,#Fisher_alpha,
                       check.names = F,stringsAsFactors = F)
  if (!is.null(tree)) {
    PD <- picante::pd(x, tree, include.root = include.root)[1]
    names(PD) <- 'PD_whole_tree'
    result <- cbind(result, PD)
  }
  rownames(result) <- rownames(x)
  result
}

