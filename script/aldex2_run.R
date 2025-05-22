#' https://mp.weixin.qq.com/s/gNxqQ7rZdorvdd2rtUD_2g
#' (a)为每个样本生成 Dirichlet 分布的蒙特卡罗样本，
#' (b)使用对数比变换转换每个实例，然后
#' (c)返回两个样本(Welch’st, Wilcoxon)
#'   或 多个样本(glm, Kruskal-Wallace)测试的测试结果。
#'   该函数还估计了两个样本分析的效应大小。
#'
#'
aldex2_run <- function(dat, sgc, opre=NULL){
  library(ALDEx2)
  if(F){
    data(selex)
    #subset only the last 400 features for efficiency
    selex.sub <- selex[1:400,]
    selex.sub[1:4,1:4]
    # 分组
    conds <- c(rep("NS", 7), rep("S", 7))
    conds
    # 目前，aldex 函数仅限于双样本检验和单因素方差分析。
    x.all <- aldex(
      selex.sub, conds, mc.samples=16, test="t",
      effect=TRUE, include.sample.summary=FALSE,
      denom="all", verbose=FALSE)
    head(x.all)
    
    dat=round(eg.phylum*30000, 0); sgc=eg.sgc2
  }

read.table(text="列名, 解释 
rab.all,（特征中所有样本的clr中值）
rab.win.VdrFecal,（NS 组样本的clr中值）
rab.win.VdrCecal,（N 组样本的clr中值）
dif.btw,（VdrFecal 和 VdrCecal 组间的clr值中位数差）
dif.win,（VdrFecal 和 VdrCecal 组内的clr值最大差异中值）
effect, （效应量中值 diff.btw /max(diff.win）
overlap,（效应量包含0的比例）
we.ep,  （Welch’s t检验的p值）
we.eBH, （Welch’s t检验校正后的p值）
wi.ep,  （等级和检验的p值） Wilcoxon Rank Sum test 
wi.eBH, （秩和检验校正后的p值）
",sep = ',', stringsAsFactors = F, header= T) -> note
  ALDEx2::aldex(
    reads = dat, conditions = sgc$sg[colnames(dat)],
    denom = 'all',effect = TRUE ) -> res
  list(
    tibble::rownames_to_column(res),
    data = tibble::rownames_to_column(dat),
    Group= tibble::enframe(sgc$sg),
    说明 = note
  )  -> odat
  
  ####
  list(
    dat = dat, sgc=sgc, res= res, odat = odat
  ) -> ol
  if(is.null(opre)){ return(ol) }else{
    dir.create(dirname(opre), showWarnings = F, recursive = T)
    writexl::write_xlsx(ol$odat, paste0(opre, '.xlsx'))
    cli::cli_alert_success('{opre}.xlsx')  
  }
}



