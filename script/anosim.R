#' anosim 相似性分析
#' 输入矩阵和分组
#' ?vegan::anosim
v_anosim <- function(
  dist, sgc, permutations=999, violin = F,
  notch=F,
  laby='Rank Dissimilarity', title='ANOSIM',
  las = 1, Between_col = 'gray',
  ptype=c('box','violin','violin2')[1],
  base_size = 17,
  opre = NULL, pw = 6, ph = 7){
  if(F){
    dist=eg.dist;sgc=eg.sgc;permutations=999
    violin = F;notch=T;laby='Rank Dissimilarity';title='ANOSIM';
    las=2;mar=c(6,4,4,2)
    Between_col = 'gray'
  }
  gpn   <- unique(sgc$sg) # 分组名
  resa <- vegan::anosim(as.dist(dist), sgc$sg[colnames(dist)],
                        permutations = permutations) # 整体分析
  resdf <- tibble::tibble(Group=paste(gpn,collapse = ' vs '),
                          permutations=resa$permutations,
                 R=resa$statistic,p.value=resa$signif)
  if(length(gpn)>2){ # 两两之间比较
    comgp <- combn(gpn,2)
    resdf <- rbind(
      resdf,
      purrr::map2_dfr(comgp[1,],comgp[2,],function(a,b){  # a='B1';b='B2'
        sgn <- sgc$sg[sgc$sg %in% c(a,b)]
        set.seed(123)
        res <- vegan::anosim(as.dist(as.matrix(dist)[names(sgn),names(sgn)]),sgn,
                             permutations = permutations)
        list(Group=paste(a,'vs',b), permutations=res$permutations,
             R=res$statistic, p.value=res$signif)
      })
    )
  }
  # par(mar=c(7,4,4,2));plot(resa,las=2)
  # p <- ggplotify::as.ggplot(function(){ par(mar=mar)
  #    plot(resa,col=c('gray',sgc$gc[levels(resa$class.vec)[-1]]),las=las,ylab="Rank Dissimilarity",xlab='',main="ANOSIM") })
  # data.frame(x=factor(resa$class.vec,levels = c('Between',gpn)),y=resa$dis.rank)
  pd <- data.frame(x=factor(resa$class.vec,levels = c('Between',unique(sgc$sg))),y=resa$dis.rank) %>%
    dplyr::mutate(width=(table(.$x)/nrow(.))[x])
  ggplot(pd,aes(x=x,y=y)) +
    geom_violin(position = position_dodge(1),size=1, scale = 'width', trim = F) +
    ggbeeswarm::geom_quasirandom(
      aes(color=x), method='quasirandom', show.legend = F) -> pv2
  ggplot(pd,aes(x=x,y=y,fill=x)) +  #geom_violin(trim=F,show.legend = F) +
    stat_ydensity(trim = F,geom = 'violin',scale = 'width',width=.7,show.legend = F) +
    geom_boxplot(width=0.2,fill='white',position=position_dodge(0.9),notch=F) -> pv #绘制箱线图
  ggplot(pd,aes(x=x,y=y,fill=x)) +    #fill=x))  空心箱线图
    stat_boxplot(geom ='errorbar',width=.5,linetype=1) +
    geom_boxplot(notch=F,show.legend = F,width=.7) -> pb
  # geom_jitter(aes(x,y))+  guides(color='none')
  # ggplot(pd,aes(x=x,y=y,fill=x)) +
  #   stat_boxplot(geom ='errorbar',linetype=1,width=table(pd$x)/nrow(pd)) +
  #   geom_boxplot(notch=notch,show.legend = F,width=table(pd$x)/nrow(pd))
  # plot(resa)
  mcols <- c(Between=Between_col, sgc$gc[levels(resa$class.vec)[-1]])
  fp <- function(p1){
      p1 + scale_fill_manual(values = mcols) +
      scale_color_manual(values = mcols) +
      labs(x='', y=laby, title = title,
           subtitle = paste('R =',round(resa$statistic,4),
                            ',P =',resa$signif)) +
      theme_bw(base_size = base_size) +
      theme(plot.title = element_text(hjust = .5),
            plot.subtitle = element_text(hjust = .5),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45,hjust = 1)) }
  resdata <- list(ANOSIM=resdf,Group=tibble::enframe(sgc$sg,'name','group'),
                  Dist=tibble::rownames_to_column(as.data.frame(dist),'.ID'))
  
  # p; plot(resa)
  if(is.null(opre)){
    return(list(
      res.ano = resa, resdf = resdf, dat=resdata,
      p = fp(pb), pv=fp(pv), pv2=fp(pv2),
      plot = function(){
        plot(resa,col=c(Between_col, sgc$gc[levels(resa$class.vec)[-1]]),
             las = las, ylab=laby, xlab='', main = title)
      }
    ))
  }else{ # 输出结果
    if(ptype == 'box'){
      pp <- fp(pb) 
    }else if(ptype =='violin'){
      pp <- fp(pv) 
    }else if(ptype =='violin2'){
      pp <- fp(pv2) 
    }else{stop('ptype = worng!!')}
    # if(violin) pp <- fp(pv) else pp <- fp(pb)
    dir.create(dirname(opre), showWarnings = F,recursive = T)
    ggsave(paste0(opre,'.pdf'),plot = pp,device = cairo_pdf,width = pw,height = ph)
    writexl::write_xlsx(resdata,path = paste0(opre,'.xlsx'))
    cat(opre,'\n')
  }
}; # v_anosim(dist=eg.dist,sgc=eg.sgc,permutations=999)
