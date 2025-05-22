
#' 差异分析 alpha多样性
diff_alpha_box <- function(
  dat.alpha, # alpha多样性表格
  sgc, # 样本分组信息 list()
  parameter=T, paired=F, exact = T,
  cname=c(Richness='Richness', Chao1='chao',
          ACE='ace', Shannon='shannon diversity',
          Simpson = 'simpson diversity',
          PD_whole_tree = 'PD_whole_tree'),
  opresingle = NULL, opre = NULL, pw = 9, ph = 6
){
  if(F){
    dat.alpha = resi$alpha; sgc=sgci; cname = alphaname; paired = opt$paired;
    parameter = T; opresingle = od1; opre = paste0(od1,'整体')
  }
  # dat.alpha <- res$alpha
  if(length(setdiff(names(sgc$sg),rownames(dat.alpha)))>0)
    stop(setdiff(names(sgc$sg),rownames(dat.alpha)),'not in dat')
  if(length(setdiff(names(cname),colnames(dat.alpha)))>0)
    stop(setdiff(names(cname),colnames(dat.alpha)),'not in dat')
  ldngp <- length(unique(sgc$sg))==2
  (testmethod <- ifelse(parameter,ifelse(ldngp,'t.test','anova'),ifelse(ldngp,'wilcox','kruskal')))
  pdat <- tibble::rownames_to_column(dat.alpha,'.ID') %>%
    dplyr::select(.ID,names(cname)) %>%
    dplyr::filter(.ID %in% names(sgc$sg)) %>%
    dplyr::mutate(.gp = sgc$sg[.ID],.gp = factor(.gp,levels = unique(sgc$sg))) # 合并分组与数据
  pval <- ggpubr::compare_means(val ~ .gp,data = tidyr::gather(pdat,name,val,-.gp,-.ID) %>%
                                  dplyr::mutate(val=ifelse(is.na(val),0,val)),
              group.by = 'name', paired=paired,exact = exact, method = testmethod) %>% dplyr::select(-.y.) %>%
             dplyr::mutate(p.signif=ifelse(p.signif=='ns','',p.signif))
  purrr::map(names(cname),function(x){
    pval2 <- dplyr::filter(pval,name == x)
    sutitle <- paste0('method =',testmethod,' p =',pval2$p.format,ifelse(pval2$p.signif!='ns',pval2$p.signif,''))
    ggplot(pdat[,c(x,'.gp')], aes_string(x='.gp', y=x, fill = '.gp'))+
      stat_boxplot(geom ='errorbar') +  geom_boxplot(show.legend = F)+
      scale_fill_manual(values=sgc$gc) +
	  # 	  #geom_point(aes(color = .gp))+
	  # geom_jitter(aes(fill=.gp),width =0.2,shape = 21,size=1)+
      labs(x='',y=cname[x],subtitle = sutitle) + theme_bw() + 
      theme(panel.grid = element_blank(),legend.title = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1)) #210812加入,旋转x轴标签45°
  }) %>% purrr::set_names(names(cname)) -> pall
  # 事后检验
  if(length(unique(sgc$sg))>2){
    reshoc <- purrr::map(names(cname),function(y){ # y='ACE'
      post_hocfun(value = pdat[[y]], gp = pdat$.gp, method = testmethod)
    }) %>% purrr::set_names(names(cname))
  }else{reshoc <- NULL}
  if(!is.null(opresingle)){
    dir.create(dirname(opresingle),showWarnings = F,recursive = T)
    purrr::walk(names(pall),function(x){
      ggsave(paste0(opresingle,x,'.pdf'),pall[[x]],device = cairo_pdf,width = 6,height = 6)
      cat(paste0(opresingle,x,'.pdf\n'))
    })
  }
  if(!is.null(opre)){
    ggsave(paste0(opre,'.pdf'),patchwork::wrap_plots(pall),device = cairo_pdf,
           width = pw,height = ph)
    writexl::write_xlsx(list(Pvalue=pval,Rawdata=pdat),paste0(opre,'.xlsx'))
    if(!is.null(reshoc)){ sink(paste0(opre,'两两比较.txt'));print(reshoc); sink() }
    cat(paste0(opre,'.pdf xlsx\n'))
  }
  return(list(p=pall,pall=patchwork::wrap_plots(pall),
              dat=list(Pvalue=pval,Rawdata=pdat),reshoc=reshoc))
}
