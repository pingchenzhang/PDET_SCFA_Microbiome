#'  alpha
#'  
diff_alpha_box <- function(
  dat.alpha, # 
  sgc, # list()
  parameter=T, paired=F, exact = T,
  cname=c(Richness='Richness', Chao1='chao',
          ACE='ace', Shannon='shannon diversity',
          Simpson = 'simpson diversity',
          PD_whole_tree = 'PD_whole_tree'),
  opre = NULL, pw = 9, ph = 6
){
  # dat.alpha <- res$alpha
  if(length(setdiff(names(sgc$sg),rownames(dat.alpha)))>0)
    stop(setdiff(names(sgc$sg),rownames(dat.alpha)),'not in dat')
  if(length(setdiff(names(cname),colnames(dat.alpha)))>0)
    stop(setdiff(names(cname),colnames(dat.alpha)),'not in dat')
  ldngp <- length(unique(sgc$sg))==2
  testmethod <- ifelse(parameter,
                       ifelse(ldngp,'t.test','anova'),
                       ifelse(ldngp,'wilcox','kruskal'))
  ####
  pdat <- tibble::rownames_to_column(dat.alpha,'.ID') %>%
    dplyr::select(.ID,names(cname)) %>%
    dplyr::filter(.ID %in% names(sgc$sg)) %>%
    dplyr::mutate(.gp = sgc$sg[.ID],.gp = factor(.gp,levels = unique(sgc$sg))) # 
  pval <- ggpubr::compare_means(val ~ .gp,data = tidyr::gather(pdat,name,val,-.gp,-.ID) %>%
                                  dplyr::mutate(val=ifelse(is.na(val),0,val)),
              group.by = 'name', paired=paired, exact = exact, method = testmethod) %>%
    dplyr::select(-.y.) %>%
             dplyr::mutate(p.signif=ifelse(p.signif=='ns','',p.signif))
  purrr::map(names(cname),function(x){
    pval2 <- dplyr::filter(pval,name == x)
    sutitle <- paste0('method = ',stringr::str_to_title( testmethod), ' P = ',
                      round(pval2$p, 3), ifelse(pval2$p.signif!='ns',pval2$p.signif,''))
    ggplot(pdat[,c(x,'.gp')], aes_string(x='.gp', y=x, fill = '.gp'))+
      stat_boxplot(geom ='errorbar') +  geom_boxplot(show.legend = F)+
      scale_fill_manual(values=sgc$gc) +
	  # 	  #geom_point(aes(color = .gp))+
	  # geom_jitter(aes(fill=.gp),width =0.2,shape = 21,size=1)+
      labs(x='',y=cname[x],subtitle = sutitle) + theme_bw() + 
      theme(panel.grid = element_blank(),legend.title = element_blank(),
            axis.text.x=element_text(angle=45, hjust=1)) #
  }) %>% purrr::set_names(names(cname)) -> pall
  # 

  if(!is.null(opre)){
    ggsave(paste0(opre,'.pdf'),patchwork::wrap_plots(pall),device = cairo_pdf,
           width = pw,height = ph)
    writexl::write_xlsx(list(Pvalue=pval,Rawdata=pdat),paste0(opre,'.xlsx'))
    cat(paste0(opre,'.pdf xlsx\n'))
  }
  return(list(p=pall,pall=patchwork::wrap_plots(pall),
              dat=list(Pvalue=pval,Rawdata=pdat)))
}
