#' froc_list
#' pROC
#' 
froc_list <- function(
   df, Y, namelist, color=NULL,
   legacy.axes = T,
   print.thres = T, # 显示 best 点
   txttype = c('auc_ci','auc_ci_pval')[1], # 显示在图片上 哪些信息
   title = NULL,    # 默认 "ROC curve (A & B)"
   line_size  = .7,
   line_alpha = 1,
   line_type   = 1,
   legend_size = 2,
   base_family = '',
   quiet = TRUE ,
   pexpand = TRUE,
   show_legend = TRUE,
   opre=NULL, ph=7, pw=7 
){
   if(F){ # test
      df <- data.frame(t(eg.phylum[c(1,7),]),y=eg.sgc2$sg[colnames(eg.phylum)])
      Y='y';col=c(1,2);legend_size = 2;base_family = ''
      line_size  = .7; line_alpha = 1;line_type   = 1;legend_size = 2; quiet = TRUE;
      pexpand = TRUE; show_legend = TRUE;print.thres=T;
      name=NULL;legacy.axes=T;title='ROC';color=NULL
      namelist<- list(a=c('Proteobacteria','Verrucomicrobia'),
                      b5454='Verrucomicrobia')
   }
   ### 
   if(is.null(color)) color <- ggsci::pal_igv()(50)[c(-4,-6)]
   color <- setNames(color[1:length(namelist)],names(namelist))
   
   df <- tibble::rownames_to_column(dplyr::rename(df,all_of(c(.y=Y))),'name') %>% # Y的列名
      dplyr::mutate(.y=as.factor(.y))
   purrr::map_dfr(names(namelist),function(type){ # type=names(namelist)[1]
      xname <- namelist[[type]] # x名
      if(length(xname)==1){
         dplyr::select(df,name,.y,value=dplyr::all_of(xname)) %>%
            dplyr::mutate(type=type)
      }else{
         model_glm <- glm(.y ~ .,data=df[,c(xname,'.y')],family='binomial')
         # df$.x <- predict(model_glm, newdata=df[,col,drop=F], type="response")
         dplyr::select(df,name,.y) %>% dplyr::mutate(value=model_glm$fitted.values,type=type)
      }
   }) %>% tidyr::pivot_wider(names_from = type,values_from = value) -> datl
   
   # res.roc <- pROC::roc(.y ~ .,ci=T,data=dplyr::select(datl,.y,all_of(names(namelist))))
   lapply(names(namelist), function(x){
      datx <- dplyr::select(datl,.y, all_of(c('.y',x=x)))
      rocx <- pROC::roc(.y ~ x, ci=T, print.thres='best', data = datx, quiet= quiet)
      if(rocx$auc<.5){ # AUC小于5的变量，调换分组的因子
         datx$.y <- forcats::fct_rev(datx$.y)
         rocx <- pROC::roc(.y ~ x, ci=T, data = datx, quiet= quiet) }
      rocx
   }) %>% setNames(names(namelist)) -> res.roc
   
   roc.pval <- function(r){
      b <- r$auc - .5
      se <- sqrt(pROC::var(r))
      z <- (b / se)
      2 * pt(-abs(z), df=Inf)  ## two-sided test
   }
   

   purrr::map_dfr(names(namelist), function(nn){ # nn=names(namelist)[1]
      if(length(namelist)==1 & length(namelist[[1]])==1 ){
         r=res.roc[[1]]
      }else{
         r=res.roc[[nn]]
      }
      auc <- round(r$auc,3);
      ci <- format(round(r$ci,3),nsmall=3)
      pROC::coords(r, x="best", transpose = TRUE,
                   input="threshold", best.method="youden") -> cutoff
      cf <- format(round(cutoff,3),nsmall=3)
      roc.pval <- roc.pval(r)
      list(name=nn,list=paste0(namelist[[nn]],collapse = ','),
           auc=auc,ci=paste0(ci[1],'-',ci[3]),
           threshold  =  cutoff['threshold'],
           specificity = cutoff['specificity'],
           sensitivity = cutoff['sensitivity'],
           #Cutoff  = paste0(cf[1],'(',paste0(cf[c(2,3)], collapse = ','),')'),
           cutoff_label = paste0(format(round(cutoff['threshold']    ,3),nsmall=3),'(',
                                 format(round(1-cutoff['specificity'],3),nsmall=3),',',
                                 format(round(cutoff['sensitivity']  ,3),nsmall=3),')'),
           P.value = roc.pval,
           P.value.format = ifelse(roc.pval<0.001, '<0.001', as.character(round(roc.pval,3))),
           color=color[nn])
   }) -> res.tab
   if(txttype == 'auc_ci_pval'){
      laba <- 'AUC(95% CI) P.value'
      dplyr::mutate(res.tab, label=paste0(name,': ',auc,'(',ci,') ',
                                          'P.value:', P.value.format #signif(P.value,digits = 3)
                                          )) -> res.tab
   }else{
      laba <- 'AUC(95% CI)'
      dplyr::mutate(res.tab, label=paste0(name,': ',auc,'(',ci,')')) -> res.tab
   }
   res.tab2 <- tibble(y=seq(0,(length(namelist)-1)*.05,.05),
                      #name=c(res.tab$name,'.'),label=c(res.tab$label,'AUC(CI)'))
                      name=c(res.tab$name),label=c(res.tab$label))
   # color['.'] <- 'black'
   if(is.null(title)) 
      title <- paste0('ROC curve (',paste0(unique(df$.y),collapse = ' & '),')')
   pROC::ggroc(res.roc, legacy.axes=legacy.axes, linetype = line_type,
               size = line_size, alpha = line_alpha) + 
      geom_abline(intercept = ifelse(legacy.axes,0,1), slope = 1,
                  color = "darkgrey", linetype = "solid") +
      geom_text(aes(x=ifelse(legacy.axes,.95,.05), y=y+.05, label=label),
                data = res.tab2, hjust=1, show.legend = F) +
      annotate('text',x=ifelse(legacy.axes,.95,.05), y=length(namelist)*.05+0.05,
               label=laba, hjust=1,color='black',fontface = "bold") +
      labs(title = title) + theme_bw(base_family = base_family) +
      theme(text = element_text(family = base_family),
            panel.grid = element_blank(),
            plot.title = element_text(hjust = .5),
            legend.title = element_blank()) +
      guides(color=guide_legend(reverse = TRUE,
                                override.aes = list(size=legend_size))) -> p
   if(print.thres) # 加上cutoff
      p + geom_point(aes(x=1-specificity,y=sensitivity),dat=res.tab,show.legend = F) +
        ggrepel::geom_text_repel(aes(x=1-specificity,y=sensitivity,
                                 label=cutoff_label), #nudge_y = -.01,
                                 segment.size = .1,
                                 dat=res.tab,show.legend = F) -> p
   if(!pexpand) 
      p <- p + scale_y_continuous(expand = expansion(mult = c(0,0))) +
               scale_x_continuous(expand = expansion(mult = c(0,0)))
   if(length(namelist)>1) p <- p + scale_color_manual(values = color)
   if(!show_legend)       p <- p + theme(legend.position = 'none')
   # plot(pROC::ci.se(res.roc),type="shape", col='blue')
   if(is.null(opre)){
      return(list(p=p,dat=list(AUC=res.tab,data=datl,rawdata=df)))
   }else{
      dir.create(dirname(opre),showWarnings = F,recursive = T)
      ggsave(paste0(opre, '.pdf'), p, device = cairo_pdf, width = pw, height = ph)
      writexl::write_xlsx(list(
         AUC=res.tab, data=datl, Rawdata=df
      ), path = paste0(opre, '.xlsx'))
      cat(opre,'.xlsx .pdf\n')
   }
}

