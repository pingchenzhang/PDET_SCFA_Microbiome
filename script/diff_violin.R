#' diff_fun 差异分析
#' 计算数据和画图
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' ggsignif::geom_signif()
#' 
diff_violin <- function(
   dat, sgc, name,
   test  = c('param','nonparam','kruskal','anova','t.test','wilcox')[1], # auto 
   pair.test= c('posthoc','pairwise','t.test','wilcox','TukeyHSD','DunnTest','NemenyiTest')[1],
   pair.adjust = c('none')[1],
   var.equal = F, paired = F, exact = NULL, # test 参数
   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
   symbols   = c("***", "**", "*", ""), # P值转signif
   ptype = c('bar','box','violin4')[1], # 画图类型
   title=NULL, laby='Relative Abundance (mean±SE)', laby.box='Abundance', # label
   sig.type = c('p.signif','p.format','p.format.signif')[1],
   p_limit=0.05, coef = 1.5, errorbar.width = .5,
   violin_trim = T, 
   bee.cex = 1, bee.alpha=1,
   priority = c("ascending","descending", "density", "random", "none")[3],
   tip_length = 0.01, tip_vjust=0.1, sig_vjust=0.01, # tip调整参数
   angle.x=0, base.size=11, size.signif=8, size.border=.5, # 字体调整
   kh=.5, kw=.5, opre=NULL, pw=6, ph=6, # 输出图形参数
   warn = F
){
   if(F){ # test #####
      dat=eg.genus;sgc=eg.sgc;name='Fusobacterium';ptype = c('bar','box','violin','bee')[2]
      # sgc$sg[1:3] <- c('B_1'='B1-B1','B_2'='B1-B1','B_3'='B1-B1');sgc$gc[1] <- c('B1-B1'='red'); dfi$val[1] <- 1
      test=c('anova','kruskal')[1];pair.test='posthoc';  bee.cex = 1; violin_trim=T
      laby='Relative Abundance (mean±SE)';laby.box='Abundance';title='plot'
      var.equal=F;paired =F;exact=F;angle.x=0; base.size=8; size.signif=8; size.border=.5;
      cutpoints = c(0, 0.001, 0.01, 0.05, 1); symbols = c("***", "**", "*", "")
      p_limit=0.05;coef = 1.5;tip_length = 0.01;sig_vjust=0.01;
      sig.type = c('p.signif','p.format','p.format.signif')[1];
      kh=.5; kw=.5;pre=NULL; pw=6; ph=6; warn=T
   }
   ol <- list(name=name) # 输出list
   if(is.null(title)) title <- name
   ptypes <-  c('bar','box','violin4')
   if(! (ptype %in% ptypes)) stop('ptype must in:',paste0(ptypes,collapse = ', '))
   ## 挑出数据 sp gp val #### dfi$val[1] <- 1
   if(length(name)>1 & name %in% rownames(dat)) stop('name 必须是长度为一的字符，必须是dat中的一个行名。')
   gpname <- unique(sgc$sg) # 分组名字
   gpnewname <- setNames(gpname,stringr::str_replace_all(gpname,'-','_'))
   dfi <- tibble(sp=names(sgc$sg),gp=sgc$sg,val=unlist(dat[name,names(sgc$sg)])) %>%
      dplyr::mutate(gp2=factor(stringr::str_replace_all(gp,'-','_'),levels = names(gpnewname))) %>%
      dplyr::mutate(gp =factor(gp,levels = gpname))
   ## 确定分析的参数: ##########
   gpnn <- length(unique(sgc$sg)) # 分组数
   sp_isequal <- (length(unique(table(sgc$sg)))==1)# 每个分组的样本数是否相等
   test.typen <- c('kruskal','anova','t.test','wilcox')
   pair.typen <-  c('t.test','wilcox','TukeyHSD','DunnTest','NemenyiTest')
   parr.typestop <- paste0('pair.type must in: pairwise, posthoc, ',paste(pair.typen,collapse = ', '))
   sig.typen = c('p.signif','p.format','p.format.signif')
   if(!sig.type %in% sig.typen) stop('sig.type must in:',paste0(sig.typen,collapse = ', '))
   if(sig.type!='p.signif'){ size.signif = max(size.signif,1); sig_vjust = sig_vjust+0.02 }
   
   # 差异计算参数 ####
   if(startsWith(test,'p')){       method <- ifelse(gpnn>2,'anova','t.test')  # 参数检验
   }else if(startsWith(test,'n')){ method <- ifelse(gpnn>2,'kruskal','wilcox')# 非参数检验
   }else if(test %in% test.typen){ method <- test
   }else{ stop('must in: ',paste(test.typen, collapse = ', ')) }
   # 两两比较参数 ####
   if(method %in% c('t.test','anova')){
      if(pair.test=='pairwise'){           pair.method <- 't.test'
      }else if(pair.test=='posthoc'){      pair.method <- 'TukeyHSD'
      }else if(pair.test %in% pair.typen){ pair.method <- pair.test
      }else{ stop(parr.typestop) }
   }else if(method %in% c('wilcox','kruskal')) {
      if(pair.test=='pairwise'){           pair.method <- 'wilcox'
      }else if(pair.test=='posthoc'){ pair.method <- ifelse(sp_isequal,'NemenyiTest','DunnTest')
      }else if(pair.test %in% pair.typen){ pair.method <- pair.test
      }else{ stop(parr.typestop) }
   }
   if(warn) cat('method:',method,'pair.method:',pair.method,'\n')
   ol$method <- method
   ol$pair.method <- pair.method
   # P值转化为星号
   pval2stars <- function(p) unclass(symnum(p,corr=FALSE,na=FALSE,cutpoints=cutpoints,symbols=symbols))
   # 长数据转宽数据
   .f_long2w <- function(dfl, scol){ # scol=c('mean','SD','n','SE')
      dplyr::select(dfl, gp, all_of(scol)) %>%
         tidyr::pivot_longer(cols=scol, names_to='name', values_to='val') %>%
         dplyr::mutate(name=paste0(name,'(',gp,')'),name=factor(name,levels = unique(name))) %>% 
         dplyr::select(-gp) %>%  tidyr::pivot_wider(names_from='name',values_from='val')
   }
   
   ## 计算P值
   .f_test <- function(m,data,x){
      if(m == 't.test'){
        t.test(val ~ gp,data=data, var.equal=var.equal, paired=paired)
      }else if(m=='anova'){
        oneway.test(val ~ gp,data=data, var.equal=var.equal)
      }else if(m=='wilcox'){
        # wilcox.test(val ~ gp, data=data, paired=paired, exact=exact)
        unique(data$gp) -> .gpn
        wilcox.test(x = data$val[data$gp == .gpn[1]],
                    y = data$val[data$gp == .gpn[2]],
                    paired=paired, exact=exact)
      }else if(m=='kruskal'){
        kruskal.test(val ~ gp,data=data)
      }
   }
   ## 两两比较
   .f_pairwise <- function(dat=dat, method=c('t.test','wilcox')[1]){
      gpcombn <- combn(unique(sgc$sg),2)
      purrr::map2_dfr(gpcombn[1,],gpcombn[2,],function(g1,g2){
         # g1=gpcombn[1,1];g2=gpcombn[2,1]
         .f_test(m = method, data = dplyr::filter(dat, gp %in% c(g1,g2)) ) %>%
            broom::tidy() %>% dplyr::mutate(group1=g1,group2=g2,.before=1)
      }) %>% dplyr::mutate(p.signif=pval2stars(p.value),.after='p.value')
   }
   ## 事后检验
   .f_tidyposthoc <- function(restest, mname){
      # print(DescTools::as.matrix.xtabs(restest))
      tibble::rownames_to_column(as.data.frame(restest[[1]])) %>% 
         tibble::as_tibble() %>%
         tidyr::separate(rowname,into = c('group1','group2'),sep = '-',remove = F) %>%
         dplyr::rename(p.value=pval) %>%
         dplyr::mutate(p.signif=pval2stars(p.value),.after='p.value') %>%
         dplyr::mutate(group1=gpnewname[group1],group2=gpnewname[group2],
                       method=stringr::str_trim(attr(restest,mname)))
   }
   ## 计算两两P值 
   .f_post_pair <- function(dat, m){
      # pacman::p_load(PMCMR,DescTools)
      if(m=='t.test'){  # 1. pairwise t.test
         # pairwise.t.test(dfi$val,dfi$gp,p.adjust.method = 'none',pool.sd = F,var.equal=var.equal,paired=paired)
         .f_pairwise(dat,method = 't.test')
      }else if(m=='wilcox'){ #2. pairwise wilcox.test
         # pairwise.wilcox.test(dfi$val,dfi$gp,p.adjust.method = 'none',paired=paired,exact =exact)
         .f_pairwise(dat,method = 'wilcox')
      }else if(m=='TukeyHSD'){ # 1. TukeyHSD
         # stats::TukeyHSD(aov(val ~ gp2, data = dat))
         DescTools::PostHocTest(aov(val ~ gp2, data = dat)) %>%
            .f_tidyposthoc(mname='method.str')
      }else if(m=='DunnTest'){   # 2. 分组内的样本数不相等
         # PMCMR::posthoc.kruskal.dunn.test(val ~ gp, data = dat,p.adjust.method='none')
         DescTools::DunnTest(val ~ gp2, data = dat, method="none")  %>%
            .f_tidyposthoc(mname = 'main')
      }else if(m=='NemenyiTest'){ # 3. 分组内的样本数相等
         # PMCMR::posthoc.kruskal.nemenyi.test(val ~ gp, data = dat)
         DescTools::NemenyiTest(val ~ gp2, data = dat) %>%
            .f_tidyposthoc(mname = 'main')
      }else{
         stop('method must in: t.test,wilcox,TukeyHSD,DunnTest,NemenyiTest')
      }
   }
   theme_bw2 <- function(){
      theme_bw(base_size = base.size) + 
         theme(panel.grid = element_blank(),
               legend.title = element_blank(),
               plot.title   = element_text(hjust = .5),
               axis.text.x = element_text(angle = angle.x, hjust = ifelse(angle.x==0,.5,1)))
   }
   ## 计算全部分组P值 ####
   pval.all <- .f_test(m = method,data = dfi) %>% broom::tidy() %>%
      dplyr::mutate(Name = name,.before=1) %>%
      dplyr::mutate(p.format=ifelse(p.value<0.001,'<0.001',as.character(round(p.value,3))),
                    p.signif=pval2stars(p.value),
                    p.format.signif=paste0(p.format,' ',p.signif),.after=p.value)
   ol$p.format.signif <- paste0(method,ifelse(pval.all$p.value<0.001,' p ',' p ='),
                               pval.all$p.format.signif)
   ol$pval.all <- pval.all
   ## 计算两两P值 OR 事后检验 ####
   if(gpnn > 2){
      pval.pair  <- .f_post_pair(dat = dfi,m = pair.method) %>%
         dplyr::mutate(Name=ol$name, .before=1) %>%
         dplyr::mutate(p.format=ifelse(p.value<0.001,'<0.001',as.character(round(p.value,3))),
                       p.format.signif=paste0(p.format,' ',p.signif), .after=p.value)
      if(pair.adjust != 'none'){
         # 矫正p值
         pval.pair <- dplyr::mutate(pval.pair, p.adjust = p.adjust(p.value,method = pair.adjust)) %>%
            dplyr::mutate(p.adjust.method = pair.adjust,
               p.adjust.format=ifelse(p.adjust<0.001,'<0.001',as.character(round(p.adjust,3))),
               p.adjust.signif=pval2stars(p.adjust),
               p.adjust.format.signif=paste0(p.adjust.format,' ',p.adjust.signif), .after=p.adjust)
         pval.pairp <- dplyr::mutate(
            pval.pair, p.value = p.adjust, p.format = p.adjust.format,
            p.format.signif = p.adjust.format.signif) %>%
            dplyr::filter(p.value <= p_limit)  # 两两有差异的 组
      }else{
         # 未矫正直接筛选原始 p.value
         pval.pairp <- dplyr::filter(pval.pair, p.value <= p_limit)  # 两两有差异的 组
      }
   }else{
      pval.pair <-  dplyr::mutate(pval.all, group1=gpname[1], group2=gpname[2]) %>%
         dplyr::mutate(Name=ol$name, .before=1) %>%
         dplyr::mutate(p.format=ifelse(p.value<0.001,'<0.001',as.character(round(p.value,3))),
                       p.format.signif=paste0(p.format,' ',p.signif),.after=p.value)
      pval.pairp <- dplyr::filter(pval.pair, p.value <= p_limit)  # 两两有差异的 组
   }
   ol$pval.pair <- pval.pair
   ## 画图数据 ###      
   dplyr::group_by(dfi, gp) %>%
      dplyr::summarise(mean=mean(val), SD=sd(val), n=dplyr::n(), SE=SD/sqrt(n),.groups = 'drop') %>%
      dplyr::mutate(Name = ol$name, gp=factor(gp, levels = unique(sgc$sg)),.before=1) -> res.meanl
   ys.mina <- max(dfi$val) # 最高点
   ys.minm <- max(res.meanl$mean + res.meanl$SE) # mean + se最高点
   gp2lev <- setNames(1:length(levels(dfi$gp)),levels(dfi$gp)) # 分组对应的水平
   # gp2lev <- setNames(as.numeric(res.meanl$gp),as.character(res.meanl$gp))
   ol$ys.mina <- ys.mina;  ol$ys.minm <- ys.minm
   if(nrow(pval.pairp)>0){
      dfp.sig <- tibble::rowid_to_column(pval.pairp,'yid') %>%
         dplyr::mutate(ynm = ys.minm*(1+yid*tip_vjust), yna = ys.mina*(1+yid*tip_vjust),.after='yid',
                       xn = purrr::map2_dbl(group1,group2,~(gp2lev[.x]+gp2lev[.y])/2))
      ol$dfp.sig <- dfp.sig
   }
   .p_add.sig <- function(
      p, yn='yna',
      ys.min = ys.mina,
      labelcol = 'p.signif'
   ){
      if(nrow(pval.pairp)>0){
         p+ geom_segment(aes(x=group1, xend=group2, y=!!sym(yn),
                             yend=!!sym(yn)),data = dfp.sig) +
            geom_segment(aes(x=group1, xend=group1, y=!!sym(yn),
                             yend=!!sym(yn)-ys.min*tip_length),data = dfp.sig) +
            geom_segment(aes(x=group2, xend=group2, y=!!sym(yn),
                             yend=!!sym(yn)-ys.min*tip_length),data = dfp.sig) +
            geom_text(aes(x=xn,y=!!sym(yn)+(ys.min*sig_vjust),
                          label=!!sym(labelcol)),
                      size=size.signif, data = dfp.sig)
      }else{ p }
   }

   
   ## 计算 均值 SD SE ####
   if(ptype %in% c('bar')){
      .f_long2w(res.meanl, scol = c('mean','SD','n','SE')) %>%
         dplyr::mutate(Name=name,.before=1) -> res.meanw
      ol$res.meanl <- res.meanl; ol$res.meanw <- res.meanw; ol$res.dat <- res.meanw
      ## 1. plot Bar mean SE pval.pair ####
      ggplot(res.meanl) +
         labs(x = "",y = laby,title = title) + theme_bw2() +
         geom_bar(aes(x=gp,y = mean,fill=gp),size=size.border,position = "dodge",
                  stat="identity",color = "black") +
         geom_errorbar(aes(x=gp,y = mean,ymin=mean-SE, ymax=mean+SE,group=gp),
                       width=.2,position=position_dodge(.9)) +
        #scale_y_continuous(expand = expansion(mult = c(0,.1))) +
        scale_y_continuous(expand = expansion(mult = c(0,.1))) +
        scale_fill_manual(values = sgc$gc)  -> pb1
      .p_add.sig(pb1,yn = 'ynm',ys.min = ys.minm, labelcol = sig.type) -> pb2
      ol$p <- pb2
      ol$dat <- list(mean_SE=res.meanl,
                     P.value=pval.all,
                     P.value_pairwise=pval.pair,
                     Rawdata=dplyr::select(dfi,-gp2))
   # }else if(ptype %in% c('box', 'violin')){
   }else{
      ## 计算 最小，4分位，中位数，最大，异常点 ----
      purrr::map_dfr(unique(sgc$sg),function(gpi){ # gpi='B';dfi$val[1]=2
         datli  <- dplyr::filter(dfi, gp==gpi)$val 
         #datout <- grDevices::boxplot.stats(datli,coef = coef)$out  # 计算异常值
         #datli2 <- setdiff(datli,datout) # 排除异常值后的数据
         # https://waterdata.usgs.gov/blog/boxplots/
         # https://ggplot2.tidyverse.org/reference/geom_boxplot.html
         quartiles <- as.numeric(quantile(datli,probs = c(0.25, 0.5, 0.75), na.rm=T))
         IQR <- diff(quartiles[c(1,3)])
         upper_dots <- datli[datli > (quartiles[3] + coef*IQR)]
         lower_dots <- datli[datli < (quartiles[1] - coef*IQR)]
         datout <- c(upper_dots,lower_dots)
         datli2 <- setdiff(datli,datout) # 排除异常值后的数据
         list(gp=gpi, min=min(datli2),
              #P25 = quantile(datli2,0.25), median = median(datli2), P75 = quantile(datli2,0.75),
              P25 = quartiles[1], median = quartiles[2], P75 = quartiles[3],
              max = max(datli2), n = length(datli),
              outlier=ifelse(length(datout) == 0, '',
                             paste0(paste0(names(datout),':',round(datout,5)),collapse = '; ')))
      }) -> res.medl
      bind_cols(.f_long2w(res.medl,scol = c('min','P25','median','P75','max','n')),
                .f_long2w(res.medl,scol = c('outlier'))) %>%
         dplyr::mutate(Name=name,.before=1) -> res.medw
      ol$res.medl <- res.medl; ol$res.medw <- res.medw; ol$res.dat <- res.medw
      ## 2. plotbox ####
      p0 <- ggplot(data = dfi, aes(x = gp, y = val)) +
         labs(x='',y=laby.box, title = title, subtitle = '') +
         scale_y_continuous(expand = expansion(mult = c(.05,.1))) +
         scale_fill_manual(values=sgc$gc) +
         scale_color_manual(values=sgc$gc) + theme_bw2()
      # 箱线图----
      p0 + stat_boxplot(aes(group=gp),geom = 'errorbar',coef = coef,width=errorbar.width) +
         geom_boxplot(aes(fill=gp),coef = coef, outlier.size = .5) -> px1# ;px1

      # 小提琴4 加箱线图 ----
      p0 + geom_violin(aes(fill=gp, color=gp), position = position_dodge(1),
                       size=.5,alpha=.5, scale = 'width', trim = violin_trim) +
         geom_boxplot(width=0.18, fill='black',
                      outlier.size = .2) +
         stat_summary(fun=function(x) boxplot.stats(x, coef = 1.5)$stats[3],
                      geom="point", color='white', size=2) -> pv4

      # 加 p.signif
      if(ptype=='box'){
         ol$p <- .p_add.sig(px1,yn = 'yna',ys.min = ys.mina, labelcol = sig.type)
      }else if(ptype=='violin4'){
         ol$p <- .p_add.sig(pv4,yn = 'yna',ys.min = ys.mina, labelcol = sig.type)
      }
      
      ol$dat <- list(mean_SE = res.medl,
                     P.value = pval.all,
                     P.value_pairwise = pval.pair,
                     Rawdata = dplyr::select(dfi,-gp2))
   };#  ol$p
   # 输出PDF return ####
   if(is.null(opre)){
      return(ol)
   }else{
      dir.create(dirname(opre),showWarnings = F,recursive = T)
      ggsave(paste0(opre,'.pdf'), ol$p, device = cairo_pdf, width = pw, height = ph)
      writexl::write_xlsx(ol$dat, paste0(opre, '.xlsx'))
      message(opre,'pdf xlsx \n')
   }
}


