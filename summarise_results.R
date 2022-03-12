library(tidyverse)
library(lubridate)
library(knitr)
library(kableExtra)

load_results_bml<-function(filename,set,k,len,dist) {
  if (k==1) {
    res <- readRDS(filename) %>%
      mutate(bml_bulkess=ifelse(bml_divergent==0&(bml_d_bulk_ess<400|bml_f_bulk_ess<400|bml_phi_bulk_ess<400),1,0),
             bml_rhat=ifelse(bml_divergent==0&(bml_d_rhat>1.1|bml_f_rhat>1.1|bml_phi_rhat>1.1),1,0)) %>% #bml_bulkess==0&
      select(bml_d,bml_f,bml_phi,bml_t,bml_divergent,bml_bulkess,bml_rhat) %>%
      rename(d1=bml_d,f1=bml_f,phi=bml_phi,time=bml_t,divergent=bml_divergent,bulkess=bml_bulkess,rhat=bml_rhat) %>%
      mutate(d2=NA,f2=NA)
  } else {
    res <- readRDS(filename) %>%
      mutate(bml_bulkess=ifelse(bml_divergent==0&(bml_d1_bulk_ess<400|bml_f1_bulk_ess<400|bml_d2_bulk_ess<400|bml_f2_bulk_ess<400|bml_phi_bulk_ess<400),1,0),
             bml_rhat=ifelse(bml_divergent==0&(bml_d1_rhat>1.1|bml_f1_rhat>1.1|bml_d2_rhat>1.1|bml_f2_rhat>1.1|bml_phi_rhat>1.1),1,0)) %>% #bml_bulkess==0&
      select(bml_d1,bml_d2,bml_f1,bml_f2,bml_phi,bml_t,bml_divergent,bml_bulkess,bml_rhat)  %>%
      rename(d1=bml_d1,f1=bml_f1,d2=bml_d2,f2=bml_f2,phi=bml_phi,time=bml_t,divergent=bml_divergent,bulkess=bml_bulkess,rhat=bml_rhat)
    
    if (mean(res$f1,na.rm=TRUE)>mean(res$f2,na.rm=TRUE)) { # then swap them
      res$temp_f1 <- res$f2
      res$temp_d1 <- res$d2
      res$f2 <- res$f1
      res$d2 <- res$d1
      res$f1 <- res$temp_f1
      res$d1 <- res$temp_d1
    }
  }
  res$k    <- k
  res$len  <- len
  res$dist <- dist
  res$set  <- set
  res$method <- 'Bayesian Exact'
  if (set=="Set 1") {
  res$true_d1  <- res$true_d2 <- 0.22
  res$true_f1  <- 0.2
  res$true_f2  <- 0.4
  res$true_phi <- 0.8
  } else {
    res$true_d1  <- 0.3
    res$true_d2 <- 0.20
    res$true_f1  <- 1/28
    res$true_f2  <- 1/7
    res$true_phi <- (-0.8)
  }
  res$sim  <- sprintf("%s %s k=%d, n=%d",dist,set,k,len)
  return(res)
}
load_results_bw<-function(filename,set,k,len,dist) {
  if (k==1) {
    res <- readRDS(filename) %>%
      mutate(bw_bulkess=ifelse(bw_divergent==0&(bw_d1_bulk_ess<400|bw_f1_bulk_ess<400|bw_phi_bulk_ess<400),1,0),
             bw_rhat=ifelse(bw_divergent==0&bw_bulkess==0&(bw_d1_rhat>1.1|bw_f1_rhat>1.1|bw_phi_rhat>1.1),1,0)) %>%
      select(bw_d1,bw_f1,bw_phi,bw_t,bw_divergent,bw_bulkess,bw_rhat) %>%
      rename(d1=bw_d1,f1=bw_f1,phi=bw_phi,time=bw_t,divergent=bw_divergent,bulkess=bw_bulkess,rhat=bw_rhat) %>%
      mutate(d2=NA,f2=NA)
  } else {
    res <- readRDS(filename) %>%
      mutate(bw_bulkess=ifelse(bw_divergent==0&(bw_d1_bulk_ess<400|bw_f1_bulk_ess<400|bw_d2_bulk_ess<400|bw_f2_bulk_ess<400|bw_phi_bulk_ess<400),1,0),
             bw_rhat=ifelse(bw_divergent==0&bw_bulkess==0&(bw_d1_rhat>1.1|bw_f1_rhat>1.1|bw_d2_rhat>1.1|bw_f2_rhat>1.1|bw_phi_rhat>1.1),1,0)) %>%
      select(bw_d1,bw_d2,bw_f1,bw_f2,bw_phi,bw_t,bw_divergent,bw_bulkess,bw_rhat) %>%
      rename(d1=bw_d1,f1=bw_f1,d2=bw_d2,f2=bw_f2,phi=bw_phi,time=bw_t,divergent=bw_divergent,bulkess=bw_bulkess,rhat=bw_rhat)
    
    if (mean(res$f1,na.rm=TRUE)>mean(res$f2,na.rm=TRUE)) { # then swap them
      res$temp_f1 <- res$f2
      res$temp_d1 <- res$d2
      res$f2 <- res$f1
      res$d2 <- res$d1
      res$f1 <- res$temp_f1
      res$d1 <- res$temp_d1
      res$temp_f1 <- res$temp_d1 <- NULL
    }
  }
  res$k    <- k
  res$len  <- len
  res$dist <- dist
  res$set  <- set
  res$method <- 'Bayesian Whittle'
  if (set=="Set 1") {
    res$true_d1  <- res$true_d2 <- 0.22
    res$true_f1  <- 0.2
    res$true_f2  <- 0.4
    res$true_phi <- 0.8
  } else {
    res$true_d1  <- 0.30
    res$true_d2  <- 0.20
    res$true_f1  <- 1/28
    res$true_f2  <- 1/7
    res$true_phi <- (-0.8)
  }
  res$sim  <- sprintf("%s %s k=%d, n=%d",dist,set,k,len)
  return(res)
}

load_results<-function(filename,set,k,len,dist,method) {
  res <- readRDS(filename)
  col_list <- colnames(res)
  newnames <- c()
  for (col in col_list) {
    newname = strsplit(col,'_',fixed=T)[[1]][2]
    if (newname=='t') newname<-'time' #else if (newname!='phi'&k==1) newname <- paste0(newname,'2')
    newnames <- c(newnames,newname)
  }
  colnames(res) <- newnames
  if (k==1) {
    res$d2 <- res$f2 <- NA
  }
  res$divergent <- res$bulkess <- res$rhat <- 0
  res$k    <- k
  res$len  <- len
  res$dist <- dist
  res$set  <- set
  res$method <- method
  if (set=="Set 1") {
    res$true_d1  <- res$true_d2 <- 0.22
    res$true_f1  <- 0.2
    res$true_f2  <- 0.4
    res$true_phi <- 0.8
  } else {
    res$true_d1  <- 0.30
    res$true_d2 <-  0.20
    res$true_f1  <- 1/28
    res$true_f2  <- 1/7
    res$true_phi <- (-0.8)
  }
  res$sim  <- sprintf("%s %s k=%d, n=%d",dist,set,k,len)
  return(res)
}

res <- rbind(load_results_bml('~/gauss_k1_500/res_gauss_500_k1_bml.RDS',set='Set 2',k=1,len=500,dist='Gauss'),
             load_results_bw('~/gauss_k1_500/res_gauss_500_k1_bw.RDS',set='Set 2',k=1,len=500,dist='Gauss'),
             load_results('~/gauss_k1_500/res_gauss_500_k1_whittle.RDS',set='Set 2',k=1,len=500,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k1_500/res_gauss_500_k1_wll.RDS',set='Set 2',k=1,len=500,dist='Gauss',method='Log Whittle'),
             
             load_results_bml('~/gauss_k1_100/res_gauss_100_k1_bml.RDS',set='Set 2',k=1,len=100,dist='Gauss'),
             load_results_bw('~/gauss_k1_100/res_gauss_100_k1_bw.RDS',set='Set 2',k=1,len=100,dist='Gauss'),
             load_results('~/gauss_k1_100/res_gauss_100_k1_whittle.RDS',set='Set 2',k=1,len=100,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k1_100/res_gauss_100_k1_wll.RDS',set='Set 2',k=1,len=100,dist='Gauss',method='Log Whittle'),
             
             load_results_bml('~/gauss_k1_1000/res_gauss_1000_k1_bml.RDS',set='Set 2',k=1,len=1000,dist='Gauss'),
             load_results_bw('~/gauss_k1_1000/res_gauss_1000_k1_bw.RDS',set='Set 2',k=1,len=1000,dist='Gauss'),
             load_results('~/gauss_k1_1000/res_gauss_1000_k1_whittle.RDS',set='Set 2',k=1,len=1000,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k1_1000/res_gauss_1000_k1_wll.RDS',set='Set 2',k=1,len=1000,dist='Gauss',method='Log Whittle'),
             
             load_results_bml('~/gauss_k2_100/res_gauss_100_k2_bml.RDS',set='Set 2',k=2,len=100,dist='Gauss'),
             load_results_bw('~/gauss_k2_100/res_gauss_100_k2_bw.RDS',set='Set 2',k=2,len=100,dist='Gauss'),
             load_results('~/gauss_k2_100/res_gauss_100_k2_whittle.RDS',set='Set 2',k=2,len=100,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k2_100/res_gauss_100_k2_wll.RDS',set='Set 2',k=2,len=100,dist='Gauss',method='Log Whittle'),
             
             load_results_bml('~/gauss_k2_500/res_gauss_500_k2_bml.RDS',set='Set 2',k=2,len=500,dist='Gauss'),
             load_results_bw('~/gauss_k2_500/res_gauss_500_k2_bw.RDS',set='Set 2',k=2,len=500,dist='Gauss'),
             load_results('~/gauss_k2_500/res_gauss_500_k2_whittle.RDS',set='Set 2',k=2,len=500,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k2_500/res_gauss_500_k2_wll.RDS',set='Set 2',k=2,len=500,dist='Gauss',method='Log Whittle'),
             load_results_bml('~/gauss_k2_1000/res_gauss_1000_k2_bml.RDS',set='Set 2',k=2,len=1000,dist='Gauss'),
             load_results_bw('~/gauss_k2_1000/res_gauss_1000_k2_bw.RDS',set='Set 2',k=2,len=1000,dist='Gauss'),
             load_results('~/gauss_k2_1000/res_gauss_1000_k2_whittle.RDS',set='Set 2',k=2,len=1000,dist='Gauss',method='Whittle'),
             load_results('~/gauss_k2_1000/res_gauss_1000_k2_wll.RDS',set='Set 2',k=2,len=1000,dist='Gauss',method='Log Whittle'),
             
             load_results_bml('~/chisq_k1_100/res_chisq_100_k1_bml.RDS',set='Set 2',k=1,len=100,dist='ChiSq'),
             load_results_bw('~/chisq_k1_100/res_chisq_100_k1_bw.RDS',set='Set 2',k=1,len=100,dist='ChiSq'),
             load_results('~/chisq_k1_100/res_chisq_100_k1_whittle.RDS',set='Set 2',k=1,len=100,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k1_100/res_chisq_100_k1_wll.RDS',set='Set 2',k=1,len=100,dist='ChiSq',method='Log Whittle'),
             
             load_results_bml('~/chisq_k1_500/res_chisq_500_k1_bml.RDS',set='Set 2',k=1,len=500,dist='ChiSq'),
             load_results_bw('~/chisq_k1_500/res_chisq_500_k1_bw.RDS',set='Set 2',k=1,len=500,dist='ChiSq'),
             load_results('~/chisq_k1_500/res_chisq_500_k1_whittle.RDS',set='Set 2',k=1,len=500,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k1_500/res_chisq_500_k1_wll.RDS',set='Set 2',k=1,len=500,dist='ChiSq',method='Log Whittle'),
             
             load_results_bml('~/chisq_k1_1000/res_chisq_1000_k1_bml.RDS',set='Set 2',k=1,len=1000,dist='ChiSq'),
             load_results_bw('~/chisq_k1_1000/res_chisq_1000_k1_bw.RDS',set='Set 2',k=1,len=1000,dist='ChiSq'),
             load_results('~/chisq_k1_1000/res_chisq_1000_k1_whittle.RDS',set='Set 2',k=1,len=1000,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k1_1000/res_chisq_1000_k1_wll.RDS',set='Set 2',k=1,len=1000,dist='ChiSq',method='Log Whittle'),

             load_results_bml('~/chisq_k2_100/res_chisq_100_k2_bml.RDS',set='Set 2',k=2,len=100,dist='ChiSq'),
             load_results_bw('~/chisq_k2_100/res_chisq_100_k2_bw.RDS',set='Set 2',k=2,len=100,dist='ChiSq'),
             load_results('~/chisq_k2_100/res_chisq_100_k2_whittle.RDS',set='Set 2',k=2,len=100,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k2_100/res_chisq_100_k2_wll.RDS',set='Set 2',k=2,len=100,dist='ChiSq',method='Log Whittle'),
             
             load_results_bml('~/chisq_k2_500/res_chisq_500_k2_bml.RDS',set='Set 2',k=2,len=500,dist='ChiSq'),
             load_results_bw('~/chisq_k2_500/res_chisq_500_k2_bw.RDS',set='Set 2',k=2,len=500,dist='ChiSq'),
             load_results('~/chisq_k2_500/res_chisq_500_k2_whittle.RDS',set='Set 2',k=2,len=500,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k2_500/res_chisq_500_k2_wll.RDS',set='Set 2',k=2,len=500,dist='ChiSq',method='Log Whittle'),
             load_results_bml('~/chisq_k2_1000/res_chisq_1000_k2_bml.RDS',set='Set 2',k=2,len=1000,dist='ChiSq'),
             load_results_bw('~/chisq_k2_1000/res_chisq_1000_k2_bw.RDS',set='Set 2',k=2,len=1000,dist='ChiSq'),
             load_results('~/chisq_k2_1000/res_chisq_1000_k2_whittle.RDS',set='Set 2',k=2,len=1000,dist='ChiSq',method='Whittle'),
             load_results('~/chisq_k2_1000/res_chisq_1000_k2_wll.RDS',set='Set 2',k=2,len=1000,dist='ChiSq',method='Log Whittle')
)

res2 <- res %>% 
  mutate(d1=d1-true_d1,
         d2=d2-true_d2,
         f1=f1-true_f1,
         f2=f2-true_f2,
         phi=phi-true_phi,
         div_count = ifelse(divergent>0,1,0),
         bulkess = ifelse(bulkess>0,1,0),
         rhat = ifelse(rhat>0,1,0)) %>%
  group_by(dist,set,k,len,method) %>%
  summarise(d1_bias=mean(d1,na.rm=TRUE),
            d2_bias=mean(d2,na.rm=TRUE),
            f1_bias=mean(f1,na.rm=TRUE),
            f2_bias=mean(f2,na.rm=TRUE),
            phi_bias=mean(phi,na.rm=TRUE),
            d1_mse=mean(d1^2,na.rm=TRUE),
            d2_mse=mean(d2^2,na.rm=TRUE),
            f1_mse=mean(f1^2,na.rm=TRUE),
            f2_mse=mean(f2^2,na.rm=TRUE),
            phi_mse=mean(phi^2,na.rm=TRUE),
            time=mean(time),
            divergent=sum(div_count,na.rm=TRUE),
            bulkess = sum(bulkess,na.rm=TRUE),
            rhat = sum(rhat,na.rm=TRUE)) %>%
  mutate(d1=ifelse(is.na(d1_bias),"",sprintf("%.4f (%.4f)",d1_bias,d1_mse)),
         f1=ifelse(is.na(f1_bias),"",sprintf("%.4f (%.4f)",f1_bias,f1_mse)),
         d2=ifelse(is.na(d2_bias),"",sprintf("%.4f (%.4f)",d2_bias,d2_mse)),
         f2=ifelse(is.na(f2_bias),"",sprintf("%.4f (%.4f)",f2_bias,f2_mse)),
         phi=ifelse(is.na(phi_bias),"",sprintf("%.4f (%.4f)",phi_bias,phi_mse)),
         time=sprintf("%.2f",time)) %>%
  ungroup() %>%
  select(dist,set,k,len,method,d1,f1,d2,f2,phi,time,divergent,bulkess,rhat) %>%
  arrange(dist,set,k,len) 

res2 %>%
  filter(k==1,set=='Set 2') %>%
  select(-k,-set,-f2,-d2) %>%
  kbl(caption="Summary of simulations",
      format = "latex",
      escape=FALSE,
      align=c('c','c','l','r','r','r','r','r','r','r','r'),
      col.names=c('Distn','n','Method',"$d_1$","$f_1$","$\\phi$",'time','N. Div.','N Bulk ESS','N Rhat')
      ) -> tbl
save_kable(tbl,file='test_k1.tex')

res2 %>%
  filter(k==2) %>%
  select(-k) %>%
  kbl(caption="Summary of simulations",
      format = "latex",
      escape=FALSE,
      align=c('c','c','r','l','r','r','r','r','r','r','r','r','r'),
      col.names=c('Distn','Set','len','method',"$d_1$","$f_1$","$d_2$","$f_2$","$\\phi$",
                  'time','N. Div.','N Bulk ESS','N Rhat')
      ) -> tbl

save_kable(tbl,file='test_k2.tex')

plot_data <- res %>%
  mutate(d1=d1-true_d1,
         d2=d2-true_d2,
         f1=f1-true_f1,
         f2=f2-true_f2,
         phi=phi-true_phi,
         method=ifelse(method=='Bayesian Exact','Bayesian\nExact',
                       ifelse(method=='Bayesian Whittle','Bayesian\nWhittle',method))
         #dist=ifelse(sim=='Gauss Set 2 k=2, n=1000'&method=='Bayesian Exact','Gauss Partial',
         #            ifelse(sim=='ChiSq Set 2 k=2, n=1000'&method=='Bayesian Exact','ChiSq Partial',dist))
         ) %>%
  mutate(conv=ifelse(divergent>0,'Divergent',ifelse(rhat>0,'Rhat > 1.1','Converged'))) %>%
  dplyr::select(dist,len,k,method,conv,d1,f1,d2,f2,phi) %>%
  pivot_longer(c('d1','d2','f1','f2','phi'),names_to='param',values_to = 'estimate') %>%
  filter(!is.na(estimate)) %>%
  mutate(param=if_else(param=='d1','d[1]',
                       if_else(param=='d2','d[2]',
                               if_else(param=='f1','f[1]',
                                       if_else(param=='f2','f[2]','phi'))))) %>%
  mutate(kn = sprintf('k=%d, n=%4d',k,len))

plot_data<- rbind(plot_data,
                  # data.frame(dist=rep(c('ChiSq','Gauss'),5),
                  #            len=rep(1000,10),
                  #            k=rep(2,10),
                  #            method=rep('Bayesian\nExact',10),
                  #            conv=rep('Rhat > 1.1',10),
                  #            param=rep(c('phi','d[1]','d[2]','f[1]','f[2]'),2),
                  #            estimate=rep(0,10),
                  #            kn=rep('k=2, n=1000',10)),
                  data.frame(dist=rep('Gauss',10),
                             len=rep(1000,10),
                             k=rep(2,10),
                             method=rep('Bayesian\nExact',10),
                             conv=c(rep('ZZDivergentInvis',5),rep('ZZRhatInvis',5)),
                             param=rep(c('phi','d[1]','d[2]','f[1]','f[2]'),2),
                             estimate=rep(0,10),
                             kn=rep('k=2, n=1000',10))
                  )

dot_data <- filter(plot_data,conv %in% c('Divergent','ZZDivergentInvis','Rhat > 1.1','ZZRhatInvis'),k==2) %>%
            mutate(conv=factor(conv,levels=c('Divergent','Rhat > 1.1','ZZDivergentInvis','ZZRhatInvis')))
ggplot(filter(plot_data,k==2,conv=='Converged'),aes(x=method,y=estimate,fill=dist)) + #
  geom_hline(yintercept=0,color='black',size=0.2) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_boxplot(size=0.2,outlier.size=0.4,outlier.color='black')+
  geom_point(data=filter(dot_data,conv %in% c('Rhat > 1.1','ZZRhatInvis')),aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  geom_point(data=filter(dot_data,conv%in%c('Divergent','ZZDivergentInvis')),aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  scale_color_manual(name='',values=c('Divergent'='red','Rhat > 1.1'='cyan',
                                      'ZZRhatInvis'=NA,'ZZDivergentInvis'=NA),
                     labels=c('Divergent'='Divergent','Rhat > 1.1'='Rhat > 1.1',
                              'ZZRhatInvis'='','ZZDivergentInvis'=''),
                     guide=guide_legend(order=2))+
  scale_fill_manual(name='', values=c('ChiSq'="#3366FF", 'Gauss'='darkgray'), labels=c(expression( chi[1]^2),'Gauss'),guide=guide_legend(order=1)) +
  facet_grid(cols=vars(param),rows=vars(kn),labeller=labeller(param=label_parsed)) +
  labs(title='Distribution of bias in parameter estimates for k=2 simulations.',
       subtitle='Red dots indicate divergent transitions for Bayesian techniques',
       x='',
       y='Bias',
       caption=expression('True values '*f[1]*'=0.0357, '*f[2]*'= 0.1429, '*d[1]*'= 0.3, '*d[2]*'= 0.2, '*phi*'= -0.8'),
       parse=TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        strip.text = element_text(size=10),
        legend.text.align = 0,
        plot.title=element_text(color='black')
        ,plot.subtitle=element_text(color='red'))

ggsave('results_k2.jpeg',width = 11,height=8,dpi=600)

# plots for f1 and f2 separately
ggplot(filter(plot_data,k==2,conv=='Converged',param%in%c("f[1]","f[2]")),aes(x=method,y=estimate,fill=dist)) + #
  geom_hline(yintercept=0,color='black',size=0.2) +
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  geom_boxplot(size=0.2,outlier.size=0.4,outlier.color='black')+
  geom_point(data=filter(dot_data,conv %in% c('Rhat > 1.1','ZZRhatInvis'),param%in%c("f[1]","f[2]")),
             aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  geom_point(data=filter(dot_data,conv%in%c('Divergent','ZZDivergentInvis'),param%in%c("f[1]","f[2]")),
             aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  scale_color_manual(name='',values=c('Divergent'='red','Rhat > 1.1'='cyan',
                                      'ZZRhatInvis'=NA,'ZZDivergentInvis'=NA),
                     labels=c('Divergent'='Divergent','Rhat > 1.1'='Rhat > 1.1',
                              'ZZRhatInvis'='','ZZDivergentInvis'=''),
                     guide=guide_legend(order=2))+
  scale_fill_manual(name='', values=c('ChiSq'="#3366FF", 'Gauss'='darkgray'), labels=c(expression( chi[1]^2),'Gauss'),guide=guide_legend(order=1)) +
  facet_grid(cols=vars(param),rows=vars(kn),labeller=labeller(param=label_parsed)) +
  labs(title='Distribution of bias in parameter estimates for k=2 simulations.',
       subtitle='Red dots indicate divergent transitions for Bayesian techniques',
       x='',
       y='Bias',
       caption=expression('True values '*f[1]*'=0.0357, '*f[2]*'= 0.1429, '*d[1]*'= 0.3, '*d[2]*'= 0.2, '*phi*'= -0.8'),
       parse=TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        strip.text = element_text(size=10),
        legend.text.align = 0,
        plot.title=element_text(color='black')
        ,plot.subtitle=element_text(color='red'))

ggsave('results_k2_f12.jpeg',width = 8,height=8,dpi=600)

# plots for k=1
plot_data<- rbind(plot_data,
                  data.frame(dist=rep(c('ChiSq','Gauss'),3),
                             len=rep(500,6),
                             k=rep(1,6),
                             method=rep('Bayesian\nExact',6),
                             conv=rep('Rhat > 1.1',6),
                             param=rep(c('phi','d[1]','f[1]'),2),
                             estimate=rep(0,6),
                             kn=rep('k=1, n= 500',6)))

ggplot(filter(plot_data,k==1,conv=='Converged'),aes(x=method,y=estimate,fill=dist)) + #
  geom_hline(yintercept=0,color='black',size=0.2) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  geom_boxplot(size=0.2,outlier.size=0.4,outlier.color='black')+
  geom_point(data=filter(plot_data,conv=='Rhat > 1.1',k==1),aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  geom_point(data=filter(plot_data,conv=='Divergent',k==1),aes(color=conv),size=0.5,position=position_dodge(width=0.69)) +
  scale_fill_manual(name='', values=c('ChiSq'="#3366FF", 'Gauss'='darkgray'), labels=c(expression( chi[1]^2),'Gauss'),guide=guide_legend(order=1)) +
  scale_color_manual(name='',values=c('Divergent'='red','Rhat > 1.1'='cyan',' '=NA),guide=guide_legend(order=2))+
  facet_grid(cols=vars(param),rows=vars(kn),labeller=labeller(param=label_parsed)) +
  labs(title='Distribution of bias in parameter estimates for k=1 simulations.',
       subtitle='Red dots indicate divergent transitions in Bayesian techniques',
       x='',
       y='Bias',
       caption=expression('True values '*f[1]*'=0.0357, '*d[1]*'= 0.3, '*phi*'= -0.8'),
       parse=TRUE) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10,angle=90,hjust=1,vjust=0.5),
        strip.text = element_text(size=10),
        legend.text.align = 0,
        plot.title=element_text(color='black')
        ,plot.subtitle=element_text(color='red'))

ggsave('results_k1.jpeg',width = 11,height=8,dpi=600)

