#######################################################################
### Pair-wise test for normal probes
rm(list=ls());
load('data/aov_probe_sig_data.Rdata');
require(foreach);
idx_celiac<-grep('^celiac$',aov_probe_sig_data$targets);
idx_control<-grep('^control$',aov_probe_sig_data$targets);
idx_FDR<-grep('^FDR$',aov_probe_sig_data$targets);
###
idx_ce_co<-c(idx_celiac,idx_control);
idx_ce_fdr<-c(idx_celiac,idx_FDR);
idx_co_fdr<-c(idx_control,idx_FDR)
pair_wise__p_value_fc__t_test<-foreach(i = 1:nrow(aov_probe_sig_data$df), .combine=rbind)%do%{
  probe_data<-aov_probe_sig_data$df[i,]
  ## Celiac vs control
  df_test<-data.frame(value=as.numeric(probe_data[idx_ce_co]),
                      targets=aov_probe_sig_data$targets[idx_ce_co],
                      stringsAsFactors = TRUE);
  fc__ce_co<-mean(df_test[df_test$targets == 'celiac',]$value) - mean(df_test[df_test$targets == 'control',]$value);
  p_value_ce_co<-t.test(value~targets,data=df_test,paired = FALSE)$p.value;
  ## Celiac vs FDR
  df_test<-data.frame(value=as.numeric(probe_data[idx_ce_fdr]),
                      targets=aov_probe_sig_data$targets[idx_ce_fdr],
                      stringsAsFactors = TRUE);
  fc__ce_fdr<-mean(df_test[df_test$targets == 'celiac',]$value) - mean(df_test[df_test$targets == 'FDR',]$value);
  p_value_ce_fdr<-t.test(value~targets,data=df_test,paired = FALSE)$p.value;
  ## Control vs FDR
  df_test<-data.frame(value=as.numeric(probe_data[idx_co_fdr]),
                      targets=aov_probe_sig_data$targets[idx_co_fdr],
                      stringsAsFactors = TRUE);
  p_value_fdr_co<-t.test(value~targets,data=df_test,paired = FALSE)$p.value;
  fc__fdr_co<-mean(df_test[df_test$targets == 'FDR',]$value) - mean(df_test[df_test$targets == 'control',]$value);
  p_val__n__fc<-data.frame('celiac_vs_control.pval'=p_value_ce_co,
                           'celiac_vs_control.fc'=fc__ce_co,
                           'celiac_vs_FDR.pval'=p_value_ce_fdr,
                           'celiac_vs_FDR.fc'=fc__ce_fdr,
                           'FDR_vs_control.pval'=p_value_fdr_co,
                           'FDR_vs_control.fc'=fc__fdr_co);
}
rownames(pair_wise__p_value_fc__t_test)<-rownames(aov_probe_sig_data$df)
save(pair_wise__p_value_fc__t_test,file="data/pair_wise__p_value_fc__t_test.RData")
##################
rm(list=ls())
load('./data/kruskal_probe_sig_data.RData')
require(foreach);
idx_celiac<-grep('^celiac$',kruskal_probe_sig_data$targets);
idx_control<-grep('^control$',kruskal_probe_sig_data$targets);
idx_FDR<-grep('^FDR$',kruskal_probe_sig_data$targets);
###
idx_ce_co<-c(idx_celiac,idx_control);
idx_ce_fdr<-c(idx_celiac,idx_FDR);
idx_co_fdr<-c(idx_control,idx_FDR)
pair_wise__p_value_fc__wilcox_test<-foreach(i = 1:nrow(kruskal_probe_sig_data$df), .combine=rbind)%do%{
  probe_data<-kruskal_probe_sig_data$df[i,]
  ## Celiac vs Celiac control
  df_test<-data.frame(value=as.numeric(probe_data[idx_ce_co]),
                      targets=kruskal_probe_sig_data$targets[idx_ce_co],
                      stringsAsFactors = TRUE);
  fc__ce_co<-median(df_test[df_test$targets == 'celiac',]$value) - median(df_test[df_test$targets == 'control',]$value);
  p_value_ce_co<-wilcox.test(value~targets,data=df_test,exact = FALSE)$p.value;
  ## Celiac vs Celiac family
  df_test<-data.frame(value=as.numeric(probe_data[idx_ce_fdr]),
                      targets=kruskal_probe_sig_data$targets[idx_ce_fdr],
                      stringsAsFactors = TRUE);
  fc__ce_fdr<-median(df_test[df_test$targets == 'celiac',]$value) - median(df_test[df_test$targets == 'FDR',]$value);
  p_value_ce_fdr<-wilcox.test(value~targets,data=df_test,exact = FALSE)$p.value;
  ## Celiac Control vs Celiac family
  df_test<-data.frame(value=as.numeric(probe_data[idx_co_fdr]),
                      targets=kruskal_probe_sig_data$targets[idx_co_fdr],
                      stringsAsFactors = TRUE);
  p_value_fdr_co<-wilcox.test(value~targets,data=df_test,exact = FALSE)$p.value;
  fc__fdr_co<-median(df_test[df_test$targets == 'FDR',]$value) - median(df_test[df_test$targets == 'control',]$value);
  p_val__n__fc<-data.frame('celiac_vs_control.pval'=p_value_ce_co,
                           'celiac_vs_control.fc'=fc__ce_co,
                           'celiac_vs_FDR.pval'=p_value_ce_fdr,
                           'celiac_vs_FDR.fc'=fc__ce_fdr,
                           'FDR_vs_control.pval'=p_value_fdr_co,
                           'FDR_vs_control.fc'=fc__fdr_co);
};
rownames(pair_wise__p_value_fc__wilcox_test)<-rownames(kruskal_probe_sig_data$df);
save(pair_wise__p_value_fc__wilcox_test,file="data/pair_wise__p_value_fc__wilcox_test.RData")
rm(list=ls())
load('./data/pair_wise__p_value_fc__t_test.RData');
load('./data/pair_wise__p_value_fc__wilcox_test.RData');
length(which(p.adjust(pair_wise__p_value_fc__t_test$celiac_vs_control.pval,method = "bonferroni") <= 0.01))
#[1] 1374
length(which(p.adjust(pair_wise__p_value_fc__t_test$celiac_vs_FDR.pval,method = "bonferroni") <= 0.01))
#[1] 2623
length(which(p.adjust(pair_wise__p_value_fc__t_test$FDR_vs_control.pval,method = "bonferroni") <= 0.01))
#[1] 2514
length(which(p.adjust(pair_wise__p_value_fc__wilcox_test$celiac_vs_FDR.pval.pval,method = "bonferroni") <= 0.01))
#[1] 0
length(which(p.adjust(pair_wise__p_value_fc__wilcox_test$FDR_vs_control.pval,method = "bonferroni") <= 0.01))
#[1] 0
length(which(p.adjust(pair_wise__p_value_fc__wilcox_test$celiac_vs_control.pval,method = "bonferroni") <= 0.01))
#[1] 0
get_p_adjust<-function(x){
  x.colnames<-colnames(x)[grep(".pval",colnames(x))];
  x.p.adjust<-apply(x[,x.colnames],2,function(y){p.adjust(y,method = 'bonferroni')});
  colnames(x.p.adjust)<-paste(x.colnames,'_adjusted',sep='');
  return(data.frame(x,x.p.adjust))
}
pair_wise__p_value_fc__t_test_adjust<-get_p_adjust(x = pair_wise__p_value_fc__t_test)
fc_pval_pair_wise_test <- pair_wise__p_value_fc__t_test_adjust
save(fc_pval_pair_wise_test,
     file = './data/fc_pval_pair_wise_test.RData')

