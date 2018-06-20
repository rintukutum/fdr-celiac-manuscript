############
### ANOVA for nomally distributed probes
rm(list=ls())
load('data/normal_probe_data.RData');
# ## Test drive
# probe_data<-normal_probe_data$df[1,];
# targets<-as.factor(normal_probe_data$targets);
# aov_test<-aov(probe_data~targets)
get_pvalue_aov_test<-function(aov_test){
  summary_values<-summary(aov_test);
  p_value<-summary_values[[1]]['targets',5];
  return(p_value)
}
require(foreach)
aov_results_pvalue<-foreach(i =1:nrow(normal_probe_data$df),.combine=c)%do%{
  probe_data<-normal_probe_data$df[i,];
  targets<-as.factor(normal_probe_data$targets);
  aov_test<-aov(probe_data~targets)
  get_pvalue_aov_test(aov_test = aov_test)
}
names(aov_results_pvalue)<-rownames(normal_probe_data$df);
save(aov_results_pvalue,file = 'data/aov_results_pvalue.Rdata')
rm(list=ls())
load('data/aov_results_pvalue.Rdata');
load('data/normal_probe_data.RData')
idx_sig<-which(aov_results_pvalue <= 0.01);
aov_probe_sig<-names(aov_results_pvalue)[idx_sig];
aov_probe_sig_data<-list(df=normal_probe_data$df[aov_probe_sig,],
                         targets=normal_probe_data$targets)
save(aov_probe_sig_data,file='data/aov_probe_sig_data.Rdata');
### Kruskal-Wallis Rank Sum Test for non-normal distributed probes
rm(list=ls())
load('data/non_normal_probe_data.Rdata');
# ## Test dirve
# probe_data<-non_normal_probe_data$df[1,];
# targets<-as.factor(non_normal_probe_data$targets)
# kruskal_wallis_rank_sum_test<-kruskal.test(x = probe_data,g = targets);
# p_value<-kruskal_wallis_rank_sum_test$p.value;
targets<-as.factor(non_normal_probe_data$targets)
require(foreach)
kruskal_results<-foreach(i =1:nrow(non_normal_probe_data$df), .combine=c)%do%{
  probe_data<-non_normal_probe_data$df[i,];
  kruskal.test(x = probe_data,g = targets)$p.value;
}
names(kruskal_results)<-rownames(non_normal_probe_data$df)
save(kruskal_results,file='data/kruskal_results.RData')
rm(list=ls());
load('data/kruskal_results.RData');
load('data/non_normal_probe_data.Rdata')
idx_sig<-which(kruskal_results <= 0.01);
kruskal_probe_sig<-names(kruskal_results)[idx_sig];
kruskal_probe_sig_data<-list(df=non_normal_probe_data$df[kruskal_probe_sig,],
                             targets=non_normal_probe_data$targets);
boxplot(kruskal_probe_sig_data$df[109,]~kruskal_probe_sig_data$targets)
save(kruskal_probe_sig_data,file='data/kruskal_probe_sig_data.RData')
