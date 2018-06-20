rm(list=ls());
load('data/preprocessing/probes_common.RData');
load('data/preprocessing/normalized_data_probe_id_details.RData');
require(foreach);
common_probe_data<-foreach(i =1:length(normalized_data_probe_id_details))%do%{
  list(E=normalized_data_probe_id_details[[i]]$data$E[probes_common,],
       name=normalized_data_probe_id_details[[i]]$name)
};
save(common_probe_data,file='data/common_probe_data.RData');
labels_targets<-list()
combined_data<-foreach(i = 1:length(common_probe_data),.combine=cbind)%do%{
  labels_targets[[i]]<-rep(common_probe_data[[i]]$name,12)
  common_probe_data[[i]]$E
}
labels_targets<-unlist(labels_targets);
png('figures/MDS_plot_normalized_data_13012_probes.png',
    width = 900,height = 800,units = 'px',res = 150)
plotMDS(x = combined_data,
        labels = '*',
        cex=1.25,
        col=rep(c("red","grey40","blue"),each=12),
        main=paste('MDS plot with',length(probes_common), 'probes',sep=' '));
legend(x = 0.85,
       y = 1.25,
       pch = '*',
       cex=.75,
       legend = c('Celiac','Control','FDR'),       
       col = rep(c("red","grey40","blue")))
dev.off()
###################
### NORMALITY CHECK
### Normality test awas performed with Shapiro Wilk test and Anderson-Darling test
rm(list=ls())
require(nortest)
load('data/common_probe_data.RData')
require(doMC);
registerDoMC(3);
normality_test<-foreach(i = 1:length(common_probe_data))%dopar%{
  expr_data<-common_probe_data[[i]]$E
  p_value<-foreach(j = 1:nrow(expr_data))%do%{
    list(shapiro_test = shapiro.test(expr_data[j,])$p.value,
         aderson_darling_test= ad.test(expr_data[j,])$p.value)
  }
  list(p_value=p_value,
       name=common_probe_data[[i]]$name)
};
save(normality_test,file='data/normality_test.RData')
rm(list=ls())
load('data/normality_test.RData')
require(foreach);
shapiro_test<-list();
aderson_darling_test<-list();
foreach(i = 1:length(normality_test))%do%{
  ### Shapiro test
  p_value<-foreach(j = 1:length(normality_test[[i]]$p_value),.combine = c)%do%{
    normality_test[[i]]$p_value[[j]]$shapiro_test
  }
  shapiro_test[[i]]<-list(p_value=p_value,
                          name=normality_test[[i]]$name);
  ### Anderson test
  p_value<-foreach(j = 1:length(normality_test[[i]]$p_value),.combine = c)%do%{
    normality_test[[i]]$p_value[[j]]$aderson_darling_test
  }
  aderson_darling_test[[i]]<-list(p_value=p_value,
                                  name=normality_test[[i]]$name);
  print(normality_test[[i]]$name);
}
### segregate the normal and non-normal probe data based on shapiro test
load('data/common_probe_data.RData');
segregated_data_shipro_test<-list()
for(i in 1:length(common_probe_data)){
  if(common_probe_data[[i]]$name ==  shapiro_test[[i]]$name){
    print('ok');
    # alpha value is 0.01
    idx_non_normal<-which(shapiro_test[[i]]$p_value <= 0.1);
    idx_normal<-(-idx_non_normal)
    non_normal_data<-common_probe_data[[i]]$E[idx_non_normal,];
    normal_data<-common_probe_data[[i]]$E[idx_normal,];
    segregated_data_shipro_test[[i]]<-list(
      normal=normal_data,
      normal_pvalue=shapiro_test[[i]]$p_value[idx_normal],
      non_normal=non_normal_data,
      non_normal_pvalue=shapiro_test[[i]]$p_value[idx_non_normal],
      name=shapiro_test[[i]]$name
      );
  }
}
### segregate the normal and non-normal probe data based on Anderson Darling test
load('data/common_probe_data.RData');
segregated_data_anderson_darling_test<-list()
for(i in 1:length(common_probe_data)){
  if(common_probe_data[[i]]$name ==  aderson_darling_test[[i]]$name){
    print('ok');
    # alpha value is 0.01
    idx_non_normal<-which(aderson_darling_test[[i]]$p_value <= 0.1);
    idx_normal<-(-idx_non_normal)
    non_normal_data<-common_probe_data[[i]]$E[idx_non_normal,];
    normal_data<-common_probe_data[[i]]$E[idx_normal,];
    segregated_data_anderson_darling_test[[i]]<-list(
      normal=normal_data,
      normal_pvalue=aderson_darling_test[[i]]$p_value[idx_normal],
      non_normal=non_normal_data,
      non_normal_pvalue=aderson_darling_test[[i]]$p_value[idx_non_normal],
      name=shapiro_test[[i]]$name);
  }
}
#########
### Find total non normal probes across three data set
non_normal_probes<-list();
for(i in 1:length(segregated_data_shipro_test)){
  non_normal_probes[[i]]<-list(shapiro_test=rownames(segregated_data_shipro_test[[i]]$non_normal),
                               anderson_darling_test=rownames(segregated_data_anderson_darling_test[[i]]$non_normal))
  
}
non_normal_probes<-unique(unlist(non_normal_probes));
save(non_normal_probes,file='data/non_normal_probes.RData');
rm(list=ls())
load('data/non_normal_probes.RData');
load('data/common_probe_data.RData');
normal_probes<-setdiff(rownames(common_probe_data[[1]]$E),non_normal_probes);
common_data_distribution_based<-list();
for(i in 1:length(common_probe_data)){
  common_data_distribution_based[[i]]<-list(normal=common_probe_data[[i]]$E[normal_probes,],
                                            non_normal=common_probe_data[[i]]$E[non_normal_probes,],
                                            name=common_probe_data[[i]]$name)
}
for(i in 1:length(common_data_distribution_based)){
  print(common_data_distribution_based[[i]]$name)
}
require(foreach);
### Normal probe data
normal_probe_df<-foreach(i = 1:length(common_data_distribution_based), .combine=cbind)%do%{
  common_data_distribution_based[[i]]$normal
}
normal_probe_data<-list(df=normal_probe_df,
                        targets=rep(c('celiac','control','FDR'),each=12));
save(normal_probe_data,file='data/normal_probe_data.RData')
### Non-normal probe data
non_normal_probe_df<-foreach(i = 1:length(common_data_distribution_based), .combine=cbind)%do%{
  common_data_distribution_based[[i]]$non_normal
}
non_normal_probe_data<-list(df=non_normal_probe_df,
                            targets=rep(c('celiac','control','FDR'),each=12));
save(non_normal_probe_data,file='data/non_normal_probe_data.Rdata')
