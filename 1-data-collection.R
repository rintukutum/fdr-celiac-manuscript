rm(list=ls());
#### Normalization and filtering
require(limma);
dir.create('data/preprocessing',showWarnings = FALSE,recursive = TRUE);
dir.create('figures',showWarnings = FALSE);

celiac.dat<-read.ilmn(
  files = "./data/raw/Celiac_without_background.txt",
  ctrlfiles = "./data/raw/Celiac_control_probes_without_background.txt",
  other.columns="Detection");
#### Background corrected
back_celiac.dat<-nec(celiac.dat)
#### Quantile normalization
norm_celiac<-normalizeBetweenArrays(object = back_celiac.dat,method = 'quantile');
celiac_data<-list(
  raw=celiac.dat, # raw data
  background_corrrected=back_celiac.dat, # background corrected
  normalized_quantile=norm_celiac
  ); # normalized using quantile 
save(celiac_data,file='data/preprocessing/celiac_data.RData')
rm(list=ls())
#### CELIAC CONTROLS
celiac_control.dat<-read.ilmn(
  files = './data/raw/CeliacControl_without background.txt',
  ctrlfiles = './data/raw/CeliacControl_control_probes_without_background.txt',
  other.columns="Detection"
);
#### Background corrected
back_celiac_control.dat<-nec(celiac_control.dat);
#### Quantile normalization
norm_celiac_control<-normalizeBetweenArrays(
  object = celiac_control.dat,
  method = 'quantile'
  );
celiac_control_data<-list(
  raw=celiac_control.dat,
  background_corrected=back_celiac_control.dat,
  normalized_quantile=norm_celiac_control);
save(celiac_control_data,file='data/preprocessing/celiac_control_data.RData')
rm(list=ls())
#### CELIAC FAMILY or FDR
celiac_family.dat<-read.ilmn(
  files = './data/raw/Family_without_background.txt',
  ctrlfiles = './data/raw/Family_control_probes_without_background.txt',
  other.columns="Detection");
#### Background correction
back_celiac_family.dat<-nec(celiac_family.dat);
#### Quantile normalization
norm_celiac_family<-normalizeBetweenArrays(object = celiac_family.dat,method = 'quantile');
celiac_family_data<-list(
  raw=celiac_family.dat,
  background_corrected=back_celiac_family.dat,
  normalized_quantile=norm_celiac_family);
save(celiac_family_data,
     file='data/preprocessing/celiac_family_data.RData')
rm(list=ls());
load('data/preprocessing/celiac_data.RData');
load('data/preprocessing/celiac_control_data.RData');
load('data/preprocessing/celiac_family_data.RData')
########################### get common probe data
get_expressed_normalized_data<-function(expression_data){
  norm_expr<-expression_data$normalized_quantile;
  expressed<-which(norm_expr$genes$Status == 'regular');
  norm_expr_expressed<-norm_expr[expressed,];
  probe_detect<-which(rowSums(norm_expr_expressed$other$Detection < 0.05) == 12);
  return(norm_expr_expressed[probe_detect,]);
}

### Celiac
celiac_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_data
);
### Celiac control / Control
celiac_control_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_control_data
);
### Celiac family / FDR
celiac_family_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_family_data
);
### Celiac
celiac_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_data
);
### Celiac control / Control
celiac_control_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_control_data
);
### Celiac family / FDR
celiac_family_expr_data<-get_expressed_normalized_data(
  expression_data = celiac_family_data
);
### common PROBES
normalized_data_probe_id_details<-list(
  list(
    probe_ids=dimnames(celiac_expr_data$E)[[1]],
    chip_ids=dimnames(celiac_expr_data$E)[[2]],
    name='celiac',
    normalization_method='quantile',
    data=celiac_expr_data),
  list(
    probe_ids=dimnames(celiac_control_expr_data$E)[[1]],
    chip_ids=dimnames(celiac_control_expr_data$E)[[2]],
    name='control',
    normalization_method='quantile',
    data=celiac_control_expr_data),
  list(
    probe_ids=dimnames(celiac_family_expr_data)[[1]],
    chip_ids=dimnames(celiac_family_expr_data)[[2]],
    name='FDR',
    normalization_method='quantile',
    data=celiac_family_expr_data)
);
save(normalized_data_probe_id_details,
     file='data/preprocessing/normalized_data_probe_id_details.RData'
)
###
all_probe_ids<-unique(c(
  normalized_data_probe_id_details[[1]]$probe_ids,
  normalized_data_probe_id_details[[2]]$probe_ids,
  normalized_data_probe_id_details[[3]]$probe_ids)
);
probes_common<-all_probe_ids
for(i in 1:length(normalized_data_probe_id_details)){
  probes_common<-intersect(probes_common,normalized_data_probe_id_details[[i]]$probe_ids);
}
save(probes_common,
     file='data/preprocessing/probes_common.RData')

