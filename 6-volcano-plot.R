rm(list=ls())
load('./data/pair_wise__p_value_fc__t_test.RData')
pw.data <- pair_wise__p_value_fc__t_test

CelVsCtrl <- data.frame(bon.p.val = p.adjust(pw.data$celiac_vs_control.pval,
                                             method = "bonferroni"),
                        log2fc = pw.data$celiac_vs_control.fc,
                        Probe = rownames(pw.data),
                        comparison="CeD/Control")

CelVsFDR <- data.frame(bon.p.val = p.adjust(pw.data$celiac_vs_FDR.pval,method = "bonferroni"),
                       log2fc = pw.data$celiac_vs_FDR.fc,
                       Probe = rownames(pw.data),
                       comparison="CeD/FDR")
FDRVsCtrl <- data.frame(bon.p.val = p.adjust(pw.data$FDR_vs_control.pval,method = "bonferroni"),
                        log2fc = pw.data$FDR_vs_control.fc,
                        Probe = rownames(pw.data),
                        comparison="FDR/Control")
df.plot <- rbind(CelVsFDR,CelVsCtrl,FDRVsCtrl)
rm(CelVsFDR)
rm(CelVsCtrl)
rm(FDRVsCtrl)
getUpDownReg <- function(p.val,fc){
  p.val_passed <- p.val <= 0.05
  fc.up <- fc >= 1
  fc.down <- fc <= -1
  status <- rep('no',length(p.val))
  status[which(p.val_passed & fc.up)] <- 'up'
  status[which(p.val_passed & fc.down)] <- 'down'
  return(status)
}
status <- getUpDownReg(p.val = df.plot$bon.p.val,
                       fc = df.plot$log2fc)
df.plot<- data.frame(df.plot,
                     diff.expr=status)
library(ggplot2)
p <- ggplot(df.plot,aes(x=log2fc, y= -log10(bon.p.val))) +
  geom_point(aes(col=diff.expr),alpha=0.5,shape=16) +
  scale_color_manual(name="Differential\nregulation",
                     labels=c("Down",
                              "No",
                              "Up"),
                     values=c('#009e73ff',
                              '#333333ff',
                              '#d55e00ff')) +
  facet_wrap(facets="comparison") +
  xlab(expression(paste(log[2],"(fold change)"))) +
  ylab(expression(paste(-log[10],"(adjusted p-value)"))) +
  theme(axis.text.x=element_text(angle = 90,
                                 vjust = 0.5,
                                 hjust = 1)
        )
pdf('./figures/2-volcano-plot.pdf',width = 5.5,height = 4)
print(p)
dev.off()