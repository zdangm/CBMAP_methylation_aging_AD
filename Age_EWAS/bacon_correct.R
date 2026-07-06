library(bacon)
library(stringr)
library(ggplot2)

qq_plot <- function(p, lim, main = "", add_text = NULL) {
  p <- p[!is.na(p) & p > 0 & p <= 1]
  
  obs <- -log10(sort(p))
  exp <- -log10(ppoints(length(obs)))
  
  plot(
    exp, obs,
    pch = 16, cex = .5,
    ylim = lim,
    xlab = expression(Expected~~-log[10](P)),
    ylab = expression(Observed~~-log[10](P)),
    main = main
  )
  
  abline(0, 1, col = "red", lwd = 2)
  
  if (!is.null(add_text)) {
    legend(
      "topleft",
      legend = add_text,
      bty = "n"
    )
  }
}

# blood 
files.path = c('/methylation/methy_AD_revision2/results/EWAS/Hannum_meth_age.RData',
               '/methylation/methy_AD_revision2/results/EWAS/Hannum_meth_age_60years.RData',
               '/methylation/methy_AD_revision2/results/EWAS/SATSA_meth_age.RData',
               '/methylation/methy_AD_revision2/results/EWAS/SATSA_meth_age_60years.RData',
               '/methylation/methy_AD_revision2/results/EWAS/PEG1_meth_age.RData',
               '/methylation/methy_AD_revision2/results/EWAS/PEG1_meth_age_60years.RData')
title = c('Hannum (all samples)','Hannum (age > 60)','SATSA (all samples)','SATSA (age > 60)','PEG1 (all samples)','PEG1 (age > 60)')
lim_list = list(c(0,120),c(0,30),c(0,50),c(0,35),c(0,60),c(0,20))
study='Blood'
png(paste0('/methylation/methy_AD_revision2/results/QQ_plot/QQ_plot_',study,'.png'),width = 12,height = 9,unit='in',res=400)  
par(mfrow = c(3, 4), mar = c(5, 5, 4, 2))
for (i in 1:6) {
  load(files.path[i])
  print(sum(res$pValue < 0.05))
  print(sum(res$fdr < 0.05))

  set.seed(2025)
  bc = bacon(effectsizes=res$Estimate, standarderrors=res$StdErr)
  inflationFactor = inflation(bc)
  biasFactor = bias(bc)
  print(paste0('inflationFactor: ',inflationFactor))
  print(paste0('biasFactor: ',biasFactor))

  res$p_bacon = pval(bc)
  res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
  post_tstat = tstat(bc)
  set.seed(2025)
  post_bc = bacon(teststatistics = post_tstat)
  post_inflationFactor = inflation(post_bc)
  print(paste0('Post inflationFactor: ',post_inflationFactor))

  qq_plot(
    p    = res$pValue,
    main = paste0(title[i],"\nuncorrected"),
    lim=lim_list[[i]],
    add_text = c(
      paste0("Inflation = ", round(inflationFactor, 3))
    )
  )

  qq_plot(
    p    = res$p_bacon,
    lim=lim_list[[i]],
    main = paste0(title[i],"\ncorrected"),
    add_text = c(
      paste0("Inflation = ", round(post_inflationFactor, 3))
  ))

  save(res, file=files.path[i])
}
dev.off()

# blood meta age EWAS
files.path = c('/methylation/methy_AD_revision2/results/EWAS/meta_meth_age_blood_60years.RData')
load(files.path)
res=na.omit(meta_final_ordered_df)
title = c('meta','meta')
lim_list = c(0,40)
study='Blood'
png(paste0('/methylation/methy_AD_revision2/results/QQ_plot/QQ_plot_meta',study,'.png'),width = 8,height = 4,unit='in',res=400)  
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

print(sum(res$pVal.final < 0.05))
print(sum(res$fdr < 0.05))

set.seed(2025)
bc = bacon(effectsizes=res$estimate, standarderrors=res$se)
inflationFactor = inflation(bc)
biasFactor = bias(bc)
print(paste0('inflationFactor: ',inflationFactor))
print(paste0('biasFactor: ',biasFactor))

res$p_bacon = pval(bc)
res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
post_tstat = tstat(bc)
set.seed(2025)
post_bc = bacon(teststatistics = post_tstat)
post_inflationFactor = inflation(post_bc)
print(paste0('Post inflationFactor: ',post_inflationFactor))

qq_plot(
  p    = res$pVal.final,
  main = paste0(title[1],"\nuncorrected"),
  lim=lim_list,
  add_text = c(
    paste0("Inflation = ", round(inflationFactor, 3))
  )
)

qq_plot(
  p    = res$p_bacon,
  lim=lim_list,
  main = paste0(title[2],"\ncorrected"),
  add_text = c(
    paste0("Inflation = ", round(post_inflationFactor, 3))
  ))
dev.off()
save(res,file='/methylation/methy_AD_revision2/results/EWAS/meta_meth_age_blood_60years_bacon.RData')

# ROSMAP age EWAS
files.path = c('/methylation/methy_AD_revision2/results/EWAS/ROSMAP_meth_age.RData',
               '/methylation/methy_AD_revision2/results/EWAS/ROSMAP_meth_age_noDementia.RData')
title = c('ROSMAP (all samples)','ROSMAP (no dementia)')
lim_list = list(c(0,70),c(0,40))
study='ROSMAP'
png(paste0('/methylation/methy_AD_revision2/results/QQ_plot/QQ_plot_',study,'.png'),width = 7,height = 7,unit='in',res=400)  
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
for (i in 1:2) {
  load(files.path[i])
  print(sum(res$pValue < 0.05))
  print(sum(res$fdr < 0.05))
  
  set.seed(2025)
  bc = bacon(effectsizes=res$Estimate, standarderrors=res$StdErr)
  inflationFactor = inflation(bc)
  biasFactor = bias(bc)
  print(paste0('inflationFactor: ',inflationFactor))
  print(paste0('biasFactor: ',biasFactor))
  
  res$p_bacon = pval(bc)
  res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
  post_tstat = tstat(bc)
  set.seed(2025)
  post_bc = bacon(teststatistics = post_tstat)
  post_inflationFactor = inflation(post_bc)
  print(paste0('Post inflationFactor: ',post_inflationFactor))
  
  qq_plot(
    p    = res$pValue,
    main = paste0(title[i],"\nuncorrected"),
    lim=lim_list[[i]],
    add_text = c(
      paste0("Inflation = ", round(inflationFactor, 3))
    )
  )
  
  qq_plot(
    p    = res$p_bacon,
    lim=lim_list[[i]],
    main = paste0(title[i],"\ncorrected"),
    add_text = c(
      paste0("Inflation = ", round(post_inflationFactor, 3))
    ))
  
  save(res, file=files.path[i])
}
dev.off()

# CBMAP age EWAS
files.path = c('/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_age.RData',
               '/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_age_no_dementia_947.RData',
               '/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_age_noDementia_60years.RData')
title = c('CBMAP (all samples)','CBMAP (no dementia)','CBMAP (no dementia & age > 60)')
lim_list = list(c(0,160),c(0,160),c(0,70))
study='CBMAP'
png(paste0('/methylation/methy_AD_revision2/results/QQ_plot/QQ_plot_',study,'.png'),width = 6,height = 9,unit='in',res=400)  
par(mfrow = c(3, 2), mar = c(5, 5, 4, 2))
for (i in 1:3) {
  load(files.path[i])
  print(sum(res$pValue < 0.05))
  print(sum(res$fdr < 0.05))
  
  set.seed(2025)
  bc = bacon(effectsizes=res$Estimate, standarderrors=res$StdErr)
  inflationFactor = inflation(bc)
  biasFactor = bias(bc)
  print(paste0('inflationFactor: ',inflationFactor))
  print(paste0('biasFactor: ',biasFactor))
  
  res$p_bacon = pval(bc)
  res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
  post_tstat = tstat(bc)
  set.seed(2025)
  post_bc = bacon(teststatistics = post_tstat)
  post_inflationFactor = inflation(post_bc)
  print(paste0('Post inflationFactor: ',post_inflationFactor))
  
  qq_plot(
    p    = res$pValue,
    main = paste0(title[i],"\nuncorrected"),
    lim=lim_list[[i]],
    add_text = c(
      paste0("Inflation = ", round(inflationFactor, 3))
    )
  )
  
  qq_plot(
    p    = res$p_bacon,
    lim=lim_list[[i]],
    main = paste0(title[i],"\ncorrected"),
    add_text = c(
      paste0("Inflation = ", round(post_inflationFactor, 3))
    ))
  
  save(res, file=files.path[i])
}
dev.off()

# AD EWAS
files.path = c('/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_hml4level.RData',
               '/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_braak.RData',
               '/methylation/methy_AD_revision2/results/EWAS/CBMAP_meth_Ascore.RData',
               '/methylation/methy_AD_revision2/results/EWAS/ROSMAP_meth_hml4level.RData',
               '/methylation/methy_AD_revision2/results/EWAS/ROSMAP_meth_braak.RData',
               '/methylation/methy_AD_revision2/results/EWAS/ROSMAP_meth_Ascore.RData')
lim_list = list(c(0,8),c(0,11))
png(paste0('/methylation/methy_AD_revision2/results/QQ_plot/QQ_plot_AD_EWAS.png'),width = 9,height = 6,unit='in',res=400)  
par(mfrow = c(2,3), mar = c(5, 5, 4, 2))
for (i in 1:6) {
  load(files.path[i])
  res = as.data.frame(t(results))
  colnames(res) = c('beta','se','t','p')
  res$fdr = p.adjust(res$p,'BH')
  if (i %in% c(1,2,3)) {study = 'CBMAP'; lim = lim_list[[1]]} else {study = 'ROSMAP'; lim = lim_list[[2]]}
  if (i %in% c(1,4)) {group = 'HML'} else if (i %in% c(2,5)) {group = 'Braak'} else {group = 'Ascore'}
  print(sum(res$p < 0.05))
  print(sum(res$fdr < 0.05))
  
  set.seed(2025)
  bc = bacon(effectsizes=res$beta, standarderrors=res$se)
  inflationFactor = inflation(bc)
  biasFactor = bias(bc)
  print(paste0('inflationFactor: ',inflationFactor))
  print(paste0('biasFactor: ',biasFactor))
  
  res$p_bacon = pval(bc)
  res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
  print(sum(res$p_bacon_fdr < 0.05))
  
  qq_plot(
    p    = res$p,
    main = paste0(study," ",group," Uncorrected"),
    lim=lim,
    add_text = c(
      paste0("Inflation = ", round(inflationFactor, 3))
    )
  )
}
dev.off()
