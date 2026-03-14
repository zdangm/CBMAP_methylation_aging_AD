library(meta)
library(rhdf5)
library(data.table)
library(dplyr)
library(future.apply)
library(tibble)
library(bacon)

results.files <- c("Hannum_meth_age.RData","PEG1_meth_age.RData","SATSA_meth_age.RData")
for(i in 1:length(results.files)){
  load(results.files[i])
  data <- res[,1:5]
  colnames(data) <- c("cpg","Estimate", "StdErr", "pValue","fdr")
  dataset <- unlist(stringr::str_split(results.files[i],"_"))[1]  %>% as.character()
  aux <- paste0(dataset,c("_estimate", "_se", "_pValue", "_fdr"))
  colnames(data) <- c("cpg", aux)
  assign(dataset,data)
}

## Create a merged final dataset 
cohort_ls <- list(Hannum = Hannum, PEG1 = PEG1, SATSA=SATSA)

### outer join input region
multi_cohorts <- Reduce(
  function(x,y, ...) merge(x, y, by = "cpg", all = TRUE, ...),
  cohort_ls
) 

## Meta analysis 
### calculate meta analysis z scores and p values
#doParallel::registerDoParallel(cores = parallel::detectCores()/4)
plan(multisession, workers = 4)
multi_matrix <- as.matrix(multi_cohorts)
meta_list <- future_apply(multi_matrix, 1, function(rowOne) {
  row_df <- as.data.frame(t(rowOne), stringsAsFactors = FALSE)
  colnames(row_df) <- colnames(multi_cohorts)
  
  est <- as.numeric(row_df[grep("estimate", colnames(row_df))])
  se <- as.numeric(row_df[grep("se", colnames(row_df))])
  direction <- paste(ifelse(is.na(est), ".", ifelse(est > 0, "+", "-")), collapse = "")
  cohort <- gsub("_se", "", grep("se", colnames(row_df), value = TRUE))
  
  meta_input <- data.frame(
    cohort = cohort,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
  
  if (all(is.na(meta_input$est)) || all(is.na(meta_input$se))) {
    return(tibble(
      cpg = row_df$cpg,
      estimate = NA,
      se = NA,
      pVal.fixed = NA,
      pVal.random = NA,
      pValQ = NA,
      direction = direction
    ))
  }
  
  f <- tryCatch({
    metagen(TE = meta_input$est, seTE = meta_input$se, data = meta_input)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(f)) {
    return(tibble(
      cpg = row_df$cpg,
      estimate = NA,
      se = NA,
      pVal.fixed = NA,
      pVal.random = NA,
      pValQ = NA,
      direction = direction
    ))
  }
  
  tibble(
    cpg = row_df$cpg,
    estimate = f$TE.fixed,
    se = f$seTE.fixed,
    pVal.fixed = f$pval.fixed,
    pVal.random = f$pval.random,
    pValQ = f$pval.Q,
    direction = direction
  )
})
meta_df <- do.call(rbind, meta_list)
meta_df$pVal.final <- ifelse(meta_df$pValQ > 0.05, meta_df$pVal.fixed, meta_df$pVal.random)
meta_df$fdr <- p.adjust(meta_df$pVal.final, method = "BH")
meta_final_df <- meta_df[, c(grep("_",colnames(meta_df),invert = TRUE),
                             grep("_",colnames(meta_df),invert = FALSE))]
meta_final_ordered_df <- meta_final_df[order(meta_final_df$pVal.final),]
meta_final_ordered_df = merge(meta_final_ordered_df, multi_cohorts, by='cpg')
meta_final_ordered_df <- meta_final_ordered_df[order(meta_final_ordered_df$pVal.final),]
save(meta_final_ordered_df,file='meta_meth_age_blood.RData')

########## Q-Q plot ###########
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

res = meta_final_ordered_df[-which(is.na(meta_final_ordered_df$pVal.final)),]
set.seed(2025)
bc = bacon(effectsizes=res$estimate, standarderrors=res$se)
inflationFactor = inflation(bc)
biasFactor = bias(bc)
print(paste0('inflationFactor: ',inflationFactor))
print(paste0('biasFactor: ',biasFactor))

res$p_bacon = pval(bc)
res$p_bacon_fdr = p.adjust(res$p_bacon, 'BH')
print(sum(res$p_bacon_fdr < 0.05)) 
png('QQ_plot_meta_blood.png',width = 9,height = 4.5,unit='in',res=400)  
par(mfrow = c(1,2), mar = c(5, 5, 4, 2))
qq_plot(
  p    = res$pVal.final,
  main = "meta\nuncorrected",
  lim=c(0,80),
  add_text = c(
    paste0("Inflation = ", round(inflationFactor, 3))
  )
)
qq_plot(
  p    = res$p_bacon,
  lim=c(0,80),
  main = "meta\ncorrected"
)
dev.off()