library(meta)
library(rhdf5)
library(data.table)
library(dplyr)
library(future.apply)
library(tibble)

results.files <- c("ROSMAP_meth_hml4level.RData","CBMAP_meth_hml4level.RData")
for(i in 1:length(results.files)){
  load(results.files[i])
  data <- as.data.frame(t(results[c(1,2,4),]),stringsAsFactors = FALSE)
  colnames(data) <- c("Estimate", "StdErr", "pValue")
  data$cpg <- colnames(results)
  data <- data[,c("cpg", "Estimate", "StdErr", "pValue")]
  data$fdr <- p.adjust(data$pValue,'BH')
  dataset <- unlist(stringr::str_split(results.files[i],"_"))[1]  %>% as.character()
  aux <- paste0(dataset,c("_estimate", "_se", "_pValue", "_fdr"))
  colnames(data) <- c("cpg", aux)
  assign(dataset,data)
}

## Create a merged final dataset 
cohort_ls <- list(CBMAP = CBMAP, ROSMAP = ROSMAP)

### outer join input region
multi_cohorts <- Reduce(
  function(x,y, ...) merge(x, y, by = "cpg", all = TRUE, ...),
  cohort_ls
) 

## Meta analysis 
### calculate meta analysis z scores and p values
#doParallel::registerDoParallel(cores = parallel::detectCores()/4)
plan(multisession, workers = 30)
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
save(meta_final_ordered_df,file='meta_hml4level2.RData')