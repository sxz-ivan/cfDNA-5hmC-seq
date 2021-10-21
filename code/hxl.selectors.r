suppressPackageStartupMessages({
  require(parallel)
  require(Rfast)
  require(caret)
  require(glmnet)
})
source('code/df.tools.r')
source('code/parallel.tools.r')
source('code/caret.pipelines.r')
source('code/sampling.r')
source('code/sxz.selectors.r')
source('code/matrix.stats.r')
hxl_glmnet = function(matrix, y, alpha=0.1, family="binomial", type.measure = "class",...,output.format = 'selected'){
  alpha.fit <- cv.glmnet(matrix, y, alpha=alpha, family=family, type.measure = type.measure,...)
  rm(matrix)
  gc()
  if(output.format == 'selected'){
    res <- as.matrix(coef(alpha.fit, s=alpha.fit$lambda.1se))
    return(rownames(res)[res[,1]!=0][-1])#remove intercept
  }else{
    return(alpha.fit)
  }
}
hxl_ensemble = function(matrix, y, n.repeats = 10, freq_cutoff = 0.95, selector = hxl_glmnet,mc.cores = 1,...){
  1:n.repeats %>% flexible.mclapply(mc.cores = mc.cores, FUN = function(i){
    selector(matrix, y,...)
  }) %>% unlist %>% table %>% as.data.frame(stringsAsFactors = F) %>% setNames(nm = c('ROWNAMES','Freq'))
}
hxlensemble_selector = function(matrix, y, n.repeats = 10, freq_cutoff = 0.95, selector = hxl_glmnet,mc.cores = 1,...){
  hxl_ensemble(matrix, y, n.repeats = n.repeats, freq_cutoff = freq_cutoff, selector = selector,mc.cores = mc.cores,...)%>%
    with({ROWNAMES[Freq> n.repeats*freq_cutoff]})
}
hxlfoldensemble_iteration = function(matrix,y,testFUN = row_t_welch,mc.cores = 1,significance_metric='fdr',logmat = F,top.cutoff=0.1,top.metric = c('abs_log2FC', 'significance'),top.method = row_L1_norm,...){
  stat = matrixtests_selector(matrix = matrix,y = y,FUN = testFUN,output.format = 'stats')%>%
    keep.rownames('ROWNAMES')
    if('abs_log2FC'%in%top.metric){stat=stat%>% within({
        if(logmat){
            abs_log2FC = abs(mean.x-mean.y)
        }else{
            abs_log2FC = abs(log2(mean.x/mean.y))
        }
    })}
  stat$significance = -log10(stat[[significance_metric]])
  selected_by_test = stat %>%
    top_stats(top = top.cutoff,metric = top.metric,method = top.method,output.format = 'stats')%>% rownames
  
  selected_by_ensemble = hxlensemble_selector(matrix = matrix[selected_by_test,] %>% t,
                                              y =y,n.repeats =10,mc.cores = mc.cores,...)
  
  selected_by_ensemble
}
hxlfoldensemble_selector = function(matrix, y, k.folds = 5, mc.cores = 1,testFUN = row_t_welch,significance_metric='fdr',logmat = F,
                                    top.cutoff=0.1,top.metric = c('abs_log2FC', 'significance'),top.method = row_L1_norm,...){
  folds = createFolds(y = y,k =k.folds,returnTrain = T)
  modeling_cores = max(floor(mc.cores/length(folds)),1)
  folds %>% flexible.mclapply(mc.cores =min(length(folds),mc.cores),FUN = function(fold){
    hxlfoldensemble_iteration(matrix = matrix[,fold],y = y[fold], FUN = testFUN,mc.cores = modeling_cores,significance_metric=significance_metric,logmat=logmat,top.top = top.cutoff,top.metric = top.metric,top.method = top.method,...)
  }) %>% unlist %>% table %>% as.data.frame(stringsAsFactors = F) %>% setNames(nm = c('ROWNAMES','Freq'))
}
# wc_L1_selector = function(matrix,y,testFUN = row_t_welch,significance_metric='pvalue',logmat = F,rank.method = row_L1_norm){
#   stat = matrixtests_selector(matrix = matrix,y = y,FUN = testFUN,output.format = 'stats')%>%
#     keep.rownames('ROWNAMES')%>% within({
#         if(logmat){
#             abs_log2FC = abs(mean.x-mean.y)
#         }else{
#             abs_log2FC = abs(log2(mean.x/mean.y))
#         }
#     })
#   stat$significance = -log10(stat[[significance_metric]])
#   selected_by_test = stat %>%
#     top_stats(metric = c('abs_log2FC', 'significance'),method = row_L1_norm,output.format = 'stats')%>% rowname
#     selected_by_test
# }
# hxlfoldensemble_iteration = function(matrix,y,tests_selector=wc_L1_selector,mc.cores = 1,logmat = F,...){
#   selected_by_test = tests_selector(matrix = matrix,y = y)
  
#   selected_by_ensemble = hxlensemble_selector(matrix = matrix[selected_by_test,] %>% t,
#                                               y =y,n.repeats =10,mc.cores = mc.cores,...)
  
#   selected_by_ensemble
# }
# hxlfoldensemble_selector = function(matrix, y, k.folds = 5, mc.cores = 1,tests_selector=wc_L1_selector,significance_metric='fdr',logmat = F, ...){
#   folds = createFolds(y = y,k =k.folds,returnTrain = T)
#   modeling_cores = max(floor(mc.cores/length(folds)),1)
#   folds %>% flexible.mclapply(mc.cores =min(length(folds),mc.cores),FUN = function(fold){
#     hxlfoldensemble_iteration(matrix = matrix[,fold],y = y[fold], mc.cores = modeling_cores,tests_selector=tests_selector,logmat=logmat,...)
#   }) %>% unlist %>% table %>% as.data.frame(stringsAsFactors = F) %>% setNames(nm = c('ROWNAMES','Freq'))
# }
hxl_trainer = function(x, y, alpha=0, family="binomial",random.seed = 50,...){
  arg_list = c(environment() %>% as.list, ...)
  set.seed(random.seed)
  alpha_predictor = do.call(glmnet, arg_list)
  alpha.fit <- do.call(cv.glmnet, arg_list)
  alpha_lambda = alpha.fit$lambda.1se
  hxl_model = list(alpha_predictor = alpha_predictor,alpha_lambda =alpha_lambda )
  if(is.factor(y)){
    hxl_model$y_levels = levels(y)
  }
  hxl_model %>%
    structure(class = c('hxl_model','list'))
}
predict.hxl_model = function(hxl_model, x){
  predicted = with(hxl_model,{
    predict(alpha_predictor, s=alpha_lambda, newx=x, type = "class") %>% as.character
  })
  if('y_levels' %in% names(hxl_model)){predicted = factor(predicted, hxl_model$y_levels)}
  return(predicted)
}

hxl_model_coef = function(hxl_model){
  with(hxl_model,{
    as.matrix(coef(alpha_predictor,s=alpha_lambda))[-1,1]
  })
}