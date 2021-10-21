suppressPackageStartupMessages({
  require(parallel)
})

flexible.mclapply = function(X, FUN,..., mc.cores =1, XRnames = T){
    args = list(...)
    if(XRnames&class(X)=='character'){Ynames = X}else{Ynames=names(X)}
  if(length(X)>0){
    if(!is.list(X)){X=as.list(X)}
      mc.cores = min(mc.cores, parallel::detectCores(), length(X))%>%max(1)
      indices = 1:length(X) %>% as.list
      if(!is.null(names(X))){indices = indices %>% setNames(nm = names(X))}
      eFUN = function(index){
        y = tryCatch(
          {y = do.call(FUN,c(list(X[[index]]),args))},
#           warning = function(cond){
#             message(paste0('Warning at element ',index))
#             message(cond)
#             return(y)
#           },
#             #https://stackoverflow.com/questions/35785142/suppress-warnings-using-trycatch-in-r
          error = function(cond){
            message(paste0('Error at element ',index))
            message(cond)
            return(NULL)
          },
          finally = {
          }
        )
        y
      }
      if(mc.cores==1){
        Y = lapply(indices,eFUN)
      }else{
        Y = mclapply(indices,eFUN, mc.cores = mc.cores)
      }
  }else{Y=list()}
    names(Y) = Ynames
  Y
}
