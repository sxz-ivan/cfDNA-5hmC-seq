suppressPackageStartupMessages({
  require(tidyverse)
#   require(fifer)
})
source('code/df.tools.r')

# bootstrap = function(factor, method = 'complete', size = NA, size.factor = 0.5,show = 'numeric',...){
#   #factor是包含分类标准为一列的表
#   #     complete, statified, at.least
#   suppressMessages({
#     df = data.frame(stringsAsFactors = F,factor= factor%>% as.character, index = 1:length(factor))
#     if(method =='by.least'){
#       df.table = df$factor %>% table %>% as.data.frame
#       sizes = df.table$Freq
#       if(is.na(size)){size = min(sizes)}
#         if(!is.na(size.factor)){size = size*size.factor}
#       sizes = pmin(size, sizes)
#       bootstrap = stratified(df,'factor',sizes,...)$index
#     }else{
#       if(is.na(size)){size =ceiling(length(factor)*size.factor)}else{size.factor =size/length(factor)}
#       if(method == 'complete'){
#         bootstrap = df$index %>% sample(size) %>% sort
#       }
#       if(method =='equal'){
#         bootstrap = stratified(df,'factor',size.factor,...)$index
#       }
#     }
#   })
#   bootstrap = bootstrap %>% sort
#   if(show=='factor'){bootstrap=factor[bootstrap]}
#   if(show=='logical'){out = rep(F,length(factor))
#   out[bootstrap] = T
#   bootstrap = out}
#   bootstrap
# }

# bootstrap.list = function(factor, n.bootstrap=10, method= 'complete', size = NA, size.factor = 0.5,...){
#   bootstrap.list = list()
#   for (i in 1:n.bootstrap){
#     bootstrap.list[[i]] = bootstrap(factor, method=method, size = size, size.factor = size.factor,...)
#   }
#   bootstrap.list %>% setNames(paste0(method,1:n.bootstrap))
# }

# createBootstraps = function(y, times =1 ,p = 0.7,method = 'equal',...){
#     bootstrap.list(factor = y, n.bootstrap = times, size.factor = p, method = method,...)
# }

# bootstrap.splits = function(factor,index,n.bootstrap = 1,n.split = 1,structure = 'shuffled',random.name = T,...){
#   package.list = bootstrap.list(factor=factor, n.bootstrap=n.bootstrap, ...)%>%
#     lapply(function(package){
#       splits = index %>% random.split(n.split = n.split,format = 'data')
#       names(splits) = stri_rand_strings(n.split,6,pattern = '[a-zA-Z]')
#       list(bootstrap = package, splits = splits)
#     })
#   names(package.list) = stri_rand_strings(n.bootstrap,6)
#   package.list
# }

step.seq.generator = function(x, step.length = 10){
	steps = floor(log(x,step.length))
	output = 1
	for(step in 0:steps){
		output = c(output, seq(step.length^step, min(x,step.length^(step+1)),step.length^step))
	}
	c(output,x) %>% unique
}

na.rm = function(index, f){
  index[!is.na(f[index])]
}

suppressPackageStartupMessages({
#   require(caret)
  require(tidyverse)
})

eqsample = function(x,size=length(x),times=floor(x/size)){
   xs = sample(x,length(x))
    size = min(size,length(x))
    
    if((length(x) %% size) > (size/2) ){
        f = rep(1:times,size) %>% sort
        xs = rep(xs,ceiling(size*times/length(x)))
        s = split(xs[1:length(f)],f)
    }else{
        xs = c(xs,xs)
        s = lapply((1:times)*ceiling(length(x)/times),function(x){
            i = x+1:size
            xs[ifelse(i > length(x), i-length(x), i)]
        })
    }
    s %>% lapply(sort)
}

eqsplit = function(f,size,times=1){
    i = 1:length(f)
    r = mapply(eqsample,split(i,f),size,MoreArgs = list(times=times),SIMPLIFY = F)
    lapply(1:times,function(t){do.call(c,r %>% lapply(function(b){b[[t]]}))})
}

# lko.index.list = function(f, k=1,format='lko.only', method = 'random.folds'){
#   #输出lko.index
#   #format='lko.only','zipped','separate'
#   #lko.only只包含了lko，没有ko
#   #separate是将lko和ko分别放两个list，但是position是对应的
#   #zipped将lko和对应的ko放在一个list element
#   #这个应该能够跟完全组合达到差不多的效果吧。。而且估计了类比比例
#   k = min(length(f), k)
#   n.folds = round(length(f)/k)
#   fold.list = createFolds(f, n.folds)
#   lko.list = (1:length(fold.list)) %>% as.list %>% lapply(function(ko.index){
#     ko = fold.list[ko.index] %>% unlist
#     lko = fold.list[-ko.index] %>% unlist
#     if(format=='zipped'){lko = list(lko = lko, ko=ko)}
#     lko
#   })
#   if(format=='separate'){lko.list = list(lko.list = lko.list, ko.list = fold.list)}
#   lko.list
# }

# createLko = function(y, times = NULL, p = NULL, k=1,...){
#     if(!is.null(p)){k = floor(y*(1-p))}
#     lko.list = lko.index.list(f = y, k = k)
#     sample(lko.list,min(length(lko.list),times),replace = F)
# }

# balanced_sample = function(y,k=NULL){
#     y_rmna = ifelse(is.na(y),'NA',y)
#     y_split = 1:length(y) %>% split(y_rmna)
#     y_split = y_split[names(y_split) %>% pop('NA')]
#     k = min(k,sapply(y_split,length))
#     lapply(y_split,function(x){sample(x,k)}) %>% unlist
# }