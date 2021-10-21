suppressPackageStartupMessages({
	require(matrixStats)
})
source('code/parallel.tools.r')

rowCVs = function(matrix){
	matrix = as.matrix(matrix)
    rowMeans(matrix)/rowSds(matrix)
}

my_jaccard = function(x, y){
    sum(x*y)/sum(ifelse((x+y)>0,1,0))
}

pairwise_jaccard = function(data,FUN.jaccard = my_jaccard,diagnal = NA,mc.cores = 1){
	samples = colnames(data)
    combinations = combn(samples, 2, FUN = NULL, simplify = F)

    jaccard_results = flexible.mclapply(combinations, function(p){FUN.jaccard(data[,p[1]],data[,p[2]])},mc.cores = mc.cores)%>%unlist

    f = sapply(combinations, function(x){x[1]})
    r = sapply(combinations, function(x){x[2]})

    jaccard_melt = rbind(data.frame(forward = f , reverse = r, value = jaccard_results),data.frame(forward =r , reverse = f, value = jaccard_results))
    jaccard_melt = rbind(jaccard_melt, data.frame(forward = samples, reverse = samples, value =diagnal))

    jaccard_matrix = dcast(jaccard_melt, formula = forward~reverse,value.var = 'value')
    rownames(jaccard_matrix) = jaccard_matrix$forward
    jaccard_matrix$forward = NULL
    jaccard_matrix = jaccard_matrix[order(rownames(jaccard_matrix)),order(colnames(jaccard_matrix))]
    jaccard_matrix
}

row_L1_norm = function(x){
  abs(rowSums(x))
}

row_L2_norm = function(x){
  apply(x, MARGIN = 2, FUN = function(x){x^2})%>% rowSums
}

rowZscore = function(x){
    x = as.matrix(x)
    (x-rowMeans(x))/rowSds(x)
}
rowRescale = function(x,to=c(0,1)){
    x= as.matrix(x)
    bottom = to[1]
    top = to[2]
    Maxs = matrix(rowMaxs(x),nrow = nrow(x), ncol = ncol(x))
    Mins = matrix(rowMins(x),nrow = nrow(x), ncol = ncol(x))
    Ranges = Maxs-Mins
    Ranges = ifelse(Ranges==0,Inf,Ranges)
    a = (top-bottom)/(Maxs-Mins)
    b = top-a*Maxs
    a*x+b
}