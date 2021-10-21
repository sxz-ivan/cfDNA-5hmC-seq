suppressPackageStartupMessages(require(edgeR))
source('code/df.tools.r')

edger.quickfit=function (counts = counts, group = group, lib.size = lib.size,...) 
{
    samples = data.frame(group=group,...)%>%droplevels
    formula = paste0('~',paste0(colnames(samples),collapse = '+'))
    d <- edgeR::DGEList(counts = counts, samples = samples, lib.size = lib.size)
    design <- model.matrix(as.formula(formula),data = samples)
    d = edgeR::calcNormFactors(d)
    d <- estimateDisp(d, design)
    fit <- glmQLFit(d, design)
    fit
}
edger.toptable = function(fit,label=NULL,logFC_direction = 1,...){
    suppressPackageStartupMessages({
        require(qvalue)
    })
    tab = glmQLFTest(fit,...)$table
        tab = tab %>% within({
                padj = p.adjust(PValue,method = 'BH')
                qvalue = qvalue(PValue)$qvalue
                label = label
                id = rownames(tab)
                logFC=logFC_direction*logFC
            })
        tab
}
combine.df = function(df1,df2){
    rownames1 = rownames(df1)
    rownames2 = rownames(df2)
    colnames1 = colnames(df1)
    for(colname in colnames1){
        col1 = df1[[colname]]
        col2 = df2[[colname]]
        if(is.factor(col1)|is.factor(col2)){
            if(is.factor(col1)){
                lv1 = levels(col1)
            }else{
                lv1 = unique(col1)
            }
            if(is.factor(col2)){
                lv2 = levels(col2)
            }else{
                lv2 = unique(col2)
            }
            lv = c(lv1,lv2) %>% unique
            df1[[colname]] = factor(col1,lv)
            df2[[colname]] = factor(col2,lv)
        }
    }
    df = bind_rows(df1,df2[,colnames1])
    rownames.df = c(rownames1,rownames2)
    rownames.df = ifelse(duplicated(rownames.df),paste0(rownames.df,'.2'),rownames.df)
    rownames(df) = rownames.df
    df
}
combine.DGEList = function(dge1,dge2){
    samples = combine.df(dge1$samples,dge2$samples)
    intersected.tags = intersect(rownames(dge1$counts),rownames(dge2$counts))
    counts = cbind(dge1$counts[intersected.tags,,drop=F],dge2$counts[intersected.tags,,drop=F])
    colnames(counts) = rownames(samples)
    genes = dge1$genes[intersected.tags,]
    DGEList(counts = counts,samples = samples,genes = genes)
}


exactMembership = function(x,dge1,dge2,
                           dispersion = "auto", rejection.region = "doubletail", big.count = 900, prior.count = 0.125){
    object = combine.DGEList(dge1%>%within({samples$group='dge1'}),dge2%>%within({samples$group='dge2'}))
    object_x = combine.DGEList(x%>%within({samples$group='x'}), object)%>%edgeR::calcNormFactors(method = "none")
    object$samples$norm.factors = object_x$samples$norm.factors[-1]
    object = edgeR::estimateCommonDisp(object)
    object = edgeR::estimateTagwiseDisp(object)
    
    rejection.region <- match.arg(rejection.region, c("doubletail", 
        "deviance", "smallp"))

    if (is.null(dispersion)) 
        dispersion <- "auto"
    if (is.character(dispersion)) {
        dispersion <- match.arg(dispersion, c("auto", "common", 
            "trended", "tagwise"))
        dispersion <- switch(dispersion, common = object$common.dispersion, 
            trended = object$trended.dispersion, tagwise = object$tagwise.dispersion, 
            auto = getDispersion(object))
    }
    ldisp <- length(dispersion)
    ntags <- nrow(object$counts)
    
    y <- object_x$counts
    lib.size <- object_x$samples$lib.size
    norm.factors <- object_x$samples$norm.factors
    if (is.null(rownames(y))) 
        rownames(y) <- paste("tag", 1:ntags, sep = ".")
    lib.size <- lib.size * norm.factors
    offset <- log(lib.size)
    lib.size.average <- exp(mean(offset))
    prior.count <- prior.count * lib.size/mean(lib.size)
    offset.aug <- log(lib.size + 2 * prior.count)
    abundance <- mglmOneGroup(y, dispersion = dispersion, offset = offset)
    e <- exp(abundance)
    
    cpm = cpm(object_x,log=T)
    sample.split = split(rownames(object_x$samples),object_x$samples$group)
    
    data.frame(row.names = rownames(y),
        p1 = get_exact.pvals(x,dge1,e,ntags,lib.size.average,dispersion,rejection.region,big.count),
    p2 = get_exact.pvals(x,dge2,e,ntags,lib.size.average,dispersion,rejection.region,big.count),
              p12 = get_exact.pvals(dge1,dge2,e,ntags,lib.size.average,dispersion,rejection.region,big.count),
              cpmx = cpm[,sample.split[['x']]],
               cpm1 = cpm[,sample.split[['dge1']]]%>%rowMeans,
               cpm2 = cpm[,sample.split[['dge2']]]%>%rowMeans)
}
get_exact.pvals =function(dge1,dge2,e,ntags,lib.size.average,dispersion,rejection.region,big.count){
    y1 = dge1$counts
    n1 <- ncol(y1)
    input.mean <- matrix(e, ntags, n1)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * dge1$samples$lib.size)
    y1 <- q2qnbinom(y1, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)

    y2 = dge2$counts
    n2 <- ncol(y2)
    input.mean <- matrix(e, ntags, n2)
    output.mean <- input.mean * lib.size.average
    input.mean <- t(t(input.mean) * dge2$samples$lib.size)
    y2 <- q2qnbinom(y2, input.mean = input.mean, output.mean = output.mean, 
        dispersion = dispersion)
    
    exact.pvals <- switch(rejection.region, doubletail = exactTestDoubleTail(y1, 
        y2, dispersion = dispersion, big.count = big.count), 
        deviance = exactTestByDeviance(y1, y2, dispersion = dispersion), 
        smallp = exactTestBySmallP(y1, y2, dispersion = dispersion))
    exact.pvals
}

medips.edger = function(counts,group,lib.size,diffnorm = "tmm"){
    if(!is.factor(group)){
        group = factor(group)
    }
    cat(paste0('DE for ',levels(group)[2],' over ',levels(group)[1],'\n'))
    d <- edgeR::DGEList(counts = counts, group = group,lib.size = lib.size)
        if (diffnorm == "tmm") {
            cat("Apply trimmed mean of M-values (TMM) for library sizes normalization...\n")
            d = edgeR::calcNormFactors(d)
        }
    if (diffnorm == "quantile" | diffnorm == "none") {
            cat("Skipping trimmed mean of M-values (TMM) library size normalization...\n")
            d = edgeR::calcNormFactors(d, method = "none")
        }
    d = edgeR::estimateCommonDisp(d)
                d = edgeR::estimateTagwiseDisp(d)
                de.com = edgeR::exactTest(d,pair = levels(group))
    rm(d)
    gc()
    de.com
}