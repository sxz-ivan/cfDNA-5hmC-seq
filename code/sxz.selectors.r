suppressPackageStartupMessages({
    require(tidyverse)
    require(matrixTests)
  require(Rfast)
})
source('code/basics.r')
source('code/matrix.stats.r')
source('code/df.tools.r')

t.selector = function(matrix, annotation ,scheme, rowname.column = 'sample', top.n = 2e6,
					rowSums.cutoff= NA, p.cutoff = 0.05, lfc.cutoff = 0.5,
					pool =F,fdr.method = NA,direction = 'both',output.format = 'selected'){
    matrix = matrix[, annotation %>% arrange_at(scheme)%>% pull(rowname.column)]
    if(!is.na(rowSums.cutoff)){matrix = matrix[rowSums(matrix)>rowSums.cutoff, ]}
	table.f = annotation%>% pull(scheme)
	if(is.factor(table.f)){table.f = table.f %>% droplevels}
    table.n = table.f %>% table %>% as.numeric
    result = matrix.t.test(matrix, n1 = table.n[1],n2 = table.n[2],pool = pool,pOnly = F)
    feature = matrix %>% rownames
    if(!is.na(fdr.method)){result$p.value = p.adjust(result$p.value,method = fdr.method)}
    if(is.na(p.cutoff)){
        p.value = result$p.value
        p.value = p.value[!is.na(p.value)]
        p.cutoff = sort(p.value)[min(length(p.value),top.n)]
    }
    p.value = result$p.value[result$p.value < p.cutoff]
    feature = feature[result$p.value < p.cutoff][order(p.value)]
    #t.stat = result$stat[result$p.value < p.cutoff][order(p.value)]
	
	lfc = log(rowMeans(matrix[feature,-(1:table.n[1])])+1e-6,2)-log(rowMeans(matrix[feature,1:table.n[1]]) + 1e-6, 2)
    
	if(direction =='both'){top.n = round(top.n/2)
    top.a = feature[lfc >lfc.cutoff] %>% head(n = top.n)
    top.b = feature[lfc <(-lfc.cutoff)] %>% head(n = top.n)
    top = c(top.a,top.b)}
    if(direction == 'L2'){top = feature[lfc >lfc.cutoff] %>% head(n = top.n)}
    if(direction == 'L1'){top = feature[lfc <(-lfc.cutoff)] %>% head(n = top.n)}
    if(output.format == 'selected'){output = top}else{output = list(selected = top, stat = result)}
    output
}
check.levels = function(f){
    if(is.factor(f)){f = f %>% droplevels}
    levels = f  %>% table %>% names
    names(levels) = c('L1','L2')
    levels
}

row_t_equalvar = function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = xx %>% col_split(f = g)
  stats = matrixTests::row_t_equalvar(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
row_t_welch= function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = x %>% col_split(f = g)
  stats = matrixTests::row_t_welch(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
row_wilcoxon_twosample= function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = x %>% col_split(f = g)
  stats = matrixTests::row_wilcoxon_twosample(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
row_t_paired= function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = x %>% col_split(f = g)
  stats = matrixTests::row_t_paired(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
row_wilcoxon_paired= function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = x %>% col_split(f = g)
  stats = matrixTests::row_wilcoxon_paired(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
row_f_var= function(x,g){
  if(is.factor(g)){g = g %>% droplevels }else{g = g %>% as.factor}
  ab_matrices = x %>% col_split(f = g)
  stats = matrixTests::row_f_var(ab_matrices[[1]], ab_matrices[[2]])
  stats
}
independence_tests = function(x, g,nlevels = NULL,method = g2tests){
  if(is.null(nlevels)){nlevels = rep(2,nrow(x))}
  ROWNAMES = rownames(x)
  if(is.factor(g)){g = g %>% droplevels %>% as.numeric}else{g = g %>% as.factor%>% as.numeric}
  g = g - 1 #vars must start from 0?
  g.nlevels = g %>% unique %>% length
  x = rbind(x,g) %>% t
  gc()
  stats = method(x, 1:(ncol(x)-1), ncol(x),c(nlevels,g.nlevels))%>%
    with({data.frame(statistic=statistic,df = df,pvalue = pchisq(statistic, df, lower.tail = FALSE),row.names = ROWNAMES)})
  stats
}
row_chi2 = function(x,g,nlevels = NULL){
  independence_tests(x,g,nlevels,chi2tests)
}
row_g2 = function(x,g,nlevels = NULL){
  independence_tests(x,g,nlevels,g2tests)
}

matrixtests_selector = function(matrix, y, top.n = NULL, p.cutoff = 0.05, p.adjust.method = 'BH',fdr.cutoff = NULL, FUN = row_wilcoxon_twosample, ..., output.format = 'selected'){
  stats = FUN(matrix, y,...)
  rm(matrix)
  gc()
  if(!is.null(p.adjust.method)){
    stats$fdr = stats[['pvalue']] %>% p.adjust(p.adjust.method)
  }
  if(output.format == 'stats'){return(stats)}
  if(!(is.null(p.adjust.method)|is.null(fdr.cutoff))){
    selected = rownames(stats)[stats$fdr < fdr.cutoff]
    fdr = stats$fdr[stats$fdr < fdr.cutoff]
    selected = selected[order(fdr)] 
  }else{
    p = stats[['pvalue']]
    selected = rownames(stats)[p<p.cutoff]
    p = p[p<p.cutoff]
    selected = selected[order(p)] 
  }
  if(!is.null(top.n)){
    selected = selected%>% head(n = top.n)
  }
  selected %>% na.omit
}
top_stats = function(stats, top = 0.01,metric = c('abs_log2FC','neglog_p'), method = row_L1_norm,output.format = 'stats'){
  suppressPackageStartupMessages({
    require(scales)
  })
  result = stats[metric] %>% apply(MARGIN = 2,FUN = scales::rescale) %>% method 
  subsetter = data.frame(pos = 1:length(result),result = result) %>% arrange(desc(result))%>%
    head(floor(ifelse(top<1,length(result)*top,top))) %>% pull(pos)
  if(output.format == 'stats'){
    return(stats[subsetter,])
  }
  subsetter
}

deseq2_stats = function(cts=NULL, y=NULL, lib.sizes=NULL, feature.lengths=NULL,output.format='result'){
    dds = build_dds(cts = cts, y = y, lib.sizes= lib.sizes,feature.lengths= feature.lengths)
    dds <- DESeq(dds)
    if(output.format=='result'){
        out = results(dds) %>% 
        .[,c('log2FoldChange','pvalue','padj')] %>% 
        setNames(nm = c('log2FoldChange','pvalue','fdr'))
    }else{out=dds}
    return(out)
}

edger_stats = function(cts=NULL, y=NULL, lib.sizes=NULL, feature.lengths=NULL,output.format='result'){
    Ann = data.frame(row.names = rownames(cts), Length = feature.lengths)
    dge <- DGEList(counts = cts, genes = Ann)

    dge <- calcNormFactors(dge,lib.size = lib.sizes)

    design <- model.matrix(~ boot$label) 

    dge <- estimateDisp(dge,design)
	fit <- glmQLFit(dge,design)
	qlf <- glmQLFTest(fit,coef=2)
    
    if(output.format=='result'){
        out = topTags(qlf,n=nrow(dge))%>% .$table %>% 
        .[rownames(dge),c('logFC','PValue','FDR')] %>%
        setNames(nm = c('log2FoldChange','pvalue','fdr'))
    }else{out=fit}
    
    rm(dge)
    rm(fit)
    rm(qlf)
    gc()
    return(out)
}

limma_stats = function(cts=NULL, y=NULL, lib.sizes=NULL, feature.lengths=NULL,output.format='result'){
    Ann = data.frame(row.names = rownames(cts), Length = feature.lengths)
    dge <- DGEList(counts = cts, genes = Ann)

    dge <- calcNormFactors(dge,lib.size = lib.sizes)

    design <- model.matrix(~ 0+y) 
    colnames(design) =c('firstlevel','secondlevel')

    dge <- voom(dge, design)
    rm(cts)
    gc()
    
    fit = lmFit(dge,design)%>%
    contrasts.fit(., makeContrasts(secondlevel-firstlevel, levels =design))%>%
    eBayes(., trend = TRUE)
    
    if(output.format=='result'){
        out = fit%>%
        topTable(., number = nrow(dge))%>%
        .[rownames(dge),c('logFC','P.Value','adj.P.Val')] %>%
        setNames(nm = c('log2FoldChange','pvalue','fdr'))
    }else{out=fit}
    
    rm(dge)
    rm(fit)
    gc()
    return(out)
}
logFoldChange = function(x,g){
    if (is.factor(g)) {
        g = g %>% droplevels
    }
    else {
        g = g %>% as.factor
    }
    ab_matrices = x %>% col_split(f = g)
    log2(rowMeans(ab_matrices[[2]])/rowMeans(ab_matrices[[1]]))
}
row_2by2 = function(x,b){
    if(!is.matrix(x)){x=as.matrix(x)}
    btab = table(b%>%droplevels)
    b1 = btab[1]
    b2 = btab[2]
    s1 = x[,b==names(b1)] %>% rowSums
    s2 = x[,b==names(b2)] %>% rowSums
    data.frame(row.names = rownames(x),a=s2,b=b2-s2,c=s1,d=b1-s1,n1=b1,n2=b2)
}

row_OR = function(x,b){
    row_2by2(x,b)%>%
    tab_OR
}

tab_OR = function(tab){
    tab%>%
    with({
        (a*d)/(b*c)
    })
}

tab_RR = function(tab){
    tab%>%
    with({
        (a*n2)/(b*n1)
    })
}

tab_fisher = function(tab, alternative="two.sided",...){
    tab %>%
    apply(MARGIN = 1,function(tab1){
       fisher.test(matrix(tab1[1:4],ncol = 2,byrow = T), alternative=alternative,...)$p.value
    }) 
}

row_fisher = function(x,b, alternative="two.sided",...){
    tab = row_2by2(x,b)
    data.frame(tab,pvalue = tab %>%tab_fisher,OR=tab%>%tab_OR,RR=tab%>%tab_RR)
}