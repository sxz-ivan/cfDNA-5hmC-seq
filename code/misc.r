ezdcast = function(df,m,v){
    f = setdiff(colnames(df),c(m,v)) %>% paste0(collapse='+')%>% paste0('~',m)
    dcast(df,formula = as.formula(f),value.var = v)
}
trim_rawcts = function(x,FUN=rowMaxs,cutoff=3){
    if(is.data.frame(x)){x=as.matrix(x)}
    x[FUN(x)>=cutoff,]
}
cachePath = function(row){
    sub(pattern = '/.RData','.RData',do.call(paste,c(list(cacheHome),row,list(sep='/','.RData'))))
}

processors = list(
    tt = function(Refs,v='qvalue',topn=10,direction='both',...){
        topn=as.numeric(topn)
    Filter = switch(direction,both = rep(T,nrow(Refs)),loss=(Refs[['logFC']]<0),gain=(Refs[['logFC']]>0))
    Refs=Refs[Filter,]
    oneside_select(Refs[['id']],Refs[[v]],topn)
},cutBy = function(Refs,cutoff=0.05,v='qvalue',direction='both',...){
        cutoff = as.numeric(cutoff)
    Filter = switch(direction,both = rep(T,nrow(Refs)),loss=(Refs[['logFC']]<0),gain=(Refs[['logFC']]>0))
    Refs[['id']][(Refs[[v]]<cutoff)&Filter]
})
processors = c(processors,list(
    ttOverlap = function(Refs,v='qvalue',topn=10,direction='both',...){
        concordant = (Refs[['cfDNADE']][['logFC']]*Refs[['tissueDE']][['logFC']]>0)
        Refs %>% lapply(function(Ref){Ref[concordant,]}) %>% lapply(processors$tt,v,topn,direction) %>% Reduce(f = intersect)
    },
    cutOverlap = function(Refs,cutoff=0.05,v='qvalue',direction='both',...){
        concordant = (Refs[['cfDNADE']][['logFC']]*Refs[['tissueDE']][['logFC']]>0)
        Refs %>% lapply(function(Ref){Ref[concordant,]}) %>% lapply(processors$cutBy,cutoff,v,direction) %>% Reduce(f = intersect)
    }))
processors = c(processors,list(
    combineTT = function(Refs,v='qvalue',topn=10,direction='both',...){
    Refs %>% lapply(processors$tt,v,topn,direction) %>% unlist %>% unique
    },
    combineCutOverlap = function(Refs,cutoff=0.05,v='qvalue',direction='both',...){
        Refs %>% lapply(processors$cutOverlap,cutoff,v,direction) %>% unlist %>% unique
    },
    combineCutBy = function(Refs,cutoff=0.05,v='qvalue',direction='both',...){
        Refs %>% lapply(processors$cutBy,cutoff,v,direction) %>% unlist %>% unique
    },
    combineTToverlap = function(Refs,v='qvalue',topn=10,direction='both',...){
        Refs %>% lapply(processors$ttOverlap,v,topn,direction) %>% unlist %>% unique
    }
))

source('code/edgeR_tweaks.r')
de.edger =function(counts,colData,coefs,sampleID='sampleID',group='y',lib.size='lib_size',extra=NULL){
    fit = edger.quickfit(counts[,colData[[sampleID]]],colData[[group]],colData[[lib.size]],colData[extra])
    de = coefs %>% names %>% lapply(function(cn){
        edger.toptable(fit,cn,1,coef = coefs[[cn]])
    })%>%setNames(names(coefs)%>% paste0('VScontrol'))
    if(length(coefs)>1){
        VS_label = paste0(rev(names(coefs)),collapse = 'VS')
        vb = edger.toptable(fit,VS_label,1,contrast = c(0,-1,1,rep(0,length(extra))))
        de = c(de,list(vb)%>%setNames(VS_label))
    }
    de
}
source('code/selectn_pc_glmnet.r')
invCV = function(x){
    if(is.data.frame(x)){x = as.matrix(x)}
    data.frame(stringsAsFactors = F,id=rownames(x),invCV = 1/(x %>% rowcvs %>% abs))
}

invSD = function(x){
    if(is.data.frame(x)){x = as.matrix(x)}
    data.frame(stringsAsFactors = F,id=rownames(x),invSD = 1/(x %>% rowSds))
}

search_geneID = function(ENSEMBL){
    bitr(ENSEMBL, fromType = "ENSEMBL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Hs.eg.db)#%>%
#     df.dedup('ENSEMBL')%>%
#     df.dedup('ENTREZID')
}

extend_geneID = function(df,fgm){
    df=merge(df,fgm,by='id')
    ids = search_geneID(df[['ENSEMBL']])
    df %>% merge(ids,by='ENSEMBL')
}

mygse=function (named_x, DB){
    gse = GSEA(gene = named_x, TERM2GENE = DB[1:2])
    gse@result = within(gse@result, {
        Description = vlookup(ID, DB, "annotID", "annotTerm")
    })
    gse
}

myenrich = function (entrezIDs, DB, universe = NULL){
    enrich = enricher(gene = entrezIDs, TERM2GENE = DB[1:2], 
        universe = universe)
    enrich@result = within(enrich@result, {
        Description = vlookup(ID, DB, "annotID", "annotTerm")
    })
    enrich
}

df.rowlistN = function(df,names = apply(df,1,paste,collapse='_')){
    df.rowlist(df,names)
}

breakLongStr = function(x,every=30,at='_',with='\n'){
    pattern = sprintf('(.{1,%s})(%s)',every,at)
    replacement = sprintf('\\1%s%s',at,with)
    gsub(pattern,replacement,x)
}

RHS = function(model,...){
    UseMethod("RHS", model)
}

RHS.default = function(model,X){
    predict(model,X,type = 'prob')[[1]]
}

model_kits = list(glmnet = list(
    model="glmnet",
    tuneGrid = expand.grid(
        .alpha = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1),
        .lambda = c(0, 0.01, 0.05, 0.1, 0.2, 0.5)),
    predictor = RHS.glmnet))

model_preprocessors = list(pca = list(preProcess = c("center", "scale", "pca")))

model_selectors = list(oneSE10=list(selectionFunction = "oneSE10"))

model_optims = list(ROC = list(sampling = "down", summaryFunction = twoClassSummary,metric = "ROC"))

base_trainControl = list(method = "repeatedcv", number = 10, repeats = 3, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE)

decideQ = function(results){results %>% filter(AUC>0.5)}

myAUC = function(predictions){
    AUCs = list()
    for(i in 1:5){
    AUCs = c(AUCs,predictions%>% df.split('QA')%>% 
             lapply(sample_n,predictions$QA %>% table %>% min) %>% 
             bind_rows %>% with(auROC(ifelse(QA==levels(model_prediction)[1],0,1),model_RHS)))
    }
    AUCs %>% unlist %>% median
}

testAUC = function(predictions){
    predictions %>%
    filter(QA!='N.A'&assignment=='testing') %>% myAUC
}

trainAUC = function(predictions){
    predictions %>%
    filter(QA!='N.A'&assignment=='training') %>% myAUC
}
MODEL_produce = function(EXTRACTION,trans_mat,subdir_trunk,noload,recreate){
    EXTRACTION %>% MODEL_fit(trans_mat,subdir_trunk,noload,recreate) %>%
    MODEL_predict(trans_mat,subdir_trunk,noload,recreate)
}

my_multi_roc = function(data){
    multi_roc(data,d_2class = 'd',m_score = 'score',grouping = 'resample',ipalette = 'black',size = 0.1)+
    theme_classic()+theme(legend.position = 'none')
}

my_jitter_plot = function(data){
    data %>% ggplot(aes(x = y.label, y = score))+
    geom_boxplot(outlier.size =0,color = "black",fill = 'white')+
    geom_jitter(aes(color=color.label),size = 1.5,width =0.2)+
    coord_flip()+theme_classic()+
    theme(axis.text = element_text(size=15),legend.position = 'none')
}

logfpkm = function(raw_mat,feature_info,rawInfo){
    raw_mat = raw_mat[feature_info$id, rawInfo$sampleID] %>% as.matrix
    log2(naive_fpkm(raw_mat[,rawInfo$sampleID]+1,
    lengths = with(feature_info[rownames(raw_mat),],{end-start}),
    lib.sizes = rawInfo$lib_size))
}
blankspace = ggplot() + geom_blank() + theme(panel.background = element_blank())
save2pdf= function (plot_list, filename = "plot.pdf", perpage = 1, width = 16, 
    height = 25) {
    require(gridExtra)
    spliter = rep(1:ceiling(length(plot_list)/perpage), rep(perpage, 
        ceiling(length(plot_list)/perpage)))
    plot_lists = c(plot_list,rep(list(blankspace),length(spliter)-length(plot_list))) %>% split(spliter)
    pdf(filename, width, height, onefile = TRUE)
    for (pl in plot_lists) {
        print(plot_grid(plotlist = pl, ncol = 1))
    }
    dev.off()
    filename
}