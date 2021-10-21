cacheHome = 'data/cacheHome'
cacheHome

source('code/session.utils2.r')

setCacheDir(cacheHome)

source('code//df.tools.r')
source('code/metric.r')
source('code/bedior.R')
source('code/misc.r')
suppressPackageStartupMessages({
    require(cowplot)
    require(ggpubr)
})

plg = function(plotlist,...){plot_grid(plotlist = plotlist,...)}

testAUC2 = function (predictions) 
{
    predictions %>% filter(QA != "N.A" & assignment == "testing") %>% 
        with(auROC(ifelse(QA == levels(model_prediction)[1], 
                0, 1), model_RHS))
}

small_round = function(x,ndigits = 2){
    nz = log(x,0.1)
    round(x,nz+ndigits)
}

suppressPackageStartupMessages({require(openblasctl)})
openblas_set_num_threads(8)
openblas_get_num_threads()

rawInfo = read.csv('data/sampleInfo.csv',stringsAsFactors = F)

cfDNA_set = rawInfo %>% filter(sampletype == 'cfDNA')

tissue_set = rawInfo %>% filter(sampletype != 'cfDNA')

label_levels = list(label = c('wbc','healthy','nonbrain','nonglioma','nongbm','gbm'))

subtypes = list(mgmt = c('Unmethylad','WeaklyMethylated','Methylated'),
   codel = c('NonCodel','Codel'),
    tert= c('TERTwt','TERTmut'),
    atrx= c('ATRXretained','ATRXpartial','ATRXlost'),
   idh= c('IDHwt','IDHmut'),
   grade= c('Grade2','Grade3','Grade4'))

nonAmbiguous = list(grade=c('Grade2','Grade4'),idh=c('IDHwt','IDHmut'),
                codel = c('NonCodel','Codel'),tert=c('TERTwt','TERTmut'),
               atrx=c('ATRXretained','ATRXlost'),mgmt=c('Unmethylad','Methylated'))

my_palette = list(
    label=c(gbm='firebrick',nongbm='lightcoral',glioma='red',nonglioma='yellow',nonbrain='purple',healthy='skyblue'),
    grade=c(Grade2='pink',Grade3='lightcoral',Grade4='firebrick'),
    idh = c(IDHwt='red',IDHmut='gold'),
    atrx = c(ATRXretained='brown',ATRXpartial='orange', ATRXlost='gold'),
    mgmt = c(Unmethylad='gold', WeaklyMethylated='lightblue', Methylated='royalblue'),
    codel=c(NonCodel='brown', Codel='gold'),
    tert=c(TERTwt='pink', TERTmut='purple')
)

recurrent_levels = c('primary','recurrent')

DEcoefs = c(list(glioma=list(glioma=2)),nonAmbiguous %>% lapply(function(Q){list(2,3) %>% setNames(Q)}))

promoter.bed = 'homer/data/genomes/hg19/annotations/basic/promoters.ann.txt' %>%
read.table1(col.names = c('id','chr','start','end','V1','V2'),sep='\t',stringsAsFactors = F) %>% 
dplyr::select_at(c('chr','start','end')) %>% filter(chr %in% paste0('chr',c(1:22))) %>%distinct_all %>%
mutate(id=paste0(chr,':',start,'-',end))%>% makeGRangesFromDataFrame(keep.extra.columns = T) %>% sort %>%
granges2bed('data/ref/promoter.bed',keep.mcols = T)

sch('promoter.counts',{
    promoter.counts = rawInfo%>%pull(sampleID)%>% 
    flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
        bed = subset(rawInfo, sampleID == sampleid)$bed_file
        cts = bed %>% routine.coverage(feature.input = promoter.bed,need.cols = 4:5) %>% setNames(sampleid)
        cts
    })%>% merge.df.list %>% trim_rawcts
    promoter.counts
})

sch('promoter.glioma_healthy.Dhm',{
    diff_colData = cfDNA_set %>% 
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })
    Dhms = de.edger(counts = promoter.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

sch('promoter.glioma_wbc.Dhm',{
    diff_colData = tissue_set %>% 
    within({
        y = ifelse(label=='wbc','wbc','glioma')%>%factor(c('wbc','glioma'))
    })
    Dhms = de.edger(counts = promoter.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']])
    Dhms %>% bind_rows
})

promoterInfo = read.table(promoter.bed,col.names = c('chr','start','end','id'),sep='\t',stringsAsFactors = F) %>% make.rownames('id')


promoter = 'homer/data/genomes/hg19/annotations/basic/promoters.ann.txt' %>%
read.table1(col.names = c('id','chr','start','end','V1','V2'),sep='\t',stringsAsFactors = F) %>% mutate(REFSEQ = sub('\\)','',sub('.*\\(','',id)))

promoter_entrez = bitr(promoter$REFSEQ,fromType = 'REFSEQ',toType = c('ENTREZID'),OrgDb = org.Hs.eg.db)

promoter_ensembl = bitr(promoter$genes,fromType = 'REFSEQ',toType = c('ENSEMBL'),OrgDb = org.Hs.eg.db)

promoterFGM = promoter %>% merge(promoter_entrez,by='REFSEQ',all.x=T) %>% merge(promoter_ensembl,by='REFSEQ',all.x=T)%>%
mutate(id=paste0(chr,':',start,'-',end))%>%filter(!(is.na(ENTREZID)&is.na(ENSEMBL))) %>%
arrange(as.numeric(ENTREZID),ENSEMBL) %>% df.dedup('id') %>% dplyr::select(id,ENTREZID,ENSEMBL)

promoter.logfpkm = log2(naive_fpkm(
    cts = promoter.counts+1,
    lengths = promoterInfo[rownames(promoter.counts),]%>% with({end-start}),
    lib.sizes = vlookup(colnames(promoter.counts),rawInfo,'sampleID','lib_size')
))

nrow(promoter.glioma_healthy.Dhm)

table(promoter.glioma_healthy.Dhm$qvalue <0.05)

source('code/fixed.pheatmap.r')

HTMAP_samples = rawInfo%>% filter(disease %in%c('GLIOMA','HEALTHY','WBC'))%>%
within({y = case_when(
    disease=='GLIOMA'&sampletype=='t'~'glioma_tumor',
    disease=='GLIOMA'&sampletype=='cfDNA'~'glioma_plasma',
    disease=='HEALTHY'~'healthy_plasma',
    TRUE~'WBC'
)})%>%make.rownames('sampleID')

HTMAP_samples$disease %>% table

HTMAP_palette = list(y=c(glioma_plasma='pink',healthy_plasma='skyblue',glioma_tumor='red',WBC='forestgreen'))

HTMAP_markers= promoter.glioma_healthy.Dhm %>% arrange(qvalue) %>% head(1000)%>%pull(id)

HTMAP1000 = promoter.logfpkm[HTMAP_markers, HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = HTMAP_samples[, 'y', drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width=5, repr.plot.height=3) 
HTMAP1000

promoter.glioma_healthy.gain.q0.05 = promoter.glioma_healthy.Dhm %>% with({id[logFC>0&qvalue<0.05]})

promoter.glioma_healthy.loss.q0.05 = promoter.glioma_healthy.Dhm %>% with({id[logFC<0&qvalue<0.05]})

promoter.glioma_wbc.gain.q0.00001 = promoter.glioma_wbc.Dhm %>% with({ id[logFC>0&qvalue<0.00001]})

promoter.glioma_wbc.loss.q0.00001 = promoter.glioma_wbc.Dhm %>% with({id[logFC<0&qvalue<0.00001]})

promoter.glioma_wbc.gain.q0.05 = promoter.glioma_wbc.Dhm %>% with({ id[logFC>0&qvalue<0.05]})

promoter.glioma_wbc.loss.q0.05 = promoter.glioma_wbc.Dhm %>% with({id[logFC<0&qvalue<0.05]})

length(promoter.glioma_healthy.gain.q0.05)

length(promoter.glioma_wbc.gain.q0.00001)

length(promoter.glioma_wbc.gain.q0.05)

length(intersect(promoter.glioma_healthy.gain.q0.05,promoter.glioma_wbc.gain.q0.00001))

length(intersect(promoter.glioma_healthy.gain.q0.05,promoter.glioma_wbc.gain.q0.05))

length(promoter.glioma_healthy.loss.q0.05)

length(promoter.glioma_wbc.loss.q0.00001)

length(promoter.glioma_wbc.loss.q0.05)

length(intersect(promoter.glioma_healthy.loss.q0.05,promoter.glioma_wbc.loss.q0.00001))

length(intersect(promoter.glioma_healthy.loss.q0.05,promoter.glioma_wbc.loss.q0.05))

length(intersect(promoter.glioma_healthy.gain.q0.05,promoter.glioma_wbc.loss.q0.00001))

length(intersect(promoter.glioma_healthy.loss.q0.05,promoter.glioma_wbc.gain.q0.00001))

tcga_Info = read_feather('data/ref/blcspg_colData_atomic.feather')

# Download TCGA RNA-seq FPKMUQ data and extract expression matrix

tcga_Info = tcga_Info %>% filter(definition=='Primary solid Tumor')%>%within({
    cancer_by_site =ifelse(primary_site=='Brain','Glioma',paste0(primary_site,' cancer'))
    cancer_by_site = cancer_by_site%>%factor(unique(c('Glioma',cancer_by_site)))
})

promoter.overlap.gain = intersect(promoter.glioma_healthy.gain.q0.05,promoter.glioma_wbc.gain.q0.00001)%>%
vlookup(promoterFGM,'id','ENSEMBL') %>% intersect(rownames(tcga_fpkmuq))

promoter.overlap.loss = intersect(promoter.glioma_healthy.loss.q0.05,promoter.glioma_wbc.loss.q0.00001)%>%
vlookup(promoterFGM,'id','ENSEMBL') %>% intersect(rownames(tcga_fpkmuq))

options(repr.plot.width=8,repr.plot.height=4)
bind_rows(
    within(tcga_Info,{direction='Hyper-hmC Genes';
                      E=tcga_fpkmuq[promoter.overlap.gain,tcga_Info$barcode]%>%rowZscore%>%colMeans}),
    within(tcga_Info,{direction='Hypo-hmC Genes';
                      E=tcga_fpkmuq[promoter.overlap.loss,tcga_Info$barcode]%>%rowZscore%>%colMeans}))%>%
ggplot()+
geom_boxplot(aes_string(x = 'cancer_by_site',y='E',fill = 'direction'),
             color='black',position = position_dodge2(0.9),outlier.size = 0.1)+
facet_wrap(facets = 'direction')+geom_hline(yintercept = 0,linetype='dashed')+
xlab('Tumor Site')+ylab('Mean Expression Z score')+
theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))

suppressPackageStartupMessages({require(regioneR)})

promoter.universe = makeGRangesFromDataFrame(promoterInfo,keep.extra.columns = T)

sch('promoter.glioma.perm',recreate=T,{
    c('gain','loss') %>% flexible.mclapply(mc.cores = 2,function(direction){
        overlapPermTest(promoter.universe[get(sprintf('promoter.glioma_healthy.%s.q0.05',direction))], 
                    promoter.universe[get(sprintf('promoter.glioma_wbc.%s.q0.00001',direction))], 
                        ntimes = 1000, genome = promoter.universe)
    })
})

options(repr.plot.width=4,repr.plot.height=3)
c('gain','loss')%>%lapply(function(direction){
direction_split = promoter.glioma.perm[[direction]]
obj = direction_split %>% with(numOverlaps)
with(obj,bind_rows(
    data.frame(stringsAsFactors = F, n_overlaps = observed, type = 'observed'),
    data.frame(stringsAsFactors = F, n_overlaps = permuted, type = 'permuted')
    )%>%mutate(direction=direction,pval = pval))})%>%bind_rows %>% df.split('type')%>%with({
        significance = ifelse(observed$pval[1]<0.001,'p < 0.001',as.character(observed$pval[1]))
        permuted %>% ggboxplot(x='direction',y='n_overlaps',color='black')+
        geom_point(data = observed,mapping = aes(x=direction,y=n_overlaps,color = direction),shape=18,stroke=3,size=2)+
        geom_text(data = observed,mapping = aes(x=direction,y=n_overlaps+0.1*max(n_overlaps),label = significance))+
        scale_color_manual(values = c(gain='red',loss='blue'))+
    scale_x_discrete(labels = c(gain='Hyper-5hmC\nPromoters',loss='Hypo-5hmC\nPromoters'))+
        ylab('Number of Overlaps')+xlab(NULL)+
        theme(axis.title.x = element_text(size = 18),
              axis.text.x = element_text(size = 14),
              strip.text.x = element_text(size = 16),legend.position = 'none')
    })

sch('promoter.glioma.reverse',recreate=T,{
    c('gain','loss') %>% flexible.mclapply(mc.cores = 2,function(direction){
        overlapPermTest(promoter.universe[get(sprintf('promoter.glioma_healthy.%s.q0.05',direction))], 
                    promoter.universe[get(sprintf('promoter.glioma_wbc.%s.q0.00001',setdiff(c('gain','loss'),direction)))], 
                        ntimes = 1000, genome = promoter.universe)
    })
})

options(repr.plot.width=4,repr.plot.height=3)
c('gain','loss')%>%lapply(function(direction){
direction_split = promoter.glioma.reverse[[direction]]
obj = direction_split %>% with(numOverlaps)
with(obj,bind_rows(
    data.frame(stringsAsFactors = F, n_overlaps = observed, type = 'observed'),
    data.frame(stringsAsFactors = F, n_overlaps = permuted, type = 'permuted')
    )%>%mutate(direction=direction,pval = pval))})%>%bind_rows %>% df.split('type')%>%with({
        significance = ifelse(observed$pval[1]<0.001,'p < 0.001',as.character(observed$pval[1]))
        permuted %>% ggboxplot(x='direction',y='n_overlaps',color='black')+
        geom_point(data = observed,mapping = aes(x=direction,y=n_overlaps,color = direction),shape=18,stroke=3,size=2)+
        geom_text(data = observed,mapping = aes(x=direction,y=n_overlaps+0.1*max(n_overlaps),label = significance))+
        scale_color_manual(values = c(gain='red',loss='blue'))+
    scale_x_discrete(labels = c(gain='Hyper-5hmC in cfDNA\nHypo-5hmC in gDNA',loss='Hypo-5hmC in cfDNA\nHyper-5hmC in gDNA'))+
        ylab('Number of Overlaps')+xlab(NULL)+
        theme(axis.title.x = element_text(size = 18),
              axis.text.x = element_text(size = 8),
              strip.text.x = element_text(size = 16),legend.position = 'none')
    })

suppressPackageStartupMessages({
    require(clusterProfiler)
    library(DOSE)
    library(msigdbr)
    require(org.Hs.eg.db)
    require(meshes)
    require(MeSH.Hsa.eg.db)
})

annotationDB_files = list.files('data/ref/DOSE_annot',full.names = T) %>% 
as.list %>%setNames(sub('.txt','',list.files('data/ref/DOSE_annot'),fixed = T))

annotationDB_files %>% names %>% cat

annotationDBs = annotationDB_files %>% lapply(read.table,colClasses = 'character',header=T,sep='\t')

FGM = promoterFGM %>% filter(id %in%rownames(promoter.counts))

id2entrez = FGM

universe = FGM %>% pull(ENTREZID) %>% unique

promoter.cfDNA_Dhm_GL = id2entrez %>% filter(id %in% with(promoter.glioma_healthy.Dhm,{id[qvalue<0.05]}))

promoter.cfDNA_Dhm_ORA = names(annotationDBs)%>%flexible.mclapply(mc.cores=10,function(annotationDB){
    ORA = promoter.cfDNA_Dhm_GL %>% with({
        ora = myenrich(ENTREZID, annotationDBs[[annotationDB]], universe)
        ora@result = within(ora@result, {
            annotationDB = rep(annotationDB, length(ID))
        })
        ora
    })
})

promoter.cfDNA_Dhm_ORAplots = lapply(promoter.cfDNA_Dhm_ORA,function(ora){
    result=ora@result %>%filter(p.adjust<0.05)
    plot = NULL
    if(nrow(result)>0){
        plot = ora%>% dotplot(showCategory=6) + ggtitle(result$annotationDB[1])
    }
    plot
})
promoter.cfDNA_Dhm_ORAplots=promoter.cfDNA_Dhm_ORAplots[!sapply(promoter.cfDNA_Dhm_ORAplots,is.null)]

options(repr.plot.width=20,repr.plot.height = 5)
promoter.cfDNA_Dhm_ORAplots%>%lapply(function(plot){
    plot+theme(text = element_text(size=20),
               axis.text.y = element_text(size=20),
               axis.text.x = element_text(size=20),
              axis.title.x = element_text(size=20))
})

custom_trainControl=trainControl(method = "repeatedcv", number = 10, repeats = 3, returnData = FALSE, classProbs = TRUE, savePredictions = FALSE, 
             sampling = "down", summaryFunction = twoClassSummary,selectionFunction = "oneSE10")

custom_tuneGrid = expand.grid(
                        .alpha = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1),
                        .lambda = c(0, 0.01, 0.05, 0.1, 0.2, 0.5))

topMarker_functions = list(
    top300=function(Dhm){
        processors$tt(Dhm,'qvalue',300)
    },
    cut0.05=function(Dhm){
        processors$cutBy(Dhm,v='qvalue',cutoff =0.05)
    }
)

topMarkers_method = 'cut0.05'

gliomaModel_colData = cfDNA_set %>% 
within({
    y = ifelse(label%in%c('gbm','nongbm'),'glioma','others')%>%factor(c('others','glioma'))
})

set.seed(1234)
heldOut = with(cfDNA_set,sampleID[createDataPartition(y = label,p = 0.2,list = F)])

sch('promoter.glioma_healthy.Dhm.train',cacheSubDir = NULL,{
    train_diff_colData = cfDNA_set %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })

    Dhms = de.edger(counts = promoter.counts,colData = train_diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

topMarkers.genbody.glioma.train = topMarker_functions[[topMarkers_method]](promoter.glioma_healthy.Dhm.train)

train_gliomaModel_colData = gliomaModel_colData %>% filter(!sampleID %in% heldOut)
test_gliomaModel_colData = gliomaModel_colData %>% filter(sampleID %in% heldOut)

sch('promoter.gliomaModel',recreate=F,cacheSubDir = NULL,{
    set.seed(1234)
    caret_trainer(matrix = t(promoter.logfpkm[topMarkers.genbody.glioma.train, train_gliomaModel_colData$sampleID]), y = train_gliomaModel_colData$y, 
                  model="glmnet",preProcess = c("center", "scale", "pca"),
                  trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
                 )
})

sch('promoter.gliomaPrediction',recreate=F,cacheSubDir = NULL,reload=T,{
      cfDNA_set %>% within({
      assignment = case_when(sampleID %in% test_gliomaModel_colData$sampleID ~ "testing", 
                             sampleID %in% train_gliomaModel_colData$sampleID ~ 'training',TRUE~'N.A')
      model_RHS = RHS.glmnet(promoter.gliomaModel, t(promoter.logfpkm[topMarkers.genbody.glioma.train, sampleID]))
      probability = predict(promoter.gliomaModel,t(promoter.logfpkm[topMarkers.genbody.glioma.train, sampleID]),type='prob')[['glioma']]
          model_prediction = predict(promoter.gliomaModel, t(promoter.logfpkm[topMarkers.genbody.glioma.train, sampleID]))
      QA = vlookup(sampleID, within(gliomaModel_colData, 
        {
          y = as.character(y)
        }), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment, model_RHS,probability, model_prediction, QA)
})

promoter.gliomaPrediction = promoter.gliomaPrediction%>% 
merge(cfDNA_set,by='sampleID') %>%
mutate(GhmC_Score=model_RHS,
       label = ifelse(QA =='glioma','glioma',label)%>%
       factor(c(label_levels$label,'glioma'))%>%droplevels)

promoter.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({confusionMatrix(model_prediction,factor(QA,levels(model_prediction)),positive = 'glioma')})

suppressPackageStartupMessages({require(pROC)})

promoter.gliomaPrediction.roc = promoter.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({roc(QA,GhmC_Score,levels = c('others','glioma'))})

ci.auc(promoter.gliomaPrediction.roc)%>%round(3)

options(repr.plot.width=5,repr.plot.height=5)
promoter.gliomaPrediction.roc %>% 
plot(print.thres=0,main='Gliomas VS Non-Gliomas with Promoters',
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

immediatePosBelow = function(q,x){
    if(all(x>q)){p=0}else{p=(1:length(x))[x<=q][which.max(x[x<=q])]}
    p
}

getCutoffs = function(y,r,q){
    y[sapply(q,immediatePosBelow,r)]
}

options(repr.plot.width=8,repr.plot.height=3.5)
promoter.gliomaPrediction %>% filter(assignment=='testing') %>%
(function(data){
    main_comparisons = list(c('nonglioma','glioma'),c('nonbrain','glioma'))
    show_cutoffs = c(0.3,0.4,0.5,0.6,0.7)
    
    cutoffs = data.frame(stringsAsFactors = F,breaks = with(data,getCutoffs(GhmC_Score,probability,show_cutoffs)),
                         note=paste0(round(show_cutoffs*100),'%'))%>%arrange(breaks)
         plot = data %>% 
        ggplot(aes(x = label, y = GhmC_Score)) + 
        geom_boxplot(outlier.size = 0, color = "black", fill = "white")+
        geom_jitter(aes(color = label), size = 1.5, width = 0.2)+theme_classic() + xlab(NULL)+
    theme(axis.text.y = element_text(size = 15),axis.title.x = element_text(size = 15),
              axis.text.x = element_text(size = 15,angle=30,hjust=1), legend.position = "none")+
        scale_color_manual(values = my_palette[['label']])+
    geom_hline(yintercept = 0,linetype='dashed')+coord_cartesian(clip="off")
    plot_list = list(
        p1=plot+ylab('GhmC-Score')+
        scale_y_continuous(breaks = cutoffs$breaks,labels = round(cutoffs$breaks,1))+
        stat_compare_means(comparisons = main_comparisons),
        p2=plot +
        scale_y_continuous(breaks = cutoffs$breaks,labels = sprintf('%s',cutoffs$note),position = 'right')+
        stat_compare_means(comparisons = main_comparisons)+
        ylab('Probability of Glioma')
    )
    plg(plot_list,nrow=1)
})

options(repr.plot.width=3,repr.plot.height=3)
promoter.gliomaPrediction %>%
filter(assignment=='testing'&!is.na(max_diameter))%>%
dplyr::select(sampleID,patientID,probability,pathologyAbbr,idh,GhmC_Score,enhancement,max_diameter)%>%
(function(data){
    plm = summary(lm(GhmC_Score~max_diameter,data))$coefficients['max_diameter','Pr(>|t|)']
    ggplot(data,aes_string(x='max_diameter',y='GhmC_Score'))+
geom_smooth(method = "lm", formula = y ~ x,color='black',linetype='dashed',alpha=0.1)+
geom_point(color='black')+
    annotate(geom = 'text',y = 0.65,x=6.5,label=sprintf('p = %.3f',plm))+
scale_color_manual(values = my_palette$grade)+theme_classic()})

cfDNA_set_idh = cfDNA_set %>% filter(idh %in% subtypes[['idh']]|label =='healthy')

idhModel_colData = cfDNA_set_idh %>% 
filter(idh %in% nonAmbiguous[['idh']]) %>% 
within({
    y = factor(idh,nonAmbiguous[['idh']])
    age_scaled=(age-mean(age))/sd(age)
})

sch('promoter.idh.tissue.Dhm',recreate=F,{
    de.edger(counts = promoter.counts,
             colData = tissue_set %>% filter(idh %in% nonAmbiguous[['idh']])%>%
             within({
                y = factor(idh,nonAmbiguous[['idh']])
            }), coefs = list(IDHmut=2))%>%bind_rows
})

promoter.idh.tissue.Dhm.q0.05 = processors$cutBy(promoter.idh.tissue.Dhm,cutoff = 0.05, v = "qvalue")

promoter.idh.tissue.Dhm.q0.05 %>% length

IDH_HTMAP_samples = rawInfo%>% filter(idh %in%c('IDHmut','IDHwt'))%>%
within({y = case_when(
    idh=='IDHmut'&sampletype=='t'~'IDHmut_tumor',
    idh=='IDHmut'&sampletype=='cfDNA'~'IDHmut_plasma',
    idh=='IDHwt'&sampletype=='cfDNA'~'IDHwt_plasma',
    idh=='IDHwt'&sampletype=='t'~'IDHwt_tumor')
        grade = Hmisc::capitalize(grade)
       })%>%make.rownames('sampleID')

IDH_HTMAP_samples$codel %>% table

IDH_HTMAP_palette = list(y=c(IDHmut_plasma='yellow',IDHwt_plasma='coral',IDHwt_tumor='firebrick',IDHmut_tumor='orange')) %>% c(my_palette)%>%
lapply(function(p){c(p,N.A='white')})

IDH_HTMAP_markers= promoter.idh.tissue.Dhm.q0.05

IDH_HTMAP = promoter.logfpkm[IDH_HTMAP_markers, IDH_HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = IDH_HTMAP_samples[, c('y'), drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = IDH_HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width=5, repr.plot.height=3) 
IDH_HTMAP

sch('promoter.idh.cfDNA.Dhm',recreate=F,reload = T,{
    train_diff_colData = cfDNA_set_idh %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(idh %in% c('healthy',nonAmbiguous[['idh']])) %>% 
    within({
        y = factor(idh,c('healthy',nonAmbiguous[['idh']]))
    })
    
    Dhms = de.edger(counts = promoter.counts,colData = train_diff_colData, 
             coefs = DEcoefs[['idh']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

promoter.idh.cfDNA.Dhm = df.split(promoter.idh.cfDNA.Dhm,'label')[c(nonAmbiguous[['idh']],paste0(nonAmbiguous[['idh']][2],'VS',nonAmbiguous[['idh']][1]))]

table(promoter.idh.cfDNA.Dhm$IDHmutVSIDHwt$qvalue <0.05 )

train_idhModel_colData = idhModel_colData %>% filter(!sampleID %in% heldOut)
test_idhModel_colData = idhModel_colData %>% filter(sampleID %in% heldOut)

promoter.idh.train_mat = t(promoter.logfpkm[promoter.idh.tissue.Dhm.q0.05, train_idhModel_colData$sampleID])

promoter.idh.train_preP = preProcess(promoter.idh.train_mat,method = c("center", "scale", "pca"))

promoter.idh.train_PCmat_age = cbind(predict(promoter.idh.train_preP,promoter.idh.train_mat),age=train_idhModel_colData$age_scaled)

promoter.idh.test_mat = t(promoter.logfpkm[promoter.idh.tissue.Dhm.q0.05, test_idhModel_colData$sampleID])

promoter.idh.test_PCmat_age = cbind(predict(promoter.idh.train_preP,promoter.idh.test_mat),age=test_idhModel_colData$age_scaled)

sch('promoter.idhModel',recreate=F,cacheSubDir = NULL,{
set.seed(1234)
    caret_trainer(
        matrix = promoter.idh.train_PCmat_age, y = train_idhModel_colData$y, 
        model="glmnet",trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
    )
})

sch('promoter.idhPrediction',recreate=F,cacheSubDir = NULL,{
   test_idhModel_colData %>% filter(idh %in% subtypes[['idh']]) %>% within({
       probability = predict(promoter.idhModel, promoter.idh.test_PCmat_age,type='prob')[[nonAmbiguous[['idh']][2]]]
      assignment = "testing"
      model_RHS = RHS.glmnet(promoter.idhModel, promoter.idh.test_PCmat_age)
      model_prediction = predict(promoter.idhModel, promoter.idh.test_PCmat_age)
      QA = vlookup(sampleID, within(idhModel_colData, {y = as.character(y)}), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment,probability, model_RHS,  model_prediction, QA)
})

promoter.idhPrediction = promoter.idhPrediction %>% 
merge(cfDNA_set,by='sampleID')

promoter.idhPrediction.roc = promoter.idhPrediction %>%
filter(assignment=='testing')%>% 
with({roc(QA,model_RHS,levels = nonAmbiguous[['idh']])})

options(repr.plot.width=5,repr.plot.height=5)
promoter.idhPrediction.roc %>% 
plot(print.thres=0,main=paste0(paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),' with Promoters'),
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

idh_age_roc = cfDNA_set_idh %>% filter(idh %in% nonAmbiguous[['idh']]) %>% 
with({roc(idh,age,levels = nonAmbiguous[['idh']])}) 

set.seed(1234)
idh_age_roc.test = roc.test(idh_age_roc,promoter.idhPrediction.roc)

options(repr.plot.width=5,repr.plot.height=5)
idh_age_roc%>%
plot(main=paste0(paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),' with Promoters'),
     print.auc=T,print.auc.x=0,print.auc.y=0.13,print.auc.pattern='AUC using Age Alone: %.3f',print.auc.adj=1)
plot(promoter.idhPrediction.roc,add=T,
     print.auc=T,print.auc.x=0,print.auc.y=0.2,print.auc.pattern='AUC using Model: %.3f',print.auc.adj=1)
text(0,0.06,sprintf('p = %.3f',idh_age_roc.test$p.value),adj = 1)
