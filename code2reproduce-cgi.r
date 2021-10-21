cacheHome = 'data/cacheHome'
cacheHome

source('code/session.utils2.r')

setCacheDir(cacheHome)

source('code/df.tools.r')
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

cgi = 'homer/data/genomes/hg19/annotations/basic/cpgIsland.ann.txt' %>%
read.table1(col.names = c('id','chr','start','end','V1','V2'),stringsAsFactors = F) %>% 
dplyr::select_at(c('chr','start','end','id')) %>% filter(chr %in% paste0('chr',c(1:22))) %>%
makeGRangesFromDataFrame(keep.extra.columns = T) %>% sort%>%
granges2bed('data/ref/cgi.bed',keep.mcols = T)

sch('cgi.counts',{
    cgi.counts = rawInfo%>%pull(sampleID)%>% 
    flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
        bed = subset(rawInfo, sampleID == sampleid)$bed_file
        cts = bed %>% routine.coverage(feature.input = cgi.bed,need.cols = 4:5) %>% setNames(sampleid)
        cts
    })%>% merge.df.list %>% trim_rawcts
    cgi.counts
})

sch('cgi.glioma_healthy.Dhm',{
    diff_colData = cfDNA_set %>% 
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })
    Dhms = de.edger(counts = cgi.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

sch('cgi.glioma_wbc.Dhm',{
    diff_colData = tissue_set %>% 
    within({
        y = ifelse(label=='wbc','wbc','glioma')%>%factor(c('wbc','glioma'))
    })
    Dhms = de.edger(counts = cgi.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']])
    Dhms %>% bind_rows
})

cgiInfo = read.table(cgi.bed,col.names = c('chr','start','end','id'),sep='\t',stringsAsFactors = F) %>% make.rownames('id')

table(cgi.glioma_healthy.Dhm$qvalue <0.05)

cgi.logfpkm = log2(naive_fpkm(
    cts = cgi.counts+1,
    lengths = cgiInfo[rownames(cgi.counts),]%>% with({end-start}),
    lib.sizes = vlookup(colnames(cgi.counts),rawInfo,'sampleID','lib_size')
))

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

HTMAP_markers= cgi.glioma_healthy.Dhm %>% arrange(qvalue) %>% head(1000)%>%pull(id)

HTMAP1000 = cgi.logfpkm[HTMAP_markers, HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = HTMAP_samples[, 'y', drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width=5, repr.plot.height=3) 
HTMAP1000

cgi.glioma_healthy.gain.q0.05 = cgi.glioma_healthy.Dhm %>% with({id[logFC>0&qvalue<0.05]})

cgi.glioma_healthy.loss.q0.05 = cgi.glioma_healthy.Dhm %>% with({id[logFC<0&qvalue<0.05]})

cgi.glioma_wbc.gain.q0.00001 = cgi.glioma_wbc.Dhm %>% with({ id[logFC>0&qvalue<0.00001]})

cgi.glioma_wbc.loss.q0.00001 = cgi.glioma_wbc.Dhm %>% with({id[logFC<0&qvalue<0.00001]})

cgi.glioma_wbc.gain.q0.05 = cgi.glioma_wbc.Dhm %>% with({ id[logFC>0&qvalue<0.05]})

cgi.glioma_wbc.loss.q0.05 = cgi.glioma_wbc.Dhm %>% with({id[logFC<0&qvalue<0.05]})

length(cgi.glioma_healthy.gain.q0.05)

length(cgi.glioma_wbc.gain.q0.00001)

length(intersect(cgi.glioma_healthy.gain.q0.05,cgi.glioma_wbc.gain.q0.00001))

length(intersect(cgi.glioma_healthy.gain.q0.05,cgi.glioma_wbc.gain.q0.05))

length(cgi.glioma_healthy.loss.q0.05)

length(cgi.glioma_wbc.loss.q0.00001)

length(intersect(cgi.glioma_healthy.loss.q0.05,cgi.glioma_wbc.loss.q0.00001))

length(intersect(cgi.glioma_healthy.loss.q0.05,cgi.glioma_wbc.loss.q0.05))

length(intersect(cgi.glioma_healthy.gain.q0.05,cgi.glioma_wbc.loss.q0.00001))

length(intersect(cgi.glioma_healthy.loss.q0.05,cgi.glioma_wbc.gain.q0.00001))

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

sch('cgi.glioma_healthy.Dhm.train',cacheSubDir = NULL,{
    train_diff_colData = cfDNA_set %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })

    Dhms = de.edger(counts = cgi.counts,colData = train_diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

topMarkers.genbody.glioma.train = topMarker_functions[[topMarkers_method]](cgi.glioma_healthy.Dhm.train)

train_gliomaModel_colData = gliomaModel_colData %>% filter(!sampleID %in% heldOut)
test_gliomaModel_colData = gliomaModel_colData %>% filter(sampleID %in% heldOut)

sch('cgi.gliomaModel',recreate=F,cacheSubDir = NULL,{
    set.seed(1234)
    caret_trainer(matrix = t(cgi.logfpkm[topMarkers.genbody.glioma.train, train_gliomaModel_colData$sampleID]), y = train_gliomaModel_colData$y, 
                  model="glmnet",preProcess = c("center", "scale", "pca"),
                  trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
                 )
})

sch('cgi.gliomaPrediction',recreate=F,cacheSubDir = NULL,reload=T,{
      cfDNA_set %>% within({
      assignment = case_when(sampleID %in% test_gliomaModel_colData$sampleID ~ "testing", 
                             sampleID %in% train_gliomaModel_colData$sampleID ~ 'training',TRUE~'N.A')
      model_RHS = RHS.glmnet(cgi.gliomaModel, t(cgi.logfpkm[topMarkers.genbody.glioma.train, sampleID]))
      probability = predict(cgi.gliomaModel,t(cgi.logfpkm[topMarkers.genbody.glioma.train, sampleID]),type='prob')[['glioma']]
          model_prediction = predict(cgi.gliomaModel, t(cgi.logfpkm[topMarkers.genbody.glioma.train, sampleID]))
      QA = vlookup(sampleID, within(gliomaModel_colData, 
        {
          y = as.character(y)
        }), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment, model_RHS,probability, model_prediction, QA)
})

cgi.gliomaPrediction = cgi.gliomaPrediction%>% 
merge(cfDNA_set,by='sampleID') %>%
mutate(GhmC_Score=model_RHS,
       label = ifelse(QA =='glioma','glioma',label)%>%
       factor(c(label_levels$label,'glioma'))%>%droplevels)

cgi.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({confusionMatrix(model_prediction,factor(QA,levels(model_prediction)),positive = 'glioma')})

suppressPackageStartupMessages({require(pROC)})

cgi.gliomaPrediction.roc = cgi.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({roc(QA,GhmC_Score,levels = c('others','glioma'))})

ci.auc(cgi.gliomaPrediction.roc)%>%round(3)

options(repr.plot.width=5,repr.plot.height=5)
cgi.gliomaPrediction.roc %>% 
plot(print.thres=0,main='Gliomas VS Non-Gliomas with CGIs',
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
cgi.gliomaPrediction %>% filter(assignment=='testing') %>%
(function(data){
    main_comparisons = list(c('nonglioma','glioma'),c('nonbrain','glioma'))
    show_cutoffs = c(0.5)
    
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

cfDNA_set_idh = cfDNA_set %>% filter(idh %in% subtypes[['idh']]|label =='healthy')

idhModel_colData = cfDNA_set_idh %>% 
filter(idh %in% nonAmbiguous[['idh']]) %>% 
within({
    y = factor(idh,nonAmbiguous[['idh']])
    age_scaled=(age-mean(age))/sd(age)
})

sch('cgi.idh.tissue.Dhm',recreate=F,{
    de.edger(counts = cgi.counts,
             colData = tissue_set %>% filter(idh %in% nonAmbiguous[['idh']])%>%
             within({
                y = factor(idh,nonAmbiguous[['idh']])
            }), coefs = list(IDHmut=2))%>%bind_rows
})

cgi.idh.tissue.Dhm.q0.05 = processors$cutBy(cgi.idh.tissue.Dhm,cutoff = 0.05, v = "qvalue")

cgi.idh.tissue.Dhm.q0.05 %>% length

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

IDH_HTMAP_markers= cgi.idh.tissue.Dhm.q0.05

IDH_HTMAP = cgi.logfpkm[IDH_HTMAP_markers, IDH_HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = IDH_HTMAP_samples[, c('y'), drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = IDH_HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

options(repr.plot.width=5, repr.plot.height=3) 
IDH_HTMAP

sch('cgi.idh.cfDNA.Dhm',recreate=F,reload = T,{
    train_diff_colData = cfDNA_set_idh %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(idh %in% c('healthy',nonAmbiguous[['idh']])) %>% 
    within({
        y = factor(idh,c('healthy',nonAmbiguous[['idh']]))
    })
    
    Dhms = de.edger(counts = cgi.counts,colData = train_diff_colData, 
             coefs = DEcoefs[['idh']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

cgi.idh.cfDNA.Dhm = df.split(cgi.idh.cfDNA.Dhm,'label')[c(nonAmbiguous[['idh']],paste0(nonAmbiguous[['idh']][2],'VS',nonAmbiguous[['idh']][1]))]

table(cgi.idh.cfDNA.Dhm$IDHmutVSIDHwt$qvalue <0.05 )

cgi.idh.cfDNA.Dhm$IDHmutVSIDHwt %>% filter(qvalue<0.05)

cgiInfo['CpG-14744',]

ggboxplot(idhModel_colData %>% within({
        CpG14744 = as.numeric(cgi.logfpkm['CpG-14744',sampleID])
    }),x = 'idh',y='CpG14744',outlier.shape = NA)+geom_jitter(aes(color=idh))+stat_compare_means()+
    xlab('CpG-14744')+ylab('cfDNA 5hmC logFPKM on CpG-14744')+coord_cartesian(clip = 'off')+theme(legend.position = 'none')

train_idhModel_colData = idhModel_colData %>% filter(!sampleID %in% heldOut)
test_idhModel_colData = idhModel_colData %>% filter(sampleID %in% heldOut)

cgi.idh.train_mat = t(cgi.logfpkm[cgi.idh.tissue.Dhm.q0.05, train_idhModel_colData$sampleID])

cgi.idh.train_preP = preProcess(cgi.idh.train_mat,method = c("center", "scale", "pca"))

cgi.idh.train_PCmat_age = cbind(predict(cgi.idh.train_preP,cgi.idh.train_mat),age=train_idhModel_colData$age_scaled)

cgi.idh.test_mat = t(cgi.logfpkm[cgi.idh.tissue.Dhm.q0.05, test_idhModel_colData$sampleID])

cgi.idh.test_PCmat_age = cbind(predict(cgi.idh.train_preP,cgi.idh.test_mat),age=test_idhModel_colData$age_scaled)

sch('cgi.idhModel',recreate=F,cacheSubDir = NULL,{
set.seed(1234)
    caret_trainer(
        matrix = cgi.idh.train_PCmat_age, y = train_idhModel_colData$y, 
        model="glmnet",trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
    )
})

sch('cgi.idhPrediction',recreate=F,cacheSubDir = NULL,{
   test_idhModel_colData %>% filter(idh %in% subtypes[['idh']]) %>% within({
       probability = predict(cgi.idhModel, cgi.idh.test_PCmat_age,type='prob')[[nonAmbiguous[['idh']][2]]]
      assignment = "testing"
      model_RHS = RHS.glmnet(cgi.idhModel, cgi.idh.test_PCmat_age)
      model_prediction = predict(cgi.idhModel, cgi.idh.test_PCmat_age)
      QA = vlookup(sampleID, within(idhModel_colData, {y = as.character(y)}), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment,probability, model_RHS,  model_prediction, QA)
})

cgi.idhPrediction = cgi.idhPrediction %>% 
merge(cfDNA_set,by='sampleID')

cgi.idhPrediction.roc = cgi.idhPrediction %>%
filter(assignment=='testing')%>% 
with({roc(QA,model_RHS,levels = nonAmbiguous[['idh']])})

options(repr.plot.width=5,repr.plot.height=5)
cgi.idhPrediction.roc %>% 
plot(print.thres=0,main=paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

idh_age_roc = cfDNA_set_idh %>% filter(idh %in% nonAmbiguous[['idh']]) %>% 
with({roc(idh,age,levels = nonAmbiguous[['idh']])}) 

set.seed(1234)
idh_age_roc.test = roc.test(idh_age_roc,cgi.idhPrediction.roc)

options(repr.plot.width=5,repr.plot.height=5)
idh_age_roc%>%
plot(main=paste0(paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),' with CGIs'),
     print.auc=T,print.auc.x=0,print.auc.y=0.13,print.auc.pattern='AUC using Age Alone: %.3f',print.auc.adj=1)
plot(cgi.idhPrediction.roc,add=T,
     print.auc=T,print.auc.x=0,print.auc.y=0.2,print.auc.pattern='AUC using Model: %.3f',print.auc.adj=1)
text(0,0.06,sprintf('p = %.3f',idh_age_roc.test$p.value),adj = 1)
