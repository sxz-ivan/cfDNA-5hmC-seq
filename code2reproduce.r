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

library(gridExtra)
library(grid)
require(gtable)

showTissue = tissue_set%>% 
within({
    grade = case_when(
        sampletype =='wbc'~'WBC',
        grade %in% subtypes$grade & pathologyFull!='low grade glioma'~grade,
        TRUE~'N.A'
    )%>%
    factor(c('N.A','WBC','Grade1',subtypes$grade))%>%droplevels
    showDisease = Hmisc::capitalize(pathologyFull)
#     showDisease = Hmisc::capitalize(ifelse(grade=='WBC','WBC',pathologyFull))
})

Info_tissue = showTissue%>% group_by(grade,showDisease) %>% 
    summarise(N = n())%>% mutate(N = ifelse(showDisease=='WBC',as.character(N),paste0(showDisease,': ',N))) %>%
    group_by(grade) %>% summarise(N = paste0(c(N,''),collapse='\n'))%>%mutate(grade=as.character(grade))%>%
setNames(c('Sample type','N'))%>%
rbind(c('Tumor',as.character(sum(showTissue$showDisease!='WBC'))))%>%.[c(1,5,2:4),]

options(repr.plot.res=120,repr.plot.height=3,repr.plot.width=4)
Info_tissue%>% tableGrob(rows=NULL,theme = ttheme_minimal(core=list(fg_params=list(vjust=1,y=0.8)))) %>% 
gtable_add_grob(
        grobs = grobTree(
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(1,"npc"),x1 = unit(1,"npc"),y1 = unit(1,"npc"),gp = gpar(lwd = 4)),
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(0,"npc"),x1 = unit(1,"npc"),y1 = unit(0,"npc"),gp = gpar(lwd = 2))
        ),t = 1, b = 1, l = 1, r = ncol(.))%>% 
gtable_add_grob(
        grobs = grobTree(
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(0,"npc"),x1 = unit(1,"npc"),y1 = unit(0,"npc"),gp = gpar(lwd = 4))
        ),t = nrow(.), b = nrow(.), l = 1, r = ncol(.))%>% 
grid.draw

showInfo = cfDNA_set %>% 
within({
    grade=sub('grade','Grade ',ifelse(is.na(grade),'-',grade))
    factor(c('N.A','Grade1',subtypes$grade))
    showGroup = ifelse(label %in% c('gbm','nongbm'),'Glioma','Non-Glioma')
    showDisease = Hmisc::capitalize(ifelse(is.na(pathologyFull),case_when(
        label == 'healthy'~'Healthy',
        disease =='HCC'~'Hepatic Cancer',
        disease =='PDAC'~'Pancreatic Cancer'
    ),pathologyFull))
})

small_round = function(x,ndigits = 2){
    nz = log(x,0.1)
    round(x,nz+ndigits)
}

Info_cfDNA = bind_rows(
    showInfo %>% group_by(showGroup) %>% summarise(
        median_age = median(age,na.rm = T),
        min_age = min(age),max_age=max(age),N=n(),male = sum(gender=='M')
    )%>%mutate(
        Age = paste0(round(median_age,1),' (',min_age,' - ',max_age,')'),
        Male = paste0(male,' (',round(100*male/N),'%)'),
        Total = paste0('N = ',N)
    )%>% dplyr::select(Total,Age,Male,showGroup) %>% make.rownames('showGroup',F) %>% t %>% keep.rownames('rowHeader'),
    data.frame(stringsAsFactors = F,grade='WHO grade',a='',b='')%>% setNames(c('rowHeader','Glioma','Non-Glioma')),
    showInfo %>% group_by(grade,showGroup,showDisease) %>% summarise(N = n())%>% mutate(N = paste0(showDisease,': ',N)) %>%
    group_by(grade,showGroup) %>% summarise(N = paste0(c(N,''),collapse='\n')) %>% dcast(grade~showGroup,value.var = 'N',fill = '-')%>% 
    setNames(c('rowHeader','Glioma','Non-Glioma'))
)%>%select_at(c('rowHeader','Glioma','Non-Glioma'))%>%
setNames(c('Plasma cfDNA','Glioma','Non-Glioma'))%>%
within({
    P=rep('',length(Glioma))
    P[2]=with(showInfo,wilcox.test(age~showGroup))$p.value %>%small_round
    P[3]=with(showInfo,chisq.test(table(gender,showGroup)))$p.value%>%small_round
})%>%setRowNames(NULL)

Info_cfDNA

options(repr.plot.res=120,repr.plot.height=6,repr.plot.width=15)
Info_cfDNA %>% tableGrob(rows=NULL,theme = ttheme_minimal(core=list(fg_params=list(vjust=1,y=0.8)))) %>% 
gtable_add_grob(
        grobs = grobTree(
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(1,"npc"),x1 = unit(1,"npc"),y1 = unit(1,"npc"),gp = gpar(lwd = 4))
        ),t = 1, b = 1, l = 1, r = ncol(.))%>% 
gtable_add_grob(
        grobs = grobTree(
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(0,"npc"),x1 = unit(1,"npc"),y1 = unit(0,"npc"),gp = gpar(lwd = 2))
        ),t = 2, b = 2, l = 1, r = ncol(.))%>%
gtable_add_grob(
        grobs = grobTree(
            segmentsGrob(x0 = unit(0,"npc"),y0 = unit(0,"npc"),x1 = unit(1,"npc"),y1 = unit(0,"npc"),gp = gpar(lwd = 4))
        ),t = nrow(.), b = nrow(.), l = 1, r = ncol(.))%>% 
grid.draw

showIDH = cfDNA_set %>% filter(idh %in% nonAmbiguous[['idh']])

Info_IDH = bind_rows(
    showIDH %>% group_by(idh) %>% summarise(
        mean_age = mean(age,na.rm = T),
        min_age = min(age),max_age=max(age),N=n(),male = sum(gender=='M')
    )%>%mutate(
        Age = paste0(round(mean_age,1),' (',min_age,' - ',max_age,')'),
        Male = paste0(male,' (',round(100*male/N),'%)'),
        Total = paste0('N = ',N)
    )%>% dplyr::select(Total,Age,Male,idh) %>% make.rownames('idh',F) %>% t %>% 
keep.rownames('rowHeader')%>% select_at(c('rowHeader','IDHmut','IDHwt')),
data.frame(stringsAsFactors = F,grade='WHO grade',a='',b='')%>% 
setNames(c('rowHeader','IDHmut','IDHwt')),
showIDH %>% group_by(grade,idh) %>% summarise(N = n()) %>% 
dcast(grade~idh,value.var = 'N',fill = '-')%>% 
setNames(c('rowHeader','IDHmut','IDHwt'))%>%
dplyr::slice(match(c(subtypes[['grade']],'N.A'),rowHeader)%>%na.omit)%>%
mutate(rowHeader=sub('Grade','Grade ',rowHeader)),
data.frame(stringsAsFactors = F,grade='TERT promoter',a='',b='')%>% 
setNames(c('rowHeader','IDHmut','IDHwt')),
showIDH %>% group_by(tert,idh) %>% summarise(N = n()) %>% 
dcast(tert~idh,value.var = 'N',fill = '-')%>% 
setNames(c('rowHeader','IDHmut','IDHwt'))%>%
dplyr::slice(match(c(subtypes[['tert']],'N.A'),rowHeader)%>%na.omit)%>%
mutate(rowHeader=c('Wildtype','Mutant','Undetermined')),
data.frame(stringsAsFactors = F,grade='1p19q Codeletion',a='',b='')%>% 
setNames(c('rowHeader','IDHmut','IDHwt')),
showIDH %>% group_by(codel,idh) %>% summarise(N = n()) %>% 
dcast(codel~idh,value.var = 'N',fill = '-')%>% 
setNames(c('rowHeader','IDHmut','IDHwt'))%>%
dplyr::slice(match(c(subtypes[['codel']],'N.A'),rowHeader)%>%na.omit)%>%
mutate(rowHeader=c('Non-codel','Codel','Undetermined')),
data.frame(stringsAsFactors = F,grade='ATRX expression',a='',b='')%>% 
setNames(c('rowHeader','IDHmut','IDHwt')),
showIDH %>% group_by(atrx,idh) %>% summarise(N = n()) %>% 
dcast(atrx~idh,value.var = 'N',fill = '-')%>% 
setNames(c('rowHeader','IDHmut','IDHwt'))%>%
dplyr::slice(match(c(subtypes[['atrx']],'N.A'),rowHeader)%>%na.omit)%>%
mutate(rowHeader=c('Expressed','Partially expressed','Loss','Undetermined')),
data.frame(stringsAsFactors = F,grade='MGMT promoter methylation',a='',b='')%>% 
setNames(c('rowHeader','IDHmut','IDHwt')),
showIDH %>% group_by(mgmt,idh) %>% summarise(N = n()) %>% 
dcast(mgmt~idh,value.var = 'N',fill = '-')%>% 
setNames(c('rowHeader','IDHmut','IDHwt'))%>%
dplyr::slice(match(c(subtypes[['mgmt']],'N.A'),rowHeader)%>%na.omit)%>%
mutate(rowHeader=c('Unmethylad','Weakly methylated','Methylated','Undetermined'))
)%>%select_at(c('rowHeader','IDHmut','IDHwt'))%>%setNames(c('Plasma cfDNA','IDHmut','IDHwt'))%>%
within({
    P=rep('',length(IDHmut))
    P[2]=with(showIDH,wilcox.test(age~idh))$p.value %>%small_round
    P[3]=with(showIDH,chisq.test(table(gender,idh)))$p.value%>%small_round
})%>%setRowNames(NULL)

Info_IDH

DEcoefs = c(list(glioma=list(glioma=2)),nonAmbiguous %>% lapply(function(Q){list(2,3) %>% setNames(Q)}))

human_simple.genome.bed = 'data/ref/human_simple.genome' %>% 
awk.filter(need.cols = '1,"\t",1,"\t",$2',OFS = "",output = 'data/ref/human_simple.genome.bed')

human_simple.autosomal.bed = human_simple.genome.bed %>% shell.io('head -n 22',output = 'data/ref/human_simple.autosomal.bed')

hg19_blacklist_simple.bed = 'data/ref/wgEncodeDacMapabilityConsensusExcludable.bed' %>% 
bedtools.sort(output = 'data/ref/wgEncodeDacMapabilityConsensusExcludable_simple.bed')

bin300.bed = bedtools.shell(hg19_blacklist_simple.bed,
               command =paste0(bedtools, " complement %s -i "),
               arguments = paste0("-g ", bedtools.genome.file), output ='pipe') %>% 
bedtools.makewindows(window_size = 300,output = 'pipe')%>%
bedtools.intersect(input.b = human_simple.autosomal.bed,arguments = '-wa',output='pipe')%>%
awk.filter(need.cols = '1,"\t",$2,"\t",$3,"\t",$1,":",$2,"-",$3',OFS = "",
           output = 'data/ref/hg19_bin300_simple.mappable.bed')

sch('bin300.counts',{
    bin300.counts = rawInfo%>%pull(sampleID)%>% 
    flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
        bed = subset(rawInfo, sampleID == sampleid)$bed_file
        cts = bed %>% routine.coverage(feature.input = bin300.bed,need.cols = 4:5) %>% setNames(sampleid)
        cts
    })%>% merge.df.list %>% trim_rawcts
    bin300.counts
})

sch('bin300.glioma_healthy.Dhm',{
    diff_colData = cfDNA_set %>% 
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })
    Dhms = de.edger(counts = bin300.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

sch('bin300.glioma_wbc.Dhm',{
    diff_colData = tissue_set %>% 
    within({
        y = ifelse(label=='wbc','wbc','glioma')%>%factor(c('wbc','glioma'))
    })
    Dhms = de.edger(counts = bin300.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']])
    Dhms %>% bind_rows
})

gene.gtf = rtracklayer::import('data/ref/gencode.v19.annotation.gtf') %>% as.data.frame(stringsAsFactors = F)%>%
filter(type=='gene') %>% make.rownames('gene_id')

protein_codingAutosomalGenes = gene.gtf %>% 
filter(seqnames %in% paste0('chr',1:22)&transcript_type=='protein_coding') %>%
select_at(c('seqnames','start','end','strand','gene_id'))

protein_codingAutosomalGenesOver100bp = protein_codingAutosomalGenes %>% filter((end-start)>=100)

forward.gr = protein_codingAutosomalGenesOver100bp %>% filter(strand=='+') %>% makeGRangesFromDataFrame(keep.extra.columns = T)

reverse.gr = protein_codingAutosomalGenesOver100bp %>% filter(strand=='-') %>% makeGRangesFromDataFrame(keep.extra.columns = T)

width(reverse.gr) %>% summary

width(forward.gr) %>% summary

max_size_each_bin = 300
n_flank = 100

forward.tiles = forward.gr %>% tile(n = 100) %>% flexible.mclapply(mc.cores = 120,function(gr){
    ntile = length(gr)
    size_each_bin = min(width(gr),max_size_each_bin)
    upstream = gr[1]
    start(upstream) = max(start(gr[1])-size_each_bin*(n_flank+1),0)
    end(upstream) = start(gr[1])
    upstream = upstream %>% tile(n=n_flank+1) %>% .[[1]]
    upstream$relative = -(n_flank:0)
    
    downstream = gr[ntile]
    end(downstream) = end(gr[ntile])+size_each_bin*(n_flank+1)
    start(downstream) = end(gr[ntile])+1
    downstream = downstream %>% tile(n=n_flank) %>% .[[1]]
    downstream$relative = ntile+1:n_flank
    
    gr$relative = 1:100
    locus = c(upstream,gr,downstream)
    locus$id = names(locus)
    locus
})

reverse.tiles = reverse.gr %>% tile(n = 100) %>% flexible.mclapply(mc.cores = 120,function(gr){
    ntile = length(gr)
    size_each_bin = min(width(gr),max_size_each_bin)
    downstream = gr[1]
    start(downstream) = max(start(gr[1])-size_each_bin*(n_flank+1),0)
    end(downstream) = start(gr[1])
    downstream = downstream %>% tile(n=n_flank) %>% .[[1]]
    downstream$relative = ntile+n_flank:1
    
    upstream = gr[ntile]
    end(upstream) = end(gr[ntile])+size_each_bin*(n_flank+1)
    start(upstream) = end(gr[ntile])+1
    upstream = upstream %>% tile(n=n_flank+1) %>% .[[1]]
    upstream$relative = -(0:n_flank)
    
    gr$relative = 100:1
    locus = c(downstream,gr,upstream)
    locus$id = names(locus)
    locus
})

gene.bins = c(reverse.tiles,forward.tiles)

gene.bins2 =as(gene.bins,"GRangesList") %>% unlist

gene.relative_positions.bed = gene.bins2 %>% setNames(NULL) %>%sort%>%granges2bed('data/cacheHome/gene.bins2.bed',keep.mcols = T)

bin300.glioma_healthy.Dhm.q0.05 = bin300.glioma_healthy.Dhm %>% with(id[qvalue<0.05])

bin300.glioma_healthy.Dhm.q0.05 %>% length

bin300.glioma_healthy.Dhm.q0.05.gr = bin300.glioma_healthy.Dhm.q0.05 %>% id2granges

start(bin300.glioma_healthy.Dhm.q0.05.gr)=start(bin300.glioma_healthy.Dhm.q0.05.gr)+1

bin300.glioma_healthy.Dhm.q0.05.bed = bin300.glioma_healthy.Dhm.q0.05.gr %>% 
granges2bed('data/cacheHome/bin300.glioma_healthy.Dhm.q0.05.bed',keep.mcols = T)

bin300.glioma_healthy.Dhm.relative = bin300.glioma_healthy.Dhm.q0.05.bed %>% 
bedtools.intersect(gene.relative_positions.bed,arguments = '-wo',output = 'pipe') %>% 
awk.filter(need.cols = c(4,8,9),output='R') %>% 
colsplit(pattern = '\t',names = c('id','relative','gene_id'))

bin300.glioma_healthy.Dhm.relative.stat = bin300.glioma_healthy.Dhm.relative %>% 
group_by(relative) %>% summarise(n = n())%>%as.data.frame(stringsAsFactors = F) %>% arrange(relative) %>%
mutate(fraction=n/length(bin300.glioma_healthy.Dhm.q0.05))

Fig.1C = bin300.glioma_healthy.Dhm.relative.stat %>% 
ggplot(aes(x = relative, y = fraction))+
geom_point(size=0.5)+geom_line()+geom_vline(xintercept = c(0,100),linetype='dashed')+
xlab('Relative Position of Genes in Percentile')+
ylab('Fraction of Differential 5hmC Sites')+coord_cartesian(clip='off')+
scale_x_continuous(limits = c(-110,220),breaks = c(-100,0,100,200),labels = c('- Gene Width/30kbp','TSS','TES','+ Gene Width/30kbp'))+
theme_classic()+ 
theme(axis.text.x=element_text(hjust=c(0.4,0.5,0.5,0.7)))

options(repr.plot.height=3,repr.plot.width=5)
Fig.1C

ggsave2('data/output/Fig.1C.svg',Fig.1C,width = 5,height = 3)

bin300_universe = bin300.bed %>% bed2granges

source('code/LOLAenrichment.r')

# hg19 annotation data files accompanying HOMER suite were filtered 
# for regions on autosomal chromosomes to prepare regionDB files as used here
# CpG shores and shelves were prepared as described in https://www.r-bloggers.com/2014/03/cpg-island-shelves/

# cgi = 'homer/data/genomes/hg19/annotations/basic/cpgIsland.ann.txt' %>%
# read.table1(col.names = c('id','chr','start','end','V1','V2'),stringsAsFactors = F) %>% 
# dplyr::select_at(c('chr','start','end','id')) %>% filter(chr %in% paste0('chr',c(1:22))) %>%
# makeGRangesFromDataFrame(keep.extra.columns = T) %>% sort

# # https://www.r-bloggers.com/2014/03/cpg-island-shelves/
# #             Extract CpG island shores
# ###############################################################
# # extract the shore defined by 2000 bp upstream of cpg islands
# shore1=flank(cgi, 2000)
# # extract the shore defined by 2000 bp downstream of cpg islands
# shore2=flank(cgi,2000,FALSE)
# # perform intersection and combine the shores where they overlap
# shore1_2=reduce(c(shore1,shore2))
# # extract the features (ranges) that are present in shores only and not in cgi (ie., shores not overlapping islands)
# cpgi_shores=setdiff(shore1_2, cgi)

# ###############################################################
# #             Extract CpG island shelves
# ###############################################################
# # extract the shore defined by 4000 bp upstream of cpg islands
# shelves1=flank(cgi, 4000)
# # extract the shore defined by 2000 bp downstream of cpg islands
# shelves2=flank(cgi,4000,FALSE)
# # perform intersection and combine the shelves where they overlap
# shelves1_2=reduce(c(shelves1,shelves2))
# # create a set of ranges consisting CpG Islands, Shores
# island_shores=c(cgi, cpgi_shores)
# # extract the features (ranges) that are present in shelves only and not in cgi  or shores(ie., shelves not overlapping islands or shores)
# cpgi_shelves=setdiff(shelves1_2, island_shores)


regionDB = loadRegionDB('data/ref/LOLADB/')

bin300.glioma_healthy.Dhm.enrichment = bin300_universe[bin300.glioma_healthy.Dhm.q0.05]%>%
smartLOLA(bin300_universe,regionDB,cores = 10)%>%
as.data.frame(stringsAsFactors = F)%>% within({
    logOR = log2(oddsRatio)
    qvalue = cut(qValue,c(0,1e-50,1e-10,1e-5,0.05),labels = paste0('q<',c(1e-50,1e-10,1e-5,0.05)))
})

description_levels = c('cgi','cpgi_shores','cpgi_shelves','cpgi_opensea','promoter','utr5','coding','introns','utr3','intergenic')

significance_alpha=c(1,0.6,0.3)%>%setNames(c('q<1e-50','q<1e-10','q<1e-05'))

Fig.1D = bin300.glioma_healthy.Dhm.enrichment %>% 
within({
    description=factor(description,description_levels)
})%>%df.split('collection') %>%
lapply(function(enrichment){
    enrichment %>% droplevels %>% ggplot()+
    geom_bar(aes_string(x = 'description',y='logOR',alpha='qvalue'),fill='gold',stat='identity',color='black')+
    scale_alpha_manual(values=significance_alpha)+xlab(NULL)+theme_classic()+
    theme(axis.title = element_text(size = 15),
          axis.text.x = element_text(size = 14,angle = 30,hjust = 1))
})%>% plg(nrow=1)

options(repr.plot.height=3,repr.plot.width=8)
Fig.1D

ggsave2('data/output/Fig.1D.svg',Fig.1D,width = 8,height = 3)

gene.gtf %>% head

genebodyInfo = gene.gtf %>% filter(seqnames %in% paste0('chr',1:22))%>%
dplyr::select(seqnames,start,end,gene_id) %>% make.rownames('gene_id') %>% 
setNames(c('chr','start','end','id'))

genebody.bed =protein_codingAutosomalGenes %>% dplyr::select(seqnames,start,end,gene_id) %>% 
makeGRangesFromDataFrame(keep.extra.columns = T)%>%
sort %>% granges2bed('data/ref/hg19_genebody_simple.bed',keep.mcols = T)

sch('genebody.counts',{
    genebody.counts = rawInfo%>%pull(sampleID)%>% 
    flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
        bed = subset(rawInfo, sampleID == sampleid)$bed_file
        cts = bed %>% routine.coverage(feature.input = genebody.bed,need.cols = 4:5) %>% setNames(sampleid)
        cts
    })%>% merge.df.list %>% trim_rawcts
    genebody.counts
})

genebody.logfpkm = log2(naive_fpkm(
    cts = genebody.counts+1,
    lengths = genebodyInfo[rownames(genebody.counts),]%>% with({end-start}),
    lib.sizes = vlookup(colnames(genebody.counts),rawInfo,'sampleID','lib_size')
))

sch('genebody.glioma_healthy.Dhm',cacheSubDir = NULL,{
    diff_colData = cfDNA_set %>% 
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })

    Dhms = de.edger(counts = genebody.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

sch('genebody.glioma_wbc.Dhm',{
    diff_colData = tissue_set %>% 
    within({
        y = ifelse(label=='wbc','wbc','glioma')%>%factor(c('wbc','glioma'))
    })
    Dhms = de.edger(counts = genebody.counts,colData = diff_colData, 
         coefs = DEcoefs[['glioma']])
    Dhms %>% bind_rows
})

table(genebody.glioma_healthy.Dhm$qvalue <0.05)

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

HTMAP_markers= genebody.glioma_healthy.Dhm %>% arrange(qvalue) %>% head(1000)%>%pull(id)

HTMAP1000 = genebody.logfpkm[HTMAP_markers, HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = HTMAP_samples[, 'y', drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

Fig.2A=HTMAP1000$gtable

options(repr.plot.width=5, repr.plot.height=3) 
grid.draw(Fig.2A)

ggsave2('data/output/Fig.2A.png',Fig.2A,width = 5,height = 3)

genebody.glioma_healthy.gain.q0.05 = genebody.glioma_healthy.Dhm %>% with({id[logFC>0&qvalue<0.05]})

genebody.glioma_healthy.loss.q0.05 = genebody.glioma_healthy.Dhm %>% with({id[logFC<0&qvalue<0.05]})

genebody.glioma_wbc.gain.q0.00001 = genebody.glioma_wbc.Dhm %>% with({ id[logFC>0&qvalue<0.00001]})

genebody.glioma_wbc.loss.q0.00001 = genebody.glioma_wbc.Dhm %>% with({id[logFC<0&qvalue<0.00001]})

length(genebody.glioma_healthy.gain.q0.05)

length(genebody.glioma_wbc.gain.q0.00001)

length(intersect(genebody.glioma_healthy.gain.q0.05,genebody.glioma_wbc.gain.q0.00001))

length(genebody.glioma_healthy.loss.q0.05)

length(genebody.glioma_wbc.loss.q0.00001)

length(intersect(genebody.glioma_healthy.loss.q0.05,genebody.glioma_wbc.loss.q0.00001))

length(intersect(genebody.glioma_healthy.gain.q0.05,genebody.glioma_wbc.loss.q0.00001))

length(intersect(genebody.glioma_healthy.loss.q0.05,genebody.glioma_wbc.gain.q0.00001))

suppressPackageStartupMessages({require(regioneR)})

genebody.universe = makeGRangesFromDataFrame(genebodyInfo,keep.extra.columns = T)

sch('genebody.glioma.perm',recreate=T,{
    c('gain','loss') %>% flexible.mclapply(mc.cores = 2,function(direction){
        overlapPermTest(genebody.universe[get(sprintf('genebody.glioma_healthy.%s.q0.05',direction))], 
                    genebody.universe[get(sprintf('genebody.glioma_wbc.%s.q0.00001',direction))], 
                        ntimes = 1000, genome = genebody.universe)
    })
})

genebody.glioma.perm

Fig.2C = c('gain','loss')%>%lapply(function(direction){
direction_split = genebody.glioma.perm[[direction]]
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
    scale_x_discrete(labels = c(gain='Hyper-5hmC\nGenes',loss='Hypo-5hmC\nGenes'))+
        ylab('Number of Overlaps')+xlab(NULL)+
        theme(axis.title.x = element_text(size = 18),
              axis.text.x = element_text(size = 14),
              strip.text.x = element_text(size = 16),legend.position = 'none')
    })

options(repr.plot.width=4,repr.plot.height=3)
Fig.2C

ggsave2('data/output/Fig.2C.svg',Fig.2C,width = 4,height = 3)



sch('genebody.glioma.reverse',recreate=T,{
    c('gain','loss') %>% flexible.mclapply(mc.cores = 2,function(direction){
        overlapPermTest(genebody.universe[get(sprintf('genebody.glioma_healthy.%s.q0.05',direction))], 
                    genebody.universe[get(sprintf('genebody.glioma_wbc.%s.q0.00001',setdiff(c('gain','loss'),direction)))], 
                        ntimes = 1000, genome = genebody.universe)
    })
})

genebody.glioma.reverse

options(repr.plot.width=4,repr.plot.height=3)
c('gain','loss')%>%lapply(function(direction){
direction_split = genebody.glioma.reverse[[direction]]
obj = direction_split %>% with(numOverlaps)
with(obj,bind_rows(
    data.frame(stringsAsFactors = F, n_overlaps = observed, type = 'observed'),
    data.frame(stringsAsFactors = F, n_overlaps = permuted, type = 'permuted')
    )%>%mutate(direction=direction,pval = pval))})%>%bind_rows %>% df.split('type')%>%with({
        significance = ifelse(observed$pval[1]<0.001,'p < 0.001',as.character(observed$pval[1]))
        permuted %>% ggboxplot(x='direction',y='n_overlaps',color='black')+
        geom_point(data = observed,mapping = aes(x=direction,y=n_overlaps,color = direction),shape=18,stroke=3,size=2)+
        geom_text(data = observed,mapping = aes(x=direction,y=500,label = significance))+
        scale_color_manual(values = c(gain='red',loss='blue'))+
    scale_x_discrete(labels = c(gain='Hyper-5hmC in cfDNA\nHypo-5hmC in gDNA',loss='Hypo-5hmC in cfDNA\nHyper-5hmC in gDNA'))+
        ylab('Number of Overlaps')+xlab(NULL)+
        theme(axis.title.x = element_text(size = 18),
              axis.text.x = element_text(size = 8),
              strip.text.x = element_text(size = 16),legend.position = 'none')
    })


tcga_Info = read_feather('data/ref/blcspg_colData_atomic.feather')

# Download TCGA RNA-seq FPKMUQ data and extract expression matrix

tcga_fpkmuq = tcga_fpkmuq[rowMedians(tcga_fpkmuq)>0,]

tcga_Info = tcga_Info %>% filter(definition=='Primary solid Tumor')%>%within({
    cancer_by_site =ifelse(primary_site=='Brain','Glioma',paste0(primary_site,' cancer'))
    cancer_by_site = cancer_by_site%>%factor(unique(c('Glioma',cancer_by_site)))
})

genebody.overlap.gain = intersect(genebody.glioma_healthy.gain.q0.05,genebody.glioma_wbc.gain.q0.00001)%>%
extract.between('^','\\.') %>% intersect(rownames(tcga_fpkmuq))

genebody.overlap.loss = intersect(genebody.glioma_healthy.loss.q0.05,genebody.glioma_wbc.loss.q0.00001)%>%
extract.between('^','\\.') %>% intersect(rownames(tcga_fpkmuq))

Fig.3A = bind_rows(
    within(tcga_Info,{direction='Hyper-5hmC Genes';
                      E=tcga_fpkmuq[genebody.overlap.gain,tcga_Info$barcode]%>%rowZscore%>%colMeans}),
    within(tcga_Info,{direction='Hypo-5hmC Genes';
                      E=tcga_fpkmuq[genebody.overlap.loss,tcga_Info$barcode]%>%rowZscore%>%colMeans}))%>%
ggplot()+
geom_boxplot(aes_string(x = 'cancer_by_site',y='E',fill = 'direction'),
             color='black',position = position_dodge2(0.9),outlier.size = 0.1)+
facet_wrap(facets = 'direction')+geom_hline(yintercept = 0,linetype='dashed')+
xlab('Tumor Site')+ylab('Mean Expression Z score')+
theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))

options(repr.plot.width=8,repr.plot.height=4)
Fig.3A

ggsave2('data/output/Fig.3A.svg',Fig.3A,width = 8,height = 4)

suppressPackageStartupMessages({
    require(clusterProfiler)
    library(DOSE)
    library(msigdbr)
    require(org.Hs.eg.db)
    require(meshes)
    require(MeSH.Hsa.eg.db)
})

# The msigdb annotation files were prepared using data from msigdbr package
# msigH = msigdbr(species = "Homo sapiens",category = 'H')
# msigH %>% dplyr::select(gs_id,entrez_gene,gs_name) %>% setNames(c('annotID','entrezID','annotTerm')) %>% distinct_all %>%
# write.table('data/ref/DOSE_annot/msigH.txt',quote = T,sep = '\t',row.names = F,col.names = T)

annotationDB_files = list.files('data/ref/DOSE_annot',full.names = T) %>% 
as.list %>%setNames(sub('.txt','',list.files('data/ref/DOSE_annot'),fixed = T))

annotationDB_files %>% names %>% cat

annotationDBs = annotationDB_files %>% lapply(read.table,colClasses = 'character',header=T,sep='\t')

FGM = data.frame(stringsAsFactors = F,id=genebodyInfo$id,ENSEMBL=extract.between(genebodyInfo$id))

id2entrez = merge(FGM , FGM$ENSEMBL %>% search_geneID, by = "ENSEMBL")

universe = id2entrez %>% pull(ENTREZID) %>% unique

genebody.cfDNA_Dhm_GL = id2entrez %>% filter(id %in% with(genebody.glioma_healthy.Dhm,{id[qvalue<0.05]}))

genebody.cfDNA_Dhm_ORA = names(annotationDBs)%>%flexible.mclapply(mc.cores=10,function(annotationDB){
    ORA = genebody.cfDNA_Dhm_GL %>% with({
        ora = myenrich(ENTREZID, annotationDBs[[annotationDB]], universe)
        ora@result = within(ora@result, {
            annotationDB = rep(annotationDB, length(ID))
        })
        ora
    })
})

genebody.cfDNA_Dhm_ORAplots = lapply(genebody.cfDNA_Dhm_ORA,function(ora){
    result=ora@result %>%filter(p.adjust<0.05)
    plot = NULL
    if(nrow(result)>0){
        plot = ora%>% dotplot(showCategory=6) + ggtitle(result$annotationDB[1])
    }
    plot
})
genebody.cfDNA_Dhm_ORAplots=genebody.cfDNA_Dhm_ORAplots[!sapply(genebody.cfDNA_Dhm_ORAplots,is.null)]

options(repr.plot.width=20,repr.plot.height = 5)
genebody.cfDNA_Dhm_ORAplots%>%lapply(function(plot){
    plot+theme(text = element_text(size=20),
               axis.text.y = element_text(size=20),
               axis.text.x = element_text(size=20),
              axis.title.x = element_text(size=20))
})

Fig.3B = genebody.cfDNA_Dhm_ORAplots[[1]]+theme(text = element_text(size=20),
               axis.text.y = element_text(size=20),
               axis.text.x = element_text(size=20),
              axis.title.x = element_text(size=20))

options(repr.plot.width=20,repr.plot.height = 5)
Fig.3B

ggsave2('data/output/Fig.3B.svg',Fig.3B,width = 20,height = 5)



genebody.cfDNA_Dhm_ORA$msigC2_KEGG@result %>% filter(p.adjust<0.05) 

Fig.3C1 = id2entrez%>%filter(SYMBOL=='APC') %>%
as.list %>% with({
    p1 = cfDNA_set %>% within({
        disease = factor(disease,c('GLIOMA','HEALTHY','MENINGIOMA','mNSCLC','HCC','PDAC'))
        H = genebody.logfpkm[id,sampleID] %>% as.numeric
    }) %>% ggplot(aes(x=disease,y=H))+
    geom_boxplot(outlier.size = 0)+geom_jitter(size=0.2)+ stat_compare_means(label.x = 2)+
    ylab(SYMBOL)+xlab('cfDNA 5hmC')+coord_cartesian(clip = 'off')+
    theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))
    p2 = tcga_Info %>% within({
        E = tcga_fpkmuq[ENSEMBL,barcode] %>% as.numeric
    }) %>% ggplot(aes(x=cancer_by_site,y=E))+coord_cartesian(clip = 'off')+
    geom_boxplot(outlier.size = 0)+geom_jitter(size=0.1)+stat_compare_means(label.x = 2)+
    ylab(SYMBOL)+xlab('TCGA RNA')+
    theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))
    plot_grid(p1,p2,nrow=1)
})

options(repr.plot.width=6,repr.plot.height = 3)
Fig.3C1

ggsave2('data/output/Fig.3C1.svg',Fig.3C1,width = 6,height = 3)

Fig.3C2 = id2entrez%>%filter(SYMBOL=='GRIA4') %>%
as.list %>% with({
    p1 = cfDNA_set %>% within({
        disease = factor(disease,c('GLIOMA','HEALTHY','MENINGIOMA','mNSCLC','HCC','PDAC'))
        H = genebody.logfpkm[id,sampleID] %>% as.numeric
    }) %>% ggplot(aes(x=disease,y=H))+
    geom_boxplot(outlier.size = 0)+geom_jitter(size=0.2,alpha=0.3)+ stat_compare_means(label.x = 2)+
    ylab(SYMBOL)+xlab('cfDNA 5hmC')+coord_cartesian(clip = 'off')+
    theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))
    p2 = tcga_Info %>% within({
        E = tcga_fpkmuq[ENSEMBL,barcode] %>% as.numeric
    }) %>% ggplot(aes(x=cancer_by_site,y=E))+coord_cartesian(clip = 'off')+
    geom_boxplot(outlier.size = 0)+geom_jitter(size=0.1,alpha=0.3)+stat_compare_means(label.x = 2)+
    ylab(SYMBOL)+xlab('TCGA RNA')+
    theme_classic()+theme(axis.text.x = element_text(angle = 30,hjust = 1))
    plot_grid(p1,p2,nrow=1)
})

options(repr.plot.width=6,repr.plot.height = 3)
Fig.3C2

ggsave2('data/output/Fig.3C2.svg',Fig.3C2,width = 6,height = 3)

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

sch('genebody.glioma_healthy.Dhm.train',cacheSubDir = NULL,{
    train_diff_colData = cfDNA_set %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })

    Dhms = de.edger(counts = genebody.counts,colData = train_diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

train_gliomaModel_colData = gliomaModel_colData %>% filter(!sampleID %in% heldOut)
test_gliomaModel_colData = gliomaModel_colData %>% filter(sampleID %in% heldOut)

train_gliomaModel_colData$y %>% table

topMarkers.genebody.glioma.train = topMarker_functions[[topMarkers_method]](genebody.glioma_healthy.Dhm.train)

sch('genebody.gliomaModel',recreate=F,cacheSubDir = NULL,{
    set.seed(1234)
    caret_trainer(matrix = t(genebody.logfpkm[topMarkers.genebody.glioma.train, train_gliomaModel_colData$sampleID]), y = train_gliomaModel_colData$y, 
                  model="glmnet",preProcess = c("center", "scale", "pca"),
                  trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
                 )
})

genebody.gliomaModel$preProcess

sch('genebody.gliomaPrediction',recreate=F,cacheSubDir = NULL,reload=T,{
      cfDNA_set %>% within({
      assignment = case_when(sampleID %in% test_gliomaModel_colData$sampleID ~ "testing", 
                             sampleID %in% train_gliomaModel_colData$sampleID ~ 'training',TRUE~'N.A')
      model_RHS = RHS.glmnet(genebody.gliomaModel, t(genebody.logfpkm[topMarkers.genebody.glioma.train, sampleID]))
      probability = predict(genebody.gliomaModel,t(genebody.logfpkm[topMarkers.genebody.glioma.train, sampleID]),type='prob')[['glioma']]
          model_prediction = predict(genebody.gliomaModel, t(genebody.logfpkm[topMarkers.genebody.glioma.train, sampleID]))
      QA = vlookup(sampleID, within(gliomaModel_colData, 
        {
          y = as.character(y)
        }), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment, model_RHS,probability, model_prediction, QA)
})

genebody.gliomaPrediction = genebody.gliomaPrediction%>% 
merge(cfDNA_set,by='sampleID') %>%
mutate(GhmC_Score=model_RHS,
       label = ifelse(QA =='glioma','glioma',label)%>%
       factor(c(label_levels$label,'glioma'))%>%droplevels)

genebody.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({confusionMatrix(model_prediction,factor(QA,levels(model_prediction)),positive = 'glioma')})

suppressPackageStartupMessages({require(pROC)})

genebody.gliomaPrediction.roc = genebody.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({roc(QA,GhmC_Score,levels = c('others','glioma'))})

ci.auc(genebody.gliomaPrediction.roc)%>%round(3)

options(repr.plot.width=5,repr.plot.height=5)
genebody.gliomaPrediction.roc %>% 
plot(print.thres=0,main='Gliomas VS Non-Gliomas',
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

options(repr.plot.width=5,repr.plot.height=5)
svg('data/output/Fig.4B.svg',width = 5,height = 5)
genebody.gliomaPrediction.roc %>% 
plot(print.thres=0,main='Gliomas VS Non-Gliomas',
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)
dev.off()

topMarkers.genebody.top300.glioma.train = topMarker_functions[['top300']](genebody.glioma_healthy.Dhm.train)

sch('genebody.top300.gliomaModel',recreate=T,cacheSubDir = NULL,{
    set.seed(1234)
    caret_trainer(matrix = t(genebody.logfpkm[topMarkers.genebody.top300.glioma.train, train_gliomaModel_colData$sampleID]), y = train_gliomaModel_colData$y, 
                  model="glmnet",preProcess = c("center", "scale"),
                  trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
                 )
})

sch('genebody.top300.gliomaPrediction',recreate=T,cacheSubDir = NULL,reload=T,{
      cfDNA_set %>% within({
      assignment = case_when(sampleID %in% test_gliomaModel_colData$sampleID ~ "testing", 
                             sampleID %in% train_gliomaModel_colData$sampleID ~ 'training',TRUE~'N.A')
      model_RHS = RHS.glmnet(genebody.top300.gliomaModel, t(genebody.logfpkm[topMarkers.genebody.top300.glioma.train, sampleID]))
      probability = predict(genebody.top300.gliomaModel,t(genebody.logfpkm[topMarkers.genebody.top300.glioma.train, sampleID]),type='prob')[['glioma']]
          model_prediction = predict(genebody.top300.gliomaModel, t(genebody.logfpkm[topMarkers.genebody.top300.glioma.train, sampleID]))
      QA = vlookup(sampleID, within(gliomaModel_colData, 
        {
          y = as.character(y)
        }), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment, model_RHS,probability, model_prediction, QA)
})

genebody.top300.gliomaPrediction = genebody.top300.gliomaPrediction%>% 
merge(cfDNA_set,by='sampleID') %>%
mutate(GhmC_Score=model_RHS,
       label = ifelse(QA =='glioma','glioma',label)%>%
       factor(c(label_levels$label,'glioma'))%>%droplevels)

options(repr.plot.width=5,repr.plot.height=5)
genebody.top300.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({confusionMatrix(model_prediction,factor(QA,levels(model_prediction)),positive = 'glioma')})

options(repr.plot.width=5,repr.plot.height=5)
genebody.top300.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({roc(QA,GhmC_Score,levels = c('others','glioma'))}) %>% 
plot(print.thres=0,main='Gliomas VS Non-Gliomas using Top 300 DhmGs',cex.main = 0.9,
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

options(repr.plot.width=5,repr.plot.height=5)
genebody.top300.gliomaPrediction %>% 
filter(assignment=='testing')%>%
with({roc(QA,GhmC_Score,levels = c('others','glioma'))}) %>% 
ci.auc %>% round(3)

immediatePosBelow = function(q,x){
    if(all(x>q)){p=0}else{p=(1:length(x))[x<=q][which.max(x[x<=q])]}
    p
}

getCutoffs = function(y,r,q){
    y[sapply(q,immediatePosBelow,r)]
}

Fig.4C = genebody.gliomaPrediction %>% filter(assignment=='testing') %>%
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

options(repr.plot.width=8,repr.plot.height=3.5)
Fig.4C

ggsave2('data/output/Fig.4C.svg',Fig.4C,height = 3.5,width = 8)

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction %>% filter(assignment=='testing')%>% 
ggboxplot(x='gender',y='GhmC_Score')+stat_compare_means()+xlab('Sex')

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction %>% filter(assignment=='testing')%>% 
ggboxplot(x='age_group',y='GhmC_Score')+stat_compare_means()+xlab('Age Group')

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction %>% filter(assignment=='testing')%>% 
filter(pathologyAbbr %in% c('GBM','DA','O2','AA3'))%>%ggboxplot(x='pathologyAbbr',y='GhmC_Score')+
stat_compare_means()+xlab('Pathology')

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction %>% filter(assignment=='testing')%>% 
filter(grade !='N.A')%>%mutate(grade=factor(grade,paste0('grade',1:4)))%>%
ggboxplot(x='grade',y='GhmC_Score')+stat_compare_means()+xlab('Grade')

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction %>% filter(assignment=='testing')%>% 
filter(idh %in% subtypes$idh)%>%ggboxplot(x='idh',y='GhmC_Score')+stat_compare_means()+xlab('IDH status')

options(repr.plot.width=5,repr.plot.height=3)
genebody.gliomaPrediction%>%
filter(assignment=='testing')%>% 
filter(!is.na(enhancement))%>% 
mutate(enhancement=ifelse(enhancement=='NonEnhancing','NonEnhancing','Enhancing'))%>%
ggboxplot(x='enhancement',y='GhmC_Score')+stat_compare_means()+xlab('Contrast Enhancement')#+
# theme(axis.text.x = element_text(angle = 15,hjust=1))

options(repr.plot.width=3,repr.plot.height=3)
genebody.gliomaPrediction %>%
filter(assignment=='testing'&!is.na(max_diameter))%>%
dplyr::select(sampleID,patientID,probability,pathologyAbbr,idh,GhmC_Score,enhancement,max_diameter)%>%
(function(data){
    plm = summary(lm(GhmC_Score~max_diameter,data))$coefficients['max_diameter','Pr(>|t|)']
    ggplot(data,aes_string(x='max_diameter',y='GhmC_Score'))+
geom_smooth(method = "lm", formula = y ~ x,color='black',linetype='dashed',alpha=0.1)+
geom_point(color='black')+
    annotate(geom = 'text',y = 0.65,x=6.5,label=sprintf('p = %.3f',plm))+
scale_color_manual(values = my_palette$grade)+theme_classic()})

# patientsWithMRI=dget('data/patientsWithMRI.R')

sch('genebody.glioma_healthy.Dhm.noMRI',recreate=T,cacheSubDir = NULL,{
    noMRI_diff_colData = cfDNA_set %>% 
    filter(!patientID %in% patientsWithMRI)%>%
    filter(label %in% c('gbm','nongbm','healthy')) %>% 
    within({
        y = ifelse(label=='healthy','healthy','glioma')%>%factor(c('healthy','glioma'))
    })

    Dhms = de.edger(counts = genebody.counts,colData = noMRI_diff_colData, 
         coefs = DEcoefs[['glioma']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

topMarkers.genbody.glioma.noMRI = topMarker_functions[[topMarkers_method]](genebody.glioma_healthy.Dhm.noMRI)

noMRI_gliomaModel_colData = gliomaModel_colData %>% filter(!patientID %in% patientsWithMRI)
withMRI_gliomaModel_colData = gliomaModel_colData %>% filter(patientID %in% patientsWithMRI)

sch('genebody.gliomaModel_MRI',recreate=T,cacheSubDir = NULL,{
    set.seed(1234)
    caret_trainer(matrix = t(genebody.logfpkm[topMarkers.genbody.glioma.noMRI, noMRI_gliomaModel_colData$sampleID]), y = noMRI_gliomaModel_colData$y, 
                  model="glmnet",preProcess = c("center", "scale", "pca"),
                  trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
                 )
})

sch('genebody.gliomaPrediction_MRI',recreate=T,cacheSubDir = NULL,reload=T,{
      cfDNA_set %>% within({
      assignment = case_when(sampleID %in% withMRI_gliomaModel_colData$sampleID ~ "testing", 
                             sampleID %in% noMRI_gliomaModel_colData$sampleID ~ 'training',TRUE~'N.A')
      model_RHS = RHS.glmnet(genebody.gliomaModel_MRI, t(genebody.logfpkm[topMarkers.genbody.glioma.noMRI, sampleID]))
      probability = predict(genebody.gliomaModel_MRI,t(genebody.logfpkm[topMarkers.genbody.glioma.noMRI, sampleID]),type='prob')[['glioma']]
          model_prediction = predict(genebody.gliomaModel_MRI, t(genebody.logfpkm[topMarkers.genbody.glioma.noMRI, sampleID]))
      QA = vlookup(sampleID, within(gliomaModel_colData, 
        {
          y = as.character(y)
        }), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment, model_RHS,probability, model_prediction, QA)
})

Fig.4D = genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing'&!is.na(max_diameter))%>% mutate(GhmC_Score=model_RHS)%>%
dplyr::select(sampleID,patientID,probability,pathologyAbbr,idh,GhmC_Score,enhancement,max_diameter)%>%
(function(data){
    plm = summary(lm(GhmC_Score~max_diameter,data))$coefficients['max_diameter','Pr(>|t|)']
    ggplot(data,aes_string(x='max_diameter',y='GhmC_Score'))+
geom_smooth(method = "lm", formula = y ~ x,linetype='dashed',alpha=0.1)+
geom_point(color='black')+
    annotate(geom = 'text',y = 0.65,x=8,label=sprintf('p = %.3f',plm))+
scale_color_manual(values = my_palette$grade)+theme_classic()})

options(repr.plot.width=3,repr.plot.height=3)
Fig.4D

ggsave2('data/output/Fig.4D.svg',Fig.4D,height = 3,width = 3)

options(repr.plot.width=3,repr.plot.height=3)
genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing'&!is.na(max_diameter)&enhancement%in%c("Weakly_Enhancing","Strongly_Enhancing"))%>%
mutate(GhmC_Score=model_RHS)%>%
dplyr::select(sampleID,patientID,probability,pathologyAbbr,idh,GhmC_Score,enhancement,max_diameter)%>%
(function(data){
    plm = summary(lm(GhmC_Score~max_diameter,data))$coefficients['max_diameter','Pr(>|t|)']
    ggplot(data,aes_string(x='max_diameter',y='GhmC_Score'))+
geom_smooth(method = "lm", formula = y ~ x,linetype='dashed',alpha=0.1)+
geom_point(color='black')+
    annotate(geom = 'text',y = 0.65,x=8,label=sprintf('p = %.3f',plm))+
scale_color_manual(values = my_palette$grade)+theme_classic()})

options(repr.plot.width=4,repr.plot.height=3)
genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing')%>%filter(!is.na(enhancement))%>% 
mutate(GhmC_Score=model_RHS)%>%
mutate(enhancement=ifelse(enhancement=='NonEnhancing','NonEnhancing','Enhancing'))%>%
ggboxplot(x='enhancement',y='GhmC_Score')+stat_compare_means()+xlab('Contrast Enhancement')#+
# theme(axis.text.x = element_text(angle = 15,hjust=1))

options(repr.plot.width=3,repr.plot.height=3)
genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing'&!is.na(max_diameter)&enhancement%in%c("NonEnhancing"))%>%
mutate(GhmC_Score=model_RHS)%>%
dplyr::select(sampleID,patientID,probability,pathologyAbbr,idh,GhmC_Score,enhancement,max_diameter)%>%
(function(data){
    plm = summary(lm(GhmC_Score~max_diameter,data))$coefficients['max_diameter','Pr(>|t|)']
    ggplot(data,aes_string(x='max_diameter',y='GhmC_Score'))+
geom_smooth(method = "lm", formula = y ~ x,linetype='dashed',alpha=0.1)+
geom_point(color='black')+
    annotate(geom = 'text',y = 0.65,x=8,label=sprintf('p = %.3f',plm))+
scale_color_manual(values = my_palette$grade)+theme_classic()})

genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing'&!is.na(max_diameter))%>%
mutate(GhmC_Score=model_RHS)%>%(function(data){
    with(data, {
        list(cor.test(GhmC_Score,max_diameter),
             summary(lm(GhmC_Score~max_diameter,data))$coefficients
            )
    })
})

genebody.gliomaPrediction_MRI %>% merge(cfDNA_set,by='sampleID') %>% 
filter(assignment=='testing'&!is.na(max_diameter)&enhancement%in%c("Weakly_Enhancing","Strongly_Enhancing"))%>%
mutate(GhmC_Score=model_RHS)%>%(function(data){
    with(data, {
        list(cor.test(GhmC_Score,max_diameter),
             summary(lm(GhmC_Score~max_diameter,data))$coefficients
            )
    })
})

cfDNA_set_idh = cfDNA_set %>% filter(idh %in% subtypes[['idh']]|label =='healthy')

idhModel_colData = cfDNA_set_idh %>% 
filter(idh %in% nonAmbiguous[['idh']]) %>% 
within({
    y = factor(idh,nonAmbiguous[['idh']])
    age_scaled=(age-mean(age))/sd(age)
})

sch('genebody.idh.tissue.Dhm',recreate=F,{
    de.edger(counts = genebody.counts,
             colData = tissue_set %>% filter(idh %in% nonAmbiguous[['idh']])%>%
             within({
                y = factor(idh,nonAmbiguous[['idh']])
            }), coefs = list(IDHmut=2))%>%bind_rows
})

genebody.idh.tissue.Dhm.q0.05 = processors$cutBy(genebody.idh.tissue.Dhm,cutoff = 0.05, v = "qvalue")

genebody.idh.tissue.Dhm.q0.00001 = processors$cutBy(genebody.idh.tissue.Dhm,cutoff = 0.00001, v = "qvalue")

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

IDH_HTMAP_markers= genebody.idh.tissue.Dhm.q0.05

IDH_HTMAP = genebody.logfpkm[IDH_HTMAP_markers, IDH_HTMAP_samples$sampleID] %>% 
rowZscore %>% mat_order %>% outlier.trimmer4heatmap(1.5) %>% 
  pheatmap(cluster_cols = F, cluster_rows = F, 
    silent = T, legend = F, #annotation_legend = F, 
    annotation_col = IDH_HTMAP_samples[, c('y'), drop = F], 
    show_rownames = F, show_colnames = F, 
    annotation_colors = IDH_HTMAP_palette, 
   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100))

Fig.5A = IDH_HTMAP$gtable

options(repr.plot.width=5, repr.plot.height=3) 
grid.draw(Fig.5A)

ggsave2('data/output/Fig.5A.png',Fig.5A,width = 5,height = 3)

sch('genebody.idh.cfDNA.Dhm',recreate=F,reload = T,{
    train_diff_colData = cfDNA_set_idh %>% 
    filter(!sampleID %in% heldOut)%>%
    filter(idh %in% c('healthy',nonAmbiguous[['idh']])) %>% 
    within({
        y = factor(idh,c('healthy',nonAmbiguous[['idh']]))
    })
    
    Dhms = de.edger(counts = genebody.counts,colData = train_diff_colData, 
             coefs = DEcoefs[['idh']], extra = c("age_group", "gender"))
    Dhms %>% bind_rows
})

genebody.idh.cfDNA.Dhm = df.split(genebody.idh.cfDNA.Dhm,'label')[c(nonAmbiguous[['idh']],paste0(nonAmbiguous[['idh']][2],'VS',nonAmbiguous[['idh']][1]))]

table(genebody.idh.cfDNA.Dhm$IDHmutVSIDHwt$qvalue <0.05 )

train_idhModel_colData = idhModel_colData %>% filter(!sampleID %in% heldOut)
test_idhModel_colData = idhModel_colData %>% filter(sampleID %in% heldOut)

genebody.idh.train_mat = t(genebody.logfpkm[genebody.idh.tissue.Dhm.q0.05, train_idhModel_colData$sampleID])

genebody.idh.train_preP = preProcess(genebody.idh.train_mat,method = c("center", "scale", "pca"))

genebody.idh.train_PCmat_age = cbind(predict(genebody.idh.train_preP,genebody.idh.train_mat),age=train_idhModel_colData$age_scaled)

genebody.idh.test_mat = t(genebody.logfpkm[genebody.idh.tissue.Dhm.q0.05, test_idhModel_colData$sampleID])

genebody.idh.test_PCmat_age = cbind(predict(genebody.idh.train_preP,genebody.idh.test_mat),age=test_idhModel_colData$age_scaled)

sch('genebody.idhModel',recreate=F,cacheSubDir = NULL,{
set.seed(1234)
    caret_trainer(
        matrix = genebody.idh.train_PCmat_age, y = train_idhModel_colData$y, 
        model="glmnet",trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
    )
})

sch('genebody.idhPrediction',recreate=F,cacheSubDir = NULL,{
   test_idhModel_colData %>% filter(idh %in% subtypes[['idh']]) %>% within({
       probability = predict(genebody.idhModel, genebody.idh.test_PCmat_age,type='prob')[[nonAmbiguous[['idh']][2]]]
      assignment = "testing"
      model_RHS = RHS.glmnet(genebody.idhModel, genebody.idh.test_PCmat_age)
      model_prediction = predict(genebody.idhModel, genebody.idh.test_PCmat_age)
      QA = vlookup(sampleID, within(idhModel_colData, {y = as.character(y)}), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment,probability, model_RHS,  model_prediction, QA)
})

genebody.idhPrediction = genebody.idhPrediction %>% 
merge(cfDNA_set,by='sampleID')

genebody.idhPrediction %>%
with({confusionMatrix(model_prediction,factor(QA,levels(model_prediction)),positive = 'IDHmut')})

genebody.idhPrediction.roc = genebody.idhPrediction %>%
filter(assignment=='testing')%>% 
with({roc(QA,model_RHS,levels = nonAmbiguous[['idh']])})

genebody.idhPrediction.roc %>% ci.auc %>% round(3)

options(repr.plot.width=5,repr.plot.height=5)
genebody.idhPrediction.roc %>% 
plot(print.thres=0,main=paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)

svg('data/output/Fig.5B.svg',width = 5,height = 5)
genebody.idhPrediction.roc %>% 
plot(print.thres=0,main=paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),
     print.thres.pattern="Cutoff = %.0f \nSpecificity = %.3f \nSensitivity = %.3f",
     print.auc=T,print.auc.x=0.3,print.auc.y=0.1)
dev.off()

idh_age_roc = cfDNA_set_idh %>% filter(idh %in% nonAmbiguous[['idh']]) %>% 
with({roc(idh,age,levels = nonAmbiguous[['idh']])}) 

set.seed(1234)
idh_age_roc.test = roc.test(idh_age_roc,genebody.idhPrediction.roc)

options(repr.plot.width=5,repr.plot.height=5)
idh_age_roc%>%
plot(main=paste0(rev(nonAmbiguous[['idh']]),collapse=' VS '),
     print.auc=T,print.auc.x=0,print.auc.y=0.13,print.auc.pattern='AUC using Age Alone: %.3f',print.auc.adj=1)
plot(genebody.idhPrediction.roc,add=T,
     print.auc=T,print.auc.x=0,print.auc.y=0.2,print.auc.pattern='AUC using Model: %.3f',print.auc.adj=1)
text(0,0.06,sprintf('p = %.3f',idh_age_roc.test$p.value),adj = 1)

caret_varImps2 = function (coefs,maxABScoef = max(abs(coefs))) 
{
    100 * abs(coefs)/maxABScoef
}

genebody.idhWX = row_wilcoxon_twosample(genebody.logfpkm[,idhModel_colData$sampleID],g=idhModel_colData$y)

genebody.idhModel_coefs = (genebody.idh.train_preP$rotation %*% coef.pc_glmnet(genebody.idhModel)[colnames(genebody.idh.train_preP$rotation)])[, 1]

genebody.idhModel_varImps = genebody.idhModel_coefs %>% caret_varImps2(max(abs(coef.pc_glmnet(genebody.idhModel))))

genebody.idhModel_terms = id2entrez %>% filter(id %in% names(genebody.idhModel_coefs))%>%
within({coef=genebody.idhModel_coefs[id];imp=genebody.idhModel_varImps[id];
        wxp=genebody.idhWX[id,'pvalue'];wxq=qvalue(wxp)$qvalue})

genebody.idhModel_terms %>%filter(wxq<0.05)%>% arrange(desc(imp)) %>% head

Fig.5C=genebody.idhModel_terms %>%filter(wxq<0.05)%>% df.rowlist %>% lapply(function(row){
    d = idhModel_colData
    id=row$id
    symbol = row$SYMBOL
     d[[id]] = as.numeric(genebody.logfpkm[id,d$sampleID])
    pp = max(d[[id]])+0.1
    q = paste0('q = ',small_round(row$wxq))
    ggboxplot(d,x = 'idh',y=id,outlier.shape = NA)+geom_jitter(aes(color=idh))+
    scale_color_manual(values = my_palette[['idh']])+geom_text(aes(x=1.5,y=pp,label=q),family='sans',fontface='plain')+
    xlab(symbol)+ylab('cfDNA 5hmC logFPKM')+coord_cartesian(clip = 'off')+theme(legend.position = 'none')
}) %>% plg(nrow=1)

options(repr.plot.height=3,repr.plot.width=5)
Fig.5C

ggsave2('data/output/Fig.5C.svg',Fig.5C,height = 3.5,width = 8)



cytobandpq.bed = 'data/ref/cytobandpq.bed'

cytobandpq.gr = bed2granges(cytobandpq.bed)

cfDNA_set_codel = cfDNA_set %>% filter(codel %in% nonAmbiguous[['codel']])

cytobandpq.cts = cfDNA_set_codel %>%
pull(sampleID)%>% 
flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
    bed = subset(cfDNA_set_codel, sampleID == sampleid)$bed_file
    cts = bed %>% routine.coverage(feature.input = cytobandpq.bed,need.cols = 4:5) %>% setNames(sampleid)
    cts
})%>% merge.df.list

cytobandpq.fpkm = naive_fpkm(cts = cytobandpq.cts, lengths = width(cytobandpq.gr),cfDNA_set_codel$lib_size)

cytobandpq.gr_1p19q=cytobandpq.gr[grepl('chr1p',cytobandpq.gr$id)|grepl('chr19q',cytobandpq.gr$id)]

options(repr.plot.width=5,repr.plot.height=3)
cytobandpq.fpkm[cytobandpq.gr_1p19q$id,] %>% keep.rownames('chr')%>% melt(id.vars = 'chr',variable.name = 'sampleID') %>% 
merge(cfDNA_set_codel,by='sampleID') %>% 
mutate(chr=factor(chr,cytobandpq.gr_1p19q$id),codel=factor(codel,nonAmbiguous[['codel']]))  %>% arrange(codel)%>%
ggboxplot('codel','value',facet.by = 'chr')+stat_compare_means(label.y=0.54)+ylim(0.4,0.55)+
ylab('5hmC FPKM')+coord_cartesian(clip='off')

mgmtPromoter.bed = 'data/ref/MGMTpromoter.bed'

cfDNA_set_mgmt=cfDNA_set %>% filter(mgmt %in% nonAmbiguous[['mgmt']] &idh =='IDHwt') %>% 
mutate(mgmtHMC = genebody.logfpkm['ENSG00000170430.9',sampleID] %>% as.numeric)

mgmtPromoter.cts = cfDNA_set_mgmt %>% 
pull(sampleID)%>%
flexible.mclapply(mc.cores = 200,FUN = function(sampleid){
    bed = subset(cfDNA_set_mgmt, sampleID == sampleid)$bed_file
    cts = bed %>% routine.coverage(feature.input = mgmtPromoter.bed,need.cols = 4:5) %>% setNames(sampleid)
    cts
})%>% merge.df.list

cfDNA_set_mgmt$mgmtpHMC = cfDNA_set_mgmt %>% with({as.numeric(mgmtPromoter.cts[sampleID])*1e6/lib_size})

options(repr.plot.width=4,repr.plot.height=3)
cfDNA_set_mgmt %>% 
ggboxplot('mgmt','mgmtHMC',outlier.shape = NA)+ 
geom_jitter(aes(x=mgmt,y=mgmtHMC,color=mgmt),width = 0.1)+
stat_compare_means()+coord_cartesian(clip='off')+ylab('MGMT Gene 5hmC')+xlab(NULL)+theme(legend.position='none')

options(repr.plot.width=4,repr.plot.height=3)
cfDNA_set_mgmt %>% 
ggboxplot('mgmt','mgmtpHMC',outlier.shape = NA)+
geom_jitter(aes(x=mgmt,y=mgmtpHMC,color=mgmt),width = 0.1)+
stat_compare_means()+coord_cartesian(clip='off')+ylab('MGMT Promoter 5hmC')+xlab(NULL)+theme(legend.position='none')



noMRI_idhModel_colData = idhModel_colData %>% filter(!patientID %in% patientsWithMRI)
withMRI_idhModel_colData = idhModel_colData %>% filter(patientID %in% patientsWithMRI)

genebody.idh.noMRI_mat = t(genebody.logfpkm[genebody.idh.tissue.Dhm.q0.05, noMRI_idhModel_colData$sampleID])

genebody.idh.noMRI_preP = preProcess(genebody.idh.noMRI_mat,method = c("center", "scale", "pca"))

genebody.idh.noMRI_PCmat_age = cbind(predict(genebody.idh.noMRI_preP,genebody.idh.noMRI_mat),age=noMRI_idhModel_colData$age_scaled)

genebody.idh.withMRI_mat = t(genebody.logfpkm[genebody.idh.tissue.Dhm.q0.05, withMRI_idhModel_colData$sampleID])

genebody.idh.withMRI_PCmat_age = cbind(predict(genebody.idh.noMRI_preP,genebody.idh.withMRI_mat),age=withMRI_idhModel_colData$age_scaled)

sch('genebody.idhModel_MRI',recreate=F,cacheSubDir = NULL,{
set.seed(1234)
    caret_trainer(
        matrix = genebody.idh.noMRI_PCmat_age, y = noMRI_idhModel_colData$y, 
        model="glmnet",trControl = custom_trainControl, tuneLength = 10,tuneGrid=custom_tuneGrid
    )
})

sch('genebody.idhPrediction_MRI',recreate=T,cacheSubDir = NULL,{
   withMRI_idhModel_colData %>% filter(idh %in% subtypes[['idh']]) %>% within({
       probability = predict(genebody.idhModel_MRI, genebody.idh.withMRI_PCmat_age,type='prob')[[nonAmbiguous[['idh']][2]]]
      assignment = "testing"
      model_RHS = RHS.glmnet(genebody.idhModel_MRI, genebody.idh.withMRI_PCmat_age)
      model_prediction = predict(genebody.idhModel_MRI, genebody.idh.withMRI_PCmat_age)
      QA = vlookup(sampleID, within(idhModel_colData, {y = as.character(y)}), "sampleID", "y", default = "N.A")
    }) %>% dplyr::select(sampleID, assignment,probability, model_RHS,  model_prediction, QA)
})
