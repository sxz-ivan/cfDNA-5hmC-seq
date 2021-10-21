suppressPackageStartupMessages({
    require(GenomicRanges)
    require(LOLA)
    require(annotatr)
})
source('code/vector.tools.r')

smartLOLA = function (userSets, userUniverse, regionDB, minOverlap = 1, cores = 1, 
    redefineUserSets = FALSE,ignoreQ=F){
    rawLOLA = runLOLA(userSets, userUniverse, regionDB, minOverlap = minOverlap, cores = cores, 
    redefineUserSets = FALSE)%>% as.data.frame(stringsAsFactors = F)
    depletionPvalue = rawLOLA %>% apply(1,function(x){
        fisher.test(matrix(x[c('support', 'b', 'c', 'd')]%>%as.numeric,ncol = 2,byrow = T), alternative = 'less')$p.value
    })
    finalLOLA = rawLOLA %>% within({
        pValueLog = ifelse(oddsRatio>1,pValueLog,log(depletionPvalue,0.1))
        logOR = log2(oddsRatio)
    })
    if(ignoreQ){
        finalLOLA = finalLOLA %>% within({qValue = NULL})
    }else{
        depletionQvalue = qvalue::qvalue(depletionPvalue)$qvalue
        finalLOLA = finalLOLA %>% within({qValue = ifelse(oddsRatio>1,qValue,depletionQvalue)})
    }
    finalLOLA
}

LOLADB_from_gr = function(named_list_grs,collection = stringi::stri_rand_strings(1,8,'[a-z]'),
                          description = names(named_list_grs),
                          cellType = NA, tissue=NA, antibody=NA, treatment=NA, dataSource=NA, filename=NA){
    suppressPackageStartupMessages({
        require(data.table)
    })
    if(is.null(description)){
        description = paste0('gr',1:length(named_list_grs))
    }
    list()%>% within({
        dbLocation = ''
        regionAnno = data.frame(
            stringsAsFactors = F,
            size = sapply(named_list_grs,length),
            dbSet=1:length(named_list_grs),
            description= description,
            collection = collection,
            cellType = cellType, 
            tissue=tissue,
            antibody=antibody,
            treatment=treatment,
            dataSource=dataSource,
            filename=filename
        )%>%data.table(keep.rownames = F,stringsAsFactors = F)
        collectionAnno = data.frame(stringsAsFactors = F,
            collectionname = collection
        )%>%data.table
        regionGRL = GRangesList(named_list_grs%>% setNames(NULL),compress = T)
    })
}

make_LOLADB = function(bed_files, DB_path, setup_method = 'ln -sf', collection = stringi::stri_rand_strings(1,8,'[a-z]'),
                          description = extract.between(bed_files,'/','.'),
                          cellType = '', dataSource=''){
    region_path = paste0(DB_path,'/',collection,'/regions')
    dir.create(region_path,recursive = T)
    file_names = bed_files %>% extract.between('/','$')
    dest_names = paste0(region_path,'/',file_names)
    setup_cmds = paste0(setup_method,' ',bed_files,' ',dest_names)
    setup_cmds %>% lapply(system)
    index_file = paste0(DB_path,'/',collection,'/index.txt')
    index_df = data.frame(
        stringsAsFactors = F,
        filename=file_names,
        description= description,
        cellType = cellType, 
        dataSource=dataSource
    )%>% write.table(index_file,row.names = F,col.names = T,quote = F,sep = '\t')
    DB_path
}

annotate_regions=function (regions, annotations, minoverlap = 1L, ignore.strand = TRUE, 
    quiet = FALSE) 
{
    if (class(regions)[1] != "GRanges") {
        stop("Error in annotate_regions(...): regions object is not GRanges.")
    }
    if (class(annotations)[1] != "GRanges") {
        stop("Error in annotate_regions(...): annotations object is not GRanges. Use build_annotations(...) to construct the annotations before calling annotate_regions(...).")
    }
    if (!quiet) {
        message("Annotating...")
    }
    intersections = GenomicRanges::findOverlaps(regions, annotations, 
        minoverlap = minoverlap, ignore.strand = ignore.strand)
#     if (length(intersections) > 0) {
        gr = regions[S4Vectors::queryHits(intersections)]
        GenomicRanges::mcols(gr)$annot = annotations[S4Vectors::subjectHits(intersections)]
        return(gr)
#     }
#     else {
#         stop("No annotations intersect the regions.")
#     }
}