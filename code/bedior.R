suppressPackageStartupMessages({
  require(dplyr)
  require(reshape2)
  require(stringi)
  require(parallel)
  require(stringr)
    require(GenomicRanges)
})

bedtools = 'bedtools'
bedtools.genome.file = 'data/ref/human_simple.genome'
bedtools.genome.bed = 'data/ref/human_simple.bed'
igvtools = 'igvtools'

source('code/vector.tools.r')

shell.io = function(input, command, output='R',pipe.in=input%>% determine.input(),execute = T,stdin = '-'){
  if(pipe.in){
    shell = paste0(input, command,' ',stdin,' ')
  }else{
    shell = paste0(command,' ',input,' ')
  }
  
  if(output =='pipe'){
    paste0(shell,' | ')
  }else{
    if(output == 'R'){
      if(execute){system(shell,intern=T)}else{shell}
    }else{
      shell = paste0(shell, ' > ', output )
      if(execute){
        system(shell)
        output
      }else{
        cat(shell)
        invisible(shell)
      }
    }
  }
}
delete = function(file,command = 'rm ',...){
  shell.io(file, command = command, ...)
}
determine.input = function(input){
  length = nchar(input)
  substr(input, length-2, length) == ' | '
}
append.process.tag = function(file ,tag = '',postfix ='.bed',replace = F){
  if(replace){replacement = paste0(tag, '')}else{replacement = paste0(tag, postfix)}
  sub(pattern = postfix, replacement = replacement,file,fixed = T)
}

awk.filter = function(input,
                      pipe.in=input%>% determine.input(),
                      filter = '',
					  condition.out = '',
					  need.cols = 0,
					  OFS = '\t',
                      execute = T,...){
  #only use double where there is single, use single in other places to avoid mis-interpretation
  need.cols = paste('$',need.cols,sep = '',collapse = ',')
  awk.shell = paste0("awk '",filter, ' { ', condition.out ,' print ',need.cols,' } ', "' OFS='",OFS,"'")
  shell.io(input, awk.shell, pipe.in = pipe.in,execute = execute,stdin = '-',...)
}

paste.bed_id = function(input, id_style = '$1,":",$2,"-",$3',...){
  awk.filter(input = input, need.cols = paste0('1,"\t",$2,"\t",$3,"\t",',id_style),OFS = "",...)
}

read.coverage = function(input, file = T, input.cols = 4:5,
                         full.fields = c('chr','start','end','id',
                                         'count','n.base.covered.A','length.A','fraction.base.covered.A'),
                         want.fields = c('id','count'),make.id.rownames = F,dedup = T
){
  full.fields = full.fields[input.cols]
  if(file){coverage = read.table(input,sep = '\t',col.names = full.fields, colClasses = c('character','numeric','numeric','character',rep('numeric',4)))
  }else{
    coverage = colsplit(string = input, pattern = '\t', names = full.fields)
  }
  coverage = coverage[,want.fields]
  if(make.id.rownames&('id' %in% want.fields)){
	  coverage = coverage[!(duplicated(coverage$id)),]
	  rownames(coverage) = coverage$id
	  coverage$id = NULL
  }else{
  if(dedup){coverage = coverage[!(duplicated(coverage$id)),]}}
  coverage
}

routine.coverage = function(bed, feature.input,need.cols = c(4:5)){
  bed %>%
    bedtools.coverage(feature.input = feature.input,output ='pipe')%>% 
    awk.filter(need.cols = need.cols,output = 'R') %>%
    read.coverage(file = F,make.id.rownames = T)
}

rmdup = function(input, output = NULL,...){
  temp = shell.io(input, command = 'uniq ', output = paste0(input,'.rmdup'),...)
  if(is.null(output)){
    file.remove(input)
    file.rename(temp, input)
    output = input
  }else{
      file.rename(temp,output)
  }
  output
}

get.libsize = function(bed , method = 'wc'){
	cmd = c(wc= 'wc -l ', sed = "sed -n '$=' ")[method]
	bed %>% 
	shell.io(cmd,output = 'R') %>% 
	strsplit.extract(position = 1, split = ' ') %>% 
	as.numeric
}

get.fragment.lengths = function(input.bed, sample.size = 1e5, ...){
	if(is.null(sample.size)){
	lengths = input.bed %>% awk.filter(need.cols = '3-$2', output = 'R',...) %>% as.numeric
	}else{
	  sample.lines(input.bed,
	               n.line = sample.size, 
	               output = 'pipe')%>%
	    lengths = awk.filter(need.cols = '3-$2', output = 'R',...) %>% as.numeric
	}
  lengths
}

bedtools.shell = function(input, 
                          command = 'cat ', arguments = ' ',
                          pipe.in=input%>% determine.input(), process.tag = '',output = input %>% append.process.tag(process.tag),
                          show.head = F,execute = T){
  bedtools.command = sprintf(command, arguments)
  shell.io(input, bedtools.command, output,pipe.in,execute,stdin = 'stdin')
}

bedtools.sort = function(input, 
                         command = paste0(bedtools, ' sort %s -i '),
                         arguments = paste0('-g ',bedtools.genome.file),
                         process.tag = '.sorted',...){
  bedtools.shell(input, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.merge = function(input, 
                          command = paste0(bedtools, ' merge %s -i '),
                          process.tag = '.merged',...){
  bedtools.shell(input, command,process.tag = process.tag, ...)
}

bedtools.slop = function(input, 
                         command = paste0(bedtools, ' slop -g ',bedtools.genome.file,' %s -i '),
                         arguments = '',
                         process.tag = '.slopped',...){
  bedtools.shell(input, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.makewindows = function(input, window_size = 2000,
                         command = paste0(bedtools, ' makewindows -w ',window_size ,' %s -b '),
                         arguments = '',
                         process.tag = '.windows',...){
  bedtools.shell(input, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.intersect = function(input.a,input.b,  
                              command = paste0(bedtools, ' intersect %s -b ', input.b, ' -a '),
                              arguments = ' ',
                              process.tag = '.intersected',...){
							  #其实只有-wa 有用，-wb好像还是会把-a写出来
  bedtools.shell(input.a, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.jaccard = function(input.b,input.a,  
                              command = paste0(bedtools, ' jaccard %s -a ', input.a, ' -b '),
                              arguments = ' ',number.only = T,output = 'R',
                              process.tag = '.jaccard',...){
  result = bedtools.shell(input.b, command, arguments = arguments ,process.tag = process.tag, ...)
  if(number.only){
	result = result %>%(function(x){strsplit(x[2], split = '\t')%>% unlist %>% as.numeric%>%.[3]})
  }
  result
}

#https://github.com/jokergoo/cotools/blob/master/R/genomic_region_correlation.R
genomicCorr.jaccard = function(query, reference, restrict = NULL) {
	if(is.null(restrict)) {
		res = sum(width(intersect(query, reference))) / sum(as.numeric(width(union(query, reference))))
	} else {
		gr1 = intersect(query, reference)
		gr1 = intersect(gr1, restrict)

		gr2 = union(query, reference)
		gr2 = intersect(gr2, restrict)
		res = sum(width(gr1)) / sum(width(gr2))
	}
	return(res)
}

cotools.jaccard = function(input.a, input.b,restrict = NULL){
    granges.a = shell.io(input.a, output = 'R',command = 'cat') %>% bed.string.2granges
    granges.b = shell.io(input.b, output = 'R',command = 'cat') %>% bed.string.2granges
    genomicCorr.jaccard(granges.a,granges.b,restrict = restrict)
}

bedtools.coverage = function(fragments.input, 
                             feature.input, 
                             command = paste0(bedtools,' coverage %s -a ',feature.input,' -b '),
                             arguments = paste0(' -counts -sorted -g ',bedtools.genome.file),process.tag = '.coverage',...){
  bedtools.shell(fragments.input, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.closest = function(input.a, input.b, 
                              command = paste0(bedtools, ' closest %s -b ', input.b, ' -a '),
                              arguments = ' ',
                              process.tag = '.closest',...){
  bedtools.shell(input.a, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.groupby = function(input, group_by = 1:3, summarise = 4, operation ='count',
                         command = paste0(bedtools, ' groupby %s -i '),
                         arguments = '',
                         process.tag = '.groupby',...){
    group_by = paste0('-g ',paste0(group_by,collapse = ','))
	summarise = paste0('-c ', summarise)
	operation = paste0('-o ', operation)
	arguments = paste(arguments, group_by, summarise, operation, collapse = ' ')
	bedtools.shell(input, command, arguments = arguments ,process.tag = process.tag, ...)
}

bedtools.genomecov = function(input, 
                              command = paste0(bedtools, ' genomecov -g ',bedtools.genome.file,' %s -i '),
                              arguments = ' -bg ',
                              process.tag = '.bdg',
                              output = input %>% append.process.tag(process.tag,replace = T),...){
  bedtools.shell(input, command, arguments = arguments,output = output ,process.tag = process.tag,...)
}

igvtools.totdf = function(input, 
                          pipe.in=input%>% determine.input(),
                          command = paste0(igvtools,' toTDF '),
                          arguments = 'hg19',
                          process.tag = '.tdf',
                          output = input %>% append.process.tag(process.tag,postfix = '.bdg',replace = T), 
#                           tmp.file.prefix = output,
                          execute = T){
  totdf.shell = command
  
    if(pipe.in){
      tmp.file = paste0(output,'.tmp.bdg')
      totdf.shell = paste(input, ' cat - > ', tmp.file,' && ', 
                             totdf.shell)
      input = tmp.file
    }
  totdf.shell = paste0(totdf.shell, input,' ', output ,' ', arguments )
  if(execute){system(totdf.shell,ignore.stderr = T)}
  if(pipe.in & execute){system(paste0('rm ',tmp.file))}
  if(execute){output}else{totdf.shell}
}

bed.2tdf = function(input, output,method = 'totdf', ...){
  output = input %>% bedtools.genomecov(output = 'pipe') %>% igvtools.totdf(output = output,...)
}

grange.2awkfilter = function(grange = NULL,chr = NULL, start = 0 , end = 10000000000,filter = ''){
  if(is.null(grange)){
    if(!is.null(chr)){
      chr = as.character(chr)
      start = as.character(start)
      end = as.character(end)
      filter = paste0('$1=="', chr, '" && $2<', end,' && $3>', start, filter)
    }
  }else{
    chr = seqnames(grange)
    start = start(grange)
    end = end(grange)
    filter = paste0('$1=="', chr, '" && $2<', end,' && $3>', start, filter)
  }
  filter
}

bed.string.2granges = function(bed.string= '', bed.fields = c('chr','start','end','id','score','strand'),col.remove =NULL,split = '\t',...){
  suppressPackageStartupMessages({require(GenomicRanges)})
    bed.fields = bed.fields[1:(str_count(bed.string[1],split)+1)]
  bed_df = colsplit(string = bed.string, pattern = split, names = bed.fields)
  if(is.null(bed_df$strand)){bed_df$strand = '*'}
  bed_df[,col.remove] = NULL
  makeGRangesFromDataFrame(df = bed_df, keep.extra.columns = T, 
                           seqnames.field = 'chr',start.field = 'start',end.field = 'end',strand.field = 'strand',...)%>%
    sort
}

bed2granges = function(bed, ...){
    bed %>% paste0(collapse = ' ')%>% shell.io('cat',pipe.in = F,output='R')%>%
    bed.string.2granges(...)
}

save.bed = function(bed, file){
	write.table(bed, file = file, quote = F,row.names = F,col.names = F,sep = '\t')
	file
}

granges2bed = function(granges, file=NULL, id = NULL, keep.strand = F ,keep.mcols = F){
    bed = granges %>%as.data.frame(stringsAsFactors = F)
    bed$width = NULL
    if(!keep.strand){bed$strand = NULL}
    if(!keep.mcols){bed = bed %>% dplyr::select(seqnames,start,end)}
    if(!is.null(id)){
        if(length(id) == nrow(bed)){
            bed = bed %>% within({id = id})
        }else{
            bed = bed %>% within({id = paste0(seqnames,":",start,'-',end)})
        }
    }
    if(is.null(file)){
        return(bed)
    }else{
        file = save.bed(bed,file)
        return(file)
    }
}
id2granges = function(id,pattern = ':|-', names = c('chr','start','end'),...){
    df = data.frame(id %>%colsplit(':|-',c('chr','start','end')),id = id,...,stringsAsFactors = F)
    df %>% makeGRangesFromDataFrame(keep.extra.columns = T)
}

filter_bdg = function(bdg,cutoff = 1,output = 'R',execute = T){
    awk.filter(bdg,filter = paste0('$4>=',cutoff),output = 'pipe')%>%
    bedtools.merge(output = output, arguments= '-d 1',execute = execute)
}
cat_bdg = function(bed, output = 'R',execute = T){
    paste0(bed,collapse = ' ')%>% 
    shell.io(command = 'cat ',output = 'pipe')%>%
    bedtools.intersect(input.b = bedtools.genome.bed,output = 'pipe',arguments = '-wa')%>%
    bedtools.sort(output = 'pipe')%>% 
    bedtools.genomecov(output = output,execute = execute)
}
consensus_regions = function(bed, output = 'R',n_intersection = 1,remove_bdg=T,execute = T){
    bdg = cat_bdg(bed,
                  output = paste0(
                      bed[1] %>% matched_sub('.*/'),
                      stringi::stri_rand_strings(1,8,pattern = '[a-z]'),
                      '.bdg'
                  ),
                  execute = execute
                 )
    
    merged = filter_bdg(bdg, cutoff = n_intersection, output = output,execute = execute)
    if(remove_bdg&execute){file.remove(bdg)}
    return(merged)
}