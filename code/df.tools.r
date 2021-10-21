suppressPackageStartupMessages({
  require(tidyverse)
})

big.head = function(df,show.col = 1:10){head(df[,show.col])}

big.read.table = function(file,nscan = 10, header = F, sep = "", quote = "\"'",stringsAsFactors = default.stringsAsFactors(),...){
  colClasses= sapply(read.table(file,nrows = nscan, 
                              header = header, sep = sep, quote = quote,
                              stringsAsFactors =stringsAsFactors),class)
  colClasses[!(colClasses %in% c('character','factor','numeric','integer'))] = 'character'
  read.table(file,header = header, sep = sep, quote = quote,
             colClasses = colClasses,
             stringsAsFactors =stringsAsFactors)
}

big.sample_n = function(big,ncol = 10,nrow = 10){sample_n(big[,sample(1:ncol(big),size =ncol)],size = nrow)}

arrange_at.keep.rownames = function(df,arrange.col,decreasing = T){
  df[order(df[[arrange.col]],decreasing = decreasing),]
}

make.rownames = function(df,rowname.column = 1,keep.rowname.column = T){
  df = as.data.frame(df)
  rownames(df) = df[[rowname.column ]]
  if(!keep.rowname.column){df[[rowname.column ]] = NULL}
  df
}

keep.rownames = function(df, id.col = 'rownames.2id'){
  df = as.data.frame(df)
  df[[id.col]] = rownames(df)
  df
}

change.colnames = function(df, colnames){
  colnames(df) = colnames
  df
}

rule.wierd.numbers = function(x,to = NA){
  x[!is.finite(x)]=to
  x
}

operate.inpipe = function(data,operations){
  environment = list()
  environment$data = data
  with(environment,operations)
  }

df.dedup = function(df,keys = colnames(df)){df[!(duplicated(df[,keys])),]}

sample_n = function(data,size,...){
  size = min(size,nrow(data))
  data %>% keep.rownames %>% dplyr::sample_n(size = size,...) %>% 
    make.rownames(rowname.column = 'rownames.2id',keep.rowname.column = F)
}

random.split = function(data,n.split,format = 'data'){
  if(is.vector(data)){size = length(data)}else{size = nrow(data)}
  f = sample(1:n.split,size = size,replace = T)
  if(format=='data'){output = data %>% split(f = f)}
  if(format=='number'){output = (1:size) %>%split(f = f)}
  if(format=='row.names' & (!is.null(rownames(data)))){output = rownames(data)%>%split(f = f)}
  output
}

df.split = function(df, f, sep = '-'){
    if(length(f)>1){
    split.f = df[,f,drop = F] %>% Reduce(f=function(x,y){paste0(x,sep,y)})
    }else{split.f = df[[f]]}
	df %>% split(f = split.f)
}

df.rowlist = function(df,names=rownames(df)){
    split(df,1:nrow(df))%>%setNames(names)
}

row_split = function(df, f){
  if(nrow(df)!=length(f)){return(df)}else{1:nrow(df) %>% split(f = f) %>% lapply(function(index){df[index,]})}
}

col_split = function(df, f){
  if(ncol(df)!=length(f)){return(df)}else{1:ncol(df) %>% split(f = f) %>% lapply(function(index){df[,index]})}
}

bind_rows.rownames = function(...){
    row.names = lapply(...,rownames) %>% unlist
    out = bind_rows(...)
    rownames(out) = row.names
    out
}
read.table1 = function(file,trial_nrows = 5,specified_classes = NULL, execute = T,nrows = -1,header = FALSE, sep = "", quote = "\"'", dec = ".", 
                       numerals = c("allow.loss", "warn.loss", "no.loss"), row.names, 
                       col.names, as.is = !stringsAsFactors, na.strings = "NA", 
                       colClasses = NA,  skip = 0, check.names = TRUE, 
                       fill = !blank.lines.skip, strip.white = FALSE, blank.lines.skip = TRUE, 
                       comment.char = "#", allowEscapes = FALSE, flush = FALSE, 
                       stringsAsFactors = default.stringsAsFactors(), fileEncoding = "", 
                       encoding = "unknown", text, skipNul = FALSE){
  colClasses = do.call(read.table, environment() %>% as.list()%>% within({nrows = trial_nrows}) %>% .[-c(2:4)]) %>% sapply(class)
  if(!is.null(specified_classes)){
    colClasses[names(specified_classes)] = specified_classes
  }
  if(!stringsAsFactors){
    colClasses[colClasses=='factor'] = 'character'
  }
  if(execute){
    out = do.call(read.table, environment() %>% as.list()%>% .[-c(2:4)])
  }else{
    out = environment() %>% as.list()%>% .[-c(2:4)]
  }
  return(out)
}
vlookup = function(x,y,key.y=1,value.y=2,default = NA){
  y = y[!duplicated(y[[key.y]]),] %>% make.rownames(rowname.column = key.y)
  p = match(x,y[[key.y]])
  ifelse(is.na(p),default,y[p,value.y])
}
use_df = function(df,...,fill=NA){
    new = data.frame(...,stringsAsFactors = F)
    for(column in setdiff(colnames(df),colnames(new))){
        new[[column]] = fill
    }
    new[,union(colnames(df),colnames(new))]
}
inflate = function(df,cols,fill=NA){
    new_cols = setdiff(cols,colnames(df))
    for(nc in new_cols){
        df[[nc]] = fill
    }
    df[,cols]
}


setRowNames = function(X,row.names){
    rownames(X)=row.names
    X
}
pastedRownames = function(df,collapse='_'){
    df %>% setRowNames(apply(.,1,paste,collapse=collapse))
}

expandGrid = function(...,sep='_', KEEP.OUT.ATTRS = TRUE, stringsAsFactors = F){
    g = expand.grid(..., KEEP.OUT.ATTRS = TRUE, stringsAsFactors = F)
#     rownames(g) = apply(g,1,paste,collapse=sep)
    g
}
extendGrid = function(grid,...){
    args = list(...)
    for(argi in 1:length(args)){
        values = args[[argi]]
        if(class(values)=='data.frame'){
            combined_colnames = c(colnames(grid),colnames(values)) %>% unique
            grid = merge(grid,values,all.x=T,by=intersect(colnames(grid),colnames(values)))[combined_colnames]
        }else{
            argn = names(args)[argi]
            values = unique(values)
            grid_list = list()
            for(v in values){
                grid[[argn]]=v
                grid_list = c(grid_list,list(grid))
            }
            grid = bind_rows(grid_list)
        }
    }
    grid
}
setupDummyCols = function(df,cols,with = NA){
    if(nrow(df)==0){
        df[1,] = with
        empty = T
    }else{empty=F}
    new_cols = setdiff(cols,colnames(df))
    for(col in new_cols){
        df[[col]] = with
    }
    if(empty){
        df = df[-1,]
    }
    return(df[,cols])
}