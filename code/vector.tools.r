suppressPackageStartupMessages({
  require(tidyverse)
})

matched_sub <- function(x,r){regmatches(x, regexpr(r, x, perl=TRUE))}
extract.between = function(x, leading = '^',trailing = '$'){
  sub(paste0(pattern = '.*',leading,' *(.*?) *',trailing,'.*'),replacement = "\\1",x)
}

strsplit.extract = function(x,position = 1, split = '|',fixed = T,...){
  sapply(strsplit(x = as.character(x),split = split,fixed = fixed,...),function(x){x[position]})
  }

extract.duplicates = function(x){
  x %in% x[duplicated(x)]
}

is_fresh = function(x){
  !duplicated(x)
}

table.na = function(x){is.na(x)%>%table}

is.empty = function(x){is.na(x)|nchar(x) == 0}

subset.by.quantile = function(x,range = c(0,100)){
	if(range[1]<range[2]){
		subsetter = (x>quantile(x,range[1])&x<quantile(x,range[2]))
	}else{
		subsetter = (x>quantile(x,range[1])|x<quantile(x,range[2]))
	}
	subsetter
}

index.of.logic = function (x, logic = T) 
{
    (1:length(x))[x == logic]
}

my.var.summary = function(x){
  x.table = table(x,useNA = 'ifany')%>% as.data.frame
  x.table$props =( x.table$Freq/length(x)*100) %>% round()%>%paste0('%')
  x.table$report = paste0(x.table$Freq,'（',x.table$props,'）')
  x.table
}

no.zero = function(x){x+min(x[x>0])}

as.list.infer.named = function(object, infer.FUN,...){
  object.list = object %>% as.list
  names(object.list) = infer.FUN(object, ...)
  object.list
}

name.list.lapply = function(name.list, ...){
  result = name.list %>%as.list %>% lapply(...)
  names(result) = name.list
}

named.object = function(x, names=stri_rand_strings(length(x),6)){
  names(x) = names
  x
}

find.index = function(subject, query){
    if(is.null(query)){return(int())}
    (1:length(subject))[subject %>% sapply(function(x){all(x == query)})]
}

pop = function(subject , pop=NULL, index = find.index(subject, pop)){
    if(length(index)==0){return(subject)}
    subject[-index]
}

find.nth = function(x, nth=1, decreasing = T){
	x[order(x,decreasing= decreasing)][nth]}

xgsub = function(x, pattern, replacement,ignore.case = FALSE, perl = FALSE, 
                fixed = FALSE, useBytes = FALSE){
  gsub(pattern, replacement,x,ignore.case,perl,fixed, useBytes)
}

xsub = function(x, pattern, replacement,ignore.case = FALSE, perl = FALSE, 
                fixed = FALSE, useBytes = FALSE){
  sub(pattern, replacement,x,ignore.case,perl,fixed, useBytes)
}

seq_frac = function(from, to = 1, frac = 0.5){
  if(frac>=1|to>from){
    return(c(to,from))
  }
  out = from
  x = floor(from*frac)
  while(x > to){
    out = c(out,x)
    x = floor(x * frac)
  }
  out =sort(c(out,to))
  return(out)
}

cut_levels = function(cuts,include.lowest= F){
    cuts = sort(cuts)
    paste0(c(ifelse(include.lowest,'(','['),rep('(',length(cuts)-1)),cuts[1:length(cuts)-1],',',cuts[-1],']')
}
ordered_num_factor = function(x,...){
    factor(as.character(x),as.character(sort(unique(x),...)))
}