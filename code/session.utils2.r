suppressPackageStartupMessages({
    require(tidyverse)
#     require(fastSave)
#     require(stringi)
    require(simpleCache)
    require(feather)
})

source('code/vector.tools.r')
source('code/parallel.tools.r')

generate_filename = function(var_name, 
                             data_home = .GlobalEnv$configs$data_home, 
                             parent=.GlobalEnv$configs$current_node
                            ){
    file_head = paste0(data_home, '/',parent, '-',var_name)
    file_tail =case_when(
        is.data.frame(.GlobalEnv[[var_name]])~'.feather',
        is.matrix(.GlobalEnv[[var_name]])~'.feather',
        TRUE~'.RData'
    )
    paste0(file_head,file_tail)
}

vsave = function(var_name,filename = generate_filename(var_name),show = .GlobalEnv$configs$vsave_show, just_update_time = F){
# 一旦保存即更新
    #updating session meta
    .GlobalEnv$session[[var_name]] = list()%>%within({
#         id = stri_rand_strings(1,length = 10,pattern = '[A-Za-z]')
        parent = .GlobalEnv$configs$current_node
        var_name = var_name
        if(is_saved(var_name) & just_update_time){
          file = filename
            system(paste0('touch ', file))
        }else{
          file = fsave(.GlobalEnv[[var_name]], filename)
        }
        time = file.info(file)[['mtime']]
    })
    .GlobalEnv$session%>% fsave(.GlobalEnv$configs$current_session)
    #updating vchildren
    .GlobalEnv$vchildren = c(.GlobalEnv$vchildren,var_name)
    
    if(show){
        return(.GlobalEnv$session[[var_name]])
    }else{
        invisible(.GlobalEnv$session[[var_name]])
    }
}

vssave = function(var_names){
    for(var_name in var_names){vsave(var_name)}
}

fsave = function(...){
    UseMethod('fsave')
}

fsave.default = function(var,filename){
    saveRDS(var, filename)
    filename
}

fsave.data.frame = function(var,filename){
  var$ROWNAMES = rownames(var)
    write_feather(var, filename)
    filename
}

fsave.matrix = function(var,filename){
    var = var %>% as.data.frame(stringsAsFactors = F)
    var$ROWNAMES = rownames(var)
    write_feather(var,filename)
    filename
}

meta_from_file = function(file){
    within(list(),{
        parent = file%>% extract.between('/','-')
        var_name = var_name %>% extract.between('-','.')
        file = file
        time = file.info(file)[['mtime']]
    })
}

meta_from_name =function(var_name){.GlobalEnv$session[[var_name]]}

vload = function(meta){
    .GlobalEnv[[meta$var_name]] = fload(meta$file)
    .GlobalEnv$session[[meta$var_name]] = meta
}
#所有用到过的数据都会被session记录下来，所以不用转么再存了，只用存session就行了

vsload = function(metas,reload = F,mc.cores = 1){
  metas = metas[metas %>% lapply(function(meta){if_load(meta$var_name,reload)}) %>% unlist]
  # print(metas)
  metas %>% flexible.mclapply(mc.cores = mc.cores,vload)
  # vsload的中的ifload一方面会对比文件是否已经加在，另一方面会对比已经加在的文件和储存文件的时间
	metas
}

fload = function(file,default = list()){
	if(file.exists(file)){
		data = tryCatch({
            read_feather(path = file)%>% as.data.frame(stringsAsFactors = F)%>%
            (function(data){
                if('ROWNAMES'%in% colnames(data)){
                  rownames(data) = data$ROWNAMES
                  data$ROWNAMES = NULL
                }
                data
            })
        },error=function(e){readRDS(file)})
	}else{data = default }
	data
}

is_old = function(var_name){
  file = .GlobalEnv$session[[var_name]]$file
  var_time = .GlobalEnv$session[[var_name]]$time 
  if(file.exists(file)){return(var_time< file.info(file)[['mtime']])}else{return(F)}
  #file至少都比meta新，因为meta一旦更新必然发生file更新（vsave），而file更新后meta不一定更新
}

if_load = function(var_name, reload=F){
    is_loaded = var_name %in% names(.GlobalEnv$session)
    is_missing = !is_loaded
    if(is_loaded){
        is_old = var_name %>% sapply(is_old)
    }else{is_old = rep(F,length(var_name))}
    is_missing | is_old | reload
}

update_session = function(session = NULL, reload = F){
    if(is.null(session)){
        .GlobalEnv$session = .GlobalEnv$session %>% vsload(reload)
    }else{invisible(session)}
}

refresh_session = function(keep_file = F){
  current_vars = names(.GlobalEnv$session) %>% setdiff(c('session','session_generations','required_pkgs'))
  for(current_var in current_vars){
    vdel(current_var,keep_file = .GlobalEnv$session[[current_var]]$parent != .GlobalEnv$configs$current_node | keep_file)
    }
  c('session_generations','required_pkgs') %>% vssave
  session %>% fsave(configs$current_session)
}

query_time = function(var_name, default = Sys.time()){
    if(var_name %in% names(.GlobalEnv$session)){
        result = .GlobalEnv$session[[var_name]]$time
    }else{
        result = default
    }
    return(result)
}

if_run = function(input=NULL, output=NULL, anyway=NULL){
    #这里就不管数据是否更新的问题了，要更新用别的命令更新
    input_is_present = input %in% names(.GlobalEnv$session)
    if(!all(input_is_present)){
        print('missing input')
        print(input[!input_is_present])
        return(F)
    }
    if(is.null(anyway)){
        output_is_present = output %in% names(.GlobalEnv$session)
        output_is_missing = !output_is_present 
        if(output_is_present){
            input_is_newer = max(input %>% lapply(query_time) %>% unlist) > query_time(output)
            #input里面最新的那个跟output比
        }else{
            input_is_newer = T
        }
        return(output_is_missing | input_is_newer)
    }else{
        return(anyway)
    }
}

is_saved = function(var_name){
    is_present = var_name %in% names(.GlobalEnv$session)
    if(is_present){is_saved = file.exists(.GlobalEnv$session[[var_name]]$file)}else{is_saved = F}
    is_saved
}

vdel = function(var_name,keep_file = T){
    if(!keep_file){file.remove(.GlobalEnv$session[[var_name]]$file)}
    .GlobalEnv$session[[var_name]] = NULL
    .GlobalEnv$session%>% fsave(.GlobalEnv$configs$current_session)
    .GlobalEnv[[var_name]] = NULL
}

load_required = function(session, requirements = .GlobalEnv$configs$require_from_parent){
  if(is.null(requirements)){
    return(session)
  }else{
    return(session[intersect(requirements,names(session))])
  }
}

split_fsave = function(obj,base,specify=NULL){
    objnames = obj %>% names
    not_specified = setdiff(objnames, specify %>% names)
    rest = not_specified %>% as.list %>% setNames(not_specified) %>% 
    lapply(function(item){
        itemclass = class(obj[[item]])
        filename = paste0(base,'_',item,'.',ifelse(itemclass %in% c('data.frame','matrix'),'feather','RData'))
        fsave(obj[[item]],filename)
        obj[[item]]=NULL
        gc()
        list(class=itemclass,filename=filename)
    })
    for(n in names(obj)){obj[[n]]=NULL}
    trunkname = paste0(base,'_TRUNK.RData')
    fsave(obj,trunkname)
    c(list(TRUNK=trunkname),specify,rest)
}

split_fload = function(recipe){
    obj = recipe[['TRUNK']] %>% fload
    for(item in names(recipe)[names(recipe)!='TRUNK']){
        obj[[item]] = as(fload(recipe[[item]]$filename),recipe[[item]]$class)
    }
    obj
}

sch = function (cacheName, instruction = NULL, 
                loadFromFile = NULL,cacheSymlink = TRUE,
                buildEnvir = NULL, reload = FALSE, 
    recreate = FALSE, noload = FALSE, cacheDir = getCacheDir(), 
    cacheSubDir = NULL, 
    assignToVariable = NULL, loadEnvir = parent.frame(), searchEnvir = getOption("SIMPLECACHE.ENV"), 
    nofail = FALSE, lifespan = NULL) 
{
    if (!"character" %in% class(cacheName)) {
        stop("simpleCache expects the cacheName variable to be a character vector.")
    }
    instruction = substitute(instruction)
    if ("character" %in% class(instruction)) {
        message("Character instruction; consider wrapping in braces.")
        parse = TRUE
    }
    else {
        parse = FALSE
    }
    if (!is.null(cacheSubDir)) {
        cacheDir = file.path(cacheDir, cacheSubDir)
    }
    if (is.null(cacheDir)) {
        message(strwrap("No cacheDir specified. You should set global option\n\t\tRCACHE.DIR with setCacheDir(), or specify a cacheDir parameter directly\n\t\tto simpleCache(). With no other option, simpleCache will use tempdir():\n\t\t", 
            initial = "", prefix = " "), tempdir())
        cacheDir = tempdir()
    }
    if (!file.exists(cacheDir)) {
        dir.create(cacheDir, recursive = TRUE)
    }
    cacheFile = file.path(cacheDir, paste0(cacheName, ".RData"))
    submitted = FALSE
    searchEnvir = append(searchEnvir, ".GlobalEnv")
    cacheExists = FALSE
    cacheWhere = NULL
    for (curEnv in searchEnvir) {
        if (!(exists(curEnv) && is.environment(get(curEnv)))) {
            warning(curEnv, " is not an environment.")
        }
        else if (exists(cacheName, where = get(curEnv))) {
            cacheExists = TRUE
            cacheWhere = curEnv
            break
        }
    }
    ret = NULL
#     if (.tooOld(cacheFile, lifespan)) {
#         message(sprintf("Stale cache: '%s' (age > %d day(s))", 
#             cacheFile, lifespan))
#         recreate = TRUE
#     }
#     print(class(instruction))
    if (cacheExists & !reload & !recreate & file.exists(cacheFile)) {
        message("::Object exists (in ", cacheWhere, ")::\t", 
            cacheName)
        ret = get(cacheName, pos = get(cacheWhere))
    }
    else if (file.exists(cacheFile) & !recreate & !noload) {
        message("::Loading cache::\t", cacheFile)
        ret = fload(cacheFile)
    }
    else if (file.exists(cacheFile) & !recreate) {
        message("::Cache exists (no load)::\t", cacheFile)
        return(NULL)
    }
    else {
        message("::Creating cache::\t", cacheFile)
        tryCatch({
            if (is.null(instruction)) {
                if(!is.null(loadFromFile)){
					ret = fload(loadFromFile)
                }  
			}
            else {
                if (is.null(buildEnvir)) {
					if (parse) {
					  ret = eval(parse(text = instruction), envir = parent.frame())
					}
					else {
					  ret = eval(instruction, envir = parent.frame())
					}
                }
                else {
                  if (exists("instruction", buildEnvir)) {
                    stop("Can't provide a variable named 'instruction' in buildEnvir")
                  }
                  buildEnvir$instruction = instruction
                  be = as.environment(buildEnvir)
                  parent.env(be) = parent.frame()
                  
                  if (parse) {
                    ret = with(be, eval(parse(text = instruction)))
                  }
                  else {
                    ret = with(be, eval(instruction))
                  }
                }
            }
        }, error = function(e) {
            if (nofail) 
                warning(e)
            else stop(e)
        })
        if (is.null(ret)) {
            message("NULL value returned, no cache created")
            return()
        }
        else {
			if(cacheSymlink & !is.null(loadFromFile)){
				system(paste0('ln -sf ',loadFromFile,' ',cacheFile))
			}else{
				fsave(ret, file = cacheFile)
			}
        }
    }
    if (noload) {
        rm(ret)
        gc()
        return()
    }
    if (is.null(assignToVariable)) {
        assignToVariable = cacheName
    }
    assign(assignToVariable, ret, envir = loadEnvir)
}
nodeSubDir = function(nodes,fromHome = T){
    if(length(nodes)==1 & fromHome){
        nodes_i = 1:match(nodes,names(cacheNodes))
    }else{
        nodes_i = match(nodes,names(cacheNodes))
    }
    paste0(unlist(cacheNodes[nodes_i]),collapse = '/')
}
schdir = function (cacheDir = NULL) 
{
    dir.create(cacheDir,showWarnings = F, recursive = TRUE)
    setCacheDir(cacheDir)
}
sch.update = function(var){
    sch(var,recreate = T)
}
lch = function (cacheName, cacheDir = getCacheDir()) 
{
    fload(paste0(cacheDir,'/',cacheName,'.RData'))
}