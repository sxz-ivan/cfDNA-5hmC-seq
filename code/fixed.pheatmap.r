suppressPackageStartupMessages({
    require(tidyverse)
    require(pheatmap)
    require(dendsort)
    require(grid)
    require(gtable)
})

# pre_transformed.abs_cor = pre_transformed %>% cor(method = 'spearman')%>% abs

# pre_transformed.abs_cor.plot_order = pre_transformed.abs_cor %>% dist %>% hclust %>% .$order
hc_ordered = function(hc){hc %>% with(labels[order])}

outlier.trimmer4heatmap = function (x, outlier.threshold = 1.5) 
{
    x = as.matrix(x)
    rowMeidan = rowMedians(x)
    offset = outlier.threshold * rowIQRs(x)
    top = matrix(rowMeidan+offset,nrow = nrow(x),ncol = ncol(x))
    bottom = matrix(rowMeidan-offset,nrow = nrow(x),ncol = ncol(x))
    x = ifelse(x>top,top,x)
    x = ifelse(x<bottom,bottom,x)
    x
}

var_order = function(mat,sort=T){
    hc = hclust(1-cor(mat%>%t)%>%as.dist,"ward.D2")
    if(sort){hc = hc %>%dendsort}
    hc %>%
    hc_ordered
}

sample_order = function(mat,sort=T){
    hc = hclust((mat%>%t)%>%dist,"ward.D2")
    if(sort){hc = hc %>%dendsort}
    hc %>%
    hc_ordered
}

mat_order = function(mat,sort=T){
    mat[var_order(mat,sort),sample_order(mat,sort)]
}

matrix_layout = function(matrix,row.layout = rownames(matrix),col.layout = colnames(matrix),...){
    matrix = matrix[row.layout,col.layout]
    row.layout = rownames(matrix)
    col.layout = colnames(matrix)
    c(environment()%>%as.list,...)%>%
    structure(class =c('matrix_layout','list'))
}

fixed_heatmap.matrix_layout = function(matrix_layout, new_matrix = NULL, ...){
    require(pheatmap)
    if(is.null(new_matrix)){new_matrix = matrix_layout$matrix}
    new_matrix[matrix_layout$row.layout, matrix_layout$col.layout] %>%
    pheatmap(cluster_rows = F,cluster_cols = F,...)
}

fixed_heatmap = function(matrix_layout, ...){
    UseMethod("fixed_heatmap", matrix_layout)
}
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {

  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis

  heatmap <- pheatmap$gtable

  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 

  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")

  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant

    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }

      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }

    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))

    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)

  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.16, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions

  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,grobs = new.flag,t = 4, l = 4)
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  pheatmap$gtable = heatmap
  # return a copy of the heatmap invisibly
  invisible(pheatmap)
}