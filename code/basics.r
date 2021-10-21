z.score = function(x){
  (x - mean(x,na.rm = T))/sd(x,na.rm = T)
}

ma.log= function(x, base = 2,offset = 1){
  log(x = x,base = base)
}

ma.exp = function(x, base = 2, offset =1){
  base^x - offset 
}
