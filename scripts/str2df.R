#get structure data for breeding and migrant runs
str2df <- function(str_in,melt=T){
  str <- readLines(str_in,warn=F)
  start <- grep("Inferred ancestry of individuals:",str)+2
  nind <- str[grep("Run parameters:",str)+1] %>% strsplit(" * ") %>% unlist() %>% .[2] %>% as.numeric()
  k <- str[grep("Run parameters:",str)+3] %>% strsplit(" * ") %>% unlist() %>% .[2] %>% as.numeric()
  str <- str[start:(start+nind-1)]
  str <- lapply(str,function(e){
    strsplit(e," * ") %>% unlist() %>% .[c(3,7:(7+k-1))]
  })
  str <- do.call(rbind.data.frame,str)
  str[,2:(2+k-1)] <- apply(str[,2:(2+k-1)],2,function(e) as.numeric(e))
  colnames(str) <- c("id",sapply(1:k,function(e) paste0("q",e)))
  if(melt){
    return(melt(str,id.vars="id"))
  } else {
    return(str)
  }
}