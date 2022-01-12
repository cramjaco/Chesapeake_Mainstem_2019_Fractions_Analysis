scientific_10_c <- function(x) {
    xout <- gsub("1e", "10^{", format(x),fixed=TRUE)
    xout <- gsub("{-0", "{-", xout,fixed=TRUE)
    xout <- gsub("{+", "{", xout,fixed=TRUE)
    xout <- gsub("{0", "{", xout,fixed=TRUE)
    xout <- paste(xout,"}",sep="")
    return(parse(text=xout))
    
}

scale_y_log10nice <- function(name=NULL,omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_y_log10(breaks=breaks10,labels=scientific_10_c(breaks10),...)
}