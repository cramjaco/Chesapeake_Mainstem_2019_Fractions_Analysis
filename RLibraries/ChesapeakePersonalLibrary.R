
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

ches_plot_options <- list(
  scale_y_log10nice() ,
    scale_x_log10(breaks = my_sizes, labels = as.character(my_sizes)) ,
  geom_point(size = 2) ,
  geom_path(aes(color = as.factor(Station))) ,
  scale_shape_manual(values = rep(21:25, 2)) ,
  scale_fill_viridis_d(option = "plasma") ,
  scale_color_viridis_d(option = "plasma")
)

reformat_sci <- function(x){
  x %>%
    str_replace("e\\+0", " x 10^") %>%
           str_replace("e\\-0", " x 10^-") %>%
           paste0("^")
}
