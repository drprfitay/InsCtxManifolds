

list.dirs.o <- list.dirs
list.dirs.o.2 <- list.dirs
list.dirs <- function(p,full.names=T,recursive=T) {np <- gsub("\\\\", "/", p); list.dirs.o(np, full.names,recursive)}

load.o <- load 
load.o.2 <- load 
load <- function (file, envir = parent.frame(), verbose = FALSE) {np <- gsub("\\\\", "/", file); load.o(np, envir, verbose)}


list.files.o  <- list.files
list.files.o.2 <- list.files
list.files <- function (path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE,
                        recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE,
                        no.. = FALSE) {
  np <- gsub("\\\\", "/", path);
  
  list.files.o(np, pattern, all.files, full.names, recursive, ignore.case, include.dirs, no..)
}


save.o <- save
save.o.2 <- save

save <- function (..., list = character(), file = stop("'file' must be specified"), 
          ascii = FALSE, version = NULL, envir = parent.frame(), compress = isTRUE(!ascii), 
          compression_level, eval.promises = TRUE, precheck = TRUE)  {
  np <- gsub("\\\\", "/", file);
  
  save.o(..., list=list, file=file, ascii=ascii, version=version,
         envir=envir, compress=compress, compression_level=compression_level, 
         eval.promises=eval.promises ,precheck=precheck)
}