options(stringsAsFactors=F)

read.tsv <- function(file, header=T, sep="\t", quote="", comment.char="", check.names=FALSE, ...) {
    read.table(file, header=header, sep=sep, quote=quote, 
               comment.char=comment.char, check.names=check.names,
               ...)     
}

write.tsv <- function(table, dir, file=NULL, gzip=FALSE, row.names=NA, col.names=NA, ...) {
    name <- deparse(substitute(table))
    table <- as.data.frame(table) 
    
    if (is.null(file)) {
        file <- file.path(dir, paste0(name, ".tsv", if (gzip) ".gz"))        
    }

    if (is.na(row.names)) {
        row.names <- is.character(attr(table, "row.names"))
    }

    if (!row.names && is.na(col.names)) {
        col.names=T
    }
    
    for (c in colnames(table)) {
        if (is.character(table[[c]])) {
            table[[c]] <- sub("#", "", table[[c]])            
        }
    }
    
    if (gzip) {
        file <- gzfile(file, "w")
    }
    write.table(table, file, quote=F,
                row.names=row.names, col.names=col.names, sep="\t")
    if (gzip) {
        close(file)
    }
}
