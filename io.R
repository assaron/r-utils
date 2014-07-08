options(stringsAsFactors=F)

read.tsv <- function(file, header=T, sep="\t", quote="", ...) {
    read.table(file, header=header, sep=sep, quote=quote, ...)     
}

write.tsv <- function(table, dir, row.names=NA, col.names=NA, ...) {
    name <- deparse(substitute(table))
    table <- as.data.frame(table) 

    if (is.na(row.names)) {
        row.names <- is.character(attr(table, "row.names"))
    }

    if (!row.names && is.na(col.names)) {
        col.names=T
    }
    write.table(table, file.path(dir, paste0(name, ".tsv")), quote=F,
                row.names=row.names, col.names=col.names, sep="\t")
}
