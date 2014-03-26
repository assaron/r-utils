write.tsv <- function(table, dir, ...) {
  name <- deparse(substitute(table))
  write.table(table, file.path(dir, paste0(name, ".tsv")), quote=F, row.names=T, col.names=NA, sep="\t")
}
