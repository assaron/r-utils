loadOrInstall <- function(package) {
    if (!do.call(require, list(package))) {        
        biocLite(package, suppressUpdates=TRUE)    
        do.call(library, package)
    }
}