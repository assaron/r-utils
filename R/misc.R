loadOrInstall <- function(package) {
    if (!do.call(require, list(package))) {        
        biocLite(package, suppressUpdates=TRUE)    
        stopifnot(do.call(require, list(package)))
    }
}

setSimilarity <- function(set1, set2){
    c1 <- unique(set1)
    c2 <- unique(set2)            
    
    a <- length(c1)
    b <- length(c2)
    ab <- length(intersect(c1, c2))                    
    if (a * b == 0) {
        0
    } else {        
        ab / sqrt(a * b)            
    }
}

pairwiseCompare <- function(FUN, list1, list2=list1, ...) {
    additionalArguments <- list(...)
    f1 <- function(...) {
        mapply(FUN=function(x, y) {    
            do.call(FUN, 
                    c(list(list1[[x]], list2[[y]]), additionalArguments))
        }, ...)    
    }
    z <- outer(seq_along(list1), seq_along(list2), FUN=f1)    
    rownames(z) <- names(list1) 
    colnames(z) <- names(list2)    
    z
}

messagef <- function(...) {
    message(sprintf(...))
}

warningf <- function(...) {
    warning(sprintf(...))
}

matrixAsColumnList <- function(m) { 
    split(m, col(m)) 
}

"%o%" <- compose
