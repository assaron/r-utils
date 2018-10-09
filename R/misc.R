#' @export
loadOrInstall <- function(package) {
    if (!do.call(require, list(package))) {        
        biocLite(package, suppressUpdates=TRUE)    
        stopifnot(do.call(require, list(package)))
    }
}

#' @export
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

#' @export
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

#' @export
messagef <- function(...) {
    message(sprintf(...))
}

#' @export
warningf <- function(...) {
    warning(sprintf(...))
}

#' @export
matrixAsColumnList <- function(m) { 
    split(m, col(m)) 
}

#' @export
"%o%" <- pryr::compose

#' @export
ulength <- pryr::compose(length, unique)

#' @export
replaceNA <- function(x, y) { x[is.na(x)] <- y; x}

#' @export
'%f%' <- function(s, a) do.call("sprintf", c(s, as.list(a))) 

#' @export
minPos <- function(x) min(x[x > 0])

#' Remove layer from ggplot object
#' from https://stackoverflow.com/questions/13407236/remove-a-layer-from-a-ggplot2-chart
#' @export
removeLayer <- function(ggplot2_object, geom=NULL, color=NULL) {
    # Delete layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
        if (!is.null(geom) && (class(x$geom)[1] != geom)) {
            return(x)
        } 
        
        if (!is.null(color) && (!("colour" %in% names(x$aes_params)) ||
                                (x$aes_params$colour != color))) {
            return(x)
        } 
        NULL
    })
    # Delete the unwanted layers.
    layers <- layers[!sapply(layers, is.null)]
    ggplot2_object$layers <- layers
    ggplot2_object
}
