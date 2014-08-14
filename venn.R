getIntersectionSets <-function(name, sets, include=unique(unlist(sets)), exclude=NULL) {
    if (length(sets) == 0) {        
        res <- list(name=setdiff(include, exclude))
        names(res) <- name
        return(res)
    }
    h.name <- names(sets)[1]
    h <- sets[[1]]
    return(c(
        getIntersectionSets(
            paste0(name, ".", h.name, "+"),
            tail(sets, n=length(sets) - 1),
            include=intersect(include, h),
            exclude=exclude
        ),
        getIntersectionSets(
            paste0(name, ".", h.name, "-"),
            tail(sets, n=length(sets) - 1),
            include=include,
            exclude=unique(c(exclude, h))
        )
    ))
}