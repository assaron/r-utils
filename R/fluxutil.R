#' @export
oneReactionDeletion <- function(model, reaction, lb=0, ub=0) {
  model.sd <- model
  model.sd <- changeBounds(model.sd, reaction, lb=lb, ub=ub)
  model.sd.sol <- optimizeProb(model.sd)  
  model.sd.sol
}

#' @export
allOneReactionDeletions <- function(model, reaction.to.compare) {
  model.sol <- optimizeProb(model)
  r.t.c.pos <- checkReactId(model, reaction.to.compare)@react_pos
  opt.flux <- model.sol@fluxdist@fluxes[r.t.c.pos]
  if (min(opt.flux) < 1e-10) {
    stop(paste0("Optimal flux via ", reaction.to.compare, " is zero"))
  }
  sapply(model@react_id, USE.NAMES=T, function(r.id) {    
    model.sd.sol <- oneReactionDeletion(model, r.id)    
    res <- model.sd.sol@fluxdist@fluxes[r.t.c.pos] / opt.flux
    if (max(res) < 0.9) {
      message(paste0(r.id, ": ", res))  
    }
    res  
  })
}

#' @export
fluxVia <- function(model, solution, reaction) {
  solution@fluxdist[checkReactId(model, reaction)@react_pos]  
}

#' @export
compareFluxes <- function(models, filter.zero=T, react=NULL) {
  m.inf <- max(sapply(models, uppbnd))
  
  models.f <- sapply(models, function(model) {
    if (!is.null(react)) {
      model <- changeObjFunc(model, react)      
    }
    
    model.sol <- optimizeProb(model, algorithm="mtf")
    #     message(show(model.sol))
    model.f <- getFluxDist(model.sol)
    names(model.f) <- model@react_id  
    model.f
  }, USE.NAMES=T)
  
  if (filter.zero) {
    models.f <- models.f[apply(abs(models.f), 1, max) > 1e-10, ]  
  }
  
  k <- length(models)
  ratios <- list()
  i <- 1
  for (j in (seq_len(k - i) + i)) {
    r <- list(models.f[, j] / models.f[, i])
    names(r) <- paste0(names(models)[j], "/", names(models[i]))
    ratios <- c(ratios, r)      
  }    
  
  ratios <- do.call(cbind, ratios)
  
  
  diffs <- list()
  i <- 1
  for (j in (seq_len(k - i) + i)) {
    d <- list(models.f[, j] - models.f[, i])
    names(d) <- paste0(names(models)[j], "-", names(models[i]))
    diffs <- c(diffs, d)      
  }    
  
  diffs <- do.call(cbind, diffs)
  
  models.f <- as.data.frame(models.f)
  models.f$zeroFlux <- apply(abs(models.f), 1, max) < 1e-10
  models.f <- cbind(models.f, ratios, diffs)  
  models.f
}

#' @export
getFluxGraph <- function(m, m.f=NULL) {
  if (is.null(m.f)) {
    m.sol <- optimizeProb(m, algorithm="mtf")
    m.f <- getFluxDist(m.sol)    
  }
  names(m.f) <- m@react_id
  
  getTables <- function(rxn, met) {
    rxn.id <- m@react_id[rxn]
    met.id <- m@met_id[met]
    
    if (met.id %in% mets2mask) {
      met.id <- paste0(met.id, "_", rxn.id)
    }
    rxn2met.flux <- m.f[rxn] * m@S[met, rxn]
    
    if (rxn2met.flux > 0) {
      from <- rxn.id
      to <- met.id
    } else {
      from <- met.id
      to <- rxn.id
    }
    
    list(
      node.table=rbind(
        data.table(id=met.id, label=met.id, nodeType="met", masked=m@met_id[met] %in% mets2mask),
        data.table(id=m@react_id[rxn], label=m@react_id[rxn], nodeType="rxn", masked=FALSE)
      ),
      edge.table=data.table(from=from, to=to, flux=abs(rxn2met.flux)))    
  }
  
  
  nz.rxns <- which(abs(m.f) > 1e-10)  
  rxn <- nz.rxns[1]
  
  tables <- lapply(nz.rxns, function(rxn) {
    tables1 <- lapply(which(m@S[, rxn] != 0), function(met) { getTables(rxn, met) })
    list(
      node.table=rbindlist(lapply(tables1, "[[", "node.table")),
      edge.table=rbindlist(lapply(tables1, "[[", "edge.table"))
    )    
  })
  
  node.table <- unique(rbindlist(lapply(tables, "[[", "node.table")))
  edge.table <- unique(rbindlist(lapply(tables, "[[", "edge.table")))
  
  res <- GAM:::graph.from.tables(
    node.table=as.data.frame(node.table),
    edge.table=as.data.frame(edge.table))
  
  attr(res, "model") <- m
  attr(res, "fluxDist") <- m.f
  
  res 
}

#' @export
changeObjFuncRel <- function(model, reacts) {
  reacts.max <- sapply(reacts, function(react) {
      optimizeProb(changeObjFunc(model, react))@lp_obj
    })
  
  changeObjFunc(model, reacts, 1 / reacts.max)
}
