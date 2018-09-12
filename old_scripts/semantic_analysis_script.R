toponyms = c("virginia","naples","bohemia","england","france","rome","peru",
             "mexico","jamaica","boston","paris","dublin","edinburgh","egypt",
             "greece","italy","china","india","mexico","oxford","cambridge")

setwd("~/Desktop/geosemantics/inst/extdata")

results = list()
for (i in 1:5) {
  load(paste("~/Desktop/geosemantics/inst/extdata/ch1_kd",i,".rda", sep=""))
  load(paste("~/Desktop/geosemantics/inst/extdata/ch1_pd",i,".rda", sep=""))
  
  kd = ppmi(kd)
  avgs = apply(kd, 1, mean)
  devs = apply(kd, 1, sd)
  lsa = irlba::irlba(kd, nv = dims)
  kds = lsa$u %*% diag(lsa$d)
  rownames(kds) = rownames(kd)
  
  results[[i]] = list()
  
  for (j in 1:length(toponyms)) {
    print(paste(i, j))
    results[[i]][[j]] = list()
    
    place = toponyms[j]
    places = place_join(place, gaz = gaz, mode = "same")
    vec = place_composition(pd, places)
    ord = which(vec > 0)
    if (length(ord) > 50) {
      dims = 50
    } else {
      dims = round(length(ord) / 5)
    }
    if (dims < 5) next()
    pavgs = apply(kd[,ord], 1, mean, na.rm = T)
    z = (pavgs - avgs) / devs
    
    lsa = irlba::irlba(kd[,ord], nv = dims)
    w = lsa$u %*% diag(lsa$d)
    rownames(w) = rownames(kd)
    
    m1 = w %*% t(w)
    
    
    m2 = kds %*% t(kds)
    sims = c()
    for (s in 1:nrow(kds)) {
      sim = cosine_similarity(m1[s,], m2[s,])
      sims = c(sims, sim)
    }
    names(sims) = rownames(w)
    
    x = z / abs(z)
    
    results[[i]][[j]]$toponym = place
    results[[i]][[j]]$doc_frequency = length(ord) / nrow(pd)
    results[[i]][[j]]$summary = summary(sims)
    results[[i]][[j]]$top_words = sort(z * sims, decreasing = T)[1:10]
    results[[i]][[j]]$distinctive_words = sort(pavgs * x * (1 - sims), decreasing = T)[1:10]
    
    semantic_results = list()
    for (k in 1:10) {
      term = names(results[[i]][[j]]$distinctive_words)[k]
      semantic_results[[k]] = list()
      semantic_results[[k]]$local = semantic_similarity(w, term)
      semantic_results[[k]]$global = semantic_similarity(kds, term)
    }
    results[[i]][[j]]$semantic_results = semantic_results
  }
}
save(results, file = "results.rda")

# kd = ppmi(kd)
# avgs = apply(kd, 1, mean)
# devs = apply(kd, 1, sd)
# lsa = irlba::irlba(kd, nv = 50)
# kds = lsa$u %*% diag(lsa$d)
# rownames(kds) = rownames(kd)
# 
# places = place_join(c("jamaica"), gaz = gaz, mode = "same")
# vec = place_composition(pd, places)
# ord = which(vec > 0)
# pavgs = apply(kd[,ord], 1, mean, na.rm = T)
# z = (pavgs - avgs) / devs
# 
# lsa = irlba::irlba(kd[,ord], nv = 50)
# w = lsa$u %*% diag(lsa$d)
# rownames(w) = rownames(kd)
# 
# m1 = w %*% t(w)
# 
# 
# m2 = kds %*% t(kds)
# sims = c()
# for (i in 1:nrow(kds)) {
#   sim = cosine_similarity(m1[i,], m2[i,])
#   sims = c(sims, sim)
# }
# names(sims) = rownames(w)
# 
# # The lexical profile of a place is top 10 words over each metric
# sort(z * sims, decreasing = T)[1:10]
# x = z / abs(z)
# sort(pavgs * x * (1 - sims), decreasing = T)[1:10]
# 
# # Then to compare uses of the terms
# semantic_similarity(w, "ocean")
# semantic_similarity(kds, "ocean")
# 
# summary(sims)
