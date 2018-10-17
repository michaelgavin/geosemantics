gaz = read.csv("~/Downloads/gaz - gaz.csv", stringsAsFactors = F)
pnames = unique(gaz$SUBJECT)
load("~/Downloads/pd1660 (1).rda")
load("~/Downloads/kd1660.rda")
pd = pd[rownames(pd) %in% pnames,]

ppmi = function(mat) {
  total = sum(mat, na.rm = T)
  pcol = apply(mat, 2, function(x) sum(x, na.rm = T) / total)
  prow = apply(mat, 1, function(x) sum(x, na.rm = T) / total)
  for (i in 1:nrow(mat)) {
    row = mat[i,]
    row = row / sum(row)
    ord = which(row < 0)
    row[ord] = 0
    row = log(row / (prow[i] * pcol))
    mat[i,] = row
  }
  mat[is.na(mat)] = 0
  mat[is.infinite(mat)] = 0
  return(mat)
}

relative_frequency = function(words, mat1, mat2) {
  
  freq1 = sum(rowSums(mat1[words,]) / sum(mat1))
  freq2 = sum(rowSums(mat2[words,]) / sum(mat2))
  
  return(freq1 / freq2)
  
}

# Get same as places
place = "carolina"
pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
vec = colSums(py[pnames,])
plot(py[place,], type = "l", lwd = 2)
lines(vec, type = "l", lwd = 2, col = "red")

# Get inside places
x = 1:226
place = "malacca"
pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% pnames))])
pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
pnames = pnames[pnames %in% rownames(py)]
pnames = unique(pnames)
if (length(pnames) > 1) {
  vec = colSums(py[pnames,])
} else {
  vec = py[place,]
}
fit = lm(vec ~ x)
plot(vec, type = "l", lwd = 2, main = paste(place, paste(fit$coefficients, collapse = " "), sep = " "))
abline(fit)

for (j in 1:ncol(py)) {
  vec = py[,j]
  vec = vec / sum(vec)
  py[,j] = vec
}

slopes = matrix(0, nrow(py), 2)
rownames(slopes) = rownames(py)
for (i in 1:nrow(py)) {
  print(i)
  place = rownames(py)[i]
  pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = pnames[pnames %in% rownames(py)]
  pnames = unique(pnames)
  if (length(pnames) > 1) {
    vec = colSums(py[pnames,])
  } else {
    vec = py[place,]
  }
  fit = lm(vec ~ x)
  slopes[i,] = fit$coefficients
}

find_matches = function(place1, place2, gazetteer = gaz) {
  place = place1
  #pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
  # pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% pnames))])
  pnames = c(place, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT == place))])
  pnames = pnames[pnames %in% rownames(py)]
  pnames1 = pnames
  
  place = place2
  
  pnames = c(place, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT == place))])
  pnames = pnames[pnames %in% rownames(py)]
  pnames2 = pnames
  
  return(intersect(pnames1, pnames2))
  
}

cos_sim = function(x, y) {
  (x %*% y)/(sqrt(x %*% x) * sqrt(y %*% y))
}

ord1 = which(rowSums(pk2) > 0)
ord2 = which(rowSums(pk3) > 0)
ord = intersect(ord1, ord2)
sims = c()
for (i in 1:length(ord)) {
  print(i)
  place = rownames(pk2)[ord[i]]
  pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$SUBJECT))])
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = pnames[pnames %in% rownames(py)]
  pnames = unique(pnames)
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% pnames))])
  pnames = unique(pnames)
  if (length(pnames) > 1) {
    vec1 = colSums(pk2[pnames,])
    vec2 = colSums(pk3[pnames,])
  } else {
    vec1 = pk2[place,]
    vec2 = pk3[place,]
  }
  
  sim = cos_sim(vec1, vec2)
  sims = c(sims, sim)
}
names(sims) = rownames(pk1)[ord]

place_composition = function(place, mat) {
  pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$SUBJECT))])
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = pnames[pnames %in% rownames(mat)]
  pnames = unique(pnames)
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% pnames))])
  pnames = unique(pnames)
  pnames = pnames[pnames %in% rownames(mat)]
  if (length(pnames) > 1) {
    vec = colSums(mat[pnames,])
  } else {
    vec = mat[place,]
  }
  return(vec)
}

ord = which(rowSums(pk3) > 0)
pmat = pk3[ord,]

for (i in 1:nrow(py)) {
  print(i)
  py[i,] = place_composition(place = rownames(py)[i], mat = py)
}


# Orthogonal Procrustes
# M = Bt(A)
# R = Ut(V)
procrustes = function(A, B) {
  M = B %*% t(A)
  M_svd = irlba(M)
  U = M_svd$u
  V = M_svd$v
  R = U %*% t(V)
  return(R)
}

R1 = procrustes(A = A, B = B)
R2 = procrustes(A = B, B = C)

geodesicDistance <- function(long1, lat1, long2, lat2) {
  deg2rad <- function(deg) return(deg*pi/180)
  long1 = deg2rad(long1)
  lat1 = deg2rad(lat1)
  long2 = deg2rad(long2)
  lat2 = deg2rad(lat2)
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

semantic_footprint = function(term, mat, coords) {
  vec = mat[,term]
  vec = vec[vec > 0]
  vec = vec[names(vec) %in% coords[,"NAME"]]
  hits = which(coords[,"NAME"] %in% names(vec))
  mean_lon = sum(coords[hits,"LON"] * vec[vec > 0]) / sum(vec)
  mean_lat = sum(coords[hits,"LAT"] * vec[vec > 0]) / sum(vec)
  distances = c()
  for (i in hits) {
    lat = coords[i,"LAT"]
    lon = coords[i,"LON"]
    d = geodesicDistance(long1 = mean_lon,
                         lat1 = mean_lat,
                         long2 = lon,
                         lat2 = lat)
    distances = c(distances,d)
  }
  mean_dist = sum(distances * vec) / sum(vec)
  sd_dist = sd(distances)
  results = list()
  results$n = length(vec)
  results$freq = sum(vec)
  results$sd = sd(vec)
  results$lon = mean_lon
  results$lat = mean_lat
  results$radius = mean_dist
  results$stdev = sd_dist
  return(results)
}

footprints = matrix(0, 7, ncol(pk))
colnames(footprints) = colnames(pk)
for (j in 1:ncol(footprints)) {
  
  print(j)
  
  res = semantic_footprint(term = colnames(footprints)[j],
                           mat = pk,
                           coords = lonlat)
  res = unlist(res)
  footprints[,j] = res
}
footprints = t(footprints)
ord = which(footprints[,1] > 10)
footprints = footprints[ord,]
colnames(footprints) = names(res)

map_term = function(mat, terms, coords) {
  if (length(terms) > 1) {
    vec = rowSums(mat[,terms])
  } else {
    vec = mat[,terms]
  }
  vec = vec[vec > 0]
  vec = vec[names(vec) %in% coords[,"NAME"]]
  hits = which(coords[,"NAME"] %in% names(vec))
  lon = coords$LON[hits]
  lat = coords$LAT[hits]
  plot(lonlat$LON, lonlat$LAT, pch=20, col="gray", main = paste(terms, collapse = ","))
  points(lon, lat, pch = 20, col = "red")
}

setwd("~/Desktop/pk")
fnames = dir()
for (f in 1:length(fnames)) {
  print(f)
  load(fnames[f])
  start = (7 * (f - 1)) + 1
  end = start + 6
  footprints = matrix(0, 7, ncol(pk))
  colnames(footprints) = colnames(pk)
  for (j in 1:ncol(pk)) {
    
    print(j)
    
    res = semantic_footprint(term = colnames(footprints)[j],
                             mat = pk,
                             coords = lonlat)
    res = unlist(res)
    footprints[,j] = res
  }
  footprints = t(footprints)
  footprints = cbind(footprints, rep(f, nrow(footprints)))
  if (f == 1) {
    geosemantic_data = footprints
  } else {
    geosemantic_data = rbind(geosemantic_data, footprints)
  }
}




#### Working with geosemantic footprint data  ####

# This is an index to geo_data, which covers all five time periods
geowords = gsub("\\.0-9", "", rownames(geosemantic_data))

# Terms are the unique words, of at least 3 characters, used in at least 25 books from each period
terms = unique(geowords)
terms = terms[nchar(terms) > 2]
hits = c()
for (i in 1:length(terms)) {
  term = terms[i]
  ord = which(geowords == term)[1:5]
  hits = c(hits, all(geosemantic_data[ord,"n"] > 25))
}
terms = terms[hits]

# Displacement is the standard deviation of latitudes and longitudes for each term from period to period
displacement = c()
for (i in 1:length(terms)) {
  print(i)
  term = terms[i]
  ord = which(geowords == term)
  displacement = c(displacement, sd(geosemantic_data[ord,"lat"] * sd(geosemantic_data[ord,"lon"], na.rm = T)))
}
names(displacement) = terms
sort(displacement, decreasing = F)[1:50]

# Semantic drift measures distance from centroid in southwest europe
semantic_drift = c()
centroid = apply(geosemantic_data[,4:5], 2, median, na.rm = T)
for (i in 1:length(terms)) {
  print(i)
  term = terms[i]
  ord = which(geowords == term)
  
  lons = geosemantic_data[ord, "lon"]
  lats = geosemantic_data[ord, "lat"]
  
  d1 = geodesicDistance(long1 = mean(lons[1:3]),
                       lat1 = mean(lats[1:3]),
                       long2 = centroid[1],
                       lat2 = centroid[2])
  
  d2 = geodesicDistance(long1 = centroid[1],
                        lat1 = centroid[2],
                        long2 = mean(lons[4:5]),
                        lat2 = mean(lats[4:5]))

  semantic_drift = c(semantic_drift, d2 - d1)
}
names(semantic_drift) = terms

# Semantic spread measures the change in the radius, post 1660 vs. pre 1660
semantic_spread = c()
for (i in 1:length(terms)) {
  print(i)
  term = terms[i]
  ord = which(geowords == term)
  
  radii = geosemantic_data[ord, "radius"]
  
  r1 = mean(radii[1:3])
  r2 = mean(radii[4:5])
  
  semantic_spread = c(semantic_spread, r2 / r1)
}
names(semantic_spread) = terms

sort(displacement, decreasing = F)[1:20]
sort(displacement, decreasing = T)[1:20]
sort(semantic_drift, decreasing = F)[1:20]
sort(semantic_drift, decreasing = T)[1:20]
sort(semantic_spread, decreasing = F)[1:20]
sort(semantic_spread, decreasing = T)[1:20]

# Check word clusters -- commodities, 
ord = which(geowords %in% names(similarity(w, "liberties")))
View(geosemantic_data[ord[order(geowords[ord])],])

load("~/Downloads/py.rda")
py = py[rownames(py) %in% gaz$SUBJECT,]
for (j in 1:ncol(py)) {
  py[,j] = py[,j] / sum(py[,j])
}
py = py[,27:ncol(py)]

timeseries = function(mat, term, compose_places = T) {
  if (compose_places == T) {
    vec = place_composition(place = term, mat = mat)
  } else {
    vec = mat[term,]
  }
  vec = vec / colSums(mat)
  x = seq(from = 0, to = 1, length.out = 100)
  fit = lm(vec[1:100] ~ x)
  slope1 = round(fit$coefficients[2], digits = 3)
  fit = lm(vec[101:200] ~ x)
  slope2 = round(fit$coefficients[2], digits = 3)
  x = seq(from = 0, to = 1, length.out = 200)
  scatter.smooth(x, vec, xlab = "", xaxt = "n", pch = 20, col = "gray",
                 main = paste(term, "pre 1600:", slope1, "post 1600", slope2, collapse = " "))
  #plot(x, vec, type = "l", lwd = 2, col = "gray", xlab = "", xaxt = "n",
  #    main = paste(term, "pre 1600:", intercept, "post 1600", slope, collapse = " "))
  axis(1, at = seq(from = 0, to = 1, length.out = 3), labels = seq(from = 1500, to = 1700, length.out = 3))
  #abline(fit$coefficients)
}

composed_py = py
for (i in 1:nrow(py)) {
  print(i)
  composed_py[i,] = place_composition(place = rownames(py)[i], mat = py)
}

slopes = matrix(0, nrow(py), 2)
colnames(slopes) = c("pre","post")
rownames(slopes) = rownames(py)
totals = colSums(composed_py)
normed_py = composed_py
for (j in 1:ncol(normed_py)) {
  normed_py[,j] = normed_py[,j] / sum(normed_py[,j])
}
for (i in 1:nrow(composed_py)) {
  print(i)
  vec = normed_py[i,]
  x = seq(from = 0, to = 1, length.out = 100)
  fit = lm(vec[1:100] ~ x)
  slope1 = fit$coefficients[2]
  fit = lm(vec[101:200] ~ x)
  slope2 = fit$coefficients[2]
  slopes[i,] = c(slope1, slope2)
}

ptotals = rowSums(composed_py)
slope_totals = sum(composed_py)
slopes_adj = slopes
for (i in 1:nrow(slopes)) {
  vec = slopes[i,] / (ptotals[i] / slope_totals)
  slopes_adj[i,] = vec
}

timeseries = function(mat = py, terms) {
  PLACE = c()
  YEAR = c()
  COUNT = c()
  for (i in 1:length(terms)) {
    term = terms[i]
    vec = place_composition(place = term, mat = mat)
    place = rep(term, length(vec))
    year = names(vec)
    count = vec
    PLACE = c(PLACE, place)
    YEAR = c(YEAR,year)
    COUNT = c(COUNT,count)
  }
  df = data.frame(PLACE, YEAR, COUNT)
  ggplot(df, aes(x = YEAR, y = COUNT, color = PLACE)) + 
    geom_point() + 
    stat_smooth(method = loess, aes(group = PLACE), se = F, fullrange = F) + 
    theme(panel.background = element_blank(),
          axis.title = element_blank()) +
    scale_x_discrete(breaks=c("1500","1600","1699"),
                     labels=c("1500", "1600", "1700"))
}

timeseries(py, continents)
timeseries(py, c("india","china","persia", "turkey"))
timeseries(py, c("egypt", "israel", "greece"))
timeseries(py, c("rome"))
timeseries(py, c("edinburgh","dublin"))
timeseries(py, c("jamaica","boston","virginia", "carolina"))
timeseries(py, c("britain"))
timeseries(py, c("peru","mexico","brazil"))

### Here we go. This worked before tweaking place_composition. ###
# Place to keyword mapping
place_join = function(places, gazetteer = gaz, mode = "all") {
  places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
  places = c(places,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
  places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
  places = unique(places)
  
  if (mode == "all") {
    places = c(places, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% places))])
    places = unique(places)
  }
  return(places)
}

place_composition = function(places, mat) {
  places = places[places %in% rownames(mat)]
  if (length(places) > 1) {
    mat = mat[places,]
    mat = as.matrix(mat)
    vec = colSums(mat)
  } else {
    vec = mat[places,]
  }
  return(vec)
}


ppmi = function(mat) {
  total = sum(mat, na.rm = T)
  pcol = apply(mat, 2, function(x) sum(x, na.rm = T) / total)
  prow = apply(mat, 1, function(x) sum(x, na.rm = T) / total)
  pmat = as.matrix(prow) %*% t(pcol)
  mat = mat / total
  mat = mat / pmat
  mat = apply(mat, 1, log)
  mat = t(mat)
  mat[mat < 0] = 0
  mat[is.na(mat)] = 0
  return(mat)
}

cos_sim = function(x, y) {
  x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
}

kd = ppmi(kd)
avgs = apply(kd, 1, mean)
devs = apply(kd, 1, sd)
lsa = irlba::irlba(kd, nv = 50)
kds = lsa$u %*% diag(lsa$d)
rownames(kds) = rownames(kd)

places = place_join(c("virginia"))
vec = place_composition(places, pd)
ord = which(vec > 0)
pavgs = apply(kd[,ord], 1, mean, na.rm = T)
z = (pavgs - avgs) / devs

lsa = irlba::irlba(kd[,ord], nv = 50)
w = lsa$u %*% diag(lsa$d)
rownames(w) = rownames(kd)

m1 = w %*% t(w)


m2 = kds %*% t(kds)
sims = c()
for (i in 1:nrow(kds)) {
  sim = cos_sim(m1[i,], m2[i,])
  sims = c(sims, sim)
}
names(sims) = rownames(w)

# The lexical profile of a place is top 10 words over each metric
sort(z * sims, decreasing = T)[1:10]
x = z / abs(z)
sort(pavgs * x * (1 - sims), decreasing = T)[1:10]

# Then to compare uses of the terms
empson::similarity(w, "ocean")
empson::similarity(kds, "ocean")

summary(sims)
