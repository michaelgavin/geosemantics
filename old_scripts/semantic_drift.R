loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

cos_sim = function(x, y) {
  x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
}

svd_loc = "C:/Users/mgavin/Desktop/svd/"
output_loc = "C:/Users/mgavin/Desktop/output/"
filepaths = paste(svd_loc,dir(svd_loc),sep="")

for (f in 1:4) {
  print(f)
  fp1 = filepaths[f]
  fp2 = filepaths[f+1]
  
  print("loading data...")
  m1 = loadRData(fp1)
  m2 = loadRData(fp2)
  
  A = m1$u %*% diag(m1$d)
  B = m2$u %*% diag(m2$d)
  print("multiplying matrices...")
  M = B %*% t(A)
  
  save(M, file = paste(output_loc, "M", f, ".rda", sep = ""))
  
}
rm(A, B)

keywords_displacement = matrix(0,length(keywords),3)
rownames(keywords_displacement) = keywords

for (f in 1:3) {
  
  m1 = loadRData(dir()[f])
  m2 = loadRData(dir()[f+1])
  
  sims = c()
  for (i in 1:nrow(m1)) {
    print(paste(f, i))
    sims = c(distances, cos_sim(m1[i,], m2[i,]))
  }
  keywords_sims[,f] = sims
  
}

distances = c()
for (i in 1:nrow(keywords_displacement)) {
  
  d = sum(keywords_displacement[i,])
  distances = c(distances, d)
}
names(distances) = keywords

# Move M files to other folder before running next step

# Now get places
for (f in 6:9) {
  print(f)
  fp1 = filepaths[f]
  fp2 = filepaths[f+1]
  
  print("loading data...")
  m1 = loadRData(fp1)
  m2 = loadRData(fp2)
  
  A = m1$u %*% diag(m1$d)
  B = m2$u %*% diag(m2$d)
  print("multiplying matrices...")
  M = A %*% t(B)
  
  save(M, file = paste(output_loc, "M", f, ".rda", sep = ""))
  
}
rm(A, B)

for (f in 11:14) {
  print(f)
  fp1 = filepaths[f]
  fp2 = filepaths[f+1]
  
  print("loading data...")
  m1 = loadRData(fp1)
  m2 = loadRData(fp2)
  
  A = m1$u %*% diag(m1$d)
  B = m2$u %*% diag(m2$d)
  print("multiplying matrices...")
  M = A %*% t(B)
  
  save(M, file = paste(output_loc, "M", f, ".rda", sep = ""))
  
}
rm(A, B)

place_displacement = matrix(0,length(places),3)
rownames(place_displacement) = places

for (f in 1:3) {
  hit = 1 + f
  
  m1 = loadRData(dir()[hit])
  m2 = loadRData(dir()[hit+1])
  
  rownames(m1) = rownames(pd)
  rownames(m2) = rownames(pd)
  
  m1 = m1[places,]
  m2 = m2[places,]
  
  m2 = m1 + m2
  
  sims = c()
  for (i in 1:nrow(m1)) {
    print(paste(i))
    place = places[i]
    pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
    pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$SUBJECT))])
    pnames = pnames[pnames %in% places]
    pnames = unique(pnames)
    if (length(pnames) > 1) {
      vec1 = colSums(m1[pnames,])
      vec2 = colSums(m2[pnames,])
    } else {
      vec1 = m1[place,]
      vec2 = m2[place,]
    }

    sims = c(sims, cos_sim(vec1, vec2))
  }
  place_displacement[,f] = sims
  
}
# ord = apply(place_displacement, 1, function(x) any(is.na(x)))
# place_displacement = place_displacement[ord == F,]
# pnames = unique(gaz$SUBJECT)
# place_displacement = place_displacement[which(rownames(place_displacement) %in% pnames == T),]

distances = c()
for (i in 1:nrow(place_displacement)) {
  
  d = sum(1 - abs(place_displacement[i,]))
  distances = c(distances, d)
}
names(distances) = rownames(place_displacement)

total_variance = c()
for (i in 1:length(distances)) {
  print(i)
  place = names(distances)[i]
  pnames = c(place,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$OBJECT))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), grep(place,gaz$SUBJECT))])
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = pnames[pnames %in% names(distances)]
  pnames = unique(pnames)
  pnames = c(pnames,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% pnames))])
  pnames = c(pnames,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% pnames))])
  pnames = unique(pnames)
  pnames = pnames[pnames %in% names(distances)]
  total_variance = c(total_variance,mean(distances[pnames], na.rm = T))
}
names(total_variance) = names(distances)

locations = 



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


### Here we go. This worked before tweaking place_composition. ###
# Place to keyword mapping
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
    row[is.na(row)] = 0
    row[is.infinite(row)] = 0
    mat[i,] = row
  }

  return(mat)
}


kd = ppmi(kd)
avgs = apply(kd, 1, mean)
devs = apply(kd, 1, sd)
lsa = irlba::irlba(kd, nv = 50)
kds = lsa$u %*% diag(lsa$d)
rownames(kds) = rownames(kd)

places = place_join(c("virginia"), pd)
vec = place_composition(places, pd)
ord = which(vec > 0)
pavgs = apply(kd[,ord], 1, mean, na.rm = T)
z = (pavgs - avgs) / devs

lsa = irlba::irlba(kd[,ord], nv = 50)
w = lsa$u %*% diag(lsa$d)
rownames(w) = rownames(kd)


sims = c()
for (i in 1:nrow(kds)) {
  sim = cos_sim(w[i,], kds[i,])
  sims = c(sims, sim)
}
names(sims) = rownames(w)

# The lexical profile of a place is top 10 words over each metric
sort(z * sims, decreasing = T)[1:10]
x = z / abs(z)
sort(pavgs * x * (1 - sims), decreasing = T)[1:10]




# Match gaz to lonlat
NAME = c()
LON = c()
LAT = c()
places = unique(gaz$SUBJECT)
continents = unique(gaz$SUBJECT[gaz$OBJECT == "continent"])
continental_places = gaz$SUBJECT[which(gaz$PREDICATE == "is in" & 
                                         gaz$OBJECT %in% continents == T)]
for (i in 1:length(places)) {
  print(i)
  place = places[i]
  
  if (place %in% lonlat$NAME) {
    
    hits = which(lonlat$NAME == place)
    lon = mean(lonlat$LON[hits])
    lat = mean(lonlat$LAT[hits])
    
    NAME = c(NAME, place)
    LON = c(LON, lon)
    LAT = c(LAT, lat)
    next()
  }
  
  pnames = place_join(pname = place, mode = "NA")
  
  if (any(pnames %in% lonlat$NAME)) {
    
    hits = which(lonlat$NAME %in% pnames)
    lon = mean(lonlat$LON[hits])
    lat = mean(lonlat$LAT[hits])
    
    NAME = c(NAME, place)
    LON = c(LON, lon)
    LAT = c(LAT, lat)
    next()
  }
  
  pnames = place_join(pname = place)
  
  if (any(pnames %in% lonlat$NAME)) {
    
    hits = which(lonlat$NAME %in% pnames)
    lon = mean(lonlat$LON[hits])
    lat = mean(lonlat$LAT[hits])
    
    NAME = c(NAME, place)
    LON = c(LON, lon)
    LAT = c(LAT, lat)
    next()
  }
  
  if (any(pnames %in% continental_places)) {
    
    pnames = gaz$OBJECT[which(gaz$SUBJECT %in% pnames == T & gaz$PREDICATE == "is in")]
    
    if (any(pnames %in% lonlat$NAME)) {
      
      hits = which(lonlat$NAME %in% pnames)
      lon = mean(lonlat$LON[hits])
      lat = mean(lonlat$LAT[hits])
      
      NAME = c(NAME, place)
      LON = c(LON, lon)
      LAT = c(LAT, lat)
      next()
    }
  }
}
matched_lonlat = data.frame(NAME, LON, LAT)
write.csv(matched_lonlat, file = "/Users/mgavin/Desktop/matched_lonlat.csv")

# Now geocorrelation
lonlat = read.csv("/Users/mgavin/Downloads/matched_lonlat.csv - matched_lonlat.csv.csv", stringsAsFactors = F)
locations = unique(lonlat$NAME)
locations = locations[locations %in% gaz$SUBJECT]


p = m_reduced$u %*% diag(m_reduced$d)
rownames(p) = pnames
locations = locations[locations %in% rownames(p)]
p = p[locations,]

lonlat = lonlat[lonlat$NAME %in% locations,]
rownames(lonlat) = lonlat$NAME
lonlat = as.matrix(lonlat[,2:3])
lonlat = lonlat[locations,]

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

dmat = matrix(0, length(locations), length(locations))
rownames(dmat) = locations
for (i in 1:nrow(dmat)) {
  print(i)
  loc = locations[i]
  d = geodesicDistance(long1 = lonlat[loc,1], 
                       lat1 = lonlat[loc,2],
                       long2 = lonlat[,1],
                       lat2 = lonlat[,2])
  dmat[i,] = d
}
dmat[is.na(dmat)] = 0

smat = matrix(0, length(locations), length(locations))
rownames(smat) = locations
for (i in 1:nrow(smat)) {
  print(i)
  loc = locations[i]
  d = 1 - similarity(p, loc, fullResults = T)
  smat[i,] = d
}
smat[is.na(smat)] = 0
dmat[is.na(dmat)] = 0


cos_sim = function(x, y) {
  x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
}
sims = c()
for (i in 1:nrow(dmat)) {
  print(i)
  sim = cos_sim(smat[i,], dmat[i,])
  sims = c(sims, sim)
}
names(sims) = locations
ord = which(is.na(sims) == F)
plot(sims[ord], dmat["paris",ord])

density = apply(dmat, 1, function(x) length(x[x < 300]) / length(x))

fit = lm(sims[ord] ~ log(d + .001))
y=predict(fit,newdata=list(x=d))
plot(d, sims[ord])
points(d,y, pch = 20, col = "red")
