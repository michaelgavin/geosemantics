#### CHAPTER ONE ANALYSIS SCRIPT ####
library(geosemantics)
library(ggplot2)
load("~/Desktop/geosemantics/data/ch1_place_year.rda")
toponyms = rownames(py)
toponyms = toponyms[toponyms %in% stops == F]
py = py[toponyms, 27:226]
totals = rowSums(py)

excluded_terms = c()
for (i in 1:length(toponyms)) {
  print(i)
  toponym = toponyms[i]
  sames = place_join(places = toponym, gaz = gaz, mode = "same")
  if (length(sames) == 1) {
    next()
  }
  if (length(sames) > 1) {
    top_place = names(sort(totals[sames], decreasing = T))[1]
    if (top_place != toponym) {
      excluded_terms = c(excluded_terms, toponym)
    }
  }
}

places = setdiff(toponyms, excluded_terms)
composed_py = matrix(0, length(places), 200)
rownames(composed_py) = places
colnames(composed_py) = colnames(py)
for (i in 1:length(places)) {
  print(i)
  place = places[[i]]
  sames = place_join(place, gaz, mode = "same")
  vec = place_composition(py, sames)
  composed_py[i,] = vec
}
totals = rowSums(composed_py)
totals = totals[names(totals) != "world"]
places = places[places != "world"]

ranges = list(1:100,101:140,141:160,161:180,181:200)

for (i in 1:length(ranges)) {
  subtotals = rowSums(composed_py[places,ranges[[i]]])
  subtotals = sort(subtotals, decreasing = T)[1:5]
  subtotals = round(subtotals / 1000, digits = 1)
  if (i == 1) {
    top_five = matrix(names(subtotals),5,1)
    top_five = cbind(top_five, subtotals)
  } else {
    top_five = cbind(top_five,names(subtotals))
    top_five = cbind(top_five, subtotals)
  }
}

## Number of places per year
vec = apply(composed_py, 2, function(x) { length(x[x>0])})
plot(vec)
plot_timeseries(composed_py, places = continents)
plot_semantic_map(composed_py, term = colnames(composed_py)[1], coords = lonlat)
plot_semantic_map(composed_py, term = colnames(composed_py)[50], coords = lonlat)
plot_semantic_map(composed_py, term = colnames(composed_py)[100], coords = lonlat)

freq_dist = table(totals)
y = as.numeric(freq_dist)
x = as.numeric(names(freq_dist))
plot(x, y, log = "x")
abline(v = median(totals), col = "red")


### Need to find slopes analysis and paste here
sort(totals[slopes[places,"post"] < 0], decreasing = T)[1:10]

## Build PK
## Need to build composed_pd matrices for each place
for (f in 1:5) {
  print(f)
  fp = paste("~/Desktop/geosemantics/inst/extdata/ch1_pd",f,".rda", sep="")
  load(fp)
  fp = paste("~/Desktop/geosemantics/inst/extdata/ch1_kd",f,".rda", sep="")
  load(fp)
  pk = pd %*% t(kd)
  if (f == 1) {
    pk = pk
  } else {
    pk = pk + pk
  }
}

pk = as.matrix(pk)
composed_pk = matrix(0, length(places), ncol(pk))
rownames(composed_pk) = places
for (i in 386:length(places)) {
  print(i)
  place = places[i]
  sames = place_join(place, gaz, mode = "same")
  vec = place_composition(pk, sames)
  if (sum(vec) == 0) { next }
  composed_pk[place,] = vec
}
save(composed_pk, file = "~/Desktop/composed_pk.rda")

totals = rowSums(composed_pk)

# Get centroid of all PK
vec = totals
coords = ch1$coordinates
vec = vec[names(vec) %in% coords[, "NAME"]]
hits = which(coords[, "NAME"] %in% names(vec))
vec = vec[coords[hits, "NAME"]]
lon = coords[hits, "LON"]
lon[lon >= 181] = 181 - lon[lon >= 181]
lat = coords[hits, "LAT"]
mean_lon = sum(lon * vec) / sum(vec)
mean_lat = sum(lat * vec) / sum(vec)

unweighted_mean_lon = mean(lon)
unweighted_mean_lat = mean(lat)

# Get places closest to centroid
d = c()
for (i in 1:length(lon)) {
  loni = lon[i]
  lati = lat[i]
  distance = great_circle_distance(lon1 = loni, lat1 = lati, lon2 = mean_lon, lat2 = mean_lat)
  d = c(d, distance)
}

# Get places closest to unweighted centroid
un_d = c()
for (i in 1:nrow(coords)) {
  loni = lon[i]
  lati = lat[i]
  distance = great_circle_distance(lon1 = loni, lat1 = lati, lon2 = unweighted_mean_lon, lat2 = unweighted_mean_lat)
  un_d = c(un_d, distance)
}

great_circle_distance(lon1 = mean_lon,
                      lat1 = mean_lat,
                      lon2 = unweighted_mean_lon, 
                      lat2 = unweighted_mean_lat)
lon = coords[,"LON"]
lat = coords[,"LAT"]
ord = which(lon < 60 & lat > 30)
plot(lon[ord], lat[ord], pch = 20, col = "gray")
points(unweighted_mean_lon, unweighted_mean_lat, pch = 18, col = "black", cex = 2)
points(mean_lon, mean_lat, pch = 18, col = "black", cex = 2)

## Get centroid, range, and variation over composed_pk
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
colnames(footprints) = names(res)

d = c()
for (i in 1:nrow(footprints)) {
  print(i)
  lon = footprints[i,"lon"]
  lat = footprints[i,"lat"]
  distance = great_circle_distance(lon1 = lon, lat1 = lat, lon2 = mean_lon, lat2 = mean_lat)
  d = c(d, distance)
}

x = d
y = footprints[,"freq"]
plot(x, y, log = "xy", main = "distance from centroid")
abline(h = mean(y), col = "red", lwd = 2)
abline(v = mean(x), col = "red", lwd = 2)

x = footprints[,"radius"]
y = footprints[,"freq"]
plot(x, y, log = "xy", main = "range")
abline(h = mean(y), col = "red", lwd = 2)
abline(v = mean(x), col = "red", lwd = 2)

### Compare spatial distributions of terms
x = d
y = footprints[,"radius"]
plot(x, y, cex = 0, main = "x = distance, y = stdev")
text(x, y, labels = names(y))
abline(h = mean(y), col = "red", lwd = 2)
abline(v = mean(x), col = "red", lwd = 2)

x = d
y = footprints[,"stdev"]
fit = lm(y ~ x)
adj_y = (fit$coefficients[2] * x) + fit$coefficients[1]
adj_y = adj_y - 50
#plot(x, y, cex = 0, main = "x = distance, y = stdev")
#text(x, y, labels = names(y))
plot(x, y, main = "x = distance, y = stdev")
abline(h = mean(y), col = "red", lwd = 2)
abline(v = mean(x), col = "red", lwd = 2)
lines(x, adj_y, pch=20, col = "red", lwd = 4)

ord = which(fit$residuals < -50)
spatial_terms = rownames(footprints)[ord]
