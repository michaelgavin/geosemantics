#### CHAPTER ONE ANALYSIS SCRIPT ####
all_toponyms = rownames(py)
toponyms = all_toponyms
unique_places = list()
totals = rowSums(py)
for (i in 1:length(all_toponyms)) {
  print(i)
  toponym = all_toponyms[i]
  if (toponym %in% toponyms == F) {
    next 
  } else {
      sames = place_join(places = toponym, gaz = gaz, mode = "same")
      if (length(sames) == 1) {
        unique_places[[i]] = toponym
        names(unique_places)[i] = toponym
      }
      if (length(sames) > 1) {
        unique_places[[i]] = sames
        names(unique_places)[i] = names(sort(totals[sames], decreasing = T))[1]
      }
    }
}
unique_places = unique_places[-which(is.na(names(unique_places)))]
toponyms = unique(names(unique_places))
places_list = unique_places[toponyms]
composed_py = matrix(0, length(toponyms), 200)
rownames(composed_py) = toponyms
colnames(composed_py) = colnames(py)
for (i in 1:length(toponyms)) {
  print(i)
  place = toponyms[i]
  vec = place_composition(py, places_list[[i]])
  composed_py[i,] = vec
}
totals = rowSums(composed_py)
totals = totals[names(totals) != "world"]

## Number of places per year
places = rownames(composed_py)

vec = apply(composed_py, 2, function(x) { length(x[x>0])})
plot_timeseries(composed_py, places = continents)
plot_semantic_map(composed_py, term = colnames(composed_py)[1], coords = lonlat)
plot_semantic_map(composed_py, term = colnames(composed_py)[50], coords = lonlat)
plot_semantic_map(composed_py, term = colnames(composed_py)[100], coords = lonlat)

y = as.numeric(freq_dist)
x = as.numeric(names(freq_dist))
plot(x, y, log = "x")
abline(v = median(totals), col = "red")

sort(totals[slopes[places,"post"] < 0], decreasing = T)[1:10]
