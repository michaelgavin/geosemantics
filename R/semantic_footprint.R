#' 
#' Measures the semantic footprint of a term (computed over place-keyword matrix)
#'
#' @export
semantic_footprint = function(mat, term, coords) {
  vec = mat[,term]
  vec = vec[vec > 0]
  vec = vec[names(vec) %in% coords[,"NAME"]]
  hits = which(coords[,"NAME"] %in% names(vec))
  vec = vec[coords[hits,"NAME"]]
  lat = coords[hits,"LAT"]
  lon = coords[hits,"LON"]
  lon[lon >= 180] = lon[lon >= 180] - 360
  mean_lon = sum(lon * vec[vec > 0]) / sum(vec)
  mean_lat = sum(lat * vec[vec > 0]) / sum(vec)
  lon = coords[hits,"LON"]
  distances = c()
  for (i in 1:length(lon)) {
    lati = lat[i]
    loni = lon[i]
    d = great_circle_distance(lon1 = mean_lon,
                         lat1 = mean_lat,
                         lon2 = loni,
                         lat2 = lati)
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
