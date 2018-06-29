#' 
#' Plots all places where a key term is used. Computed over place-keyword matrix.
#'
#' @export
plot_semantic_map = function(mat, terms, coords) {
  LON = coords$LON
  LAT = coords$LAT
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
  lon[lon >= 181] = 181 - lon[lon >= 181]
  LON[LON >= 181] = 181 - LON[LON >= 181]
  plot(LON, LAT, pch=20, col="gray", main = paste(terms, collapse = ","))
  points(lon, lat, pch = 20, col = "red")
}
