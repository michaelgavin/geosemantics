#' Compute similarity over a matrix
#' 
#' @param mat A matrix
#' @param term A word (must be in the rownames of mat)
#' @param output "top" or "all"
#' @return similarity scores
#' 
#' @export
semantic_similarity = function(mat, term, output = "top") {
  cos_sim = function(x, y) {
    x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  }
  if (term %in% rownames(mat) == F) {
    stop("The term is not among the rownames of your matrix. This can happen
         for several reasons: 1) the term is just not included in your
         data, 2) the matrix does not have rownames, or 3) you are trying
         to calculate over the columns, in which case you need to transpose
         the matrix with t(mat).")
  }
  vec = mat[term,]
  results = unlist(apply(mat, 1, cos_sim, vec))
  if (output == "top") {
    results = sort(results, decreasing = T)[1:12]
  }
  return(results)
}
