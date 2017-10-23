infer.missing.values <- function(data.array) {
  out <- data.array
  data.dim <- dim(data.array)[-1]
  template <- array(0, data.dim)
  #slice.index gives the indexing for each element for slicing along the given dimension,
  # e.g. entries in slice.array(data.array, 1) will match the row number
  for (j in 1:prod(data.dim)) {
    arr_ind <- vec_to_arr_index(j, data.dim)
    weight.array <- array(0, data.dim)
    for (dim.num in seq_along(data.dim)) {
      weight.array <- weight.array - (slice.index(template, dim.num) - arr_ind[dim.num])^2 / 2
    }
    weight.array <- exp(weight.array)
    for (i in 1:dim(data.array)[1]) {
      current.point <- array(data.array[slice.index(data.array, 1) == i], data.dim)
      if (is.na(current.point[matrix(arr_ind, 1)])) {
        out[matrix(c(i, arr_ind), 1)] <- sum(weight.array * current.point, na.rm=TRUE) / sum(as.numeric(!is.na(current.point)) * weight.array)
      }
    }
  }
  return(out)
}
