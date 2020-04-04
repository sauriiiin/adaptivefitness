isoutlier <- function (data) {
  m <- median(data, na.rm = T)
  md <- mad(data, na.rm = T)
  ll <- m - 3 * md
  ul <- m + 3 * md
  output <- data < ll | data > ul
  return(output)
}