isoutlier <- function (data, n) {
  m <- median(data, na.rm = T)
  md <- mad(data, na.rm = T)
  ll <- m - n * md
  ul <- m + n * md
  output <- data < ll | data > ul
  return(output)
}