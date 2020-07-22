

#' @export
#' @useDynLib HiddenSemiMarkov add_
addjbdjkbsdbwk <- function(x, y) {
  .C(add_, x, y, numeric(1))[[3]]
}
