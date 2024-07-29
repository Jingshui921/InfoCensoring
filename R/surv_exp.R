#' Survival function for Exp(rate) at T=t
#'
#' @param t
#' @param rate
#'
#' @return
#' @export
#'
#' @examples
surv_exp <- function(t, rate) {
  pexp(t, rate = rate, lower.tail = FALSE)
}
