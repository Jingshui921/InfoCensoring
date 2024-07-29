#' The quantile function for piecewise exponential distribution
#' T ~ Exp(lambda_1),         0 < t < t_1,
#' T ~ Exp(lambda_2),       t_1 < t < t_1 + t_2,
#' T ~ Exp(lambda_2), t_1 + t_2 < t < t_1 + t_2 + t_3,
#' ....
#' @param P vector of probabilities
#' @param duration vector of durations (t_1, t_2, t_3, ...) where T ~ Exp(lambda_k), t_1 + ... + t_{k-1} < t < t_1 + ... + t_k
#' @param rate vector of rates (lambda_1, lambda_2, lambda_3, ...)
#' @param lower.tail logical; if TRUE (default), probabilities are CDF P(T<=q), otherwise, survival P(T>q)
#'
#' @return t where P(T <= t) = q if lower.tail =TRUE; or t where P(T > t) = q if lower.tail =FALSE;
#' @export
#'
#' @examples
# qpwexp(p = c(0,
#        pexp(2, rate = log(2)/10),
#        seq(0.1, 1, 0.1)),
#        duration = c(4,4,100),
#        rate = c(log(2)/10, 0.6*log(2)/10, 0.8*log(2)/10))


qpwexp <- function(p, duration, rate, lower.tail =TRUE){
  # quantile qexp(pexp(6, rate = 1), rate = 1)=6

  if (lower.tail == FALSE) {p = 1 - p}

  pw_cum_haz <- c(0, cumsum(duration*rate))
  cum_haz <- -log(1-p)

  index <- findInterval(cum_haz, pw_cum_haz, all.inside = TRUE)
  cut_points <- c(0, cumsum(duration))

  x <- cut_points[index] + ((cum_haz - pw_cum_haz[index])/rate[index])

  return(x)
}
