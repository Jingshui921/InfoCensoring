#' The density function for piecewise exponential distribution
#' T ~ Exp(lambda_1),         0 < t < t_1,
#' T ~ Exp(lambda_2),       t_1 < t < t_1 + t_2,
#' T ~ Exp(lambda_2), t_1 + t_2 < t < t_1 + t_2 + t_3,
#' ....
#' @param x vector of quantiles.
#' @param duration (t_1, t_2, t_3, ...) where T ~ Exp(lambda_k), t_1 + ... + t_{k-1} < t < t_1 + ... + t_k
#' @param rate vector of rates (lambda_1, lambda_2, lambda_3, ...)
#'
#' @return the density function of T at t = x
#' @export
#'
#' @examples
#' dpwexp(x = c(0,4,6,8,12,Inf),
#'        duration = c(4,4,100),
#'        rate = c(log(2)/10, 0.6*log(2)/10, 0.8*log(2)/10))

dpwexp <- function(x, duration, rate){
  # density f(t) = h(t)*S(t)

  cut_points <- c(0, cumsum(duration))
  index <- findInterval(x, cut_points, all.inside = TRUE)

  return(rate[index]*ppwexp(x, duration, rate, lower.tail = FALSE))
}
