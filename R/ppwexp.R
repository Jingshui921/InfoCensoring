#' The distribution function for piecewise exponential distribution
#' T ~ Exp(lambda_1),         0 < t < t_1,
#' T ~ Exp(lambda_2),       t_1 < t < t_1 + t_2,
#' T ~ Exp(lambda_2), t_1 + t_2 < t < t_1 + t_2 + t_3,
#' ....
#' @param q vector of quantiles.
#' @param duration (t_1, t_2, t_3, ...) where T ~ Exp(lambda_k), t_1 + ... + t_{k-1} < t < t_1 + ... + t_k
#' @param rate vector of rates (lambda_1, lambda_2, lambda_3, ...)
#' @param lower.tail logical; if TRUE (default), probabilities are CDF P(T<=q), otherwise, survival P(T>q)
#'
#' @return the cumulative distribution function of T at T=q if lower.tail = TRUE;
#' the survival function of T at T=q if lower.tail = FALSE
#' @export
#'
#' @examples
#' ppwexp(q = c(0,4,6,8,12,Inf),
#'        duration = c(4,4,100),
#'        rate = c(log(2)/10, 0.6*log(2)/10, 0.8*log(2)/10))

ppwexp <- function(q, duration, rate, lower.tail = TRUE){
  # distribution F(t) = 1-S(t) = 1-exp(-H(t))

  cut_points <- c(0, cumsum(duration))
  index <- findInterval(q, cut_points, all.inside = TRUE) # cdf is right-continuous

  pw_cum_haz <- c(0, cumsum(duration*rate))
  cum_haz <- pw_cum_haz[index]+rate[index]*(q-cut_points[index])

  surv <- exp(-cum_haz)

  if(lower.tail) {
    return(1-surv)
  } else{
    return(surv)
  }

}
