#' Title
#'
#' @param cutpoint
#' @param prob.dpo
#'
#' @return
#' @export
#'
#' @examples
pw_dpo <- function(cutpoint = c(0, 6, 12, 18, 24, 30, 36),
                   prob.dpo = c(0, .12, .18, .20, .30, 0, 0)) {

  # T - time to dropout(dpo) follows piecewise-exponential distribution
  # where T|t_k-1 < t < t_k ~ Exp(lambda_k), k = 0,1,2,..., t_0=0

  # rate parameter was determined by the survival probabilities at chosen cutpoints
  # P(T <= t_k) = p_k

  # t_k = cutpoint[k]
  # p_k = prob.dpo[k]

  lambda <- sapply(2:length(cutpoint),
                   function(k) ifelse(prob.dpo[k]==0, Inf,
                                      -log(1- (prob.dpo[k] - prob.dpo[k-1]))/(cutpoint[k] - cutpoint[k-1])))

  return(data.frame(cutpoint = cutpoint[-1],
                    prob.dpo = prob.dpo[-1],
                    lambda))
}
