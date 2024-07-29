#' Title
#'
#' @param delay_time
#' @param cut_date
#' @param rate_0
#' @param rate_1
#'
#' @return
#' @export
#'
#' @examples
dat_pw_surv <- function(delay_time = 4, cut_date = 28.6,
                        rate_0 = log(2) / c(10, 10),
                        rate_1 = log(2) / c(10, 10/0.6)){
  # Inputs: T ~ Exp(rate_0) for t<delay_time and Exp(rate_1) for t>delay_time
  # Return: the survival function for T between 0 and cut_date

  t1 <- seq(0, delay_time, 0.1)
  t2 <- seq(max(t1)+0.1, cut_date, 0.1)

  surv_0 <- c(surv_exp(t1, rate = rate_0[1]), surv_exp(delay_time, rate = rate_0[1])*surv_exp(t2-delay_time, rate = rate_0[2]))
  # all.equal(surv_0, surv_exp(c(t1,t2), rate = rate_0[1]))
  surv_1 <- c(surv_exp(t1, rate = rate_1[1]), surv_exp(delay_time, rate = rate_1[1])*surv_exp(t2-delay_time, rate = rate_1[2]))

  dat <- data.frame(time = rep(c(t1,t2), 2),
                    trt = c(rep("Control", length(c(t1,t2))), rep("Experimental", length(c(t1,t2)))),
                    surv = c(surv_0, surv_1)) %>%
    mutate(trt = factor(trt))
  return(dat)
}
