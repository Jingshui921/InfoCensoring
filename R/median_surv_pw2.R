#' Calculate the median(T) for T ~ Exp(rate.phase1) for t<delay and Exp(rate.phase2) for t>delay.
#' Say delay = 4, and med.T > 4, then P(T<med.T)=P(T<4)+P(T>=4) * P(4 <= T < med.T)=P(T<4)+P(T>=4)*P(0 <= T-4 < med.T-4)
#'
#' @param delay.time T ~ Exp(rate.phase1) for t < delay.time
#' @param rate.phase1 T ~ Exp(rate.phase1) for t < delay.time
#' @param rate.phase2 T-delay.time ~ Exp(rate.phase2) for t >= delay.time
#'
#' @return med.T so that P(T<med.T) = 0.5
#' @export
#'
#' @examples
#'  median_surv_pw2(rate.phase1 = log(2) / 10,
#'                  rate.phase2 = 0.6 * log(2) / 10,
#'                  delay.time = 4)
#'  # The median survival time is 14 for T ~ Exp(log(2)/10) for t < 4 and T-4 ~ Exp(0.6*log(2)/10) for t>=4.
#'
#'  median_surv_pw2(rate.phase1 = log(2) / 10,
#'                  rate.phase2 = 0.6 * log(2) / 10,
#'                  delay.time = 15)
#'  # The median survival is 10 for T ~ Exp(log(2)/10) for t < 15 and T-15 ~ Exp(0.6*log(2)/10) for t>=15.
#'  The warning message reminds you that the median PFS is witin the delay window, med.T < delay.time.

median_surv_pw2 <- function(rate.phase1 = log(2)/c(4, 10, 20),
                            rate.phase2 = 0.6*log(2)/c(4, 10, 20),
                            delay.time = c(2, 4, 6)){

  prob.ev.phase1 <- pexp(delay.time, rate = rate.phase1)

  if(sum(0.5-prob.ev.phase1 < 0) > 0){
    warning("median PFS less than delay.time")
    med.T <- qexp(p = 0.5, rate = rate.phase1)
  } else {
    time2med.phase2 <- qexp((0.5-prob.ev.phase1)/(1-prob.ev.phase1), rate = rate.phase2)
    med.T <- delay.time+time2med.phase2
  }
  return(med.T)
}
