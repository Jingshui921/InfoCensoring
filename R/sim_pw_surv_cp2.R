#' Title
#' @importFrom tibble tibble
#' @importFrom simtrial randomize_by_fixed_block rpw_enroll rpwexpinvRcpp
#' @importFrom copula rotCopula gumbelCopula claytonCopula
#' @import dplyr
#'
#' @param n Number of observations.
#'   If length(n) > 1, the length is taken to be the number required.
#' @param stratum A tibble with stratum specified in `stratum`,
#'   probability (incidence) of each stratum in `p`.
#' @param block Vector of treatments to be included in each block.
#' @param enroll_rate Enrollment rates; see details and examples.
#' @param fail_rate Failure rates; see details and examples;
#'   note that treatments need to be the same as input in block.
#' @param dropout_rate Dropout rates; see details and examples;
#'   note that treatments need to be the same as input in block.
#' @param copula The copula model to use. Currently only accept inputs from c('clayton', 'gumbel')
#' @param theta_trt The theta parameter that characterizes the correlation between failure time and dropout time in the experimental arm in the copula model.
#' @param theta_ctr The theta parameter that characterizes the correlation between failure time and dropout time in the control arm in the copula model.
#'
#' @return A `tibble` with the following variables for each observation:
#' - `stratum`.
#' - `enroll_time`: Enrollment time for the observation.
#' - `Treatment`: Treatment group; this will be one of the values
#'   in the input `block`.
#' - `fail_time`: Failure time generated using [rpwexp()].
#' - `dropout_time`: Dropout time generated using [rpwexp()].
#' - `cte`: Calendar time of enrollment plot the minimum of
#'   failure time and dropout time.
#' - `fail`: Indicator that `cte` was set using failure time;
#'   i.e., 1 is a failure, 0 is a dropout.
#' @export
#'
#' @examples
#'
#' # Example 1: unstratified design
#' n <- 700
#' enroll_duration <- c(2, 2, 8)
#' enroll_relative_speed <- c(1, 2, 4)
#' q <- 700 / sum(enroll_duration*enroll_relative_speed)
#' enroll_rate <- tibble(stratum = "All",
#'                       duration = enroll_duration, rate = q*enroll_relative_speed)
#' block <- c(rep("control", 2), rep("experimental", 2))
#' stratum <- tibble(stratum = "All", p = 1)
#' fail_rate <- tibble(stratum = "All",
#'                     period = rep(1:2, 2),
#'                     treatment = c(rep("control", 2), rep("experimental", 2)),
#'                     duration = rep(c(4, 999), 2),
#'                     rate = log(2) / c(10, 10, 10, 10/0.6))
#' dropout_rate <- tibble(stratum = "All",
#'                        period = rep(1:2, 2),
#'                        treatment = c(rep("control", 2), rep("experimental", 2)),
#'                        duration = rep(c(12, 999),2),
#'                        rate = c(-log(1-0.25)/12, 0.00000001, -log(1-0.20)/12, 0.00000001))
#' # Example 1a: non-informative censoring in both groups
#' sim_pw_surv_cp2(n = n,
#'                stratum = stratum,
#'                block = block,
#'                enroll_rate = enroll_rate,
#'                fail_rate = fail_rate,
#'                copula = 'clayton',
#'                theta_trt = 0,
#'                theta_ctr = 0)
#'
#' # Example 1b: positively correlated in both groups
#' sim_pw_surv_cp2(n = n,
#'                stratum = stratum,
#'                block = block,
#'                enroll_rate = enroll_rate,
#'                fail_rate = fail_rate,
#'                dropout_rate = dropout_rate,
#'                copula = 'clayton',
#'                theta_trt = 1,
#'                theta_ctr = 1)
#'
#' # Example 2: stratified design
#' n <- 700
#' enroll_rate <- tibble(stratum = rep(c("Biomarker positive", "Biomarker negative"), each = 3),
#'                       duration = rep(c(2, 2, 8), 2), rate = c(1, 2, 4, 2, 1, 4))
#' block <- c(rep("control", 2), rep("experimental", 2))
#' stratum <- tibble(stratum = c("Biomarker positive", "Biomarker negative"), p = c(0.6, 0.4))
#' fail_rate <- tibble(stratum = rep(c("Biomarker positive", "Biomarker negative"), each = 4),
#'                     period = rep(1:2, 4),
#'                     treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
#'                     duration = c(rep(c(4, 999), 2), rep(c(6, 999), 2)),
#'                     rate = log(2) / c(10, 10, 10, 10/0.6, 10, 10, 10, 10/0.8))
#' dropout_rate <- tibble(stratum = rep(c("Biomarker positive", "Biomarker negative"), each = 4),
#'                        period = rep(1:2, 4),
#'                        treatment = rep(c(rep("control", 2), rep("experimental", 2)), 2),
#'                        duration = rep(c(12, 999),4),
#'                        rate = c(-log(1-0.25)/12, 1e-8, -log(1-0.20)/12, 1e-8,
#'                                 -log(1-0.35)/12, 1e-8, -log(1-0.30)/12, 1e-8))
#' sim_pw_surv_cp2(n = n,
#'                stratum = stratum,
#'                block = block,
#'                enroll_rate = enroll_rate,
#'                fail_rate = fail_rate,
#'                dropout_rate = dropout_rate,
#'                copula = 'clayton',
#'                theta_trt = 1,
#'                theta_ctr = 1)

sim_pw_surv_cp2 <- function(n,
                           stratum,
                           block,
                           enroll_rate,
                           fail_rate,
                           dropout_rate,
                           copula = c('clayton', 'gumbel'),
                           theta_trt,
                           theta_ctr) {

  copula <- match.arg(copula)

  # Start tibble by generating stratum and enrollment times
  x <- tibble(
    stratum = sample(x = stratum$stratum,
                     size = n,
                     replace = TRUE,
                     prob = stratum$p)
  ) %>%
    mutate(enroll_time = rpw_enroll(n, enroll_rate)) %>%
    group_by(stratum) %>%
    # Assign treatment
    mutate(treatment = randomize_by_fixed_block(n = n(), block = block)) %>%
    # Generate time to failure and time to dropout
    group_by(stratum, treatment)

  unique_stratum <- unique(x$stratum)
  unique_treatment <- unique(x$treatment)

  x$fail_time <- 0
  x$dropout_time <- 0

  for (sr in unique_stratum) {
    # sr <- "All"

    for (tr in unique_treatment) {
      # tr <- "experimental"

      # [Among entire sample] The row index for patients in current stratum and current treatment group.
      indx <- x$stratum == sr & x$treatment == tr
      n_srtr <- sum(indx)

      fail_rate_srtr <- fail_rate[fail_rate$stratum == sr & fail_rate$treatment == tr, c('duration', 'rate'), drop = FALSE]
      drop_rate_srtr <- dropout_rate[dropout_rate$stratum == sr & dropout_rate$treatment == tr, c('duration', 'rate'), drop = FALSE]

      # Disclaimer:
      # this is a quick fixed for allowing arm-specific correlation across all stratum.
      # treatment must be either "control" or "experimental"
      # the more ideal way to do this is to add theta into failure_rate and drop_rate tibble
      # it could require input structure change that over-complicates the problem at this moment.

      # YZ: alternative is save the theta input as
      # theta = tibble(treatment = c("experimental", "control"), theta = c(0, 0))
      theta <- ifelse(tr == "control", theta_ctr,
                      ifelse(tr == "experimental", theta_trt, NA))

      if(copula == 'gumbel'){
        marg <- mvdc(rotCopula(gumbelCopula(param = theta,
                                            dim = 2,
                                            use.indepC = "FALSE"),
                               flip = rep(TRUE, 2)), #survival copula
                     margins = c('pwexp', 'pwexp'),
                     paramMargins = list(list(rate = fail_rate_srtr$rate, duration = fail_rate_srtr$duration),
                                         list(rate = drop_rate_srtr$rate, duration = drop_rate_srtr$duration)))
      }

      if(copula == 'clayton'){
        marg <- mvdc(rotCopula(claytonCopula(param = theta,
                                             dim = 2,
                                             use.indepC = "FALSE"),
                               flip = rep(TRUE, 2)), #survival copula
                     margins = c('pwexp', 'pwexp'),
                     paramMargins = list(list(rate = fail_rate_srtr$rate, duration = fail_rate_srtr$duration),
                                         list(rate = drop_rate_srtr$rate, duration = drop_rate_srtr$duration)))
      }

      time_cp <- copula::rMvdc(n_srtr, mvdc = marg)
      x$fail_time[indx] <- time_cp[,1]
      x$dropout_time[indx] <- time_cp[,2]

    }
  }
  # Set calendar time-to-event and failure indicator
  ans <- x %>%
    mutate(
      cte = pmin(dropout_time, fail_time) + enroll_time,
      fail = (fail_time <= dropout_time) * 1
    )
  ans
}
