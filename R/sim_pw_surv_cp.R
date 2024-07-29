#' Title
#' @importFrom tibble tibble
#' @importFrom simtrial randomize_by_fixed_block rpw_enroll rpwexpinvRcpp
#' @importFrom copula rotCopula gumbelCopula claytonCopula
#' @import dplyr
#'
#' @param n
#' @param stratum
#' @param block
#' @param enroll_rate
#' @param fail_rate
#' @param dropout_rate
#' @param copula
#' @param theta_trt
#' @param theta_ctr
#'
#' @return
#' @export
#'
#' @examples
sim_pw_surv_cp <- function(n,
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

    x$fail_time_cp <- 0
    x$drop_time_cp <- 0
    x$ttfe_time_cp  <- 0
    x$fail_flag_cp <- 0

    for (sr in unique_stratum) {

      for (tr in unique_treatment) {

        # [Among entire sample] The row index for patients in current stratum and current treatment group.
        indx <- x$stratum == sr & x$treatment == tr
        n_srtr <- sum(indx)

        # Simulate failure time for indx
        x$fail_time[indx] <- rpwexpinvRcpp(
          n = n_srtr,
          fail_rate = fail_rate[fail_rate$stratum == sr & fail_rate$treatment == tr, , drop = FALSE]
        )

        # Simulate dropout time for indx
        x$dropout_time[indx] <- rpwexpinvRcpp(
          n = n_srtr,
          fail_rate = dropout_rate[dropout_rate$stratum == sr & dropout_rate$treatment == tr, , drop = FALSE]
        )

        ################################# Beginning of rgrumbel_pwexp #############################################
        # After QC, the code that generates copula with pw exp margins can be packaged into a function
        # Input: fail_rate_srtr, drop_rate_srtr
        # Output: ans_fail, ans_drop, ans_ttfe, ans_fail_flag

        # [Notice!] I have not tested the code below for single rate fail or dropout.
        # do not enter duration of 0 which makes no sense.

        # Simulate failure and dropout time simultaneously from a copula model

        fail_rate_srtr <- fail_rate[fail_rate$stratum == sr & fail_rate$treatment == tr, c('duration', 'rate'), drop = FALSE]
        drop_rate_srtr <- dropout_rate[dropout_rate$stratum == sr & dropout_rate$treatment == tr, c('duration', 'rate'), drop = FALSE]

        fail_rate_srtr <- collapse(fail_rate_srtr)
        drop_rate_srtr <- collapse(drop_rate_srtr)

        fail_end_time <- cumsum(fail_rate_srtr$duration)
        drop_end_time <- cumsum(drop_rate_srtr$duration)

        end_max <- max(c(fail_end_time, drop_end_time))
        fail_end_time[length(fail_end_time)] <- end_max
        drop_end_time[length(drop_end_time)] <- end_max

        joint_end_time <- c(sort(unique(c(fail_end_time, drop_end_time))))
        joint_start_time <- c(0, joint_end_time[-length(joint_end_time)])
        joint_duration <- joint_end_time - joint_start_time

        fail_rep <- diff(c(0, sapply(1:length(fail_end_time), function(x) sum(fail_end_time[x] >= joint_end_time))))
        drop_rep <- diff(c(0, sapply(1:length(drop_end_time), function(x) sum(drop_end_time[x] >= joint_end_time))))

        n_joint <- length(joint_end_time) # total number of windows

        fail_add <- n_joint-sum(fail_rep) # fail_add might be redundant, and always equal to 0
        drop_add <- n_joint-sum(drop_rep) # drop_add might be redundant, and always equal to 0

        fail_rep[length(fail_rep)] <- fail_rep[length(fail_rep)] + fail_add
        drop_rep[length(drop_rep)] <- drop_rep[length(drop_rep)] + drop_add

        joint_fail_rate <- rep(fail_rate_srtr$rate, times = fail_rep)
        joint_drop_rate <- rep(drop_rate_srtr$rate, times = drop_rep)

        # (joint_rate <- data.frame(
        #   window = 1:n_joint,
        #   joint_start_time,
        #   joint_end_time,
        #   joint_duration,
        #   joint_fail_rate,
        #   joint_drop_rate))

        # [Among current sr tr] # how many tte to generate
        indx_srtr <- rep(TRUE, n_srtr)
        ans_fail <- rep(0, n_srtr)
        ans_drop <- rep(0, n_srtr)
        ans_ttfe <- rep(0, n_srtr)
        ans_fail_flag <- rep(0, n_srtr)

        for (i in 1:n_joint) { # i for interval time window.

          # i <- 1
          # (i <- i+1)
          n_event_left <- sum(indx_srtr) # More intuitively, it should be #patients left, n_pat_left

          if (n_event_left == 0) {
            break
          }

          # Set failure time (or dropout time) to Inf for interval i if failure rate = 0
          # And the joint copula of dropout time (or failure time) reduce to marginal distribution of dropout.
          if (joint_fail_rate[i] == 0 | joint_drop_rate[i] == 0){

            if (joint_fail_rate[i] == 0 & joint_drop_rate[i] != 0) {
              ans_fail[indx_srtr] <- joint_start_time[i] + rep(Inf, n_event_left)
              ans_drop[indx_srtr] <- joint_start_time[i] + stats::rexp(n_event_left, joint_drop_rate[i])
            }

            if (joint_fail_rate[i] != 0 & joint_drop_rate[i] == 0) {
              ans_fail[indx_srtr] <- joint_start_time[i] + stats::rexp(n_event_left, joint_fail_rate[i])
              ans_drop[indx_srtr] <- joint_start_time[i] + rep(Inf, n_event_left)
            }

            if (joint_fail_rate[i] == 0 & joint_drop_rate[i] == 0) {
              ans_fail[indx_srtr] <- joint_start_time[i] + rep(Inf, n_event_left)
              ans_drop[indx_srtr] <- joint_start_time[i] + rep(Inf, n_event_left)
            }

          } else { # joint_fail_rate[i] != 0 & joint_drop_rate[i] != 0

            # Disclaimer:
            # this is a quick fixed for allowing arm-specific correlation across all stratum.
            # treatment must be either "control" or "experimental"
            # the more ideal way to do this is to add theta into failure_rate and drop_rate
            # it could require input structure change that over-complicates the problem at this moment.

            theta <- ifelse(tr == "control", theta_ctr,
                            ifelse(tr == "experimental", theta_trt, NA))

            if(copula == 'gumbel'){
              marg <- mvdc(rotCopula(gumbelCopula(param = theta,
                                                  dim = 2,
                                                  use.indepC = "FALSE"),
                                     flip = rep(TRUE, 2)), #survival copula
                           margins = c('exp', 'exp'),
                           paramMargins = list(list(rate = joint_fail_rate[i]),
                                               list(rate = joint_drop_rate[i])))
            }

            if(copula == 'clayton'){
              marg <- mvdc(rotCopula(claytonCopula(param = theta,
                                                   dim = 2,
                                                   use.indepC = "FALSE"),
                                     flip = rep(TRUE, 2)), #survival copula
                           margins = c('exp', 'exp'),
                           paramMargins = list(list(rate = joint_fail_rate[i]),
                                               list(rate = joint_drop_rate[i])))
            }


            time_cp <-  joint_start_time[i] + copula::rMvdc(n_event_left, mvdc = marg)
            # plot(time_cp)

            ans_fail[indx_srtr] <- time_cp[,1]
            ans_drop[indx_srtr] <- time_cp[,2]
            ans_ttfe[indx_srtr] <- pmin(time_cp[,1], time_cp[,2])
            ans_fail_flag[indx_srtr] <- ifelse(time_cp[,1] <= time_cp[,2],1,0)
            # "Fail or Dropout" 1 for event, 0 for dropout

            # cbind(time_cp,ans_ttfe[indx_srtr],ans_fail_flag[indx_srtr])
            # cbind(ans_fail, ans_drop, ans_ttfe, ans_fail_flag)

            if (i < n_joint) { # if not the last interval
              indx_srtr <- (ans_ttfe > joint_end_time[i])
              # indx_srtr = TRUE for those whose both failure and dropout time fall outside the current interval
            }
          }
        } # END OF for (i in 1:n_joint)
        # data.frame(ans_fail, ans_drop, ans_ttfe, ans_fail_flag)
        ################################## End of of rgrumbel_pwexp ###################################################

        x$fail_time_cp[indx] <- ans_fail
        x$drop_time_cp[indx] <- ans_drop
        x$ttfe_time_cp[indx] <- ans_ttfe
        x$fail_flag_cp[indx] <- ans_fail_flag

      } # END OF for (tr in unique_treatment)
    } # END OF for (sr in unique_stratum)

    # Set calendar time-to-event and failure indicator
    ans <- x %>%
      mutate(
        cte = pmin(dropout_time, fail_time) + enroll_time,
        fail = (fail_time <= dropout_time) * 1,
        cte_cp = ttfe_time_cp + enroll_time
      ) %>%
      rename(fail_cp = fail_flag_cp) %>%
      select(stratum, enroll_time, treatment,
             fail_time, dropout_time, cte, fail,
             fail_time_cp, drop_time_cp, cte_cp, fail_cp, ttfe_time_cp)

    ans
  }
