library(gsDesign2)
library(tibble)
library(dplyr)

 
sapply(list.files(path = '~/R', full.names = TRUE), source)

fixed_nph <- NULL

for (beta in c(0.1, 0.2)) {
  for(med.pfs.0 in c(4, 10, 20)){
    #  Median PFS in the control group (month)
    for(marg.dpo.risk.0 in c(.05, .15, .25)){
      # Annual dropout rate in the control group

      # Target type 1 and type 2 error
      alpha <- 0.025

      # Constant HR after delay time
      HR <- 0.6

      # Fixed dropout rate ratio ()
      DR <- 0.8

      # Delayed duration in month (i.e. no treatment effect prior to xx month)
      # delay_time <- med.pfs.0/4
      delay_time <- ifelse(med.pfs.0 == 4, 2,
                           ifelse(med.pfs.0 == 10, 4,
                                  ifelse(med.pfs.0 == 20, 6, NA)))


      # Median PFS and  Annual dropout rate in the experimental group (month)
      # med.pfs.1 <- median_surv_pw2(rate.phase1 = log(2)/med.pfs.0,
      #                              rate.phase2 = HR*log(2)/med.pfs.0,
      #                              delay.time = delay_time)
      med.pfs.1 <- qpwexp(p=0.5, duration = c(delay_time,100), rate=c(log(2)/med.pfs.0, HR*log(2)/med.pfs.0))
      marg.dpo.risk.1 <- DR*marg.dpo.risk.0


      # Enrollment rate decrease over time.
      (enroll_rate <- tibble(stratum = "All", duration = c(2, 2, 8), rate = c(1, 2, 4)))

      # Minimum follow up
      min_fu <- max(med.pfs.0, med.pfs.1)

      # planned trial duration in months
      cut_date <- sum(enroll_rate$duration) + min_fu

      fail_rate <- tibble(
        stratum = "All",
        duration = c(delay_time, cut_date-delay_time),
        fail_rate = log(2)/c(med.pfs.0,med.pfs.0),
        hr = c(1, HR),
        dropout_rate = rep(-log(1-max(marg.dpo.risk.0, marg.dpo.risk.1))/12, 2)
        # use the max annual dropout rate for a more conservative sample size estimates
      )

      res <- gs_design_ahr(
        enroll_rate = enroll_rate,
        fail_rate = fail_rate,
        analysis_time = cut_date, # total study duration
        ratio = 1,                # randomization ratio (exp:control)
        alpha = alpha,            # type I error
        beta  = beta,             # 1-beta: targeted power
        # efficacy bound
        upper = gs_b,             # fixed boundary
        upar = -qnorm(0.025),
        # futility bound
        lower = gs_b,             # fixed boundary
        lpar = -Inf,              # no futility bound, only efficacy, one-sided design
      )

      # res %>%  summary()  %>%  as_gt()

      fixed_nph <- bind_rows(fixed_nph ,
                             res$analysis %>%
                               mutate(alpha, power = 1-beta, HR, delay_time, enroll_time= sum(enroll_rate$duration), min_fu, med.pfs.0, med.pfs.1, marg.dpo.risk.0) %>%
                               select(alpha, power, HR, delay_time, med.pfs.0, med.pfs.1, marg.dpo.risk.0, enroll_time, min_fu,
                                      time, n, event, ahr))

    }
  }
}

View(fixed_nph%>% mutate(enroll_permonth = n/enroll_time))

fixed_nph%>% mutate(enroll_permonth = n/enroll_time) %>% group_by(power) %>% gt::gt()
