library(gsDesign2)
library(simtrial)

######### Trial Design #########

# median PFS in the control
med.pfs.0 <- 10

# marginal dropout risk in both arms
marg.dpo.risk <- 0.15

# Constant HR after delay time
HR <- 0.6

# Delayed window (i.e. no treatment effect prior to xx month)
delay_time <- 4

# Median PFS and  Annual dropout rate in the experimental group (month)
med.pfs.1 <- qpwexp(p=0.5, duration = c(delay_time,100), rate=c(log(2)/med.pfs.0, HR*log(2)/med.pfs.0))

# Minimum follow up
min_fu <- max(med.pfs.0, med.pfs.1)

enroll_duration <- c(2, 2, 8)
enroll_relative_speed <- c(1, 2, 4)
q <- 700 / sum(enroll_duration*enroll_relative_speed)
enroll_rate <- tibble(duration = enroll_duration,
                      rate = q*enroll_relative_speed,
                      stratum = "All")

# planned trial duration in months
cut_date <- sum(enroll_duration) + min_fu

fail_rate <- tibble(
  stratum = "All",
  duration = c(4, cut_date-delay_time),
  fail_rate = log(2)/c(med.pfs.0,med.pfs.0),
  hr = c(1, HR),
  dropout_rate = rep(-log(1-marg.dpo.risk)/cut_date, 2)
)


total_duration <- 27
ratio <- 1
event <- 450
analysis_time <- 27
study_duration <- 27


#########  Average hazard  #########

total_infomation <- ahr(enroll_rate = enroll_rate,
                        fail_rate = fail_rate,
                        total_duration = total_duration,
                        ratio=ratio)
total_infomation
ahr <- total_infomation$ahr

pw_information <- pw_info(enroll_rate = enroll_rate,
                          fail_rate = fail_rate,
                          total_duration = total_duration,
                          ratio=ratio
)
pw_information
sum(pw_information$info)

######### Expected number of events #########
q_e <- ratio / (1 + ratio)
q_c <- 1 - q_e
ans <- NULL

strata <- unique(enroll_rate$stratum)
s <- strata
event <- NULL
# subset to stratum
enroll <- enroll_rate %>% filter(stratum == s)
fail <- fail_rate %>% filter(stratum == s)
# update enrollment rates
enroll_c <- enroll %>% mutate(rate = rate * q_c)
enroll_e <- enroll %>% mutate(rate = rate * q_e)
# update failure rates
fail_c <- fail
fail_e <- fail %>% mutate(fail_rate = fail_rate * hr)
# compute expected number of events
event_c <- expected_event(
  enroll_rate = enroll_c,
  fail_rate = fail_c,
  total_duration = total_duration,
  simple = FALSE
)
event_e <- expected_event(
  enroll_rate = enroll_e,
  fail_rate = fail_e,
  total_duration = total_duration,
  simple = FALSE
)

# Combine control and experimental; by period recompute HR, events, information
event <- rbind(
  event_c %>% mutate(treatment = "control"),
  event_e %>% mutate(treatment = "experimental")
) %>%
  arrange(t, treatment) %>%
  ungroup() %>%
  # recompute HR, events, info by period
  group_by(t) %>%
  summarize(
    stratum = s,
    info = (sum(1 / event))^(-1),
    event = sum(event),
    hr = last(fail_rate) / first(fail_rate)
  ) %>%
  rbind(event)

event_c
event_e
event

c(1/(1/event_c$event[1] + 1/event_e$event[1]),
  1/(1/event_c$event[2] + 1/event_e$event[2]),
  1/(1/event_c$event[1] + 1/event_e$event[1])+1/(1/event_c$event[2] + 1/event_e$event[2]))


