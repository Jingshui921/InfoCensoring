library(survival)
library(simtrial)
library(copula)
library(survminer)
library(tibble)
library(dplyr)
library(tidyverse)

library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(tictoc)
library(mvtnorm)


remove(list = ls())

sapply(list.files(path = '~/R', full.names = TRUE), source)

# Parellism Computing
detectCores()
cl <- makeCluster(12)
registerDoParallel(cl)

################## Fixed Parameters #####################
# Iteration = the number of simulated dataset you generate
k <- 10000
# k <- 100000 # ten times for estimating empirical type I error.

# Random Seed
set.seed(20230605)
seed.sample <- c(sample(1:20000000,size=k))  #3259734 16331171

# Non-stratified
stratum <- tibble(stratum = "All", p = 1)

# Treatment block
# For every 4 patients enrolled, first 2 goes to control and second 2 goes to treatment.
block <- c(rep("control", 2), rep("experimental", 2))

# Constant HR after delay time
HR <- 0.6

# Flemming-Harrington type weights
# unweighted (0,0)
# early      (1,0)
# middle     (1,1)
# late       (0,1)
rho_gamma <- tibble(rho = c(0, 1, 1, 0),
                    gamma = c(0, 0, 1, 1))

# Total planned sample size and target event number
n <- 700
n_ev <- 450

# Delayed enrollment
# Enrollment rate increase over time.
enroll_duration <- c(2, 2, 8)
enroll_relative_speed <- c(1, 2, 4)
q <- n / sum(enroll_duration*enroll_relative_speed)
enroll_rate <- tibble(stratum = "All",
                      duration = enroll_duration, rate = q*enroll_relative_speed)


################## Varying Parameters #####################

# Median PFS in the control group (month)
# (med.pfs.0 <- 4)
# (med.pfs.0 <- 20)
med.pfs.0 <- 10

# Delayed duration in month (i.e. no treatment effect prior to xx month)
# delay_time <- 6
delay_time <- ifelse(med.pfs.0 == 4, 2, ifelse(med.pfs.0 == 10, 4, 6))


# Median PFS in the experimental group (month)
med.pfs.1 <- qpwexp(p=0.5, duration = c(delay_time,100), rate=c(log(2)/med.pfs.0, HR*log(2)/med.pfs.0))


# Derive average median PFS that defines dropout window
(med.pfs.avg <- (med.pfs.0+med.pfs.1)/2)

# Minimum follow up
min_fu <- max(med.pfs.0, med.pfs.1)

# planned trial duration in months
cut_date <- sum(enroll_rate$duration) + min_fu

# Failure Time
# NPH - Delayed Treatment Effect
fail_rate <- tibble(
  stratum = "All",
  period = rep(1:2, 2),
  treatment = c(rep("control", 2), rep("experimental", 2)),
  duration = rep(c(delay_time, 999), 2),
  rate = log(2) / c(med.pfs.0, med.pfs.0, med.pfs.0, med.pfs.0/HR)
)

# Marginal dropout(dpo) risk in the control group
marg.dpo.risk.0 <- 0.25

# Fixed marginal dropout risk ratio or annual dropout rate ratio
DR <- 0.8
marg.dpo.risk.1 <- DR*marg.dpo.risk.0

# Dropout Time

# # # # Scheme 2 fixed marginal ANNUAL dropout risk over the allowed time window # # # #

case <- 2

pattern <-
  switch(case, 'early_annual', 'late_annual', 'overall_annual')
dropout_rate <-
  switch(case,
         tibble(
           stratum = "All",
           period = rep(1:2, 2),
           treatment = c(rep("control", 2),
                         rep("experimental", 2)),
           duration = rep(c(med.pfs.avg, 999),2),
           rate = c(-log(1-marg.dpo.risk.0)/12, 0.00000001,
                    -log(1-marg.dpo.risk.1)/12, 0.00000001)),
         tibble(
           stratum = "All",
           period = rep(1:2, 2),
           treatment = c(rep("control", 2),
                         rep("experimental", 2)),
           duration = rep(c(med.pfs.avg, 999),2),
           rate = c(0.00000001, -log(1-marg.dpo.risk.0)/12,
                    0.00000001, -log(1-marg.dpo.risk.1)/12)),
         tibble(
           stratum = "All",
           period = rep(1, 2),
           treatment = c(rep("control", 1),
                         rep("experimental", 1)),
           duration = rep(c(999),2),
           rate = c(-log(1-marg.dpo.risk.0)/12,
                    -log(1-marg.dpo.risk.1)/12))
         )


# choice of copula
# copula <- 'gumbel'
copula <- 'clayton'

# concordance parameter
## Both negative:
theta_trt <- -0.5
theta_ctr <- -0.5

## Both independent: 0
theta_trt <- 0
theta_ctr <- 0

## Both positive:
theta_trt <- 1
theta_ctr <- 1

# theta_trt <- 5
# theta_ctr <- 5

## Trt Negative, Ctr Positive
theta_trt <- -0.5
theta_ctr <- 1

# theta_trt <- -0.5
# theta_ctr <- 5


## Trt Independent, Ctr Positive
theta_trt <- 0
theta_ctr <- 1

# theta_trt <- 0
# theta_ctr <- 5

paste(copula,'_med',med.pfs.0, '_delay', delay_time,
      '_ev',n_ev,
      '_',pattern,
      '_theta', theta_trt,'_',theta_ctr,
      '.csv', sep = '')

###### Simulation ######
# Result Shell
df <- NULL
tic()
df <- foreach (i = 1:k,
               .combine=rbind,
               .packages = c('simtrial','dplyr', 'tidyverse', 'tibble', 'survival', 'mvtnorm', 'survminer', 'copula'),
               .errorhandling = "remove") %dopar% {

                 setwd("~/Jingyi/2023intern")
                 devtools::load_all()

                 seed <- seed.sample[i]
                 set.seed(seed)
                 # reproducibility notes:
                 # when the seed is fixed,
                 # same enroll rate leads to same enroll time,
                 # and failure rate to fail time,
                 # and dropout rate to dropout time.

                 dat <- sim_pw_surv_cp2(n = n,
                                       stratum = stratum,
                                       block = block,
                                       enroll_rate = enroll_rate,
                                       fail_rate = fail_rate,
                                       dropout_rate = dropout_rate,
                                       copula = copula,
                                       theta_trt = theta_trt,
                                       theta_ctr = theta_ctr)
                 # options(digits = 3)
                 # View(dat)

                 te_date <- dat %>% get_cut_date_by_event(event=n_ev)

                 dat <- dat %>%
                   mutate(tte = ifelse(te_date < enroll_time, NA, pmin(cte, te_date) - enroll_time),
                          event = fail * (cte <= te_date)) %>%
                   mutate(censored = 1-event,
                          censor_type = ifelse(te_date < enroll_time ,'Not Enrolled',
                                               ifelse(te_date < cte , 'Administrative',
                                                      ifelse(cte <= te_date & fail == 0, 'Dropout',
                                                             ifelse(cte <= te_date & fail == 1, 'Event', NA)))))
                 # # head(dat)
                 # dev.new()
                 # plot(survfit(Surv(tte, censored) ~ treatment, data=dat))
                 # plot(survfit(Surv(tte, event) ~ treatment, data=dat))
                 # # sanity check
                 # sum(is.na(dat$tte))
                 # table(dat$censor_type, dat$treatment)
                 # table(dat$treatment)
                 # te_date - cut_date

                 # Aggregate KM plot for both dropout and administrative censoring
                 fit.1 <- survfit(Surv(tte, censored) ~ treatment, data=dat %>% filter(treatment == 'experimental'))
                 fit.0 <- survfit(Surv(tte, censored) ~ treatment, data=dat %>% filter(treatment == 'control'))

                 cut_points <- seq(0, cut_date+4, 0.5)
                 prob.1 <- as.data.frame(do.call(cbind, summary(fit.1, times = cut_points)[c('time','surv')]))
                 prob.0 <- as.data.frame(do.call(cbind, summary(fit.0, times = cut_points)[c('time','surv')]))

                 rkm.1 <- left_join(data.frame(time = cut_points),
                                    prob.1, by='time') %>% fill(surv,.direction = 'down') %>%
                   pivot_wider(names_from = time, values_from = surv, names_prefix = 'rkm_trt')
                 rkm.0 <- left_join(data.frame(time = cut_points),
                                    prob.0, by='time') %>% fill(surv,.direction = 'down') %>%
                   pivot_wider(names_from = time, values_from = surv, names_prefix = 'rkm_ctr')

                 # Aggregate KM plot for PFS
                 fit.1 <- survfit(Surv(tte, event) ~ treatment, data=dat %>% filter(treatment == 'experimental'))
                 fit.0 <- survfit(Surv(tte, event) ~ treatment, data=dat %>% filter(treatment == 'control'))

                 cut_points <- seq(0, cut_date+4, 0.5)
                 prob.1 <- as.data.frame(do.call(cbind, summary(fit.1, times = cut_points)[c('time','surv')]))
                 prob.0 <- as.data.frame(do.call(cbind, summary(fit.0, times = cut_points)[c('time','surv')]))

                 km.1 <- left_join(data.frame(time = cut_points),
                                   prob.1, by='time') %>% fill(surv,.direction = 'down') %>%
                   pivot_wider(names_from = time, values_from = surv, names_prefix = 'km_trt')
                 km.0 <- left_join(data.frame(time = cut_points),
                                   prob.0, by='time') %>% fill(surv,.direction = 'down') %>%
                   pivot_wider(names_from = time, values_from = surv, names_prefix = 'km_ctr')

                 # Median PFS
                 df.med <- as.data.frame(do.call(cbind,
                                                 c(survminer::surv_median(fit.1)[2],
                                                   survminer::surv_median(fit.0)[2])))
                 colnames(df.med) <-  c('med_trt', #'med_trt_low','med_trt_up',
                                        'med_ctr'  #,'med_ctr_low','med_ctr_up'
                 )

                 # Overall censoring by types
                 ov <- data.frame(ad_ctr = sum(dat$censor_type == 'Administrative' & dat$treatment == 'control'),
                                  do_ctr = sum(dat$censor_type == 'Dropout' & dat$treatment == 'control'),
                                  ev_ctr = sum(dat$censor_type == 'Event' & dat$treatment == 'control'),
                                  ad_trt = sum(dat$censor_type == 'Administrative' & dat$treatment == 'experimental'),
                                  do_trt = sum(dat$censor_type == 'Dropout' & dat$treatment == 'experimental'),
                                  ev_trt = sum(dat$censor_type == 'Event' & dat$treatment == 'experimental'))

                 # Generalized log-rank tests
                 fit.wlr <- doAnalysis(dat, rho_gamma=rho_gamma, n_stratum=nrow(stratum))
                 n.event <- fit.wlr[1,'event']
                 ln.hr <- fit.wlr[1,'ln_hr']
                 fit.wlr <- fit.wlr[,3:ncol(fit.wlr)]
                 z.stat <- fit.wlr$z

                 # Compute p-value for MaxCombo
                 p.max <- 1 - pmvnorm(
                   lower = rep(min(fit.wlr$z), nrow(fit.wlr)),
                   corr = data.matrix(select(fit.wlr, -c(rho, gamma, z))),
                   algorithm = GenzBretz(maxpts = 50000, abseps = 0.00001)
                 )[1]


                 res <- cbind(te_date = as.numeric(te_date),
                              rkm.1, rkm.0, km.1, km.0,
                              ov,
                              n.event, ln.hr,
                              z.00 = z.stat[1],
                              z.10 = z.stat[2],
                              z.11 = z.stat[3],
                              z.01 = z.stat[4],
                              p.max, df.med)

                 return(c(i,res))

               }
toc()
# 20 : 5254.473 sec elapsed


paste(copula,'_med',med.pfs.0, '_delay', delay_time,
      '_ev',n_ev,
      '_',pattern,
      '_theta', theta_trt,'_',theta_ctr,
      '.csv', sep = '')

###### Saving Output ######
# Transform the output into a dataframe
colnames(df)[1] <- c('Iteration')
res <- as.data.frame(matrix(unlist(df), nrow = dim(df)[1]))
colnames(res) <- colnames(df)

head(res)

# derive the p-values of WLR tests
res$p.00 <- pnorm(res$z.00)
res$p.10 <- pnorm(res$z.10)
res$p.11 <- pnorm(res$z.11)
res$p.01 <- pnorm(res$z.01)

# Adding fixed parameters
res$copula <- copula
res$pattern <- pattern
res$theta_trt <- theta_trt
res$theta_ctr <- theta_ctr
res$n <- n
res$marg.dpo.risk.0 <- marg.dpo.risk.0
res$DR <- DR
res$delay_time <- delay_time
res$med.pfs.0 <- med.pfs.0
res$HR <- HR
res$cut_date <- cut_date #276 cols

if(sum(complete.cases(res)) != k){
  res[!complete.cases(res),]
}

if (nrow(res) != k) {
  cbind(iter = which(!(1:k %in% res$Iteration)),
        seed = seed.sample[which(!(1:k %in% res$Iteration))])
}

getwd()
(title <- paste(copula,'_med',med.pfs.0, '_delay', delay_time,
                '_ev',n_ev,
                '_',pattern,
                '_theta', theta_trt,'_',theta_ctr,
                '.csv', sep = ''))

# save raw simulated dataset
write.table(res,
            file = file.path('vignette','fixtures',
                             'Dependent dropout via Copula',
                             'NPH',title),
            sep = ',',
            append = FALSE,
            col.names =  !FALSE,
            row.names = FALSE,
            quote=FALSE)

# save rkm

rkm <- matrix(c(pattern, theta_trt, theta_ctr,
                med.pfs.0 = med.pfs.0, delay_time,
                marg.dpo.risk.0 = marg.dpo.risk.0, DR,
                cut_date = cut_date,
                res %>% select(starts_with('rkm_')) %>% apply(MARGIN=2, mean, na.rm= TRUE)),
              nrow = 1)
colnames(rkm) <- c('pattern', 'theta_trt', 'theta_ctr', 'med.pfs.0', 'delay_time',
                   'marg.dpo.risk.0', 'DR', 'cut_date',
                   res %>% select(starts_with('rkm_')) %>% colnames())
rkm <- as.data.frame(rkm)
rkm$copula <- copula
(title.rkm <- paste('Copula_RKM_med',med.pfs.0, '_delay', delay_time,'.csv', sep = ''))

write.table(rkm,
            file = file.path('vignette','fixtures',
                             'Dependent dropout via Copula',
                             'NPH',title.rkm),
            sep = ',',
            append = TRUE,
            col.names = !file.exists(file.path('vignette','fixtures',
                                               'Dependent dropout via Copula',
                                               'NPH',title.rkm)),
            row.names = FALSE,
            quote=FALSE)

# save km
km <- matrix(c(pattern, theta_trt, theta_ctr,
               med.pfs.0 = med.pfs.0, delay_time,
               marg.dpo.risk.0 = marg.dpo.risk.0, DR,
               cut_date = cut_date,
               res %>% select(starts_with('km_')) %>% apply(MARGIN=2, mean, na.rm= TRUE)),
             nrow = 1)
colnames(km) <- c('pattern', 'theta_trt', 'theta_ctr', 'med.pfs.0', 'delay_time',
                  'marg.dpo.risk.0', 'DR','cut_date',
                  res %>% select(starts_with('km_')) %>% colnames())
km <- as.data.frame(km)
km$copula <- copula
(title.km <- paste('Copula_KM_med',med.pfs.0, '_delay', delay_time,'.csv', sep = ''))

write.table(km,
            file = file.path('vignette','fixtures',
                             'Dependent dropout via Copula',
                             'NPH',title.km),
            sep = ',',
            append = TRUE,
            col.names = !file.exists(file.path('vignette','fixtures',
                                               'Dependent dropout via Copula',
                                               'NPH',title.km)),
            row.names = FALSE,
            quote=FALSE)


# save ov
ov <- matrix(c(copula, pattern, theta_trt, theta_ctr, n,
               med.pfs.0 = med.pfs.0, delay_time,
               marg.dpo.risk.0 = marg.dpo.risk.0,
               DR,
               cut_date = round(cut_date,4),
               res %>% select(starts_with(c('ad_','do_','ev_'))) %>%
                 apply(MARGIN=2, mean, na.rm= TRUE)),
             nrow = 1)
colnames(ov) <- c('copula', 'pattern','theta_trt', 'theta_ctr', 'n', 'med.pfs.0', 'delay_time',
                  'marg.dpo.risk.0', 'DR', 'cut_date',
                  res %>% select(starts_with(c('ad_','do_','ev_'))) %>% colnames())

write.table(ov,
            file = file.path('vignette','fixtures',
                             'Dependent dropout via Copula',
                             'NPH','Copula overall counts.csv'),
            sep = ',',
            append = TRUE,
            col.names =  !file.exists(file.path('vignette','fixtures',
                                                'Dependent dropout via Copula',
                                                'NPH','Copula overall counts.csv')),
            row.names = FALSE,
            quote=FALSE)

# save est and power
power <- data.frame(n, n_ev, copula, pattern, theta_trt, theta_ctr,
                    marg.dpo.risk.0, DR,
                    delay_time, med.pfs.0, med.pfs.1, HR, cut_date,
                    te_date = mean(res$te_date,, na.rm= TRUE),
                    HR.cox = mean(exp(res$ln.hr), na.rm= TRUE),
                    HR.cox.se = sd(exp(res$ln.hr), na.rm= TRUE),
                    med.trt    = mean(res$med_trt, na.rm= TRUE),
                    med.trt.se = sd(res$med_trt, na.rm= TRUE),
                    med.ctr    = mean(res$med_ctr, na.rm= TRUE),
                    med.ctr.se = sd(res$med_ctr, na.rm= TRUE),
                    power.max = mean(res$p.max<= 0.025, na.rm = TRUE),
                    power.00  = mean(res$p.00 <= 0.025, na.rm = TRUE),
                    power.10  = mean(res$p.10 <= 0.025, na.rm = TRUE),
                    power.11  = mean(res$p.11 <= 0.025, na.rm = TRUE),
                    power.01  = mean(res$p.01 <= 0.025, na.rm = TRUE))

write.table(power,
            file = file.path('vignette','fixtures',
                             'Dependent dropout via Copula',
                             'NPH','Copula Power.csv'),
            sep = ',',
            append = TRUE,
            col.names =  !file.exists(file.path('vignette','fixtures',
                                                'Dependent dropout via Copula',
                                                'NPH','Copula Power.csv')),
            row.names = FALSE,
            quote=FALSE)

# mcse.power
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8086
max(sqrt(power$power.max*(1-power$power.max)/nrow(res)),
    sqrt(power$power.00 *(1-power$power.00 )/nrow(res)),
    sqrt(power$power.10 *(1-power$power.10 )/nrow(res)),
    sqrt(power$power.11 *(1-power$power.11 )/nrow(res)),
    sqrt(power$power.01 *(1-power$power.01 )/nrow(res)))

# take a peek at the rkm and km
rkm <- data.frame(rkm) %>%
  pivot_longer(cols = starts_with('rkm_'), names_to = 'grp_time', values_to = 'surv',
               names_prefix = 'rkm_') %>%
  separate_wider_position(grp_time, widths = c(trt=3, time=4),too_few = "align_start") %>%
  mutate(time = as.numeric(time),
         surv = as.numeric(surv)) %>%
  mutate(trt = factor(trt, levels = c("ctr", "trt"),
                      labels = c('Control','Experimental'))
  )

ggplot(rkm, aes(x = time, y = surv, group = trt, color = trt, linetype = trt)) +
  geom_line(linewidth = 0.8, show.legend = TRUE) +
  geom_vline(xintercept = delay_time, linetype = "twodash", colour = "gray") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks=seq(0,30,3)) +
  xlab('Months') +
  ylab('') +
  theme_bw() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_colour_manual(values=c("#EFC000FF","#0073C2FF"))

km <- data.frame(km) %>%
  pivot_longer(cols = starts_with('km_'), names_to = 'grp_time', values_to = 'surv',
               names_prefix = 'km_') %>%
  separate_wider_position(grp_time, widths = c(trt=3, time=4),too_few = "align_start") %>%
  mutate(time = as.numeric(time),
         surv = as.numeric(surv)) %>%
  mutate(trt = factor(trt, levels = c("ctr", "trt"),
                      labels = c('Control','Experimental'))
  )

ggplot(km, aes(x = time, y = surv, group = trt, color = trt, linetype = trt)) +
  geom_line(linewidth = 0.8, show.legend = TRUE) +
  geom_vline(xintercept = delay_time, linetype = "twodash", colour = "gray") +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(breaks=seq(0,30,3)) +
  xlab('Months') +
  ylab('') +
  theme_bw() +
  theme(legend.position="top", legend.title = element_blank()) +
  scale_colour_manual(values=c("#EFC000FF","#0073C2FF"))
