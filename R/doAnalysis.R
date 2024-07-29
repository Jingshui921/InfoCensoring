#' Title
#'
#' @param d is a data.frame simulated from sim_pw_surv(). It must contain the following columns:
#' arm = c("experimental", "control"), rho_gamma = the FH-type weights,
#' @param rho_gamma
#' @param n_stratum
#'
#' @return
#' @export
#'
#' @examples
doAnalysis <- function(d, rho_gamma, n_stratum) {
  if (nrow(rho_gamma) == 1) {
    z <- tibble(z = (d %>% counting_process(arm = "experimental") %>% wlr(rho_gamma = rho_gamma))$z)
  } else {
    z <- d %>%
      counting_process(arm = "experimental") %>%
      wlr(rho_gamma = rho_gamma, return_corr = TRUE)
  }

  ans <- tibble(
    event = sum(d$event),
    ln_hr = ifelse(n_stratum > 1,
                   survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental") + survival::strata(stratum), data = d)$coefficients,
                   survival::coxph(survival::Surv(tte, event) ~ (treatment == "experimental"), data = d)$coefficients
    ) %>% as.numeric()
  )

  ans <- cbind(ans, z)
  return(ans)
}
