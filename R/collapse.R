#' Title
#'
#' @param rate_srtr
#'
#' @return
#' @export
#'
#' @examples
#'
collapse <- function(rate_srtr){
  # if the rate is the same in two consecutive intervals, collapse into one interval
  if(nrow(rate_srtr) == 1) {
    return(rate_srtr)
  } else{
    for(j in 1:(nrow(rate_srtr)-1)){
      if(rate_srtr$rate[j]==rate_srtr$rate[j+1]){
        rate_srtr$duration[j] <- rate_srtr$duration[j]+rate_srtr$duration[j+1]
        rate_srtr$duration[j+1] <- 0
      }
    }
    rate_srtr <- rate_srtr[rate_srtr$duration>0,,drop = FALSE]
    return(rate_srtr)
  }
}
