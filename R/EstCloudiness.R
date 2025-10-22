#' Cloudiness
#'
#' @param Tx Daily Maximum temperature in degree Celcius
#' @param Tn Daily Minimum temperature in degree Celcius
#' @param trans transmissivity (default=NULL--> use transmissivity function)
#' @param transMin Minimum transmissivity (default=0.15)
#' @param transMax Maximum transmissivity (default=0.75)
#' @param opt Cloudiness estimation methods: Black or linear (default=linear)
#'
#' @return Cloudiness in fraction
#' @export
#'
#' @examples
EstCloudiness <- function (Tx=(-999), Tn=(-999), trans=NULL, transMin = 0.15, transMax = 0.75, opt = "linear")
{
  suppressWarnings(if (length(which(Tx< -999) > 0) && is.null(trans)){
    print("Error: Please enter either Max&Min temp or transmissivity")
  } else {
    if (is.null(trans))	trans <- transmissivity(Tx, Tn)
    if (opt=="Black") {
      cl <- (0.34 - sqrt(0.34^2 + 4*0.458*(0.803-trans)))/(-2*0.458)
      cl[which(trans > 0.803)] <- 0
    } else {
      cl <- 1 - (trans-transMin) / (transMax-transMin)
    }
    cl[which(cl > 1)] <- 1
    cl[which(cl < 0)] <- 0
    return(cl)
  } )
}

