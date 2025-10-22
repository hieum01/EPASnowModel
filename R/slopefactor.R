#' Slope factor
#'
#' Adjust solar radiation for land slope and aspect relative to the sun, 1=level ground
#'
#' @param lat Latitde [rad]
#' @param Jday Julian date or day of the year [day]
#' @param slope Slope of the ground [rad]
#' @param aspect Ground aspect [rad from north]
#'
#' @return Slope factor [0,1]
#' @export
#'
#' @examples
slopefactor <- function(lat,Jday,slope,aspect){

  SolAsp <- rep(pi, length(Jday))  # Average Solar aspect is binary - either north (0) or south (pi) for the day
  SolAsp[which(lat - declination(Jday) < 0)] <- 0   #
  SF <- cos(slope) - sin(slope)*cos(aspect-(pi-SolAsp))/tan(solarangle(lat,Jday))
  SF[which(SF < 0)] <- 0  ## Slope factors less than zero are completely shaded

  return( SF )
}

