#' Solar Angle
#'
#' Angle of solar inclination from horizontal at solar noon [rad]
#'
#' @param lat latitdue [rad]
#' @param Jday Julian date or day of the year [day]
#'
#' @return Solar angle in radian
#' @export
#'
#' @examples
solarangle <-
  function(lat,Jday){

   # solar declination [rad]
    dec<-declination(Jday)

    return(asin(sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(0)))
  }

