#' solar declination
#'
#' @param Jday Julian date or day of the year [day]
#'
#' @return solar declination [rad]
#' @export
#'
#' @examples
declination <-
  function(Jday){
    return(0.4102*sin(pi*(Jday-80)/180))
  }

