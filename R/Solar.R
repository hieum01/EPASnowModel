#' Solar radiation
#'
#' Estimate solar radiation
#'
#' @param lat Latitude: degree or rad
#' @param Jday Julian date or day of the year [day]
#' @param Tx Maximum temperature in degree Celcius
#' @param Tn Minimum temperature in degree Celcius
#' @param albedo Surface albedo (default=0.2)
#' @param forest Forest ratio
#' @param slope Slope
#' @param aspect Ground aspect [rad from north]
#' @param units Radiation unit: Wm2 or kJm2d (default=kJm2d)
#' @param latUnits Latitude unit (default=unknown)
#' @param printWarn option for printing warning (default=TRUE)
#'
#' @return Estimated solar radiation
#' @export
#'
#' @examples
Solar <-  function (lat, Jday, Tx, Tn, albedo=0.2, forest=0, slope=0, aspect = 0, units="kJm2d", latUnits = "unknown", printWarn=TRUE) {

  if ((abs(lat) > pi/2 & latUnits == "unknown") | latUnits == "degrees" ){
    if (printWarn==TRUE)
      lat <- lat*pi/180
  } else if (latUnits == "unknown"){
    if (printWarn==TRUE) warning("In Solar(): Input latitude units are not specified and assumed to be radians")
  }

  if (units == "kJm2d") convert <- 1 else convert <- 86.4  # can convert to W/m2
  return( signif((1 - albedo) * (1 - forest) * transmissivity(Tx, Tn) *
                   PotentialSolar(lat, Jday) * slopefactor(lat, Jday, slope, aspect) / convert , 5 ))
}
