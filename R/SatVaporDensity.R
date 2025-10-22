#' Saturated Vapor Density
#'
#' @param T_C Temperature in degree Celsius
#'
#' @return  Saturated Vapor Density [kg/m3]
#' @export
#'
#' @examples SatVaporDensity(5)
SatVaporDensity <- function(T_C){
  #	T_C	= Temperature [C]
  VP <- 0.611 * exp((17.3*T_C)/(237.2+T_C))  # saturated water pressure in kPa
  return(round(VP/(0.462 * (T_C+273.15)), 5))
}

