#' Longwave
#'
#' Incoming longwave radiation based on the Sephan-Boltzman equation
#'
#' @param emissivity Terrestrial Emissivity
#' @param temp temperature of the emitting body in degree Celsius
#' @param TimeStep Time step in hour: Ex) 1 hour=1, daily=24 (default)
#'
#' @return Incoming longwave radiation [kJ m-2/TimeStep]
#' @export
#'
#' @examples
Longwave <-
  function(emissivity,temp,TimeStep=24){
    #emissivity: [-]
    #temp: temperature of the emitting body [C]
    sigma<-5.670374419/(10^8)  #  Stefan–Boltzmann constant [W/m2/K4]
    SBconstant<-sigma*3600*TimeStep/10^3

#    SBconstant<-0.00000490 #[kJ m-2 K-4 d-1]

    tempK<-temp+273.15 #[degrees K]

    return(emissivity*SBconstant*tempK^4)
  }

