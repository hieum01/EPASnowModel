#' PBSM with S-curve and Sublimation
#'
#' Original PBSM has been Modified by incorporating a S-shape curve method that separates rain and snow by a S-curve
#' Added ice rain in snowpack if snowpack temperature is below 0 and rain is on the snowpack.
#' Added Sublimation physics (Mass transfer based on Latent Heat Flux).
#'
#' Added two arguments (Train and Tsnow) that are minimum temperature (Train) for all precipitation and maximum temperature (Tsnow) for all snow.
#' Author: Hyung Eum, AEP
#' Date: March 2021 (Modified for Sublimation Feb 2026)
#'
#' @param Date Date (YYYY-MM-DD)
#' @param precip_mm Precipitation in mm/day
#' @param Tmax_C Maximum temperature in degree Celsius
#' @param Tmin_C Minimum temperature in degree Celsius
#' @param lat_deg Latitude in degree
#' @param slope Slope factor
#' @param aspect Ground aspect [rad from north]
#' @param tempHt Height of temperature station
#' @param windHt Height of wind station
#' @param groundAlbedo Ground albedo (default=0.25)
#' @param SurfEmissiv Surface Emissivity
#' @param windSp Wind speed in m/sec (default=2.0)
#' @param forest Forest ratio
#' @param startingSnowDepth_m Initial value of snow depth in meter
#' @param startingSnowDensity_kg_m3 Initial value of snow density in kg/m3
#' @param Train Upper temperature of S curve: precipitation is all rain if temperature >= Train
#' @param Tsnow Lower temperature of S curve: precipitation is all snow if temperature <= Tsnow
#' @param Max_SnowDensity_kg_m3 Maximum snow density in kg/m3, default=550 from Voordendag et al.(2021)
#'
#' @return Time series of snow simulation including sublimation
#' @export
SnowMelt_S_Curve_Sublim <- function(Date, precip_mm, Tmax_C, Tmin_C, lat_deg, slope=0, aspect=0, tempHt=1, windHt=2,
                             groundAlbedo=0.25, SurfEmissiv=0.95, windSp=2, forest=0, startingSnowDepth_m=0,
                             startingSnowDensity_kg_m3=550, Train=10, Tsnow=0, Max_SnowDensity_kg_m3=550){

  ## Constants :
  WaterDens <- 1000          # kg/m3
  IceDens   <- 917           # kg/m3
  lambda    <- 3.35*10^5     # latent heat of fusion (kJ/m3) - Note: standard is ~334 kJ/kg
  lambdaV   <- 2500          # (kJ/kg) latent heat of vaporization (Water)

  # [ADDED] Latent Heat of Sublimation
  # This constant is physically distinct from vaporization. It represents the energy required
  # to go directly from solid (snow) to gas (vapor).
  lambdaS   <- LatHeatFreez+lambdaV         # (kJ/kg) latent heat of sublimation (Snow)

  SnowHeatCap <- 2.1         # kJ/kg/C
  LatHeatFreez <- 333.3      # kJ/kg
  Cw        <- 4.2*10^3      # Heat Capacity of Water (kJ/m3/C)

  ##	Converted Inputs :
  Tav <- (Tmax_C+Tmin_C)/2       # degrees C
  precip_m <- precip_mm*0.001    # precip in m
  R_m <- precip_m                # (m) depth of rain

  #==========================================
  # Added by Hyung Eum
  R_m_Snow <- array(0, dim=c(length(R_m)))
  Ice_Snow_Depth <- R_m_Snow

  #------------------------------------------------------------------
  # Rain/Snow Separation (S-Curve)
  # Stefan W. Kienzle (2008). A new temperature based method to separate rain and snow.
  #----------------------------------------------------------------------
  Tt <- 2                  # Temperature(degrees C) where snow ratio of 50 % occurs
  Tr <- Train - Tsnow      # A range of temperature that intermediate-phase occurs

  R_m[which(Tav <= Tsnow)] <- 0  # ASSUMES ALL SNOW at < Tsnow

  Upper.Intermediate.Phase <- which(Tav > Tt & Tav < Train)
  Lower.Intermediate.Phase <- which(Tav <= Tt & Tav > Tsnow)

  F_S_Upp <- (Tav[Upper.Intermediate.Phase]-Tt)/(1.4*Tr)
  F_S_Low <- (Tav[Lower.Intermediate.Phase]-Tt)/(1.4*Tr)

  Ratio.Upper <- 5*F_S_Upp^3 - 6.76*F_S_Upp^2 + 3.19*F_S_Upp + 0.5
  Ratio.Upper[which(Ratio.Upper>1)] <- 1
  Ratio.Upper[which(Ratio.Upper<0)] <- 0

  Ratio.Lower <- 5*F_S_Low^3 + 6.76*F_S_Low^2 + 3.19*F_S_Low + 0.5
  Ratio.Lower[which(Ratio.Lower>1)] <- 1
  Ratio.Lower[which(Ratio.Lower<0)] <- 0

  R_m[Upper.Intermediate.Phase] <- Ratio.Upper * precip_m[Upper.Intermediate.Phase]
  R_m[Lower.Intermediate.Phase] <- Ratio.Lower * precip_m[Lower.Intermediate.Phase]

  #=============================================================
  # Snow Density Calculation
  # Lafaysse et al. (2017)
  NewSnowDensity <- 50 + (1.7 * ((Tav+15)^1.5))
  NewSnowDensity[which(NewSnowDensity < 50)] <- 50
  NewSnowDensity[which(is.na(NewSnowDensity))] <- 50

  NewSnowWatEq <- precip_m
  NewSnowWatEq[which(Tav >= Train)] <- 0
  NewSnowWatEq[Upper.Intermediate.Phase] <- (1-Ratio.Upper)*precip_m[Upper.Intermediate.Phase]
  NewSnowWatEq[Lower.Intermediate.Phase] <- (1-Ratio.Lower)*precip_m[Lower.Intermediate.Phase]

  #----------------------------------------------------------------------------------------------
  NewSnow <- NewSnowWatEq * WaterDens / NewSnowDensity  # m
  JDay <- strptime(Date, format="%Y-%m-%d")$yday+1
  lat <- lat_deg * pi / 180

  # Thermal Resistance (day/m)
  rh <- log((windHt+0.001)/0.001) * log((tempHt+0.0002)/0.0002) / (0.41*0.41*windSp*86400)
  if (length(windSp)==1) rh <- rep(rh, length(precip_mm))

  #-----------------------------------------------------------------------------------------
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  AE <- AtmosphericEmissivity(Tav, cloudiness)

  #-----------------------------------------------------
  #  New Variables :
  SnowTemp      <- rep(0, length(precip_m))
  rhos          <- SatVaporDensity(SnowTemp)    # vapor density at surface (kg/m3)
  rhoa          <- SatVaporDensity(Tmin_C)      # vapor density of atmosphere (kg/m3)

  SnowWaterEq   <- vector(length=length(precip_mm))
  TE            <- rep(SurfEmissiv, length(precip_mm))
  DCoef         <- rep(0, length(precip_mm))
  SnowDensity   <- rep(450, length(precip_mm))
  SnowDepth     <- vector(length=length(precip_mm))
  SnowMelt      <- rep(0, length(precip_mm))
  Albedo        <- rep(groundAlbedo, length(precip_mm))

  # [ADDED] Sublimation Variable
  # We need a new vector to store daily sublimation mass (m) so it can be exported.
  Sublimation_m <- rep(0, length(precip_mm))

  ##	Energy Terms
  H       <- vector(length=length(precip_mm))
  E       <- vector(length=length(precip_mm))
  S       <- vector(length=length(precip_mm))
  La      <- Longwave(AE, Tav)
  Lt      <- vector(length=length(precip_mm))
  G       <- 173
  P       <- Cw * R_m * Tmin_C
  Energy  <- vector(length=length(precip_mm))

  ##  Initial Values (Day 1)
  SnowWaterEq[1] <- startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens
  SnowDepth[1]   <- startingSnowDepth_m
  Albedo[1]      <- ifelse(NewSnow[1] > 0, 0.98-(0.98-0.50)*exp(-4*NewSnow[1]*10),
                           ifelse(startingSnowDepth_m == 0, groundAlbedo, max(groundAlbedo, 0.5+(groundAlbedo-0.85)/10)))

  S[1] <- Solar(lat=lat, Jday=JDay[1], Tx=Tmax_C[1], Tn=Tmin_C[1], albedo=Albedo[1], forest=forest, aspect=aspect, slope=slope)

  H[1] <- 1.29 * (Tav[1] - SnowTemp[1]) / rh[1]
  if(startingSnowDepth_m > 0) TE[1] <- 0.97
  Lt[1] <- Longwave(TE[1], SnowTemp[1])

  # [MODIFIED] Dynamic Latent Heat Selection (Day 1)
  # If snow is present, use Sublimation (lambdaS). If bare ground, use Vaporization (lambdaV).
  CurrentLatentHeat <- lambdaV
  if(startingSnowDepth_m > 0 || NewSnow[1] > 0) CurrentLatentHeat <- lambdaS

  # [MODIFIED] Energy Flux Calculation
  # Uses the selected Latent Heat constant
  E[1] <- CurrentLatentHeat * (rhoa[1] - rhos[1]) / rh[1]

  # [ADDED] Mass Flux Calculation (Sublimation)
  # Convert Energy Flux (kJ/m2/d) to Mass Flux (m/day)
  # Negative = Sublimation (Loss), Positive = Condensation/Deposition (Gain)
  VaporMassFlux_m <- E[1] / CurrentLatentHeat / WaterDens

  # [ADDED] Mass Limiter
  # Calculate available snow (Yesterday's SWE + New Snow)
  AvailableSWE <- SnowWaterEq[1] + NewSnowWatEq[1]

  # If sublimation (loss) is greater than available snow, cap it.
  if (VaporMassFlux_m < 0 && abs(VaporMassFlux_m) > AvailableSWE) {
    VaporMassFlux_m <- -AvailableSWE
    E[1] <- VaporMassFlux_m * WaterDens * CurrentLatentHeat # Recalculate Energy limit to match mass limit
  }
  Sublimation_m[1] <- VaporMassFlux_m

  Energy[1] <- S[1] + La[1] - Lt[1] + H[1] + E[1] + G + P[1]

  SnowDensity[1] <- ifelse((startingSnowDepth_m+NewSnow[1])>0,
                           min(Max_SnowDensity_kg_m3, (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity[1]*NewSnow[1])/(startingSnowDepth_m+NewSnow[1])),
                           Max_SnowDensity_kg_m3)

  # [MODIFIED] Update Melt Calculation
  # We must reduce the available SWE by the amount sublimated BEFORE calculating melt.
  SWE_After_Sublimation <- max(0, AvailableSWE + Sublimation_m[1])

  SnowMelt[1] <- max(0, min(SWE_After_Sublimation,
                            (Energy[1] - SnowHeatCap * SWE_After_Sublimation * WaterDens * (0-SnowTemp[1])) / (LatHeatFreez * WaterDens)))

  # [MODIFIED] Update Depth and SWE
  # Sublimation_m is added (it is negative for loss) to the mass balance
  SnowDepth[1]   <- max(0, (SnowWaterEq[1] + NewSnowWatEq[1] + Sublimation_m[1] - SnowMelt[1]) * WaterDens / SnowDensity[1])
  SnowWaterEq[1] <- max(0, SnowWaterEq[1] + NewSnowWatEq[1] + Sublimation_m[1] - SnowMelt[1])

  ##  Snow Melt Loop
  for (i in 2:length(precip_m)){
    # Albedo updates
    if (NewSnow[i] > 0){
      Albedo[i] <- 0.98-(0.98-Albedo[i-1])*exp(-4*NewSnow[i]*10)
    } else if (SnowDepth[i-1] < 0.1 | Albedo[i-1]< 0.3){
      Albedo[i] <- max(groundAlbedo, Albedo[i-1]+(groundAlbedo-0.85)/10)
    } else {
      Albedo[i] <- 0.35-(0.35-0.98)*exp(-1*(0.177+(log((-0.3+0.98)/(Albedo[i-1]-0.3)))^2.16)^0.46)
    }

    S[i] <- Solar(lat=lat, Jday=JDay[i], Tx=Tmax_C[i], Tn=Tmin_C[i], albedo=Albedo[i-1], forest=forest, aspect=aspect, slope=slope, printWarn=FALSE)

    if(SnowDepth[i-1] > 0) TE[i] <- 0.97

    if(SnowWaterEq[i-1] > 0 | NewSnowWatEq[i] > 0) {
      DCoef[i] <- 6.2
      if(SnowMelt[i-1] == 0){
        SnowTemp[i] <- max(min(0,Tmin_C[i]), min(0, (SnowTemp[i-1] + min(-SnowTemp[i-1], Energy[i-1]/((SnowDensity[i-1]*SnowDepth[i-1] + NewSnow[i]*NewSnowDensity[i])*SnowHeatCap*1000)))))
      }
    }

    # Rain on snow freezing
    if(SnowWaterEq[i-1] > 0 | R_m[i]>0) {
      if (SnowTemp[i]<0) {
        R_m_Snow[i] <- R_m[i]
        R_m[i] <- 0
        Ice_Snow_Depth[i] <- R_m_Snow[i]*WaterDens/IceDens
      }
    }

    rhos[i] <- SatVaporDensity(SnowTemp[i])
    H[i]    <- 1.29*(Tav[i]-SnowTemp[i])/rh[i]
    Lt[i]   <- Longwave(TE[i], SnowTemp[i])

    # [MODIFIED] Dynamic Latent Heat Selection
    # Choose Sublimation constant if snow exists, otherwise Vaporization
    CurrentLatentHeat <- lambdaV
    if(SnowWaterEq[i-1] > 0 | NewSnowWatEq[i] > 0) CurrentLatentHeat <- lambdaS

    # [MODIFIED] Calculate Vapor Energy using correct Latent Heat
    E[i] <- CurrentLatentHeat * (rhoa[i] - rhos[i]) / rh[i]

    # [ADDED] Mass Flux Calculation
    # Convert Energy Flux to Mass Flux (m). If E<0 (loss), this is negative.
    VaporMassFlux_m <- E[i] / CurrentLatentHeat / WaterDens

    # [ADDED] Mass Limiter
    # Check total available snow (Yesterday SWE + New Snow + Frozen Rain)
    AvailableSWE <- SnowWaterEq[i-1] + NewSnowWatEq[i] + R_m_Snow[i]

    if (VaporMassFlux_m < 0 && abs(VaporMassFlux_m) > AvailableSWE) {
      VaporMassFlux_m <- -AvailableSWE
      E[i] <- VaporMassFlux_m * WaterDens * CurrentLatentHeat # Recalculate Energy limit
    }
    Sublimation_m[i] <- VaporMassFlux_m

    # Energy Balance
    Energy[i] <- S[i] + La[i] - Lt[i] + H[i] + E[i] + G + P[i]

    # Snow Density Update
    if (Energy[i]>0) k <- 2 else k <- 1
    SnowDensity[i] <- ifelse((SnowDepth[i-1]+NewSnow[i])>0, min(Max_SnowDensity_kg_m3,
                                                                ((SnowDensity[i-1]+k*30*(Max_SnowDensity_kg_m3-SnowDensity[i-1])*exp(-DCoef[i]))*SnowDepth[i-1] +
                                                                   NewSnowDensity[i]*NewSnow[i]+Ice_Snow_Depth[i]*IceDens)/
                                                                  (SnowDepth[i-1]+NewSnow[i]+Ice_Snow_Depth[i])), Max_SnowDensity_kg_m3)

    # [MODIFIED] Snow Melt Calculation
    # Account for mass lost to sublimation before calculating melt
    SWE_After_Sublimation <- max(0, AvailableSWE + Sublimation_m[i])

    SnowMelt[i] <- max(0, min(SWE_After_Sublimation,
                              (Energy[i]-SnowHeatCap*(SWE_After_Sublimation)*WaterDens*(0-SnowTemp[i]))/(LatHeatFreez*WaterDens) ) )

    # [MODIFIED] Update Depth and SWE
    # Include Sublimation_m[i] in the mass balance
    SnowDepth[i]   <- max(0, (SnowWaterEq[i-1] + NewSnowWatEq[i] + R_m_Snow[i] + Sublimation_m[i] - SnowMelt[i]) * WaterDens / SnowDensity[i])
    SnowWaterEq[i] <- max(0,  SnowWaterEq[i-1] + NewSnowWatEq[i] + R_m_Snow[i] + Sublimation_m[i] - SnowMelt[i])
  }

  # [MODIFIED] Results Output
  # Added "Sublimation_mm" column
  Results <- data.frame(Date, Tmax_C, Tmin_C, precip_mm, R_m*1000, SnowDensity, R_m_Snow*1000,
                        NewSnowWatEq*1000, SnowMelt*1000, Sublimation_m*1000, NewSnow, SnowDepth, SnowWaterEq*1000)
  colnames(Results) <- c("Date", "MaxT_C", "MinT_C", "Precip_mm", "Rain_mm", 'Snow_Density', 'IceRain_mm',
                         "SnowfallWatEq_mm", "SnowMelt_mm", "Sublimation_mm", "NewSnow_m", "SnowDepth_m", "SnowWaterEq_mm")
  return(Results)
}
