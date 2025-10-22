#' EPA PBSM snow model
#' @description
#' Process-based snow model with six hydroclimate variables (prcp, air temperature, wind speed, relative humidity and short/longwave radiations) at a sub-daily time step (e.g., hourly data)
#'
#' @details
#' A modified version of PBSM which simulates snow processes with seven hydroclimate variables at sub-daily time step. All of hydroclimate variables are usually from weather forecast data (e.g., GEPS and CFS).
#'
#' @param Date Date(YYYY-MM-DD HH:MM)
#' @param time.step Time.step in hours (default=1)
#' @param precip_mm Precipitation in mm/day
#' @param Tair Air temperature in degree Celsius
#' @param windSp_ms Wind speed in m/sec
#' @param Relhum_frac Relative humidity in fraction
#' @param swrad_wm2 Incoming shortwave (solar) radiation in W/m2
#' @param lwrad_wm2 Incoming longwave radiation in W/m2
#' @param lat_deg Latitude in degree
#' @param tempHt Height of temperature station
#' @param windHt Height of wind station
#' @param groundAlbedo Ground albedo (default=0.25)
#' @param SurfEmissiv Surface Emissivity
#' @param forest Forest ratio
#' @param startingSnowDepth_m Initial value of snow depth in meter
#' @param startingSnowDensity_kg_m3 Initial value of snow density in kg/m3
#' @param Max_SnowDensity_kg_m3 Maximum snow density in kg/m3, default=550 from Voordendag et al.(2021)
#'
#' @return Time series of snow simulation
#' @export
#'
#' @examples
SnowMelt_EPA_Subdaily<-function(Date, time.step=1,precip_mm, Tair, windSp_ms,Relhum_frac,swrad_wm2,lwrad_wm2,
                       lat_deg, tempHt=1, windHt=2,groundAlbedo=0.25,SurfEmissiv=0.95,
                       forest=0, startingSnowDepth_m=0,startingSnowDensity_kg_m3=550,
                       Max_SnowDensity_kg_m3=550){
  ## Constants :
  WaterDens <- 1000			# kg/m3
  IceDens<- 917          # kg/m3
  lambda <- 3.35*10^5			# latent heat of fusion (kJ/m3)
  lambdaV <- 2500				# (kJ/kg) latent heat of vaporization
  SnowHeatCap <- 2.1			# kJ/kg/C
  LatHeatFreez <- 333.3		# kJ/kg
  Cw <- 4.2*10^3				# Heat Capacity of Water (kJ/m3/C)
  #Max_SnowDensity_kg_m3<-550 #Maximum snow density (Voordendag et al., 2021)

  ##	Converted Inputs :
  Tav <- Tair		# degrees C
  precip_m <- precip_mm*0.001	 	# precip in m
  R_m <- precip_m					# (m) depth of rain
  R_m[which(Tav < 0)] <- 0		# ASSUMES ALL SNOW at < 0 C
  #==========================================
  # Added by Hyung Eum
  R_m_Snow<-array(0,dim=c(length(R_m)))
  Ice_Snow_Depth<-R_m_Snow
  #=============================================================
  # Modified by A76=rho_min+max(kp(Ta+lp)^1.5) by Hyung Eum (Nov.2022)
  # Reference: Lafaysse et al. (2017)
  # Lafaysse, M., Cluzet, B., Dumont, M., Lejeune, Y., Vionnet, V., Morin, S., 2017.
  # A multiphysical ensemble system of numerical snow modelling.
  # The Cryosphere 11, 1173-1198. https://doi.org/10.5194/tc-11-1173-2017
  NewSnowDensity <- 50+(1.7*((Tav+15)^1.5))		# kg/m3
  NewSnowDensity[which(NewSnowDensity < 50)] <- 50
  NewSnowDensity[which(is.na(NewSnowDensity))] <- 50

  NewSnowWatEq <- precip_m				# m
  NewSnowWatEq[which(Tav >= 0)] <- 0			# No new snow if average temp above or equals 0 C
  NewSnow <- NewSnowWatEq*WaterDens/NewSnowDensity		# m
  #JDay <- strptime(Date, format="%Y-%m-%d")$yday+1
  lat <- lat_deg*pi/180		#	latitude in radians
  rh 	<- log((windHt+0.001)/0.001)*log((tempHt+0.0002)/0.0002)/(0.41*0.41*windSp_ms*(3600*time.step))	# (time.step/m) Thermal Resistance
  if (length(windSp_ms)==1) rh <- rep(rh,length(precip_mm))									##	creates a vector of rh values
  #-----------------------------------------------------------------------------------------
  #cloudiness 		<- EstCloudiness(Tmax_C,Tmin_C) # Cloudiness
  #-----------------------------------------------------------------------------------------
  #AE 				<- AtmosphericEmissivity(Tav, cloudiness)	# (-) Atmospheric Emissivity function
  #-----------------------------------------------------
  #  New Variables	:
  SnowTemp 		<- rep(0,length(precip_m)) 		# Degrees C
  #-----------------------------------------------------------------------------------------
  rhos 			<- SatVaporDensity(SnowTemp)	# 	vapor density at surface (kg/m3)
  rhoa 			<- SatVaporDensity(Tav)*Relhum_frac		#	vapor density of atmoshpere (kg/m3) Modified by Hyung ()
  #rhoa 			<- SatVaporDensity(Tmin_C)		#	vapor density of atmoshpere (kg/m3)
  #-----------------------------------------------------------------------------------------
  SnowWaterEq 	<- vector(length=length(precip_mm))		#  (m) Equiv depth of water
  TE 				<- rep(SurfEmissiv,length(precip_mm))	#	(-) Terrestrial Emissivity
  DCoef 			<- rep(0,length(precip_mm))				#   Density Coefficient (-) (Simplified version)
  SnowDensity 	<- rep(450,length(precip_mm))			#  (kg/m3)  Max density is 450
  SnowDepth 		<- vector(length=length(precip_mm))		#  (m)
  SnowMelt 		<- rep(0,length(precip_mm))				#  (m)
  Albedo 			<- rep(groundAlbedo,length(precip_mm)) 	#  (-) This will change for days with snow

  ##	Energy Terms
  H 		<- vector(length=length(precip_mm))	#	Sensible Heat exchanged (kJ/m2/d)
  E 		<- vector(length=length(precip_mm))	#	Vapor Energy	(kJ/m2/d)
  S_incoming 		<- swrad_wm2*3600*time.step*10^-3		#	Incoming Solar Radiation (kJ/m2) (convert W/m2 to kJ/m2) (1W/m2=1*3600*time.step*10^-3)
  S<-S_incoming
  S[1]<-(1-Albedo[1])*(1-forest)*S_incoming[1] # Net incident solar radiation
  #--------------------------------------------------------------------------------------
  La 		<- lwrad_wm2*3600*time.step*10^-3					#	Atmospheric Longwave Radiation (convert W/m2 to kJ/m2/d) (1W/m2=1*3600*24*10^-3)
  #--------------------------------------------------------------------------------------
  Lt 		<- vector(length=length(precip_mm))	#	Terrestrial Longwave Radiation (kJ/m2/d)
  G 		<- 2.0*3600*time.step*10^-3					#	Ground Conduction (kJ/m2/time.step) Assuming a constant value of 2.0 J/m2/second, U.S. Army Corps of Engineers, Runoff from Snowmelt (EM 1110-2-1406), Washington, DC, 1998
  P 		<- Cw * R_m * Tav					# 	Precipitation Heat (kJ/m2)
  Energy 	<- vector(length=length(precip_mm))	# Net Energy (kJ/m2)

  ##  Initial Values.
  SnowWaterEq[1] 	<- startingSnowDepth_m * startingSnowDensity_kg_m3 / WaterDens
  SnowDepth[1] 	<- startingSnowDepth_m
  Albedo[1] <- ifelse(NewSnow[1] > 0, 0.98-(0.98-0.50)*exp(-4*NewSnow[1]*10),ifelse(startingSnowDepth_m == 0, groundAlbedo, max(groundAlbedo, 0.5+(groundAlbedo-0.85)/10)))  # If snow on the ground or new snow, assume Albedo at previous time step was 0.5
  #--------------------------------------------------------------------------------------
  H[1] <- 1.29*(Tav[1]-SnowTemp[1])/rh[1]
  E[1] <- lambdaV*(rhoa[1]-rhos[1])/rh[1]
  if(startingSnowDepth_m>0) TE[1] <- 0.97
  Lt[1] <- Longwave(emissivity=TE[1],temp=SnowTemp[1],TimeStep=time.step)
  Energy[1] <- S[1] + La[1] - Lt[1] + H[1] + E[1] + G + P[1]
  SnowDensity[1] <- ifelse((startingSnowDepth_m+NewSnow[1])>0, min(Max_SnowDensity_kg_m3, (startingSnowDensity_kg_m3*startingSnowDepth_m + NewSnowDensity[1]*NewSnow[1])/(startingSnowDepth_m+NewSnow[1])), Max_SnowDensity_kg_m3)
  SnowMelt[1] <- max(0,	min((startingSnowDepth_m/10+NewSnowWatEq[1]),  # yesterday on ground + today new
                            (Energy[1]-SnowHeatCap*(startingSnowDepth_m/10+NewSnowWatEq[1])*WaterDens*(0-SnowTemp[1]))/(LatHeatFreez*WaterDens) ) )
  SnowDepth[1] <- max(0,(SnowDepth[1] + NewSnowWatEq[1]-SnowMelt[1])*WaterDens/SnowDensity[1])
  SnowWaterEq[1] <- max(0,SnowWaterEq[1]-SnowMelt[1]+NewSnowWatEq[1])

  ##  Snow Melt Loop
  for (i in 2:length(precip_m)){
    if (NewSnow[i] > 0){
      Albedo[i] <- 0.98-(0.98-Albedo[i-1])*exp(-4*NewSnow[i]*10)
    } else if (SnowDepth[i-1] < 0.1 | Albedo[i-1]< 0.3){  # Added a condition of "Albedo[i-1]< 0.3" to avoid NaN
      Albedo[i] <- max(groundAlbedo, Albedo[i-1]+(groundAlbedo-0.85)/10)
    } else Albedo[i] <- 0.35-(0.35-0.98)*exp(-1*(0.177+(log((-0.3+0.98)/(Albedo[i-1]-0.3)))^2.16)^0.46)

    #S[i] <- Solar(lat=lat,Jday=JDay[i], Tx=Tmax_C[i], Tn=Tmin_C[i], albedo=Albedo[i-1], forest=forest, aspect=aspect, slope=slope, printWarn=FALSE)
    S[i]<-(1-Albedo[i])*(1-forest)*S_incoming[i]

    if(SnowDepth[i-1] > 0) TE[i] <- 0.97 	#	(-) Terrestrial Emissivity
    if(SnowWaterEq[i-1] > 0 | NewSnowWatEq[i] > 0) {
      DCoef[i] <- 6.2
      if(SnowMelt[i-1] == 0){
        SnowTemp[i] <- max(min(0,Tav[i]),min(0,(SnowTemp[i-1]+min(-SnowTemp[i-1],Energy[i-1]/((SnowDensity[i-1]*
                                                                                                    SnowDepth[i-1]+NewSnow[i]*NewSnowDensity[i])*SnowHeatCap*1000)))))
      }

    }
    #=========================================
    # Rain on snow
    if(SnowWaterEq[i-1] > 0 | R_m[i]>0) {
      if (SnowTemp[i]<0) { # rain freezing in snowpack
        R_m_Snow[i]<-R_m[i]
        R_m[i]<-0
        Ice_Snow_Depth[i]<-R_m_Snow[i]*WaterDens/IceDens		# m (ice depth in snowpack)
      }
    }
    #=========================================
    rhos[i] <- SatVaporDensity(SnowTemp[i])
    H[i] <- 1.29*(Tav[i]-SnowTemp[i])/rh[i]
    E[i] <- lambdaV*(rhoa[i]-rhos[i])/rh[i]
    Lt[i] <- Longwave(emissivity=TE[i],temp=SnowTemp[i],TimeStep=time.step)
    Energy[i] <- S[i] + La[i] - Lt[i] + H[i] + E[i] + G + P[i]

    if (Energy[i]>0) k <- 2 else k <- 1

    #SnowDensity[i] <- ifelse((SnowDepth[i-1]+NewSnow[i])>0, min(Max_SnowDensity_kg_m3,
    #	((SnowDensity[i-1]+k*30*(Max_SnowDensity_kg_m3-SnowDensity[i-1])*exp(-DCoef[i]))*SnowDepth[i-1] + NewSnowDensity[i]*NewSnow[i])/(SnowDepth[i-1]+NewSnow[i])), Max_SnowDensity_kg_m3)
    SnowDensity[i] <- ifelse((SnowDepth[i-1]+NewSnow[i])>0, min(Max_SnowDensity_kg_m3,
                                                                ((SnowDensity[i-1]+k*30*(Max_SnowDensity_kg_m3-SnowDensity[i-1])*exp(-DCoef[i]))*SnowDepth[i-1] +
                                                                   NewSnowDensity[i]*NewSnow[i]+Ice_Snow_Depth[i]*IceDens)/
                                                                  (SnowDepth[i-1]+NewSnow[i]+Ice_Snow_Depth[i])), Max_SnowDensity_kg_m3)

    SnowMelt[i] <- max(0,	min( (SnowWaterEq[i-1]+NewSnowWatEq[i]),  # yesterday on ground + today new
                               (Energy[i]-SnowHeatCap*(SnowWaterEq[i-1]+NewSnowWatEq[i])*WaterDens*(0-SnowTemp[i]))/(LatHeatFreez*WaterDens) )  )

    SnowDepth[i] <- max(0,(SnowWaterEq[i-1]+NewSnowWatEq[i]-SnowMelt[i]+R_m_Snow[i])*WaterDens/SnowDensity[i])

    SnowWaterEq[i] <- max(0,SnowWaterEq[i-1]-SnowMelt[i]+NewSnowWatEq[i]+R_m_Snow[i])	# (m) Equiv depth of water
  }

  Results<-data.frame(Date, T_C, precip_mm, R_m*1000, SnowDensity,R_m_Snow*1000,NewSnowWatEq*1000,SnowMelt*1000, NewSnow, SnowDepth, SnowWaterEq*1000)
  colnames(Results)<-c("Date", "Tair_C", "Precip_mm", "Rain_mm",'Snow_Density','IceRain_mm',"SnowfallWatEq_mm", "SnowMelt_mm", "NewSnow_m", "SnowDepth_m", "SnowWaterEq_mm")
  return(Results)
}
