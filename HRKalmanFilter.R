## ----setup, include=FALSE----------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(Ranadu)
source('chunks/getPlatform.R')  ## reads "Platform" or "platform_identifier"


## ----include=TRUE, echo = TRUE-----------------------------------------------------
VarList <- c('LAT', 'LON', 'VEW', 'VNS', 'PITCH', 'ROLL', 'THDG',
             'LATKF', 'LONKF', 'ALTKF', 'VEWKF', 'VNSKF', 'ROCKF', 'ATX',
             'AKKF', 'SSKF', 'ADIFR', 'BDIFR', 'QCF', 'TASX', 'GGVSPD', 'VSPD',
             'GGVEW', 'GGVNS',
             'PITCHKF', 'ROLLKF', 'THDGKF', 'PSXC', 'GGLAT', 'GGALT', 'ACINS')
fname1 <- setFileName('WCR-TEST', 'tf02KF')
D1 <- getNetCDF(fname1, VarList)
Rate <- 1
D1$Grav <- Gravity (D1$LAT, D1$GGALT)
Cradeg <- pi/180
Data <- D1  # ROC uses data.frame Data
source('chunks/ROC.R')
Cor <- rep(0, nrow(D1)*9)
dim(Cor) <- c(nrow(D1), 9)
D1 <- Data
Cor[, 1] = with(D1, LATKF - LAT)
Cor[, 2] = with(D1, LONKF - LON)
Cor[, 3] = with(D1, ALTKF - ZROC)
Cor[, 4] = with(D1, VEWKF - VEW)
Cor[, 5] = with(D1, VNSKF - VNS)
Cor[, 6] = with(D1, ROCKF - ROC)
Cor[, 7] = with(D1, PITCHKF - PITCH)
Cor[, 8] = with(D1, ROLLKF - ROLL)
Cor[, 9] = with(D1, THDGKF - THDG)
rm(Data)
## Get the needed variables from the 25-Hz file
fname25 <- setFileName('WCR-TEST', 'tf02h')
FI <- DataFileInfo(fname25)
VarList <- VarList[-which(!(VarList %in% FI$Variables))]
D25 <- getNetCDF(fname25, VarList)
D25$Grav <- Gravity (D25$LAT, D25$GGALT)
Data <- D25
Rate <- 25
source('chunks/ROC.R')
D25 <- Data
rm(Data)


## ----interp, include = TRUE--------------------------------------------------------
## interpolate the results over the full data.frame:
IntFilter <- function (X, inRate, outRate) {
  if (inRate == outRate) {return (X)}
  ratio <- as.integer(outRate/inRate)    ## expected to be an integer
  x <- 0:(length(X)-1)
  A <- stats::approx (x, X, n=length(X)*ratio-ratio+1)
  T <- A$y
  T <- signal::filter(signal::sgolay(4,75),T)
  ## now shift to match Rate:
  n <- as.integer (ratio / 2)
  NL = length(T)
  T <- c(rep(T[1],n), T, rep(T[NL],ratio-n-1))  ## OK, even or odd ratio
  return (T)
}
## find time overlap
SE1 <- getStartEnd(D1)
SE25 <- getStartEnd(D25)
nsec <- function(T) {
  as.integer(T %% 100 + (T %/% 100) * 60 + (T %/% 10000) * 3600)
}
if(SE1[1] > SE25[1]) {
  gap <- nsec(SE1[1]) - nsec(SE25[1])
  for (n in 1:gap) {
    Cor <- rbind(Cor[1,], Cor)
  }
} else if (SE1[1] < SE25[1]) {
  gap <- nsec(SE25[1]) - nsec(SE1[1])
  Cor <- Cor[(gap+1):nrow(Cor), ]
}
if(SE1[2] > SE25[2]) {
  gap <- nsec(SE1[2]) - nsec(SE25[2])
  Cor <- Cor[1:(nrow(Cor) - gap), ]
} else if (SE1[2] < SE25[2]) {
  gap <- nsec(SE25[2]) - nsec(SE1[2])
  for (n in 1:gap) {
    Cor <- rbind(Cor, Cor[nrow(Cor), ])
  }
}
C25 <- vector('numeric', nrow(D25) * 9)
dim (C25) <- c(nrow(D25), 9)
for (j in 1:9) {
  C25[, j] <- IntFilter (Cor[, j], 1, 25)
}


## ----newvars, include = TRUE-------------------------------------------------------
D25$LATKF <- D25$LAT + C25[, 1]
D25$LONKF <- D25$LON + C25[, 2]
D25$ALTKF <- D25$ZROC + C25[, 3]
D25$VEWKF <- D25$VEW + C25[, 4]
D25$VNSKF <- D25$VNS + C25[, 5]
D25$ROCKF <- D25$ROC + C25[, 6]
D25$PITCHKF <- D25$PITCH + C25[, 7]
D25$ROLLKF  <- D25$ROLL + C25[, 8]
D25$THDGKF  <- D25$THDG + C25[, 9]


## ----radome, include = TRUE--------------------------------------------------------
Platform <- ifelse (grepl('130', getPlatform(fname1)), "C130", "GV")
cffn <- 19.70547
cff <- 21.481
cfs <- c(4.525341674, 19.933222011, -0.001960992)
if (Platform == "C130") { # C-130 WECAN values
  cff <- 10.3123
  cfs <- c(5.6885, 14.0452, -0.00461)
}
D25$QR <- D25$ADIFR / D25$QCF
D25$QR[D25$QCF < 20] <- NA
D25$QR[is.infinite(D25$QR)] <- NA
CutoffFreq <- 600 * Rate
D25$QRS <- zoo::na.approx (as.vector(D25$QR), maxgap=1000*Rate, na.rm = FALSE, rule=2)
D25$QRS[is.na(D25$QRS)] <- 0
D25$QRS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D25$QRS)
D25$QRF <-  D25$QR - D25$QRS
D25$QCFS <- zoo::na.approx (as.vector(D25$QCF), maxgap=1000*Rate, na.rm = FALSE)
D25$QCFS[is.na(D25$QCFS)] <- 0
D25$QCFS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D25$QCFS)
D25$AKKF <- cff * D25$QRF + cfs[1] + cfs[2] * D25$QRS + cfs[3] * D25$QCFS
D25$AKKF[D25$QCF < 10] <- NA
D25$SSKF <- 0.008 + 22.301 * D25$BDIFR / D25$QCF
if (grepl('130', getPlatform(fname25))) { # C-130 WECAN values; use with heading offset +0.76
  D25$SSKF <- 0.85 + 12.6582 * D25$BDIFR / D25$QCF
}
D25$SSKF[D25$QCF < 10] <- NA


## ----getWind, include = TRUE-------------------------------------------------------
## This follows "chunks/wind-calculation.R" from the KalmanFilter directory
DataW <- D25
DataW$ATTACK <- DataW$AKKF  ## These are the values used in WindProcessor()
DataW$SSLIP <- DataW$SSKF
DataN <- WindProcessor (DataW, AC = Platform)  ## This function is in Ranadu
D25$WICC <- DataN$WIN
D25$WDCC <- DataN$WDN
D25$WSCC <- DataN$WSN
DataW$PITCH <- D25$PITCHKF
DataW$ROLL <- D25$ROLLKF
DataW$THDG <- D25$THDGKF
DataW$VEW <- D25$VEWKF
DataW$VNS <- D25$VNSKF
DataW$GGVSPD <- D25$ROCKF
DataN <- WindProcessor(DataW, AC = Platform, CompF=FALSE)    ## suppress comp filter and GPS lever arm)
D25$WDKF <- DataN$WDN
D25$WSKF <- DataN$WSN
D25$WIKF <- DataN$WIN
## add longitudinal and lateral components analogous to UXC and VYC:
.hdg <- D25$THDGKF * Cradeg
.wd <- D25$WDKF * Cradeg + pi
D25$UXKF <- D25$WSKF * (sin(.hdg)*sin(.wd) + cos(.hdg)*cos(.wd))
.hdg <- .hdg - pi/2
D25$VYKF <- D25$WSKF * (sin(.hdg)*sin(.wd) + cos(.hdg)*cos(.wd))
DataW$GGVSPD <- D25[, 'GGVSPD']   ## std script has code to substitute for GGVSPD
DataN <- WindProcessor(DataW, AC = Platform, CompF=FALSE)    ## suppress comp filter and GPS lever arm)
D25$WIKFG <- DataN$WIN


## ----newCDF, include = TRUE--------------------------------------------------------
## create-new-netcdf.R

fnew <- sub ('h.nc', 'hKF.nc', fname25)
## beware: overwrites without warning!!
Z <- file.copy (fname25, fnew, overwrite=TRUE)  ## BEWARE: overwrites without warning!!

# function to copy attributes from old variable (e.g., PITCH) to new one (e.g., PITCHKF)
copy_attributes <- function (atv, v, nfile) {
  for (i in 1:length(atv)) {
    aname <- names(atv[i])
    if (grepl ('name', aname)) {next}  # skips long and standard names
    if (grepl ('units', aname)) {next}
    if (grepl ('Dependencies', aname)) {next}
    if (grepl ('actual_range', aname)) {next}
    if (is.numeric (atv[[i]])) {
      ncatt_put (nfile, v, attname=aname, attval=as.numeric(atv[[i]]))
    } else {
      ncatt_put (nfile, v, attname=aname, attval=as.character (atv[[i]]))
    }
  }
}


## ----writeCDF, include = FALSE-----------------------------------------------------
## modify-new-netcdf.R

## variables needed for attributes:

VarListRef <- VarList
FI <- DataFileInfo (fname25)

if (!('GGVSPD' %in% FI$Variables)) {
  if ('GGVSPDB' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'GGVSPDB'
  } else if ('VSPD_A' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'VSPD_A'
  } else if ('VSPD_G' %in% FI$Variables) {
    VarList [which (VarList == 'GGVSPD')] <- 'VSPD_G'
  } else {
    print ('ERROR: no VSPD variable found')
    exit()
  }
}
for (Var in VarList) {
  if (!(Var %in% FI$Variables)) {
    print (sprintf (' required variable %s not found in file %s; skipping...', Var, fname))
    exit()
  }
}

netCDFfile <- nc_open (fnew, write=TRUE) 
Rate <- 25
Dimensions <- attr (D25, "Dimensions")
Dim <- Dimensions[["Time"]]
if ("sps25" %in% names (Dimensions)) {
  Rate <- 25
  Dim <- list(Dimensions[["sps25"]], Dimensions[["Time"]])
}
if ("sps50" %in% names (Dimensions)) {
  Rate <- 50
  Dim <- list(Dimensions[["sps50"]], Dimensions[["Time"]])
}
DATT <- D25  ## save to ensure that attributes are preserved

## variables to add to the netCDF file:
VarNew <- c('LATKF', 'LONKF', 'ALTKF', 'VEWKF', 'VNSKF', 'ROCKF', 'PITCHKF', 'ROLLKF', 'THDGKF',
            'WDKF', 'WSKF', 'WIKF', 'WICC', 'WDCC', 'WSCC', 'UXKF', 'VYKF', 
            'AKKF', 'SSKF')
VarOld <- c('LAT', 'LON', 'ALT', 'VEW', 'VNS', 'VSPD', 'PITCH', 'ROLL', 'THDG', 
            'WD', 'WS', 'WI', 'WIC', 'WDC', 'WSC', 'UXC', 'VYC', 'AKRD', 'SSRD')
VarUnits <- c('degrees', 'degrees', 'm', 'm/s', 'm/s', 'm/s', 'degrees', 'degrees', 'degrees',
              'degrees', 'm/s', 'm/s', 'm/s', 'degrees', 'm/s', 'm/s', 'm/s',
              'degrees', 'degrees')
VarStdName <- c('latitude, KF', 'longitude, KF', 'altitude MSL, KF',
                 'eastward groundspeed, KF', 'northward groundspeed, KF', 'rate of climb, KF',
                 'pitch, KF', 'roll, KF', 'heading, KF',
                 'wind direction, KF', 'wind speed, KF', 'vertical wind, KF', 
                 'vertical wind, KF and GGVSPD', 'WIC recalc', 'WDC recalc', 'WSC recalc',
                 'longitudinal wind, KF', 'lateral wind, KF', 
                 'angle of attack, KF', 'sideslip, KF')
VarLongName <- c('INS latitude, Kalman-filter-corrected',
                'INS longitude, Kalman-filter-corrected',
                'INS altitude, Kalman-filter-corrected',
                'INS eastward ground speed, Kalman-filter-corrected',
                'INS northward ground speed, Kalman-filter-corrected',
                'INS rate of climb, Kalman-filter-corrected',
                'INS aircraft pitch angle, Kalman-filter-corrected',
                'INS aircraft roll angle, Kalman-filter-corrected',
                'INS aircraft true heading angle, Kalman-filter-corrected',
                'horizontal wind direction, Kalman-filter-corrected',
                'horizontal wind speed, Kalman-filter-corrected',
                'vertical wind speed, Kalman-filter-corrected',
                'original WIC, recalculated',
                'original WDC, recalculated',
                'original WSC, recalculated',
                'longitudinal component of the horizontal wind, Kalman-filter-corrected',
                'lateral component of the horizontal wind, Kalman-filter-corrected',
                'angle of attack, from Kalman-filter routine', 'sideslip, from Kalman-filter routine')

## create the new variables
varCDF <- list ()
for (i in 1:length(VarNew)) {
  print (sprintf ('new-netcdf %d%% done', as.integer(100*(i-1)/length(VarNew))))
  varCDF[[i]] <- ncvar_def (VarNew[i],  
                            units=VarUnits[i], 
                            dim=Dim, 
                            missval=as.single(-32767.), prec='float', 
                            longname=VarLongName[i])
  if (i == 1) {
    newfile <- ncvar_add (netCDFfile, varCDF[[i]])
  } else {
    newfile <- ncvar_add (newfile, varCDF[[i]])
  }
  ATV <- ncatt_get (netCDFfile, VarOld[i])
  copy_attributes (ATV, VarNew[i], newfile)
  ncatt_put (newfile, VarNew[i], attname="standard_name", 
             attval=VarStdName[i])
  ncvar_put (newfile, varCDF[[i]], D25[, VarNew[i]], count=c(25, nrow(D25)/25))
}
nc_close (newfile)


