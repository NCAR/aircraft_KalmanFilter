## Running the KalmanFilter.R script:
## This script can be run in three ways:
## 1. As an R script using a command in a terminal window following this pattern:
##    Rscript KalmanFilter.R Project Flight UpdateAKRD[y/N] UpdateSSRD[y,N] Simple[y/N] Interval[10]")
##    Example: Rscript KalmanFilter.R CSET rf01 y y n 15
##    You can also append a directory for output. The default is DataDirectory()/Project/KF/,
##    e.g., /scr/raf_data/WCR-TEST/KF/
## 2. Interactively from an R console or RStudio console:
##    Set working directory to the KalmanFilter directory.
##    Using R:
##      Start R;
##      Enter the command "source('KalmanFilter.R')" (without double-quote marks)
##      Answer the questions interactively to select the project/flight, etc.
##    Using RStudio:
##      Set the working directory to the KalmanFilter directory.
##      Load the KalmanFilter.R file.
##      Click "Source" and answer questions interactively.
##      (The above instructions for "Using R" also work from the RStudio console.)
## 3. Use the KalmanFilter Shiny app:
##   Start the Shiny app on datavis using this address:
##     https://datavis.eol.ucar.edu:3838/KalmanFilter
##   Optionally, click the 'help' button.
##   Enter the appropriate parameters for processing in the window that appears.
##   Click the "Click here ..." button.
##
##
## The run argument form "Flight" can be entered in several different ways:
##   a). To process rf05, enter either "5" or "rf05".
##   b). To process rf05h, enter "rf05h".
##   c). To process tf02, enter "tf02". ("tf" or "ff" is required for test/ferry flights)
##   d). To process the next unprocessed low-rate flight, enter "NEXT" (without the double-quote marks)
##   e). To process all unprocessed low-rate flights in a given project, enter "ALL".
##   f). To process the next unprocessed high-rate flight, enter "NEXTHR".
##   g). To process all unprocessed high-rate flights, enter "ALLHR".
## The action differs for these choices if a processed file already exists:
##   For a), b) or c) the processed file will be overwritten without warning.
##   For d), e) f) or g) only files not previously processed will be processed.

print (sprintf ('Kalman Processor is starting -- %s', Sys.time()))
thisFileName <- "KalmanFilter"
library(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
library(numDeriv, quietly = TRUE, warn.conflicts=FALSE)    ## needed for the jacobian() function
# library(signal)      ## used for filtering -- but ::signal avoids error message
library(reshape2)    ## used with ggplot facet plots
library(grid)
options(stringsAsFactors=FALSE)
source('chunks/getPlatform.R')  ## reads "Platform" or "platform_identifier"

## Data location: use $DAT if it names an existing directory; otherwise
## Ranadu's DataDirectory() falls back to $DATA_DIR, then /scr/raf_data/, etc.
## $DAT may be /scr/raf/local_productiondata if running as nimbus or include the
## project (/scr/raf_data/TI3GER-2) if not running as nimbus; the project part
## is removed because this script and Ranadu's setFileName() both append
## Project/ themselves.
useDAT <- function(Project) {
  dat <- sub('/+$', '', Sys.getenv('DAT'))
  if (!nzchar(dat)) return(invisible(NULL))
  if (!dir.exists(dat)) {
    print(sprintf('DAT=%s is not an existing directory; ignoring it', dat))
    return(invisible(NULL))
  }
  if (basename(dat) == Project) dat <- dirname(dat)
  Sys.setenv(DATA_DIR = dat)
}

## See if this is being run as user "shiny":
SHINY <- ifelse(Sys.info()['user'] == 'shiny', TRUE, FALSE)

# setwd ('~/RStudio/KalmanFilter')

##----------------------------------------------------------
## These are the run options to set via command-line or UI:
GeneratePlots <- TRUE
ShowPlots <- FALSE  ## leave FALSE for batch runs or shiny-app runs
UpdateAKRD <- FALSE
UpdateSSRD <- FALSE
SimpleOnly <- FALSE   ## set TRUE to use CorrectPitch/CorrectHeading w/o KF
ALL <- FALSE
ALLHR <- FALSE
NEXT <- FALSE
NEXTHR <- FALSE
NSTEP <- 15      ## update interval
HighRate <- FALSE

Directory <- DataDirectory ()
KFDirectory <- ifelse(SHINY, 'KFoutput', Directory)
##----------------------------------------------------------

Project <- 'WCR-TEST'
ProjectDir <- Project
Flight <- 1
## Functions to supply next flight for which processing is needed:
getNext <- function(ProjectDir, Project) {
  Fl <- sort (list.files (sprintf ("%s%s/", Directory, ProjectDir),
                          sprintf ("%srf..KF.nc", Project)), decreasing = TRUE)[1]
  ## Consider also if a processed version exists in KF or KFoutput:
  Fl <- sort (c(Fl, 
                list.files (sprintf ("%s%s/KF/", Directory, ProjectDir),
                            sprintf ("%srf..KF.nc", Project)),
                list.files ("KFoutput", sprintf ("%srf..KF.nc", Project))),
              decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 1
  } else {
    Flight <- sub (".*rf", '',  sub ("KF.nc", '', Fl))
    Flight <- as.numeric(Flight)+1
  }
  return (Flight)
}

getNextHR <- function(ProjectDir, Project) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf..hKF.nc", Project)), decreasing = TRUE)[1]
  ## Consider also if a processed version exists in KF or KFoutput:
  Fl <- sort (c(Fl, 
                list.files (sprintf ("%s%s/KF/", DataDirectory (), ProjectDir),
                            sprintf ("%srf..hKF.nc", Project)),
                list.files ("KFoutput", sprintf ("%srf..hKF.nc", Project))),
              decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 1
  } else {
    Flight <- sub (".*rf", '',  sub ("hKF.nc", '', Fl))
    Flight <- as.numeric(Flight)+1
  }
  return (Flight)
}

if (!interactive()) {  ## can run interactively or via Rscript; this is for Rscript:
  run.args <- commandArgs (TRUE)
  if (length (run.args) > 0) {
    if (nchar(run.args[1]) > 1) {
      Project <- run.args[1]
      ProjectDir <- Project
    }
    useDAT(Project)
    Directory <- DataDirectory ()
    KFDirectory <- ifelse(SHINY, 'KFoutput', Directory)
    print(sprintf('Data directory: %s', Directory))
  } else {
    print ("Usage: Rscript KalmanFilter.R Project Flight UpdateAKRD[y/N] UpdateSSRD[y,N] Simple[y/N] Interval[10]")
    print ("Example: Rscript KalmanFilter.R CSET rf01 y y n 15")
    stop("exiting...")
  }
  ## Flight
  if (length (run.args) > 1) {
    if (run.args[2] == 'ALL') {
      flight <- run.args[2]
      ALL <- TRUE
    } else if (run.args[2] == 'ALLHR') {
      flight <- run.args[2]
      ALLHR <- TRUE
      HighRate <- TRUE
    } else if (run.args[2] == 'NEXT') {
      ## Find max rf in data directory,
      ## Use as default if none supplied via command line:
      flight <- sprintf('rf%02d', getNext(ProjectDir, Project))
    } else if (run.args[2] == 'NEXTHR') {
      HighRate <- TRUE
      ## Find max rf in data directory,
      ## Use as default if none supplied via command line:
      flight <- sprintf('rf%02d', getNextHR(ProjectDir, Project))
      HighRate <- TRUE
    } else {
      ## Possible assignments are 1,2,..., rf01,..., rf01h,..., tf01,..., ff01..., tf01h,..., ff01h,...
      if (!grepl('f', run.args[2])) { ## only number is supplied; interpret as rf
        fno <- as.numeric(run.args[2])
        flight <- paste0('rf', sprintf('%02d', fno))
      } else {
        flight <- run.args[2]
        if (grepl('h', run.args[2])) {
          HighRate <- TRUE
        }
      }
    }  
  }
  if (length (run.args) > 2) {
    UpdateAKRD <- run.args[3]
    if (UpdateAKRD == 'y' || UpdateAKRD == 'Y' || UpdateAKRD == 'TRUE') {
      UpdateAKRD <- TRUE
    } else {
      UpdateAKRD <- FALSE
    }
  }
  if (length (run.args) > 3) {
    UpdateSSRD <- run.args[4]
    if (UpdateSSRD == 'y' || UpdateSSRD == 'Y' || UpdateSSRD == 'TRUE') {
      UpdateSSRD <- TRUE
    } else {
      UpdateSSRD <- FALSE
    }
  }
  if (length (run.args) > 4) {
    SimpleOnly <- run.args[5]
    if (SimpleOnly == 'y' || SimpleOnly == 'Y' || SimpleOnly == 'TRUE') {
      SimpleOnly <- TRUE
    } else {
      SimpleOnly <- FALSE
    }
  }
  if (length (run.args) > 5) {
    NSTEP <- as.numeric (run.args[6])
  }
  if (length (run.args) > 6) {
    outputDirectory <- run.args[7]
  } else {
    outputDirectory <- paste0(DataDirectory(), Project, '/KF/')
  }
  if (length (run.args) > 7) {
    GeneratePlots <- run.args[8]
    if (GeneratePlots == 'y' || GeneratePlots == 'Y' || GeneratePlots == 'TRUE') {
      GeneratePlots <- TRUE
    } else {
      GeneratePlots <- FALSE
    }
  }
} else {  ## This is the interactive section
  x <- readline (sprintf ("Project is %s; CR to accept or enter new project name: ", Project))
  if (nchar(x) > 1) {
    Project <- x
    if (grepl ('HIPPO', Project)) {
      ProjectDir <- 'HIPPO'
    } else {
      ProjectDir <- Project
    }
  }
  x <- readline (sprintf ("Flight is %s; CR to accept, number or 'ALL' or 'NEXT' for new flight name: ", Flight))
  if (x == 'ALL') {
    ALL <- TRUE
    flight <- 'ALL'
  } else if (x == 'ALLHR') {
    ALLHR <- TRUE
    flight <- 'ALLHR'
    HighRate <- TRUE
  } else if (x == 'NEXT') {
    flight <- sprintf('rf%02d', getNext(ProjectDir, Project))
  } else if (x == 'NEXTHR') {
    flight <- sprintf('rf%02d', getNextHR(ProjectDir, Project))
    HighRate <- TRUE
  } else if (nchar(x) > 0) {
    if (!is.na(as.numeric(x))) {
      flight <- sprintf('rf%02d', as.numeric(x))
    } else {
      flight <- x
      if (grepl('h', flight)) {
        HighRate <- TRUE
      }
    }
  }
  x <- readline (sprintf ("newAK is %s; CR to accept, Y or T to enable, N or F to disable: ", UpdateAKRD))
  if (nchar(x) > 0) 
    UpdateAKRD <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("newSS is %s; CR to accept, Y or T to enable, N or F to disable: ", UpdateSSRD))
  if (nchar(x) > 0) 
    UpdateSSRD <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("simple is %s; CR to accept, Y or T to enable, N or F to disable: ", SimpleOnly))
  if (nchar(x) > 0) 
    SimpleOnly <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- readline (sprintf ("NSTEP is %d, CR to accept or number to change: ", NSTEP))
  if (nchar(x) > 0 && !is.na(as.numeric(x))) {
    NSTEP <- as.numeric(x)
  }
  outputDirectory <- paste0(DataDirectory(), Project, '/KF/')
  x <- readline (sprintf ("Output directory is %s; CR to accept, full path entry to change", outputDirectory))
  if (nchar(x) > 0) {
    outputDirectory <- x
  }
  x <- readline (sprintf ("Generate Plots is %s; CR to accept, Y or T to enable, N or F to disable: ", GeneratePlots))
  if (nchar(x) > 0) 
    GeneratePlots <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
}

print (sprintf ('run controls:  Project: %s;  Flight: %s;  UpdateAKRD: %s;  UpdateSSRD: %s;  SimpleOnly: %s;  Time increment: %d',
                Project, flight, UpdateAKRD, UpdateSSRD, SimpleOnly, NSTEP))
if (!SHINY) {
  print (sprintf ('The output directory will be %s', outputDirectory))
}
## Here is the "ALL" loop:
if (ALL) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf...nc", Project)))
  FlKF <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                            sprintf ("%srf..KF.nc", Project)))
  FlKF <- sort (c(FlKF, 
                  list.files (sprintf ("%s%s/KF/", DataDirectory (), ProjectDir),
                              sprintf ("%srf..KF.nc", Project)),
                  list.files ('KFoutput', 
                                    sprintf ("%srf..KF.nc", Project))))
  if (is.na (Fl[1])) {
    print (sprintf ('no files to process for project %s', Project))
    exit()
  } else {
    ip <- pmatch(sub('.nc', '', Fl), sub('KF.nc', '', FlKF))
    if (all(!is.na(ip))) {
      print (sprintf ('no unprocessed files in project directory for %s', Project))
      exit()
    } else if (any(!is.na(ip))) {
      Fl <- Fl[-ip[!is.na(ip)]]
    }
    ## now process each flight in Fl
    Fl <- sub('.nc', '', sub(Project, '', Fl))
    print(' ALL loop: These files will be processed:')
    print(Fl)
  }
} else if(ALLHR) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf..h.nc", Project)))
  FlKF <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                            sprintf ("%srf..hKF.nc", Project)))
  FlKF <- sort (c(FlKF, 
                  list.files (sprintf ("%s%s/KF/", DataDirectory (), ProjectDir),
                              sprintf ("%srf..hKF.nc", Project)),
                  list.files ('KFoutput', 
                                    sprintf ("%srf..hKF.nc", Project))))
  if (is.na (Fl[1])) {
    print (sprintf ('no files to process for project %s', Project))
    exit()
  } else {
    ip <- pmatch(sub('.nc', '', Fl), sub('KF.nc', '', FlKF))
    if (all(!is.na(ip))) {
      print (sprintf ('no unprocessed files in project directory for %s', Project))
      exit()
    } else if (any(!is.na(ip))) {
      Fl <- Fl[-ip[!is.na(ip)]]
    }
    ## now process each flight in Fl
    Fl <- sub('.nc', '', sub(Project, '', Fl))
    print(' ALLHR loop: These files will be processed:')
    print(Fl)
    HighRate <- TRUE
  }
} else {
  Fl <- flight
}
for (flight in Fl) {
  if(HighRate) {
    fname = sprintf("%s%s/%s%s.nc", Directory, ProjectDir, Project, flight)
    ## Check if the LR KF file exists; if not, branch out of loop
    fcheck <- sub('h.nc', 'KF.nc', fname)
    fcheck2 <- paste0('KFoutput/', sub('.*/', '', sub('h.nc$', 'KF.nc', fname)))
    fcheck3 <- sub(Project, paste0(Project, '/KF'), fcheck)
    if (!file.exists(fcheck) && !file.exists(fcheck2) && !file.exists(fcheck3)) {
      print (sprintf ('LR processed KF file %s or %s not present; skipping file %s', 
                      fcheck, fcheck2, fname))
      next
    }
  } else {
    fname = sprintf("%s%s/%s%s.nc", Directory, ProjectDir, Project, flight)
  }
  
  # function to copy attributes from old variable (e.g., PITCH) to new one (e.g., PITCHKF)
  # (placed here because it is used for both LR and HR processing)
  copy_attributes <- function (atv, v, nfile) {
    for (i in 1:length(atv)) {
      aname <- names(atv[i])
      if (grepl ('name', aname)) {next}  # skips long and standard names
      if (grepl ('units', aname)) {next}
      if (grepl ('Dependencies', aname)) {next}
      if (grepl ('actual_range', aname)) {next}
      if (aname == '_FillValue') {next}  # already set by ncvar_def(missval=); type must match
      if (is.numeric (atv[[i]])) {
        ncatt_put (nfile, v, attname=aname, attval=as.numeric(atv[[i]]))
      } else {
        ncatt_put (nfile, v, attname=aname, attval=as.character (atv[[i]]))
      }
    }
  }
  
  ## There are two main segments, the first for HR and the second for LR.
  ## That permits both to be merged into one script.
  if (HighRate) {
    print (sprintf ('reading LR data -- %s', Sys.time()))
    VarList <- c('LAT', 'LON', 'VEW', 'VNS', 'PITCH', 'ROLL', 'THDG',
                 'LATKF', 'LONKF', 'ALTKF', 'VEWKF', 'VNSKF', 'ROCKF', 'ATX',
                 'AKKF', 'SSKF', 'ADIFR', 'BDIFR', 'QCF', 'TASX', 'GGVSPD', 'VSPD',
                 'GGVEW', 'GGVNS',
                 'PITCHKF', 'ROLLKF', 'THDGKF', 'PSXC', 'GGLAT', 'GGALT', 'ACINS')
    fname1 <- setFileName(Project, paste0(sub('h$', '', flight), 'KF'))
    fname1b <- paste0('KFoutput/', Project, sub('h$', '', flight), 'KF.nc')
    fname1c <- sub(sprintf('%s/', Project), sprintf('%s/KF/', Project), fname1)
    if(file.exists(fname1b)) {
      fname1 <- fname1b
    }
    if(file.exists(fname1c)) {
      fname1 <- fname1c
    }
    D1 <- getNetCDF(fname1, VarList)
    Rate <- 1
    D1$Grav <- Gravity (D1$LAT, D1$GGALT)
    Cradeg <- pi/180
    Data <- D1  # ROC uses data.frame Data
    print (sprintf ('Kalman-loop 5%% done'))
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
    print (sprintf ('reading HR data -- %s', Sys.time()))
    print (sprintf ('Kalman-loop 10%% done'))
    ## Get the needed variables from the 25-Hz file
    fname25 <- setFileName(Project, flight)
    FI <- DataFileInfo(fname25)
    VarList <- VarList[-which(!(VarList %in% FI$Variables))]
    D25 <- getNetCDF(fname25, VarList)
    D25$Grav <- Gravity (D25$LAT, D25$GGALT)
    Data <- D25
    Rate <- 25
    print (sprintf ('ROC HR  -- %s', Sys.time()))
    print (sprintf ('Kalman-loop 15%% done'))
    source('chunks/ROC.R')
    D25 <- Data
    rm(Data)
    
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
      if (gap > 0) {
        for (n in 1:gap) {
          Cor <- rbind(Cor, Cor[nrow(Cor), ])
        }
      }
    }
    print (sprintf ('interpolation -- %s', Sys.time()))
    print (sprintf ('Kalman-loop 25%% done'))
    C25 <- vector('numeric', nrow(D25) * 9)
    dim (C25) <- c(nrow(D25), 9)
    for (j in 1:9) {
      # print (sprintf ('j=%d, nrow(D25)=%d', j, nrow(D25)))
      # print (str(IntFilter (Cor[, j], 1, 25)))
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
    print (sprintf ('wind calculation -- %s', Sys.time()))
    print (sprintf ('Kalman-loop 50%% done'))
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
    print (sprintf ('Kalman-loop 75%% done'))
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
    print (sprintf ('making new netCDF 25%%'))
    
    ## ----newCDF, include = TRUE--------------------------------------------------------
    ## create-new-netcdf.R
    print (sprintf ('create new netCDF file -- %s', Sys.time()))
    
    fnew <- sub ('h.nc', 'hKF.nc', fname25)
    fnew <- paste0(outputDirectory, sub('.*/', '', fnew))
    if(SHINY) {  ## for the datavis shiny server, only this output is allowed:
      fnew <- paste0('KFoutput/', sub('.*/', '', fnew))
    }
    print (sprintf ('making new netCDF file named %s', fnew))
    ## Make the directory KF if needed:
    if (grepl ('/KF/', fnew)) {
      if (!dir.exists(outputDirectory)) {
        dr <- dir.create(outputDirectory)
        print (sprintf ('creating the output directory named %s; return is %s',
                        outputDirectory, dr))
      }
    }
    ## beware: overwrites without warning!!
    Z <- file.copy (fname25, fnew, overwrite=TRUE)  ## BEWARE: overwrites without warning!!
    
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
  } else {  ## The LR segment begins here:
    Cradeg <- pi/180
    OmegaE <- StandardConstant ('Omega')  ## Earth's rotation rate
    OmegaE <- 15*Cradeg/3600              ## better match to INS?
    Ree <- 6378137                        ## for radii of curvature
    Ecc <- 0.08181919
    
    source ('chunks/RotationCorrection.R')
    source ('chunks/STMFF.R')
    
    ## adjust GPS velocity components for GPS antenna location
    FI_KF <- DataFileInfo (fname, LLrange=FALSE)
    FI_KF$Platform <- getPlatform(fname)  ## also used in chunks/AdjustCal.R
    LG <- ifelse (grepl('130', FI_KF$Platform), -9.88, -4.30)
    MaxGap <- 1000
    .span <- 25    
    
    
    ## These very small adjustments prevent gradual ramping during the flight.
    # D1$BROLLR <- D1$BROLLR + 0.00026
    # D1$BPITCHR <- D1$BPITCHR + 0.00026
    ## ----new-data, include=TRUE, cache=FALSE---------------------------------
    
    SaveRData2 <- sprintf("%s2.Rdata", thisFileName)
    SaveRData <- sprintf("%s.Rdata", thisFileName)
    
    print ('reading netCDF file')
    source ('chunks/AcquireData.R')
    ## add ROC variable
    source ('chunks/ROC.R')
    
    # VV <- c('BLONGA', 'BLATA', 'BNORMA')
    # .shift <- c(-50,-50,-50)
    # names(.shift) <- VV
    # for (V in VV) {
    #   Data[, V] <- ShiftInTime (Data[, V], Rate, .shift[V])
    # }
    # # Data[, 'VSPD'] <- ShiftInTime (Data[, 'VSPD'], Rate, 40)
    # # Data[, 'VEW'] <- ShiftInTime (Data[, 'VEW'], Rate, 60)
    # # Data[, 'VNS'] <- ShiftInTime (Data[, 'VNS'], Rate, 60)
    # Data[, 'PITCH'] <- ShiftInTime (Data[, 'PITCH'], Rate, 20)
    # Data[, 'ROLL'] <- ShiftInTime (Data[, 'ROLL'], Rate, 20)
    .shift = 60
    
    D1 <- Data  ## make adjustments to a copy; avoid changing original
    ## adjustments:
    source ('chunks/AdjustCal.R')
    
    ## GPS shifted vs INS in AdjustCal.R
    
    ## interpolate and protect against missing values
    VV <- c('LAT', 'LON', 'ZROC', 'VEW', 'VNS', 'ROC', 'PITCH', 'ROLL', 'THDG',
            'BPITCHR', 'BROLLR', 'BYAWR')
    DSave <- D1
    for (V in VV) {
      D1[, V] <- zoo::na.approx (as.vector (D1[, V]), maxgap=1000*Rate, na.rm=FALSE, rule=2)
      D1[is.na(D1[, V]), V] <- 0
    }
    D1 <- transferAttributes(DSave, D1)
    attr(D1, 'Platform') <- FI_KF$Platform  ## used in chunks/Kalman-setup.R
    
    ## transform to the a-frame for comparison to the IRU:
    VL <- matrix(c(D1$VEW, D1$VNS, D1$VSPD), ncol=3) 
    LA <- matrix (c(D1$vedot, D1$vndot, -D1$vudot - D1$Grav), ncol=3) + RotationCorrection (D1, VL)
    AA <- XformLA (D1, LA, .inverse=TRUE)
    AA[,3] <- AA[,3] - D1$Grav
    fa1 <- lm(D1$BLONGA ~ AA[, 1])
    fa2 <- lm(D1$BLATA ~ AA[, 2])
    fa3 <- lm(D1$BNORMA ~ AA[, 3])  
    fm1 <- lm (D1$vedot ~ D1$LACCX)
    fm2 <- lm (D1$vndot ~ D1$LACCY)
    fm3 <- lm (D1$vudot ~ D1$LACCZ)
    cfa1 <- coef(fa1); cfa2 <- coef(fa2); cfa3 <- coef(fa3)
    # save(cfa1, cfa2, cfa3, file="./BodyAccCal.Rdata")  ## Don't over-write; data from new file.
    
    ## ----Kalman-setup, include=TRUE, cache=CACHE-----------------------------
    
    ## the chunk is 'sourced' here so the same code can be used in KalmanFilter.R
    # source ('chunks/Kalman-setup.R')
    ## Kalman-setup
    ## initialize matrices needed by the Kalman filter and load the starting-point
    ## for the error-state vector.
    
    if (!SimpleOnly) {
      source ('chunks/Kalman-setup.R')
    } 
    PC <- CorrectPitch(D1, .span=901)
    D1$PC <- PC[, 1]
    D1$RC <- PC[, 2]
    D1$HC <- CorrectHeading(D1, .plotfile='KFplots/HCPlot.pdf')
    
    ## ----Kalman-loop, include=TRUE, cache=CACHE------------------------------
    
    if (!SimpleOnly) {
      print (sprintf ('main loop is starting -- %s', Sys.time()))
      source ('chunks/Kalman-loop.R')
      print (sprintf ('main loop is done -- %s', Sys.time()))
    } 
    
    ## ----plot-Kalman-position, include=TRUE, 
    
    ## Avoid some bad measurements at the start and end
    ## (e.g., from two-pass filtering of LAT/LON)
    minT <- D1$Time[1]- as.integer (D1$Time[1]) %% 60 + 120
    maxT <- D1$Time[nrow(D1)]- as.integer (D1$Time[nrow(D1)]) %% 60 - 120
    r1 <- max(which(D1$Time == minT)[1], which(D1$TASX > 90)[1])
    r2 <- which(D1$TASX > 90); r2 <- r2[length(r2)]
    r2 <- min(which(D1$Time == maxT)[1], r2)
    r <- r1:r2
    DP <- D1[r, ]
    DP <- transferAttributes(D1, DP)
    if (GeneratePlots && !SimpleOnly) {
      print ("generating plots")
      png (file='KFplots/Position.png', width=900, height=600, res=150)
      print (ggplotWAC (with (DP,
                              data.frame(Time, DLONM, CLONM, DLATM, CLATM, DALT, CALT)),
                        ylab=expression (paste ('difference in position [km]')),
                        panels=3, 
                        labelL=c('KF-GPS', 'correction'),
                        labelP=c('east', 'north', 'up'),
                        legend.position=c(0.2,0.68), theme.version=1,
                        gtitle=sprintf ('%s Flight %s', Project, flight)))
      invisible(dev.off())
      if (ShowPlots) {
        print (ggplotWAC (with (DP,
                                data.frame(Time, DLONM, CLONM, DLATM, CLATM, DALT, CALT)),
                          ylab=expression (paste ('difference in position [km]')),
                          panels=3, 
                          labelL=c('KF-GPS', 'correction'),
                          labelP=c('east', 'north', 'up'),
                          legend.position=c(0.2,0.68), theme.version=1))
      }
      pct <- as.integer(100/8)
      print (sprintf ('figures %d%% done', pct))
      sdZ <- sd(DP$DALT, na.rm=TRUE)
      sdvew <- sd(DP$GGVEW - DP$VEWKF, na.rm=TRUE)
      sdvns <- sd(DP$GGVNS - DP$VNSKF, na.rm=TRUE)
      meanvz <- mean(DP$ROCKF - DP[, VROC], na.rm=TRUE)
      sdvz <- sd(DP[, VROC] - DP$ROCKF, na.rm=TRUE)
      
      
      ## ----plot-Kalman-velocity
      
      png(file='KFplots/Velocity.png', width=900, height=600, res=150)
      print (ggplotWAC (with (DP,
                              data.frame(Time, DVEW, CVEW, DVNS, CVNS, DROC, CROC)),
                        ylab=expression (paste ('velocity component [m ',s^-1,']')),
                        panels=3, 
                        labelL=c('KF-GPS', 'correction'),
                        labelP=c('east', 'north', 'up'),
                        legend.position=c(0.2, 0.97), theme.version=1,
                        gtitle=sprintf ('%s Flight %s', Project, flight)))
      invisible(dev.off())
      if (ShowPlots) {
        print (ggplotWAC (with (DP,
                                data.frame(Time, DVEW, CVEW, DVNS, CVNS, DROC, CROC)),
                          ylab=expression (paste ('velocity component [m ',s^-1,']')),
                          panels=3, 
                          labelL=c('KF-GPS', 'correction'),
                          labelP=c('east', 'north', 'up'),
                          legend.position=c(0.2, 0.97), theme.version=1))
      }
      pct <- as.integer(100*2/8)
      print (sprintf ('figures %d%% done', pct))
      
      
      ## ----errors-in-pitch,
      ## must construct using many of the elements of ggplotWAC but in this order:
      attr(DP$CPLF, 'filtered') <- TRUE
      attr(DP$CRLF, 'filtered') <- TRUE
      attributes(DP$SDCPL) <- NULL
      attributes(DP$SDCRL) <- NULL
      d1 <- with(DP, data.frame(Time, CPL, CPLF, CRL, CRLF))
      lines_per_panel <- 2; panels <- 2
      labelL=c('KF value', 'smoothed')
      labelP=c('pitch', 'roll')
      DLP <- nrow(d1)
      VarGroup <- rep (gl (lines_per_panel, DLP, labels=labelL), panels)
      PanelGroup <- gl (panels, lines_per_panel*DLP, labels=labelP)
      dd <- data.frame(reshape2::melt(d1, 1), VarGroup, PanelGroup)
      lvl <- levels(dd$VarGroup)
      d2 <- with(DP, data.frame(Time, "plo"=CPLF-2*SDCPL, "rlo"=CRLF-2*SDCRL,
                                "phi"=CPLF+2*SDCPL, "rhi"=CRLF+2*SDCRL))
      attributes(d2$plo) <- NULL  ## (Don't know why this makes the ribbon-plot work...)
      ## note the required order below:
      de <- data.frame (reshape2::melt(d2, 1, c(2,4,3,5), factorsAsStrings=TRUE), VarGroup, PanelGroup)
      g7 <- ggplot(dd, aes(x=Time)) + facet_grid(PanelGroup ~ .) 
      g7 <- g7 + geom_ribbon(aes(x=Time, ymin=value, ymax=value), data=de, alpha=0.25)  ##alpha=.25 
      g7 <- g7 + geom_path(aes(x=Time, y=value, colour=VarGroup, linetype=VarGroup, linewidth=VarGroup), data=dd)
      g7 <- g7 + ylim(-0.05,0.05)
      g7 <- g7 + scale_linetype_manual ('', labels=lvl, breaks=lvl, values = c(1,1))
      g7 <- g7 + scale_colour_manual('', labels = lvl, breaks=lvl, values = c('blue', 'red'))
      g7 <- g7 + scale_linewidth_manual('', labels=lvl, breaks=lvl, values = c(0.5, 1.5))
      g7 <- g7 + theme_WAC(1) + theme (legend.position=c(0.5,0.95))
      g7 <- g7 + theme(axis.text.x = element_text (size=11.5, margin=margin(15,0,0,0)))
      g7 <- g7 + theme(axis.title.x = element_text (size=12))
      g7 <- g7 + labs (x='Time [UTC]', y=expression (paste ('l-frame error [',degree,']')))
      g7 <- g7 + ggtitle (sprintf ('%s Flight %s', Project, flight))
      ## I'm not sure why this is necessary; printing g7 in the usual way did not leave
      ## the ribbon visible, apparently because the resulting pdf file did not include it.
      ## Now the plot is generated here, but placed in the document in the LaTeX code
      ## preceding this chunk. This out-of-order sequence requires run
      # png (filename='Fig7.png', width=600, height=480, res=150)
      png(file='KFplots/AAlframe.png', width=900, height=600, res=150)
      print (g7)
      invisible(dev.off())
      if (ShowPlots) {print (g7)}
      # invisible (dev.off())
      pct <- as.integer(100*3/8)
      print (sprintf ('figures %d%% done', pct))
      
      
      ## ----a-frame-errors
      
      attr(DP$CRAF, 'filtered') <- TRUE
      attr(DP$CPAF, 'filtered') <- TRUE
      d1 <- with(DP, data.frame(Time, CPAF, PC, CRAF, RC))
      lines_per_panel <- 2; panels <- 2
      labelL=c('smoothed KF value', 'PC')
      labelP=c('pitch', 'roll')
      DLP <- nrow(d1)
      VarGroup <- rep (gl (lines_per_panel, DLP, labels=labelL), panels)
      PanelGroup <- gl (panels, lines_per_panel*DLP, labels=labelP)
      dd <- data.frame(reshape2::melt(d1, 1), VarGroup, PanelGroup)
      lvl <- levels(dd$VarGroup)
      d2 <- with(DP, data.frame(Time, "plo"=CPAF-2*SDCPA, "rlo"=CRAF-2*SDCRA,
                                "phi"=CPAF+2*SDCPA, "rhi"=CRAF+2*SDCRA))
      attributes(d2$plo) <- NULL
      ## note the required order below:
      de <- data.frame (reshape2::melt(d2, 1, c(2,4,3,5)), VarGroup, PanelGroup)
      g9 <- ggplot(dd, aes(x=Time)) + facet_grid(PanelGroup ~ .) 
      g9 <- g9 + geom_ribbon(aes(x=Time, ymin=value, ymax=value), data=de, alpha=0.25) 
      g9 <- g9 + geom_line(aes(x=Time, y=value, colour=VarGroup, linetype=VarGroup), data=dd)
      g9 <- g9 + ylim(-0.05,0.05)
      g9 <- g9 + scale_linetype_manual ('', labels=lvl, breaks=lvl, values = c(1,1))
      g9 <- g9 + scale_colour_manual('', labels = lvl, breaks=lvl, values = c('blue', 'red'))
      g9 <- g9 + theme_WAC(1) + theme (legend.position=c(0.5,0.95))
      g9 <- g9 + theme(axis.text.x = element_text (size=11.5, margin=margin(15,0,0,0)))
      g9 <- g9 + theme(axis.title.x = element_text (size=12))
      g9 <- g9 + labs (x='Time [UTC]', y=expression (paste ('a-frame error [',degree,']')))
      g9 <- g9 + ggtitle (sprintf ('%s Flight %s', Project, flight))
      ## I'm not sure why this is necessary; printing g9 in the usual way did not leave
      ## the ribbon visible, apparently because the resulting pdf file did not include it.
      ## Now the plot is generated here, but placed in the document in the LaTeX code
      ## preceding this chunk. This out-of-order sequence requires run
      # png (filename='Fig9.png', width=600, height=480, res=150)
      png(file='KFplots/AAaframe.png', width=900, height=600, res=150)
      print (g9)
      invisible(dev.off())
      if (ShowPlots) {print (g9)}
      # invisible (dev.off())
      pct <- as.integer(100*4/8)
      print (sprintf ('figures %d%% done', pct))
      
      pct <- as.integer(100*5/8)
      print (sprintf ('figures %d%% done', pct))
    }
    
    if (SimpleOnly) { ## Skip now because these were calculated earlier
      ## these are the errors, negative of the corrections
      # PC <- CorrectPitch(DP, .span=901)
      # DP$PC <- PC[, 1]
      # DP$RC <- PC[, 2]
      # DP$HC <- CorrectHeading (DP, .plotfile='KFplots/HCPlot.pdf')
    } else {
      DP$dP <- DP$deltaPsi/Cradeg
      HA <- with(DP, sqrt(LACCX^2+LACCY^2))
      DP$dP[HA < 1] <- NA
      sddP <- sd(DP$dP, na.rm=TRUE)
    }
    
    ## ----plot-Kalman-heading
    
    # plotWAC(subset(DP,,c(Time, dP)))
    if (GeneratePlots && !SimpleOnly) {
      png(file='KFplots/HDG.png', width=900, height=600, res=150)
      grid.newpage()
      suppressWarnings (
        with(DP,
             ggplotWAC (data.frame(Time, dP),
                        ylim=c(-1.0, 1.0),
                        ylab=expression(paste ('error in heading [',degree,']')),
                        legend.position=NA, position=c(2,2), theme.version=1,
                        gtitle=sprintf ('%s Flight %s', Project, flight))
        )
      )
      DP$dP <- DP$deltaPsi/Cradeg
      HA <- with(DP, sqrt(LACCX^2+LACCY^2))
      DP$dP[HA < 1] <- NA
      DP$SDH <- sqrt(VCor[r,9]) / (Cradeg * 30)
      ## minus sign to change back to error from correction:
      DP$CCTHDG <- -SmoothInterp (DP$CTHDG, .Length=181)
      DP$CCCTHDG <- DP$CCTHDG
      DP$CCCTHDG[DP$SDH > 0.02] <- NA
      suppressWarnings (
        with (DP,
              ggplotWAC(data.frame(Time, 'Kalman'=CCTHDG, 'Ranadu'=HC, 'spline'=-HCS, 'KF'=CCCTHDG), 
                        ylim=c(-1.0, 1.0), lwd=c(0.7,1,1,1.5), lty=c(3,1,1,1),
                        ylab=expression (paste ('error in heading [',degree,']')),
                        col=c('blue', 'forestgreen', 'red', 'blue'), position=c(1,2),
                        legend.position=c(0.5,0.95), theme.version=1
              )
        )
      )
      invisible(dev.off())
      DP$dP <- DP$deltaPsi/Cradeg
      HA <- with(DP, sqrt(LACCX^2+LACCY^2))
      DP$dP[HA < 1] <- NA
      sddP <- sd(DP$dP, na.rm=TRUE)
      if (ShowPlots) {
        grid.newpage()
        suppressWarnings (
          with(DP,
               ggplotWAC (data.frame(Time, dP),
                          ylim=c(-0.4,0.4),
                          ylab=expression (paste ('error in heading [',degree,']')),
                          legend.position=NA, position=c(2,2), theme.version=1,
                          gtitle=sprintf ('%s Flight %s', Project, flight))
          )
        )
        DP$dP <- DP$deltaPsi/Cradeg
        HA <- with(DP, sqrt(LACCX^2+LACCY^2))
        DP$dP[HA < 1] <- NA
        DP$SDH <- sqrt(VCor[r,9]) / (Cradeg * 30)
        DP$CCTHDG <- SmoothInterp (DP$CTHDG, .Length=181)
        DP$CCCTHDG <- DP$CCTHDG
        DP$CCCTHDG[DP$SDH > 0.02] <- NA
        suppressWarnings (
          with (DP,
                ggplotWAC(data.frame(Time, 'Kalman'=CCTHDG, 'Ranadu'=HC, 'spline'=-HCS, 'KF'=CCCTHDG), 
                          ylim=c(-0.4,0.4), lwd=c(0.7,1,1,1.5), lty=c(3,1,1,1),
                          ylab=expression (paste ('error in heading [',degree,']')),
                          col=c('blue', 'forestgreen', 'red', 'blue'), position=c(1,2),
                          legend.position=c(0.5,0.95), theme.version=1
                )
          )
        )
      }
      pct <- as.integer(100*6/8)
      print (sprintf ('figures %d%% done', pct))
    }
    
    ## ----new-AKRD
    # load (file='~/RStudio/Reprocessing/AKRD-fit-coef.Rdata')
    cffn <- 19.70547
    cff <- 21.481
    cfs <- c(4.525341674, 19.933222011, -0.001960992)
    if (grepl('130', FI_KF$Platform)) { # C-130 WECAN values
      cff <- 10.3123
      cfs <- c(5.6885, 14.0452, -0.00461)
    }
    if (SimpleOnly) {  # Skip because these were calculated previously
      # PC <- CorrectPitch(D1, .span=901)
      # D1$PC <- PC[, 1]
      # D1$RC <- PC[, 2]
      # D1$HC <- CorrectHeading (D1)
    }
    D1$QR <- D1$ADIFR / D1$QCF
    D1$QR[D1$QCF < 20] <- NA
    D1$QR[is.infinite(D1$QR)] <- NA
    CutoffFreq <- 600 * Rate
    D1$QRS <- zoo::na.approx (as.vector(D1$QR), maxgap=1000*Rate, na.rm = FALSE, rule=2)
    D1$QRS[is.na(D1$QRS)] <- 0
    D1$QRS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D1$QRS)
    D1$QRF <-  D1$QR - D1$QRS
    D1$QCFS <- zoo::na.approx (as.vector(D1$QCF), maxgap=1000*Rate, na.rm = FALSE)
    D1$QCFS[is.na(D1$QCFS)] <- 0
    D1$QCFS <- signal::filtfilt (signal::butter (3, 2/CutoffFreq), D1$QCFS)
    D1$AKKF <- cff * D1$QRF + cfs[1] + cfs[2] * D1$QRS + cfs[3] * D1$QCFS
    D1$AKKF[D1$QCF < 10] <- NA
    D1$SSKF <- 0.008 + 22.301 * D1$BDIFR / D1$QCF
    if (grepl('130', FI_KF$Platform)) { # C-130 WECAN values; use with heading offset +0.76
      D1$SSKF <- 0.85 + 12.6582 * D1$BDIFR / D1$QCF
    }
    D1$SSKF[D1$QCF < 10] <- NA
    if (UpdateAKRD) {
      D1$ATTACK <- D1$AKRD <- D1$AKKF
    }
    if (UpdateSSRD) {
      D1$SSLIP <- D1$SSRD <- D1$SSKF
    }
    
    ## ----wind-calculation, include=TRUE, cache=CACHE-------------------------
    
    if (SimpleOnly) {
      ## get the wind variables:
      DataW <- D1
      DataN <- WindProcessor (DataW)
      D1$WICC <- DataN$WIN
      D1$WDCC <- DataN$WDN
      D1$WSCC <- DataN$WSN
      DataW$PITCH <- D1$PITCH - PC[,1]
      DataW$ROLL <- D1$ROLL - PC[,2]
      DataW$THDG <- D1$THDG - HC
      DataW$VEW <- D1$VEWC
      DataW$VNS <- D1$VNSC
      DataW$GGVSPD <- D1$ROC
      DataN <- WindProcessor(DataW, LG=0, CompF=TRUE)    ## suppress GPS lever arm)
      D1$WDSC <- DataN$WDN
      D1$WSSC <- DataN$WSN
      D1$WISC <- DataN$WIN
      if (GeneratePlots) {
        png(file='KFplots/AAaframe.png', width=900, height=600, res=150)
        with (D1[r,], plotWAC (data.frame (Time, PC[r,1], PC[r,2])))
        invisible(dev.off())
        if (ShowPlots) {
          with (D1[r,], plotWAC (data.frame (Time, PC[r,1], PC[r,2])))
        }
      }
    } else {
      source ('chunks/wind-calculation.R')
    }
    
    if (GeneratePlots) {
      if (SimpleOnly) {
        png(file='KFplots/Wind.png', width=900, height=600, res=150)
        with (D1[r, ], plotWAC (data.frame(Time, WICC, WISC, WICC-WISC), ylim=c(-5,5)))
        invisible (dev.off())
        if (ShowPlots) {
          with (D1[r, ], plotWAC (data.frame(Time, WICC, WISC, WICC-WISC), ylim=c(-5,5)))
        }
        png(file='KFplots/Wind2.png', width=900, height=600, res=150)
        with (D1[r, ], plotWAC (data.frame (Time, WSCC-WSSC, WDCC-WDSC), ylim=c(-10,10)))
        invisible (dev.off())
        if (ShowPlots) {
          with (D1[r, ], plotWAC (data.frame (Time, WSCC-WSSC, WDCC-WDSC), ylim=c(-10,10)))
        }
      } else {
        png(file='KFplots/Wind.png', width=900, height=600, res=150)
        print (ggplotWAC(
          with(D1[r, ], 
               data.frame (Time, 'Kalman'=WIKF, 'original'=WICC, "DWI"=WIKF-WICC, 'SKIP'=WIKF*0)),
          col=c('blue', 'red'), lwd=c(1.4,0.8,1,0), lty=c(1,42),
          ylab=expression(paste('vertical wind [m ', s^-1, ']')),
          panels=2,
          labelL=c('Kalman', 'original'),
          labelP=c('variables', 'difference'),
          legend.position=c(0.8,0.94), theme.version=1,
          gtitle=sprintf ('%s Flight %s', Project, flight))
        )
        invisible(dev.off())
        pct <- as.integer(100*7/8)
        print (sprintf ('figures %d%% done', pct))
        ## ----hw-plot, include=TRUE,
        D1$WDDIFF <- D1$WDKF - D1$WDCC
        D1$WDDIFF[is.na(D1$WDDIFF)] <- 0
        D1$WDDIFF[D1$WDDIFF > 180] <- D1$WDDIFF[D1$WDDIFF > 180] - 360
        D1$WDDIFF[D1$WDDIFF < -180] <- D1$WDDIFF[D1$WDDIFF < -180] + 360
        D1$WDDIFF[D1$WDDIFF == 0] <- NA
        png(file='KFplots/Wind2.png', width=900, height=600, res=150)
        print (ggplotWAC(
          with(D1[r, ], 
               data.frame (Time, 'direction'=WDDIFF, 'speed'=WSKF-WSCC)),
          ylab=expression(paste('KF correction [', degree, ' (top) or m ', s^-1, ' (bottom)]')),
          panels=2,
          labelL=c('Kalman correction'),
          labelP=c('            direction', '               speed'),
          legend.position=c(0.8,0.94), theme.version=1,
          gtitle=sprintf ('%s Flight %s', Project, flight)) + ylim(-5,5)
        )
        invisible(dev.off())
      }
      if (ShowPlots) {
        print (ggplotWAC(
          with(D1[r, ], 
               data.frame (Time, 'Kalman'=WIKF, 'original'=WICC, "DWI"=WIKF-WICC, 'SKIP'=WIKF*0)),
          col=c('blue', 'red'), lwd=c(1.4,0.8,1,0), lty=c(1,42),
          ylab=expression(paste('vertical wind [m ', s^-1, ']')),
          panels=2,
          labelL=c('Kalman', 'original'),
          labelP=c('variables', 'difference'),
          legend.position=c(0.8,0.94), theme.version=1)
        )
        print (ggplotWAC(
          with(D1[r, ], 
               data.frame (Time, 'direction'=(WDKF-WDCC + 540) %% 360 - 180, 'speed'=WSKF-WSCC)),
          ylab=expression(paste('KF correction [', degree, ' (top) or m ', s^-1, ' (bottom)]')),
          panels=2,
          labelL=c('Kalman correction'),
          labelP=c('            direction', '               speed'),
          legend.position=c(0.8,0.94), theme.version=1) + ylim(-5,5)
        )
      }
    }
    if (file.exists ('KFplots/HCPlot.pdf')) {system ('convert KFplots/HCPlot.pdf KFplots/HCPlot.png')}
    cmd <- sprintf ('cd KFplots;tar cvfz %sPlots%s.tgz *png', Project, flight)
    system(cmd)
    print (sprintf ('plots generated -- %s', Sys.time()))
    
    ## ----create-new-netcdf, cache=FALSE--------------------------------------
    
    print (sprintf("making new netCDF file -- %s", Sys.time()))
    # source ('chunks/create-new-netcdf.R')
    fnew <- sub ('.nc', 'KF.nc', fname)
    fnew <- paste0(outputDirectory, sub('.*/', '', fnew))
    if(SHINY) {  ## for the datavis shiny server, only this output is allowed:
      fnew <- paste0('KFoutput/', sub('.*/', '', fnew))
    }
    print (sprintf ('making new netCDF file named %s', fnew))
    ## Make the directory KF if needed:
    if (grepl ('/KF/', fnew)) {
      if (!dir.exists(outputDirectory)) {
        dr <- dir.create(outputDirectory)
        print (sprintf ('creating the output directory named %s; return is %s',
                        outputDirectory, dr))
      }
    }
    ## beware: overwrites without warning!!
    Z <- file.copy (fname, fnew, overwrite=TRUE)  ## BEWARE: overwrites without warning!!
    
    ## ----modify-new-netcdf, include=TRUE-------------------------------------
    
    source ('chunks/modify-new-netcdf.R')
  }
}
## modify-new-netcdf.R
print (sprintf ("Kalman Processor is finished -- %s", Sys.time()))

