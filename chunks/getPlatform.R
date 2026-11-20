## getPlatform.R
## Return the aircraft tail number (e.g., "N130AR" or "N677F") from the global
## attributes of a netCDF file. Older RAF files use "Platform"; newer files
## use "platform_identifier". Stops rather than guessing, because the C-130
## and GV use different calibration constants.
getPlatform <- function(fname) {
  nc <- ncdf4::nc_open(fname)
  on.exit(ncdf4::nc_close(nc))
  for (att in c('Platform', 'platform_identifier')) {
    a <- ncdf4::ncatt_get(nc, 0, att)
    if (a$hasatt && nzchar(a$value)) return(a$value)
  }
  stop(sprintf("no 'Platform' or 'platform_identifier' global attribute in %s", fname))
}
