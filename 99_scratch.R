library(ncdf4)

thetao_file <- ncdf4::nc_open(sst_path2)

names(thetao_file$var)
names(thetao_file$dim)

thetao_raw <- ncdf4::ncvar_get(thetao_file, "thetao")
thetao_time <- ncdf4::ncvar_get(thetao_file, "time")

thetao_file$dim$time$units

ncdf4::nc_close(thetao_file)
dim(thetao_raw)
summary(as.numeric(thetao_raw))
