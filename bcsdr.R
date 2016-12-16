# setwd("~/Documents/bcsd-python")
source('functions.R')

# 1. read in reanalysis data
gcm_ncin = nc_open("./data/merra_example.nc")
gcm_lon = ncvar_get(gcm_ncin,"lon")
gcm_lat = ncvar_get(gcm_ncin,"lat")
gcm_time= ncvar_get(gcm_ncin,"time")
gcm_var="PRECTOTLAND"
X=ncvar_get(gcm_ncin,gcm_var) # 9    11 13263

# 2. read in observation data
obs_ncin=nc_open("./data/prism_example.nc")
obs_lon = ncvar_get(obs_ncin,"lon")
obs_lat = ncvar_get(obs_ncin,"lat")
obs_time= ncvar_get(obs_ncin,"time")
obs_var="ppt"
Y=ncvar_get(obs_ncin,obs_var) # 30    30 12418





