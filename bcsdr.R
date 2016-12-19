# setwd("~/Documents/bcsd-python")
# setwd("D:/works/bcsd-python")
source('functions.R')
library(RNetCDF)
#library(ncdf4)
library(Rmpi)
library(snow)
# ################### Read data with ncdf4 package ################### 
# # school cluster (discovery) cannot install this package
# # 1. read in reanalysis data
# gcm_ncin = nc_open("./data/merra_example.nc") # 1980.1.1 - 2015.12.31 (13149 days)
# gcm_lon = ncvar_get(gcm_ncin,"lon")
# gcm_lat = ncvar_get(gcm_ncin,"lat")
# gcm_time= ncvar_get(gcm_ncin,"time") # ! there are duplicate times, for example X0[,,732]==X0[,,733]
# gcm_var="PRECTOTLAND"
# X0=ncvar_get(gcm_ncin,gcm_var) # 9    11 13263
# 
# # remove duplicate time (the later one)
# Nt=length(gcm_time)
# id_sel=c(TRUE,!diff(gcm_time)==0) # remove duplicate id 
# id_sel=(366<1:Nt)&(1:Nt<=(Nt-365))&id_sel # remove 1980 & 2015 to fit observation data
# gcm_time = gcm_time[id_sel]
# X0=X0[,,id_sel]
# 
# # 2. read in observation data
# obs_ncin=nc_open("./data/prism_example.nc") # 1981.1.1 - 2014.12.31 (12418 days)
# obs_lon = ncvar_get(obs_ncin,"lon")
# obs_lat = ncvar_get(obs_ncin,"lat")
# obs_time= ncvar_get(obs_ncin,"time")
# obs_var="ppt"
# Y0=ncvar_get(obs_ncin,obs_var) # 30    30 12418

################### Read data with RNetCDF package ################### 
# 1. read in reanalysis data

# It seems that RNetCDF cannot read in netcdf4 file(version 4). 
# So here I use python xarray library to resave the data into version 3, with the following commands:
# # http://stackoverflow.com/questions/40645472/how-to-read-netcdf4-data-if-i-only-have-netcdf3-tool/40645667#40645667
# import xarray as xr
# ds=xr.open_dataset("merra_example.nc")
# ds.to_netcdf("merra_example_v3.nc",format='NETCDF3_CLASSIC')

gcm_ncin = open.nc("./data/merra_example_v3.nc")
gcm_lon = var.get.nc(gcm_ncin,"lon")
gcm_lat = var.get.nc(gcm_ncin,"lat")
gcm_time= var.get.nc(gcm_ncin,"time") # ! there are duplicate times, for example X0[,,732]==X0[,,733]
gcm_var="PRECTOTLAND"
X0=var.get.nc(gcm_ncin,gcm_var) # 9    11 13263

# remove duplicate time (the later one)
Nt=length(gcm_time)
id_sel=c(TRUE,!diff(gcm_time)==0) # remove duplicate id 
id_sel=(366<1:Nt)&(1:Nt<=(Nt-365))&id_sel # remove 1980 & 2015 to fit observation data
gcm_time = gcm_time[id_sel]
X0=X0[,,id_sel]

# 2. read in observation data
obs_ncin=open.nc("./data/prism_example.nc") 
obs_lon=var.get.nc(obs_ncin,"lon")
obs_lat=var.get.nc(obs_ncin,"lat")
obs_time=var.get.nc(obs_ncin,"time")
obs_var="ppt"
Y0=var.get.nc(obs_ncin,obs_var) # 30    30 12418

################### BCSD ################### 
# Assuming equal spaced for GCM, get delta_lon and delta_lat to find the grid.
# And perform analysis grid by grid
gcm_delta_lon=gcm_lon[2]-gcm_lon[1] 
gcm_delta_lat=gcm_lat[2]-gcm_lat[1] 
gcm_n_lon=length(gcm_lon)
gcm_n_lat=length(gcm_lat)

Y_bcsd=Y_cbcsd=Y0

for(ilon in 1:(gcm_n_lon-1)){ 
  for(jlat in 2:(gcm_n_lat)){ 
    i=gcm_lon[ilon] # i is current GCM lon
    j=gcm_lat[jlat] # j is current GCM lat

    cat("GCM grid: lon",i,"lat",j,"...\n")
    
    # 0.1 Select data: current GCM grid (GCM point corresponds to the upper left corner)
    X=X0[ilon,jlat,] 
    
    delta_lon = obs_lon - i # select obs > gcm : gcm grid point on the west (left)
    obs_selected_lon=(0<=delta_lon)&(delta_lon < gcm_delta_lon)
    delta_lat = j - obs_lat # select obs < gcm : gcm grid point on the north (top) 
    obs_selected_lat=(0<=delta_lat)&(delta_lat < gcm_delta_lat)
    
    if(sum(obs_selected_lon)*sum(obs_selected_lat)==0){
      cat("No observation in this GCM grid. \n")
      next
    }
    
    # 0.2 Select data: current observation in the GCM grid
    Y=Y0[obs_selected_lon,obs_selected_lat,] 
    Ybar=apply(Y,3,mean, na.rm = TRUE) # upscaled observation
    
    # 1. Bias correction: X and Ybar
    cat("Bias correction ...\n")
    Q=quantile(Ybar,ecdf(X)(X)) #Qi_test=quantile(Yim,ecdf(Xim)(Xim_test))
    # Bias correction (copula): X and Ybar
    cat("Bias correction (copula) ...\n")
    CQ=X
    ntime=length(X)
    # for(k in 1:ntime){
    #   cat("correcting time ",k,"\n")
    #   CQ[k]=max_cond(X[k],X,Ybar)
    # }
    
    cl <- makeCluster(30, type = "MPI") 
    ml=clusterApplyLB(cl, 1:ntime, function(x) {cat("correcting time ",x,"\n");max_cond(X[x],X,Ybar)})
    for(k in 1:ntime){
      CQ[k]=ml[[k]]
    }
    stopCluster(cl)
    
    # 2. Scaling factor 
    cat("Spatial Disaggregation ... \n")
    M=apply(Y,c(1,2),mean, na.rm = TRUE) # Local *Historical* mean
    Y_Q=Y_CQ=Y
    for(k in 1:ntime){
      Y_Q[,,k]=Q[k]*M/mean(Q)
      Y_CQ[,,k]=CQ[k]*M/mean(CQ)
    }
    
    # 3. Results
    Y_bcsd[obs_selected_lon,obs_selected_lat,]=Y_Q
    Y_cbcsd[obs_selected_lon,obs_selected_lat,]=Y_CQ
    
    #mseij=mean((Y_Q-Y)^2)
    #cat("mse:",mseij,"\n")
  }
}

save(Y_bcsd, Y_cbcsd, Y, file="bcsd_results.RData")






