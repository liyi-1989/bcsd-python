#library(ncdf4)
library(MASS)

# transform lon*lat*time to time*space
xyt2ts=function(P){
  nx=dim(P)[1]
  ny=dim(P)[2]
  nt=dim(P)[3]
  Q=matrix(0,nrow=nt,ncol=nx*ny)
  for(k in 1:nt){
    for(i in 1:nx){
      for(j in 1:ny){
        count=(i-1)*ny+j
        Q[k,count]=P[i,j,k]
      }
    }
  }
  return(Q)
}
# transform time*space to lon*lat*time
ts2xyt=function(Q,nx,ny){
  nt=dim(Q)[1]
  if(dim(Q)[2]!=(nx*ny)){
    stop("Dimension mismatch!")
  }
  P=array(0,c(nx,ny,nt))
  for(k in 1:nt){
    P[,,k]=matrix(Q[k,],nrow=nx,ncol=ny,byrow = T)
  }
  return(P)
}

# daily to monthly for P: lon*lat*365
d2m=function(P,method="mean"){
  
  #---------------------------------
  MON=1:365
  names(MON)=c(rep("JAN",31),rep("FEB",28),rep("MAR",31),rep("APR",30),rep("MAY",31),rep("JUN",30),
               rep("JUL",31),rep("AUG",31),rep("SEP",30),rep("OCT",31),rep("NOV",30),rep("DEC",31))
  Mon=data.frame(mon=as.character(unique(names(MON))),num=1:12,days=c(31,28,31,30,31,30,31,31,30,31,30,31))
  Mon["Days"]=cumsum(Mon["days"])
  #----------------------------------
  
  nx=dim(P)[1]
  ny=dim(P)[2]
  nt=dim(P)[3]
  nyear=nt/365
  Q=array(0,c(nx,ny,12*nyear))
  L=rep(names(MON),nyear)
  
  if(method=="mean"){
    for(i in 1:nx){
      for(j in 1:ny){
        for(k in 1:12){
          for(l in 1:nyear){
            count=(l-1)*12+k
            if(k==1){
              s_day=(l-1)*365+1
            }else{
              s_day=(l-1)*365+Mon[k-1,"Days"]+1
            }
            e_day=(l-1)*365+Mon[k,"Days"]
            Q[i,j,count]=mean(P[i,j,s_day:e_day])
          }
        }
      }
    }
  }else if(method=="sum"){
    for(i in 1:nx){
      for(j in 1:ny){
        for(k in 1:12){
          for(l in 1:nyear){
            count=(l-1)*12+k
            if(k==1){
              s_day=(l-1)*365+1
            }else{
              s_day=(l-1)*365+Mon[k-1,"Days"]+1
            }
            e_day=(l-1)*365+Mon[k,"Days"]
            Q[i,j,count]=sum(P[i,j,s_day:e_day])
          }
        }
      }
    }
  }
  
  return(Q)
}


get_gcm_all=function(){
  files=list.files("./gcm_train")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm_train/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    Pr1=d2m(Pr)
    pr=rbind(pr,xyt2ts(Pr1))
  }
  return(pr)
}

# Load the GCM data
get_gcm_train=function(mylon,mylat){
  cat("In position: ",mylon,mylat,"\n")
  files=list.files("./gcm_train")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm_train/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    pr=c(pr,Pr[((mylon-0.1)<lon)&(lon<(mylon+0.1)),((mylat-0.1)<lat)&(lat<(mylat+0.1)),])
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}

get_gcm_test=function(mylon,mylat){
  cat("In position: ",mylon,mylat,"\n")
  files=list.files("./gcm_test")
  pr=NULL
  for(i in files){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("./gcm_test/",i))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"pr")
    pr=c(pr,Pr[((mylon-0.1)<lon)&(lon<(mylon+0.1)),((mylat-0.1)<lat)&(lat<(mylat+0.1)),])
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}
# Load the observation data, unusful, use matlab do this.
get_obs=function(mylon,mylat){
  cat("In position: [",mylon,"]*[",mylat,"]\n")
  pr=NULL
  for(i in 1979:2003){
    cat("Read in file ",i,"\n")
    ncin = nc_open(paste0("/media/liyi/新加卷/downscaling/pr_",i,".nc"))
    lon = ncvar_get(ncin,"lon")
    lat = ncvar_get(ncin,"lat")
    Pr=ncvar_get(ncin,"precipitation_amount") #[day,lon,lat] 
    tmp=Pr[,((mylon[1]-0.01)<lon)&(lon<(mylon[2]+0.01)),((mylat[1]-0.01)<lat)&(lat<(mylat[2]+0.01))]
    tmp=apply(tmp,1,mean)
    if(length(tmp)==366){
      tmp[31+28]=mean(tmp[(31+28):(31+29)])
      tmp=tmp[-(31+29)]
    }
    pr=c(pr,tmp)
    rm(Pr)
  }
  prm=matrix(pr,365,length(pr)/365)
  return(prm)
}


# Y=get_obs(mylon,mylat)

## Density estimation

# plot(rank(Xim)/25,rank(Yim)/25)
# plot(density(Xim))
# x0=2
# density(Xim, from=x0, to=x0, n=1)$x
# density(Xim, from=x0, to=x0, n=1)$y
# 
# f2=kde2d(rank(Xim)/25,rank(Yim)/25)
# image(f2)
# persp(f2, phi = 30, theta = 20, d = 5)
# contour(f2)
# 
# kde2d(rank(Xim)/25,rank(Yim)/25,lims=c(0.5,0.5,0.3,0.3),n=1)



kde1d_point=function(x0,x_given){
  density(x_given, from=x0, to=x0, n=1)$y
}

kde2d_point=function(x0,y0,x_given,y_given){
  kde2d(x_given,y_given,lims=c(x0,x0,y0,y0),n=1)$z
}

max_cond=function(x_fix,x_given,y_given){
  ygrid=seq(from=min(y_given),to=max(y_given),length=500)
  yt=ygrid
  n_given=length(x_given)
  u=rank(x_given)/n_given
  v=rank(y_given)/n_given
  for(i in 1:length(ygrid)){
    #yt[i]=kde2d_point(rank(c(x_fix,x_given))[1]/26,rank(c(ygrid[i],y_given))[1]/26,u,v)*kde1d_point(x_fix,x_given)
    yt[i]=kde2d_point(rank(c(x_fix,x_given))[1]/(n_given+1),rank(c(ygrid[i],y_given))[1]/(n_given+1),u,v)*kde1d_point(ygrid[i],y_given)
  }
  #plot(ygrid,yt,type="o")
  #ygrid[yt==max(yt)]
  ygrid[which.max(yt)]
}



#======================== Covariance Estimation  ==============================
regularization_by_thresholding = function (S, MP=1) {
  p=ncol(S)
  n=nrow(S)
  return(S*(abs(S)>MP*sqrt(log(p)/n)))
}