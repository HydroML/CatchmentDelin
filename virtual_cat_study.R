#####################################################################################################################################################3
#### satellite rainfall with virtual catchments track years of data, same/different unit hydrographs, level of spatial cor, max NSE ########3
######################################################################################################################################################
rm(list = ls())
library(ncdf4)
library(ncdf4.helpers)
library(lubridate)
library(infotheo)
library(sf)
library(sp)
library(jtools)

#prepare and save data
shapename <- read_sf('C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/CAMELS_DE/catchments/CAMELS_DE_catchments.shp')
shapename$geometry<-st_transform(shapename$geometry, 4326)
gids<-shapename$gauge_id
precip_cors<-rep(0,length(gids))
precip_rcors<-rep(0,length(gids))
precip_means<-rep(0,length(gids))
precip_numoc<-rep(0,length(gids))
precip_autocor<-rep(0,length(gids))
for(i in 1:length(gids)){
  catchment=gids[i]
  
  x_repr = shapename$geometry[which(shapename$gauge_id==catchment)]
  x_repr= x_repr[[1]][[1]]
  #plot(x_repr,type="l")
  
  maxlat<-ceiling((max(x_repr[,2])+0.5)*10)/10
  maxlon<- ceiling((max(x_repr[,1])+0.5)*10)/10
  minlon<- floor((min(x_repr[,1])-0.5)*10)/10
  minlat<- floor((min(x_repr[,2])-0.5)*10)/10
  
  
  importyears<-1983:2020
  currentWeatherdata<-nc_open(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/CAMELS_DE/precip/persiann_cdr/CDR_Germany_2025-03-19085141pm_",importyears[1],".nc"))
  londata<-ncvar_get(currentWeatherdata,varid = "lon")
  latdata<-ncvar_get(currentWeatherdata,varid = "lat")
  
  whichlon<-which(londata>= minlon & londata<=maxlon)
  whichlat<-which(latdata>=minlat & latdata<=maxlat)
  lonvec<-londata[whichlon]
  latvec<-latdata[whichlat]
  posinfo<-expand.grid(latvec,lonvec)
  colnames(posinfo)<-c("lat","lon")
  posinfo$whichpos_lat<-NA
  posinfo$whichpos_lon<-NA
  for(j in 1:nrow(posinfo)){
    posinfo$whichpos_lat[j]<-which(latdata==posinfo$lat[j])
    posinfo$whichpos_lon[j]<-which(londata==posinfo$lon[j])
  }
  nc_close(currentWeatherdata)
  
  p_dates<-seq(from=as.Date("1983-01-01"),to=as.Date("2020-12-31"),by="day")
  precip_mat<-matrix(NA,nrow = nrow(posinfo),ncol = length(p_dates))
  curday<-1
  for(k in 1:length(importyears)){
    currentWeatherdata<-nc_open(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/CAMELS_DE/precip/persiann_cdr/CDR_Germany_2025-03-19085141pm_",importyears[k],".nc"))
    londata<-ncvar_get(currentWeatherdata,varid = "lon")
    latdata<-ncvar_get(currentWeatherdata,varid = "lat")
    curdat<-ncvar_get(currentWeatherdata,varid = "precip")
    nc_close(currentWeatherdata)
    for(j in 1:nrow(precip_mat)){
      precip_mat[j,curday:(curday+dim(curdat)[3]-1)]<-curdat[posinfo$whichpos_lon[j],posinfo$whichpos_lat[j],]
    }
    curday<-curday+dim(curdat)[3]
  }
  
  precip_mat[precip_mat>175]<-0
  precip_mat[precip_mat<0]<-0
  
  goodpoints<-which(rowSums(precip_mat)>0)
  precip_mat<-precip_mat[goodpoints,]
  posinfo<-posinfo[goodpoints,]
  
  
  posinfo$incat<-NA
  for(j in 1:nrow(posinfo)){
    posinfo$incat[j]<- point.in.polygon(posinfo$lon[j], posinfo$lat[j], x_repr[,1], x_repr[,2])
  }
  if(sum(posinfo$incat==0)) posinfo$incat[which.min((posinfo$lon -mean(x_repr[,1]))^2 + (posinfo$lat - mean(x_repr[,2]))^2)]<-1
  #plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[posinfo$incat+1],cex=6,
  #     ylim=c(min(posinfo$lat)-0.2,max(posinfo$lat)+0.2), xlim=c(min(posinfo$lon)-0.05,max(posinfo$lon)+0.05))
  #lines(x_repr,col="black",lwd=3)
  precip_means[i]<-mean(precip_mat[posinfo$incat==1,])
  precip_numoc[i]<-mean(rowMeans(matrix(precip_mat[posinfo$incat==1,],nrow = sum(posinfo$incat==1))>0))
  
  #detect if grid has all neighbors
  posinfo$num_neigh<-0
  posinfo$cor_neigh<-0
  posinfo$rcor_neigh<-0
  posinfo$auto_cor<-0
  for(j in 1:nrow(posinfo)){
    curpos_lat<-posinfo$whichpos_lat[j]
    curpos_lon<-posinfo$whichpos_lon[j]
    posinfo$auto_cor[j]<-acf(precip_mat[j,],plot = F)$acf[2]
    for(k in setdiff(1:nrow(posinfo),j)){
      if(curpos_lat==posinfo$whichpos_lat[k]+1 & curpos_lon==posinfo$whichpos_lon[k]){
        posinfo$num_neigh[j]<-posinfo$num_neigh[j]+1
        posinfo$cor_neigh[j]<-posinfo$cor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,])
        #posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,],method="spearman")
        posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ mutinformation(discretize(precip_mat[j,]),discretize(precip_mat[k,]))
      } else if(curpos_lat==posinfo$whichpos_lat[k]-1 & curpos_lon==posinfo$whichpos_lon[k]){
        posinfo$num_neigh[j]<-posinfo$num_neigh[j]+1
        posinfo$cor_neigh[j]<-posinfo$cor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,])
        #posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,],method="spearman")
        posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ mutinformation(discretize(precip_mat[j,]),discretize(precip_mat[k,]))
      } else if(curpos_lat==posinfo$whichpos_lat[k] & curpos_lon==posinfo$whichpos_lon[k]+1){
        posinfo$num_neigh[j]<-posinfo$num_neigh[j]+1
        posinfo$cor_neigh[j]<-posinfo$cor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,])
        #posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,],method="spearman")
        posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ mutinformation(discretize(precip_mat[j,]),discretize(precip_mat[k,]))
      } else if(curpos_lat==posinfo$whichpos_lat[k] & curpos_lon==posinfo$whichpos_lon[k]-1){
        posinfo$num_neigh[j]<-posinfo$num_neigh[j]+1
        posinfo$cor_neigh[j]<-posinfo$cor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,])
        #posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ cor(precip_mat[j,],precip_mat[k,],method="spearman")
        posinfo$rcor_neigh[j]<-posinfo$rcor_neigh[j]+ mutinformation(discretize(precip_mat[j,]),discretize(precip_mat[k,]))
      } 
    }
  }
  
  saveRDS(list(posinfo=posinfo,precip_mat=precip_mat), 
          file = paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",catchment,".rds"))
  
  precip_cors[i]<-mean(posinfo$cor_neigh[posinfo$num_neigh>0]/posinfo$num_neigh[posinfo$num_neigh>0])
  precip_rcors[i]<-mean(posinfo$rcor_neigh[posinfo$num_neigh>0]/posinfo$num_neigh[posinfo$num_neigh>0])
  precip_autocor[i]<-mean(posinfo$auto_cor)
  print(i)
}


shapename$cor<-precip_cors
shapename$rcor<-precip_rcors
shapename$pmean<-precip_means
shapename$poc<-precip_numoc
shapename$autocor<-precip_autocor


saveRDS(shapename, 
        file = paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/shpfile_withp.rds"))

shapename<-readRDS("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/shpfile_withp.rds")

plot(shapename["rcor"],axes = TRUE,main=NULL)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/precip_mi.jpeg",width = 5, height = 7, units = 'in',res = 400)
plot(shapename["rcor"],axes = TRUE,pal=hcl.colors(13,palette = "Blue-Yellow",rev=T),main=NULL)
dev.off()

plot(shapename["cor"],axes = TRUE,main=NULL)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/precip_cor.jpeg",width = 5, height = 7, units = 'in',res = 400)
plot(shapename["cor"],axes = TRUE,pal=hcl.colors(7,palette = "Blue-Yellow",rev=T),main=NULL)
dev.off()

plot(shapename["pmean"],axes = TRUE,main=NULL)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/pmean.jpeg",width = 5, height = 7, units = 'in',res = 400)
plot(shapename["pmean"],axes = TRUE,pal=hcl.colors(15,palette = "Blue-Yellow",rev=T),main=NULL)
dev.off()

plot(shapename["poc"],axes = TRUE,main=NULL)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/poc.jpeg",width = 5, height = 7, units = 'in',res = 400)
plot(shapename["poc"],axes = TRUE,pal=hcl.colors(15,palette = "Blue-Yellow",rev=T),main=NULL)
dev.off()

plot(shapename["autocor"],axes = TRUE,main=NULL)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/autocor.jpeg",width = 5, height = 7, units = 'in',res = 400)
plot(shapename["autocor"],axes = TRUE,pal=hcl.colors(13,palette = "Blue-Yellow",rev=T),main=NULL)
dev.off()



####################################################################################################################################################
#############################################   run the different methods on virtual catchmnet data   ##########################
##################################################################################################################################################
rm(list = ls())

get_cat_delin<-function(p,q,md){
  if(md=="mi"){
    return(get_mi(p=p,q=q))
  } else if(md=="te"){
    return(get_te(p=p,q=q))
  } else if(md=="cause"){
    return(get_cause(p=p,q=q))
  } else{
    return(-99)
  }
}

get_mi<-function(p,q){
  in_cat_pred<-rep(1,nrow(p))
  nlag<-3
  Y_t=q[(1+nlag):length(q)]
  Y_t= discretize(Y_t)
  fulldat<-matrix(0,nrow = nrow(Y_t),ncol = nlag+1)
  for(gridsquare in 1:nrow(p)){
    for(lag in 0:nlag){
      fulldat[,lag+1]<- p[gridsquare,(nlag- lag+1):(length(q)-lag)]
    }
    in_cat_pred[gridsquare]<-mutinformation(discretize(fulldat),Y_t)
  }
  
  bs_mi_res<-rep(0,20)
  for(b in 1:length(bs_mi_res)){
    for(lag in 0:nlag){
      fulldat[,lag+1]<-p[sample(1:nrow(p),1),sample(1:ncol(p),nrow(Y_t),replace = T)]
    }
    bs_mi_res[b]<-mutinformation(discretize(fulldat),Y_t)
  }
  thresh<-mean(bs_mi_res)+ 3*sd(bs_mi_res)
  in_cat_pred<-as.numeric(in_cat_pred > thresh)
  return(in_cat_pred)
}


get_te<-function(p,q){
  in_cat_pred<-sample(c(0,1),nrow(p),replace = T)
  nlag<-3
  Y_t=q[(1+nlag):length(q)]
  Y_t= discretize(Y_t)
  fulldat1<-data.frame(prev_q= q[(nlag):(length(q)-1)])
  partdat_mi<-mutinformation(discretize(fulldat1),Y_t)
  for(gridsquare in 1:nrow(p)){
    fulldat<-fulldat1
    for(lag in 0:nlag){
      fulldat<- cbind(fulldat,p[gridsquare,(nlag- lag+1):(length(q)-lag)])
    }
    fulldat<-discretize(fulldat)
    fulldat_mi<-mutinformation(fulldat,Y_t)
    in_cat_pred[gridsquare]<-fulldat_mi- partdat_mi
  }
  
  bs_mi_res<-rep(0,20)
  for(b in 1:length(bs_mi_res)){
    fulldat<-fulldat1
    for(lag in 0:nlag){
      fulldat<- cbind(fulldat,p[sample(1:nrow(p),1),sample(1:ncol(p),nrow(Y_t),replace = T)])
    }
    fulldat<-discretize(fulldat)
    fulldat_mi<-mutinformation(fulldat,Y_t)
    bs_mi_res[b]<-fulldat_mi- partdat_mi
  }
  thresh<-mean(bs_mi_res)+ 3*sd(bs_mi_res)
  in_cat_pred<-as.numeric(in_cat_pred > thresh)
  return(in_cat_pred)
}




get_cause<-function(p,q){
  nlag<-5
  Y_t=q[(nlag):length(q)]
  fulldat<-as.data.frame(matrix(0,nrow = length(Y_t),ncol = nlag*nrow(p)+1))
  fulldat$V1<-Y_t
  counter<-2
  for(i in 1:nrow(p)){
    for(j in 1:nlag){
      fulldat[,counter]<- p[i,(nlag- j+1):(ncol(p)-j+1)]
      counter=counter+1
    }
  }
  fit<-lm(V1~.,data = fulldat)
  in_cat_pred<-coef(fit)[-1]
  thresh=as.numeric(3*summary(fit)$coefficients[-1, 2])
  

  in_cat_pred<-as.numeric(in_cat_pred > thresh)[seq(from=2,by=nlag,length.out=nrow(p))]
  
  
  return(in_cat_pred)
}



shapename<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/shpfile_withp.rds"))

gids<-shapename$gauge_id
gids<-gids[1:1500]

nreps<-1
maxR2s<-c(0.6,0.75,0.9)
yearsOfData<-c(10,24,38)
signaldiff<-1 #seq(from=0,to=1,length.out=3)
meth<-c("mi","te","cause")



missclass_array<-array(NA,dim=c(length(gids),nreps,length(signaldiff), length(maxR2s), length(yearsOfData),length(meth)))
fp_array<-array(NA,dim=c(length(gids),nreps,length(signaldiff), length(maxR2s), length(yearsOfData),length(meth)))
fn_array<-array(NA,dim=c(length(gids),nreps,length(signaldiff), length(maxR2s), length(yearsOfData),length(meth)))


for(h in 1:length(gids)){
  curcat<-gids[h]
  p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
  posinfo<-p_data$posinfo
  precip_mat<-p_data$precip_mat
  tot_incat<-sum(posinfo$incat)
  
  for(i in 1:nreps){
    
    for(j in 1:length(signaldiff)){
      cursignaldiff<-signaldiff[j]
      alphadiffs<-runif(tot_incat,min = 0,max=cursignaldiff)
      thetadiffs<-runif(tot_incat,min = 0,max=cursignaldiff)
      alphas<-rep(1,tot_incat)+alphadiffs
      thetas<-rep(2,tot_incat)-thetadiffs
      #hist(rgamma(10000,shape = alphas[1],scale = thetas[1]),breaks=100,xlim = c(0,15))
      gt_streamflow<-rep(0, ncol(precip_mat))
      for(grid in 1:tot_incat){
        uh<- dgamma(0:10 +0.5,shape = alphas[grid],scale = thetas[grid])
        uh<-uh/sum(uh)
        what<-convolve(precip_mat[which(posinfo$incat==1)[grid],],uh,type="open",conj = F)
        gt_streamflow<-gt_streamflow + what[(length(uh)):length(what)]
      }
      
      for(k in 1:length(maxR2s)){
        curmaxR2<-maxR2s[k]
        noiseVar<-(var(gt_streamflow)-curmaxR2*var(gt_streamflow))/curmaxR2
        noise<-as.numeric(arima.sim(list(order=c(1,0,0),ar=0.01),n=length(gt_streamflow),sd=1))
        noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
        streamflow=gt_streamflow+noise
        
        
        for(l in 1:length(yearsOfData)){
          curyears<-yearsOfData[l]
          useable_days<- length(seq(from=as.Date("1983-01-01"),to=as.Date(paste0(1983+curyears-1,"-12-31")),by="day"))
          useable_p<-precip_mat[,1:useable_days]
          useable_streamflow<-streamflow[1:useable_days]
          
          for(m in 1:length(meth)){
            curmeth<-meth[m]
            results<-get_cat_delin(p=useable_p,q=useable_streamflow,md=curmeth)
            missclass_array[h,i,j,k,l,m]<-mean(results!= posinfo$incat)
            fp_array[h,i,j,k,l,m]<-mean(results[results==1]!= posinfo$incat[results==1])
            fn_array[h,i,j,k,l,m]<-mean(results[results==0]!= posinfo$incat[results==0])
          }
        }
      }
    }
  }
  print(h)
}

saveRDS(list(missclass=missclass_array,fp=fp_array,fn=fn_array), 
        file = paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/results_list.rds"))








what=readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/results_list.rds"))
res_to_plot<-what$missclass[1:1500,,,,,3,drop = F]
n<-length(res_to_plot)
rest_to_plot_dat<-data.frame(missclass=rep(0,n),catchment=rep(0,n),rep=rep(0,n),signaldiff=rep(0,n),maxR2=rep(0,n),nyears=rep(0,n),sc=rep(0,n),npix=rep(0,n))
index<-1
for(i in 1:dim(res_to_plot)[1]){
  curcat<-gids[i]
  p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
  posinfo<-p_data$posinfo
  precip_mat<-p_data$precip_mat
  
  for(j in 1:dim(res_to_plot)[2]){
    for(k in 1:dim(res_to_plot)[3]){
      for(l in 1:dim(res_to_plot)[4]){
        for(m in 1:dim(res_to_plot)[5]){
          rest_to_plot_dat$catchment[index]<-i
          rest_to_plot_dat$rep[index]<-j
          rest_to_plot_dat$signaldiff[index]<-c("Low","Medium","High")[k]
          rest_to_plot_dat$maxR2[index]<-c("Low","Medium","High")[l]
          rest_to_plot_dat$nyears[index]<-c("Low","Medium","High")[m]
          rest_to_plot_dat$sc[index]<-shapename$cor[i]
          rest_to_plot_dat$missclass[index]<-res_to_plot[i,j,k,l,m,1]
          rest_to_plot_dat$npix[index]<-nrow(posinfo)
          index<-index+1
        }
      }
    }
  }
  print(i)
}

rest_to_plot_dat$signaldiff <- factor(rest_to_plot_dat$signaldiff , levels=c("Low","Medium","High"))
rest_to_plot_dat$maxR2 <- factor(rest_to_plot_dat$maxR2 , levels=c("Low","Medium","High"))
rest_to_plot_dat$nyears <- factor(rest_to_plot_dat$nyears , levels=c("Low","Medium","High"))
rest_to_plot_dat$sc_class<-"Low"
rest_to_plot_dat$sc_class[rest_to_plot_dat$sc> quantile(rest_to_plot_dat$sc,1/3)]<-"Medium"
rest_to_plot_dat$sc_class[rest_to_plot_dat$sc> quantile(rest_to_plot_dat$sc,2/3)]<-"High"
rest_to_plot_dat$sc_class <- factor(rest_to_plot_dat$sc_class , levels=c("Low","Medium","High"))

rest_to_plot_dat$npix_class<-"Low"
rest_to_plot_dat$npix<-rest_to_plot_dat$npix + rnorm(nrow(rest_to_plot_dat),sd=1e-6)
rest_to_plot_dat$npix_class[rest_to_plot_dat$npix>= quantile(rest_to_plot_dat$npix,1/3)]<-"Medium"
rest_to_plot_dat$npix_class[rest_to_plot_dat$npix>= quantile(rest_to_plot_dat$npix,2/3)]<-"High"
rest_to_plot_dat$npix_class <- factor(rest_to_plot_dat$npix_class , levels=c("Low","Medium","High"))



boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$signaldiff, 
        ylab="Classification Error" , xlab="Unit hydrograph variability",col="steelblue")
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/diff_vs_error.jpeg",width = 5, height = 5, units = 'in',res = 400)
boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$signaldiff, 
        ylab="Classification Error" , xlab="Unit hydrograph variability",col="steelblue")
dev.off()


boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$maxR2, 
        ylab="Classification Error" , xlab="Signal-to-noise ratio",col="steelblue")
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/stn_vs_error.jpeg",width = 5, height = 5, units = 'in',res = 400)
boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$maxR2, 
        ylab="Classification Error" , xlab="Signal-to-noise ratio",col="steelblue")
dev.off()


boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$nyears, 
        ylab="Classification Error" , xlab="Years of available data",col="steelblue")
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/years_vs_error.jpeg",width = 5, height = 5, units = 'in',res = 400)
boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$nyears, 
        ylab="Classification Error" , xlab="Length of available time series",col="steelblue")
dev.off()


boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$sc_class, 
        ylab="Classification Error" , xlab="Precipitation spatial autocorrelation",col="steelblue")
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/cor_vs_error.jpeg",width = 5, height = 5, units = 'in',res = 400)
boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$sc_class, 
        ylab="Classification Error" , xlab="Precipitation spatial autocorrelation",col="steelblue")
dev.off()


boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$npix_class, 
        ylab="Classification Error" , xlab="Number of pixels",col="steelblue")
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/size_vs_error.jpeg",width = 5, height = 5, units = 'in',res = 400)
boxplot(rest_to_plot_dat$missclass ~ rest_to_plot_dat$npix_class, 
        ylab="Classification Error" , xlab="Number of pixels",col="steelblue")
dev.off()


fit<-lm(missclass~ maxR2 + nyears + sc_class + npix_class,data = rest_to_plot_dat)
summary(fit)

plot_summs(fit, plot.distributions = TRUE, inner_ci_level = .99)


##################################################################################################################################
###################################################   get average results   #####################################################
##################################################################################################################################

what=readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/results_list.rds"))
ncat<-150

mean(what$missclass[1:ncat,,,,,3,drop = F])
sd(what$missclass[1:ncat,,,,,3,drop = F])


mean(what$missclass[1:ncat,,,,,2,drop = F])
sd(what$missclass[1:ncat,,,,,2,drop = F])


mean(what$missclass[1:ncat,,,,,1,drop = F])
sd(what$missclass[1:ncat,,,,,1,drop = F])



mean(what$fp[1:ncat,,,,,3,drop = F],na.rm = T)
sd(what$fp[1:ncat,,,,,3,drop = F],na.rm = T)


mean(what$fp[1:ncat,,,,,2,drop = F],na.rm = T)
sd(what$fp[1:ncat,,,,,2,drop = F],na.rm = T)


mean(what$fp[1:ncat,,,,,1,drop = F],na.rm = T)
sd(what$fp[1:ncat,,,,,1,drop = F],na.rm = T)



mean(what$fn[1:ncat,,,,,3,drop = F],na.rm = T)
sd(what$fn[1:ncat,,,,,3,drop = F],na.rm = T)


mean(what$fn[1:ncat,,,,,2,drop = F],na.rm = T)
sd(what$fn[1:ncat,,,,,2,drop = F],na.rm = T)


mean(what$fn[1:ncat,,,,,1,drop = F],na.rm = T)
sd(what$fn[1:ncat,,,,,1,drop = F],na.rm = T)







####################################################################################################################################################
#############################################   continuous changes in scenarios analysis   ##########################
##################################################################################################################################################
rm(list = ls())
library(caret)
library(ranger)
library(ncdf4)
library(ncdf4.helpers)
library(lubridate)
library(infotheo)
library(sf)
library(sp)
library(jtools)

get_cat_delin<-function(p,q,md){
  if(md=="mi"){
    rep(0,nrow(p))#return(get_mi(p=p,q=q))
  } else if(md=="te"){
    rep(0,nrow(p))#return(get_te(p=p,q=q))
  } else if(md=="cause"){
    return(get_cause(p=p,q=q))
  } else{
    return(-99)
  }
}

get_mi<-function(p,q){
  in_cat_pred<-rep(1,nrow(p))
  nlag<-3
  Y_t=q[(1+nlag):length(q)]
  fulldat<-as.data.frame(matrix(0,nrow = length(Y_t),ncol = nlag+1))
  for(gridsquare in 1:nrow(p)){
    for(lag in 0:nlag){
      fulldat[,lag+1]<- p[gridsquare,(nlag- lag+1):(length(q)-lag)]
    }
    mod<-ranger(x=fulldat,y=Y_t,num.trees = 50,write.forest = F,splitrule = "extratrees",replace = F,sample.fraction = 0.4)
    in_cat_pred[gridsquare]<-mod$r.squared
    #in_cat_pred[gridsquare]<-mutinformation(discretize(fulldat),Y_t,method = "mm")
  }
  
  bs_mi_res<-rep(0,20)
  for(b in 1:length(bs_mi_res)){
    for(lag in 0:nlag){
      fulldat[,lag+1]<-p[sample(1:nrow(p),1),sample(1:ncol(p),length(Y_t),replace = T)]
    }
    mod<-ranger(x=fulldat,y=Y_t,num.trees = 50,write.forest = F,splitrule = "extratrees",replace = F,sample.fraction = 0.4)
    bs_mi_res[b]<-mod$r.squared
    #bs_mi_res[b]<-mutinformation(discretize(fulldat),Y_t,method = "mm")
  }
  thresh<-3*sd(bs_mi_res)
  in_cat_pred<-as.numeric(in_cat_pred > thresh)
  return(in_cat_pred)
}

# curtime<-Sys.time()
# results<-get_cat_delin(p=useable_p,q=useable_streamflow,md=curmeth)
# Sys.time()-curtime


get_te<-function(p,q){
  in_cat_pred<-sample(c(0,1),nrow(p),replace = T)
  nlag<-3
  Y_t=q[(1+nlag):length(q)]
  #Y_t= discretize(Y_t)
  fulldat1<-data.frame(prev_q1= q[(nlag):(length(q)-1)],prev_q2= q[(nlag-1):(length(q)-2)],prev_q3= q[(nlag-2):(length(q)-3)])
  #partdat_mi<-mutinformation(discretize(fulldat1),Y_t,method = "mm")
  mod<-ranger(x=fulldat1,y=Y_t,num.trees = 50,write.forest = F,splitrule = "extratrees",replace = F,sample.fraction = 0.4)
  partdat_mi<-mod$r.squared
  for(gridsquare in 1:nrow(p)){
    fulldat<-fulldat1
    for(lag in 0:nlag){
      fulldat<- cbind(fulldat,p[gridsquare,(nlag- lag+1):(length(q)-lag)])
    }
    #fulldat<-discretize(fulldat)
    #fulldat_mi<-mutinformation(fulldat,Y_t,method = "mm")
    mod<-ranger(x=fulldat,y=Y_t,num.trees = 50,write.forest = F,splitrule = "extratrees",replace = F,sample.fraction = 0.4)
    fulldat_mi<-mod$r.squared
    in_cat_pred[gridsquare]<-fulldat_mi- partdat_mi
  }
  
  bs_mi_res<-rep(0,20)
  for(b in 1:length(bs_mi_res)){
    fulldat<-fulldat1
    for(lag in 0:nlag){
      fulldat<- cbind(fulldat,p[sample(1:nrow(p),1),sample(1:ncol(p),length(Y_t),replace = T)])
    }
    #fulldat<-discretize(fulldat)
    #fulldat_mi<-mutinformation(fulldat,Y_t,method = "mm")
    mod<-ranger(x=fulldat,y=Y_t,num.trees = 50,write.forest = F,splitrule = "extratrees",replace = F,sample.fraction = 0.4)
    fulldat_mi<-mod$r.squared
    bs_mi_res[b]<-fulldat_mi- partdat_mi
  }
  thresh<-3*sd(bs_mi_res)
  in_cat_pred<-as.numeric(in_cat_pred > thresh)
  return(in_cat_pred)
}




get_cause<-function(p,q){
  nlag<-10
  Y_t=q[(nlag):length(q)]
  fulldat<-as.data.frame(matrix(0,nrow = length(Y_t),ncol = nlag*nrow(p)+1))
  fulldat$V1<-Y_t
  counter<-2
  for(i in 1:nrow(p)){
    for(j in 1:nlag){
      fulldat[,counter]<- p[i,(nlag- j+1):(ncol(p)-j+1)]
      counter=counter+1
    }
  }
  fit<-lm(V1~.,data = fulldat)
  in_cat_pred<-coef(fit)[-1]
  thresh=as.numeric(3*summary(fit)$coefficients[-1, 2])
  covmat<-vcov(fit)[-1,-1]
  
  pixel_contributions<-rep(0,nrow(p))
  pixel_threshes<-rep(0,nrow(p))
  for(j in 1:nrow(p)){
    pixel_contributions[j]<-sum(in_cat_pred[(1+(j-1)*nlag):(j*nlag)])
    pixel_threshes[j]<-sqrt(sum(covmat[(1+(j-1)*nlag):(j*nlag),(1+(j-1)*nlag):(j*nlag)]))
  }
  
  #in_cat_pred<-as.numeric(in_cat_pred > thresh)[seq(from=2,by=nlag,length.out=nrow(p))]
  in_cat_pred<-as.numeric(pixel_contributions>(3*pixel_threshes))
  
  return(in_cat_pred)
}



shapename<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/shpfile_withp.rds"))

gids<-shapename$gauge_id
gids<-gids#[1:1555]

nreps<-5
meth<-c("mi","te","cause")

rest_to_plot_dat<-data.frame(catchment=rep("",length(gids)*nreps*length(meth)))
rest_to_plot_dat$classAcc<-0
rest_to_plot_dat$true_pos<-0
rest_to_plot_dat$true_neg<-0
rest_to_plot_dat$signaldiff<-0
rest_to_plot_dat$maxR2<-0
rest_to_plot_dat$nyears<-0
rest_to_plot_dat$spatial_cor<-0
rest_to_plot_dat$npix<-0
rest_to_plot_dat$autocor<-0
rest_to_plot_dat$pmean<-0
rest_to_plot_dat$method<-"mi"
curloc<-1


for(h in 1:length(gids)){
  curcat<-gids[h]
  p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
  posinfo<-p_data$posinfo
  precip_mat<-p_data$precip_mat
  tot_incat<-sum(posinfo$incat)
  
  for(i in 1:nreps){
    
  
    thetas<-runif(tot_incat,min = 1,max=5)
    alphas<-3/thetas
    #hist(rgamma(10000,shape = alphas[1],scale = thetas[1]),breaks=100,xlim = c(0,15))
    gt_streamflow<-rep(0, ncol(precip_mat))
    for(grid in 1:tot_incat){
      uh<- dgamma(0:49 +0.5,shape = alphas[grid],scale = thetas[grid])
      uh<-uh/sum(uh)
      what<-convolve(precip_mat[which(posinfo$incat==1)[grid],],uh,type="open",conj = F)
      gt_streamflow<-gt_streamflow + what[(length(uh)):length(what)]
    }
    
    curmaxR2<-runif(1,min=0.5,max=0.9)
    noiseVar<-(var(gt_streamflow)-curmaxR2*var(gt_streamflow))/curmaxR2
    noise<-as.numeric(arima.sim(list(order=c(1,0,0),ar=0.1),n=length(gt_streamflow),sd=1))
    noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
    streamflow=gt_streamflow+noise
    
    curyears<-sample(5:38,1)
    useable_days<- length(seq(from=as.Date("1983-01-01"),to=as.Date(paste0(1983+curyears-1,"-12-31")),by="day"))
    useable_p<-precip_mat[,1:useable_days]
    useable_streamflow<-streamflow[1:useable_days]
    
    
    for(m in 1:length(meth)){
      if(tot_incat>1) rest_to_plot_dat$signaldiff[curloc]<-mean(alphas*thetas^2)
      curmeth<-meth[m]
      results<-get_cat_delin(p=useable_p,q=useable_streamflow,md=curmeth)
      rest_to_plot_dat$catchment[curloc]<-curcat
      rest_to_plot_dat$classAcc[curloc]<-mean(results== posinfo$incat)
      rest_to_plot_dat$true_pos[curloc]<-mean(results[posinfo$incat==1]== posinfo$incat[posinfo$incat==1])
      rest_to_plot_dat$true_neg[curloc]<-mean(results[posinfo$incat==0]== posinfo$incat[posinfo$incat==0])
      rest_to_plot_dat$maxR2[curloc]<-curmaxR2
      rest_to_plot_dat$nyears[curloc]<-curyears
      rest_to_plot_dat$npix[curloc]<-nrow(posinfo)
      rest_to_plot_dat$spatial_cor[curloc]<-shapename$cor[h]
      rest_to_plot_dat$method[curloc]<-curmeth
      rest_to_plot_dat$autocor[curloc]<-mean(posinfo$auto_cor)
      rest_to_plot_dat$pmean[curloc]<-shapename$pmean[h]
      curloc<-curloc+1
    }
    
  }
  print(h)
}

saveRDS(rest_to_plot_dat, 
        file = paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/results_dat_contin.rds"))



rest_to_plot_dat<-readRDS("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/results_dat_contin.rds")




rest_to_plot_dat$nyears<-rest_to_plot_dat$nyears/100
rest_to_plot_dat$npix<-rest_to_plot_dat$npix/10
rest_to_plot_dat$spatial_cor<-rest_to_plot_dat$spatial_cor*10
rest_to_plot_dat$autocor<-rest_to_plot_dat$autocor*10



fit<-lm(classAcc~signaldiff + maxR2+ nyears + spatial_cor + npix,data = rest_to_plot_dat[seq(from=3,by=3,to=nrow(rest_to_plot_dat)),])
summary(fit)

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12))

plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/regressionCoefs.jpeg",width = 9, height = 5, units = 'in',res = 400)
plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

dev.off()


#try dry vs wet regions:


fit<-lm(classAcc~signaldiff + maxR2+ nyears + spatial_cor + npix,data = rest_to_plot_dat[rest_to_plot_dat$method=="cause" & rest_to_plot_dat$pmean>2.6,])
summary(fit)

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12))

plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/regressionCoefs_wet.jpeg",width = 9, height = 5, units = 'in',res = 400)
plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

dev.off()



fit<-lm(classAcc~signaldiff + maxR2+ nyears + spatial_cor + npix,data = rest_to_plot_dat[rest_to_plot_dat$method=="cause" & rest_to_plot_dat$pmean<2.6,])
summary(fit)

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Helvetica'),
        legend.title=element_blank(), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size = 12))

plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/regressionCoefs_dry.jpeg",width = 9, height = 5, units = 'in',res = 400)
plot_summs(fit, inner_ci_level = .99,legend.title = "",
           coefs=c("Spread of unit hydrograph" = "signaldiff", "Signal-to-noise (NSE)" = "maxR2","Length of time series"="nyears", "Spatial autocorrelation"="spatial_cor", "Temporal autocorrelation"="autocor", "# of rainfall pixels"="npix")) + apatheme + labs(x="Impact on classification accuracy (Regression coefficents)",y="") + xlim(-0.025, 0.05)

dev.off()








rest_to_plot_dat$classAcc_perf<-0
rest_to_plot_dat$classAcc_perf[rest_to_plot_dat$classAcc==1]<-1
mylogit <- glm(classAcc_perf ~ signaldiff + maxR2+ nyears + spatial_cor + npix,data = rest_to_plot_dat[seq(from=3,by=3,to=nrow(rest_to_plot_dat)),], family = "binomial")
summary(mylogit)

confusionMatrix(data = as.factor(as.numeric(mylogit$fitted.values>0.5)), reference = as.factor(mylogit$data$classAcc_perf))

plot_summs(mylogit, plot.distributions = TRUE, inner_ci_level = .99)



mean(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="mi"])
mean(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="mi"])
mean(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="mi"])


mean(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="te"])
mean(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="te"])
mean(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="te"])


mean(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="cause"])
mean(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="cause"])
mean(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="cause"])


sd(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="mi"])
sd(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="mi"])
sd(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="mi"])


sd(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="te"])
sd(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="te"])
sd(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="te"])


sd(rest_to_plot_dat$classAcc[rest_to_plot_dat$method=="cause"])
sd(rest_to_plot_dat$true_pos[rest_to_plot_dat$method=="cause"])
sd(rest_to_plot_dat$true_neg[rest_to_plot_dat$method=="cause"])



########################################## plot results for biggest catchment  #################################################

curcat<-"DEF13090"  #"DE213740"
p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
posinfo<-p_data$posinfo
precip_mat<-p_data$precip_mat
tot_incat<-sum(posinfo$incat)
x_repr = shapename$geometry[which(shapename$gauge_id==curcat)]
x_repr= x_repr[[1]][[1]]


plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual/Topographic Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/true_bigCat.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual/Topographic Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


alphas<-runif(tot_incat,min = 1,max=2)
thetas<-runif(tot_incat,min = 1,max=2)
#hist(rgamma(10000,shape = alphas[1],scale = thetas[1]),breaks=100,xlim = c(0,15))
gt_streamflow<-rep(0, ncol(precip_mat))
for(grid in 1:tot_incat){
  uh<- dgamma(0:49 +0.5,shape = alphas[grid],scale = thetas[grid])
  uh<-uh/sum(uh)
  what<-convolve(precip_mat[which(posinfo$incat==1)[grid],],uh,type="open",conj = F)
  gt_streamflow<-gt_streamflow + what[(length(uh)):length(what)]
}

curmaxR2<-0.8
noiseVar<-(var(gt_streamflow)-curmaxR2*var(gt_streamflow))/curmaxR2
noise<-as.numeric(arima.sim(list(order=c(1,0,0),ar=0.1),n=length(gt_streamflow),sd=1))
noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
streamflow=gt_streamflow+noise

curyears<-38
useable_days<- length(seq(from=as.Date("1983-01-01"),to=as.Date(paste0(1983+curyears-1,"-12-31")),by="day"))
useable_p<-precip_mat[,1:useable_days]
useable_streamflow<-streamflow[1:useable_days]

results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="mi")

plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Mutual Information Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/mi_bigCat.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Mutual Information Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="te")
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Transfer Entropy Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/te_bigCat.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Transfer Entropy Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="cause")
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/causal_bigCat.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()









curcat<-"DEA11040"  #"DE213740"
p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
posinfo<-p_data$posinfo
precip_mat<-p_data$precip_mat
tot_incat<-sum(posinfo$incat)
x_repr = shapename$geometry[which(shapename$gauge_id==curcat)]
x_repr= x_repr[[1]][[1]]


plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual/Topographic Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/true_bigCat2.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual/Topographic Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


alphas<-runif(tot_incat,min = 1,max=2)
thetas<-runif(tot_incat,min = 1,max=2)
#hist(rgamma(10000,shape = alphas[1],scale = thetas[1]),breaks=100,xlim = c(0,15))
gt_streamflow<-rep(0, ncol(precip_mat))
for(grid in 1:tot_incat){
  uh<- dgamma(0:49 +0.5,shape = alphas[grid],scale = thetas[grid])
  uh<-uh/sum(uh)
  what<-convolve(precip_mat[which(posinfo$incat==1)[grid],],uh,type="open",conj = F)
  gt_streamflow<-gt_streamflow + what[(length(uh)):length(what)]
}

curmaxR2<-0.85
noiseVar<-(var(gt_streamflow)-curmaxR2*var(gt_streamflow))/curmaxR2
noise<-as.numeric(arima.sim(list(order=c(1,0,0),ar=0.1),n=length(gt_streamflow),sd=1))
noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
streamflow=gt_streamflow+noise

curyears<-20
useable_days<- length(seq(from=as.Date("1983-01-01"),to=as.Date(paste0(1983+curyears-1,"-12-31")),by="day"))
useable_p<-precip_mat[,1:useable_days]
useable_streamflow<-streamflow[1:useable_days]

results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="mi")

plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Mutual Information Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/mi_bigCat2.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Mutual Information Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="te")
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Transfer Entropy Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/te_bigCat2.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Transfer Entropy Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()


results<-get_cat_delin(p=useable_p,q=useable_streamflow,md="cause")
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/causal_bigCat2.jpeg",width = 4, height = 4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()




########################################## compare with topography #################################################

get_cause<-function(p,q){
  nlag<-12
  Y_t=q[(nlag):length(q)]
  fulldat<-as.data.frame(matrix(0,nrow = length(Y_t),ncol = nlag*nrow(p)+1))
  fulldat$V1<-Y_t
  counter<-2
  for(i in 1:nrow(p)){
    for(j in 1:nlag){
      fulldat[,counter]<- p[i,(nlag- j+1):(ncol(p)-j+1)]
      counter=counter+1
    }
  }
  fit<-lm(V1~.,data = fulldat)
  in_cat_pred<-coef(fit)[-1]
  thresh=as.numeric(3*summary(fit)$coefficients[-1, 2])
  covmat<-vcov(fit)[-1,-1]
  
  pixel_contributions<-rep(0,nrow(p))
  pixel_threshes<-rep(0,nrow(p))
  for(j in 1:nrow(p)){
    pixel_contributions[j]<-sum(in_cat_pred[(1+(j-1)*nlag):(j*nlag)])
    pixel_threshes[j]<-sqrt(sum(covmat[(1+(j-1)*nlag):(j*nlag),(1+(j-1)*nlag):(j*nlag)]))
  }
  
  #in_cat_pred<-as.numeric(in_cat_pred > thresh)[seq(from=2,by=nlag,length.out=nrow(p))]
  in_cat_pred<-as.numeric(pixel_contributions>(3*pixel_threshes))
  
  return(in_cat_pred)
}

curcat<-"DE811440"  #"DE213740"
p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
posinfo<-p_data$posinfo
precip_mat<-p_data$precip_mat
tot_incat<-sum(posinfo$incat)
x_repr = shapename$geometry[which(shapename$gauge_id==curcat)]
x_repr= x_repr[[1]][[1]]
posinfo$incat[posinfo$lat==53.75 & posinfo$lon==13]<-1
posinfo$incat[posinfo$lon==13.75 & posinfo$lat==53.5]<-1
posinfo$incat[posinfo$lon==13.75 & posinfo$lat==53.25]<-1
posinfo$incat[posinfo$lon==13.5 & posinfo$lat==53.25]<-1

posinfo_true<-posinfo


plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue","grey")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/true_Cat_topCompare.jpeg",width = 4.35, height = 4.4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue","grey")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="True Virtual Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()



p_data<-readRDS(paste0("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/data/precipStuff_",curcat,".rds"))
posinfo<-p_data$posinfo
precip_mat<-p_data$precip_mat
tot_incat<-sum(posinfo$incat)
x_repr = shapename$geometry[which(shapename$gauge_id==curcat)]
x_repr= x_repr[[1]][[1]]
posinfo$incat[posinfo$lat==53.75 & posinfo$lon==13]<-1
#posinfo$incat[posinfo$lon==13.75 & posinfo$lat==53.5]<-1
#posinfo$incat[posinfo$lon==13.75 & posinfo$lat==53.25]<-1
#posinfo$incat[posinfo$lon==13.5 & posinfo$lat==53.25]<-1

posinfo_topo<-posinfo

plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue","grey")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="Topographic Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/topo_Cat_topCompare.jpeg",width = 4.35, height = 4.4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue","grey")[posinfo$incat+1],cex=6,xlab="longitude",ylab="latitude",main="Topographic Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()



alpha=3
theta=1
fastker<-hist(rgamma(900000,shape = alpha,scale = theta),breaks=100,xlim = c(0,50))
alpha=2
theta=8
slowker<-hist(rgamma(900000,shape = alpha,scale = theta),breaks=100,xlim = c(0,50))

plot(fastker$mids,fastker$density,type="l",xlim=c(0,49),lwd=3,ylab="Response weight",xlab="lag (days)")
lines(slowker$mids,slowker$density,col="red",lwd=3)
legend("topright",legend = c("Within-boundary unit hydrograph","Cross-boundary unit hydrograph"),lwd=5,col = c("black","red"))

jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/unithydrographs_lastexp.jpeg",width = 4.35, height = 4.4, units = 'in',res = 400)
par(mar = c(4, 4, 1, 2))
plot(fastker$mids,fastker$density,type="l",xlim=c(0,49),lwd=3,ylab="Response magnitude",xlab="lag (days)")
lines(slowker$mids,slowker$density,col="red",lwd=3)
legend("topright",legend = c("Within-boundary unit hydrograph","Cross-boundary unit hydrograph"),cex=0.85,lwd=5,col = c("black","red"))
dev.off()

tot_incat<-sum(posinfo_true$incat)
gt_streamflow<-rep(0, ncol(precip_mat))
for(grid in 1:tot_incat){
  curgrid<-which(posinfo_true$incat==1)[grid]
  if(posinfo_topo$incat[curgrid]==1 & posinfo_true$incat[curgrid]==1){
    alpha=3
    theta=1
    #hist(rgamma(10000,shape = alpha,scale = theta),breaks=100,xlim = c(0,50))
  } else{
    alpha=2
    theta=8
    #hist(rgamma(10000,shape = alpha,scale = theta),breaks=100,xlim = c(0,50))
  }
  uh<- dgamma(0:49 +0.5,shape = alpha,scale = theta)
  uh<-uh/sum(uh)
  what<-convolve(precip_mat[which(posinfo_true$incat==1)[grid],],uh,type="open",conj = F)
  gt_streamflow<-gt_streamflow + what[(length(uh)):length(what)]
}

curmaxR2<-0.9
noiseVar<-(var(gt_streamflow)-curmaxR2*var(gt_streamflow))/curmaxR2
noise<-as.numeric(arima.sim(list(order=c(1,0,0),ar=0.01),n=length(gt_streamflow),sd=1))
noise<- (noise - mean(noise))* sqrt(noiseVar)/sd(noise)
streamflow=gt_streamflow+noise

curyears<-38
useable_days<- length(seq(from=as.Date("1983-01-01"),to=as.Date(paste0(1983+curyears-1,"-12-31")),by="day"))
useable_p<-precip_mat[,1:useable_days]
useable_streamflow<-streamflow[1:useable_days]


results<-get_cause(p=useable_p,q=useable_streamflow)
plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
jpeg("C:/Users/joeja/Desktop/researchMasters/CatchmentDelin/Paper/figures/causal_Cat_topCompare.jpeg",width = 4.35, height = 4.4, units = 'in',res = 400)
plot(posinfo$lon,posinfo$lat-0.025,pch=15,col=c("lightsalmon","steelblue")[results+1],cex=6,xlab="longitude",ylab="latitude",main="Causal Inference Catchment")
lines(x_repr,col="black",lwd=3)
dev.off()














x=rnorm(10000)
ytm1<-x+rnorm(10000)
yt<-x+ytm1+rnorm(10000)

summary(lm(yt~x+ytm1))



