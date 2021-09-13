#' generate a simulation of PTR-TOF-MS raw data from exhaled breath around a nominal mass
#' from a list of temporal profile
#' @param mzProfileExp
#' @param mzAmbiantProfile
#' @param typeProfile could be expiration, ambient or constant
#' @param mz mass nominal of the simulation
#' @param mzLarge width of the simulation
#' @param npeak number of peaks
#' @param ppmSep separations between the peaks
#' @param heigth intensity of the peak
#' @param ratio ratio of the neighbor peaks
#' @param backnoiseRate lambda of the poison noise
#' @param pulseMean mean of intensity of the gaussian single ions pulse mean
#' @param pulseSd sd of the gaussian single ions pulse
#' @param dbaseLine
#' @param coefbaseline
#' @param plot should the simulation be plotted
#' @import rhdf5, signal
#' readtsvxpirationProfile)
#' data(AmbiantProfile)
#' simulate(expirationProfile,AmbiantProfile)
simulate<-function(mzProfileExp,mzAmbiantProfile,typeProfile=c("expiration","ambiant","constant")[1:2],
                   mz=80.05,mzLarge=0.5,
                   ppmSep=250,backnoiseRate=0.1,pulseMean=11,pulseSd=3,
                   dbaseLine=0,coefbaseline=0,
                   heigth=10^3,
                   ratio=c(1.5),
                   plot=TRUE,npeak=length(typeProfile)){

  i <- sample(seq_along(mzProfile),size = 1)
  if(length(typeProfile) != npeak) stop("npeak != length(typeProfile)")

  file <- names(mZexpiration)[i]
  filefullNames <- fulleNameslistPeak[ file == basename(fulleNameslistPeak) ]

  profile<-list()
  if(typeProfile=="ambiant") listProfile<-mzAmbiantProfile else listProfile<-mzProfileExp
  for(j in 1:npeak){
    profile[[j]] <- eval(parse(text =typeProfile[j]))(i,listProfile)
  }

  time<-as.numeric(Reduce(c,lapply(strsplit(colnames(mzProfile[[i]]),"-"),function(x) utils::tail(x,1))))
  temporalEvolution<- do.call(rbind,lapply(profile,function(x) x$profil))


  mzLarge<-mzLarge
  mzNom<-round(mz)


  if(npeak>=2){
    for(j in 2:npeak){
      mz[j]<- mz[j-1]+ppmSep*mzNom/10^6}
  }

  #asymetrie<- runif(n = 1,min = 0.35,max = 0.65)
  #res<-rnorm(npeak,mean = 5000,sd = 100)
  parameter.1<- mzNom/(rnorm(npeak,mean = 5000,sd = 500)*2)
  #parameter.2 <- parameter.1*asymetrie
  parameter.2<- mzNom/(rnorm(npeak,mean = 5000,sd = 500)*2)

  #parameter.1<-(mzNom/res)*asymetrie
  #parameter.2<-(mzNom/res)*(1-asymetrie)
  parameterPeak<-data.frame(Mz=mz,
                            parameter.1=parameter.1,
                            parameter.2=parameter.2)

  mzVn<-rhdf5::h5read(filefullNames, "/FullSpectra/MassAxis", bit64conversion='bit64')

  mz<- mzVn[mzVn > mzNom-mzLarge & mzVn < mzNom+mzLarge]


  #simulation
  ratiox<-1

  if(npeak>=2) ratiox<-c(ratiox,rep(ratio,npeak-1))


  backnoiseMat<-rep(0,length(time)*length(mz))

  for(c in seq_along(backnoiseMat)){
    k<-rpois(1,lambda = backnoiseRate )
    phd <- rnorm(k,mean = pulseMean , sd = pulseSd )
    phd[phd<0]<-0
    backnoiseMat[c]<-sum(phd)
  }

  #backnoiseMat<-round(rexp(length(time)*length(mz),rate =backnoise),4)
  #backnoiseMat<-rpois(length(time)*length(mz),lambda = backnoise)
  #backnoiseMat<-rpois(length(time)*length(mz),lambda = 0.1)*14.5
  #law<-hist(spnoise,breaks = 500,freq = FALSE)
  #backnoiseMat<- sample(x = law$breaks[-length(law$breaks)], size = length(time),
  #                      prob = law$density,replace = TRUE)

  baseline <- coefbaseline * mz^dbaseLine

  simulatedPeak<-lapply(seq(1,npeak),function(p){
    rawSimulate<-matrix(0,ncol=length(time),nrow=length(mz))
    peakSimulate<-rawSimulate
    for(j in seq_along(time)){

      height<- trunc(temporalEvolution[p,j] * heigth * ratiox[p])
      if(height<0) height<-0

      law <- sech2(p = parameterPeak$Mz[p],
                             lf =parameterPeak$parameter.1[p],
                             lr=parameterPeak$parameter.2[p],
                             h=1 ,
                             x=c(mz[1],mz+diff(mz)[1]))

      tirage <- hist(sample(x = c(mz[1],mz+diff(mz)[1]),
                            size = as.integer(height),
                            prob = law/sum(law),replace = TRUE),
                     breaks =  c(mz[1],mz+diff(mz)[1]),
                     plot = FALSE)

      peaks<-tirage$counts
      peakSimulate[,j]<-peaks

      #peakSimulate[,j]<-law[-length(law)] * height
      # peaks<- peakSimulate[,j] + rnorm(mz,sd=50)* ( peakSimulate[,j] >0.05*height)
      # peaks[peaks<0]<-0


      rawSimulate[,j] <- peaks + baseline


    }

    rawSimulate<-round(rawSimulate,1)
    return(list(rawSimulate,peakSimulate))
  })


  rawSimulate<-simulatedPeak[[1]][[1]] + backnoiseMat
  if(npeak>=2){
    for(k in 2:npeak){
      rawSimulate<-rawSimulate + simulatedPeak[[k]][[1]]
    }
  }

  peakSimulate<-Reduce(cbind,lapply(simulatedPeak,function(x) colSums(x[[1]])))
  peakSimulate<-matrix(peakSimulate,ncol = npeak)
  rownames(rawSimulate)<-mz
  colnames(rawSimulate)<-time

  param<-list(file=file,
              profilmz=profile,
              time=length(time),
              ppmSep=ppmSep,
              backnoise=backnoiseRate,
              timeProcess=NA,
              ratio=ratio,
              heigth=heigth,
              Mz= parameterPeak$Mz)

  mMean<-mean(c(parameterPeak$Mz))

  if(plot){
    image(z=t(rawSimulate[mz > mMean-mzNom*500/10^6 & mz < mMean+mzNom*500/10^6 ,]),x=
            time,
          y=mz[mz > mMean-mzNom*500/10^6 & mz < mMean+mzNom*500/10^6 ],
          main=paste("simulated data, ppmsep=",ppmSep,
                     ", count - ratio =",heigth,paste(ratio,collapse = " ")),xlab="time",ylab="mz")
  }

  return(list(rawSimulate=rawSimulate,
              peakSimulate=peakSimulate,
              param=param))
}


createEntireH5file<-function(simuList,file.origine,name,dirSimulate,fulleNameslistPeak){

  filefullNames <- fulleNameslistPeak[ file.origine == basename(fulleNameslistPeak)]
  file.new <- file.path(dirSimulate,paste0(name,".h5"))

  mzVn<-rhdf5::h5read(filefullNames, "/FullSpectra/MassAxis", bit64conversion='bit64')

  mzVnNew<-mzVn
  time<-as.numeric(Reduce(c,lapply(strsplit(colnames(mzProfile[[file.origine]]),"-"),function(x) utils::tail(x,1))))

  rawM<-matrix(0,nrow=length(mzVnNew),ncol=length(time))
  for( m in seq_along(simuList)){
    rawM[mzVnNew > m - 0.5 & mzVnNew < m+0.5,] <- simuList[[m]]$rawSimulate
  }

  rawM[rawM<0]<-0
  #rawM<-round(rawM)
  NbrWrite <- ceiling(length(time)/10)


  rawM_New<-matrix(0,nrow= nrow(rawM)*NbrWrite*10)
  dim(rawM_New)<- c(nrow(rawM),1,10,NbrWrite)

  for(j in 1:(NbrWrite-1)){
    rawM_New[,1,,j] <- rawM[,(10*(j-1)+1):min((j*10),length(time))]
  }
  j<-NbrWrite


  if(NbrWrite*10 == length(time)){
    rawM_New[,1,,j] <- rawM[,(10*(j-1)+1):min((j*10),length(time))]
  }else {
    lastWrite<-c(rawM[,(10*(j-1)+1):min((j*10),length(time))])
    lastWrite[(length(lastWrite)+1):(length(lastWrite)+length((length(time)+1):(NbrWrite*10))*nrow(rawM))]<-0
    rawM_New[,1,,NbrWrite]<- matrix(lastWrite,ncol=10)
  }

  rawM_New<-array(rawM_New, dim= c(length(mzVnNew),1,10,NbrWrite))

  file.origine<-filefullNames
  file.copy(filefullNames, file.new,overwrite = TRUE)
  #rawM_New[rawM_New<0.1]<-0
  #createH5file
  #
  # rhdf5::h5createFile(file.new)
  # read<-rhdf5::h5read(file.origine, "/", bit64conversion='bit64')
  #
  # ## craete group and subgroup
  # listGroup <- rhdf5::h5ls(file.origine)
  # group<-listGroup[listGroup$otype ==  "H5I_GROUP",]
  # for( i in 1:dim(group)[1])
  #   rhdf5::h5createGroup(file.new,paste(group$group[i],group$name[i],sep="/"))
  #
  # ## write data
  # data<-listGroup[listGroup$otype ==  "H5I_DATASET",]
  #
  # ## No change
  # DATASET<-paste(data$group,data$name,sep="/")
  # DATASETNochange<-  c( "/AcquisitionLog/Log",
  #                       "/AddTraces/PTR-Instrument/TwInfo" ,
  #                       "/AddTraces/PTR-Misc/TwInfo",
  #                       "/AddTraces/PTR-Reaction/TwInfo",
  #                       "/PTR-Concentration/TwInfo",
  #                       "/PTR-Peaktable/Info",
  #                       "/PTR-Peaktable/Data",
  #                       "/PTR-PrimaryIonSettings/Data",
  #                       "/PTR-PrimaryIonSettings/Info" ,
  #                       "/PTR-Transmission/Data",
  #                       "/PTR-Transmission/Info",
  #                       "/PeakData/PeakTable",
  #                       "/RawData/Ndigo5G/Temperature/TwInfo",
  #                       "/TPS2/TwInfo",
  #                       "/FullSpectra/MassCalibration",
  #                       "/PTR-Concentration/TwData",
  #                       "/PeakData/PeakData",
  #                       "/RawData/Ndigo5G/Temperature/TwData",
  #                       "/TPS2/TwData",
  #                       "/TimingData/BufTimes")
  #
  # for(i in 1:dim(data)[1]){
  #   harchi <- paste(data$group[i],data$name[i],sep="/")
  #   if(harchi %in% DATASETNochange){
  #     harchi_sep <- strsplit(harchi,"/")[[1]][-1]
  #     if(length(harchi_sep)==3)
  #       rhdf5::h5write(read[[harchi_sep[1]]][[harchi_sep[2]]][[harchi_sep[3]]],file.new,harchi)
  #     if(length(harchi_sep)==2)
  #       rhdf5::h5write(read[[harchi_sep[1]]][[harchi_sep[2]]],file.new,harchi)
  #   }
  # }
  #

  #rhdf5::h5write(mzVnNew,file.new, "/FullSpectra/MassAxis")
  #rhdf5::h5createDataset(file.new,"/FullSpectra/TofData", dim(rawM_New), chunk = c(dim(rawM_New)[1],1,10,1),level = 9)
  #raw_M_origine<- rhdf5::h5read(file.origine,"/FullSpectra/TofData")

  rhdf5::h5write(rawM_New,file.new, "/FullSpectra/TofData")
  rhdf5::h5write(rowSums(rawM_New),file.new, "/FullSpectra/SumSpectrum")



}


#normalize a expiration profile
expiration<-function(i,mzProfile){
  mzExp<-sample(seq_len(nrow(mzProfile[[i]])),size = 1)
  profil<-as.numeric(mzProfile[[i]][mzExp,])
  profil<-signal::sgolayfilt(profil,n = OptimalWindowsSG(profil,noiseacf = 0.1))
  profil<-profil/sum(profil)*length(profil)
  return(list(profil=profil,mz=mzExp))
}

#normalize a ambiant air VOC temporal profile
ambiant <-function(i,mzAmbiantProfile){
  mzAmbiant<-sample(seq_len(nrow(mzAmbiantProfile[[i]])),size = 1)
  profil<-as.numeric(mzAmbiantProfile[[i]][mzAmbiant,])
  profil<-signal::sgolayfilt(profil,n = OptimalWindowsSG(profil,noiseacf = 0.1))
  profil<-profil/sum(profil)*length(profil)
  return(list(profil=profil,mz=mzAmbiant))
}

#generate a contant vector
constant<- function(i,mzProfile){
  n<-length(mzProfile[[i]][1,])
  profil<-rep(1,n)
  return(list(profil=profil))
}


sech2 <- function(p, lf, lr, h, x, l.shape = NULL) {
  h/(cosh((log(sqrt(2) + 1)/lf) * (x - p))^2 * (x <= p) + cosh((log(sqrt(2) + 1)/lr) *
                                                                 (x - p))^2 * (x > p))
}

OptimalWindowsSG <- function(sp, noiseacf, d = 3) {
  n = 5
  spf <- signal::sgolayfilt(sp, p = d, n = n)
  res <- sp - spf
  acf_res0 <- stats::acf(res, plot = FALSE)[1]$acf[1]
  n = n + 2
  repeat {
    spf <- signal::sgolayfilt(sp, p = d, n = n)
    res <- sp - spf
    acf_res1 <- stats::acf(res, plot = FALSE)[1]$acf[1]
    if (abs(acf_res0 - noiseacf) < abs(acf_res1 - noiseacf))
      break
    # si res 0 est plus proche que res1
    n = n + 2
    if (n >= length(sp))
      break
    acf_res0 <- acf_res1
  }
  n
}
