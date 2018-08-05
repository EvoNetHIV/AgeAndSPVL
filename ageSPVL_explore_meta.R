
nruns <- 10
popatts <- popsumm <- list()

agecoef <- iSPVL <- ageinf <- prev <- numinc <- agematch <- rep(NA, nruns)

for (i in 1:nruns) {
  if(i<10) filler <- "0" else filler <- ""
  obj <- get(paste("ageSPVL_m",filler,i,sep=""))
  popatts[[i]] <- obj$pop[[1]]
  popsumm[[i]] <- get(paste("ageSPVL_m",filler,i,sep=""))$popsumm[[1]]
  agecoef[i] <- lm(log(popatts[[i]]$SetPoint[popatts[[i]]$Time_Inf>0],10)~
               popatts[[i]]$age_infection[popatts[[i]]$Time_Inf>0])$coef[[2]]
  iSPVL[i] <- mean(log10(popatts[[i]]$SetPoint[popatts[[i]]$Time_Inf>0]),na.rm=TRUE)
  ageinf[i] <- mean(popatts[[i]]$age_infection,na.rm=TRUE)
  prev[i] <- tail(popsumm[[i]]$prevalence,1)
  numinc[i] <- sum(popatts[[i]]$Time_Inf>0, na.rm=TRUE)
  agematch[i] <- obj$nwparam[[1]]$coef.form['absdiff.sqrt_age']
}

plot(numinc, prev)
plot(prev, iSPVL)
plot(prev, iSPVL, type='b')
iSPVL
prev
plot(agecoef); abline(h=0)
plot(agematch)


