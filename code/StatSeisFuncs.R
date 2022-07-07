####### A collection of code for  statistical seismology
######  Mostly exploratory data analysis with some simulation and data handling for specific catalogs


##### CORSSA FMD/Mc stuff
#### From http://www.corssa.org/en/articles/theme_4/ Completeness magnitude in earthquake catalogs
#### mbs method is altered slightly to fix error in bave (old code doesn't actually average) and to increase the number of steps considered

### Calculate frequency-magnitude distribution, called by other Mc functions
## Inputs: (vector of) magnitudes, binsize
## Outputs: dataframe with magnitude bins (mi), cumulative (cum) and non-cumulative (noncum) frequency-magnitude
fmd <- function(mag,mbin){
  mi <- seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm <- length(mi)
  cumnbmag <- numeric(nbm)
  nbmag <- numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] <- length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp <- c(cumnbmag,0)
  nbmag <- abs(diff(cumnbmagtmp))
  res <- list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}

### Calculate Mc by method of maximum curvature
## Inputs: (vector of) magnitudes, binsize
## Outputs: Mc value
maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}


### Calculate Mc by method of Goodness-of-fit test (GFT) [Wiemer & Wyss, 2000]
## Inputs: (vector of) magnitudes, binsize
## Outputs: Mc value, whether this is at 90 or 95% confidence, sequence of trialled Mc values (Mco) and associated R
gft <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  McBound <- maxc(mag,mbin)$Mc
  Mco <- McBound-0.4+(seq(15)-1)/10
  R <- numeric(15)
  for(i in 1:15){
    indmag <- which(mag > Mco[i]-mbin/2)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    FMDcum_model <- 10^(a-b*FMD$m)
    indmi <- which(FMD$m >= Mco[i])
    R[i] <- sum(abs(FMD$cum[indmi]-FMDcum_model[indmi]))/sum(FMD$cum[indmi])*100
    #in Wiemer&Wyss [2000]: 100-R
  }
  indGFT <- which(R <= 5)         #95% confidence
  if(length(indGFT) != 0){
    Mc <- Mco[indGFT[1]]
    best <- "95%"
  } else{
    indGFT <- which(R <= 10)    #90% confidence
    if(length(indGFT) != 0){
      Mc <- Mco[indGFT[1]]
      best <- "90%"
    } else{
      Mc <- McBound
      best <- "MAXC"
    }
  }
  return(list(Mc=Mc, best=best, Mco=Mco, R=R))
}

### Mc estimation by b-val Stability (MBS) [Cao & Gao, 2002], modification with Shi & Bolt [1982] uncertainty [Woesner & Wiemer, 2005]
## Inputs: (vector of) magnitudes, binsize
## Outputs: Mc estimate, list of trial Mc values and associated b value at each Mc, uncertainty in b from Shi and Bolt (1982) and 5-bin average b-value (bave)
mbs <- function(mag,mbin){
  McBound <- maxc(mag,mbin)$Mc
  ## Remember that mbin will affect how far your fixed bins stretch
  #Mco <- McBound-0.7+(seq(30)-1)/10
  Mco <- McBound+(seq(40)-1)/10
  #Mco <- McBound+1 +(seq(30)-1)/10
  #print(Mco)
  bi <- numeric(40); unc <- numeric(40)
  for(i in 1:40){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    bi[i] <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    unc[i] <- 2.3*bi[i]^2*sqrt(sum((mag[indmag]-mean(mag[indmag]))^2)/(nbev*(nbev-1)))
  }
  bave <- numeric(35)
  for(i in 1:35) bave[i] <- mean(bi[i:(i+5)])
  
  dbi_old <- abs(diff(bi))
  indMBS_old <- which(dbi_old <= 0.03)
  
  dbi <- abs(bave[1:35]-bi[1:35])
  indMBS <- which(dbi <= unc[1:35])
  Mc <- Mco[indMBS[1]]
  return(list(Mc=Mc, Mco=Mco, bi=bi, unc=unc, bave=bave))
}

###Sometimes you need more bins to see where (if) b stabilises, especially if bin size is small and the catalog spans a really wide range of magnitudes. 
####### Same as above, but also supply number of bins to consider as n
mbs_mod <- function(mag,mbin, n, calc_inc_size = 0.1){
  ### If n = number of bins to consider, n-5 is number of averages
  nl <- n-5
  McBound <- maxc(mag,mbin)$Mc
  
  Mco <- McBound+(seq(n)-1)/(1/calc_inc_size)
  
  bi <- numeric(n); unc <- numeric(n)
  for(i in 1:n){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    bi[i] <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    unc[i] <- 2.3*bi[i]^2*sqrt(sum((mag[indmag]-mean(mag[indmag]))^2)/(nbev*(nbev-1)))
  }
  bave <- numeric(nl)
  for(i in 1:nl) bave[i] <- mean(bi[i:(i+5)])
  
  dbi_old <- abs(diff(bi))
  indMBS_old <- which(dbi_old <= 0.03)
  
  dbi <- abs(bave[1:nl]-bi[1:nl])
  indMBS <- which(dbi <= unc[1:nl])
  Mc <- Mco[indMBS[1]]
  return(list(Mc=Mc, Mco=Mco, bi=bi, unc=unc, bave=bave))
}

#### Again change bounds on this - in this case start at M0 and extent to Max mag to avoid issues like having a kink in the GR!
#### Set increment inc over which to test Mc values, defaults to 0.1
## Inputs same as gft() but with added increment. 
gft_mod <- function(mag,mbin, inc=0.1){
  FMD <- fmd(mag,mbin)
  #McBound <- maxc(mag,mbin)$Mc
  #Mco <- McBound-0.4+(seq(n)-1)/10
  Mco <- seq(min(mag), max(mag), by=inc)
  n <- length(Mco)
  R <- numeric(n)
  for(i in 1:n){
    indmag <- which(mag > Mco[i]-mbin/2)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    FMDcum_model <- 10^(a-b*FMD$m)
    indmi <- which(FMD$m >= Mco[i])
    R[i] <- sum(abs(FMD$cum[indmi]-FMDcum_model[indmi]))/sum(FMD$cum[indmi])*100
    #in Wiemer&Wyss [2000]: 100-R
  }
  indGFT <- which(R <= 5)         #95% confidence
  if(length(indGFT) != 0){
    Mc <- Mco[indGFT[1]]
    best <- "95%"
  } else{
    indGFT <- which(R <= 10)    #90% confidence
    if(length(indGFT) != 0){
      Mc <- Mco[indGFT[1]]
      best <- "90%"
    } else{
      Mc <- McBound
      best <- "MAXC"
    }
  }
  return(list(Mc=Mc, best=best, Mco=Mco, R=R))
}


######################### Added by Kirsty (so more likely to be buggy :P)

##### Function to plot b-value stability graph
### Uses modified mbs so can provide number of bins (defaults to 40)
### Inputs: Only magnitudes and bin size required, but optionally the increment size, a different title for the plot, specified legend location, limits and legend flag
### Outputs: A plot
plot_mbs <- function(mag, mbin, n=40, calc_inc_size=0.1, main="", leg_loc="bottomright", ylim = c(0.5, 2), legend=TRUE){
  nl <- n-5
  Mbs_out <- mbs_mod(mag, mbin, n)
  plot(Mbs_out$Mco, Mbs_out$bi, pch=16, ylim=ylim, xlab="Mc", ylab="b-value", main=main)
  lines(Mbs_out$Mco, Mbs_out$bi, col="black", lty=1)
  abline(v=Mbs_out$Mc, lty=2)
  points(Mbs_out$Mco, Mbs_out$bi + Mbs_out$unc, col="red")
  lines(Mbs_out$Mco, Mbs_out$bi + Mbs_out$unc, col="red", lty=2)
  points(Mbs_out$Mco, Mbs_out$bi - Mbs_out$unc, col="red")
  lines(Mbs_out$Mco, Mbs_out$bi - Mbs_out$unc, col="red", lty=2)
  lines(Mbs_out$Mco[1:nl], Mbs_out$bave, col="gray")
  if(legend == TRUE){
    legend(leg_loc, lty=c(1,2,1, 2), col=c("black", "red", "gray", "black"), legend=c("b", "b uncertainty", "5-bin average b", "Mc"))
  }
}


##### Function to Plot GR with ML b-value and poisson confidence interval
## Inputs: Magnitudes, Mc estimate, bin size
## Outputs: A plot of (non-cumulative) GR distribution
plot_GR_F <- function(mag, Mc, mbin){
  Z <- hist(mag, plot=FALSE, breaks=seq(min(mag), max(mag)+2, mbin))
  ind <- which(mag >= Mc - mbin/2)
  b <- log10(exp(1))/(mean(mag[ind])-(Mc-mbin/2))
  #print(paste0("Mc estimate (maxc): ", Mc, " ML b-value estimate: ", b))
  rate.mle <- log(10)*b
  m <- which( Z$counts>0 )
  plot(Z$mids[m], log10(Z$counts[m]) , pch=2, col=2, cex=0.4, xlab="Magnitude",ylab="log10(Count)", xlim=c(min(Z$mids),max(Z$mids)), main=paste0("Mc = ", Mc, " , b = ", b))
  ### Because you basically need an idea of how big a "complete" dataset is before you can rescale
  a <- log10(length(ind)) + b*(Mc - mbin/2)
  N <- 10^(a-(b*min(abs(Z$mids))))
  #print(paste0("a estimate: ", a))
  #print(N)
  C <- dexp(Z$mids, rate.mle)*mbin*N
  #print(C)
  lines(Z$mids, log10(C), lty=1, col="gray")
  qhi <- 0.975
  qlo <- 0.025
  nhi <- qpois(qhi,C)
  nlo <- qpois(qlo,C)
  points(Z$mids, log10(nhi), lty=2, type="l")
  points(Z$mids, log10(nlo), lty=2, type="l")
  abline(v=Mc, lty=2, col=adjustcolor("blue", alpha=1))
}

####Function to calculate GR parameters a and b
## Inputs: Magnitudes (vector), Mc estimate, bin size
## Outputs: list containing maximum likelihood a and b-value estimates for GR FMD 
calc_GR_params <- function(mag, Mc, mbin) {
  indmag <- which(mag > Mc-mbin/2)
  b <- log10(exp(1))/(mean(mag[indmag])-(Mc-mbin/2))
  a <- log10(length(indmag))+b*(Mc - mbin/2)
  params <- c(a, b)
  return(params)
}

###### Convert datatimes from INGV format to datetimes format
## Input: INGV time column (date and time separated by T)
## Output: INGV datetimes as POSIX datetime in UTC
INGV_datatimes <- function(INGV_time_col){
  DT <- as.data.frame(matrix(unlist(strsplit(INGV_time_col, split="T")), nrow=length(INGV_time_col), byrow=T))
  options(digits.secs = 3)
  INGV_datetime <- as.POSIXct(paste0(DT$V1, " ", DT$V2), tz="UTC")
}

##### Function to Plot GR with ML b-value and poisson confidence interval
## Inputs: Magnitudes, Mc estimate, bin size
## Outputs: A plot of (cumulative) GR distribution
plot_GR_cumulative <- function(mag, Mc, mbin, legend=TRUE, leg_loc = "topright"){
  pars <- calc_GR_params(mag, Mc, mbin)
  fmd <- fmd(mag, mbin)
  plot(fmd$m, log10(fmd$cum), xlab="magnitude", ylab= "log10(event count)", pch=2, col=2, cex=0.4)
  abline(a=pars[1], b=-pars[2], lty=1, col="gray")
  ct <- 10^(pars[1] - fmd$m*pars[2])
  #lines(fmd$m, ct,lty=1, col="gray" )
  qhi <- 0.975
  qlo <- 0.025
  nhi <- qpois(qhi,ct)
  nlo <- qpois(qlo,ct)
  points(fmd$m, log10(nhi), lty=2, type="l")
  points(fmd$m, log10(nlo), lty=2, type="l")
  abline(v=Mc, lty=2, col=adjustcolor("blue", alpha=1))
  if(legend == TRUE){
    legend(leg_loc, lty=c(1, 2,2), col=c("gray", "black", "blue"), legend=c(paste0("ML b estimate = ", round(pars[2], digits=2)), "Poisson confidence interval", paste0("Mc = ", Mc)))
    
  }
}

#### Choose Mc with Nick's method
## Roberts et al 2015 (https://www.sciencedirect.com/science/article/abs/pii/S0377027315003509) decsribe a method for choosing Mc based on the three Mc methods above
## This function basically implements the workflow in Fig 9 to recommend a suitable Mc estimate given the data
## Input: magnitudes, bin size
## Output: Prints
Choose_Mc <- function(m, mbin){
  maxc_mc <- maxc(m, mbin)$Mc
  gft_mc <- gft(m, mbin)$Mc
  mbs_mc <- mbs(m, mbin)$Mc
  
  diff_gft_mbs <- abs(gft_mc - mbs_mc)
  diff_gft_maxc <- abs(gft_mc - maxc_mc)
  if (diff_gft_maxc <= 0.1 & diff_gft_mbs <= 0.1){
    
    ind <- which(m >= maxc_mc - mbin/2)
    nbev <- length(ind)
    maxc_b <- log10(exp(1))/(mean(m[ind])-(maxc_mc-mbin/2))
    maxc_unc <- 2.3*maxc_b^2*sqrt(sum((m[ind]-mean(m[ind]))^2)/(nbev*(nbev-1)))
    if (maxc_unc < 0.25){
      print("All similar - use maxc")
      mc <- maxc_mc
    }
    else{
      ind <- which(m >= mbs_mc - mbin/2)
      nbev <- length(ind)
      mbs_b <- log10(exp(1))/(mean(m[ind])-(mbs_mc-mbin/2))
      mbs_unc <- 2.3*mbs_b^2*sqrt(sum((m[ind]-mean(m[ind]))^2)/(nbev*(nbev-1)))
      
      if(mbs_unc < 0.25){
        mc <- mbs_mc
        print("All similar - use mbs based on uncertainty")
      }
      
      else{
        ind <- which(m >= gft_mc - mbin/2)
        nbev <- length(ind)
        gft_b <- log10(exp(1))/(mean(m[ind])-(gft_mc-mbin/2))
        gft_unc <- 2.3*mbs_b^2*sqrt(sum((m[ind]-mean(m[ind]))^2)/(nbev*(nbev-1)))
        
        if(mbs_unc < 0.25){
          mc <- gft_mc
          print("All similar - use gft based on uncertainty")
        }
        else{
          print("GR model not appropriate")
        }
      }
    }
  }
  else{
    ind <- which(m >= mbs_mc - mbin/2)
    nbev <- length(ind)
    mbs_b <- log10(exp(1))/(mean(m[ind])-(mbs_mc-mbin/2))
    mbs_unc <- 2.3*mbs_b^2*sqrt(sum((m[ind]-mean(m[ind]))^2)/(nbev*(nbev-1)))
    
    if(mbs_unc < 0.25){
      mc <- mbs_mc
      print("Mc not similar but mbs has low uncertainty - use mbs")
    }
    else {
      print("GR model not appropriate")
    }
  }
   return(mc) 
}

########################### Sarah Touati's ETAS simulation code
########### Code for ETAS simulation from Sarah's thesis
####### You will need all of these functions for S-T ETAS simulation and to call the code with something like:
##### params=c(0.01, 5, 1, 0.01, 1.2, 2, 1.5)
#### coords2 <- c(-120, -115, 32, 36)
#### ETASevents1 <- etas.sim.spatial(params, coords2, 10000, seed=5, bvalue=1) 
####

etas.sim.spatial <- function (params, coords, max.events, seed = 5, 
                              bvalue = 1) 
{
  set.seed(seed)
  
  
  # work out the temporal length we should run each aftershock sequence
  fraction <- 0.05 ## this is the fraction of aftershocks we want
  ## each sequence to be short by, on average
  cc <- params[4]
  theta <- params[5] - 1
  sequence.length <- cc * fraction^(-1/theta)
  
  
  # work out the spatial length each sequence should be on average
  fraction <- 0.05 ## this is the fraction of aftershocks we want
  ## each sequence to be short by, on average
  dd <- params[6]
  qq <- params[7]
  sequence.spatial.length <- dd * sqrt(fraction^(-1/(qq-1)) - 1)
  
  
  # create background sequence to pass into simulation, estimating the
  # required number of background events conservatively
  bg.times <- c()
  bg.magnitudes <- c()
  bg.xs <- c()
  bg.ys <- c()
  bg.gas.indices <- c()
  mu <- params[1]
  if(length(coords)==2)
  {
    centre <- coords[1]
    radius <- coords[2]
  }
  else
  {
    x1 <- coords[1]
    x2 <- coords[2]
    y1 <- coords[3]
    y2 <- coords[4]
  }
  index <- 1
  ti <- 0
  n <- calc.n(params,bvalue)
  repeat
  {
    # increment ti
    tau <- rexp(1, rate = mu)
    ti <- ti + tau
    
    # magnitude
    magnitude <- rexp(1, bvalue * log(10))
    while (magnitude  < 0.1){magnitude <- rexp(1, bvalue * log(10))}
    #print(magnitude)
    #if (magnitude > 0.1)
    #  break
    
    # coords
    reject <- FALSE
    if(length(coords)==2)
    {
      x <- runif(1, (centre-radius-2*sequence.spatial.length), 
                 (centre+radius+2*sequence.spatial.length))
      y <- runif(1, (centre-radius-2*sequence.spatial.length), 
                 (centre+radius+2*sequence.spatial.length))
      if(sqrt((x-centre)^2 + (y-centre)^2) > 
         (radius + 2*sequence.spatial.length)) reject <- TRUE
    }
    else
    {
      
      y <- runif(1, (y1 - ((2*sequence.spatial.length) / 6378) * (180 / pi)), 
                 (y2 + ((2*sequence.spatial.length) / 6378) * (180 / pi)))
      x <- runif(1, (x1 - ((2*sequence.spatial.length) /6378) * (180 / pi) / cos(y * pi/180)), 
                 (x2 + ((2*sequence.spatial.length) /6378) * (180 / pi)/ cos(y * pi/180)))
    }
    
    if(!reject)
    {
      # add the event to bg.events
      bg.times <- c(bg.times,ti)
      bg.magnitudes <- c(bg.magnitudes,magnitude)
      bg.xs <- c(bg.xs,x)
      bg.ys <- c(bg.ys,y)
      bg.gas.indices <- c(bg.gas.indices,as.character(index))
      
      index <- index + 1
      
      # check if we have (more than) enough independent events
      if(length(bg.times) >= max.events/(1+n)) break
    }
  }
  
  
  # set the background rate (mu) to 0
  params[1] <- 0
  
  
  # create the aftershocks
  all.events <- create.aftershocks.spatial(times=bg.times, 
                                           magnitudes=bg.magnitudes, xs=bg.xs, ys=bg.ys, 
                                           gas.indices=bg.gas.indices, params=params, coords=coords, 
                                           bvalue=bvalue, sequence.length=sequence.length, start.time=0, 
                                           max.events=max.events)
  
  
  # extract the data
  times <- all.events$times
  magnitudes <- all.events$magnitudes
  xs <- all.events$xs
  ys <- all.events$ys
  gas.indices <- all.events$gas.indices
  remove(all.events)
  
  
  # order the events chronologically
  ii <- order(times)
  times <- times[ii]
  magnitudes <- magnitudes[ii]
  xs <- xs[ii]
  ys <- ys[ii]
  gas.indices <- gas.indices[ii]
  
  
  # truncate data to required length
  max.time <- Inf
  if(length(times) > max.events) max.time <- times[max.events]
  use <- times <= max.time
  times <- times[use]
  magnitudes <- magnitudes[use]
  xs <- xs[use]
  ys <- ys[use]
  gas.indices <- gas.indices[use]
  
  
  return.object <- list(time=times, magnitude=magnitudes, x=xs, y=ys, 
                        gas.indices=gas.indices)
  sim <- as.data.frame(return.object)
  return(sim)
}

calc.n <- function(params, bvalue, m.max=NULL, m0=NULL)
{
  beta <- bvalue*log(10)
  A <- params[2]
  alpha <- params[3]
  c <- params[4]
  p <- params[5]
  
  n <- A*c/(p-1) * beta/(beta-alpha)
  if(!is.null(m.max))
  {
    n <- n * 
      (1-exp(-(beta-alpha)*(m.max-m0)))/(1-exp(-beta*(m.max-m0)))
  }
  
  return(n)
}

create.aftershocks.spatial <- function(times, magnitudes, xs, ys, 
                                       gas.indices, params, coords, bvalue, sequence.length, start.time, 
                                       max.events)
{
  if(length(coords)==2)
  {
    centre <- coords[1]
    radius <- coords[2]
  }
  else
  {
    x1 <- coords[1]
    x2 <- coords[2]
    y1 <- coords[3]
    y2 <- coords[4]
  }
  
  # objects to eventually hold all events (those passed in to the 
  # function, and aftershocks created for them)
  new.times <- c()
  new.magnitudes <- c()
  new.xs <- c()
  new.ys <- c()
  new.gas.indices <- c()
  
  num.mainshocks <- 0
  for(i in 1:length(times))
  {
    # select an event
    parent.time <- times[i]
    parent.mag <- magnitudes[i]
    parent.x <- xs[i]
    parent.y <- ys[i]
    parent.index <- gas.indices[i]
    
    # make its time the starting time
    ti <- parent.time
    
    # if this is a b/g event (no dot in index), check if the number of
    # events collected so far is enough yet
    if(i > 1 & length(grep("[0-9]+[.][0-9].*",parent.index))==0)
    {
      if(length(coords)==2)
      {
        new.rs <- sqrt((new.xs-centre)^2 + (new.ys-centre)^2)
        use <- new.rs <= radius
      }
      else use <- 
          new.xs < x2 & new.xs >= x1 & new.ys < y2 & new.ys >= y1
      num.aftershocks <- length(new.times[use])
      if((num.aftershocks + num.mainshocks) > max.events)
      {
        # we have enough events, so we can finish here.
        end.time <- times[i-1]
        
        # throw away unused mainshocks
        times <- times[1:(i-1)]
        magnitudes <- magnitudes[1:(i-1)]
        xs <- xs[1:(i-1)]
        ys <- ys[1:(i-1)]
        gas.indices <- gas.indices[1:(i-1)]
        
        # remove events outside boundaries
        new.times <- new.times[use]
        new.magnitudes <- new.magnitudes[use]
        new.xs <- new.xs[use]
        new.ys <- new.ys[use]
        new.gas.indices <- new.gas.indices[use]
        if(length(coords)==2)
        {
          rs <- sqrt((xs-centre)^2 + (ys-centre)^2)
          use <- rs <= radius
        }
        else use <- xs < x2 & xs >= x1 & ys < y2 & ys >= y1
        times <- times[use]
        magnitudes <- magnitudes[use]
        xs <- xs[use]
        ys <- ys[use]
        gas.indices <- gas.indices[use]
        
        # wrap aftershocks occurring after the end time around to 
        # the beginning, and reset the first-generation part of 
        # their index to break the connection with the later ones
        event.indices <- which(new.times > end.time)
        diffs <- new.times[event.indices] - end.time
        sim.length <- end.time - start.time
        wraps <- floor(diffs/sim.length)
        all.parent.indices <- as.numeric(sub("[.].*", "", 
                                             new.gas.indices[event.indices]))
        to.add <- wraps * max(all.parent.indices)
        new.parent.indices <- all.parent.indices + to.add
        for(j in 1:length(new.parent.indices))
        {
          new.gas.indices[event.indices[j]] <- sub("^[0-9]+", 
                                                   new.parent.indices[j], 
                                                   new.gas.indices[event.indices[j]])
        }
        new.times[event.indices] <- new.times[event.indices] - 
          (wraps * sim.length)
        
        # include mainshocks in returned events
        new.times <- c(times, new.times)
        new.magnitudes <- c(magnitudes, new.magnitudes)
        new.xs <- c(xs, new.xs)
        new.ys <- c(ys, new.ys)
        new.gas.indices <- c(gas.indices, new.gas.indices)
        
        final.events <- list(times=new.times, 
                             magnitudes=new.magnitudes, xs=new.xs, ys=new.ys, 
                             gas.indices=new.gas.indices)
        return(final.events)
      }
    }
    
    # increment the number of mainshocks if this mainshock is within 
    # the spatial boundaries
    if(length(coords)==2)
      if(sqrt((parent.x-centre)^2 + (parent.y-centre)^2) <= radius) 
        num.mainshocks <- num.mainshocks + 1
    else if(length(coords)==4)
      if(parent.x < x2 & parent.x >= x1 & parent.y < y2 & 
         parent.y >= y1) 
        num.mainshocks <- num.mainshocks + 1
    
    # objects to hold its aftershocks
    aftershock.times <- c()
    aftershock.magnitudes <- c()
    aftershock.xs <- c()
    aftershock.ys <- c()
    aftershock.gas.indices <- c()
    
    # calculate the initial rate, at the time of this event
    Rmax <- conditional.intensity(data.mag=parent.mag, 
                                  data.time=parent.time, eval.time=ti, params=params)
    repeat
    {
      # increment the time by an appropriate amount, tau
      if(Rmax > 0) tau <- rexp(1, rate = Rmax)
      else tau <- Inf
      ti <- ti + tau
      
      # check if this sequence is long enough yet
      if ((ti-parent.time) > sequence.length) break
      
      # calculate the new rate
      rate <- conditional.intensity(data.mag=parent.mag, 
                                    data.time=parent.time, eval.time=ti, params=params)
      #print(rate)
      #print(Rmax)
      # decide whether to create an aftershock at this time 
      # (thinning method)
      if (runif(1, 0, 1) <= rate/Rmax)
      {
        # select a magnitude
        new.mag <- rexp(1, bvalue * log(10))
        while (new.mag  < 0.1){new.mag <- rexp(1, bvalue * log(10))}
        #if (new.mag > 0.1){break} 
        
        # select a distance
        dd <- params[6]
        qq <- params[7]
        u <- runif(1, 0, 1)
        root <- sqrt((1 - u)^(-1/(qq-1)) - 1)
        distance <- dd*root
        
        
        # select an orientation and work out the event coordinates
        angle <- runif(1, 0, 2*pi)
        #((2*sequence.spatial.length*1000) /6378) * (180 / pi)
        
        new.y <- parent.y + (distance*cos(angle)/6378)*(180/pi)
        new.x <- parent.x + (distance*sin(angle)/6378)*(180/pi)/cos(new.y*pi/180)
        # make the index the same as the parent's but with ".i" 
        # appended, where i is the aftershock number within this 
        # sequence
        new.index <- paste(parent.index, ".", 
                           length(aftershock.times), sep="")
        
        # add the aftershock to the data objects
        aftershock.times <- c(aftershock.times, ti)
        aftershock.magnitudes <- c(aftershock.magnitudes, new.mag)
        aftershock.xs <- c(aftershock.xs, new.x)
        aftershock.ys <- c(aftershock.ys, new.y)
        aftershock.gas.indices <- 
          c(aftershock.gas.indices, new.index)
      }
      
      # re-set the initial rate for the next iteration
      Rmax <- rate 
    }
    
    
    if(length(aftershock.times) > 0) 
    {
      #counter <<- counter + 1
      #print(counter)
      # create aftershock sequences for each of these aftershocks
      aftershocks <- 
        create.aftershocks.spatial(times=aftershock.times, 
                                   magnitudes=aftershock.magnitudes, xs=aftershock.xs, 
                                   ys=aftershock.ys, gas.indices=aftershock.gas.indices, 
                                   params=params, coords=coords, bvalue=bvalue, 
                                   sequence.length=sequence.length, start.time=start.time, 
                                   max.events=max.events)
      remove(aftershock.times)
      remove(aftershock.magnitudes)
      remove(aftershock.xs)
      remove(aftershock.ys)
      remove(aftershock.gas.indices)
      
      # add the created aftershocks to our collection of events
      new.times <- c(new.times, aftershocks$time)
      new.magnitudes <- c(new.magnitudes, aftershocks$magnitude)
      new.xs <- c(new.xs, aftershocks$xs)
      new.ys <- c(new.ys, aftershocks$ys)
      new.gas.indices <- c(new.gas.indices, aftershocks$gas.indices)
    }
  }
  
  # include parent events in returned events
  new.times <- c(times, new.times)
  new.magnitudes <- c(magnitudes, new.magnitudes)
  new.xs <- c(xs, new.xs)
  new.ys <- c(ys, new.ys)
  new.gas.indices <- c(gas.indices, new.gas.indices)
  
  events <- list(times=new.times, magnitudes=new.magnitudes, xs=new.xs, 
                 ys=new.ys, gas.indices=new.gas.indices)
  return(events)
}

conditional.intensity <- function(data.mag, data.time, eval.time, params) 
{
  mu <- params[1]
  A <- params[2]
  alpha <- params[3]
  CC <- params[4]
  P <- params[5]
  
  if (length(data.time) > 0) 
  {
    triggering.ci <- A * sum(exp(params[3] * data.mag) * 
                               (1 + (eval.time - data.time)/CC)^(-P))
  }
  else triggering.ci <- 0
  ci <- mu + triggering.ci
  return(ci)
}

