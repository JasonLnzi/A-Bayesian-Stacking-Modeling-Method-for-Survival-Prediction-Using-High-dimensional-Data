
bs.dat.t <- function(cross.fwrd.dat, t) {
  
  ys = cross.fwrd.dat$time; ds = cross.fwrd.dat$status
  
  Zi <- 1 * {ys > t}														# Response in the weighted LS problem
  Di <- 1 * {{ds == 1} | {ys >= t}}             # delta i
  
  # Create survival curves
  # Adjusted Delta Function at time point t
  modG <- survfit(Surv(ys, I(1-ds)) ~ 1, type = "kaplan-meier")
  tmpSurv <- rep(modG$surv, modG$n.event + modG$n.censor)	# Compute censoring survival at each observed time point (allowing for possible ties)
  
  Gc <- unlist(sapply(1:length(ys), function(x){                                                    	# Define function to evaluate censoring survival at each time point relative to t
    
    if(mean(ys <= t) != 1){ tmp1 <- {ys[x] <= t} * tmpSurv[x] + {ys[x] > t} * tmpSurv[ys > t][1] }  	# Calculate censoring survival for observations occurring prior to time t
    if(mean(ys <= t) == 1){ tmp1 <- tmpSurv[x] }                                                    	# Calculate censoring survival for observations occurring after time t
    
    return(tmp1) }))
  
  Wi <- Di / Gc		# Compute weights for the weighted LS problem
  
  # Remove any observations with weight 0 or infinity...
  booRemove <- ! Wi %in% c(0, NaN, Inf)   
  times = ys[booRemove==TRUE]; 
  times[times > t] = t; 
  

  rslt <- cbind(W = Wi[booRemove], Z = Zi[booRemove], T = times)	# Retain weights, responses, times and survival esitmates
  rslt = data.frame(rslt); names(rslt) = c("weight", "Zind", "tsubr")
  
  out <- list(rslt = rslt, booRemove = booRemove)
  # return weight, Zind, and Evaluate survival for cv.grps at time t
  return(out)
  
}
