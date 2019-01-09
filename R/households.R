#' calculate the probability density for a given incubation period
#' 
#' @param d value of incubation at which to calculate probability density.
#' @param p named vector of parameter values.  Contains element "m.incub", the mean
#' incubation period
#' @return the probability density of the incubation period at d given model
#' parameters
#' @export
f.incub <- function(d, p) {
  # We assume that the incubation period has an exponential distribution
  return( exp(-d / p[["m.incub"]]) / p[["m.incub"]] )
}

#' calculate the hazard of infection at time t
#' 
#' @param t time at which to calculate the hazard of infection.
#' @param p named vector of parameter values.  Contains elements "beta", the
#' transmission parameter; and "m.inf", the mean infectious period
#' @return the hazard of infection at time t given model parameters
#' @export
h <- function(t,p) {
  return( p[["beta"]]*exp(-t/p[["m.inf"]])/p[["m.inf"]] )
}

#' calculate the cumulative hazard of infection at time t
#' 
#' @param t time at which to calculate the hazard of infection.
#' @param p named vector of parameter values.  Contains elements "beta", the
#' transmission parameter; and "m.inf", the mean infectious period
#' @return the cumulative hazard of infection at time t given model parameters
#' @export
H <- function(t,p) {
  return( p[["beta"]]*(1-exp(-t/p[["m.inf"]])) )
}

#' calculate the probability of having survived infection at time t
#' 
#' @param t time at which to calculate the hazard of infection.
#' @param p named vector of parameter values.  Contains elements "beta", the
#' transmission parameter; and "m.inf", the mean infectious period
#' @return the probability of having survived infection at time t given model parameters
#' @export
S <-function(t,p) {
  return( exp(-H(t,p)) )
}

#' calculate the contribution of a cluster to the log likelihood, for the model
#' of incubation only
#' 
#' @param p named vector of parameter values.  Contains element "m.incub", the mean
#' incubation period
#' @param data data frame with a single row, containing the column "time.onset.donor"
#' which is the time of symptom onset for the donor, i.e. the incubation period
#' for the donor.
#' @return the contribution of the cluster to the log likelihood
#' @export
calculateLogLikCluster_incubation<-function(p, data)	{
  LL<-log(f.incub(data[["time.onset.donor"]], p))
  return(LL)
}

#' calculate log likelihood for the model of incubation only, summed across clusters
#' 
#' @param p named vector of parameter values.  Contains element "m.incub", the mean
#' incubation period
#' @param data data frame where each row represents a donor-recipient pair.
#' It contains the column "time.onset.donor" which is the time of symptom onset 
#' for the donor, i.e. the incubation period for the donor.
#' @return the log likelihood
#' @export
calculateLogLik_incubation<-function(p, data)	{
  logLik.cluster<-apply(data, 1, calculateLogLikCluster_incubation,p=p)
  return( sum(logLik.cluster) )
}

#' calculate the contribution of a cluster to the log likelihood, for the model
#' of incubation and transmission
#' 
#' @param p named vector of parameter values.  Contains elements "m.incub", the mean
#' incubation period;"beta", the transmission parameter; and "m.inf", the mean 
#' infectious period
#' @param data data frame with a single row, containing the columns 
#' "time.onset.donor": the observed time of symptom onset for the donor;
#' "recipient.infected": 1 if the recipient was infected, 0 if not;
#' "time.infection.recipient": the inferred time at which the recipient was infected;
#' "time.onset.recipient": the observed time of symptom onset for the recipient
#' @return the contribution of the cluster to the log likelihood
#' @export
calculateLogLikCluster_transmission <-function(p, data) {
  # contribution to log likelihood from the observed incubation period of the donor
  LL<-log(f.incub(data[["time.onset.donor"]],p))
  if(data[["recipient.infected"]]==1) {	#for recipients that were infected
    # contribution to log likelihood from surviving until the time of infection
    LL<-LL+log(S(data[["time.infection.recipient"]],p)) +
      # contribution to log likelihood from becoming infected at the time of infection
      log(h(data[["time.infection.recipient"]],p)) +
      # contribution to log likelihood from the inferred incubation period of 
      # the recipient
      log(f.incub(data[["time.onset.recipient"]]-data[["time.infection.recipient"]],p))
  } else  {
    # contribution to log likelihood from surviving until the end of the experiment
    LL<-LL+log(S(1000,p))
  }
  return(LL)
}

#' run MCMC algorithm to fit parameters to data
#' 
#' @param fit_parnames vector of names of parameter(s) to be fitted
#' @param nbIteration number of iterations for which to run MCMC
#' @param randomWalkSD vector containing standard deviations of proposal
#' distributions for each fitted parameter
#' @param currentParameter named vector of parameter values at which to start fitting
#' @param data data for fitting, where each row is a cluster
#' @return list with the elements
#' storedParameter: matrix of fitted parameter values for each iteration
#' logLik: log likelihood for each iteration
#' acceptance_rate: proportion of proposed values that were accepted, for each
#' fitted parameter
#' @export
run_MCMC <- function(fit_parnames, 
                     nbIteration,
                     randomWalkSD,
                     currentParameter,
                     data) {
  
  nbParameter <- length(fit_parnames)
  storedParameter<-matrix(0,ncol=nbParameter,nrow=nbIteration)	#parameter values at the different iterations of the MCMC are stored in a matrix
  colnames(storedParameter) <- fit_parnames
  logLik<-rep(0,nbIteration)	#log likelihood at the different iterations of the MCMC is stored in a vector
  
  nbAccepted<-rep(0,nbParameter)		#for each parameter, number of updates that are accepted
  names(nbAccepted) <- fit_parnames
  
  currentLogLik<-calculateLogLik_incubation(currentParameter, data)
  
  storedParameter[1,]<-currentParameter	#store current parameter value
  logLik[1]<-currentLogLik
  
  par_ind <- match(fit_parnames, names(currentParameter))
  
  for(iter in 2:nbIteration) {
    #-----updating each parameter independently
    for(parID in par_ind) {
      oldValue<-currentParameter[parID]
      newValue<-oldValue*exp(randomWalkSD[parID]*rnorm(1))
      currentParameter[parID]<-newValue
      newLogLik<-calculateLogLik_incubation(currentParameter, data)
      
      logRatioProposal<-log(newValue)-log(oldValue)
      differenceLogLik<-newLogLik-currentLogLik
      if(log(runif(1))<differenceLogLik+logRatioProposal)  {	#candidate value is accepted
        currentLogLik<-newLogLik
        nbAccepted[parID]<-nbAccepted[parID]+1
      }
      else	{
        currentParameter[parID]<-oldValue
      }	#candidate value is not accepted 
    }
    
    #-----store data 
    storedParameter[iter,]<-currentParameter
    logLik[iter]<-currentLogLik
  }
  acceptance_rate <- nbAccepted/nbIteration	#acceptance rate	
  
  return(list(storedParameter = storedParameter,
              logLik = logLik,
              acceptance_rate = acceptance_rate))
}

#' run MCMC with data augmentation to fit parameters to data
#' 
#' @param fit_parnames vector of names of parameter(s) to be fitted
#' @param nbIteration number of iterations for which to run MCMC
#' @param randomWalkSD vector containing standard deviations of proposal
#' distributions for each fitted parameter
#' @param currentParameter named vector of parameter values at which to start fitting
#' @param data data for fitting, where each row is a cluster
#' @return list with the elements
#' storedParameter: matrix of fitted parameter values for each iteration
#' logLik: log likelihood for each iteration
#' acceptance_rate_parameters: proportion of proposed values that were accepted,
#' for each fitted parameter
#' acceptance_rate_DA: proportion of proposed data augmentation values that were accepted
#' @export
run_MCMC_DA <- function(fit_parnames,
                        nbIteration,
                        randomWalkSD,
                        currentParameter,
                        data) {
  
  n.cluster <- nrow(data)
  nbParameter <- length(fit_parnames)
  
  storedParameter<-matrix(0,ncol=nbParameter,nrow=nbIteration)	#parameter values at the different iterations of the MCMC are stored in a matrix
  colnames(storedParameter) <- fit_parnames
  logLik<-rep(0,nbIteration)	#log likelihood at the different iterations of the MCMC is stored in a vector
  
  nbAccepted<-rep(0,nbParameter)		#for each parameter, number of updates that are accepted
  names(nbAccepted) <- fit_parnames
  nAcceptedUpdateTime<-0			#number of updates of the time of infection that are accepted
  
  data <- t(apply(data, 1, augment_data))
  currentLogLik.cluster<-apply(data,1, calculateLogLikCluster_transmission,
                               p = currentParameter) #vector of log-likelihood for each cluster
  currentLogLik<-sum(currentLogLik.cluster)
  
  storedParameter[1,]<-currentParameter	#store current parameter value
  logLik[1]<-currentLogLik
  
  par_ind <- match(fit_parnames, names(currentParameter))
  
  for(iter in 2:nbIteration) {
    #-----updating each parameter independently
    for(parID in par_ind) {
      oldValue<-currentParameter[parID]
      newValue<-oldValue*exp(randomWalkSD[parID]*rnorm(1))
      currentParameter[parID]<-newValue
      newLogLik.cluster<-apply(data,1, calculateLogLikCluster_transmission,
                              p = currentParameter)
      newLogLik<-sum(newLogLik.cluster)
      
      logRatioProposal<-log(newValue)-log(oldValue)
      differenceLogLik<-newLogLik-currentLogLik
      if(log(runif(1))<differenceLogLik+logRatioProposal) {	#candidate value is accepted
        currentLogLik<-newLogLik
        currentLogLik.cluster<-newLogLik.cluster
        nbAccepted[parID]<-nbAccepted[parID]+1
      }
      else	{
        currentParameter[parID]<-oldValue
      }	#candidate value is not accepted
    }
    
    #-----update the time of infection in each infected cluster
    
    DA_flag <- apply(data, 1, DA_condition)
    
    for(cluster in 1:n.cluster)    {
      if(DA_flag[cluster]) {
        oldValues <- data[cluster,]
        newValues <- augment_data(data[cluster,])
        data[cluster,] <- newValues
        newLogLikOfCluster<-calculateLogLikCluster_transmission(p=currentParameter,
                                                      data = data[cluster,])
        
        logRatioProposal<- logRatioProposal_DA(newValues, oldValues)
        differenceLogLik<-newLogLikOfCluster-currentLogLik.cluster[cluster]
        if(log(runif(1))<differenceLogLik+logRatioProposal) {	#candidate value is accepted
          currentLogLik<-currentLogLik+differenceLogLik
          currentLogLik.cluster[cluster]<-newLogLikOfCluster
          nAcceptedUpdateTime<-nAcceptedUpdateTime+1
        }
        else	{
          data[cluster, ]<-oldValues
        }	#candidate value is not accepted
      }
    }
    
    #--------store data in matrice
    storedParameter[iter,]<-currentParameter
    logLik[iter]<-currentLogLik
  }
  
  acceptance_rate_parameters <- nbAccepted/nbIteration	#acceptance rate	
  acceptance_rate_DA <- 
    nAcceptedUpdateTime / nbIteration / sum(DA_flag)
  
  return(list(storedParameter = storedParameter,
              logLik = logLik,
              acceptance_rate_parameters = acceptance_rate_parameters,
              acceptance_rate_DA = acceptance_rate_DA))
}

#' function to augment data
#' 
#' @param data data frame with a single row, containing the columns 
#' "time.onset.donor": the observed time of symptom onset for the donor;
#' "recipient.infected": 1 if the recipient was infected, 0 if not;
#' "time.infection.recipient": the inferred time at which the recipient was infected;
#' "time.onset.recipient": the observed time of symptom onset for the recipient
#' @return the same data frame but with a newly generated value for time.infection.recipient
#' @export
augment_data <- function(data) {
  newValues <- data
  newValues[["time.infection.recipient"]] <-
    runif(1,
          data[["time.infection.donor"]],
          data[["time.onset.recipient"]])
  newValues
}

#' function to decide whether data should be augmented for a cluster
#' 
#' @param data data frame with a single row, containing the columns 
#' "time.onset.donor": the observed time of symptom onset for the donor;
#' "recipient.infected": 1 if the recipient was infected, 0 if not;
#' "time.infection.recipient": the inferred time at which the recipient was infected;
#' "time.onset.recipient": the observed time of symptom onset for the recipient
#' @return TRUE if the data should be augmented, FALSE otherwise
#' @export
DA_condition <- function(data) {
  data[["recipient.infected"]] == 1
}

#' function to return the log ratio of the proposal distribution at the new and
#' old values, for the data augmentation
#' @param newValues proposed values for augmented data
#' @param oldValues old values for augmented data
#' @return the log ratio of the proposal distribution at the new and
#' old values
#' @export
logRatioProposal_DA <- function(newValues, oldValues) {
  return(0)
}