network.inits <- function(network, n.chains){

  response <- network$response
  
  inits <- if(response == "multinomial"){
    multinomial.inits(network, n.chains)
  } else if(response == "binomial"){
    binomial.inits(network, n.chains)
  } else if(response == "normal"){
    normal.inits(network, n.chains)
  }
  return(inits)
}

normal.inits <- function(network, n.chains){

  with(network,{
    Eta <- Outcomes[b.id]
    se.Eta <- SE[b.id]
    delta <- Outcomes - rep(Eta, times = na)
    delta <- delta[!b.id,] #eliminate base-arm
 
    inits <- make.inits(network, n.chains, delta, Eta, se.Eta)
    return(inits)
  })
}

binomial.inits <- function(network, n.chains){

  with(network,{
  
    Outcomes <- Outcomes + 0.5 # ensure ratios are always defined
    N <- N + 1
    p <- Outcomes/N
    logits <- log(p/(1-p))
    se.logits <- sqrt(1/Outcomes + 1/(N - Outcomes))

    Eta <- logits[b.id]
    se.Eta <- se.logits[b.id]
    delta <- logits - rep(Eta, times = na)
    delta <- delta[!b.id,]
    
    inits = make.inits(network, n.chains, delta, Eta, se.Eta)
    return(inits)  
  })
}

make.inits <- function(network, n.chains, delta, Eta, se.Eta){

  with(network,{
    
  # dependent variable for regression
  y <- delta

  # design matrix
  base.tx <- Treat[b.id]    # base treatment for N studies
  end.Study <- c(0, cumsum(na))  # end row number of each trial
  rows <- end.Study - seq(0, nstudy)   # end number of each trial not including base treatment arms
  design.mat <- matrix(0, sum(na) - nstudy, ntreat) # no. non-base arms x #txs
  for (i in seq(nstudy)){
    studytx <- Treat[(end.Study[i]+1):end.Study[i+1]]  #treatments in ith Study
    nonbase.tx <- studytx[studytx!=base.tx[i]]    #non-baseline treatments for ith Study
    design.mat[(1+rows[i]):rows[i+1],base.tx[i]] <- -1
    for (j in seq(length(nonbase.tx)))
      design.mat[j+rows[i],nonbase.tx[j]] <- 1
  }
  design.mat <- design.mat[,-1,drop=F]

  fit <- summary(lm(y ~ design.mat - 1))
  print(summary(fit))
  d <- se.d <- rep(NA, ntreat)
  d[-1] <- coef(fit)[,1]
  se.d[-1] <- coef(fit)[,2]
  resid.var <- fit$sigma^2

  # covariate
  if(!is.null(covariate)) {
    x.cen = matrix(0, nrow = sum(na), ncol = dim(covariate)[2])
    for(i in 1:dim(covariate)[2]){
      x.cen[,i] <- rep(covariate[,i], times = na)
    }
    x.cen <- x.cen[-seq(dim(x.cen)[1])[b.id],,drop=F]
    x.cen <- scale(x.cen, scale = FALSE)

    slope <- se.slope <- array(NA, c(ntreat, dim(covariate)[2]))
    for(i in 1:dim(covariate)[2]){
      fit2 <- if(covariate.model == "common" || covariate.model == "exchangeable"){
        summary(lm(y ~ x.cen[,i] -1))
      } else if(covariate.model == "independent"){
        summary(lm(y ~ design.mat:x.cen[,i] - 1))
      }
      slope[-1,i] <- coef(fit2)[,1]
      se.slope[-1,i] <- coef(fit2)[,2]
    }
  }

  # baseline
  if(baseline != "none"){
    baseline.cen <- rep(Eta, na)
    baseline.cen <- baseline.cen[-seq(length(baseline.cen))[b.id]]
    baseline.cen <- scale(baseline.cen, scale = FALSE)

    baseline.slope <- baseline.se.slope <- rep(NA, ntreat)
    fit3 <- if(baseline == "common" || baseline == "exchangeable"){
      summary(lm(y ~ baseline.cen -1))
    } else if(baseline == "independent"){
      summary(lm(y ~ design.mat:baseline.cen - 1))
    }
    baseline.slope[-1] <- coef(fit3)[,1]
    baseline.se.slope[-1] <- coef(fit3)[,2]
  }

  ############# Generate initial values
  initial.values = list()
  for(i in 1:n.chains){
    initial.values[[i]] = list()
  }
  for(i in 1:n.chains){
    random.Eta <- rnorm(length(Eta))
    initial.values[[i]][["Eta"]] <- Eta + se.Eta * random.Eta
  }

  if(!is.nan(fit$fstat[1])){
    for(i in 1:n.chains){
      random.d = rnorm(length(d))
      initial.values[[i]][["d"]] <- d + se.d * random.d

      if(type == "random"){

        df <- fit$df[2]
        random.ISigma <- rchisq(1, df)
        sigma2 <- resid.var * df/random.ISigma

        if(hy.prior[[1]] == "dunif"){
          if(sqrt(sigma2) > network$prior.data$hy.prior.2){
            stop("data has more variability than your prior does")
          }
        }

        if(hy.prior[[1]] == "dgamma"){
          initial.values[[i]][["prec"]] <- 1/sigma2
        } else if(hy.prior[[1]] == "dunif" || hy.prior[[1]] == "dhnorm"){
          initial.values[[i]][["sd"]] <- sqrt(sigma2)
        }

        # generate values for delta
        delta = matrix(NA, nrow = nrow(t), ncol = ncol(t))
        for(j in 2:ncol(delta)){
          diff_d <- ifelse(is.na(d[t[,1]]), d[t[,j]], d[t[,j]] - d[t[,1]])
          for(ii in 1:nrow(delta)){
            if(!is.na(diff_d[ii])) delta[ii,j] = rnorm(1, mean = diff_d[ii], sd = sqrt(sigma2))
          }
        }
        initial.values[[i]][["delta"]] <- delta
      }
    }
  }

  if (!is.null(covariate)) {
    if(!is.nan(fit2$fstat[1])){
      for(i in 1:n.chains){
        random.slope <- matrix(rnorm(dim(slope)[1]*dim(slope)[2]),dim(slope))
        for(j in 1:dim(covariate)[2]){
          initial.values[[i]][[paste("beta", j, sep = "")]] = slope[,j] + se.slope[,j] * random.slope[,j]
        }
      }
    }
  }

  if(baseline != "none"){
    if(!is.nan(fit3$fstat[1])){
      for(i in 1:n.chains){
        random.baseline = rnorm(length(baseline.slope))
        initial.values[[i]][["b_bl"]] = baseline.slope + baseline.se.slope * random.baseline
      }
    }
  }
  return(initial.values)
  })
  
  
}

############################################ multinomial inits functions

multinomial.inits <- function(network, n.chains)
{
  with(network,{
    
  if (length(miss.patterns[[1]])!= 1){
    Dimputed <- multi.impute.data(network)
  } else{
    Dimputed <- Outcomes
  }
  Dimputed = Dimputed + 0.5

  logits <- as.matrix(log(Dimputed[, -1]) - log(Dimputed[, 1]))
  se.logits <- as.matrix(sqrt(1/Dimputed[, -1] + 1/Dimputed[, 1]))

  Eta <- se.Eta <- matrix(NA, nstudy, ncat)
  Eta[,2:ncat] <- logits[b.id,]
  se.Eta[,2:ncat] <- se.logits[b.id,]

  delta <- logits - apply(as.matrix(Eta[, -1]), 2, rep, times = na)
  rows.of.basetreat <- seq(dim(as.matrix(delta))[1])*as.numeric(b.id)
  delta <- delta[-rows.of.basetreat,,drop=F]   # Eliminate base treatment arms

  ###################### Using delta, Eta, and se.Eta make initial values

  y <- delta            # dependent variable for regression (part of Delta)
  d <- se.d <- matrix(NA, length(unique(Treat)), ncat - 1)
  resid.var <- rep(NA, ncat -1)
  base.tx <- Treat[b.id]    # base treatment for N studies
  end.Study <- c(0, cumsum(na))  # end row number of each trial
  rows <- end.Study - seq(0, nstudy)   # end number of each trial not including base treatment arms
  design.mat <- matrix(0, sum(na) - nstudy, ntreat) # no. non-base arms x #txs

  for (i in seq(nstudy)){
    studytx <- Treat[(end.Study[i]+1):end.Study[i+1]]  #treatments in ith Study
    nonbase.tx <- studytx[studytx!=base.tx[i]]    #non-baseline treatments for ith Study
    design.mat[(1+rows[i]):rows[i+1],base.tx[i]] <- -1
    for (j in seq(length(nonbase.tx)))
      design.mat[j+rows[i],nonbase.tx[j]] <- 1
  }
  design.mat <- design.mat[,-1,drop=F]

  for(k in 1:(ncat - 1)){
    fit <- summary(lm(y[,k] ~ design.mat - 1))
    d[-1,k] <- coef(fit)[1:(ntreat-1), 1]
    se.d[-1,k] <- coef(fit)[1:(ntreat-1), 2]
    resid.var[k] <- fit$sigma^2
  }

  # covariate
  if(!is.null(covariate)){
    x.cen <- matrix(0, nrow = sum(na), ncol = dim(covariate)[2])
    for(i in 1:dim(covariate)[2]){
      x.cen[,i] <- rep(covariate[,i], na)
    }
    x.cen <- x.cen[-seq(dim(x.cen)[1])[b.id],,drop=F]
    x.cen <- scale(x.cen, scale = FALSE)
    slope <- se.slope <- array(NA, c(ntreat, dim(covariate)[2], ncat - 1))

    for(i in 1:dim(covariate)[2]){
      for(k in 1:(ncat-1)){

        fit2 <- if(covariate.model == "independent" || covariate.model == "exchangeable"){
          summary(lm(y[,k] ~ design.mat:x.cen[,i] - 1))
        } else if(covariate.model == "common"){
          summary(lm(y[,k] ~ x.cen[,i] - 1))
        }
        slope[-1,i,k] <- coef(fit2)[,1]
        se.slope[-1,i,k] <- coef(fit2)[,2]
      }
    }
  }

  # baseline
  if(baseline != "none"){
    baseline.cen <- apply(as.matrix(Eta[, -1]), 2, rep, times = na)
    baseline.cen <- baseline.cen[-seq(dim(baseline.cen)[1])[b.id],]
    baseline.cen <- scale(baseline.cen, scale = FALSE)

    baseline.slope <- baseline.se.slope <- matrix(nrow = ntreat, ncol = ncat -1)

    for(k in 1:(ncat -1)){
      fit3 <- if(baseline == "common" || baseline == "exchangeable"){
        summary(lm(y[,k] ~ baseline.cen[,k] -1))
      } else if(baseline == "independent"){
        summary(lm(y[,k] ~ design.mat:baseline.cen[,k] - 1))
      }
      baseline.slope[-1, k] <- coef(fit3)[,1]
      baseline.se.slope[-1, k] <- coef(fit3)[,2]
    }
  }

  ################################################
  initial.values = list()
  for(i in 1:n.chains){
    initial.values[[i]] = list()
  }

  for(i in 1:n.chains){
    random.Eta <- matrix(rnorm(dim(Eta)[1]*dim(Eta)[2]),dim(Eta)[1],dim(Eta)[2])
    initial.values[[i]][["Eta"]] <- Eta + se.Eta * random.Eta
  }

  if(!is.nan(fit$fstat[1])){
    for(i in 1:n.chains){
      random.d = matrix(rnorm(dim(d)[1]*dim(d)[2]),dim(d)[1],dim(d)[2])
      initial.values[[i]][["d"]] = d + se.d * random.d

      if(type == "random"){
        df <- fit$df[2]
        random.ISigma <- rchisq(1, df)
        sigma2 <- resid.var * df/random.ISigma

        initial.values[[i]][["prec"]] <- 1/sigma2 * diag(ncat - 1)

        if(max(na) == 2){
          delta <- array(NA, dim = c(nstudy, max(na), ncat))
          for(j in 2:max(na)){
            for(m in 1:(ncat-1)){
              diff_d <- ifelse(is.na(d[t[,1],m]), d[t[,j],m], d[t[,j],m] - d[t[,1],m])
              for(ii in 1:nstudy){
                if(!is.na(diff_d[ii])) delta[ii,j,m+1] <- rnorm(1, mean = diff_d[ii], sd = sqrt(sigma2))
              }
            }
          }
          initial.values[[i]][["delta"]] <- delta
        }
      }
    }
  }

  if (!is.null(covariate)) {
    if(!is.nan(fit2$fstat[1])){
      for(i in 1:n.chains){
        random.slope <- array(rnorm(dim(slope)[1]*dim(slope)[2]*dim(slope)[3]),dim(slope))
        for(j in 1:dim(covariate)[2]){
          initial.values[[i]][[paste("beta", j, sep = "")]] = slope[,j,] + se.slope[,j,] * random.slope[,j,]
        }
      }
    }
  }

  if(baseline != "none"){
    if(!is.nan(fit3$fstat[1])){
      for(i in 1:n.chains){
        random.baseline = matrix(rnorm(dim(baseline.slope)[1]*dim(baseline.slope)[2]),dim(baseline.slope))
        initial.values[[i]][["b_bl"]] = baseline.slope + baseline.se.slope * random.baseline
      }
    }
  }

  return(initial.values)
  })
}


multi.impute.data <- function(network)
{
  #Take partial sums by study and allocate them to missing outcomes according to either defined category probabilities or to probabilities computed empirically from available data. Empirical scheme first estimates allocation probabilities based on complete data and then updates by successive missing data patterns.
  #
  #1. Fill in all data for all outcomes with complete information
  #2. Pull off summed outcome columns and back out the known data (e.g. if one type of death count known, subtract this from total deaths)
  #3. Renormalize imputation probabilities among outcomes with missing values
  #4. Split summed outcome categories by imputation probabilities for each sum
  #5. For each outcome category, average the imputed values gotten from each partial sum
  #6. Apply correction factor to ensure that sum of imputed values add up to total to be imputed
  #

  with(network,{
  
  rows.all = vector("list", length = npattern)
  for(i in seq(npattern)){
    rows.all[[i]] = seq(nrow)[pattern == levels(pattern)[i]]
  }

  Dimputed = matrix(NA,dim(D)[1],ncat)
  count = 0
  imputed.prop = rep(1/ncat,ncat)

  for (i in seq(length(miss.patterns[[1]]))) {
    rows = rows.all[[i]]        #data rows in missing data pattern
    cols.data = miss.patterns[[1]][[i]][[2]]   #data columns in first combo of missing data pattern
    is.complete.cols = cols.data %in% seq(ncat)   #which data columns are complete
    if (any(is.complete.cols)) {
      complete.cols = cols.data[is.complete.cols] #col no. of complete cols
      incomplete.cols = cols.data[!is.complete.cols] #col nos. of incomplete cols
      Dimputed[rows, complete.cols] = D[rows, complete.cols] #Put in complete data
    }
    else
      incomplete.cols = cols.data
    if (!all(is.complete.cols)) { #If some columns with missing data
      pmat = miss.patterns[[2]][incomplete.cols,,drop=F] #Parameters corresponding to incomplete cols
      if (any(is.complete.cols)) {
        sums.to.split = D[rows, incomplete.cols, drop=F] - D[rows, complete.cols, drop=F]%*%t(pmat[, complete.cols,drop=F]) #back out known data
        pmat[,complete.cols] = 0  #set backed out columns to zero
        imputed.prop[complete.cols] = 0 #set imputation probabilities for complete data cols to zero
      }
      else
        sums.to.split = D[rows, incomplete.cols, drop=F]
      imputed.prop = imputed.prop/sum(imputed.prop)   #renormalize
      for (j in seq(length(rows))) {
        x0 = matrix(rep(sums.to.split [j,], each=ncat),ncol=length(incomplete.cols))*t(pmat)
        x1 = imputed.prop*t(pmat)
        x2 = x0*x1/rep(apply(x1,2,sum),each=ncat,ncol=dim(pmat)[1])
        x2[x2==0] = NA
        x3 = apply(x2, 1, mean, na.rm=T) # average across potential imputed values
        x5 = (N[rows[j]]- sum(Dimputed[rows[j],], na.rm=T))/sum(x3, na.rm=T)  #Factor to adjust imputations
        x6 = round(x3*x5) # Apply factor to imputations
        if (any(is.complete.cols))
          Dimputed[rows[j],seq(ncat)[-complete.cols]] = x6[!is.na(x6)]
        else
          Dimputed[rows[j],seq(ncat)] = x6[!is.na(x6)]
        Dimputed[rows[j],1] = Dimputed[rows[j],1] + N[rows[j]] - sum(Dimputed[rows[j],])  #Correction for rounding so totals add
      }
    }

    running.total = apply(Dimputed,2,sum,na.rm=T)
    imputed.prop = running.total/sum(running.total) # Proportion of events in each category
  }
  return(Dimputed)
  })
}


