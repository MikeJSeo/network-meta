#' Make a network object for the unrelated mean effects model (inconsistency model) containing data, priors, and a JAGS model file
#' 
#' This is similar to the function \code{\link{network.data}}, except this is used for the unrelated mean effects model.
#'
#' @param Outcomes Arm-level outcomes. If it is a multinomial response, the matrix would be arms (row) by multinomial categories (column). If it is binomial or normal, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm. Treatments should have positive integer values starting from 1 to total number of treatments.
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses.
#' @param SE A vector of standard error for each arm. Used only for normal response.
#' @param response Specification of the outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial".
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param rank.preference Set it equal to "higher" if higher values are preferred (i.e. assumes events are good). Set it equal to "lower" if lower values are preferred (i.e. assumes events are bad). Default is "higher".
#' @param mean.d Prior mean for the relative effect
#' @param prec.d Prior precision for the relative effect
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal. It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter.
#' @param dic This is an indicator for whether user wants to calculate DIC. Model stores less information if you set it to FALSE.
#' @references S. Dias, N.J. Welton, A.J. Sutton, D.M. Caldwell, G. Lu, and A.E. Ades (2013), \emph{Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials}, Medical Decision Making 33(5):641-656. [\url{https://doi.org/10.1177/0272989X12455847}]
#' @export

ume.network.data <- function(Outcomes, Study, Treat, N = NULL, SE = NULL, response = NULL, type = "random",
                             mean.d = 0, prec.d = 0.0001, hy.prior = list("dunif", 0, 5), dic = TRUE){
  
  if(missing(Study) || missing(Treat) || missing(Outcomes)){
    stop("Study, Treat, and Outcomes have to be all specified")
  }
  
  if(response == "multinomial" || response == "binomial"){
    if(is.null(N)){
      stop("If the response is multinomial or binomial, N has to be specified")
    }
  } else if (response == "normal"){
    if(is.null(SE)){
      stop("If the response is normal, SE has to be specified")
    }
  }
  
  if(!type %in% c("fixed", "random")){
    stop("type has to be either fixed or random")
  }
  
  na <- rle(Study)$lengths
  if(any(na == 1)) stop("study cannot have only 1 arm or arms have to be next to each other in each study")
  
  nstudy <- length(unique(Study))
  ntreat <- unique(as.vector(Treat))
  ntreat <- length(ntreat[!is.na(ntreat)])
  Outcomes <- as.matrix(Outcomes)
  
  ends <- cumsum(na) # End row of trials
  starts <- c(1, ends[-length(ends)] + 1) # Start row of trials
  b.Treat <- rep(NA, length(na))
  b.id <- rep(F, sum(na))
  for (i in 1:length(na)){
    limits <- starts[i]:ends[i]
    b.Treat[i] <- min(Treat[limits])
    b.id[limits[b.Treat[i] == Treat[limits]]] <- T
  }
  
  # make r, t, se(or n) that has dimensions suitable for rjags coding.
  t <- make.byStudy.matrix(Treat, Study)
  if(response == "binomial" || response == "multinomial"){
    n <- make.byStudy.matrix(N, Study)
  } else if(response == "normal"){
    se <- make.byStudy.matrix(SE, Study)
  }
  r <- make.byStudy.Outcome(Outcomes, Study, nstudy, na)
  if(response != "multinomial"){
    r <- r[,,1]
  }
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, response = response, hy.prior = hy.prior, mean.d = mean.d, prec.d = prec.d, dic = dic)
  
  if(response == "binomial"){
    network$N = N
    network$n = n
  } else if (response == "normal"){
    netowkr$SE = SE
    network$se = se
  }
  
  code <- ume.network.rjags(network)
  network$code <- code
  
  class(network) <- "ume.network.data"
  return(network)
}

ume.network.rjags <- function(network){
  
  network <- with(network, {
    if(response == "binomial"){
      ume.binomial.rjags(network)  
    }
  })
  
}


ume.binomial.rjags <- function(network){
  
  with(network, {
    
    code <- paste0("model\n{",
                   "\n\tfor(i in 1:", nstudy, ") {",
                   "\n\t\tdelta[i,1] <- 0",
                   "\n\t\tmu[i] ~ dnorm(0,.0001)",
                   "\n\t\tfor(k in 1:na[i]) {",
                   "\n\t\t\tr[i,k] ~ dbin(p[i,k], n[i,k])")
    
    
    if(type == "fixed"){
      code <- paste0(code, "\n\t\t\tlogit(p[i,k]) <- mu[i] + d[t[i,1], t[i,k]]")
    } else if(type == "random"){
      code <- paste0(code, "\n\t\t\tlogit(p[i,k]) <- mu[i] + delta[i,k]")
    }
                   
    code <- paste0(code,          
                   "\n\t\t\trhat[i,k] <- p[i,k] * n[i,k]",
                   "\n\t\t\tdev[i,k] <- 2 * (r[i,k] * (log(r[i,k])- log(rhat[i,k])) + (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k])))",
                   "\n\t\t}",
                   "\n\t\tresdev[i] <- sum(dev[i,1:na[i]])"
                   )
    
    if(type == "random"){
      code <- paste0(code, "\n\t\tfor (k in 2:na[i]) {",
                           "\n\t\t\tdelta[i,k] ~ dnorm(d[t[i,1],t[i,k]], tau)",
                           "\n\t\t}")
    }
    
    code <- paste0(code, 
                   "\n\t}",
                   "\n\ttotresdev <- sum(resdev[])")
    
    if(type == "fixed"){
      code <- paste0(code, "\n\tfor(k in 1:", ntreat, ") {",
                     "\n\t\td[k,k] <- 0",
                     "\n\t}")
    }
                   
    code <- paste0(code,
                   "\n\tfor(c in 1:", ntreat -1, ") {",
                   "\n\t\tfor(k in (c+1):", ntreat, ") {",
                   "\n\t\t\td[c,k] ~ dnorm(", mean.d, ", ", prec.d, ")",
                   "\n\t\t}",
                   "\n\t}")
                   
    if(type == "random"){
      code <- paste0(code, ume.hy.prior.rjags(hy.prior))
    }
    
    code <- paste0(code, "\n}")
    return(code)
  })
}



ume.hy.prior.rjags <- function(hy.prior){
  
  code <- ""
  distr <- hy.prior[[1]]
  if (distr == "dunif") {
    code <- paste0(code,
                   "\n\tsd ~ dunif(hy.prior.1, hy.prior.2)",
                   "\n\ttau <- pow(sd,-2)")
  } else if(distr == "dgamma"){
    code <- paste0(code,
                   "\n\tsd <- pow(tau, -0.5)",
                   "\n\ttau ~ dgamma(hy.prior.1, hy.prior.2)")
  } else if(distr == "dhnorm"){
    code <- paste0(code,
                   "\n\tsd ~ dnorm(hy.prior.1, hy.prior.2)T(0,)",
                   "\n\ttau <- pow(sd, -2)")
  }
  return(code)
}


#' Run the model using the network object
#' 
#' This is similar to the function \code{\link{network.run}}, except this is used for the unrelated mean effects model.
#'
#' @param network network object created from \code{\link{ume.network.data}} function
#' @param inits Initial values for the parameters being sampled. If left unspecified, program will generate reasonable initial values.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that user is willing to run. If the algorithm is not converging, it will run up to \code{max.run} iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks. If the algorithm converges fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run Final number of iterations that the user wants to store. If after the algorithm converges, user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run extra iterations to reach the specified number of runs
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence diagnostic. Point estimate is used to test convergence of parameters for study effect (eta), relative effect (d), and heterogeneity (log variance (logvar)).
#' @param extra.pars.save Parameters that user wants to save besides the default parameters saved. See code using \code{cat(network$code)} to see which parameters can be saved.
#' @return
#' \item{data_rjags}{Data that is put into rjags function jags.model}
#' \item{inits}{Initial values that are either specified by the user or generated as a default}
#' \item{pars.save}{Parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' \item{deviance}{Contains deviance statistics such as pD (effective number of parameters) and DIC (Deviance Information Criterion)}
#' \item{rank.tx}{Rank probability calculated for each treatments. \code{rank.preference} parameter in \code{\link{ume.network.data}} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' network <- with(thrombolytic, {
#'  ume.network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- ume.network.run(network)
#' @export

ume.network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                                 conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "ume.network.data")) {
    stop('Given network is not ume.network.data. Run ume.network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- list(r = r, t = t, na = na)
    
    if(response == "binomial"){
      data$n <- n
    }
    
    if(type == "random"){
      data$hy.prior.1 <- hy.prior[[2]]
      data$hy.prior.2 <- hy.prior[[3]]
    }
    
    pars.save <- c("d")
    
    if(type == "random"){
      pars.save <- c(pars.save, "sd", "delta")  
    }
    
    if(dic == TRUE){
      pars.save <- c(pars.save, "totresdev", "resdev")
      if(response == "binomial"){
        pars.save <- c(pars.save, "rhat", "dev")
      }
    }
    
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save, pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    if(is.null(inits)){
      inits <- ume.network.inits(network, n.chains)
    }
    samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
    result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
    result <- c(result, samples)
    
    result$deviance <- calculate.deviance(result)
    
#    result$rank.tx <- rank.tx(result)
    class(result) <- "ume.network.result"
    return(result)
  })
}


ume.network.inits <- function(network, n.chains){
  
  response <- network$response
  
  inits <- if(response == "multinomial"){
    ume.multinomial.inits(network, n.chains)
  } else if(response == "binomial"){
    ume.binomial.inits(network, n.chains)
  } else if(response == "normal"){
    ume.normal.inits(network, n.chains)
  }
  return(inits)
}

ume.binomial.inits <- function(network, n.chains){
  
  with(network,{
    
    Outcomes <- Outcomes + 0.5 # ensure ratios are always defined
    N <- N + 1
    p <- Outcomes/N
    logits <- log(p/(1-p))
    se.logits <- sqrt(1/Outcomes + 1/(N - Outcomes))
    
    mu <- logits[b.id]
    se.mu <- se.logits[b.id]
    delta <- logits - rep(mu, times = na)
    delta <- delta[!b.id,]
    
    inits = make.inits(network, n.chains, delta, mu, se.mu)
    return(inits)  
  })
  
}



make.inits <- function(network, n.chains, delta, mu, se.mu){
  
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
    d <- se.d <- rep(NA, ntreat)
    d[-1] <- coef(fit)[,1]
    se.d[-1] <- coef(fit)[,2]
    resid.var <- fit$sigma^2
    
    
    ############# Generate initial values
    initial.values = list()
    for(i in 1:n.chains){
      initial.values[[i]] = list()
    }
    for(i in 1:n.chains){
      random.mu <- rnorm(length(mu))
      initial.values[[i]][["mu"]] <- mu + se.mu * random.mu
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
            if(sqrt(sigma2) > hy.prior[[3]]){
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
    return(initial.values)
  })
}


pick.summary.variables.ume <- function(result, extra.pars = NULL, only.pars = NULL){
  samples <- result[["samples"]]
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  varnames.split <- gsub("[[:digit:]]","",varnames.split)
  
  if(!is.null(only.pars)){
    if(!all(only.pars %in% varnames.split)){
      stop(paste0(only.pars, "was not sampled"))
    }
  }
  if(is.null(only.pars)){
    pars <- c("d", "sd", "sigma")
  } else{
    pars <- only.pars
  }
  if(!is.null(extra.pars)){
    if(!extra.pars %in% varnames.split){
      stop(paste0(extra.pars, " is not saved in result"))
    }
    pars <- c(pars, extra.pars)
  }
  summary.samples <- lapply(samples, function(x){x[,varnames.split %in% pars, drop = F]})
  summary.samples <- coda::mcmc.list(summary.samples)
  summary.samples
}



#' Summarize result run by \code{\link{ume.network.run}}
#'
#' This function uses summary function in coda package to summarize mcmc.list object. Monte carlo error (Time-series SE) is also obtained using the coda package and is printed in the summary as a default.
#'
#' @param object Result object created by \code{\link{ume.network.run}} function
#' @param ... Additional arguments affecting the summary produced
#' @examples
#' network <- with(smoking, {
#'  ume.network.data(Outcomes, Study, Treat, N = N, response = "binomial", type = "random")
#' })
#' result <- ume.network.run(network) 
#' summary(result)
#' #summary(result, only.pars = "sd")
#' #summary(result, extra.pars = c("delta"))
#' @export

summary.ume.network.result <- function(object, ...){
  
  if(!inherits(object, "ume.network.result")) {
    stop('This is not the output from ume.network.run. Need to run ume.network.run function first')
  }
  summary.samples <- pick.summary.variables.ume(object, ...)
  
  rval <- list("summary.samples"= summary(summary.samples),
               "deviance" = unlist(object$deviance[1:3]),
               "total_n" = sum(object$network$na))
  class(rval) <- 'summary.ume.network.result'
  rval
}




#' Plot traceplot and posterior density of the result using contrast data
#'
#' This function uses plotting function in coda package to plot mcmc.list object
#'
#' @param x Result object created by \code{\link{ume.network.run}} function
#' @examples
#' network <- with(smoking, {
#'  ume.network.data(Outcomes, Study, Treat, N = N, response = "binomial", type = "random")
#' })
#' result <- ume.network.run(network, only.pars = "sd")
#' plot(result)
#' @export

plot.ume.network.result <- function(x, ...) {
  
  if(!inherits(x, "ume.network.result")) {
    stop('This is not the output from ume.network.run. Need to run ume.network.run function first')
  }
  summary.samples <- pick.summary.variables.ume(x, ...)
  plot(summary.samples)
}


