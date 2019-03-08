#' Assessing consistency using node-splitting model
#'
#' This function is used to make node-splitting inconsistency model. The structure is similar to the function \code{\link{network.data}}
#'
#' @param Outcomes Arm-level outcomes. If it is a multinomial response, the matrix would be arms (row) by multinomial categories (column). If it is binomial or normal, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm. Treatments should have positive integer values starting from 1 to total number of treatments. In a study, lowest number is taken as the baseline treatment. Also, in each study first arm should be the baseline treatment and lower number treatment arms should come before the higher ones. (i.e. treatment arm should go 1, 2, 4 and not 1, 4, 2). See example data for clarification.
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses.
#' @param SE A vector of standard error for each arm. Used only for normal response.
#' @param response Specification of the outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial".
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @references S. Dias, N.J. Welton, D.M. Caldwellb, A.E. Ades (2010), \emph{Checking consistency in mixed treatment}, Statistics in Medicine 29(7-8, Sp. Iss. SI): 932-944. [\url{https://doi.org/10.1002/sim.3767}]
#' @export

nodesplit.network.data <- function(Outcomes, Study, Treat,  N = NULL, SE = NULL, response = NULL, type = "random"){
  
  orig <- network.data(Outcomes, Study, Treat, N = N, response = response, type = type)
  
  network <- with(orig, {  list(data = data, Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, response = response)})
    
  #  list(Outcomes = orig$Outcomes, Study = orig$Study, Treat = orig$Treat, r = orig$r, t = orig$t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, response = response)
  
  if(response == "binomial" || response == "multinomial"){
    network$N = N
    network$n = n
  } else if (response == "normal"){
    network$SE = SE
    network$se = se
  }
  
  code <- nodesplit.network.rjags(network)
  network$code <- code
  
  class(network) <- "nodesplit.network.data"
  return(network)
}

nodesplit.network.rjags <- function(network){
  
  network <- with(network, {
    if(response == "binomial"){
      nodesplit.binomial.rjags(network)  
    } else if(response == "normal"){
      nodesplit.normal.rjags(network)
    } else if(response == "multinomial"){
      nodesplit.multinomial.rjags(network)
    }
  })
}


nodesplit.binomial.rjags <- function(network){
  
  with(network, {
    
    code <- paste0("model{",
                   "\n\tfor(i in 1:", nstudy, ") {",
                   "\n\t\tw[i,1] <- 0",
                   "\n\t\tj[i,1] <- 0",
                   "\n\t\tdelta[i,1] <- 0",
                   "\n\t\tmu[i] ~ dnorm(0,.0001)",
                   "\n\t\tfor (k in 1:na[i]) {",
                   "\n\t\t\tr[i,k] ~ dbin(p[i,t[i,k]], n[i,k])",
                   "\n\t\t\tlogit(p[i,t[i,k]]) <- mu[i] + delta[i,t[i,k]]",
                   "\n\t\t\tindex[i,k] <- split[i] * (equals(t[i,k], pair[1]) + equals(t[i,k], pair[2]))",
                   "\n\t\t\trhat[i,k] <- p[i,t[i,k]] * n[i,k]",
                   "\n\t\t\tdev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k]))
                          + (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k])))",
                   "\n\t\t}",
                   "\n\t\tresdev[i] <- sum(dev[i,1:na[i]])",
                   "\n\t\tfor(k in 2:na[i]){",
                   "\n\t\t\tdelta[i,t[i,k]] ~ dnorm(md[i,t[i,k]], taud[i,t[i,k]])",
                   "\n\t\t\t# mean of LOR distritbutions, split into direct and indirect",
                   "\n\t\t\tmd[i,t[i,k]] <- (d[t[i,k]] - d[t[i,1]]+sw[i,k]) * (1- index[i,k]) + direct * index[i,k]",
                   "\n\t\t\t# adjustment for multi-arm RCTs with correction for arms removed to split node",
                   "\n\t\t\tj[i,k] <- k - (equals(1, split[i]) * step(k-3))",
                   "\n\t\t\ttaud[i, t[i,k]] <- tau * 2 * (j[i,k] -1)/ j[i,k]",
                   "\n\t\t\tw[i,k] <- (delta[i, t[i,k]] - d[t[i,k]] + d[t[i,1]]) * (1- index[i,k])",
                   "\n\t\t\tsw[i,k] <- sum(w[i,1:(k-1)])/ (j[i,k] - 1)",
                   "\n\t\t}",
                   "\n\t}",
                   "\n\td[1] <- 0",
                   "\n\tdirect ~ dnorm(0, .0001)",
                   "\n\tfor (k in 2:nt) { d[k] ~ dnorm(0, .0001)}",
                   "\n\tsd ~ dunif(0,10)",
                   "\n\tvar <- pow(sd, 2)",
                   "\n\ttau <- 1/var",
                   "\n\ttotresdev <- sum(resdev[])",
                   "\n\tfor (c in 1:", ntreat-1, ") {",
                   "\n\t\tfor(k in (c+1):nt) {",
                   "\n\t\t\tor[c,k] <- exp(d[k] - d[c])",
                   "\n\t\t\tlor[c,k] <- d[k] - d[c]",
                   "\n\t\t}",
                   "\n\t}",
                   "\n\tdiff <- direct - lor[pair[1], pair[2]]",
                   "\n\t#calculate p-value",
                   "\n\tprob <- step(diff)")

      
    code <- paste0(code, "\n}")
    return(code)
  })
}



#' Run the model using the network object
#' 
#' This is similar to the function \code{\link{network.run}}, except this is used to test inconsistency.
#'
#' @param network network object created from \code{\link{nodesplit.network.data}} function
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
#' \item{rank.tx}{Rank probability calculated for each treatments. \code{rank.preference} parameter in \code{\link{nodesplit.network.data}} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' network <- with(smoking, {
#'  nodesplit.network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- nodesplit.network.run(network)
#' @export

nodesplit.network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                            conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "nodesplit.network.data")) {
    stop('Given network is not nodesplit.network.data. Run nodesplit.network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- list(r = r, t = t, na = na)
    
    if(response == "binomial" || response == "multinomial"){
      data$n <- n
    } else if(response == "normal"){
      data$se <- se
    }
    
    # if(type == "random"){
    #   data$hy.prior.1 <- hy.prior[[2]]
    #   data$hy.prior.2 <- hy.prior[[3]]
    # }
    
    data$mean.d = mean.d
    data$prec.d = prec.d
    data$mean.mu = mean.mu
    data$prec.mu = prec.mu
    
    pars.save <- c("d")
    
    if(type == "random"){
      pars.save <- c(pars.save, "delta")
      if(response %in% c("normal", "binoimal")){
        pars.save <- c(pars.save, "sd")
      } else if (response == "multinomial"){
        pars.save <- c(pars.save, "sigma")
      }
    }
    
    # if(dic == TRUE){
    #   pars.save <- c(pars.save, "totresdev", "resdev", "dev")
    #   if(response == "binomial" || response == "multinomial"){
    #     pars.save <- c(pars.save, "rhat")
    #   } else if(response == "normal"){
    #     pars.save <- c(pars.save, "theta")
    #   }
    # }
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save, pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    # if(is.null(inits)){
    #   inits <- ume.network.inits(network, n.chains)
    # }
    samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
    result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
    result <- c(result, samples)
    
    result$deviance <- calculate.deviance(result)
    
    class(result) <- "nodesplit.network.result"
    return(result)
  })
}