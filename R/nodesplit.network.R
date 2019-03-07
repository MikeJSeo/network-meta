#' Assessing consistency using node-splitting model
#'
#' This function is used to make node-splitting inconsistency model. The structure is similar to the function \code{\link{network.data}}
#'
#' @param Outcomes Arm-level outcomes. If it is a multinomial response, the matrix would be arms (row) by multinomial categories (column). If it is binomial or normal, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm. Treatments should have positive integer values starting from 1 to total number of treatments. In a study, lowest number is taken as the baseline treatment.
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses.
#' @param SE A vector of standard error for each arm. Used only for normal response.
#' @param response Specification of the outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial".
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @references S. Dias, N.J. Welton, D.M. Caldwellb, A.E. Ades (2010), \emph{Checking consistency in mixed treatment}, Statistics in Medicine 29(7-8, Sp. Iss. SI): 932-944. [\url{https://doi.org/10.1002/sim.3767}]
#' @export

nodesplit.network.data <- function(Outcomes, Study, Treat,  N = NULL, SE = NULL, response = NULL, type = "random"){
  
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
  
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, b.Treat = b.Treat response = response)
  
  if(response == "binomial" || response == "multinomial"){
    network$N = N
    network$n = n
  } else if (response == "normal"){
    network$SE = SE
    network$se = se
  }
  
  #code <- nodesplit.network.rjags(network)
  #network$code <- code
  
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
    
    code <- paste0("model\n",
                   "\n\tfor(i in 1:", nstudy, ") {",
                   "w[i,1] <- 0",
                   "j[i,1] <- 0",
                   "delta[i,bi[i]] <- 0"
                   
                   
                   # w[i,1] <-0
                   # j[i,1] <-0
                   # delta[i,bi[i]] <- 0
                   # mu[i] ~ dnorm(0,.0001)                                    # vague priors for 24 trial baselines
                   # for (k in 1:na[i])  {
                   
                   
                   )
    
    code <- paste0("model\n{",
                   "\n\tfor(i in 1:", nstudy, ") {",
                   "\n\t\tdelta[i,1] <- 0",
                   "\n\t\tmu[i] ~ dnorm(mean.mu,prec.mu)",
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
                     "\n\t\t\tdelta[i,k] ~ dnorm(d[t[i,1],t[i,k]], prec)",
                     "\n\t\t}")
    }
    
    code <- paste0(code, 
                   "\n\t}",
                   "\n\ttotresdev <- sum(resdev[])")
    
    code <- paste0(code,
                   "\n\tfor(c in 1:", ntreat -1, ") {",
                   "\n\t\tfor(k in (c+1):", ntreat, ") {",
                   "\n\t\t\td[c,k] ~ dnorm(mean.d, prec.d)",
                   "\n\t\t}",
                   "\n\t}")
    
    if(type == "random"){
      code <- paste0(code, ume.hy.prior.rjags(hy.prior, 0))
    }
    
    code <- paste0(code, "\n}")
    return(code)
  })
}