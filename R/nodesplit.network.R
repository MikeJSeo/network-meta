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
  
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, response = response)
  
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