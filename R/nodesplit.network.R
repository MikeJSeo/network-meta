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
  
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, type = type, rank.preference = NULL, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, response = response, hy.prior = hy.prior)
  
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
      ume.binomial.rjags(network)  
    } else if(response == "normal"){
      ume.normal.rjags(network)
    } else if(response == "multinomial"){
      ume.multinomial.rjags(network)
    }
  })
}