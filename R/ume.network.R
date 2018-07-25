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
#' @references S. Dias, N.J. Welton, A.J. Sutton, D.M. Caldwell, G. Lu, and A.E. Ades (2013), \emph{Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials}, Medical Decision Making 33(5):641-656. [\url{https://doi.org/10.1177/0272989X12455847}]
#' @export

ume.network.data <- function(Outcomes, Study, Treat, N = NULL, SE = NULL, response = NULL, type = "random",
                             mean.d = NULL, prec.d = NULL, hy.prior = NULL){
  
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
  
  ntreat <- unique(as.vector(Treat))
  ntreat <- length(ntreat[!is.na(ntreat)])
  
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
  network <- list(Outcomes = Outcomes, Study = Study, Treat = Treat, r = r, t = t, n = n, se = se, type = type, rank.preference = rank.preference, miss.matrix = miss.matrix, nrow = nrow, ncol = ncol, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, t = t, r = r, response = response, baseline = baseline, baseline.risk = baseline.risk, covariate = covariate, covariate.model = covariate.model, dic = dic)
  
}



ume.network.rjags <- function(network){
  
  with(network, {
    
    code <- paste0("model\n{",
                   "\n\tfor(i in 1:", nstudy, ") {",
                   "\n\t\tdelta[i,1] <- 0",
                   "\n\t\tmu[i] ~ dnorm(0,.0001)",
                   "\n\t\tfor(k in 1:na[i]) {",
                   "\n\t\t\tr[i,k] ~ dbin(p[i,k], n[i,k])",
                   "\n\t\t\tlogit(p[i,k]) <- mu[i] + delta[i,k]",
                   "\n\t\t\trhat[i,k] <- p[i,k] * n[i,k]",
                   "\n\t\t\tdev[i,k] <- 2 * (r[i,k] * (log(r[i,k])- log(rhat[i,k])) + (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k])))",
                   "\n\t\t}",
                   "\n\t\tfor (k in 2:na[i]) {",
                   "\n\t\t\tdelta[i,k] ~ dnorm(md[i,k], taud[i,k])",
                   "\n\t\t\tmd[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k]",
                   "\n\t\t\ttaud[i,k] <- tau * 2 * (k-1) / k",
                   "\n\t\t\tw[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])",
                   "\n\t\t\tsw[i,k] <- sum(w[i,1:(k-1)])/ (k-1)",
                   "\n\t\t}",
                   "\n\t}",
                   "\n\ttotresdev <- sum(resdev[])",
                   "\n\td[1] <- 0",
                   "\n\tfor (k in 2:", ntreat, "){",
                   "\n\t\td[k] ~ dnorm(mean.d, prec.d)",
                   "\n\t\tsd ~ dunif(0,5)",
                   "\n\t\ttau <- pow(sd, -2)",
                   "\n\t}")
                   
    return(code)
  })
}
