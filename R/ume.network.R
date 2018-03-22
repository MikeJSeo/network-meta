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


r[,1] n[,1] r[,2] n[,2] r[,3] n[,3] t[,1] t[,2] t[,3] na[]
9 140 23 140 10 138 1 3 4 3 # trial 1 ACD
11 78 12 85 29 170 2 3 4 3 # trial 2 BCD
75 731 363 714 NA 1 1 3 NA 2 # 3
2 106 9 205 NA 1 1 3 NA 2 # 4
58 549 237 1561 NA 1 1 3 NA 2 # 5
0 33 9 48 NA 1 1 3 NA 2 # 6
3 100 31 98 NA 1 1 3 NA 2 # 7
1 31 26 95 NA 1 1 3 NA 2 # 8
6 39 17 77 NA 1 1 3 NA 2 # 9
79 702 77 694 NA 1 1 2 NA 2 # 10
18 671 21 535 NA 1 1 2 NA 2 # 11
64 642 107 761 NA 1 1 3 NA 2 # 12
5 62 8 90 NA 1 1 3 NA 2 # 13
20 234 34 237 NA 1 1 3 NA 2 # 14
0 20 9 20 NA 1 1 4 NA 2 # 15
8 116 19 149 NA 1 1 2 NA 2 # 16
95 1107 143 1031 NA 1 1 3 NA 2 # 17
15 187 36 504 NA 1 1 3 NA 2 # 18
78 584 73 675 NA 1 1 3 NA 2 # 19
69 1177 54 888 NA 1 1 3 NA 2 # 20
20 49 16 43 NA 1 2 3 NA 2 # 21
7 66 32 127 NA 1 2 4 NA 2 # 22
12 76 20 74 NA 1 3 4 NA 2 # 23
9 55 3 26 NA 1 3 4 NA 2 # 24
  
  