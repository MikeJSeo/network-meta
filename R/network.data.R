#' Make a network object containing data, priors, and a jags model file
#'
#' @param Outcomes Arm-level outcomes. If it is multinomial response, the matrix would have arms for the row and categories for the column. Otherwise, it would be a vector.
#' @param Study A vector of study indicator for each arm
#' @param Treat A vector of treatment indicator for each arm
#' @param N A vector of total number of observations in each arm. Used for binomial and multinomial responses
#' @param SE A vector of standard error for each arm. Used only for normal response.
#' @param response Specification for Outcomes type. Must specify one of the following: "normal", "binomial", or "multinomial".
#' @param Treat.order This specifies the treatment order and therefore how treatments are compared. The first treatment that is specified is considered as a base treatment. Default order would be alphabetical. This would hold true for numbers. If the treatments are coded 1,2,etc, then treatment with a value of 1 would be assigned as baseline treatment.
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. The default is "random".
#' @param rank.preference Set it equal to "higher" if higher values are preferred (i.e. assumes events are good). Set it equal to "lower if lower values are preferred (i.e. assumes events are bad).
#' @param miss.matrix This is parameter for running incomplete multinomial. Still under revision.
#' @param baseline "independent" allows baseline effects for each treatment while "common" restricts same baseline effect for all treatment. Default is "none" which doesn't incorporate baseline effect.
#' @param covariate A covariate matrix with each row representing each trial and column representing each covariate. Covariate information is needed for each study, and so the user doesn't need to repeatly specify covariates for each arm.
#' @param covariate.type Should be a vector indicating the type of the covariate. Covariate can be either "continuous" or "discrete". If it continuous, covariates are centered. If the covariate is discrete it is not centered and it has to be in a dummy integer format (i.e. 0,1,2,...). The code doesn't factor the covariates for the user, so user needs to specify dummy variables if factor is needed.
#' @param covariate.model "independent" allows covariate effects for each treatment. "common" restricts same covariate effect for all treatment. Lastly, "exchangeable"  assumes that the covariate effects are different but related and strength is borrowed across them. We set "common" to be default.
#' @param dic This is an indicator for whether you want to calculate DIC. It stores less information if you set it to FALSE.
#' @param mean.d Prior mean for the relative effect
#' @param prec.d Prior precision for the relative effect
#' @param mean.Eta Prior mean for the study effect
#' @param prec.Eta Prior precision for the study effect
#' @param mean.bl Prior mean for the baseline effect
#' @param prec.bl Prior precision for the baseline effect
#' @param hy.prior.bl Between treatment heterogeneity in baseline effect (for exchangeable regression coefficient only). Format of the data is same as hy.prior.
#' @param mean.cov Prior mean for the covariate effect
#' @param prec.cov Prior precision for the covariate effect
#' @param hy.prior.cov Between treatment heterogeneity in covariate effect (for exchangeable regression coefficient only). Format of the data is same as hy.prior. Default is set to be dunif(0, 5) for binary, dunif(0, 100) for normal, and wishart with identity scale matrix and (# of categories - 1) degrees of freedom.
#' @param hy.prior Prior for the heterogeneity parameter. Supports uniform, gamma, and half normal for normal and binomial response and wishart for multinomial response. It should be a list of length 3, where first element should be the distribution (one of dunif, dgamma, dhnorm, dwish) and the next two are the parameters associated with the distribution. For example, list("dunif", 0, 5) give uniform prior with lower bound 0 and upper bound 5 for the heterogeneity parameter. For wishart distribution, the last two parameter would be the scale matrix and the degrees of freedom.
#' @return Creates list of variables that are used to run the model.
#' \item{data}{Data combining all the input data. User can check this to insure the data is correctly specified. For modelling purposes, any character valued study or treatment variable is changed to numeric values based on alphabetical order}
#' \item{nrow}{Total number of arms in the meta-analysis}
#' \item{ncol}{Number of columns in the Outcomes. Will equal 1 for binary and normal and number of categories for multinomial}
#' \item{nstudy}{Number of study}
#' \item{na}{Number of arms for each study}
#' \item{ntreat}{Number of treatment}
#' \item{b.id}{Indicator in sequence of all treatments for which treatment is base treatment in Study}
#' \item{t}{Treat made into a matrix which has dimensions number of study by max number of arms in studies}
#' \item{r}{Outcomes made into array that is suitable for use in our rjags code. For multinomial, it has 3 dimensions: number of study by max number of arms in studies by number of categories.}
#' \item{mx}{If the continuous covariate is included, it calculates the mean of the covariates which is used to center the covariates. The numeric indicator after mx refers to column number of the covariates if there is more than one covariates included. For discrete covariate, covariates are not centered.}
#' \item{mx_bl}{If the baseline effect is specified, it also calculates the mean baseline risk}
#' \item{prior.data}{Prior data created using the user inputs or default values. If no user input is specifies for the prior, it uses default values}
#' \item{code}{Rjags model file code that is generated using information provided by the user. To view model file inside R, use \code{cat(network$code)}}
#' @examples
#' ###Blocker data example
#' blocker
#' Outcomes <- blocker[["Outcomes"]]
#' Study <- blocker[["Study"]]
#' Treat <- blocker[["Treat"]]
#' N <- blocker[["N"]]
#' network <- network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' str(network)
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

network.data <- function(Outcomes, Study, Treat, N = NULL, SE = NULL, response = "multinomial", Treat.order = NULL, type = "random", rank.preference = "higher",
                         miss.matrix = NULL, baseline = "none", covariate = NULL, covariate.type = NULL, covariate.model = NULL, dic = TRUE,
                         mean.d = NULL, prec.d = NULL, mean.Eta = NULL, prec.Eta = NULL, mean.bl = NULL, prec.bl = NULL, hy.prior.base = NULL,
                         mean.cov = NULL, prec.cov = NULL, hy.prior.cov = NULL, hy.prior = NULL) {

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

  if(is.null(baseline)){
    stop("baseline has to be none, independent, exchangeable, or common")
  } else if(!baseline %in% c("none", "independent", "exchangeable", "common")){
    stop("baseline has to be none, independent, exchangeable, or common")
  }

  if(!type %in% c("fixed", "random")){
    stop("type has to be either fixed or random")
  }

  if(!is.null(covariate)){
    if(!is.matrix(covariate) && !is.vector(covariate)){
      stop("covariate has to be vector if there is only one and matrix if there is more than one")
    }
    if(is.null(covariate.model)){
      covariate.model <- "common"
    }
    covariate <- as.matrix(covariate)
    if(is.null(covariate.type)){
      covariate.type = rep("continuous", dim(covariate)[2])
    }
    if(!all(apply(covariate, 2, is.numeric))){
      stop("covariate has to be numeric type")
    }
  }

  if(!is.null(Treat.order)){
    check.Treat.order(Treat, Treat.order)
  } else{
    Treat.order <- sort(unique(Treat))
  }

  Orig_Study <- Study
  Orig_Treat <- Treat

  #relabel Treatment and Study names into to a numeric sequence (1:nstudy)
  na <- rle(Study)$lengths
  print(na)
  print(Study)
  if(any(na == 1)) stop("study cannot have only 1 arm or arms have to be next to each other in each study")
  Study <- rep(1:length(unique(Study)), times = na)
  Treat <- relabel.vec(Treat, Treat.order)

  data <- if(response == "multinomial" || response == "binomial"){
    cbind(Outcomes, N, Study, Orig_Study, Treat, Orig_Treat)
  } else if(response == "normal"){
    cbind(Outcomes, SE, Study, Orig_Study, Treat, Orig_Treat)
  }
  data <- as.data.frame(data)

  nrow <- dim(data)[1]
  Outcomes <- as.matrix(Outcomes)
  ncol <- dim(Outcomes)[2]
  nstudy <- length(unique(Study))
  ntreat <- length(unique(Treat))

  # permute the data so that base treatment arm is always listed first in each study
  data <- data[order(as.numeric(as.character(data$Study)), data$Treat),,drop=FALSE]
  data[!names(data) %in% c("Orig_Study", "Orig_Treat")] <- suppressWarnings(sapply(data[!names(data) %in% c("Orig_Study", "Orig_Treat")], function(x) {as.numeric(paste(x))}))
  rownames(data) <- seq(nrow(data))
  Outcomes <- as.matrix(data[,1:ncol])
  Treat <- data$Treat
  Study <- data$Study

  if(response == "binomial" || response == "multinomial"){
    N <- data$N
  } else if(response == "normal"){
    SE <- data$SE
  }

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
  network <- list(data = data, Outcomes = Outcomes, Study = Study, Treat = Treat, Treat.order = Treat.order, type = type, rank.preference = rank.preference, miss.matrix = miss.matrix, nrow = nrow, ncol = ncol, nstudy = nstudy, na = na, ntreat = ntreat, b.id = b.id, t = t, r = r, response = response, baseline = baseline, covariate = covariate, covariate.model = covariate.model, dic = dic)

  ######################### missing pattern ##############################
  # check miss.matrix is specified correctly (multinomial only)
  ncat <- ifelse(is.null(miss.matrix), ncol, dim(miss.matrix)[2])
  if(!is.null(miss.matrix)){
    check.miss.matrix(miss.matrix, Outcomes, nrow, ncol, ncat)
  }
  if(response == "multinomial"){
    pattern <- find.pattern(Outcomes, nrow, ncol)
    npattern <- length(unique(pattern))
  }

  if(response == "multinomial"){
    # calculate a list called miss.patterns which is used in later function calls
    miss.patterns <- if(!is.null(miss.matrix)){
      create.miss.patterns(pattern, miss.matrix, Outcomes, nrow, ncat, ncol, Study)
    } else{
      list(list(list(Study,seq(ncat))),diag(ncat))
    }
    network$miss.patterns = miss.patterns
    network$ncat = ncat
    network$N = N

    if(!is.null(miss.matrix)){
      network$npattern = npattern
      network$pattern = pattern
      rr <- divide.r(r, miss.patterns, npattern)
      nn <- divide.n(n, miss.patterns, npattern)
      network <- c(network, rr, nn)
      network[["r"]] <- NULL
    } else{
      network$n <- n
    }
  } else if(response == "binomial"){
    network$N = N
    network$n = n
  } else if (response == "normal"){
    network$SE = SE
    network$se = se
  }

  ###### Calculate baseline log odds
  if(baseline %in% c("independent", "common", "exchangeable")){
    if(response == "normal"){
      network$mx_bl = mean(r[,1][t[,1]==1], na.rm = TRUE)
    }
    if(response == "binomial"){
      rdummy = r[,1]
      rdummy[r[,1] == 0] = 0.5
      ndummy = n[,1]
      ndummy[r[,1] == 0] = ndummy[r[,1] == 0] + 1

      #take only the non-active control group (treatment A) to calculate the observed mean log odds
      p = (rdummy/ndummy)[t[,1]==1]
      lodds = log(p/(1-p))
      network$mx_bl = mean(lodds, na.rm = TRUE)
    }
    if(response == "multinomial"){
      #the first response in the multinomial is the reference, we will call it J
      category = which(apply((miss.patterns[[2]]), 1, sum) == 1)
      J = category[1]

      #take only the control group
      P_J = (r[,1,J]/n[,1])[t[,1]==1]
      P_j = matrix(0, ncol = ncat-1, nrow = sum(t[,1]==1))
      for(j in 2:ncat){
        P_j[,j-1] = (r[,1,category[j]]/n[,1])[t[,1]==1]
      }
      lodds = log(P_j/P_J)
      network$mx_bl = apply(lodds, 2, mean, na.rm = TRUE)
    }
  }
  # calculate mean of covariate
  if(!is.null(covariate)){
    for(i in 1:dim(covariate)[2]){
      nam <- paste("mx",i, sep = "")
      nam <- assign(nam, mean(covariate[,i], na.rm = TRUE))

      network[[paste("mx",i, sep = "")]] <- ifelse(covariate.type[i] == "continuous", nam, 0)
      nam2 <- paste("x", i, sep = "")
      nam2 <- assign(nam2, covariate[,i])

      network[[paste("x", i, sep = "")]] <- nam2
    }
  }

  # heterogeneity parameter
  if(!is.null(hy.prior)){
    check.hy.prior(hy.prior, response)
  } else{
    hy.prior <- hy.prior.default(network)
  }

  #between treatment heterogeneity in covariate effect
  if(!is.null(covariate)){
    if(covariate.model == "exchangeable"){
      if(!is.null(hy.prior.cov)){
        check.hy.prior(hy.prior.cov, response)
      } else{
        hy.prior.cov <- hy.prior.default(network)
      }
    }
  }

  if(baseline == "exchangeable"){
    if(!is.null(hy.prior.base)){
      check.hy.prior(hy.prior.base, response)
    } else{
      hy.prior.base <- hy.prior.default(network)
    }
  }

  prior.data <- network.prior.default(network, mean.d, prec.d, mean.Eta, prec.Eta, mean.bl, prec.bl, hy.prior.base, mean.cov, prec.cov, hy.prior.cov, hy.prior)
  if(type == "fixed"){
    prior.data <- prior.data[!names(prior.data) %in% c("hy.prior.1", "hy.prior.2")]
  }

  network$prior.data <- prior.data
  network$hy.prior <- hy.prior
  network$hy.prior.cov <- hy.prior.cov
  network$hy.prior.base <- hy.prior.base

  code <- network.rjags(network)
  network$code <- code

  class(network) <- "network.data"
  network
}

make.byStudy.Outcome = function(Outcomes, Study, nstudy, na){
  r = structure(.Data = rep(NA, nstudy*max(na)* dim(Outcomes)[2]), .Dim = c(nstudy, max(na), dim(Outcomes)[2]))

  arms_index = NULL
  for(i in 1:length(na)){
    arms_index = c(arms_index, seq(na[i]))
  }
  for(i in 1:dim(Outcomes)[1]){
    r[Study[i],arms_index[i],] = Outcomes[i,]
  }
  return(r)
}


make.byStudy.matrix = function(Treat, Study){

  #make the by-arms vector into a by-study matrix
  nstudy = length(unique(Study))
  na = rle(Study)$lengths
  Study = rep(1:nstudy, times = na)

  t = matrix(NA, nrow = nstudy, ncol = max(na))

  arms_index = NULL
  for(i in 1:length(na)){
    arms_index = c(arms_index, seq(na[i]))
  }

  for(i in 1:length(Treat)){
    t[Study[i], arms_index[i]] = Treat[i]
  }
  return(t)
}

check.Treat.order <- function(Treat, Treat.order){
  if(!is.null(Treat.order)){
    if(length(unique(Treat.order)) != length(unique(Treat))){
      stop("Need to specify treatment order for all treatments")
    }
    if(!all(Treat.order %in% Treat)){
      stop("Treat.order names have to be specified correctly.")
    }
  }
}

create.miss.patterns = function(pattern, miss.matrix, Outcomes, nrow, ncat, ncol, Study){
  #### create miss.patterns matrix #######
  num.pattern = length(levels(pattern))
  rows.all = vector("list", length = num.pattern)
  for(i in seq(num.pattern)){
    rows.all[[i]] = Study[pattern == levels(pattern)[i]]
  }

  missing = as.matrix(Outcomes)
  missing = replace(missing, !is.na(missing), 0)
  for(i in 1:ncol){
    missing[,i][missing[,i]==0] = i
  }

  miss.matrix2 = miss.matrix
  for(j in 1:dim(miss.matrix)[2]){
    miss.matrix2[,j][miss.matrix2[,j] == 1] = j
  }

  cols.data.pre = unique(missing)
  for(i in seq(num.pattern)){
    for(j in (ncat+1):ncol){
      if(any(!is.na(cols.data.pre[i,miss.matrix2[j,]]))){
        cols.data.pre[i,j] = NA
      }
    }
  }
  cols.data.all = vector("list", length = num.pattern)
  for(i in seq(num.pattern)){
    cols.data.all[[i]] = as.vector(cols.data.pre[i,][!is.na(cols.data.pre[i,])])
  }

  miss.patterns = vector("list")
  miss.patterns[[2]] = miss.matrix
  for(i in seq(num.pattern)){
    miss.patterns[[1]][[i]] = list(rows.all[[i]], cols.data.all[[i]])
  }
  return(miss.patterns)
}

check.miss.matrix <- function(miss.matrix, Outcomes, nrow, ncol, ncat){

  if(any(is.na(Outcomes))){
    if(is.null(miss.matrix)){
      stop("Data contains missing data. Need missing matrix specification")
    }
  }
  # Check if user entered miss.matrix correctly
  if(!is.null(miss.matrix)){
    if(dim(miss.matrix)[1] != ncol){
      stop("row length of miss.matrix should equal column length of Outcomes")
    }
    miss.matrix2 = miss.matrix
    for(i in 1:ncat){
      miss.matrix2[,i][miss.matrix2[,i] == 1] = i
    }
    for(i in 1:nrow){
      for(j in 1:ncol){
        if(!is.na(Outcomes[i,j]) && !is.na(sum(Outcomes[i, miss.matrix2[j,]]))){
          if(Outcomes[i,j] != sum(Outcomes[i, miss.matrix2[j,]])){
            print(Outcomes[i,j])
            print(Outcomes[i, miss.matrix2[j,]])
            stop("collapsed set of outcome do not add up; either your Outcomes input is wrong or you mispecified miss.matrix")
          }
        }
      }
    }
  }
}

find.pattern <- function(Outcomes, nrow, ncol){
  #Find missing pattern if any
  if(any(is.na(Outcomes))){
    missing = as.matrix(Outcomes)
    missing = replace(missing, !is.na(missing), 0)
    for(i in 1:ncol){
      missing[,i][missing[,i]==0] = i
    }
    pattern = vector("integer")
    for(i in seq(nrow)){
      pattern[i] = paste(missing[i,], collapse = "")
    }
    pattern = as.factor(pattern)
    levels(pattern) = 1:length(levels(pattern))
  } else{
    pattern = as.factor(rep(1, nrow))
    missing = NULL
  }
  return(pattern)
}

relabel.vec <- function(x, order)
{
  old.x <- x
  x <- rep(NA, length(old.x))
  for (i in seq(length(order))) x[old.x == order[i]] <- i #relabel studies in numerical order starting with one
  return(x)
}

divide.r <- function(r, miss.patterns, npattern){
  rr = vector("list", npattern)
  for(i in seq(npattern)){
    rows = unique(miss.patterns[[1]][[i]][[1]])
    rr[[i]] = r[rows,,miss.patterns[[1]][[i]][[2]]]
    names(rr)[i] = paste("r",i,sep = "")
  }
  return(rr)
}

divide.n <- function(n, miss.patterns, npattern){
  nn = vector("list", npattern)
  for(i in seq(npattern)){
    rows = unique(miss.patterns[[1]][[i]][[1]])
    nn[[i]] = n[rows,]
    names(nn)[i] = paste("n",i,sep="")
  }
  return(nn)
}

check.hy.prior <- function(hy.prior, response){

  if(length(hy.prior) != 3) stop("length of the hy.prior has to be 3 (Distribution, and two parameters for the distribution)")

  distr = hy.prior[[1]]
  stopifnot(distr %in% c("dunif", "dgamma", "dhnorm", "dwish"))

  if(response == "normal" || response == "binomial"){
    stopifnot(distr %in% c("dunif", "dgamma", "dhnorm"))
  } else if(response == "multinomial"){
    stopifnot(distr %in% c("dwish"))
  }
}
