#' Make a network object for contrast-level data containing data, priors, and a JAGS model file
#' 
#' This is similar to the function \code{\link{network.data}}, except it uses contrast-level data instead of arms-level data. Contrast-level format uses treatment differences relative to the control arm.
#' Note that in two arm trials there is only one contrast value per trial, but in three arm trials there are two contrast values relative to the control arm.
#'
#' @param Outcomes Contrast-level outcomes. Outcome is assumed to be normally distributed. Outcome should be a matrix with dimensions number of studies by maximum number of contrasts. If the maximum number of arms in a trial is three, then there should be two columns. See parkinsons_contrast data for an example. All the missing value in the matrix would be denoted as NA.
#' @param Treat A matrix of treatment for each arm. This will be a matrix with dimensions number of study by maximum number of arms. If the maximum arms in a trial is three, then the matrix should have three columns. All the missing value in the matrix should be denoted as NA. Treatments should have positive integer values starting from 1 to total number of treatments.
#' @param SE A matrix of standard error for each contrasts. The matrix would be same dimensions as Outcomes.
#' @param na A vector of number of arms in each study.
#' @param V Needed if you have multi-arm trials. Length of this vector should be number of studies. If the study is multi-arm trial, need to specify variance of the baseline treatment in that trial. Denote it with NA if the study only has two-arm trials.
#' @param type Type of model fitted: either "random" for random effects model or "fixed" for fixed effects model. Default is "random".
#' @param rank.preference Set it equal to "higher" if higher values are preferred (i.e. assumes events are good). Set it equal to "lower" if lower values are preferred (i.e. assumes events are bad). Default is "higher".
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}] 
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

contrast.network.data <- function(Outcomes, Treat, SE, na, V = NULL, type = "random", rank.preference = "higher"){

  if(missing(Outcomes) || missing(Treat) || missing(SE) || missing(na)){
    stop("Outcomes, Treat, SE, and na have to be all specified")
  }
  
  if(any(na == 1)) stop("study cannot have only 1 arm")
  if(is.unsorted(na)) stop("please sort the studies so that studies with higher number of arms are at the end (in a increasing order)")
  
  if(max(na) >= 3 & is.null(V)){
    stop("Need to specify variance of the baseline treatment if you have multi-arm trials")
  }
  
  if(!type %in% c("fixed", "random")){
    stop("type has to be either fixed or random")
  }
  
  if(!rank.preference %in% c("higher", "lower")){
    stop("rank preference has to be either higher or lower")
  }

  # Attach NA column for the first column
  Outcomes <- cbind(NA, Outcomes)
  SE <- cbind(NA, SE)
  
  na_count <- as.vector(table(na))
  ntreat <- unique(as.vector(parkinsons_contrast$Treat))
  ntreat <- length(ntreat[!is.na(ntreat)])
  
  network <- list(Outcomes = Outcomes, Treat = Treat, SE = SE, na = na, na_count = na_count, ntreat = ntreat)
  
  if(!is.null(V)){
    network$V <- V
  }
  
  code <- contrast.network.rjags(network)
  network$code <- code
  
  class(network) <- "contrast.network.data"
  return(network)

}

contrast.network.rjags <- function(network){
  
  with(network, {
    
    code <- paste0("model\n{",
                   "\n\tfor(i in 1:", na_count[1], ") {",
                   "\n\t\ty[i,2] ~ dnorm(delta[i,2], prec[i,2])",
                   "\n\t\tresdev[i] <- (y[i,2] - delta[i,2]) * (y[i,2] - delta[i,2]) * prec[i,2]",
                   "\n\t}")
    
    if(length(na_count) > 1){
    
      for(i in 2:length(na_count)){
        
        code <- paste0(code, "\n\tfor(i in ", cumsum(na_count)[i-1] + 1, ":", cumsum(na_count)[i], ") {", 
                       "\n\t\tfor(k in 1:(na[i]-1)) {",
                       "\n\t\t\tfor(j in 1:(na[i]-1)) {",
                       "\n\t\t\t\tSigma[i,j,k] <- V[i]*(1-equals(j,k)) + Var[i,k+1] * equals(j,k)",
                       "\n\t\t\t}",
                       "\n\t\t}",
                       "\n\t\tOmega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])",
                       "\n\t\ty[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]], Omega[i, 1:(na[i]-1), 1:(na[i]-1)])",
                       "\n\t\tfor(k in 1:(na[i]-1)){",
                       "\n\t\t\tydiff[i,k] <- y[i,(k+1)] - delta[i,(k+1)]",
                       "\n\t\t\tz[i,k] <- inprod(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])",
                       "\n\t\t}",
                       "\n\t\tresdev[i] <- inprod(ydiff[i,1:(na[i]-1)], z[i, 1:(na[i]-1)])",
                       "\n\t}")
      }  
    }
    
    code <- paste0(code, "\n\tfor(i in 1:", sum(na_count), ") {",
                   "\n\t\tw[i,1] <- 0",
                   "\n\t\tdelta[i,1] <- 0",
                   "\n\t\tfor(k in 2:na[i]) {",
                   "\n\t\t\tVar[i,k] <- pow(se[i,k], 2)",
                   "\n\t\t\tprec[i,k] <- 1/Var[i,k]",
                   "\n\t\t}",
                   "\n\t\tfor(k in 2:na[i]) {",
                   "\n\t\t\tdelta[i,k] ~ dnorm(md[i,k], taud[i,k])",
                   "\n\t\t\tmd[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k]",
                   "\n\t\t\ttaud[i,k] <- tau * 2 * (k-1)/k",
                   "\n\t\t\tw[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]])",
                   "\n\t\t\tsw[i,k] <- sum(w[i,1:(k-1)])/ (k-1)",
                   "\n\t\t}",
                   "\n\t}",
                   "\n\ttotresdev <- sum(resdev[])",
                   "\n\td[1] <- 0",
                   "\n\tfor(k in 2:", ntreat, ") {",
                   "\n\t\td[k] ~ dnorm(0,.0001)",
                   "\n\t}",
                   "\n\tsd ~ dunif(0,5)",
                   "\n\ttau <- pow(sd,-2)",
                   "\n}")
    return(code)
  })
}

#' Run the model using the network object
#' 
#' This is similar to the function \code{\link{network.data}}, except it uses contrast-level data instead of arms-level data.
#'
#' @param network contrast level network object created from \code{\link{contrast.network.data}} function
#' @param inits Initial values for the parameters being sampled. If left unspecified, program will generate reasonable initial values.
#' @param n.chains Number of chains to run
#' @param max.run Maximum number of iterations that user is willing to run. If the algorithm is not converging, it will run up to \code{max.run} iterations before printing a message that it did not converge
#' @param setsize Number of iterations that are run between convergence checks. If the algorithm converges fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run Final number of iterations that the user wants to store. If after the algorithm converges, user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run extra iterations to reach the specified number of runs
#' @param conv.limit Convergence limit for Gelman and Rubin's convergence diagnostic. Point estimate is used to test convergence of parameters for study effect (eta), relative effect (d), and heterogeneity (log variance (logvar)).
#' @param extra.pars.save Parameters that user wants to save besides the default parameters saved. See code using \code{cat(network$code)} to see which parameters can be saved.
#' @return
#' \item{data_rjags}{Data that is put into rjags function \code{\link{contrast.jags.model}}}
#' \item{inits}{Initial values that are either specified by the user or generated as a default}
#' \item{pars.save}{Parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{burnin}{Half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{MCMC samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{Maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' \item{deviance}{Contains deviance statistics such as pD (effective number of parameters) and DIC (Deviance Information Criterion)}
#' \item{rank.tx}{Rank probability calculated for each treatments. \code{rank.preference} parameter in \code{\link{network.data}} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' @export

contrast.network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                        conv.limit = 1.05, extra.pars.save = NULL){
  
  if (!inherits(network, "contrast.network.data")) {
    stop('Given network is not contrast.network.data. Run contrast.network.data function first')
  }
  
  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }
  
  with(network, {
    
    data <- list(y = Outcomes, t = Treat, se = SE, na = na)
    
    if(!is.null(V)){
      data$V <- V
    }
    
    pars.save <- c("d", "sd")
    
    pars.save <- c(pars.save, "totresdev", "delta")
    
    if(!is.null(extra.pars.save)) {
      extra.pars.save.check(extra.pars.save)
      pars.save <- c(pars.save, extra.pars.save)
    }
    
    # if(is.null(inits)){
    #   inits <- network.inits(network, n.chains)
    # }
    samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)
    result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
    result <- c(result, samples)
    
    # if(dic == TRUE){
    #   result$deviance <- calculate.deviance(result)
    # }
    # result$rank.tx <- rank.tx(result)
    class(result) <- "contrast.network.result"
    return(result)
  })
}

#' Summarize result run by \code{\link{contrast.network.run}}
#'
#' This function uses summary function in coda package to summarize mcmc.list object. Monte carlo error (Time-series SE) is also obtained using the coda package and is printed in the summary as a default.
#'
#' @param object Result object created by \code{\link{contrast.network.run}} function
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network) 
#' summary(result)
#' @export

summary.contrast.network.result <- function(object){
  
  if(!inherits(object, "contrast.network.result")) {
    stop('This is not the output from contrast.network.run. Need to run contrast.network.run function first')
  }
  rval <- list("summary.samples"= summary(object$samples),
               "deviance" = unlist(object$deviance[1:3]),
               "total_n" = sum(object$network$na))
  class(rval) <- 'summary.contrast.network.result'
  rval
}


#' Plot traceplot and posterior density of the result using contrast data
#'
#' This function uses plotting function in coda package to plot mcmc.list object
#'
#' @param x Result object created by \code{\link{contrast.network.run}} function
#' @examples
#' network <- with(parkinsons_contrast, {
#'  contrast.network.data(Outcomes, Treat, SE, na, V)
#' })
#' result <- contrast.network.run(network)
#' plot(result)
#' @export

plot.contrast.network.result <- function(x) {

  if(!inherits(x, "contrast.network.result")) {
    stop('This is not the output from contrast.network.run. Need to run contrast.network.run function first')
  }
  plot(x$samples)
}


#' Find deviance statistics such as DIC and pD.
#'
#' Calculates deviance statistics. This function automatically called in \code{\link{contrast.network.run}} and the deviance statistics are stored after sampling is finished.
#'
#' @param result Object created by \code{\link{contrast.network.run}} function
#' @return
#' \item{Dbar}{Overall residual deviance}
#' \item{pD}{Sum of leverage_arm (i.e. total leverage)}
#' \item{DIC}{Deviance information criteria (sum of Dbar and pD)}
#' @examples
#' #parkinsons
#' network <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' result <- network.run(network)
#' calculate.deviance(result)
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

calculate.contrast.deviance <- function(result){

  
  network <- result$network
  samples <- result$samples
  
  totresdev <- lapply(samples, function(x){ x[,"totresdev"]})
  Dbar <- mean(unlist(totresdev))
  
  ybar <- lapply(samples, function(x){ x[,grep("delta\\[", dimnames(samples[[1]])[[2]])] })
  ybar <- do.call(rbind, ybar)
  ybar <- apply(ybar, 2, mean)
  ybar_arm <- devtilda_arm <- matrix(NA, nrow = sum(network$na_count), ncol = max(network$na))  
  
  with(network, {
    # Find Omega..
    Sigma <- array(NA, c(sum(na_count), max(na) -1, max(na) -1 ))
    for(i in 1:sum(na_count)){
      for(k in 1:(na[i] - 1)){
        for(j in 1:(na[i] - 1)){
          Sigma[i,j,k] <- V[i] * (1- (j ==k)) + SE[i, k+1] * (j==k)
        }
      }
    }
    Omega <- Sigma
    for(i in 1:sum(na_count)){
      Omega[i,,] <- solve(Sigma[i,,])  
    }
    
    # 2 arm
    for(i in 1:na_count[1]){
      for(j in 2:na[i]){
        r_value <- Outcomes[i,j]
        se_value <- SE[i,j]
        ybar_arm[i,j] <- ybar[which(paste("delta[", i, ",", j, "]", sep = "") == names(ybar))]
        devtilda_arm[i,j] <- ifelse(se_value != 0, (r_value - ybar_arm[i,j])^2 / se_value^2, 0)
      }
    }
    
    # 3 or more
    if(length(na_count) > 1){
      for(ii in 2:length(na_count)){
        
        for(i in cumsum(na_count)[ii-1]: cumsum(na_count)[ii]){
          for(j in 2:na[i]){
            r_value <- Outcomes[i,j]
            omega_value <- Omega[i,,]
            ybar_arm[i,j] <- ybar[which(paste("delta[", i, ",", j, "]", sep = "") == names(ybar))]
            devtilda_arm[i,j] <- ifelse(all(se_value != 0), (r_value - ybar_arm[i,j]) %*% omega_value %*% (r_value - ybar_arm[i,j]), 0)
          }
        }   
      }  
    }
    return(devtilda_arm)
    
  })
  
  
  #for()
  
    
    if(length(na_count) > 1){
      
      for(i in 2:length(na_count)){
        
        code <- paste0(code, "\n\tfor(i in ", cumsum(na_count)[i-1] + 1, ":", cumsum(na_count)[i], ") {", 
                       "\n\t\tfor(k in 1:(na[i]-1)) {",
                       "\n\t\t\tfor(j in 1:(na[i]-1)) {",
                       "\n\t\t\t\tSigma[i,j,k] <- V[i]*(1-equals(j,k)) + Var[i,k+1] * equals(j,k)",
                       "\n\t\t\t}",
                       "\n\t\t}",
                       "\n\t\tOmega[i,1:(na[i]-1),1:(na[i]-1)] <- inverse(Sigma[i,,])",
                       "\n\t\ty[i,2:na[i]] ~ dmnorm(delta[i,2:na[i]], Omega[i, 1:(na[i]-1), 1:(na[i]-1)])",
                       "\n\t\tfor(k in 1:(na[i]-1)){",
                       "\n\t\t\tydiff[i,k] <- y[i,(k+1)] - delta[i,(k+1)]",
                       "\n\t\t\tz[i,k] <- inprod(Omega[i,k,1:(na[i]-1)], ydiff[i,1:(na[i]-1)])",
                       "\n\t\t}",
                       "\n\t\tresdev[i] <- inprod(ydiff[i,1:(na[i]-1)], z[i, 1:(na[i]-1)])",
                       "\n\t}")
        
    
  
  for(i in 1:sum(network$na_count)){
    for(j in 2:network$na[i]){
      r_value <- network$Outcomes[i,j]
      se_value <- network$SE[i,j]
      ybar_arm[i,j] <- ybar[which(paste("delta[", i, ",", j, "]", sep = "") == names(ybar))]
      devtilda_arm[i,j] <- ifelse(se_value != 0, (r_value - ybar_arm[i,j])^2 / se_value^2, 0)
    }
  }

}