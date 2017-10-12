#' Run the model using the network object
#'
#' @param network network object created from \code{network.data} function
#' @param inits initial values for the parameters being sampled. If left unspecified, program will generate reasonable initial values.
#' @param n.chains number of chains to run
#' @param max.run maximum number of iterations that user is willing to run. If the algorithm is not converging, it will run up to max.run iterations before printing a message that it did not converge
#' @param setsize number of iterations that are run between convergence checks. If the algorithm converges fast, user wouldn't need a big setsize. The number that is printed between each convergence checks is the gelman-rubin diagnostics and we would want that to be below the conv.limit the user specifies.
#' @param n.run final number of iterations that the user wants to store. If after the algorithm converges, user wants less number of iterations, we thin the sequence. If the user wants more iterations, we run extra iterations to reach the specified number of runs
#' @param conv.limit convergence limit for Gelman and Rubin's convergence diagnostic. Point estimate is used to test convergence of parameters for study effect (eta), relative effect (d), and heterogeneity (log variance (logvar)).
#' @param extra.pars.save parameters that user wants to save besides the default parameters saved.
#' @return
#' \item{data_rjags}{data that is put into rjags function \code{jags.model}}
#' \item{inits}{initial values that are either specified by the user or generated as a default}
#' \item{pars.save}{parameters that are saved. Add more parameters in extra.pars.save if other variables are desired}
#' \item{burnin}{half of the converged sequence is thrown out as a burnin}
#' \item{n.thin}{If the number of iterations user wants (n.run) is less than the number of converged sequence after burnin, we thin the sequence and store the thinning interval}
#' \item{samples}{mcmc samples stored using jags. The returned samples have the form of mcmc.list and can be directly applied to coda functions}
#' \item{max.gelman}{maximum Gelman and Rubin's convergence diagnostic calculated for the final sample}
#' \item{deviance}{contains deviance statistics such as pD (effective number of parameters) and DIC (Deviance Information Criterion)}
#' \item{rank.tx}{rank probability calculated for each treatments. Rank.preference parameter in \code{network.data} is used to define whether higher or lower value is preferred. The numbers are probabilities that a given treatment has been in certain rank in the sequence.}
#' @examples
#' #parkinson's example (normal)
#' parkinsons
#' network <- with(parkinsons,{
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' result <- network.run(network)
#' @export

network.run <- function(network, inits = NULL, n.chains = 3, max.run = 100000, setsize = 10000, n.run = 50000,
                        conv.limit = 1.05, extra.pars.save = NULL){

  if (!inherits(network, "network.data")) {
    stop('Given network is not network.data. Run network.data function first')
  }

  if(max.run < setsize){
    stop("setsize should be smaller than max.run")
  }

  data <- with(network,
    if(response == "binomial"){
      list(na = na, t = t, r = r, n = n)
    } else if(response == "normal"){
      list(na = na, t = t, r = r, se = se)
    } else if(response == "multinomial"){
      list(na = na, t = t)
    })

  if(network$response == "multinomial"){
    if(is.null(network$miss.matrix)){
      data$r <- network$r
      data$n <- network$n
    }
    for(i in 1:length(network$miss.patterns[[1]])){
      data[[paste("r",i,sep="")]] <- network[[paste("r",i,sep="")]]
      data[[paste("n",i,sep="")]] <- network[[paste("n",i,sep="")]]
    }
  }
  data <- c(data, network$prior.data)

  # add covariate info
  if(!is.null(network$covariate)){
    for(i in seq(dim(network$covariate)[2])){
      data[[paste("mx",i, sep = "")]] = network[[paste("mx",i, sep = "")]]
      data[[paste("x",i, sep = "")]] = network[[paste("x",i, sep = "")]]
    }
  }

  # add baseline info
  if(network$baseline != "none"){
    data$mx_bl = network$mx_bl
  }

  ########## parameters to save in the model
  pars.save <- with(network,
    if(network$response == "binomial" || network$response == "normal"){
      c("Eta", "d", "sd", "logvar","prob")
    } else if(network$response == "multinomial"){
      c("Eta", "d", "sigma", "sigma_transformed","prob")
    })
  if(network$type == "fixed"){
    pars.save <- pars.save[!pars.save %in% c("sd", "sigma", "logvar", "sigma_transformed")]
  }

  if(!is.null(extra.pars.save)) {
    extra.pars.save.check(extra.pars.save)
    pars.save <- c(pars.save, extra.pars.save)
  }
  if(network$dic == TRUE){
    pars.save <- c(pars.save, "totresdev")
    if(network$response == "binomial"){
      pars.save <- c(pars.save, "rhat", "dev")
    } else if(network$response == "normal"){
      pars.save <- c(pars.save, "theta", "dev")
    } else if(network$response == "multinomial"){
      if(is.null(network$miss.matrix)){
        pars.save <- c(pars.save, "rhat", "dev")
      } else{
        for(i in 1:network$npattern){
          pars.save <- c(pars.save, paste0("rhat", i), paste0("dev", i))
        }
      }
    }
  }
  if(network$baseline != "none"){
    pars.save = c(pars.save, "b_bl")
    if(network$baseline %in% c("common", "exchangeable")){
      pars.save <- c(pars.save, "B")
    }
  }
  if(!is.null(network$covariate)){
    for(i in seq(dim(network$covariate)[2])){
      pars.save = c(pars.save, paste("beta",i,sep = ""))
    }
  }
  pars.save <- unique(pars.save)

  if(is.null(inits)){
    inits <- network.inits(network, n.chains)
  }
  
  samples <- jags.fit(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit)

  result <- list(network = network, data.rjags = data, inits = inits, pars.save = pars.save)
  result <- c(result, samples)

  if(network$dic == TRUE){
    result$deviance <- calculate.deviance(result)
  }
  result$rank.tx <- rank.tx(result)

  class(result) <- "network.result"
  result
}


jags.fit <- function(network, data, pars.save, inits, n.chains, max.run, setsize, n.run, conv.limit) {

  mod = rjags::jags.model(textConnection(network$code), data = data, inits = inits, n.chains = n.chains, n.adapt = 0)

  adapted <- FALSE
  count <- 0
  while(!adapted){
    adapted <- rjags::adapt(mod, setsize, end.adaptation = FALSE)
    count <- count + 1
    if(count == 100){
      stop("algorithm has not adapted")
    }
  }
  
  conv.save <- if(network$response == "multinomial"){
    c("d", "Eta", "sigma_transformed")
  } else if(network$response == "binomial" || network$response == "normal"){
    c("d", "Eta", "logvar")
  }
  if(network$type == "fixed"){
    conv.save <- conv.save[!conv.save %in% c("logvar", "sigma_transformed")]
  }

  samples <- rjags::coda.samples(model = mod, variable.names = pars.save, n.iter = setsize)
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  conv.save.variables <- varnames.split %in% conv.save

  max.gelman <- find.max.gelman(samples, conv.save.variables)
  print(max.gelman)
  check <- max.gelman > conv.limit

  if(check) {
    count <- 1
    while (check & count < max.run/setsize) {
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = setsize)
      samples <- add.mcmc(samples, samples2)

      count <- count + 1

      max.gelman <- find.max.gelman(samples, conv.save.variables)
      check <- max.gelman > conv.limit
      print(max.gelman)
    }
  }

  start <- mcpar(samples[[1]])[1]
  end <- mcpar(samples[[1]])[2]
  mid <- (end + start-1)/2
  burnin <- ceiling(end - mid)
  samples <- window(samples, mid+1, end, 1) #keep the last half of the converged sequence
  samples <- new.mcmc(samples)

  n.thin <- 1
  if(check == TRUE){
    print("code didn't converge according to gelman-rubin diagnostics")
  } else if(n.run < burnin){
    n.thin <- ceiling(burnin/n.run)
    extra.run <- n.run * n.thin - burnin
    if(extra.run != 0){
      samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
      samples <- add.mcmc(samples, samples2)
    }
    samples <- window(samples, 1, dim(samples[[1]])[1], n.thin)
  } else if(n.run > burnin){
    extra.run <- n.run - burnin
    samples2 <- rjags::coda.samples(mod, variable.names = pars.save, n.iter = extra.run)
    samples <- add.mcmc(samples, samples2)
  }
  max.gelman <- find.max.gelman(samples, conv.save.variables)
  print(max.gelman)

  out <-list(burnin = burnin, n.thin = n.thin, samples = samples, max.gelman = max.gelman)
  return(out)
}

extra.pars.save.check <- function(extra.pars.save){
  if(!is.atomic(extra.pars.save) || !is.vector(extra.pars.save)) stop("extra pars should be a vector of strings")
  for(i in 1:length(extra.pars.save)){
    if(!is.character(extra.pars.save[i])) stop("extra pars should be a vector of strings")
  }
}

new.mcmc <- function(x){
  n.chains <- length(x)
  n.var <- nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- x[[i]]
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

add.mcmc <- function(x, y){

  n.chains <- length(x)
  n.var <- nvar(x)
  newobjects <- vector("list", length = n.chains)

  for(i in 1:n.chains){
    newobjects[[i]] <- matrix(NA, nrow = 0, ncol = n.var, dimnames = list(NULL, dimnames(x[[1]])[[2]]))
    newobjects[[i]] <- rbind(x[[i]], y[[i]])
    newobjects[[i]] <- mcmc(newobjects[[i]])
  }
  mcmc.list(newobjects)
}

find.max.gelman <- function(samples, index){

  samples2 <- lapply(samples, function(x){ x[,index]})
  samples2 <- lapply(samples2, function(x) { x[,colSums(abs(x)) != 0] })

  max(gelman.diag(samples2, multivariate = FALSE)$psrf[,1]) #look at point estimate instead of 95% C.I.
}

