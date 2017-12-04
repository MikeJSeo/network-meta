pick.summary.variables <- function(result, extra.pars = NULL, only.pars = NULL){
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
    pars <- c("d", "sd", "sigma", "b_bl", "beta", "B", "sdB")
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

#' summarize result run by \code{network.run}
#'
#' Use summary function in \code{coda} to summarize mcmc.list object
#'
#' @param object result object created by \code{network.run} function
#' @param ... additional arguments affecting the summary produced
#' #' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' result <- network.run(network)
#' summary(result)
#' @export

summary.network.result <- function(object, ...){

  if(!inherits(object, "network.result")) {
    stop('This is not the output from network.run. Need to run network.run function first')
  }
  summary.samples <- pick.summary.variables(object, ...)

  rval <- list("summary.samples"= summary(summary.samples),
               "Treat.order" =  object$network$Treat.order,
               "deviance" = unlist(object$deviance[1:3]),
               "total_n" = sum(object$network$na))
  class(rval) <- 'summary.network.result'
  rval
}

#' plot traceplot and posterior density of the result
#'
#' Use plotting function in \code{coda} to plot mcmc.list object
#'
#' @param x result object created by \code{network.run} function
#' @param ... additional arguments affecting the plot produced
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' result <- network.run(network)
#' plot(result, only.pars = "sd")
#' @export

plot.network.result <- function(x, ...) {
  summary.samples <- pick.summary.variables(x, ...)
  plot(summary.samples)
}


#' Use coda package to plot gelman-diagnostic plot
#'
#' @param result object created by \code{network.run} function
#' @param extra.pars extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars parameters that user wants to plot only
#' @examples
#' #blocker
#' network <- with(blocker,{
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- network.run(network)
#' network.gelman.plot(result, only.pars = "d")
#' @export

network.gelman.plot <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  for(v in 1:nvar(summary.samples)){
    gelman.plot(summary.samples[,v,drop=FALSE])
  }
}

#' Use coda package to find gelman-diagnostic diagnostics
#'
#' @param result object created by \code{network.run} function
#' @param extra.pars extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars parameters that user wants to plot only
#' @examples
#' network <- with(statins, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial",
#'  Treat.order = c("Placebo", "Statin"), covariate = covariate, covariate.type = "discrete")
#' })
#' result <- network.run(network)
#' network.gelman.diag(result, extra.pars = "Eta")
#' @export

network.gelman.diag <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  gelman.diag(summary.samples, multivariate = FALSE)
}

#' Use coda package to find autocorrelation diagnostics
#'
#' @param result object created by \code{network.run} function
#' @param lags a vector of lags at which to calculate the autocorrelation
#' @param extra.pars extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars parameters that user wants to plot only
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' 
#' result <- network.run(network)
#' network.autocorr.diag(result, only.pars = "d")
#' @export

network.autocorr.diag <- function(result, lags = c(0,1,5,10,50), extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))

  autocorr.diag(summary.samples, lags = lags)
}

#' Use coda package to plot autocorrelation plot
#'
#' @param result object created by \code{network.run} function
#' @param extra.pars extra parameters that the user wants to plot other than the default parameters.
#' @param only.pars parameters that user wants to plot only
#' @examples
#' #cardiovascular
#' Study <- cardiovascular[["Study"]]
#' Treat <- cardiovascular[["Treat"]]
#' Outcomes <- cardiovascular[["Outcomes"]]
#' N <- cardiovascular[["N"]]
#' network <- with(cardiovascular, {
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' result <- network.run(network)
#' network.autocorr.plot(result, only.pars = "d")
#' @export

network.autocorr.plot <- function(result, extra.pars = NULL, only.pars = NULL){
  summary.samples <- pick.summary.variables(result, extra.pars, only.pars)
  summary.samples <- mcmc.list(lapply(summary.samples, function(x) { x[,colSums(abs(x)) != 0] }))
  autocorr.plot(summary.samples)
}

#' Find relative effects for different base treatment and comparison treatments
#'
#' @param result object created by \code{network.run} function
#' @param base.treatment base treatment user wants for the relative effects. Base treatment is initially set by \code{Treat.order} parameter in \code{network.data} (first one in the list). If set to null, default is to use base treatment.
#' @param comparison.treatments treatments that user wants to compare against base treatment. If set to null, all the treatments besides base treatment is considered as comparison treatments.
#' @param base.category base category user wants for the relative effects. This is only used for multinomial data.
#' @param comparison.categories category that user wants to compare against base.category
#' @param covariate covariate value at which to compute relative effects.
#' @return
#' This returns a mcmc.list sample of relative effects for the base treatment specified. This allows user to obtain relative effects of different base.treatment after the sampling has been done.
#' For a simple summary use relative.effects.table.
#' @examples
#' #We can fit two different models with different base treatment and we can
#' #obtain same relative effects estimate using this function
#' #parkinsons
#' network <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' result <- network.run(network) 
#' summary(result)
#'
#' network2 <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal",
#'  Treat.order = c(2,1,3,4,5))
#' })
#' result2 <- network.run(network2)
#'
#' summary(result)
#' summary(relative.effects(result2, base.treatment = 1))
#'
#' #This also works for comparing different base.category for multinomial.
#' #We fit two different models and compare the estimates again.
#' #cardiovascular
#' 
#' network3 <- with(cardiovascular, {
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' result3 <- network.run(network3)
#' 
#' network4 <- with(cardiovascular, {
#'  network.data(Outcomes[,c(2,1,3)], Study, Treat, N, response = "multinomial")
#' })
#' result4 <- network.run(network4)
#'
#' summary(result3)
#' summary(relative.effects(result4, base.category = 2))
#' @seealso \code{\link{relative.effects.table}}
#' @export

relative.effects <- function(result, base.treatment = NULL, comparison.treatments = NULL, base.category = NULL, comparison.categories = NULL, covariate = NULL){

  network <- result$network

  if(!is.null(covariate)){
    stopifnot(length(covariate) == dim(network$covariate)[2])
  }

  Treat.order <- network$Treat.order
  if(!is.null(base.treatment)){
    stopifnot(base.treatment %in% Treat.order)
  } else{
    base.treatment <- Treat.order[1]
  }
  if(!is.null(comparison.treatments)){
    stopifnot(comparison.treatments %in% Treat.order)
    stopifnot(!comparison.treatments %in% base.treatment)
  } else{
    comparison.treatments <- Treat.order[-which(Treat.order == base.treatment)]
  }
  if(!is.null(covariate)){
    summary.samples <- pick.summary.variables(result, only.pars = c("d", "beta"))
  } else{
    summary.samples <- pick.summary.variables(result, only.pars = c("d"))
  }
  vars <- dimnames(summary.samples[[1]])[[2]]

  if(network$response != "multinomial"){
    effects <- matrix(0, nrow = network$ntreat, ncol = length(comparison.treatments))
    effects[which(Treat.order == base.treatment),] = -1

    col_name = NULL
    for(i in 1:ncol(effects)){
      effects[which(comparison.treatments[i] == Treat.order),i] = 1
      col_name <- c(col_name, paste0("d_treatment", base.treatment, comparison.treatments[i]))
    }

    if(!is.null(covariate)){
      cov_matrix <-  covariate_centerered  <- NULL
      for(i in 1:length(covariate)){
        cov <- effects
        covariate_centered <- covariate[i] - network[[paste0("mx",i)]]
        cov <- cov * covariate_centered
        cov_matrix <- rbind(cov_matrix, cov)
      }
      effects <- rbind(cov_matrix, effects)
    }
    colnames(effects) <- col_name
    rownames(effects) <- vars

    samples <- as.mcmc.list(lapply(summary.samples, function(chain){
      samples <- chain %*% effects
      colnames(samples) <- colnames(effects)
      mcmc(samples, start = start(chain), end = end(chain), thin = thin(chain))
    }))
  } else{
    vars_d <- vars[grep("d\\[", vars)]
    categories_row <- as.numeric(substr(vars_d, nchar(vars_d[1])-1, nchar(vars_d[1])-1))
    categories_row <- categories_row+1
    ncat <- network$ncat

    if(!is.null(base.category)){
      stopifnot(base.category %in% 1:ncat)
    } else{
      base.category <- 1
    }
    if(!is.null(comparison.categories)){
      stopifnot(comparison.categories %in% 1:ncat)
      stopifnot(!comparison.categories %in% base.category)
    } else{
      comparison.categories <- (1:ncat)[-base.category]
    }

    effects <- matrix(0, nrow = network$ntreat*(network$ncat-1), length(vars), ncol = length(comparison.treatments) * length(comparison.categories))
    categories_column <- rep(comparison.categories, each = length(comparison.treatments))

    effects[which(rep(Treat.order, ncat-1) == base.treatment),] <- -1
    col_name <- NULL
    for(i in 1:ncol(effects)){
      effects[which(rep(Treat.order, ncat-1) == rep(comparison.treatments, length(comparison.categories))[i]),i] <- 1
      col_name <- c(col_name, paste0("d_treatment", base.treatment, rep(comparison.treatments, length(comparison.categories))[i]))
    }
    colnames(effects) <- col_name

    for(i in 1:ncol(effects)){
      effects[which(categories_row == base.category),i] <- -effects[which(categories_row == base.category),i]
      effects[which(categories_row != base.category & categories_row != rep(comparison.categories, each = length(comparison.treatments))[i]),i] <- 0
      colnames(effects)[i] <- paste0(colnames(effects)[i], "_category", base.category, rep(comparison.categories, each = length(comparison.treatments))[i])
    }

    if(!is.null(covariate)){
      cov_matrix <-  covariate_centerered  <- NULL
      for(i in 1:length(covariate)){
        cov <- effects
        covariate_centered <- covariate[i] - network[[paste0("mx",i)]]
        cov <- cov * covariate_centered
        cov_matrix <- rbind(cov_matrix, cov)
      }
      effects <- rbind(cov_matrix, effects)
    }
    rownames(effects) <- vars

    samples <- as.mcmc.list(lapply(summary.samples, function(chain){
      samples <- chain %*% effects
      colnames(samples) <- colnames(effects)
      mcmc(samples, start = start(chain), end = end(chain), thin = thin(chain))
    }))
  }
  samples
}

#' Make a summary table for relative effects
#' 
#' Relative effects in units of log odds ratio for binomial and multinomial data and real number scale for normal data.
#'
#' @param result object created by \code{network.run} function
#' @param summary_stat specifies what type of statistics user wants. Options are: "mean", "ci", "quantile", "sd", "p-value".
#' ci gives 95% confidence interval (0.025, 0.5, 0.975) and quantile gives specific quantile specified in probs parameter. 
#' P-value is the probability relative effect (in binomial, log odds ratio) is less than 0.
#' @param base.category specifies for which base category user wants for the summary. Used only for multinoimal.
#' @examples
#' #cardiovascular
#' network <- with(cardiovascular,{
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' result <- network.run(network)
#' exp(relative.effects.table(result)) #look at odds ratio instead of log odds ratio
#' @seealso \code{\link{relative.effects}}
#' @export

relative.effects.table <- function(result, summary_stat = "mean", probs = NULL, base.category = NULL){

  stopifnot(summary_stat %in% c("mean", "quantile", "sd", "p-value"))

  if(!is.null(probs)){
    if(length(probs) == 1){
      stop("length of probs should be 1")
    }
  }
  
  Treat.order <- result$network$Treat.order

  ts <- 1:length(Treat.order)
  comps <- combn(ts, 2)

  if(result$network$response != "multinomial"){
    tbl <- matrix(NA, nrow = length(ts), ncol = length(ts), dimnames = list(Treat.order, Treat.order))

    for (i in 1:ncol(comps)) {
      comp <- comps[, i]
      samples <- as.matrix(relative.effects(result, base.treatment = Treat.order[comp[1]], comparison.treatments = Treat.order[comp[2]]))

      if(summary_stat == "mean"){
        tbl[comp[1], comp[2]] <- mean(samples)
        tbl[comp[2], comp[1]] <- -tbl[comp[1], comp[2]]
      } else if(summary_stat == "ci"){
        q <- round(quantile(samples, probs = c(0.025, 0.5, 0.975)), 6)
        tbl[comp[1], comp[2]] <- paste0("[", q[1], ",", q[2], ",", q[3], "]")
        tbl[comp[2], comp[1]] <- paste0("[", -q[3], ",", -q[2], ",", -q[1], "]")
      } else if(summary_stat == "quantile"){
        tbl[comp[1], comp[2]] <- round(quantile(samples, probs = probs), 6)
        tbl[comp[2], comp[1]] <- -tbl[comp[1], comp[2]]
      } else if(summary_stat == "sd"){
        tbl[comp[1], comp[2]] <- tbl[comp[2], comp[1]] <- sd(samples)
      } else if(summary_stat == "p-value"){
        tbl[comp[1], comp[2]] <- sum(samples < 0)/ dim(samples)[1]
        tbl[comp[2], comp[1]] <- 1 - tbl[comp[1], comp[2]]
      }
    }
  } else if(result$network$response == "multinomial"){
    ncat <- result$network$ncat
    tbl <- array(NA, dim = c(length(ts), length(ts), ncat -1), dimnames = list(Treat.order, Treat.order, NULL))

    for (i in 1:ncol(comps)) {
      comp <- comps[, i]
      samples <- as.matrix(relative.effects(result, base.treatment = Treat.order[comp[1]], comparison.treatments = Treat.order[comp[2]], base.category = base.category))

      if(summary_stat == "mean"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, mean)
        tbl[comp[2], comp[1],]  <- -tbl[comp[1], comp[2],]
      } else if(summary_stat == "ci"){
        q <- round(apply(samples, 2, quantile, probs = c(0.025, 0.5, 0.975)), 6)
        q1 <- apply(q, 2, function(x){ paste0("[", x[1], ",", x[2], ",", x[3], "]")})
        q2 <- apply(q, 2, function(x){ paste0("[", -x[3], ",", -x[2], ",", -x[1], "]")})
        tbl[comp[1], comp[2],] <- q1
        tbl[comp[2], comp[1],] <- q2
      } else if(summary_stat == "quantile"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, quantile, probs = probs)
        tbl[comp[2], comp[1],] <- -tbl[comp[1], comp[2],]
      } else if(summary_stat == "sd"){
        tbl[comp[1], comp[2],] <- tbl[comp[2], comp[1],] <- apply(samples, 2, sd)
      } else if(summary_stat == "p-value"){
        tbl[comp[1], comp[2],] <- apply(samples, 2, function(x){ sum(x <0) / length(x)})
        tbl[comp[2], comp[1],] <- 1 - tbl[comp[1], comp[2],]
      }
    }
  }
  tbl
}

#' Creates a treatment rank table
#'
#' @param result object created by \code{network.run} function
#' @return
#' This makes a table of ranking for each treament. Each number in the cell represents a probability certain treatment was in such rank.
#' This table is also stored as an output from \code{network.run}.
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- network.run(network)
#' rank.tx(result)
#' @seealso \code{\link{network.rank.tx.plot}}
#' @export

rank.tx <- function(result){
  samples <- result[["samples"]]
  varnames <- dimnames(samples[[1]])[[2]]
  varnames.split <- sapply(strsplit(varnames, "\\["), '[[', 1)
  varnames.split <- gsub("[[:digit:]]","",varnames.split)

  rank.samples <- lapply(samples, function(x){x[,varnames.split %in% "prob"]})

  Treat.order <- result$network$Treat.order
  response <- result$network$response

  if(response != "multinomial"){
    prob.matrix <- matrix(NA, nrow = length(Treat.order), ncol = length(Treat.order), dimnames = list(paste0("rank ", 1:length(Treat.order)), paste0("treatment ", Treat.order)))
    for(i in 1:nrow(prob.matrix)){
      for(j in 1:ncol(prob.matrix)){
        prob.matrix[i,j] <- mean(unlist(lapply(rank.samples, function(x){ x[,paste0("prob[", i, ",", j, "]")]})))
      }
    }
  } else if(response == "multinomial"){
    ncat <- result$network$ncat
    prob.matrix <- array(NA, dim = c(length(Treat.order), length(Treat.order), ncat-1), dimnames = list(paste0("rank ", 1:length(Treat.order)), paste0("treatment ", Treat.order),  paste0("Category ", 1:(ncat-1))))
    for(i in 1:nrow(prob.matrix)){
      for(j in 1:ncol(prob.matrix)){
        for(k in 1:(ncat-1)){
          prob.matrix[i,j,k] <- mean(unlist(lapply(rank.samples, function(x){ x[,paste0("prob[", i, ",", j, ",", k, "]")]})))
        }
      }
    }
  }
  return(prob.matrix)
}

#' Creates a treatment rank plot
#'
#' @param result object created by \code{network.run} function
#' @param txnames treatment names used in creating legend
#' @param catnames category names. Only used in multinomial.
#' @param legend.position x,y position of the legend
#' @examples
#' network <-with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- network.run(network)
#' network.rank.tx.plot(result, txnames = c("a", "b"))
#' @seealso \code{\link{rank.tx}}
#' @export

network.rank.tx.plot <- function(result, txnames = NULL, catnames = NULL, legend.position = c(1,1)){

  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),pty="s",yaxt="n",ylab="Probability",xlab="Rank")
    axis(side=1,at=seq(ntreat))
    axis(side=2,at=seq(0,1,by=0.2))
    for (i in seq(ntreat)) {
      points(seq(ntreat), rank.table[,i],type="b",lty=i,col=i,pch=i)
    }
    legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))
    for (j in seq(ncat)) {
      plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),pty="s",yaxt="n",ylab="Probability",xlab="Rank")
      axis(side=1,at=seq(ntreat))
      axis(side=2,at=seq(0,1,by=0.2))
      title(catnames[j])
      for (i in seq(ntreat)) {
        points(seq(ntreat), rank.table[,i,j],type="b",lty=i,col=i,pch=i)
      }
      legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
    }
  }
}

#' Creates a treatment cumulative rank plot
#'
#' @param result object created by \code{network.run} function
#' @param txnames treatment names used in creating legend
#' @param catnames category names. Only used in multinomial.
#' @param legend.position x,y position of the legend
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- network.run(network)
#' network.cumrank.tx.plot(result, txnames = c("control", "beta blocker"))
#' @seealso \code{\link{rank.tx}}
#' @export

network.cumrank.tx.plot <- function(result, txnames = NULL, catnames = NULL, legend.position = c(1,1)){
  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    x <- apply(rank.table,2,cumsum)
    plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),yaxt="n",ylab="Cumulative Probability",xlab="Rank")
    axis(side=1,at=seq(ntreat))
    axis(side=2,at=seq(0,1,by=0.2))
    for (j in seq(ntreat))
      points(seq(ntreat), x[,j],type="l",lty=j,col=j)
    legend(legend.position[1], legend.position[2], txnames,lty=1:(ntreat),bty="n",cex=.75,col=1:(ntreat))
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))

    for (i in seq(ncat))  {
      x = apply(rank.table[,,i],2,cumsum)
      plot(seq(ntreat),seq(ntreat),type="n",xaxt="n",ylim=c(0,1),yaxt="n",ylab="Cumulative Probability",xlab="Rank")
      axis(side=1,at=seq(ntreat))
      axis(side=2,at=seq(0,1,by=0.2))
      title(catnames[i])
      for (j in seq(ntreat))
        points(seq(ntreat), x[,j],type="l",lty=j,col=j)
      legend(legend.position[1], legend.position[2],txnames,lty=1:ntreat,bty="n",cex=.75,col=1:ntreat)
    }
  }
}

#' Creates a treatment rank plot
#'
#' SUCRA is the surface under the cumulative ranking distribution defined in Salanti et al. (2011)
#'
#' @param result object created by \code{network.run} function
#' @param txnames treatment names used in creating legend
#' @param catnames category names. Only used in multinomial.
#' @examples
#' ########### certolizumab (with baseline risk)
#' network <- with(certolizumab, {
#'  network.data(Outcomes, Study, Treat, N=N, response = "binomial", Treat.order,
#'  baseline = "common", hy.prior = list("dhnorm", 0, 9.77))
#' })
#' result <- network.run(network)
#' sucra(result)
#' @seealso \code{\link{rank.tx}}
#' @references G. Salanti, A.E. Ades, J.P.A. Ioannidisa (2011), \emph{Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial}, Journal of Clinical Epidemiology 64(2):163-71. [\url{https://doi.org/10.1016/j.jclinepi.2010.03.016}]
#' @export

sucra = function(result, txnames = NULL, catnames = NULL)
{
  rank.table <- rank.tx(result)
  ntreat = dim(rank.table)[1]
  if (is.null(txnames)) txnames <- paste("Treatment", result$network$Treat.order)

  if(result$network$response != "multinomial"){
    if(ntreat ==2){
      x <- rank.table[-ntreat,]
    } else{
      x <- apply(apply(rank.table[-ntreat,],2,cumsum),2,sum)/(ntreat-1)
    }
    names(x) <- txnames
  } else if(result$network$response == "multinomial"){
    ncat <- dim(rank.table)[3]
    if (is.null(catnames)) catnames <- paste("Outcome Category with base 1 and comparison", 1+seq(ncat))
    x <- array(NA,dim(rank.table)[2:3])
    for (i in seq(ncat)){
      if(ntreat ==2){
        x[,i] <-   rank.table[-ntreat,,i]
      } else{
        x[,i] <- apply(apply(rank.table[-ntreat,,i],2,cumsum),2,sum)/(ntreat-1)
      }
      dimnames(x) <- list(txnames,catnames)
    }
  }
  return(x)
}

#################### Deviance calculation and plots

#' Find deviance statistics such as DIC and pD.
#'
#' Calculates deviance statistics. This function is automatically called in \code{network.run} and the deviance statistics are stored after sampling is finished.
#'
#' @param result object created by \code{network.run} function
#' @return
#' \item{Dbar}{overall residual deviance}
#' \item{pD}{sum of leverage_arm (i.e. total leverage)}
#' \item{DIC}{deviance information criteria (sum of Dbar and pD)}
#' \item{data.points}{total number of arms in the meta analysis}
#' \item{dev_arm}{posterior mean of the residual deviance in each trial arm}
#' \item{devtilda_arm}{deviance at the posterior mean of the fitted values}
#' \item{leverage_arm}{dev_arm - devtilda_arm for each trial}
#' \item{rtilda_arm}{posterior mean of the fitted value for binomial and multinomial}
#' \item{ybar_arm}{posterior mean of the fitted value for normal}
#' @examples
#' #parkinsons
#' network <- with(parkinsons, {
#'  network.data(Outcomes, Study, Treat, SE = SE, response = "normal")
#' })
#' result <- network.run(network)
#' calculate.deviance(result)
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

calculate.deviance <- function(result){

  if(result$network$dic == FALSE){
    stop("need to set dic = TRUE to be able to use this function")
  }

  network <- result$network
  samples <- result$samples

  totresdev <- lapply(samples, function(x){ x[,"totresdev"]})
  Dbar <- mean(unlist(totresdev))

  ###### find residual deviance by arm
  if(network$response == "multinomial" & !is.null(network$miss.matrix)){
    dev <- list()
    for(ii in seq(network$npattern)){
      dev_each <- lapply(samples, function(x) { x[,grep(paste0("dev", ii, "\\["), dimnames(samples[[1]])[[2]])]})
      dev_each <- do.call(rbind, dev_each)
      dev_each <- apply(dev_each, 2, mean)

      n_value <- network[[paste0("n", ii)]]
      dev_matrix <- matrix(NA, nrow = dim(n_value)[1], ncol = dim(n_value)[2])

      for(i in 1:dim(dev_matrix)[1]){
        for(j in 1:dim(dev_matrix)[2]){
          ind <- which(paste0("dev", ii, "[", i, ",", j, "]") == names(dev_each))
          if(length(ind) != 0){
            dev_matrix[i,j] <- dev_each[ind]
          }
        }
      }
      dev[[paste0("dev", ii)]] <- dev_matrix
    }
    dev_arm <- do.call(rbind, dev)
  } else{
    dev <- lapply(samples, function(x) { x[,grep("dev\\[", dimnames(samples[[1]])[[2]])]})
    dev <- do.call(rbind, dev)
    dev <- apply(dev, 2, mean)

    dev_matrix <- matrix(NA, nrow =  network$nstudy, ncol = max(network$na))
    for(i in 1:dim(dev_matrix)[1]){
      for(j in 1:dim(dev_matrix)[2]){
        ind <- which(paste("dev[", i, ",", j, "]", sep = "") == names(dev))
        if(length(ind) != 0){
          dev_matrix[i,j] <- dev[ind]
        }
      }
    }
    dev_arm <- dev_matrix
  }

  ############find leverage
  if(network$response == "binomial"){
    rtilda <- lapply(samples, function(x){ x[,grep("rhat\\[", dimnames(samples[[1]])[[2]])] })
    rtilda <- do.call(rbind, rtilda)
    rtilda <- apply(rtilda, 2, mean)

    rtilda_arm <- devtilda_arm <- matrix(NA, nrow = network$nstudy, ncol = max(network$na))
    for(i in 1:network$nstudy){
      for(j in 1:network$na[i]){
        r_value <- network$r[i,j]
        n_value <- network$n[i,j]
        rtilda_arm[i,j] <- rtilda[which(paste("rhat[", i, ",", j, "]", sep = "") == names(rtilda))]

        devtilda_arm[i,j] <- ifelse(r_value != 0, 2 * r_value * (log(r_value)-log(rtilda_arm[i,j])), 0)
        devtilda_arm[i,j] <- devtilda_arm[i,j] + ifelse((n_value - r_value) != 0, 2 * (n_value-r_value) *(log(n_value-r_value) - log(n_value- rtilda_arm[i,j])), 0)
      }
    }
  } else if(network$response == "normal"){
    ybar <- lapply(samples, function(x){ x[,grep("theta\\[", dimnames(samples[[1]])[[2]])] })
    ybar <- do.call(rbind, ybar)
    ybar <- apply(ybar, 2, mean)

    ybar_arm <- devtilda_arm <- matrix(NA, nrow = network$nstudy, ncol = max(network$na))

    for(i in 1:network$nstudy){
      for(j in 1:network$na[i]){
        r_value <- network$r[i,j]
        se_value <- network$se[i,j]
        ybar_arm[i,j] <- ybar[which(paste("theta[", i, ",", j, "]", sep = "") == names(ybar))]
        devtilda_arm[i,j] <- ifelse(se_value != 0, (r_value - ybar_arm[i,j])^2 / se_value^2, 0)
      }
    }
  } else if(network$response == "multinomial"){
    if(is.null(network$miss.matrix)){ #complete dataset
      rtilda <- lapply(samples, function(x){ x[,grep("rhat\\[", dimnames(samples[[1]])[[2]])]})
      rtilda <- do.call(rbind, rtilda)
      rtilda <- apply(rtilda, 2, mean)

      rtilda_arm <- devtilda_category <- array(NA, dim = c(network$nstudy, max(network$na), network$ncat))
      for(i in 1:network$nstudy){
        for(j in 1:network$na[i]){
          for(k in 1:network$ncat){
            r_value <- network$r[i,j,k]
            rtilda_arm[i,j,k] <- rtilda[which(paste("rhat[", i, ",", j, ",", k, "]", sep = "") == names(rtilda))]
            devtilda_category[i,j,k] <- ifelse(r_value != 0,  2 * r_value * log(r_value/rtilda_arm[i,j,k]), 0)
          }
        }
      }
      devtilda_arm <- apply(devtilda_category, 1:2, sum)
    } else{ #incomplete datacase
      devtilda_value <- rtilda_arm <- list()
      for(ii in seq(network$npattern)){
        r_values <- network[[paste0("r",ii)]]
        devtilda_category <- rtilda_matrix <- array(NA, dim = dim(r_values))

        rtilda <- lapply(samples, function(x){ x[,grep(paste0("rhat", ii, "\\["), dimnames(samples[[1]])[[2]])] })
        rtilda <- do.call(rbind, rtilda)
        rtilda <- apply(rtilda, 2, mean)

        for(i in 1:dim(r_values)[1]){
          for(j in 1:dim(r_values)[2]){
            for(k in 1:dim(r_values)[3]){
              found <- which(paste("rhat", ii, "[", i, ",", j, ",", k, "]", sep = "") == names(rtilda))
              r_value <- r_values[i,j,k]
              if(!is.na(r_value) & length(found) != 0){
                rtilda_matrix [i,j,k] <- rtilda[found]
                devtilda_category[i,j,k] <- ifelse(r_value != 0,  2 * r_value * log(r_value/rtilda_matrix[i,j,k]), 0)
              }
            }
          }
        }
        devtilda_matrix <- apply(devtilda_category, 1:2, sum)
        rtilda_arm[[ii]] <- rtilda_matrix
        devtilda_value[[ii]] <- devtilda_matrix
      }
      devtilda_arm <- do.call(rbind, devtilda_value)
    }
  }
  leverage_arm <- dev_arm - devtilda_arm
  pD <- sum(leverage_arm, na.rm = TRUE)
  DIC <- Dbar + pD

  out <- list(Dbar = Dbar, pD = pD, DIC = DIC, data.points = sum(network$na), dev_arm = dev_arm, devtilda_arm = devtilda_arm, leverage_arm = leverage_arm)
  if(network$response == "binomial" || network$response == "multinomial"){
    out$rtilda_arm = rtilda_arm
  } else if(network$response == "normal"){
    out$ybar_arm = ybar_arm
  }
  return(out)
}

#' make a deviance plot
#'
#' This makes a deviance plot which plots residual deviance (dev_arm) vs. all the arms for each study.
#' @param result object created by \code{network.run} function
#' @examples
#' network <- with(blocker, {
#'  network.data(Outcomes, Study, Treat, N = N, response = "binomial")
#' })
#' result <- network.run(network)
#' network.deviance.plot(result)
#' @export

network.deviance.plot <- function(result){
  deviance <- result$deviance
  dev_vector <- c(t(deviance$dev_arm))
  dev_vector <- dev_vector[!is.na(dev_vector)]
  plot(seq(sum(result$network$na)), dev_vector, xlab = "Arm", ylab = "Residual Deviance", main = "Per-arm residual deviance")
}

#' make a leverage plot
#'
#' Make a leverage vs. square root of residual deviance plot
#' @param result object created by \code{network.run} function
#' @export

network.leverage.plot <- function(result){
  deviance <- result$deviance
  dev <- sqrt(apply(deviance$dev_arm, 1, mean, na.rm = TRUE))
  leverage <- apply(deviance$leverage, 1, mean, na.rm = TRUE)
  plot(dev, leverage, xlim = c(0, max(c(dev, 2.5))), ylim = c(0, max(c(leverage,4))),
       xlab = "Square root of residual deviance", ylab = "Leverage", main = "Leverage versus residual deviance")
  mtext("Per-study mean per-datapoint contribution")
}

#' Make a covariate plot
#'
#' Make a covariate plot of how the relative effect changes as the covariate value changes. Plot is created for each one of the covariate.
#' User needs to specify one base treatment and one comparison treatment to make this plot (base category and comparison category is needed for multinomial).
#' It then uses the \code{\link{relative.effects}} to calculate the correct relative effect.
#' 2.5\% and 97.5\% C.I. are drawn along with the median value.
#'
#' @param result object created by \code{network.run} function
#' @param base.treatment base treatment for relative effect
#' @param comparison.treatment treatment comparing against base treatment
#' @param base.category base category for multinomial data
#' @param comparison.category comparison category for multinomial data
#' @param covariate.name A vector of covariate names naming of the covariate that goes into x axis label
#' @examples
#' ########### certolizumab (with covariate)
#' network <- with(certolizumab, {
#'  network.data(Outcomes, Study, Treat, N=N, response="binomial", Treat.order,
#'  covariate = covariate, hy.prior = list("dhnorm", 0, 9.77))
#' })
#' result <- network.run(network)
#' network.covariate.plot(result, base.treatment = "Placebo", comparison.treatment = "CZP",
#' covariate.name = "Disease Duration")
#' @export

network.covariate.plot <- function(result, base.treatment = NULL, comparison.treatment= NULL, base.category = NULL, comparison.category = NULL, covariate.name = NULL){

  if(is.null(network$covariate)){
    stop("need to provide covariate information to make this plot")
  }
  if(result$network$response != "multinomial"){
    if(is.null(base.treatment) || is.null(comparison.treatment)){
      stop("need to specify both base.treatment and comparison.treatment")
    }
  } else{
    if(is.null(base.treatment) || is.null(comparison.treatment) || is.null(base.category) || is.null(comparison.category)){
      stop("need to specify all base.treatment, comparison.treatment, base.category, and comparison.category")
    }
  }

  network <- result$network
  observed <- network$covariate
  xvals <- matrix(NA, nrow = dim(network$covariate)[2], ncol = 7)
  xlim <- matrix(NA, nrow = dim(network$covariate)[2], ncol = 2)
  covariate_mx <- NULL
  for(i in 1:dim(network$covariate)[2]){
    xlim[i,] <- c(min(observed[,i], na.rm = TRUE), max(observed[,i], na.rm = TRUE))
    xvals[i,] <- seq(xlim[i,1], xlim[i,2], length.out = 7)
    covariate_mx <- c(covariate_mx, network[[paste0("mx",i)]])
  }

  for(i in 1:dim(network$covariate)[2]){
    res <- lapply(xvals[i,], function(xval) {
      covariate <- covariate_mx
      covariate[i] <- xval
      if(network$response != "multinomial"){
        samples <- relative.effects(result, base.treatment, comparison.treatment, covariate = covariate)
      } else{
        samples <- relative.effects(result, base.treatment, comparison.treatment, base.category, comparison.category, covariate = covariate)
      }

      samples <- as.matrix(samples)
      stats <- t(apply(samples, 2, quantile, probs = c(0.025, 0.5, 0.975)))
      data.frame(median = stats[,"50%"], lower = stats[,"2.5%"], upper = stats[,"97.5%"])
    })
    res <- do.call(rbind,res)

    dim_names <- if(network$response != "multinomial"){
      dimnames(as.matrix(relative.effects(result, base.treatment, comparison.treatment)))[[2]]
    } else{
      dimnames(as.matrix(relative.effects(result, base.treatment, comparison.treatment, base.category, comparison.category)))[[2]]
    }

    ylim <- c(min(res), max(res))
    xlab_name <- ifelse(is.null(covariate.name), paste0("covariate ", i), covariate.name[i])

    plot(xvals[i,], res$median, type = "l", xlim = xlim[i,], ylim = ylim, main = "Treatment effect vs. covariate", xlab = xlab_name, ylab = dim_names)
    lines(xvals[i,], res$lower, lty = 2)
    lines(xvals[i,], res$upper, lty = 2)
  }
}

#' Calculates correlation matrix for multinomial heterogeneity parameter.
#'
#' Calculates correlation matrix from the variance matrix for heterogeneity parameter. Only used for multinomial.
#' @param result object created by \code{network.run} function
#' @examples
#' #cardiovascular
#' network <- with(cardiovascular, {
#'  network.data(Outcomes, Study, Treat, N, response = "multinomial")
#' })
#' result <- network.run(network)
#' variance.tx.effects(result)
#' @export

variance.tx.effects = function(result)
{
  if(result$network$response != "multinomial"){
    stop("this function is used only for multinomial response")
  }
  samples_sigma <- pick.summary.variables(result, only.pars = c("sigma"))
  samples_sigma <- do.call(rbind, samples_sigma)
  samples_sigma <- apply(samples_sigma, 2, mean)

  sigma_matrix <- matrix(samples_sigma, nrow = result$network$ncat-1)
  cor_matrix <- sigma_matrix/outer(sqrt(diag(sigma_matrix)),sqrt(diag(sigma_matrix)))

  return(list(sigma_matrix = sigma_matrix, cor_matrix = cor_matrix))
}


#' Draws forest plot
#'
#' Draws forest plot of pooled treatment effect. Reports odds ratio for binomial and multinomial outcomes and continuous scale for normal outcomes.
#' @param result object created by \code{network.run} function
#' @param level confidence interval level (default is 0.95)
#' @param xlim horizontal limits of the plot region
#' @param alim x-axis limit on the forest plot
#' @param ylim vertical limits of the plot region. 
#' @param at position of the x-axis tick marks. If left unspecified, the function tries to set it at sensible values
#' @param steps the number of tick marks for the x-axis (the default is 5). Ignored when the user specifies the positions via the at argument.
#' @param refline value at which a vertical 'refrence' line should be drawn (the default is 0). The line can be suppressed by setting this argument to NA
#' @param treat_lab optional vector with labels for the k studies. If unspecified, simple labels are created within the function. To suppress labels, set this argument to NA.
#' @references W. Viechtbauer (2010), \emph{Conducting meta-analyses in R with the metafor package}, Journal of Statistical Software, 36(3):1-48. [\url{https://doi.org/10.18637/jss.v036.i03}]
#' @export

network.forest.plot2 <- function(result, level = 0.95){
  
  mean_store <- relative.effects.table(result)
  y <- mean_store[result$network$Treat.order[1],]
  y <- y[!is.na(y)]
  
  if(result$network$response %in% c("binomial", "multinomial")){
    y <- exp(y)
  } 
  
  Treat.order <- result$network$Treat.order
  
  ts <- 1:length(Treat.order)
  comps <- combn(ts, 2)
  
  # if(result$network$response != "multinomial"){
  #   tbl <- matrix(NA, nrow = length(ts), ncol = length(ts), dimnames = list(Treat.order, Treat.order))
  #   
  #   for (i in 1:ncol(comps)) {
  #     comp <- comps[, i]
  #     samples <- as.matrix(relative.effects(result, base.treatment = Treat.order[comp[1]], comparison.treatments = Treat.order[comp[2]]))
  #     
  
  
  odds_ratio <- matrix(NA, nrow = length(result.list), ncol = 3)
  
  for(i in 1:length(result.list)){
    result <- result.list[[i]]
    samples <- do.call(rbind, result$samples)
    odds_ratio[i,] <- exp(quantile(samples[,grep("beta", colnames(samples))], c((1 -level)/2, 0.5, 1 - (1 -level)/2)))
  }
  
  odds <- as.data.frame(odds_ratio)
  names(odds) <- c("lower", "OR", "upper")
  
  if(is.null(result.name)){
    odds$vars <- row.names(odds)  
  } else{
    if(length(result.name) != length(result.list)){
      stop("result.name should have same length as result.list")
    }
    odds$vars <- result.name
  }
  ticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
  ggplot(odds, aes(y = OR, x = factor(vars))) + 
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
    scale_y_log10(breaks = ticks, labels = ticks) +
    geom_hline(yintercept = 1, linetype = 2) +
    coord_flip() +
    labs(x = "Variables", y = "Odds Ratio", title = title) +
    theme_bw()  
}




#' Draws forest plot
#' 
#' Draws forest plot of pooled treatment effect. Reports odds ratio for binomial and multinomial outcomes and continuous scale for normal outcomes.
#' @param result object created by \code{network.run} function
#' @level numerical value between 0 and 1 specifying the confidence level percentage (default is at 0.95)
#' @xlim horizontal limits of the plot region
#' @alim x-axis limit on the forest plot
#' @ylim vertical limits of the plot region. 
#' @at position of the x-axis tick marks. If left unspecified, the function tries to set it at sensible values
#' @steps the number of tick marks for the x-axis (the default is 5). Ignored when the user specifies the positions via the at argument.
#' @refline value at which a vertical 'refrence' line should be drawn (the default is 0). The line can be suppressed by setting this argument to NA
#' @treat_lab optional vector with labels for the k studies. If unspecified, simple labels are created within the function. To suppress labels, set this argument to NA.
#' @references W. Viechtbauer (2010), \emph{Conducting meta-analyses in R with the metafor package}, Journal of Statistical Software, 36(3):1-48. [\url{https://doi.org/10.18637/jss.v036.i03}]
#' @export

network.forest.plot <- function(result, level = 0.95, xlim = NULL, alim = NULL, ylim = NULL, at = NULL, steps = 5, refline = 0, treat_lab = NULL){
  
  mean_store <- relative.effects.table(result)
  y <- mean_store[result$network$Treat.order[1],]
  y <- y[!is.na(y)]
  
  se_store <- relative.effects.table(result, summary_stat = "sd")
  se <- se_store[result$network$Treat.order[1],]
  se <- se[!is.na(se)]
  variance <- se^2
    
  if(result$network$response %in% c("binomial", "multinomial")){
    y <- exp(y)
    se <- exp(se)
    variance <- exp(variance)
  } 
  
  ci.lb <- y - qnorm((1 - level)/2, lower.tail = FALSE) * se
  ci.ub <- y + qnorm((1 - level)/2, lower.tail = FALSE) * se
  
  if(is.null(treat_lab)){
    treat_lab <- paste("Treatment", names(y))
  }
  
  k <- length(y)
  
  if(k != length(treat_lab)){
    stop("Number of treatment label does not match the number of outcomes")
  }
  
  rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
  if(is.null(xlim)){
    xlim <- c(min(y, na.rm = TRUE) - rng, max(y, na.rm = TRUE) + rng)
  }
  
  if(is.null(ylim)){
    ylim <- c(0.5, k + 3)
  }
  
  if(is.null(alim)){
    alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, na.rm = TRUE)), n = steps - 1))
  }
  
  if(is.null(at)){
    at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub, na.rm = TRUE)), n = steps - 1)
  }
  
  # par.mar <- par("mar")
  # par.mar.adj <- par.mar - c(0,3,1,1)
  # par.mar.adj[par.mar.adj < 0] <- 0
  # par(mar = par.mar.adj)
  # on.exit(par(mar = par.mar))
    
  plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", yaxt = "n", xaxt = "n", xaxs = "i", bty = "n")
  abline(h = ylim[2] -2)
  
  if(is.numeric(refline)){
    segments(refline, ylim[1] - 5, refline, ylim[2] - 2, lty = "dotted", col = "black")
  }
  
  axis(side = 1, at = at)
  
  if(is.null(xlab)){
    if(result$network$response %in% c("binomial", "multinomial")){
      xlab <- "Odds Ratio"
    } else{
      xlab <- "Observed Outcome"
    }
  }
  mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2, line = par("mgp")[1] - 0.5)
  
  weight <- 1/sqrt(variance)
  psize <- weight/ sum(weight, na.rm = TRUE)
  psize <- (psize - min(psize, na.rm = TRUE)) / (max(psize, na.rm = TRUE) - min(psize, na.rm = TRUE))
  psize <- psize + 0.5
  
  rows <- k:1
  for(i in 1:k){
    if(!is.na(y[i])){
      segments(ci.lb[i], rows[i], ci.ub[i], rows[i])
      points(y[i], rows[i], cex = psize[i], pch = 0)
    }
  }
  text(xlim[1], rows, treat_lab, pos = 4)
  
  annotext <- cbind(y, ci.lb, ci.ub)
  annotext <- formatC(annotext, format = "f", digits = 2)
  annotext <- cbind(annotext[,1], " [", annotext[,2], ", ", annotext[,3], "]")
  annotext <- apply(annotext, 1, paste, collapse = "")
  text(x = xlim[2], rows, labels = annotext, pos = 2)
  
  text(x = alim[1], ylim[2]-1, pos = 4, "Comparison: other vs 'placebo'", cex = 1.5)
}
