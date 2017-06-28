### set prior to default values if not specified

network.prior.default <- function(network, mean.d, prec.d, mean.Eta, prec.Eta, mean.bl, prec.bl, hy.prior.base, mean.cov, prec.cov, hy.prior.cov, hy.prior) {

  response <- network$response
  ncat <- network$ncat
  baseline <- network$baseline
  covariate <- network$covariate
  covariate.model <- network$covariate.model

  if(response == "binomial" || response == "normal"){
    if(is.null(mean.d)){
      mean.d = 0
    }
    if(is.null(prec.d)){
      prec.d = 0.0001
    }
    if(is.null(mean.Eta)){
      mean.Eta = 0
    }
    if(is.null(prec.Eta)){
      prec.Eta = 0.0001
    }
    hy.prior.1 <- hy.prior[[2]]
    hy.prior.2 <- hy.prior[[3]]

    prior.data = list(mean.d = mean.d, prec.d = prec.d, mean.Eta = mean.Eta, prec.Eta = prec.Eta, hy.prior.1 = hy.prior.1, hy.prior.2 = hy.prior.2)

    if(baseline != "none"){
      if(is.null(mean.bl)){
        mean.bl = 0
      }
      if(is.null(prec.bl)){
        prec.bl = 0.0001
      }
      prior.data$mean.bl = mean.bl
      prior.data$prec.bl = prec.bl
      if(baseline == "exchangeable"){
        prior.data$hy.prior.base.1 <- hy.prior.base[[2]]
        prior.data$hy.prior.base.2 <- hy.prior.base[[3]]
      }
    }

    if(!is.null(covariate)){
      if(is.null(mean.cov)){
        mean.cov = 0
      }
      if(is.null(prec.cov)){
        prec.cov = 0.0001
      }
      prior.data$mean.cov = mean.cov
      prior.data$prec.cov = prec.cov

      if(covariate.model == "exchangeable"){
        prior.data$hy.prior.cov.1 <- hy.prior.cov[[2]]
        prior.data$hy.prior.cov.2 <- hy.prior.cov[[3]]
      }
    }
    return(prior.data)
  } else if(response == "multinomial"){
    if(is.null(mean.d)){
      mean.d = rep(0, ncat - 1)
    }
    if(is.null(prec.d)){
      prec.d = diag(0.0001, ncat - 1)
    }
    if(is.null(mean.Eta)){
      mean.Eta = rep(0, ncat - 1)
    }
    if(is.null(prec.Eta)){
      prec.Eta = diag(0.25, ncat - 1)
    }
    hy.prior.1 <- hy.prior[[2]]
    hy.prior.2 <- hy.prior[[3]]
    prior.data = list(mean.d = mean.d, prec.d = prec.d, mean.Eta = mean.Eta, prec.Eta = prec.Eta, hy.prior.1 = hy.prior.1, hy.prior.2 = hy.prior.2)

    if(network$baseline != "none"){
      if(is.null(mean.bl)){
        mean.bl = rep(0, network$ncat - 1)
      }
      if(is.null(prec.bl)){
        prec.bl = diag(0.0001, network$ncat - 1)
      }
      prior.data$mean.bl = mean.bl
      prior.data$prec.bl = prec.bl
      if(baseline == "exchangeable"){
        prior.data$hy.prior.base.1 <- hy.prior.base[[2]]
        prior.data$hy.prior.base.2 <- hy.prior.base[[3]]
      }
    }
    if(!is.null(network$covariate)){
      if(is.null(mean.cov)){
        mean.cov = rep(0, network$ncat - 1)
      }
      if(is.null(prec.cov)){
        prec.cov = diag(0.0001, network$ncat - 1)
      }
      prior.data$mean.cov = mean.cov
      prior.data$prec.cov = prec.cov

      if(covariate.model == "exchangeable"){
        prior.data$hy.prior.cov.1 <- hy.prior.cov[[2]]
        prior.data$hy.prior.cov.2 <- hy.prior.cov[[3]]
      }
    }
  }
  return(prior.data)
}

hy.prior.default <- function(network){
  response <- network$response
  ncat <- network$ncat

  hy.prior <- if(response == "binomial"){
    list("dunif", 0, 5)
  } else if(response == "normal"){
    list("dunif", 0, 100)
  } else if(response == "multinomial"){
    list("dwish", diag(ncat - 1), ncat -1)
  }
  return(hy.prior)
}
