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

contrast.network.data <- function(Outcomes, Treat, SE, na, V = NULL, type = "random", rank.preference = "higher"
                                  
                                  ){

  if(missing(Outcomes) || missing(Treat) || missing(SE) || missing(na)){
    stop("OUtcomes, Treat, SE, and na have to be all specified")
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
  
  network <- list(Outcomes = Outcomes, Treat = Treat, SE = SE, na = na)
  
  if(is.null(V)){
    network$V <- V
  }
  
  class(network) <- "contrast.network.data"
  return(network)

}

