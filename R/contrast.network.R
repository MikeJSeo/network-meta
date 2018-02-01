#' Make a network object for contrast-level data containing data, priors, and a JAGS model file
#' 
#' This is very similar function as \code{\link{network.data}}. Except it uses contrast-level data instead of arms-level data. Contrast-level format is treatment differences relative to the control arm.
#' Note that in two arm trials there is only one contrast value per trial, but in three arm trials there are two contrast values relative to the control arm.
#'
#' @param Outcomes Contrast-level outcomes. Outcome is treated as normally distributed. 
#' @param Study A vector of study indicator for each contrasts. If there are only two arm trials, this will be a sequence from 1 to total number of studies.
#' @param Treat A vector of treatment indicator for each arm
#' @param SE A vector of standard error for each contrasts.
#' @param V Needed if you have multi-arm trials. Length of this vector should be number of studies. If the study is multi-arm trial, need to specify variance of the baseline treatment in that trial. Denote it with NA if the study only has two-arm trials.
#' @references A.J. Franchini, S. Dias, A.E. Ades, J.P. Jansen, N.J. Welton (2012), \emph{Accounting for correlation in network meta-analysis with multi-arm trials}, Research Synthesis Methods 3(2):142-160. [\url{https://doi.org/10.1002/jrsm.1049}] 
#' @references S. Dias, A.J. Sutton, A.E. Ades, and N.J. Welton (2013a), \emph{A Generalized Linear Modeling Framework for Pairwise and Network Meta-analysis of Randomized Controlled Trials}, Medical Decision Making 33(5):607-617. [\url{https://doi.org/10.1177/0272989X12458724}]
#' @export

contrast.network.data <- function(Outcomes, Study, Treat, SE, V){

  if(missing(SE) || missing(Study) || missing(Treat) || missing(Outcomes)){
    stop("Study, Treat, Outcomes, and SE have to be all specified")
  }
  
  na <- rle(Study)$lengths
  if(any(na == 1)) stop("study cannot have only 1 arm or arms have to be next to each other in each study")
  Study <- rep(1:length(unique(Study)), times = na)
  Treat <- relabel.vec(Treat, Treat.order)
  
  
  Outcomes
  
  na <- rle(Study)$lengths
  if(any(na == 1)) stop("study cannot have only 1 arm or arms have to be next to each other in each study")
  Study <- rep(1:length(unique(Study)), times = na)
  Treat <- relabel.vec(Treat, Treat.order)
  
  
}