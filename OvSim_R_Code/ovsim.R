#' OvSim
#' 
#' @param NF Starting number of follicles.
#' @param ND Number of 'days' or steps of the simulation.
#' @param IGP Initial growing pool of follicles.
#' @param phold Probability of 'holding' at the current number of cells.
#' This is a fixed probability in the basic simulation, but can be customized.
#' As an example, if `phold` is the character string 'custom1', the code
#' shows a simple decline in `phold` from 0.995 down to 0.99 depending on
#' the number of growth-activiated follicles at a given point in time being
#' less than 100 (i.e. in older stages of the lifecycle). 
#' Following this example, the `OvSim` function may be modified easily for
#' other proposed variations of `phold`.
#' @param cond.pdub Conditional probability of 'doubling', given that we
#' didn't 'hold' initially.
#' @param pcelllive The probability of a cell
#' surviving a doubling step, because they sometimes don't make it.
#' @param ejectnum After reaching this number of cells, a healthy follicle
#' is ejected?
#' @param cyclength The length of the ovulatory cycle in days.
#' @param puberty Whether an initial `IGP` (300 by defualt)
#' should have a head start.
#' @examples
#' ans <- ovsim()
#' ans2 <- ovsim(phold="custom1")
#' @export
ovsim <- function(NF = 3000,
                  ND = 420,
                  IGP = 300,
                  phold = 0.995,
                  cond.pdub = 0.9,
                  pcelllive = 0.8,
                  ejectnum = 50000,
                  cyclength = 4,
                  puberty = TRUE)
{
  ## NOTE
  ##
  ## This should eventually be made for a formal package.  Currently,
  ## just as an example, some plots are produced within this function.
  ## That's unnecessary and was just quick-and-dirty for our convenience.
  ## However, at the moment this serves to document the result of our
  ## paper. - Jay Emerson
  
  # Initialize, randomly drawing 1, 2 or 3 cells (uniformly for now)
  x <- sample(1:3, NF, replace=TRUE)
  
  ### Want phold to possibly fall over time because fewer cells are in
  ### the growth state (and producing AMH in this growth state).  An
  ### example is added below.
  
  phold.orig <- phold
  if (phold.orig=="custom1") phold <- 0.995
  phold.new <- phold
  pdub <- (1-phold) * cond.pdub        # Unconditional doubling probability
  pdie <- (1-phold) * (1-cond.pdub)    # Unconditional death probability

  # These are the unconditional transition probabilities at the start.
  # JJJ: 2^16 is what gets us over the current 50,000 ejection number,
  # for example.
  probs <- matrix(c(phold, pdub, pdie), NF, 3, byrow=TRUE)
  JJJ <- min(which(2^(1:25) > ejectnum))
  if (puberty) {
    x[1:IGP] <- 2^sample(1:JJJ, IGP, replace=TRUE,
                         dexp(1:JJJ, rate=0.3))
    probs[1:IGP,] <- matrix(c(0, cond.pdub, 1-cond.pdub),
                            nrow=IGP, ncol=3, byrow=TRUE)
  }
  
  z <- matrix(0, NF, ND)
  nga.byday <- rep(0, ND)
  for (i in 1:ND) {
    
    if (i > 1) {
      if (phold.orig=="custom1") {
        # Here calculate phold based on the number growth-activated
        # at i-1: sum(z[,i-1]>0 & probs[,i-1]==0)
        val <- 0.99    # 0.995 for constant, 0.99 for adjusting based
        # on the 0.995 upper value in this example customization
        phold.new <- val + (0.995 - val) *
          min(100, sum(z[,i-1]>0 & probs[,1]==0)) / 100
      }
      # The following only change if phold.new changed, above, but
      # are needed anyway so I just recalculate here:
      pdub.new <- (1-phold.new) * cond.pdub        # Unconditional doubling prob
      pdie.new <- (1-phold.new) * (1-cond.pdub)    # Unconditional death prob
    }
    for (j in 1:NF) {
      if (i>1 && probs[j,1]!=0 && probs[j,1]!=1) {
        probs[j,] <- c(phold.new, pdub.new, pdie.new) 
      }
      
      res <- sample(1:3, 1, prob=probs[j,]) # Status result
      
      # If "double" then we never "hold" in future steps and the probs change:
      if (res==2) { probs[j,] <- c(0, cond.pdub, 1-cond.pdub) }
      
      # If not "hold" then either double (2) or die (3)
      if (res > 1) {
        x[j] <- x[j] * 2 * (res==2)   # i.e. res==3 is die trigger here
        if (x[j] > 0) x[j] <- rbinom(1, x[j], pcelllive)
      }
      
      if (x[j] > ejectnum) {
        probs[j,] <- c(1, 0, 0)
        x[j] <- ejectnum
      } # Fol. ejected, stop
      
    } # End of 1:NF loop (j)
    z[,i] <- x
    nga.byday[i] <- sum(probs[,1]!=0 &
                          (z[,i]==1 | z[,i]==2 | z[,i]==3))
  } # End of 1:ND loop (i)
  
  z[z>ejectnum] <- ejectnum  # Probably not necessary, but... can't hurt
  
  return( list(all=z,
               dead.ND=sum(z[,ncol(z)] == 0),
               nongrowing.ND=sum(z[,ncol(z)]>0 & z[,ncol(z)]<ejectnum & probs[,1]!=0),
               growing.ND=sum(z[,ncol(z)]>0 & probs[,1]==0),
               ejected.ND=sum(z[,ncol(z)] == ejectnum),
               nongrowthactivated.byday=nga.byday,
               NF=NF, ND=ND, IGP=IGP, phold=phold.orig, cond.pdub=cond.pdub,
               pcelllive=pcelllive, puberty=puberty,
               ejectnum=ejectnum, cyclength=cyclength) )
}
