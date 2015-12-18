#' Josh's Follicle Simulation
#' 
#' @param NF Starting number of follicles.
#' @param ND Number of 'days' or steps of the simulation.
#' @param IGP Initial growing pool of follicles.
#' @param phold Probability of 'holding' at the current number of cells.
#' @param cond.pdub Conditional probability of 'doubling', given that we
#' didn't 'hold' initially.
#' @param pcelllive Josh's secret sauce of the probability of a cell
#' surviving a doubling step, because they sometimes don't make it.
#' @param ejectnum After reaching this number of cells, a healthy follicle
#' is ejected?
#' @param cyclength The length of the ovulatory cycle in days.
#' @param puberty Whether an initial 750 should have a head start.
#' @param verbose Print the progress along the way?
#' @param pdfname If a PDF file name is given, a plot will be produced.
josh <- function(NF = 3000,
                 ND = 420,
                 IGP = 300,
                 phold = 0.995,
                 cond.pdub = 0.9,
                 pcelllive = 0.8,
                 ejectnum = 50000,
                 cyclength = 4,
                 puberty = TRUE,
                 verbose = TRUE,
                 pdfname = NA)
{
                 
  # Initialize, randomly drawing 1, 2 or 3 cells
  x <- sample(1:3, NF, replace=TRUE)

  pdub <- (1-phold) * cond.pdub        # Unconditional doubling probability
  pdie <- (1-phold) * (1-cond.pdub)    # Unconditional death probability

  # These are the unconditional transition probabilities at the start:
  # JAY: Why 300 and 19?  2^19 > 500000
  # Dec 1: Changing 300 to 750 and reducing 19 to 16 because
  # we think 2^16 is what gets us over the new 50,000 ejection number.
  probs <- matrix(c(phold, pdub, pdie), NF, 3, byrow=TRUE)
  if (puberty) {
    x[1:IGP] <- 2^sample(1:16, IGP, replace=TRUE,
                         dexp(1:16, rate=0.3))
    probs[1:IGP,] <- matrix(c(0, cond.pdub, 1-cond.pdub),
                            nrow=IGP, ncol=3, byrow=TRUE)
  }

  z <- matrix(0, NF, ND)
  for (i in 1:ND) {
    if (verbose && (i/10 == floor(i/10))) cat("Step", i, "\n")
    for (j in 1:NF) {
      res <- sample(1:3, 1, prob=probs[j,]) # Status result
      
      # If "double" then we never "hold" in future steps and the probs change:
      if (res==2) { probs[j,] <- c(0, cond.pdub, 1-cond.pdub) }
      
      # If not "hold" then either double (2) or die (3)
      
      if (res > 1) {
        x[j] <- x[j] * 2 * (res==2)
        if (x[j] > 0) x[j] <- rbinom(1, x[j], pcelllive)
      }

      if (x[j] >= ejectnum) { probs[j,] <- c(1, 0, 0) } # Fol. ejected, stop
    } # End of 1:NF loop
    z[,i] <- x
  } # End of 1:ND loop

  z[z>ejectnum] <- ejectnum
  
  if (verbose) {
    cat("Dead follicles:\t\t", sum(z[,ncol(z)] == 0), "\n", sep="\t")
    cat("Ejected eggs:\t", sum(z[,ncol(z)] == ejectnum), "\n", sep="\t")
    cat("Follicles underdeveloped:", sum(z[,ncol(z)]>0 & z[,ncol(z)]<ejectnum),
        "\n", sep="\t")
  }

  if (!is.na(pdfname)) {
    pdf(pdfname, width=7.5, height=10, )
      par(mfrow=c(3,1))
      par(mar=c(5,7,4,2), ylbias=0.95)
    
      # Last one that we want ********************************
      pfn <- read.csv("extrapoints.csv", as.is=TRUE)
      b <- read.csv("confsim.csv", as.is=TRUE)
      b <- apply(b, 2, quantile, probs=c(0.01, 0.99))
      foo <- apply(z==1 | z==2 | z==3, 2, sum)
      plot(1:ND, foo, xlab="Time (Days)",
           ylab="Number of Primordial Follicles (1st and 99th Percentiles)",
           pch="", ylim=range(pretty(c(0, foo))),
           main="Decline of Primordial Follicle Pool")
      polygon(c(1:ND, ND:1), c(b[1,], b[2,ND:1]), col="grey",
              density=60)
      points(pfn$days, pfn$pfn)
      wantfoo <- foo
      
      if (FALSE) {
      # The ones that make it, leading to an egg ejected from the ovary:
      plot(1:ND, seq(1, 500000, length.out=ND), 
           pch="", main=paste("Cell Counts of Ejected Eggs (", sum(z[,ncol(z)] == ejectnum),
                              ")", sep=""),
           xlab="Day", ylab="",
           las=1, yaxt="n")
      axis(2, at=c(1, 1e5, 2e5, 3e5, 4e5, 5e5),
           labels=format(c(1, 1e5, 2e5, 3e5, 4e5, 5e5), scientific=FALSE), las=1)
      these <- which(apply(z==ejectnum, 1, any))
      for (i in these) {
        end <- min(which(z[i,]==ejectnum))
        lines(1:end, z[i,1:end], lty=3)
      }
      }
      
      # The ones that start to develop but die:     ****************
      plot(1:ND, seq(1, 50000, length.out=ND), 
           pch="", main="Growth of Follicles That Die Via Atresia",
           xlab="Time (Days)", ylab="Number of Granulosa Cells",
           las=1, yaxt="n")
      #main=paste("Cell Counts of Dead Follicles (", sum(z[,ncol(z)] == 0),
      #           ")", sep=""),
      axis(2, at=c(1, 1e4, 2e4, 3e4, 4e4, 5e4),
           labels=format(c(1, 1e4, 2e4, 3e4, 4e4, 5e4), scientific=FALSE), las=1)
      these <- which(apply(z==0, 1, any))
      for (i in these) {
        end <- max(1, min(which(z[i,]==0)) - 1)
        lines(1:end, z[i,1:end], lty=3)
        text(which.max(z[i,]), max(z[i,]), "D", cex=0.5)
      }
      
      if (FALSE) {
      # Partly developed or "holds":
      # The ones that start to develop but die:
      plot(1:ND, seq(1, 500000, length.out=ND), 
           pch="", main=paste("Cell Counts of Incomplete Follicles (",
                              sum(z[,ncol(z)]>0 & z[,ncol(z)]<ejectnum),
                              ")", sep=""),
           xlab="Day", ylab="",
           las=1, yaxt="n")
      axis(2, at=c(1, 1e5, 2e5, 3e5, 4e5, 5e5),
           labels=format(c(1, 1e5, 2e5, 3e5, 4e5, 5e5), scientific=FALSE), las=1)
      for (i in 1:nrow(z)) {
        if (z[i,ncol(z)] > 0 && z[i,ncol(z)] < ejectnum) {
          lines(1:ND, z[i,], lty=3)
          text(ND+1, max(z[i,]), "U") 
        }
      }
      
      foo <- apply(!(z==1 | z==2 | z==3), 1, function(a) min(c(420, which(a))) )
      hist(foo, main="Distribution of Times of First Growth",
           xlab="Time of First Growth")
      }
        
      
      # More breaks and add the lines ****************************
      foo <- apply(z==ejectnum, 1,
                   function(a) {
                     if (any(a)) return(min(which(a)))
                     else return(NA)
                   })
      hist(foo, main="Distribution of Ovulatory Follicle Development",
           xlab="Time (Days)", 
           ylab=paste("Number of Ovulatory Follicles per", cyclength,
                      "Day Estrus Cycle"),
           breaks=seq(0, 420, by=cyclength))
      abline(v=seq(1, 420, by=28), col="black", lty=2, lwd=1)
      a <- table(cut(foo, breaks=seq(0, 420, by=cyclength)))
      aa <- data.frame(a)
      print(aa)    
      
    dev.off()
  }
  
  return( list(all=z, dead=sum(z[,ncol(z)] == 0),
               underdeveloped=sum(z[,ncol(z)]>0 & z[,ncol(z)]<500000),
               ejected=sum(z[,ncol(z)] == 0),
               wantfoo=wantfoo) )
}

# Want to run 1000 times and get the info for the first
# plot, and collect $wantfoo for "error bands" around
# the curve.
NNN <- 1
foo <- matrix(0, NNN, 420)
for (i in 1:NNN) {
  set.seed(i)
  ans <- josh(pdfname=paste("test", i, ".pdf", sep=""))
  foo[i,] <- ans$wantfoo
}
if (NNN>1) write.table(foo, "confsim.csv", sep=",", row.names=F, col.names=T)


write.csv(ans$all, "sample.csv", row.names=F)

# another <- josh(phold = 0.99,
#                 cond.pdub = 0.93,
#                 pcelllive = 0.85,
#                 pdfname="test2.pdf")

# Primordial, Primary, Secondary, Pre-Antrial, Egg
# 1-3, 4-70, 71-5000, 5001-499999, 500k+

# Spiffy up plot titles and file name to reflect inputs?


