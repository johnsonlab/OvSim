#' OvSim
#' 
#' @param NF Starting number of follicles.
#' @param ND Number of 'days' or steps of the simulation.
#' @param IGP Initial growing pool of follicles.
#' @param phold Probability of 'holding' at the current number of cells.
#' @param cond.pdub Conditional probability of 'doubling', given that we
#' didn't 'hold' initially.
#' @param pcelllive The probability of a cell
#' surviving a doubling step, because they sometimes don't make it.
#' @param ejectnum After reaching this number of cells, a healthy follicle
#' is ejected?
#' @param cyclength The length of the ovulatory cycle in days.
#' @param puberty Whether an initial 750 should have a head start.
#' @param verbose Print the progress along the way?
#' @param pdfname If a PDF file name is given, a plot will be produced.
hovsim <- function(NF = 50000,
                   Y = 0.5,
                  ND = 12775,
                  IGP = 300,
                  phold = 0.9995,
                  cond.pdub = 0.9,
                  pcelllive = 0.75,
                  ejectnum = 500000,
                  cyclength = 30,
                  puberty = FALSE,
                  verbose = TRUE,
                  pdfname = NA)
{
  ## NOTE
  ##
  ## This should eventually be made for a formal package.  Currently,
  ## just as an example, some plots are produced within this function.
  ## That's unnecessary and was just quick-and-dirty for our convenience.
  ## However, at the moment this serves to document the result of our
  ## paper. - Jay Emerson
  
  # Initialize, randomly drawing 1, 2 or 3 cells
  x <- sample(1:3, NF, replace=TRUE)
  
  pdub <- (1-phold) * cond.pdub        # Unconditional doubling probability
  pdie <- (1-phold) * (1-cond.pdub)    # Unconditional death probability
  
  # These are the unconditional transition probabilities at the start.
  # JJJ: 2^16 is what gets us over the new 50,000 ejection number.
  probs <- matrix(c(phold, pdub, pdie), NF, 3, byrow=TRUE)
  JJJ <- 16
  if (puberty) {
    x[1:IGP] <- 2^sample(1:JJJ, IGP, replace=TRUE,
                         dexp(1:JJJ, rate=0.3))
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
  
  wantfoo <- NULL
  if (!is.na(pdfname)) {
    pdf(pdfname, width=7.5, height=10, )
    par(mfrow=c(3,1))
    par(mar=c(5,7,4,2), ylbias=0.95)
    
    # For the paper, relies on extrapoints.csv and confsim.csv
    # being available in the current directory
    # ********************************
    #pfn <- read.csv("extrapoints.csv", as.is=TRUE)
    b <- read.csv("confsim.csv", as.is=TRUE)
    b <- apply(b, 2, quantile, probs=c(0.01, 0.99))
    foo <- apply(z==1 | z==2 | z==3, 2, sum)
    plot(1:ND, foo, xlab="Time (Days)",
         ylab="Number of Primordial Follicles (1st and 99th Percentiles)",
         pch="", ylim=range(pretty(c(0, foo))),
         main="Decline of Primordial Follicle Pool")
    polygon(c(1:ND, ND:1), c(b[1,], b[2,ND:1]), col="grey",
            density=60)
    #points(pfn$days, pfn$pfn)
    wantfoo <- foo
    
    # The ones that start to develop but die:     ****************
    plot(1:ND, seq(1, 50000, length.out=ND), 
         pch="", main="Growth of Follicles That Die Via Atresia",
         xlab="Time (Days)", ylab="Number of Granulosa Cells",
         las=1, yaxt="n")
   # axis(2, at=c(1, 1e4, 2e4, 3e4, 4e4, 5e4),
     #    labels=format(c(1, 1e4, 2e4, 3e4, 4e4, 5e4), 
         axis(2, at=c(1, 1e5, 2e5, 3e5, 4e5, 5e5),
         labels=format(c(1, 1e5, 2e5, 3e5, 4e5, 5e5), 
         scientific=FALSE), las=1)
    these <- which(apply(z==0, 1, any))
    for (i in these) {
      end <- max(1, min(which(z[i,]==0)) - 1)
      lines(1:end, z[i,1:end], lty=3)
      text(which.max(z[i,]), max(z[i,]), "D", cex=0.5)
    }
    
    # Last one ****************************
    foo <- apply(z==ejectnum, 1,
                 function(a) {
                   if (any(a)) return(min(which(a)))
                   else return(NA)
                 })
    hist(foo, main="Distribution of Ovulatory Follicle Development",
         xlab="Time (Days)", 
         ylab=paste("Number of Ovulatory Follicles per", cyclength,
                    "Day Menstrual Cycle"),
         breaks=seq(0, 12775, by=cyclength))                    
         #breaks=seq(0, 365, by=cyclength))
    abline(v=seq(1, 12775, by=30), col="black", lty=2, lwd=1)
    #abline(v=seq(1, 420, by=30), col="black", lty=2, lwd=1)
    a <- table(cut(foo, breaks=seq(0, 12775, by=cyclength)))
    #a <- table(cut(foo, breaks=seq(0, 420, by=cyclength)))
    aa <- data.frame(a)
    print(aa)    
    
    dev.off()
  }
  
  return( list(all=z, dead=sum(z[,ncol(z)] == 0),
               underdeveloped=sum(z[,ncol(z)]>0 & z[,ncol(z)]<ejectnum),
               ejected=sum(z[,ncol(z)] == 0),
               wantfoo=wantfoo) )
}
