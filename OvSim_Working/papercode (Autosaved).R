
source("ovsim.R")

if (FALSE) {
  
  # Run this twice commenting/uncommenting to produce the two needed
  # CSV files.
  NNN <- 10 #1000
  foo <- matrix(0, NNN, 420)
  for (i in 1:NNN) {
    cat("\n\nStep", i, "\n")
    set.seed(i)
    ans <- ovsim()
    #ans <- ovsim(phold="custom1")
    foo[i,] <- ans$nongrowthactivated.byday
  }
  #if (NNN>1) write.table(foo, "confsim_99.csv", sep=",", row.names=F, col.names=T)
  if (NNN>1) write.table(foo, "confsim_constant.csv", sep=",", row.names=F, col.names=T)
  
}

# Now do a single run and produce the plots.

set.seed(1)
ans <- ovsim()      # Use constant phold for the single run in the paper

pdfname <- "newplots.pdf"
z <- ans$all
ND <- ans$ND
ejectnum <- ans$ejectnum
cyclength <- ans$cyclength

wantfoo <- NULL
if (!is.na(pdfname)) {
  pdf(pdfname, width=7.5, height=10)
  par(mfrow=c(3,1))
  par(mar=c(5,7,4,2), ylbias=0.95)
  
  # For the paper, relies on extrapoints.csv and confsim*.csv
  # being available in the current directory
  # ********************************
  pfn <- read.csv("extrapoints_moredata.csv", as.is=TRUE)
  names(pfn) <- c("days", "pfn")
  bc <- read.csv("confsim_constant.csv", as.is=TRUE)
  bc <- apply(bc, 2, quantile, probs=c(0.01, 0.99))
  b99 <- read.csv("confsim_99.csv", as.is=TRUE)
  b99 <- apply(b99, 2, quantile, probs=c(0.01, 0.99))
  foo <- ans$nongrowthactivated.byday  #apply(z==1 | z==2 | z==3, 2, sum)  # foo --> Non-growth activated XXXXXXXX
  plot(1:ND, foo, xlab="Time (Days)",
       ylab="Number of Primordial Follicles (1st and 99th Percentiles)",
       pch="", ylim=range(pretty(c(0, foo))),
       main="Decline of Primordial Follicle Pool")
  polygon(c(1:ND, ND:1), c(bc[1,], bc[2,ND:1]), col="grey",
          density=60)
  polygon(c(1:ND, ND:1), c(b99[1,], b99[2,ND:1]), col="dimgrey",
          density=60)
  points(pfn$days, pfn$pfn)
  legend("topright", legend=c("Actual Follicle Counts (C57Bl/6)",
                              "Constant Probability of Growth Arrest",
                              "Declining Probability of Growth Arrest"),
         pt.cex=2, fill=c("white", "grey", "dimgrey"), border=c(0,1,1),
         density=c(NA, 60, 60) #, pch=c(1, NA, NA),
  )
  points(291, 2920, cex=2)
  
  wantfoo <- foo                    # Reconsider, not needed now.
  
  # The ones that start to develop but die:     ****************
  plot(1:ND, seq(1, 50000, length.out=ND), 
       pch="", main="Growth of Follicles That Die Via Atresia",
       xlab="Time (Days)", ylab="Number of Granulosa Cells",
       las=1, yaxt="n")
  axis(2, at=c(1, 1e4, 2e4, 3e4, 4e4, 5e4),
       labels=format(c(1, 1e4, 2e4, 3e4, 4e4, 5e4), scientific=FALSE), las=1)
  these <- which(apply(z==0, 1, any))
  for (i in these) {
    end <- max(1, min(which(z[i,]==0)) - 1)
    lines(1:end, z[i,1:end], lty=3)
    text(which.max(z[i,]), max(z[i,]), "D", cex=0.5)
  }
  
  # Last one ****************************  # Change this foo to something...
  foo <- apply(z==ejectnum, 1,
               function(a) {
                 if (any(a)) return(min(which(a)))
                 else return(NA)
               })
  hist(foo, main="Distribution of Ovulatory Follicle Development",
       xlab="Time (Days)", 
       ylab=paste("Number of Ovulatory Follicles per", cyclength,
                  "Day Estrus Cycle"),
       breaks=seq(0, ND, by=cyclength))
  abline(v=seq(1, ND, by=28), col="black", lty=2, lwd=1)
  a <- table(cut(foo, breaks=seq(0, ND, by=cyclength)))
  aa <- data.frame(a)
  print(aa)    
  
  dev.off()
}

  