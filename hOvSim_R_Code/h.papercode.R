# Want to run 5 times 5 years and get the info for the first
# plot, and collect $wantfoo for "error bands" around
# the curve.
NNN <- 2

#Y = 1
#ND = 1825=5yrs
#ND = 10950=30yrs
#ND = 12775=35yrs

foo <- matrix(0, NNN, 12775)
#foo <- matrix(0, NNN, ND)

for (i in 1:NNN) {
  cat("Step", i, "\n")
  set.seed(i)
  ans <- hovsim(pdfname=paste("test", i, ".pdf", sep=""))
  foo[i,] <- ans$wantfoo
}
if (NNN>1) write.table(foo, "confsim.csv", sep=",", row.names=F, col.names=T)

# Use a single run for the paper (WLOG the first of the 1000)
set.seed(1)
ans <- hovsim(pdfname="final.pdf")
write.csv(ans$all, "finalsample.csv", row.names=F)
