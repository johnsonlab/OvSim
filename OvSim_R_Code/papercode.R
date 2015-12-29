# Want to run 1000 times and get the info for the first
# plot, and collect $wantfoo for "error bands" around
# the curve.
NNN <- 1000
foo <- matrix(0, NNN, 420)
for (i in 1:NNN) {
  cat("Step", i, "\n")
  set.seed(i)
  ans <- ovsim(pdfname=paste("test", i, ".pdf", sep=""))
  foo[i,] <- ans$wantfoo
}
if (NNN>1) write.table(foo, "confsim.csv", sep=",", row.names=F, col.names=T)

# Use a single run for the paper (WLOG the first of the 1000)
set.seed(1)
ans <- ovsim(pdfname="final.pdf")
write.csv(ans$all, "finalsample.csv", row.names=F)
