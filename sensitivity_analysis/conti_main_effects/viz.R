load('./sensitivity_analysis/ikss_me12.RData')

nvar <- 7
G <- 30
grid <- seq(0, 1, length=G)

m <- me2$m
q1 <- me2$q1
q2 <- me2$q2

plot(0, xlab="normalized range of parameters", ylab="main effect", xlim=c(0,1), ylim=range(c(q1,q2)), type="n")
for(j in 1:nvar) {
  lines(grid, m[,j], col=j, lwd=2)
}
legend("bottomright", paste0("x",1:nvar), fill=1:nvar, horiz=TRUE, cex=0.75)
