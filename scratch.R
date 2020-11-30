## libraries and custom functions
library(plgp)

mylhs <- function(n, m) {
  # generate the Latin hypercube 
  l <- seq(-(n - 1)/2, (n - 1)/2, by = 1)
  L <- matrix(NA, nrow = n, ncol = m)
  for(j in 1:m) {
    L[, j] <- sample(l, n)
  }
  
  # draw the random uniforms and turn the hypercube into a sample
  U <- matrix(runif(n*m), ncol = m)
  X <- (L + (n - 1)/2 + U)/n
  colnames(X) <- paste0("x", 1:m)
  
  # return the design and the grid it lives on for visualization
  return(list(X = X, g = c((l + (n - 1)/2)/n,1)))
}

mymaxmin <- function(n, m, T = 100000) {
  X <- matrix(runif(n*m), ncol = m)
  d <- distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  
  for(t in 1:T) {
    row_idx <- sample(1:n, 1)
    xold <- X[row_idx, ]
    X[row_idx] <- runif(m)
    
    d <- distance(X)
    d <- d[upper.tri(d)]
    mdprime <- min(d)
    
    if(mdprime > md) {md <- mdprime} 
    else {X[row, ] <- xold}
  }
}

mymaxmin_seq <- function(n, m, T = 100000, Xorig = NULL) {
  X <- matrix(runif(n*m), ncol = m)
  d <- distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  
  if(!is.null(Xorig)) {
    md2 <- min(distance(X, Xorig))
    if(md2 < md) {md <- md2}
  }
  
  for (t in 1:T) {
    row_idx <- sample(1:n, 1)
    xold <- X[row_idx, ]
    X[row, ] <- runif(m)
    
    d <- distance(X)
    d <- d[upper.tri(d)]
    midprime <- min(d)
    
    if(!is.null(Xorig)) {
      mdprime2 <- min(distance(X, Xorig))
      if(mdprime2 < mdprime) {mdprime <- mdprime2}
    }
    
    if(mdprime > md) {f} 
    else {X[row, ] <- xold}
  }
  
  return(X)
}


## Latin hypercube sample -----
# uniform
m <- 2
n <- 10
X <- matrix(runif(n*m), ncol = m)
colnames(X) <- paste0('x', 1:m)

plot(X, xlim = c(0, 1), ylim = c(0, 1))

# LHS
l <- seq(-(n-1)/2, (n-1)/2, by = 1)
L <- matrix(NA, nrow = n, ncol = m)
for (j in 1:m) {L[, j] <- sample(l, n)}
U <- matrix(runif(n*m), ncol = 2)
X <- (L + (n-1)/2 + U)/n

plot(X, xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2')
abline(h = c((l + (n-1)/2)/n, 1), col = 'grey', lty = 'dashed')
abline(v = c((l + (n-1)/2)/n, 1), col = 'grey', lty = 'dashed')


## maxmin designs -----
m <- 2
n <- 10

X1 <- matrix(runif(n*m), ncol = m)
dX1 <- distance(X1)
dX1 <- dX1[upper.tri(dX1)]
md1 <- min(dX1)

X2 <- matrix(runif(n*m), ncol = m)
dX2 <- distance(X2)
dX2 <- dX2[upper.tri(dX2)]
md2 <- min(dX2)

T <- 100000
for (t in 1:T) {
  row <- sample(1:n, 1)
  xold <- X1[row, ]
  X1[row, ] <- runif(m)
  
  d <- distance(X1)
  d <- d[upper.tri(d)]
  mdprime <- min(d)
  
  if(mdprime > md1) {
    md1 <- mdprime
  } else {
    X1[row, ] <- xold
  }
}

plot(X1, xlim = c(0, 1), ylim = c(0, 1), xlab = 'x1', ylab = 'x2')

