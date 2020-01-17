# Based on Susan Homes' function at
# https://web.stanford.edu/class/stats366/Gibbs.html
gibbsNimate <- function(NSimul = 50, rho = 0.5)
{
  ##NSimul=Number of simulations
  ##rho is the correlation
  
  SEQ <- seq(from = -3, to = 3, length.out = 1000)
  grid_vals <- expand.grid(x1 = SEQ, x2 = SEQ)
  sigma <- sqrt(1 - rho^2)
  d <- with(grid_vals, 1 / (2 * pi * sigma) * exp(-0.5 * (x1^2 + x2^2 - 2 * rho * x1 * x2)))
  contour(x = SEQ, y = SEQ, z = matrix(d, length(SEQ), length(SEQ)),
          main = "Bivariate Normal", col = "green", las = 1, 
          xlab = expression(x[1]), ylab = expression(x[2]), bg = "transparent")
  dev.hold()

  ## initialization
  x <- rep(NA_real_, NSimul)
  y <- rep(NA_real_, NSimul)
  x[1] <- rnorm(1)
  y[1] <- rnorm(1)
  
  ## Gibbs sampler and plot points
  
  for(b in 2:NSimul) {
    contour(x = SEQ, y = SEQ, z = matrix(d, length(SEQ), length(SEQ)),
            main = "Bivariate Normal", col = "green", las = 1, 
            xlab = expression(x[1]), ylab = expression(x[2]), bg = "transparent")
    points(x[1:(b-1)], y[1:(b-1)])
    dev.hold()
    ## sample the x direction conditional on y
    x[b] <- rnorm(1, mean = rho * y[b - 1], sd = sigma)
    lines(c(x[b-1], x[b]), c(y[b - 1], y[b - 1]))
    # ani.pause()
    ## sample the y direction conditional on x
    y[b] <- rnorm(1, mean = rho * x[b], sd = sigma)
    lines(c(x[b], x[b]), c(y[b - 1], y[b]))
    points(x[b], y[b],col="red")
    # ani.pause()
  }
  
  ## return the draws
  return(cbind(x1 = x, x2 = y))
}

# library(animation)
# oopt <- ani.options(interval = 2)
draws <- gibbsNimate(25, rho = 0.5)
