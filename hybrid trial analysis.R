bayes.two.prop <- function(y0, n0, y1, n1, prior0 = c(0.5, 0.5), 
                           prior1 = c(0.5, 0.5)){
  post0 <- c(y0 + prior0[1], n0 - y0 + prior0[2])
  post1 <- c(y1 + prior1[1], n1 - y1 + prior1[2])
  theta0 <- rbeta(10000, post0[1], post0[2])
  theta1 <- rbeta(10000, post1[1], post1[2])
  
  list("est" = post1[1]/sum(post1) - post0[1]/sum(post0),
       "ci" = quantile(theta1 - theta0, probs = c(0.025, 0.975)),
       "pp" = beta.ineq(a = post1[1], b = post1[2], 
                        c = post0[1], d = post0[2], delta = 0),
       "n0" = sum(post0))
}

## -------------------------------------------------------------------
## numerical comparison of randomly distributed beta variables
## -------------------------------------------------------------------
##x~beta(a,b), y~beta(c,d)
##P(X>y+delta)
beta.ineq <- function(a, b, c, d, delta)
{ 
  if (a <=0 | b <= 0 | c <=0 | d <= 0) 
    stop("paramters has to be positive")
  if (a <0.01 | b < 0.01 | c <0.01 | d < 0.01) 
    stop("paramters are to close to 0")
  if (delta>1) 
    stop("delta>1!")
  if (delta<0) 
    stop("delta<0!")
  
  integrand <- function(x) { dbeta(x, a, b)*pbeta(x-delta, c, d) }
  integrate(integrand, delta, 1, rel.tol=1e-4)$value
}
