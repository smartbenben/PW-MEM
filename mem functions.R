## ---------------------------------------------------
## functions for MEM
## ---------------------------------------------------
require("partitions")

# R: the number of trials to consider
# max_cl: maximum number of blocks
get.part <- function(R, max_cl){
  ## generate all the possible partitions 
  ## and store them in a matrix 
  part_mat <- t(setparts(R))
  part_mat <- part_mat[apply(part_mat, 1, function(x){
    length(unique(x))<=max_cl
  }), ]
  part_mat <- data.frame(part_mat)
  names(part_mat) <- LETTERS[1:R]
  return(part_mat)
}


set.mem.prior <- function(num_study, delta){
  part <- get.part(R = num_study, max_cl = num_study)
  K <- nrow(part)
  # number of blocks/unique response rate in each partition
  n_bk <- apply(part, 1, function(x){
    length(unique(x))
  })
  prior <- n_bk^(delta)/sum(n_bk^(delta))
  
  list("part" = part, "prior" = prior)
}

##x: vector of responses in each basket
##n: vector of sample sizes in each basket
##beta prior for response rate: a0, b0 [default: a0 = b0 = 1]
##max_cl: maximum number of clusters a basket trial should be partitioned to
update.part.bin <- function(x, n, prior_part, part, a0 = 1, b0 = 1){
  R <- length(x)
  K <- nrow(part)
  
  p <- foreach(k = 1:K,.combine = "c")%do%{
    grp <- unlist(part[k, ])
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #calculate marginal probs m(s_j)
    prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_part[k]
  }
  post_part <- p/sum(p)
  
  idx <- which.max(post_part)
  
  sim_mat <- part[, 1:(R-1)]==part[, R]
  
  if(is.null(nrow(sim_mat))){
    post_sim <- sum(sim_mat*post_part)
  }else{
    post_sim <- colSums(sim_mat*post_part)
  }
  
  return(list("part_hat" = part[idx, ], "phat" = post_part[idx], 
              "post_part" = post_part,
              "post_sim" = post_sim))
}

gaussian.mar <- function(y, source, wt){
  n <- ybar <- s <- NULL
  for(i in 0:1){
    n[i+1] <- sum(wt[source==i])
    ybar[i+1] <- sum(y[source==i]*wt[source==i])/n[i+1]
    s[i+1] <- sqrt((wt[source==i]%*%c(y[source==i]-ybar[i+1])^2)/(n[i+1]))
  }
  
  y0 <- y[source==0]
  y1 <- y[source==1]
  wt0 <- wt[source==0]
  wt1 <- wt[source==1]
  
  d <- NULL
  
  a <- 0.5*sum(n/s^2)
  b <- sum(n*ybar/s^2)
  d[1] <- sqrt(pi/a)*exp(b^2/(4*a))
  
  a <- 0.5*(n/s^2)
  b <- (n*ybar/s^2)
  d[2] <- prod(sqrt(pi/a)*exp(b^2/(4*a)))
  
  list("mden" = d, "ybar" = ybar, "sd" = s, "n" = n)
}


##alex's method
calc.mar <- function(y, source, wt){
  n <- ybar <- s <- yss <- mden <- v <- NULL
  for(i in 0:1){
    n[i+1] <- sum(wt[source==i])
    ybar[i+1] <- sum(y[source==i]*wt[source==i])/n[i+1]
    yss[i+1] <- sum(y[source==i]^2*wt[source==i])
    s[i+1] <- sqrt((wt[source==i]%*%c(y[source==i]-ybar[i+1])^2)/(n[i+1]-1))
    v[i+1] <- s[i+1]^2/n[i+1]
  }
  ## marginal probability 
  mden[1] <- sqrt(2*pi)^1 / sqrt((1/v[1] + 1/v[2])) * exp(-0.5 * ((ybar[1]-ybar[2])^2/(v[1] + v[2])) )
  mden[2] <- sqrt(2*pi)^2 / sqrt(1/(v[1]*v[2]))
  #num <- n^(-0.5)
  #den <- (s^(n-1))*(2*pi)^(0.5*(n-1))
  #(num/den)*exp( (-yss+n*ybar^2)/(2*s^2) )
  #m1 <- calc.mar(y)
  #m2 <- calc.mar(y[source==1])*calc.mar(y[source==0])
  #c(m1, m2)*fit0$prior/sum(c(m1, m2)*fit0$prior)
  list("mden" = mden, "ybar" = ybar, "sd" = s, "n" = n)
}
