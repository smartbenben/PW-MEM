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

