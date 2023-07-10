

trial.analysis.binary <- function(sim_dta, method, ps_fm = NULL, y_fm = NULL, 
                           fit0 = NULL){
  out <- sapply(sim_dta, function(tmp){
    ## summarize internal data
    IN_y0 <- sum(tmp$response[(tmp$trt==0)&(tmp$source==1)])
    IN_N0 <- sum((tmp$trt==0)&(tmp$source==1))
    
    IN_y1 <- sum(tmp$response[(tmp$trt==1)&(tmp$source==1)])
    IN_N1 <- sum((tmp$trt==1)&(tmp$source==1))
    
    if(method == "IN_only"){
      tst <- bayes.two.prop(y0 = IN_y0, n0 = IN_N0, y1 = IN_y1, n1 = IN_N1)
    }
    if(method == "pull"){
      tst <- bayes.two.prop(y0 = sum(tmp$response[tmp$trt==0]), n0 = sum(tmp$trt==0), 
                            y1 = sum(tmp$response[tmp$trt==1]), n1 = sum(tmp$trt==1))
    }
    if(method == "msm"){
      #msm_fit <- msm.boots(ps_mod_fm = ps_fm, y_mod_fm = y_fm, 
      #                     dat = tmp, B = 1000)
      ## msm:
      #svyglm always returns 'model-robust' standard errors; 
      #the Horvitz-Thompson-type standard errors used everywhere in 
      # the survey package are a generalisation of the model-robust 
      #'sandwich' estimators. 
      
      ps_mod <- glm(ps_fm, data = tmp[tmp$trt==0, ], family = binomial)
      ps <- predict(ps_mod, newdata = tmp, type = "response")
      wt <- ifelse(tmp$source == 1, 1, ps/(1-ps))
      svyD <- svydesign(~ 1, weights = ~ wt, data = tmp)
      msm <- svyglm(y_fm, design = svyD, 
                    family = quasibinomial(link = "logit"))
      #tst <- unlist(msm_fit$te)
      #aggregate(wt, by = list(tmp$source), sum)
      tst <- c(msm$coef[2], confint.default(msm,"trt"), summary(msm)$coef[2, 3])
    }
    if(method == "mem"){
      ps_mod <- glm(ps_fm, data = tmp[tmp$trt==0, ], family = binomial)
      ps <- predict(ps_mod, newdata = tmp, type = "response")
      wt <- ifelse(tmp$source == 1, 1, ps/(1-ps))
      # These weights are used to reweight outcomes of patients in both groups 
      # to produce pseudo-samples of control and experimental patients with 
      # (approximately) common distribution of pre-treatment characteristics
      EX_wy <- sum(wt*tmp$response*(1-tmp$source))
      EX_wn <- sum(wt*(1-tmp$source))
      # assess exchangeability between the weighted pseudo control samples
      # and the internal controls
      fit <- update.part.bin(x = c(IN_y0, EX_wy), 
                             n = c(IN_N0, EX_wn), 
                             prior_part = fit0$prior, 
                             part = fit0$part)
      # post_part[1]: the power parameter in power prior
      # (0.5, 0.5): the initial prior
      a0_star <- EX_wy*fit$post_part[1]
      b0_star <- (EX_wn-EX_wy)*fit$post_part[1]
      #r <- min(1, target_n/(a0_star+b0_star))
      tst <- bayes.two.prop(y0 = IN_y0+a0_star, n0 = IN_N0+EX_wn*fit$post_part[1], 
                            y1 = IN_y1, n1 = IN_N1,
                            prior0 = c(0.5, 0.5), 
                            prior1 = c(0.5, 0.5))
    }
    unlist(tst)
  })
  t(out)
}

trial.analysis.normal <- function(sim_dta, method, ps_fm = NULL, y_fm = NULL, 
                                  fit0 = NULL, alt = "greater"){
  out <- sapply(sim_dta, function(tmp){
    ## summarize internal data
    IN_y0 <- sum(tmp$response[(tmp$trt==0)&(tmp$source==1)])
    IN_N0 <- sum((tmp$trt==0)&(tmp$source==1))
    
    IN_y1 <- sum(tmp$response[(tmp$trt==1)&(tmp$source==1)])
    IN_N1 <- sum((tmp$trt==1)&(tmp$source==1))
    
    if(method == "IN_only"){
      mod <- lm(response ~ trt, data = tmp[tmp$source==1, ])
      tst <- c(mod$coef[2], confint(mod)[2, ], summary(mod)$coef[2, 3])
    }
    if(method == "pull"){
      mod <-  lm(response ~ trt, data = tmp)
      tst <- c(mod$coef[2], confint(mod)[2, ], summary(mod)$coef[2, 3])
    }
    if(method == "msm"){
      ps_mod <- glm(ps_fm, data = tmp[tmp$trt==0, ], family = binomial)
      ps <- predict(ps_mod, newdata = tmp, type = "response")
      tmp$wt <- ifelse(tmp$source == 1, 1, ps/(1-ps))
      
      svyD <- svydesign(~ 1, weights = ~ wt, data = tmp)
      msm <- svyglm(y_fm, design = svyD)
      tst <- c(msm$coef[2], confint(msm)[2, ], summary(msm)$coef[2, 3])
      
      ## boostrap takes more time...
      if(0>1){
      msm_boot <- sapply(1:1000, function(x){
        idx <- sample(1:nrow(tmp), nrow(tmp), replace = T)
        bdta <- tmp[idx, ]
        
        ps_mod <- glm(ps_fm, data = bdta, family = binomial)
        ps <- predict(ps_mod, type = "response")
        bdta$wt <- ifelse(bdta$source == 1, 1, ps/(1-ps))
        
        y_mod <- lm(y_fm, weights = wt, data = bdta)
        y_mod$coef[2]
      })}
      #tst <- c(mean(msm_boot), quantile(msm_boot, c(0.025, 0.975)),
      #mean(msm_boot)/sd(msm_boot))
    }
    if(method == "mem"){
      ps_mod <- glm(ps_fm, data = tmp[tmp$trt==0, ], family = binomial)
      ps <- predict(ps_mod, newdata = tmp, type = "response")
      tmp$wt <- ifelse(tmp$source == 1, 1, ps/(1-ps))
    
      mod <- gaussian.mar(y = tmp$response[tmp$trt==0], 
                          source = tmp$source[tmp$trt==0],
                          wt = tmp$wt[tmp$trt==0])
      post_part <- mod$mden*fit0$prior/sum(mod$mden*fit0$prior)
      
      #r <- post_part[1]*min(1, target_n/(mod$n[1]*post_part[1]))
      r <- post_part[1]
      Ap <- mod$n[2]/mod$sd[2]^2+r*mod$n[1]/mod$sd[1]^2
      Bp <- mod$n[2]*mod$ybar[2]/mod$sd[2]^2 + r*mod$n[1]*mod$ybar[1]/mod$sd[1]^2
      mu_ctr <- Bp/Ap
      sigma_ctr <- sqrt(1/Ap)
      ## assume the variance for the treatment group is known(sample variance)
      ## assume the mean of the treatment group has N(0, sd = 100)
      s_trt <- sd(tmp$response[tmp$trt==1])
      ybar_trt <- mean(tmp$response[tmp$trt==1])
      prec_trt <- 1/100^2 + IN_N1/s_trt^2
      sigma_trt <- prec_trt^(-0.5)
      mu_trt <-(ybar_trt*IN_N1/s_trt^2)/prec_trt
      
      ## their difference
      qt <- qnorm(c(0.5, 0.025, 0.975), mean = mu_trt-mu_ctr, 
            sd = sqrt(sigma_trt^2 + sigma_ctr^2))
      pp <- 1-pnorm(0, mean = mu_trt-mu_ctr, 
            sd = sqrt(sigma_trt^2 + sigma_ctr^2))
      tst <- c(qt, pp)
    }
    unlist(tst)
  })
  t(out)
}

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

summary.simu <- function(sim_out, TE, gam){
  sim_out$rej <- sim_out$pp>gam
  sim_out$cover <- (sim_out$lci<TE)&(sim_out$uci>TE)
  sim_out$bias <- sim_out$est - TE
  
  pow <- aggregate(rej ~ method + scenario, data = sim_out, mean)
  
  mse <- aggregate(bias^2 ~ method + scenario, data = sim_out, mean)
  
  cp <- aggregate(cover ~ method + scenario, data = sim_out, mean)
  
  #en0 <- aggregate(n0 ~ method + scenario, data =sim_out, mean)
  
  list("pow" = pow, "cp" = cp, "mse" = mse, "out" = sim_out)
}


trial.analysis.map <- function(sim_dta, binary, map_rob){
  out <- sapply(sim_dta, function(tmp){
    ## summarize internal data
    IN_y0 <- sum(tmp$response[(tmp$trt==0)&(tmp$source==1)])
    IN_N0 <- sum((tmp$trt==0)&(tmp$source==1))
    
    IN_y1 <- sum(tmp$response[(tmp$trt==1)&(tmp$source==1)])
    IN_N1 <- sum((tmp$trt==1)&(tmp$source==1))
    
    if(binary){
      m0 <- postmix(map_rob, n = IN_N0, r = IN_y0)
      m1 <- mixbeta(c(1, 1+IN_y1, 1+IN_N1-IN_y1))
    }else{
      IN_s0 <- sd(tmp$response[(tmp$trt==0)&(tmp$source==1)])
      IN_s1 <- sd(tmp$response[(tmp$trt==1)&(tmp$source==1)])
      m0 <- postmix(map_rob, m = IN_y0/IN_N0, 
                    se = IN_s0/sqrt(IN_N0))
      weak_prior <- mixnorm(c(1, 0, 10), param = "ms")
      m1 <- postmix(weak_prior, m = IN_y1/IN_N1, 
                    se = IN_s1/sqrt(IN_N1))
    }
   c(qmixdiff(m1, m0, c(0.5, 0.025, 0.975)),
         pmixdiff(m1, m0, 0, FALSE))
  })
  t(out)
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
