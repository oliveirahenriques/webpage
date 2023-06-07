QMLE.bootstrap <- function(fit, data, n.bootfit = 999, n.ahead = 1, tau = 0.05, alpha = 0.05){
  
  # initial paarameters
  
  # fit = garchx(data, order = c(1,1))
  
  # n.bootfit = 1000
  # 
  # tau = 0.05
  
  n.bootpred = n.ahead 
  
  N = length(data)
  
  set.seed(1)
  
  # ---------------------------------------
  
  # generate paths of equal length to data based on empirical re-sampling of z
  # Pascual, Romo and Ruiz (2006) p.2296 equation (5)
  
  fz = as.numeric(residuals(fit))
  
  empz = matrix(0, ncol = n.bootfit, nrow = N)
  
  empz = apply(as.data.frame(1:n.bootfit), 1, FUN=function(i){
    
    sample(fz, N, replace = TRUE)
    
  })
  
  # presigma uses the same starting values as the original fit
  # in paper they use alternatively the unconditional long run sigma 
  # Pascual, Romo and Ruiz (2006) p.2296 equation (5) (P.2296 paragraph 2  "...marginal variance..."
  
  coef = as.numeric(coef(fit))
  
  spec =  
    ugarchspec(
      mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
      variance.model = list(model = 'sGARCH', garchOrder = c(1,1)),  
      fixed.pars = list(
        mu = 0, # our mu (intercept)
        ar1 = 0, # our phi_1 (AR(1) parameter of mu_t)
        ma1 = 0, # our theta_1 (MA(1) parameter of mu_t)
        omega = coef[1], # our alpha_0 (intercept)
        alpha1 = coef[2], # our alpha_1 (ARCH(1) parameter of sigma_t^2)
        beta1 = coef[3])) # our beta_1 (GARCH(1) parameter of sigma_t^2)
  
  presigma = tail(sqrt(fitted(fit)),1)
  
  prereturns = tail(data, 1)
  
  preresiduals = tail(fz, 1) 
  
  paths = ugarchpath(spec, 
                     n.sim = N, 
                     m.sim = n.bootfit, 
                     presigma = presigma,
                     prereturns = prereturns,
                     preresiduals = preresiduals,
                     n.start = 0, 
                     custom.dist = list(name = "sample", distfit = as.matrix(empz)))
  
  fitlist = vector(mode="list", length = n.bootfit)
  
  simseries = fitted(paths)
  
  nx = NCOL(simseries)
  
  # generate path based forecast values
  # for each path we generate n.bootpred vectors of resampled data of length n.ahead
  # Equation (6) in the PRR paper (again using z from original fit)
  #-------------------------------------------------------------------------
  
  fitlist = lapply(as.list(1:nx), FUN = function(i){
    
    fit.boot = garchx(y = as.numeric(simseries[,i]), order = c(1,1))
    
    theta = coef(fit.boot)
    
    df = tibble('ID' = i,
                'omega' = theta[1],
                'alpha' = theta[2], 
                'beta' = theta[3],
                'epsilon' = data[length(data):1],
                'eta' = c(as.numeric(residuals(fit.boot)),NA), 
                'j' = seq(0,length(data)-1,1),
                'sum' = beta^j*(lead(epsilon)^2 - omega/(1-alpha-beta))) %>% 
      drop_na() %>% 
      reframe(
        ID = unique(ID),
        omega = unique(omega),
        alpha = unique(alpha),
        beta = unique(beta),
        epsilon = first(epsilon), 
        sum = sum(sum),
        sigma2.hat = omega + alpha*epsilon^2 + beta*(omega/(1-alpha - beta) + alpha*sum),
        extremile.hat = sqrt(sigma2.hat)*extremile(eta, probs = tau),
        expectile.hat = sqrt(sigma2.hat)*expectile(eta, probs = tau),
        VaR.hat = sqrt(sigma2.hat)*quantile(eta, probs = tau),
        ES.hat = sqrt(sigma2.hat)*mean(if_else(eta < quantile(eta, probs = tau),eta,NA),na.rm = TRUE))
    
    return(df)
    
  })
  
  confidence.interval = 
    fitlist %>% 
    bind_rows() %>% 
    select(ID,contains('hat')) %>% 
    pivot_longer(-ID, names_to = 'measure', values_to = 'estimate') %>% 
    group_by(measure) %>% 
    summarise_at(vars(estimate), 
                 .funs = list(lower_bound = ~ quantile(., probs = alpha/2),
                              upper_bound = ~ quantile(., probs = 1-alpha/2))) %>% 
    left_join(
      tibble(eta = as.numeric(residuals(fit)),
             sigma2.hat = as.numeric(predict(fit,n.ahead = 1))) %>%
        reframe(
          ID = 1,
          sigma2.hat = unique(sigma2.hat),
          extremile.hat = sqrt(sigma2.hat)*extremile(eta, probs = tau),
          expectile.hat = sqrt(sigma2.hat)*expectile(eta, probs = tau),
          VaR.hat = sqrt(sigma2.hat)*quantile(eta, probs = tau),
          ES.hat = sqrt(sigma2.hat)*mean(if_else(eta < quantile(eta, probs = tau),eta,NA),na.rm = TRUE)) %>% 
        pivot_longer(-ID,names_to = 'measure', values_to = 'estimate') %>% 
        select(-ID)) %>% 
    relocate(measure,lower_bound,estimate,upper_bound)
  
  rm(fitlist)
  
  rm(paths)
  
  gc(verbose = FALSE)
  
  return(confidence.interval)
  
}
