bootstrap = function(fitORspec, data = NULL, n.ahead = 10,
                     n.bootfit = 100, n.bootpred = 500, rseed = NA, 
                     solver = "solnp", solver.control = list(), fit.control = list(), 
                     external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
                     mexsimdata = NULL, vexsimdata = NULL, method = 'MLE'){
  
  require(xts)
  
  if(method == 'MLE'){
  
  fit = fitORspec
  
  model = fit@model
  
  vmodel = fit@model$modeldesc$vmodel
  
  m = model$maxOrder
  
  data = model$modeldata$data
  
  N = model$modeldata$T
  
  ns = model$n.start
  
  if(is.na(rseed[1])){
    sseed1 = NA
    sseed2 = NA
  } else{
    if(length(rseed) < n.bootpred){
      stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
    } else {
      sseed = rseed
      sseed1 = sseed[1:n.bootfit]
      sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
    }
  }
  
  # generate paths of equal length to data based on empirical re-sampling of z
  # Pascual, Romo and Ruiz (2006) p.2296 equation (5)
  
  fz = fit@fit$z
  
  empz = matrix(sample(fz, N, replace = TRUE), ncol = n.bootfit, nrow = N)
  
  # presigma uses the same starting values as the original fit
  # in paper they use alternatively the unconditional long run sigma 
  # Pascual, Romo and Ruiz (2006) p.2296 equation (5) (P.2296 paragraph 2  "...marginal variance..."
  
  paths = ugarchsim(fit, n.sim = N, m.sim = n.bootfit, 
                    presigma = tail(fit@fit$sigma, m), 
                    prereturns = tail(model$modeldata$data[1:N], m), 
                    preresiduals = tail(residuals(fit), m), 
                    startMethod = "sample", 
                    custom.dist = list(name = "sample", distfit = as.matrix(empz)),
                    rseed = sseed1, mexsimdata = mexsimdata, vexsimdata = vexsimdata)
  
  fitlist = vector(mode = "list", length = n.bootfit)
  
  simseries = fitted(paths)
  
  spec = getspec(fit)
  
  # help the optimization with good starting parameters
  
  spec@model$start.pars = as.list(coef(fit))
  
  nx = NCOL(simseries)
  
  # get the distribution of the parameters (by fitting to the new path data)
  #-------------------------------------------------------------------------
  
  fitlist = 
    ugarchfit(
      spec = spec,
      data = xts(as.numeric(simseries), as.Date(1:NROW(data), origin = "1970-01-01")),
      solver = solver,
      fit.control = fit.control, 
      solver.control = solver.control)
  
  boot.std.residuals = as.data.frame(residuals(fitlist, standardize = TRUE))
  
  ans = list(boot.std.residuals = boot.std.residuals)
  
  return(ans)
  
  }else{

    if(method == 'QMLE'){

      fit = 
        ugarchfit(
        spec =
          ugarchspec(
            mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
            variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
            distribution.model = 'norm'), data = data, solver = 'hybrid')
      
      model = fit@model
      
      vmodel = fit@model$modeldesc$vmodel
      
      m = model$maxOrder
      
      data = model$modeldata$data
      
      N = model$modeldata$T
      
      ns = model$n.start
      
      if(is.na(rseed[1])){
        
        sseed1 = NA
        
        sseed2 = NA
        
        }else{
          
          if(length(rseed) < n.bootpred){
            
            stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
            
            }else{
              
              sseed = rseed
              
              sseed1 = sseed[1:n.bootfit]
              
              sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
            }
          
      }
      
      # fitORspec = garchx::garchx(teste, order =  c(1,1))
      
      coefs = coef(fitORspec)
      
      spec = 
        ugarchspec(
          mean.model = list(armaOrder = c(0,0), include.mean = FALSE), # ARMA order
          variance.model = list(model = "sGARCH", garchOrder = c(1,1)), # GARCH order
          distribution.model = 'norm',
          fixed.pars =
            list(
              mu = 0, # our mu (intercept)
              ar1 = 0, # our phi_1 (AR(1) parameter of mu_t)
              ma1 = 0, # our theta_1 (MA(1) parameter of mu_t)
              omega = (20^2/252)*(1-coefs[2]-coefs[3]), # our alpha_0 (intercept)
              alpha1 = coefs[2], # our alpha_1 (ARCH(1) parameter of sigma_t^2)
              beta1 = coefs[3])) # our beta_1 (GARCH(1) parameter of sigma_t^2)

      # generate paths of equal length to data based on empirical re-sampling of z
      # Pascual, Romo and Ruiz (2006) p.2296 equation (5)
      
      # fz = fit@fit$z
      
      fz = residuals(fitORspec)
      
      empz = matrix(sample(fz, N, replace = TRUE), ncol = n.bootfit, nrow = N)
      
      # presigma uses the same starting values as the original fit
      # in paper they use alternatively the unconditional long run sigma 
      # Pascual, Romo and Ruiz (2006) p.2296 equation (5) (P.2296 paragraph 2  "...marginal variance..."
      
      paths = ugarchpath(spec, n.sim = N, m.sim = n.bootfit, 
                         presigma = tail(fitted(fitORspec), m), 
                         prereturns = tail(xdata[1:N], m), 
                         preresiduals = tail(residuals(fitORspec), m), 
                         rseed = sseed1, n.start = 0, 
                         custom.dist = list(name = "sample", distfit = as.matrix(empz)), 
                         mexsimdata = mexsimdata, vexsimdata = vexsimdata)
      
      fitlist = vector(mode = "list", length = n.bootfit)
      
      simseries = fitted(paths)
      
      spec = getspec(fit)
      
      # help the optimization with good starting parameters
      
      spec@model$start.pars = as.list(coef(fit))
      
      nx = NCOL(simseries)
      
      # get the distribution of the parameters (by fitting to the new path data)
      #-------------------------------------------------------------------------
      
      fitlist = 
        garchx(
          data = xts(as.numeric(simseries), as.Date(1:NROW(data), origin = "1970-01-01")),
          order = c(1,1))
        
      boot.std.residuals = as.data.frame(residuals(fitlist))
      
      ans = list(boot.std.residuals = boot.std.residuals)
      
      return(ans)

    }

  }
  
}