
getARL <- function(L, lambda, m, nu, n, ubc, tmin, tmax, state = 'steady') {
 
  if (state == 'steady') {
    
    pracma::integral2(fun = integrandSteady, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
                L = L, lambda = lambda, mm = m, nu = nu, ubc = ubc, tmin = tmin, tmax = tmax)$Q
    
  } else if (state == 'zero') {
   
    pracma::integral2(fun = integrandZero, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
              L = L, lambda = lambda, mm = m, nu = nu, nn = n, ubc = ubc, tmin = tmin, tmax = tmax)$Q
     
  }
  
}


getCC.CUC <- function(ARL0, interval = c(1, 4), lambda, m, nu = m - 1, n = 5, ubCons = 1, t.interval = c(20, 50), state = 'steady', tol = 1e-2) {
  
  rootFinding <- function(L, ARL0, lambda, m, nu, n, ubc, tmin, tmax, state = 'steady') {
    
    ARL <- getARL(L, lambda, m, nu, n, ubc, tmin, tmax, state)
    
    cat('L:', L,'and ARLin:', ARL, '\n')
    
    ARL0 - ARL
    
  }
  
  tmin <- t.interval[1]
  tmax <- t.interval[2]
 
  ubc <- ubCons
 
  uniroot(rootFinding, interval = interval, ARL0 = ARL0, lambda = lambda, m = m, nu = nu, n = n, ubc = ubc, 
            tmin = tmin, tmax = tmax, state = state, tol = tol)$root
  
}

c4.f <- function(nu) sqrt(2 / nu) * 1 / beta(nu / 2, 1 / 2) * sqrt(pi) 


pCARLEWMA <- function(q, L, lambda, m, nu, n, ubc, tmin, tmax, state = 'steady') {
  
  if (state == 'steady') {
    
    pracma::integral2(fun = qintegrandSteady, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
              q= q, L = L, lambda = lambda, mm = m, nu = nu, ubc = ubc, tmin = tmin, tmax = tmax)$Q
    
  } else if (state == 'zero') {
    
    pracma::integral2(fun = qintegrandZero, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
              q = q, L = L, lambda = lambda, mm = m, nu = nu, nn = n, ubc = ubc, tmin = tmin, tmax = tmax)$Q
    
  }
  
}



getCC.EPC <- function(p0, interval = c(1, 4), ARL0, epstilda, lambda, m, nu = m - 1, n = 5, ubCons = 1, t.interval = c(20, 50), state = 'steady', tol = 1e-2) {
  
  rootFinding <- function(L, p0, ARLb, lambda, m, nu, n, ubc, tmin, tmax, state) {
    
    p <- pCARLEWMA(ARLb, L, lambda, m, nu, n, ubc, tmin, tmax, state)
    cat('L:', L, 'and p:', p, '\n')
    p0 - p
    
  }
 
  tmin <- t.interval[1]
  tmax <- t.interval[2]
 
  ubc <- ubCons
 
  eps <- epstilda
 
  ARLb <- (1 - eps) * ARL0
  
  uniroot(rootFinding, interval = interval, p0 = p0, ARLb = ARLb, lambda = lambda, m = m, nu = nu, n = n, ubc = ubc, 
          tmin = tmin, tmax = tmax, state = state, tol = tol)$root
  
}


getCC <- function(
          m,
          nu = m - 1,
		  n = 5,
          ARL0 = 370,
          interval = c(1, 4),
		  lambda = 0.1,
          EPC.p0 = 0.05,
          EPC.epstilda = 0,
          cc.option = c('EPC'),
		  ubCons = 1, 
		  t.interval = c(20, 50),
		  state = 'steady',
		  tol = 1e-2) 
{

  if (cc.option == 'CUC') {

    cc <- getCC.CUC(
      ARL0 = ARL0, 
	  interval = interval, 
	  lambda = lambda, 
	  m = m, 
	  nu = nu, 
	  n = n, 
	  ubCons = ubCons, 
	  t.interval = t.interval, 
	  state = state, 
	  tol = tol
    )

  } else if (cc.option == 'EPC') {

    cc <- getCC.EPC(
	  p0 = EPC.p0, 
	  interval = interval, 
	  ARL0 = ARL0, 
	  epstilda = EPC.epstilda, 
	  lambda = lambda, 
	  m = m, 
	  nu = nu, 
	  n = n, 
	  ubCons = ubCons, 
	  t.interval = t.interval, 
	  state = state, 
	  tol = tol
	)
	
  }

  return(cc)

}



PH2EWMA <- function(
  X2,
  X1,
  cc = NULL,
  ARL0 = 370,
  interval = c(1, 4),
  lambda = 0.1,
  EPC.p0 = 0.05,
  EPC.epstilda = 0,
  cc.option = c('CUC', 'EPC'),
  t.interval = c(20, 50),
  state = 'steady',
  tol = 1e-2,
  ubCons.option = TRUE,
  plot.option = TRUE) 
{

  m <- dim(X1)[1]
  nu <- m - 1
  n <- dim(X1)[2]
  
  if (ubCons.option == TRUE) {

    ubCons = c4.f(nu)

  } else {

    ubCons = 1

  }

  m2 <- dim(X2)[1]

  X1bar <- rowMeans(X1)
  X1barbar <- mean(X1bar)
  X1Var <- var(X1bar)

  X2bar <- rowMeans(X2)

  if (is.null(cc)) {

    cc <- rep(NA, 2)
    txt1 <- rep(NA, 2)
	txt2 <- rep(NA, 2)
	
    lower.limits <- rep(NA, 2)
    upper.limits <- lower.limits

    if (length(which(cc.option == 'CUC')) > 0) {

      cc[1] <- getCC.CUC(
					ARL0 = ARL0, 
					interval = interval, 
					lambda = lambda, 
					m = m, 
					nu = nu, 
					n = n, 
					ubCons = ubCons, 
					t.interval = t.interval, 
					state = state, 
					tol = tol
				)
			  
	   txt1[1] <- 'CUC'
	   txt2[1] <- 'CUC'

    }

    if (length(which(cc.option == 'EPC')) > 0) {

      cc[2] <- getCC.EPC(
					p0 = EPC.p0, 
					interval = interval, 
					ARL0 = ARL0, 
					epstilda = EPC.epstilda, 
					lambda = lambda, 
					m = m, 
					nu = nu, 
					n = n, 
					ubCons = ubCons, 
					t.interval = t.interval, 
					state = state, 
					tol = tol
				)
	    
		txt1[2] <- 'EPC'
		txt2[2] <- 'EPC'

    }

  } else {
  
		cc.num <- length(cc)
		txt1 <- rep(NA, cc.num)
		txt2 <- rep(NA, cc.num)
		
		for (i in 1:cc.num) {
			txt1[i] <- paste('LCL', i, sep = '')
			txt2[i] <- paste('UCL', i, sep = '')
		}
  }

  cc <- cc[!is.na(cc)]

  cc.num <- length(cc)

  lower.limits <- X1barbar - cc * sqrt(X1Var) / ubCons
  upper.limits <- X1barbar + cc * sqrt(X1Var) / ubCons

  if (plot.option == TRUE) {

    plot(c(1, m2), c(min(X2bar, lower.limits), max(X2bar, upper.limits)), type = 'n',
          xlab = 'Subgroup', ylab = 'Sample Mean')
    points(1:m2, X2bar, type = 'o', lty = 1)

    for (i in 1:cc.num) {

      abline(h = lower.limits[i], lty = i)
      text(round(m2 * 0.8), lower.limits[i], paste(txt1[i], ' = ', round(lower.limits[i], 4)), pos = 3)

      abline(h = upper.limits[i], lty = i)
      text(round(m2 * 0.8), upper.limits[i], paste(txt2[i], ' = ', round(upper.limits[i], 4)), pos = 1)

    }

  }

  out <- list(
            CL = X1barbar,
            sigma = sqrt(X1Var) / ubCons,
            PH2.cc = cc,
            LCL = lower.limits,
            UCL = upper.limits,
            CS = X2bar)

  return(out)

}