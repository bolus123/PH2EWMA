
getARL <- function(L, lambda, m, nu, n, ubc, tmin, tmax, state = 'steady') {
 
  if (state == 'steady') {
    
    pracma::integral2(fun = integrandSteady, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
                L = L, lambda = lambda, mm = m, nu = nu, ubc = ubc, tmin = tmin, tmax = tmax)$Q
    
  } else if (state == 'zero') {
   
    pracma::integral2(fun = integrandZero, xmin = 0, xmax = 1, ymin = 0, ymax = 1, singular = TRUE, vectorized = FALSE, 
              L = L, lambda = lambda, mm = m, nu = nu, nn = n, ubc = ubc, tmin = tmin, tmax = tmax)$Q
     
  }
  
}


getCCEWMAUncond <- function(ARL0, interval, lambda, m, nu, n, ubc, t.interval = c(20, 50), state = 'steady', tol = 1e-2) {
  
  rootFinding <- function(L, ARL0, lambda, m, nu, n, ubc, tmin, tmax, state = 'steady') {
    
    ARL <- getARL(L, lambda, m, nu, n, ubc, tmin, tmax, state)
    
    cat('L:', L,'and ARLin:', ARL, '\n')
    
    ARL0 - ARL
    
  }
  
  tmin <- t.interval[1]
  tmax <- t.interval[2]
 
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



getCCEWMACond <- function(p0, interval = c(1, 4), ARL0, eps, lambda, m, nu, n, ubc, t.interval = c(20, 50), state = 'steady', tol = 1e-4) {
  
  rootFinding <- function(L, p0, ARLb, lambda, m, nu, n, ubc, tmin, tmax, state) {
    
    p <- pCARLEWMA(ARLb, L, lambda, m, nu, n, ubc, tmin, tmax, state)
    cat('L:', L, 'and p:', p, '\n')
    p0 - p
    
  }
 
  tmin <- t.interval[1]
  tmax <- t.interval[2]
 
  ARLb <- (1 - eps) * ARL0
  
  uniroot(rootFinding, interval = interval, p0 = p0, ARLb = ARLb, lambda = lambda, m = m, nu = nu, n = n, ubc = ubc, 
          tmin = tmin, tmax = tmax, state = state, tol = tol)$root
  
}
