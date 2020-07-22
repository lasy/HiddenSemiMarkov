
.fitnbinom <- function(eta) {
  shiftthresh=1e-20
  j = which(eta > shiftthresh)
  maxshift = min(j)
  Mtmp = max(j)

  optimize_nbinom_par = function(shift, eta, U){
    shifted_U = U-shift
    m <- weighted.mean(shifted_U,eta[U])
    v <- cov.wt(data.frame(shifted_U),wt=eta[U])$cov %>% as.numeric()
    size <- ifelse(v > m, m^2/(v - m), 100)
    densfun <- function(par) sum(dnbinom(shifted_U, size=par[1], mu=par[2], log=TRUE)*eta[U])
    opt_res = optim(c(size,m),densfun,control=list(fnscale=-1))
    c(value = opt_res$value, size = opt_res$par[1], mu = opt_res$par[2])
  }

  if(Mtmp>maxshift){
    optimal_shift_and_par = sapply(1:maxshift,function(shift) suppressWarnings(optimize_nbinom_par(shift, eta = eta, U = maxshift:Mtmp)))
    shift = which.max(optimal_shift_and_par[1,]) %>%  as.numeric()
    size = optimal_shift_and_par[2,shift] %>%  as.numeric()
    mu = optimal_shift_and_par[3,shift] %>%  as.numeric()
  }else{shift = maxshift; size = 1; mu = 0}

  c(shift = shift,size=size,mu=mu,prob=size/(size+mu))
}


.dnbinom.hsmm.sojourn <- function(x,size,prob=NULL,shift,mu=NULL,log=FALSE) {
  if(shift<0) stop(".dnbinom.hsmm.sojourn: shift must be > 0")
  if(is.null(mu)){
    if(log) dnbinom(x-shift,size,prob,log=TRUE)
    else dnbinom(x-shift,size,prob)
  }
  else {
    if(log) dnbinom(x-shift,size=size,mu=mu,log=TRUE)
    else dnbinom(x-shift,size=size,mu=mu)
  }
}

.rnbinom.hsmm.sojourn <- function(n,size,prob,shift) {
  if(shift<0) stop(".dnbinom.hsmm.sojourn: shift must be > 0")
  rnbinom(n,size,prob) + shift
}



#' estimates gamma distribution parameters using method of moments
#' @export
gammafit <- function(x,wt=NULL) {
  tol = 1e-08

  if(is.null(wt)) wt = rep(1,length(x))

  if(sum(wt != 0)==1){
    # if only one value in wt is different than 0, the cov will be NaN
    xhat = which(wt != 0)
    shape = 1000*xhat^2
  }else{
    tmp = cov.wt(data.frame(x),wt=wt)
    xhat = tmp$center
    xs = sqrt(tmp$cov)
    s = log(xhat) - mean(weighted.mean(log(x),wt))
    aold = (xhat/xs)^2
    a = Inf
    while(abs(a-aold)>tol) {
      a = aold - (log(aold) - digamma(aold) - s)/((1/aold) - trigamma(aold))
      aold=a
    }
    shape = a
  }
  list(shape=shape,scale=xhat/shape)
}
