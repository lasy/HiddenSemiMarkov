

#' Provides the list of currently supported sojourn distributions and the list of their parameters.
#'
#' @return A data.frame with the following columns: \code{distribution_type} and \code{parameters}
#'
#' @export
#' @examples
#' available_sojourn_dist()

available_sojourn_dist = function(){
  available_sojourn = rbind(
    data.frame(distribution_type = "nonparametric", parameters = "d", stringsAsFactors = FALSE),
    data.frame(distribution_type = "ksmoothed_nonparametric", parameters = "d", stringsAsFactors = FALSE),
    data.frame(distribution_type = "gamma", parameters = "shape, scale", stringsAsFactors = FALSE),
    data.frame(distribution_type = "poisson", parameters = "shift, lambda", stringsAsFactors = FALSE),
    data.frame(distribution_type = "lnorm", parameters = "meanlog, sdlog", stringsAsFactors = FALSE),
    data.frame(distribution_type = "logarithmic", parameters = "shape", stringsAsFactors = FALSE),
    data.frame(distribution_type = "nbinom", parameters = "size, mu or prob, shift", stringsAsFactors = FALSE)
  )
  available_sojourn
}





.get_longest_sojourn = function(model){

  # check sojourn
  model$sojourn = .check_sojourn(model$sojourn, model$J, model$state_names)

  p = 0.99
  Ms = purrr::map_int(
    .x = model$sojourn,
    .f = .get_longest_sojourn_for_state_j,
    p = p
  )
  M = max(Ms)
  M
}


.get_longest_sojourn_for_state_j = function(state_sojourn, p){
  # shortcut
  sojourn_distribution = state_sojourn$type

  # finding M
  if(sojourn_distribution == "nonparametric" | sojourn_distribution == "ksmoothed_nonparametric" ){
    M = max(which(state_sojourn$d > 0))
  }else if(sojourn_distribution == "poisson"){
    M = max(state_sojourn$shift + qpois(p = p, lambda = state_sojourn$lambda))
  }else if(sojourn_distribution == "lnorm"){
    M = max(qlnorm(p = p, meanlog = state_sojourn$meanlog , sdlog = state_sojourn$sdlog))
  }else if(sojourn_distribution == "gamma"){
    M = max(qgamma(p = p, shape = state_sojourn$shape , scale = state_sojourn$scale))
  }else if(sojourn_distribution == "logarithmic"){
    M = 10
    while(all(.dlog(M,state_sojourn$shape) > 0.01)){M = M*2}
  }else if(sojourn_distribution %in% c("nbinom")){
    if(is.null(state_sojourn$mu)) M = max(qnbinom(p = p, prob = state_sojourn$prob))
    if(is.null(state_sojourn$mu))  M = max(qnbinom(p = p, mu = state_sojourn$mu))
  }else{
    stop("This sojourn distribution is currently not supported.")
  }

  ceiling(M) %>% as.integer()
}

#' @export
.build_d_from_sojourn_dist <- function(model,M) {
  # shortcuts
  J = model$J

  d = purrr::map_dfc(
    .x = 1:J,
    .f = function(j){
      sojourn_this_state = model$sojourn[[j]]

      if(sojourn_this_state$type %in% c("nonparametric","ksmoothed_nonparametric")){
        this_state_d = sojourn_this_state$d
        if(length(this_state_d)<M) this_state_d = c(this_state_d, rep(0, M - length(this_state_d)))
        if(length(this_state_d)>M) this_state_d = head(this_state_d, M)
      }

      if(sojourn_this_state$type == "gamma")
        this_state_d = dgamma(1:M, shape = sojourn_this_state$shape, scale = sojourn_this_state$scale)

      if(sojourn_this_state$type == "poisson"){
        if(is.null(sojourn_this_state$shift)) sojourn_this_state$shift = 1
        this_state_d = .dpois.hsmm.sojourn(1:M,sojourn_this_state$lambda,sojourn_this_state$shift)
      }

      if(sojourn_this_state$type == "lnorm")
        this_state_d = dlnorm(1:M, meanlog = sojourn_this_state$meanlog, sdlog = sojourn_this_state$sdlog)

      if(sojourn_this_state$type == "logarithmic")
        this_state_d = .dlog(1:M, sojourn_this_state$shape)

      if(sojourn_this_state$type == "nbinom"){
        if(is.null(sojourn_this_state$shift)) sojourn_this_state$shift = 1
        if(is.null(sojourn_this_state$mu))
          this_state_d = .dnbinom.hsmm.sojourn(1:M, size = sojourn_this_state$size, prob = sojourn_this_state$prob, shift = sojourn_this_state$shift)
        if(is.null(sojourn_this_state$prob))
          this_state_d = .dnbinom.hsmm.sojourn(1:M, size = sojourn_this_state$size, mu = sojourn_this_state$mu, shift = sojourn_this_state$shift)
      }

      if(sum(this_state_d) != 0) this_state_d = this_state_d/sum(this_state_d) else this_state_d = 0*this_state_d
      this_state_d = data.frame(this_state_d) %>% magrittr::set_colnames(j)
    }
  )
  d %>% as.matrix()
}





############# SOJOURN DIST UTIL FUNCTIONS


.rlog <- function(n,p,M=10000) sample(1:M,n,TRUE,.dlog(1:M,p))

.dlog <- function(x,p) (-1 / log(1-p)) * p^x / x


.logdistrfit <- function(wt) {
  xbar = sum(wt*(1:length(wt)))
  fn <- function(p) xbar + p/((1-p)*log(1-p))
  uniroot(fn,c(1e-10,1-1e-10))$root
}


# rng for a shifted Poisson distribution
.rpois.hsmm.sojourn <- function(n,lambda,shift) rpois(n,lambda)+shift

# density function for a shifted Poisson distribution
.dpois.hsmm.sojourn <- function(x=NULL,lambda,shift,log=FALSE)   {
  if(shift<0) stop(".dpois.hsmm.sojourn: shift must be > 0")
  if(log) dpois(x-shift,lambda,log=TRUE)
  else dpois(x-shift,lambda)
}




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



#' Estimates gamma distribution parameters using method of moments
#' @param x a vector of values that need to be described by a Gamma distribution.
#' @param wt (optional) a vector of weights with values in [0,1] and of the same length of \code{x}.
#' @export
#'
gammafit = function(x , wt=NULL ) {
  tol = 1e-08

  if(is.null(wt)) wt = rep(1,length(x))

  if(sum(wt != 0)==1){
    # if only one value in wt is different than 0, the cov will be NaN
    xhat = x[which(wt != 0)] # was xhat = which(wt != 0)
    shape = 1000*xhat^2
  }else{
    tmp = cov.wt(data.frame(x), wt = wt)
    xhat = tmp$center
    xs = sqrt(tmp$cov)
    s = log(xhat) - mean(weighted.mean(log(x), wt))
    aold = (xhat/xs)^2
    a = Inf
    while(abs(a-aold) > tol) {
      a = aold - (log(aold) - digamma(aold) - s)/((1/aold) - trigamma(aold))
      aold = a
    }
    shape = a
  }
  list(shape = shape, scale = xhat/shape)
}





# .get_minimum_sojourn = function(model = model){
#   sojourn.distribution = model$sojourn$type
#   if(sojourn.distribution %in% c("nonparametric","ksmoothed_nonparametric")){
#     M = max(apply(model$sojourn$d, 2, function(x) max(which(x > 0))))
#   }else if(sojourn.distribution == "gamma"){
#     M = max(model$sojourn$shape * model$sojourn$scale + 3 * model$sojourn$shape^0.5 * model$sojourn$scale) # mean + 3 sd
#   }else if(sojourn.distribution == "poisson"){
#     M = max(model$sojourn$shift + model$sojourn$lambda + 3 * model$sojourn$lambda^0.5 ) # shift +  mean + 3 sd
#   }else if(sojourn.distribution == "lnorm"){
#     M = max(model$sojourn$meanlog + 3*model$sojourn$sdlog ) #  mean + 3 sd #CHECK THIS ONE
#   }else if(sojourn.distribution == "logarithmic"){
#     M = max(1 + 3*10*model$sojourn$shape^2 ) #  mean + 3 sd  #CHECK THIS ONE
#   }else if(sojourn.distribution == "nbinom"){
#     M = max(model$sojourn$shift + model$sojourn$size) #  shift + size
#   }
#   M
# }
