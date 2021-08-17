fitBBlogReg = function(G_count, A_count, mdl.mtx, dispersion=NA) {
  
  bb.negll <- function(par){
    
    if (is.na(dispersion)) {
      covariates.vector = par[1:length(par)-1]
      rho = par[length(par)]
    } else {
      covariates.vector = par[1:length(par)]
      rho = dispersion
      # rho = 1/((1/dispersion)+1)
    }
    
    linear.predictor = mdl.mtx%*%covariates.vector
    
    mu = 1/(1 + exp(-linear.predictor))
    
    mu[mu==0] = 1e-9
    mu[mu==1] = 0.99999999
    mu[is.na(mu)] = 1e-9
    
    precision = (1/rho) - 1
    

    alpha = mu * precision
    beta = (1-mu) * precision
    
    LL = -sum( lgamma(alpha+beta) - lgamma(alpha) - lgamma(beta) - lgamma(alpha+beta+n) + lgamma(alpha+k) + lgamma(n-k+beta))
    
    return(LL)
  }
  
  k = G_count
  n = G_count + A_count
  
  if (is.na(dispersion)) {
    num.params = dim(mdl.mtx)[2]
    mle.lin.model = nlminb(objective = bb.negll, start = c(rep(0, num.params), 0.333),
                           lower = c(rep(-Inf, num.params), 1/Inf), upper = c(rep(Inf, num.params), 0.99999),
                           hessian=FALSE)
    return(mle.lin.model$par)
    
  } else {
    num.params = dim(mdl.mtx)[2]
    mle.lin.model = ucminf(fn = bb.negll, par=rep(0, num.params),
                           hessian=2)
    reg.params = mle.lin.model$par
    fisher.info.mtx = mle.lin.model$invhessian
    return(list(reg.params, fisher.info.mtx))
  }
  
}

BBwaldTest = function(reg.params, fisher.info.mtx) {
  SEs = sqrt(diag(fisher.info.mtx))
  wald.statistic = reg.params[2]/SEs[2]
  p.value = 2*pnorm(-abs(wald.statistic))
  log.odds.ratio = reg.params[2]  
  return(c(log.odds.ratio, p.value))
}

getBBdisperiosns = function(G.counts, A.counts, mdl.mtx) {
  
  num.params = dim(mdl.mtx)[2]
  
  sites = row.names(G.counts)
  
  dispersions = c()
  
  for (site in sites) {
    G_test = as.vector(unlist(G.counts[site,])) 
    A_test = as.vector(unlist(A.counts[site,]))
    dispersions = c(dispersions, fitBBlogReg(G_test, A_test, mdl.mtx)[num.params+1])
  }
  
  return(dispersions)
} 

estMAPdisps = function(G.counts, A.counts, dispersions) {
  total_counts = A.counts+G.counts
  log_mean_total_counts = log10(apply(total_counts, 1, mean))
  local.df = data.frame(log_mean_total_counts, dispersion=dispersions)
  
  trend.df = local.df[local.df$dispersion > 1e-4,]
  
  beta.trend = betareg(dispersion ~ log_mean_total_counts, data = trend.df, link = 'logit')
  dispersions.fit = predict(beta.trend, local.df)
  
  local.df$dispersion.fit = dispersions.fit
  
  beta.map = function(params, i) {
    a = params[1]
    b = params[2]
    L1 = dbeta(dispersions[i], a, b, log = TRUE)
    L2 = dbeta(dispersions.fit[i], a, b, log = TRUE)
    
    return(-sum(L1+L2))
  }
  
  dispersions.map = c()
  
  for (i in seq(dim(local.df)[1])) {
    tmp.params = optim(par=c(1,1), fn=beta.map, i=i)$par
    a.tmp = tmp.params[1]
    b.tmp = tmp.params[2]
    dispersion.map = a.tmp/(a.tmp + b.tmp)
    dispersions.map = c(dispersions.map, dispersion.map)
  }
  
  local.df$dispersion.map = dispersions.map
  
  return(local.df)
}

getFinalResTable = function(G.counts, A.counts, mdl.mtx, dispersions.map) {
  
  sites = row.names(G.counts)
  
  log.odds.ratios = c()
  p.values = c()
  
  for (site in sites) {
    G_test = as.vector(unlist(G.counts[site,])) 
    A_test = as.vector(unlist(A.counts[site,]))
    dispersion = dispersions.map[site,]$dispersion.map
    
    params.and.cov.mtx = fitBBlogReg(G_test, A_test, mdl.mtx, dispersion)
    test.res = BBwaldTest(params.and.cov.mtx[[1]], params.and.cov.mtx[[2]])
    log.odds.ratios = c(log.odds.ratios, test.res[1])
    p.values = c(p.values, test.res[2])
  }
  
  
  diff_editing_res = data.frame(site=sites, log.odds.ratio=log.odds.ratios, p.value=p.values)
  

  diff_editing_res$adj.p.value = p.adjust(diff_editing_res$p.value, method="BH")
  
  diff_editing_res = diff_editing_res[order(diff_editing_res$adj.p.value),]
  
  return(diff_editing_res)
}


