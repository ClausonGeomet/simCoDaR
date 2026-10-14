# Independent behavioral oracle: public API and analytical densities only.
ind_emp_data <- function(z) data.frame(A = plogis(sqrt(2)*z), B = 1-plogis(sqrt(2)*z))
ind_emp_density <- function(fit, x) {
  h <- fit$groups[[1]]$patterns[[1]]; m <- h$model
  z <- log(as.matrix(x[fit$parts])) %*% t(h$basis)
  switch(fit$model,
    dirichlet = dbeta(x[[1]], m$alpha[1], m$alpha[2], log=TRUE) + log(2)/2 + rowSums(log(x)),
    logistic_normal = dnorm(z[,1], m$mu[1], sqrt(m$sigma[1,1]), log=TRUE),
    ilr_student_t = dt((z[,1]-m$mu[1])/sqrt(m$sigma[1,1]), m$df, log=TRUE)-log(m$sigma[1,1])/2,
    ilr_mixture = log(Reduce(`+`, lapply(seq_along(m$weights), function(j)
      m$weights[j]*dnorm(z[,1],m$means[j,1],sqrt(m$covariances[[j]][1,1]))))))
}
ind_emp_metrics <- function(x,y,b,ref) {
  x <- as.matrix(x); y <- as.matrix(y); ref <- as.matrix(ref)
  x <- x/rowSums(x); y <- y/rowSums(y); ref <- ref/rowSums(ref)
  z <- log(x)%*%t(b); w <- log(y)%*%t(b)
  cross <- outer(z[,1], w[,1], function(a,b) abs(a-b))
  lx <- log(x[,1]/x[,2]); ly <- log(y[,1]/y[,2])
  qx <- quantile(lx,c(.05,.25,.5,.75,.95),type=7,names=FALSE)
  qy <- quantile(ly,c(.05,.25,.5,.75,.95),type=7,names=FALSE)
  bounds <- apply(ref,2,quantile,probs=c(.01,.99),type=7,names=FALSE)
  list(ilr_energy=2*mean(cross)-mean(abs(outer(z[,1],z[,1],`-`)))-mean(abs(outer(w[,1],w[,1],`-`))),
    ilr_covariance_relative_error=abs(var(z[,1])-var(w[,1]))/max(abs(var(z[,1])),1e-10),
    mean_absolute_logratio_quantile_error=mean(abs(qx-qy)), logratio_quantile=qy,
    logratio_lower_tail=mean(ly<qx[1]),logratio_upper_tail=mean(ly>qx[5]),
    reference_outside_cell_fraction=mean(sweep(y,2,bounds[1,],`<`)|sweep(y,2,bounds[2,],`>`)))
}
test_that('independent common-measure densities, metrics, MCSE and RNG hold', {
  train <- ind_emp_data(seq(-2,2,length.out=80)); held <- ind_emp_data(c(-1,.2,1.4))
  for (model in c('dirichlet','logistic_normal','ilr_student_t','ilr_mixture')) {
    fit <- fit_assay_simulator(train,c('A','B'),model,seed=4105,student_df=5,
      mixture_components=1:2,measurement_basis='abstract_composition')
    frozen <- serialize(fit,NULL)
    x <- data.frame(A=c(.1,.3,.5,.8)); x$B <- 1-x$A
    v <- validate_assay_simulator(fit,x,train,seed=4300)
    expect_lte(max(abs(v$scores$log_score-ind_emp_density(fit,x))),1e-9)
    expect_lte(abs(v$score_summary$mean_log_score-mean(ind_emp_density(fit,x))),1e-12)
    expect_identical(v$score_summary$measure,'ilr')
    set.seed(7401); before <- .Random.seed
    a <- validate_assay_simulator(fit,held,train,seed=4300)
    expect_identical(.Random.seed,before)
    expect_identical(a,validate_assay_simulator(fit,held,train,seed=4300))
    expect_identical(serialize(fit,NULL),frozen)
    expect_equal(nrow(a$failures),0L)
    for (r in 1:20) {
      y <- simulate_assays(fit,71,seed=4300+r)
      oracle <- ind_emp_metrics(held,y[c('A','B')],fit$groups[[1]]$patterns[[1]]$basis,train)
      rows <- a$replicate_metrics[a$replicate_metrics$replicate==r,]
      for (key in names(oracle)) expect_lte(max(abs(rows$generated[rows$metric==key]-oracle[[key]])),1e-9)
      expect_lte(max(abs(rowSums(y[c('A','B')])-1)),1e-12)
    }
    for (i in seq_len(nrow(a$mc_summary))) {
      s <- a$mc_summary[i,]; rows <- a$replicate_metrics
      same <- rows$metric==s$metric
      for (key in c('part_a','part_b','probability')) same <- same & if(is.na(s[[key]])) is.na(rows[[key]]) else !is.na(rows[[key]]) & rows[[key]]==s[[key]]
      values <- rows[[s$statistic]][same]; values <- values[is.finite(values)]
      expect_equal(s$n_success,length(values)); expect_equal(s$n_failure,20-length(values))
      if(length(values)>=2) expect_lte(abs(s$mcse-sd(values)/sqrt(length(values))),1e-12)
    }
    extra <- held; extra$confidential_marker <- 'DO_NOT_RETAIN_abc123'
    expect_identical(a,validate_assay_simulator(fit,extra,train,seed=4300))
    null <- validate_assay_simulator(fit,held,seed=4300)
    expect_identical(a$scores,null$scores)
    expect_true(all(is.na(subset(null$replicate_metrics,metric=='reference_outside_cell_fraction')$generated)))
    other <- validate_assay_simulator(fit,held,train[1:40,],seed=4300)
    keep <- a$replicate_metrics$metric!='reference_outside_cell_fraction'
    expect_identical(a$replicate_metrics[keep,],other$replicate_metrics[keep,])
    expect_identical(a$scores,validate_assay_simulator(fit,held,train,n_draws=12,replicates=2,seed=19)$scores)
  }
})
test_that('validation invalid inputs restore RNG; no initial RNG is created', {
  train <- ind_emp_data(seq(-2,2,length.out=80))
  fit <- fit_assay_simulator(train,c('A','B'),measurement_basis='abstract_composition')
  set.seed(701); before <- .Random.seed
  for(key in c('n_draws','replicates')) for(value in list(1,2.5,NA,Inf)) {
    args <- list(fit=fit,heldout=train); args[[key]] <- value
    expect_error(do.call(validate_assay_simulator,args)); expect_identical(.Random.seed,before)
  }
  for(seed in list(NULL,-1,1.5,NA,Inf,.Machine$integer.max,.Machine$integer.max-10)) {
    expect_error(validate_assay_simulator(fit,train,seed=seed)); expect_identical(.Random.seed,before)
  }
  bad <- list(train[FALSE,],train[1,,drop=FALSE],train['A'])
  for(value in c(0,-1,NA,Inf)) { x<-train; x$A[1]<-value; bad[[length(bad)+1]]<-x }
  x<-train; names(x)<-c('A','A'); bad[[length(bad)+1]]<-x
  x<-train; names(x)<-c('A',''); bad[[length(bad)+1]]<-x
  x<-train; names(x)<-c('A',NA); bad[[length(bad)+1]]<-x
  x<-train; x$A<-as.character(x$A); bad[[length(bad)+1]]<-x
  for(x in bad) { expect_error(validate_assay_simulator(fit,x)); expect_error(validate_assay_simulator(fit,train,reference=x)); expect_identical(.Random.seed,before) }
  expect_error(validate_assay_simulator(list(),train))
  grouped<-train; grouped$group<-rep(c('g','h'),each=40)
  gf<-fit_assay_simulator(grouped,c('A','B'),group_col='group')
  expect_error(validate_assay_simulator(gf,train),'ungrouped')
  face<-train; face$B<-0
  sf<-fit_assay_simulator(face,c('A','B'),zero_mode='structural',structural_parts='B'); expect_error(validate_assay_simulator(sf,train),'full positive')
  if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)
  a<-validate_assay_simulator(fit,train,seed=4300)
  expect_false(exists('.Random.seed',.GlobalEnv,inherits=FALSE))
  expect_identical(a,validate_assay_simulator(fit,train,seed=4300))
  expect_error(validate_assay_simulator(fit,train,seed=-1))
  expect_false(exists('.Random.seed',.GlobalEnv,inherits=FALSE))
})
