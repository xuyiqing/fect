make_binary_panel <- function(seed=410) {
    set.seed(seed)
    N <- 20L; TT <- 22L
    d <- expand.grid(time=seq_len(TT), id=seq_len(N))
    d$D <- as.integer(d$id > 15 & d$time >= 15 + (d$id %% 3))
    d$X <- rnorm(nrow(d), sd=0.3)
    u <- rnorm(N, sd=0.25); t <- rnorm(TT, sd=0.2)
    d$Y <- rbinom(nrow(d),1,pnorm(d$X+u[d$id]+t[d$time]+0.3*d$D))
    d
}

binary_test_fit <- function(d, ...) {
    suppressMessages(fect(Y~D+X,data=d,index=c("id","time"),binary=TRUE,
        force="two-way",QR=TRUE,parallel=FALSE,tol=1e-3,...))
}

test_that("binary initializer and fixed-rank fits handle all additive FE choices", {
    d <- make_binary_panel()
    for (force in c("none","unit","time","two-way")) {
        for (QR in c(FALSE,TRUE)) {
            fit <- suppressMessages(fect(Y~D+X,data=d,index=c("id","time"),
                binary=TRUE,force=force,QR=QR,parallel=FALSE,CV=FALSE,r=0,se=FALSE))
            expect_true(all(is.finite(fit$Y.ct.full)))
            expect_true(all(fit$Y.ct.full >= 0 & fit$Y.ct.full <= 1))
        }
    }
})

test_that("rolling buffers leave the requested usable history", {
    II <- matrix(1L,15,10); D <- matrix(0L,15,10)
    folds <- .build_cv_mask_rolling(II,D,20,3,2,0.2,5,seed=12)
    for(f in folds) {
        train<-II; train[f$cv.id]<-0L
        expect_true(all(colSums(train)>=5))
        expect_length(f$est.id,6)
        expect_true(all(f$est.id %in% f$cv.id))
        expect_gt(length(f$cv.id),length(f$est.id))
    }
    expect_error(.build_cv_mask_rolling(matrix(1L,8,10),matrix(0L,8,10),
        2,3,1,.2,5),"no eligible units")
})

test_that("excluded outcomes cannot affect binary initialization or fitting", {
    d<-make_binary_panel(); Y<-matrix(d$Y,22,20); D<-matrix(d$D,22,20)
    X<-array(d$X,c(22,20,1)); II<-1-D
    f<-.build_cv_mask_rolling(II,D,1,2,1,.2,5,seed=5)[[1]]
    mask<-II; mask[f$cv.id]<-0L
    fit<-.fect_binary_estimate(Y,X,mask,1,3,TRUE,1e-3)
    altered<-Y; altered[mask==0]<-1-altered[mask==0]
    again<-.fect_binary_estimate(altered,X,mask,1,3,TRUE,1e-3)
    expect_equal(fit$fit,again$fit,tolerance=1e-12)
})

test_that("binary CV scores probabilities only on designated holdouts", {
    d<-make_binary_panel()
    fit<-binary_test_fit(d,CV=TRUE,r=c(0,1),k=2,cv.nobs=2,
        cv.buffer=1,cv.prop=.2,min.T0=5,seed=42,cv.rule="min",se=FALSE)
    expect_identical(fit$cv.loss,"probability_mspe")
    expect_equal(fit$cv.counts,c(8L,8L))
    Y<-fit$Y.dat; D<-fit$D.dat
    X<-array(d$X,c(22,20,1))
    for (b in 1:2) {
        f<-fit$cv.folds[[b]]; mask<-1-D; mask[f$cv.id]<-0L
        z<-.fect_binary_estimate(Y,X,mask,0,3,TRUE,1e-3)
        expected<-mean((Y[f$est.id]-pnorm(z$fit[f$est.id]))^2)
        expect_equal(fit$cv.loss.per.fold[1,b],expected,tolerance=1e-10)
    }
    expect_equal(fit$CV.out[,"MSPE"],rowMeans(fit$cv.loss.per.fold))
    expect_equal(fit$r.cv,unname(fit$CV.out[which.min(fit$CV.out[,"MSPE"]),"r"]))
    again<-binary_test_fit(d,CV=TRUE,r=c(0,1),k=2,cv.nobs=2,
        cv.buffer=1,cv.prop=.2,min.T0=5,seed=42,cv.rule="min",se=FALSE)
    expect_identical(again$cv.folds,fit$cv.folds)
    expect_equal(again$CV.out,fit$CV.out)
})

test_that("binary inference restrictions are enforced", {
    d<-make_binary_panel()
    expect_error(binary_test_fit(d,CV=FALSE,r=0,se=TRUE,vartype="parametric"),
                 "only nonparametric bootstrap or jackknife")
    ## treatment reversals are supported (fect 1.1.x behaviour)
    d$D[d$id==20 & d$time==22]<-0
    fit <- binary_test_fit(d,CV=FALSE,r=0,se=FALSE)
    expect_true(is.numeric(fit$att.avg))
    expect_true(length(fit$att.off) >= 1)
})

test_that("binary CV settings survive jackknife and standalone dispatch", {
    withr::local_envvar(TESTTHAT="true")
    d<-make_binary_panel()
    common<-list(d=d,CV=TRUE,r=c(0,0),k=2,cv.nobs=2,cv.prop=.2,
                 cv.buffer=2,min.T0=6,seed=18,cv.rule="min")
    plain<-do.call(binary_test_fit,c(common,list(se=FALSE)))
    jack<-do.call(binary_test_fit,c(common,list(se=TRUE,vartype="jackknife")))
    expect_equal(jack$cv.settings,plain$cv.settings)
    expect_equal(jack$CV.out,plain$CV.out)
    expect_equal(jack$att.avg,plain$att.avg,tolerance=1e-10)
    expect_true(all(is.finite(jack$est.avg)))
    separate<-suppressMessages(r.cv.rolling(Y~D+X,data=d,index=c("id","time"),
        binary=TRUE,method="ife",force="two-way",QR=TRUE,r.max=0,k=2,
        cv.nobs=2,cv.buffer=2,cv.prop=.2,min.T0=6,seed=18,cv.rule="min",
        parallel=FALSE,tol=1e-3))
    expect_equal(separate$mspe.per.fold,plain$cv.loss.per.fold)
})

test_that("binary block CV and nonparametric bootstrap use the same point fit", {
    withr::local_envvar(TESTTHAT="true")
    d<-make_binary_panel()
    common<-list(d=d,CV=TRUE,r=c(0,0),k=2,cv.method="block",cv.nobs=3,
                 cv.donut=1,cv.prop=.1,min.T0=5,seed=33,cv.rule="min")
    plain<-do.call(binary_test_fit,c(common,list(se=FALSE)))
    boot<-do.call(binary_test_fit,c(common,list(se=TRUE,vartype="bootstrap",nboots=8)))
    expect_equal(boot$CV.out,plain$CV.out)
    expect_equal(boot$cv.folds,plain$cv.folds)
    expect_equal(boot$att.avg,plain$att.avg,tolerance=1e-10)
    expect_true(all(is.finite(boot$est.avg)))
    expect_true(all(plain$cv.counts > 0))
    expect_true(all(vapply(plain$cv.folds,function(f)
        length(f$est.id) < length(f$cv.id),logical(1))))
})

test_that("failed binary fits cannot win CV by dropping difficult folds", {
    d<-make_binary_panel()
    expect_error(binary_test_fit(d,CV=TRUE,r=c(0,1),k=2,cv.nobs=2,
        seed=11,max.iteration=1,se=FALSE),"All binary CV candidates failed")
})

test_that("positive-rank binary fits work without covariates on incomplete panels", {
    d <- make_binary_panel()
    d <- d[-c(3,37,88,122,250,363),]
    for (force in c("none","unit","time","two-way")) for (QR in c(FALSE,TRUE)) {
        fit <- suppressMessages(fect(Y~D,data=d,index=c("id","time"),
            binary=TRUE,force=force,QR=QR,CV=FALSE,r=1,se=FALSE,
            parallel=FALSE,tol=1e-3))
        expect_true(is.finite(fit$att.avg))
        expect_true(all(is.finite(fit$Y.ct.full)))
        expect_true(all(fit$Y.ct.full >= 0 & fit$Y.ct.full <= 1))
    }
})

test_that("positive-rank binary inference keeps the selected model", {
    withr::local_envvar(TESTTHAT="true")
    d <- make_binary_panel()
    plain <- binary_test_fit(d,CV=FALSE,r=1,se=FALSE)
    for (vartype in c("bootstrap","jackknife")) {
        fit <- binary_test_fit(d,CV=FALSE,r=1,se=TRUE,vartype=vartype,
            nboots=8,seed=19)
        expect_equal(fit$r.cv,1)
        expect_equal(fit$att.avg,plain$att.avg,tolerance=1e-10)
        expect_true(all(is.finite(fit$est.avg)))
    }
})
