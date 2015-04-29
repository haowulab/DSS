############################################################################
## Wrapper function for Fitting multiple factor DM model.
## Take a bsseq object and a design matrix
############################################################################
DMLfit.multiFactor <- function(BSobj, design, formula) {
    ## some checking to make sure input data is correct
    if(length(sampleNames(BSobj)) != nrow(design))
        stop("Dimension of data and design don't match. ")

    ## make design matrix out of formula
    X <- model.matrix(formula, design)
    if(nrow(X) <= ncol(X))
        stop("No enough degree of freedom to fit the linear model. Drop some terms in formula.")

    ## take counts from BSobj
    N0 <- getBSseq(BSobj, "Cov")
    Y0 <- getBSseq(BSobj, "M")

    ## fit the model
    fit <- DMLfit.multiFactor.engine(Y0, N0, X)

    ## return
    list(gr=getBSseq(BSobj, "gr"), design=design, formula=formula, fit=fit)

}


############################################################################
## Engine function for Fitting multiple factor DM model.
############################################################################
DMLfit.multiFactor.engine <- function(Y0, N0, X0) {
    if( !is.matrix(Y0) | !is.matrix(N0) )
        stop("Y and N need to be matrices.\n")

    ## get dimensions
    p = NCOL(X0)
    n = NROW(X0)
    C = nrow(Y0)

    ## loop over CpG sites
    result = vector("list", C)
    beta = se.beta = matrix(NA, nrow=C, ncol=p)
    cat("Fitting DML model for CpG site: ")
    for(i in 1:C) {
        if(i %% 1e5 == 0)
            cat(i, ", ")
        ## take counts for current CpG
        Y = Y0[i,]
        N = N0[i,]
        X = X0
        ## fit model
        tmp = DMLfit.oneCG(Y, N, X, n, p)
        if(is.null(tmp)) next
        ## save point estimates and SE for this CpG
        beta[i,] = tmp$beta0
        se.beta[i,] = tmp$se.beta0
    }

    list(beta=beta, se.beta=se.beta)
}

##############################################################
## DML model fitting for one CpG
## This is the "core" function and all methods are in here!!
##############################################################

DMLfit.oneCG <- function(Y, N, X, n, p) {
    ## small constants to bound p and phi
    c0 = 0.1
    c1 = 0.001

    ## check to make sure data is complete
    ix <- N > 0
    if(mean(ix) < 1) { ## has missing entries
        X <- X[ix,,drop=FALSE]
        Y <- Y[ix]
        N <- N[ix]
        ## check design
        if(nrow(X) < ncol(X) + 1) ## not enough df for regression
            return(NULL)
        if(any(abs(svd(X)$d) <1e-8)) ## design is not of full rank because of missing. Skip
            return(NULL)
    }

    ## Transform the methylation levels. Add a small constant to bound away from 0/1.
    Z = asin(2*(Y+c0)/(N+2*c0) - 1)

    ## First round of weighted least square.
    ## QR decomposition has to be performed for each CpG because it's WLS!!!
    XTVinv = t(X * N)
    beta0 = solve(XTVinv %*% X) %*% (XTVinv %*% Z) ### this parenthesis is also helpful to speed up

    ## get dispersion estimates, and restrict a bit to bound away from 0/1.
    phiHat = (sum( (Z - X %*% beta0)^2 * N) - (n - p)) * n / (n - p) / sum(N-1)
    phiHat = min(max(c1, phiHat),1-c1)

    ## Shrinkage phiHat a bit --- how to do this???

    ## second round of regression.
    XTVinv = t(X * (N/(1+(N-1)*phiHat)))  ###t(X)%*%VInv
    XTVinvX.inv = solve(XTVinv %*% X)
    beta0 = solve(XTVinv %*% X) %*% (XTVinv %*% Z)
    se.beta0 = sqrt(diag(XTVinvX.inv))

    ## return
    list(beta0=beta0, se.beta0=se.beta0)
}

##############################################################
### hypothesis testing function
##############################################################
DMLtest.multiFactor <- function(DMLfit, coef=2) {
    ## figure out index of the factor to be tested
    ## coef = find.factor(DMLfit$design, DMLfit$formula, factor)

    ## hypothesis testing
    fit = DMLfit$fit
    tmpRes = DMLtest.multiFactor.engine(fit, coef)

    ## return a data frame
    gr = DMLfit$gr
    res = data.frame(chr=seqnames(gr), pos=start(gr), tmpRes)
    invisible(res)
}

##############################################################
## engine function for testing
##############################################################
DMLtest.multiFactor.engine <- function(testResult, coef) {
    ## obtain point estimates and SE
    betas = testResult$beta[,coef]
    ses = testResult$se.beta[,coef]

    ## Wald test, get p-values and FDR
    stat = betas / ses
    pvals = 2*(1- pnorm(abs(stat)))
    fdrs = p.adjust(pvals, method="BH")

    ## return a matrix
    res = cbind(stat, pvals, fdrs)
    invisible(res)
}

##############################################################
## find index of factor to be tested.
## Below works only if there's no -1 (model HAS intercept).
## If there's no intercept this will become tricky.
## I need to think how to do that.
##############################################################
find.factor <- function(design, formula, factor) {
    formula.terms <- terms(formula)
    formula.labels <- attr(formula.terms, "term.labels")
    ix <- grep(factor, formula.labels)
    ix+1
}

