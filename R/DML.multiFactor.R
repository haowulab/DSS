############################################################################
## Wrapper function for Fitting multiple factor DM model.
## Take a bsseq object and a design matrix
############################################################################
DMLfit.multiFactor <- function(BSobj, design, formula, smoothing=FALSE, smoothing.span=500) {
    ## some checking to make sure input data is correct
    if(length(sampleNames(BSobj)) != nrow(design))
        stop("Dimension of data and design don't match. ")

    ## make design matrix out of formula
    X <- model.matrix(formula, design)
    if(nrow(X) <= ncol(X)+1)
        stop("No enough degree of freedom to fit the linear model. Drop some terms in formula.")

    ## take counts from BSobj
    N0 <- getBSseq(BSobj, "Cov")
    Y0 <- getBSseq(BSobj, "M")

    ## compute the response variable Z, which is transformed methylation level
    c0 = 0.1
    if(smoothing) { ## with smoothing. The mean methylation levels are smoothed
        allchr <- as.character(seqnames(BSobj))
        allpos <- start(BSobj)
        N0.sm = N0; Y0.sm = Y0
        for(i in 1:ncol(N0)) {
            N0.sm[,i] <- round(smooth.chr(as.double(N0[,i]), smoothing.span, allchr, allpos, "avg"))
            Y0.sm[,i] <- round(smooth.chr(as.double(Y0[,i]), smoothing.span, allchr, allpos, "avg"))
        }
        Z0 = asin(2*(Y0.sm+c0)/(N0.sm+2*c0) - 1)

    } else { ## no smoothing
        Z0 = asin(2*(Y0+c0)/(N0+2*c0) - 1)
    }

    ## fit the model
    fit <- DMLfit.multiFactor.engine(as.array(Y0), as.array(N0), X, as.array(Z0))

    ## return
    list(gr=getBSseq(BSobj, "gr"), design=design, formula=formula, X=X, fit=fit)

}


############################################################################
## Engine function for Fitting multiple factor DM model.
############################################################################
DMLfit.multiFactor.engine <- function(Y0, N0, X0, Z0) {

    if( (!is.matrix(Y0) | !is.matrix(N0)) )
        stop("Y and N need to be matrices.\n")

    ## get dimensions
    p = NCOL(X0)
    n = NROW(X0)
    C = nrow(Y0)

    ## loop over CpG sites
    beta = matrix(NA, nrow=C, ncol=p)
    var.beta = matrix(NA, nrow=C, ncol=p*p)
    phi = numeric(C)

    cat("Fitting DML model for CpG site: ")
    for( i in 1:C ) {
        if(i %% 1e5 == 0)
            cat(i, ", ")
        ## take counts for current CpG and fit model
        tmp = DMLfit.oneCG(Y0[i,], N0[i,], X0, Z0[i,], n, p)
        if(is.null(tmp)) next
        ## save point estimates and SE for this CpG
        beta[i,] = tmp$beta0
        var.beta[i,] = tmp$var.beta0
        phi[i] = tmp$phi
    }

    list(beta=beta, var.beta=var.beta, phi=phi)
}


##############################################################
## DML model fitting for one CpG
## This is the "core" function and all methods are in here!!
##############################################################

DMLfit.oneCG <- function(Y, N, X, Z, n, p) {
    ## small constants to bound p and phi
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
        Z <- Z[ix]
    }

    ## Transform the methylation levels. Add a small constant to bound away from 0/1.
    # Z = asin(2*(Y+c0)/(N+2*c0) - 1)

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

    ## return. I'll flatten the var/cov matrix for easy storing.
    list(beta0=beta0, se.beta0=se.beta0, var.beta0 = as.vector(XTVinvX.inv), phi=phiHat)
}

##############################################################
### hypothesis testing function
##############################################################
DMLtest.multiFactor <- function(DMLfit, coef=2, term, Contrast) {
    ## figure out index of the factor to be tested
    ## coef = find.factor(DMLfit$design, DMLfit$formula, factor)
    ## check inputs
    flag = 0
    if(!missing(coef))
        flag = flag + 1
    if(!missing(term))
        flag = flag + 1
    if(!missing(Contrast))
        flag = flag + 1

    if(flag == 0)
        stop("Must specify one of the following parameter for testing: coef, term, or Contrast.\n")
    if(flag > 1)
        stop("You can only specify one of the following parameter for testing: coef, term, or Contrast.\n")

    if(!missing(coef)) { # specified coef
        res = DMLtest.multiFactor.coef(DMLfit, coef)
    } else if(!missing(term)) { # specify term
        ## create a contrast matrix for testing the term
        Contrast = makeContrast(DMLfit, term)
        ## testing
        res = DMLtest.multiFactor.Contrast(DMLfit, Contrast)

    } else if(!missing(Contrast)) { # specify contrast
        ## check contrast matrix
        if( nrow(Contrast) != ncol(DMLfit$X) )
            stop("Input Contrast matrix has wrong dimension: its number of rows must match the number of columns of the design matrix.\n")

        ## testing
        res = DMLtest.multiFactor.Contrast(DMLfit, Contrast)

    }

    class(res)[2] = "DMLtest.multiFactor"
    invisible(res)

}

##############################################################
## Hypothesis testing when specify a coef for testing.
## This only tests one column in the design matrix.
## Wald test will be used.
##############################################################

DMLtest.multiFactor.coef <- function(DMLfit, coef) {
    if(is.character(coef)) {
        tmp = which(colnames(DMLfit$X) == coef)
        if(length(tmp) == 0)
            stop(paste0("Can't find terms to be tested: ", coef,
                        ". Make sure it matches a column name in design matrix."))
        coef = tmp
    }

    ## hypothesis testing
    p = ncol(DMLfit$X)
    fit = DMLfit$fit
    betas = fit$beta[,coef]
    ## take out SE estimates from var/cov matrix
    tmp = t(apply(fit$var.beta, 1, function(x) diag(matrix(x, ncol=p))))
    ses = sqrt(tmp[,coef])

    ## Wald test, get p-values and FDR
    stat = betas / ses
    pvals = 2*pnorm(-abs(stat))  #2*(1- pnorm(abs(stat)))
    fdrs = p.adjust(pvals, method="BH")

    ## return a data frame
    gr = DMLfit$gr
    res = data.frame(chr=seqnames(gr), pos=start(gr), stat, pvals, fdrs)
    invisible(res)
}

##############################################################
## Hypothesis testing when specify a contrast matrix.
## This tests multiple columns in the design matrix.
## F-test will be used.
##############################################################

DMLtest.multiFactor.Contrast <- function(DMLfit, Contrast) {
    p = ncol(DMLfit$X)
    fit = DMLfit$fit
    betas = fit$beta
    ## A^T * beta
    Abeta = betas %*% Contrast

    ## loop through CpG sites -- have to do this since the var/cov matrices of the beta estimates
    ## are different for each site
    stat = rep( NA, nrow(betas) )
    for( i in 1:nrow(betas) ) {
        Sigma = matrix(fit$var.beta[i,], ncol=p)
        tmp = solve(t(Contrast) %*% Sigma %*% Contrast)
        thisAbeta = Abeta[i,,drop=FALSE]
        stat[i] = thisAbeta %*% tmp %*% t(thisAbeta)
    }

    ## get the sign of the contrast if there's only one contrast.
    ## This is to be added to test statistics
    ## When Contrast has multiple rows, there won't be a sign for test statistics.
    if(nrow(Contrast) == 1)
        signs = sign(betas %*% Contrast)
    else signs = 1

    ## get p-values. Stat follows F_{r, R}
    ## I found that using F distribution, the p-values are pretty large.
    ## Use sqrt(f) and normal gives much smaller p-values,
    ## and this is consistent with the Wald test in two-group comparison.
    r = ncol(Contrast)
##     R = nrow(DMLfit$X) - ncol(DMLfit$X)
##     stat = stat / r
##     pvals = 1 - pf(stat, r, R)
    stat = sqrt(stat / r) * signs
    pvals = 2*pnorm(-abs(stat))
    fdrs = p.adjust(pvals, method="BH")

    ## return a data frame
    gr = DMLfit$gr
    res = data.frame(chr=seqnames(gr), pos=start(gr), stat, pvals, fdrs)
    attr(res, "Contrast") = Contrast
    invisible(res)

}


##############################################################
## make contrast matrix given a model and a term to be tested
## Input term can be a vector (testing multiple terms)
##############################################################
makeContrast <- function(fit, term) {
    formula.terms = attr(terms(fit$formula), "term.labels")
    ix = match(term, formula.terms)
    if( any(is.na(ix)) )
        stop("Some term(s) to be tested can't be found in the formula.\n")

    ## make contrast matrix. All columns in the design matrix related to
    ## the provided term (including interactions) should be tested.
    allcolnam = colnames(fit$X)
    ix.term = NULL
    ixcol = attr(fit$X, "assign")
    for( t in term ) {
        iii = which(formula.terms == t)
        ix.term = c( ix.term, which(ixcol==iii) )
    }

    ## make matrix.
    L = matrix(0, ncol=ncol(fit$X), nrow=length(ix.term))
    for(i in 1:nrow(L))
        L[i, ix.term[i]] = 1

    return(t(L))
}


