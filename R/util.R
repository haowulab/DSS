###############################################
### some utility functions
###############################################

### generate negative binomial rv, given mu and phi (over-dispersion)
rnegbinom <- function (n, mu =1, phi=0.01){
    rpois(n, rgamma(n, shape=1/phi,scale=mu*phi))
}


### Smoothing function. Smooth by chr
smooth.chr <- function(x, ws, allchr, allpos, method=c("avg", "sum")) {
    method <- match.arg(method)
    if(method == "avg")
        flag = 1
    else
      flag = 0

    ## remove NA's
    ix = !is.na(x)
    x2 = x[ix]
    allchr2 = allchr[ix]
    allpos2 = allpos[ix]

    n0=length(x2)
    idx=split(1:n0, allchr2)
    res2=rep(0, n0)
    for(i in seq(along=idx)) {
        res2[idx[[i]]]=.Call("windowFilter", x2[idx[[i]]], as.integer(allpos2[idx[[i]]]), as.integer(ws), as.integer(flag))
    }

    if(sum(ix) > 0) {
        res = rep(NA, length(x))
        res[ix] = res2
    } else res = res2

    return(res)
}

