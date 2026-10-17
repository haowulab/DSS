# Six samples. The unweighted design is full rank. Site 1 coverage makes
# the weighted cross-product singular; site 2 is an ordinary CpG.
# Z is the arcsine transform DMLfit.multiFactor passes to the engine.

local_sites <- function() {
  X <- cbind(1, c(0, 0, 0, 1, 1, 1))
  N <- rbind(
    c(1e16, 1e16, 1e16, 1, 1, 1),
    rep(20, 6)
  )
  Y <- rbind(
    c(1e15, 2e15, 3e15, 0, 1, 1),
    c(2, 4, 6, 8, 10, 12)
  )
  c0 <- 0.1
  Z <- asin(2 * (Y + c0) / (N + 2 * c0) - 1)
  list(X = X, Y = Y, N = N, Z = Z, n = nrow(X), p = ncol(X))
}

test_that("solve() still throws on the singular weighted cross-product", {
  dat <- local_sites()
  X <- dat$X
  N1 <- dat$N[1, ]

  expect_equal(qr(X)$rank, ncol(X))
  expect_true(is.matrix(solve(t(X) %*% X)))

  weighted <- t(X * N1) %*% X
  expect_error(solve(weighted), "singular")

  expect_true(all(dat$Y >= 0 & dat$Y <= dat$N))
  expect_true(is.matrix(solve(t(X * dat$N[2, ]) %*% X)))
})

test_that("a singular CpG is kept as a missing test", {
  dat <- local_sites()
  X <- dat$X
  p <- dat$p

  one_bad <- DSS:::DMLfit.oneCG(
    dat$Y[1, ], dat$N[1, ], X, dat$Z[1, ], dat$n, p
  )
  expect_false(is.null(one_bad))
  expect_identical(dim(one_bad$beta0), c(p, 1L))
  expect_length(one_bad$se.beta0, p)
  expect_length(one_bad$var.beta0, p * p)
  expect_identical(one_bad$phi, NA_real_)
  expect_true(all(!is.finite(one_bad$beta0)))
  expect_true(all(!is.finite(one_bad$se.beta0)))
  expect_true(all(!is.finite(one_bad$var.beta0)))

  fit <- expect_output(
    DSS:::DMLfit.multiFactor.engine(dat$Y, dat$N, X, dat$Z),
    "Fitting DML model"
  )
  expect_equal(nrow(fit$beta), 2L)
  expect_equal(nrow(fit$var.beta), 2L)
  expect_equal(length(fit$phi), 2L)
  expect_true(all(!is.finite(fit$beta[1, ])))
  expect_true(all(!is.finite(fit$var.beta[1, ])))
  expect_identical(fit$phi[1], NA_real_)

  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1L, 2L),
    end = c(1L, 2L)
  ))
  tested <- DSS::DMLtest.multiFactor(list(gr = gr, X = X, fit = fit), coef = 2)
  expect_equal(nrow(tested), 2L)
  expect_true(!is.finite(tested$pvals[1]))
  expect_true(is.finite(tested$pvals[2]))

  one_ok <- DSS:::DMLfit.oneCG(
    dat$Y[2, ], dat$N[2, ], X, dat$Z[2, ], dat$n, p
  )
  expect_equal(as.numeric(fit$beta[2, ]), as.numeric(one_ok$beta0), tolerance = 1e-8)
  expect_equal(as.numeric(fit$var.beta[2, ]), as.numeric(one_ok$var.beta0), tolerance = 1e-8)
  expect_equal(fit$phi[2], one_ok$phi, tolerance = 1e-8)
  se_engine <- sqrt(diag(matrix(fit$var.beta[2, ], ncol = p)))
  expect_equal(se_engine, as.numeric(one_ok$se.beta0), tolerance = 1e-8)
})

test_that("a non-singular oneCG failure still aborts the fit", {
  dat <- local_sites()
  Nbad <- dat$N
  Nbad[1, 1] <- NA_real_

  err <- expect_output(
    tryCatch(
      DSS:::DMLfit.multiFactor.engine(dat$Y, Nbad, dat$X, dat$Z),
      error = function(e) e
    ),
    "Fitting DML model"
  )
  expect_s3_class(err, "error")
  expect_no_match(conditionMessage(err), "singular")
  ref <- tryCatch(if (NA) 1, error = function(e) e)
  expect_identical(conditionMessage(err), conditionMessage(ref))
  expect_identical(class(err), class(ref))
})
