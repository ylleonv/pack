summary.pcglm <- function(object, ...) {
  coef <- object$coefficients
  se   <- object$stderr
  tval <- coef/se

  object$coefficients <- cbind("Estimate"     = coef,
                               "Std. Error" = se,
                               "z value"    = tval,
                               "Pr(>|z|)"   = 2*pnorm(-abs(tval)))
  colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(object$coefficients, P.values=TRUE, has.Pvalue=TRUE, ...)
  # cf src/stats/R/lm.R and case with no weights and an intercept
  # f <- object$fitted.values
  # r <- object$residuals
  #mss <- sum((f - mean(f))^2)
  # mss <- if (object$intercept) sum((f - mean(f))^2) else sum(f^2)
  # rss <- sum(r^2)
  #
  # object$r.squared <- mss/(mss + rss)
  # df.int <- if (object$intercept) 1L else 0L
  # n <- length(f)
  # rdf <- object$df
  # object$adj.r.squared <- 1 - (1 - object$r.squared) * ((n - df.int)/rdf)
  class(object) <- "summary.pcglm"
  # object
}


