function(x, w, tol=1e-12) {
   # Does maximum likelihood estimation of negative binomial
   # parameters for sample 'x', with weights 'w'.

   x <- x[!is.na(x)]

   n <- length(x)
   m <- mean(x)
   tab <- tabulate(x+1L)
   u <- which(tab > 0) - 1L
   tab <- tab[tab > 0]
   x. <- sum(x)
   a <- 1
   b <- 1

   f <- sum(tab*digamma(a+u)) - n*digamma(a) - sum(log(1+w*b))
   g <- x./b - sum((a+x)*w/(1+w*b))
   for (i in 1:1000) {
      if ((oldgrad <- f*f + g*g) < tol) break
      dfda <- sum(tab*trigamma(a+u)) - n*trigamma(a)
      dfdb <- dgda <- -sum(w/(1+w*b))
      dgdb <- -x./b^2 + sum((a+x)*w^2/(1+w*b)^2)
      denom <- dfdb*dgda - dfda*dgdb
      da = (f*dgdb - g*dfdb) / denom
      db = (g*dfda - f*dgda) / denom
      f <- sum(tab*digamma(a+da+u)) - n*digamma(a+da) - sum(log(1+w*(b+db)))
      g <- x./(b+db) - sum((a+da+x)*w/(1+w*(b+db)))
      for (i in 1:20) {
         if (f*f + g*g < oldgrad) break
         da <- da/2
         db <- db/2
         f <- sum(tab*digamma(a+da+u)) - n*digamma(a+da) -
            sum(log(1+w*(b+db)))
         g <- x./(b+db) - sum((a+da+x)*w/(1+w*(b+db)))
      }
      a <- a+da
      b <- b+db
   }

   return (sum(tab*lgamma(a+u)) - n*lgamma(a) -
      sum((a+x)*log(1+w*b) - x*log(w*b)))

}
