mlNB <- function(x, tol=1e-6) {
   # Does maximum likelihood estimation of negative binomial
   # parameters for sample 'x'.

   n <- length(x)
   m <- mean(x)
   tab <- tabulate(x+1L)
   u <- 0:(length(tab)-1L)
   a <- 1
   new.a <- a + 2*tol

   while (abs(new.a - a) > tol) {
      a <- new.a
      num <- sum(tab*digamma(a+u)) - n*digamma(a) + n*log(a/(m+a))
      denom <- sum(tab*trigamma(a+u)) - n*trigamma(a) + n*m/(a*(m+a))
      new.a <- a - num/denom
   }

   return(c(a=new.a, p=a/(m+a)))

}
