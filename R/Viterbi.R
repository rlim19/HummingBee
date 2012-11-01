Viterbi <- function (Q, initialProb, emissionProb, blockSizes = NULL) {

   n <- nrow(emissionProb)
   m <- nrow(Q)
   ViterbiPath <- integer(n)
   if(is.null(blockSizes))
      blockSizes <- n

   emissionProb <- log(emissionProb)
   Q <- log(Q)

   emissionProb[is.infinite(emissionProb)] <- -325
   Q[is.infinite(Q)] <- -325
   # lower than the smallest value R can display

   cumulative.n <- 0

   for (n.i in blockSizes) {

      viterbi <- .Fortran("vit", n.i, m, log(initialProb),
         emissionProb[(cumulative.n + 1):(cumulative.n + n.i),], Q,
         integer(n.i), matrix(double(n.i * m), nrow = n.i), package = "HummingBee")

      ViterbiPath[(cumulative.n + 1):(cumulative.n + n.i)] <- viterbi[[6]]
      cumulative.n <- cumulative.n + n.i

   }

   return(ViterbiPath)

}
