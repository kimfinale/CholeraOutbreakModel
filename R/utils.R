# final epidemic size
# the simplest SIR model
# R = 1 - exp(-R0*R) where R is the final epidemic size (or R(\infty) for the SIR model)
final_epidemic_size <- function(R0 = 2) {
  y = function(x) x - 1 + exp(-R0*x)
  final_size <- uniroot(y, interval=c(1e-6,1-1e-6))$root

  return(final_size)

}
