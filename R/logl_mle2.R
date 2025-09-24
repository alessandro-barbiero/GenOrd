#' @importFrom cubature cuhre
#' @export

# another version of the logL function accepting rho and df as distinct scalar parameters
logL_mle2 <- function(rho, df, x) logL(c(rho, df), x)

