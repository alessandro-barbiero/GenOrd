#' @importFrom stats pt integrate
#' @export
p_rect_t <- function(ax, bx, ay, by, rho, df, rel.tol = 1e-4) {

  tdiff <- function(lo, hi, df) {
    p_hi <- ifelse(is.infinite(hi) & hi > 0, 1, pt(hi, df))
    p_lo <- ifelse(is.infinite(lo) & lo < 0, 0, pt(lo, df))
    pmax(p_hi - p_lo, 0)
  }

  ua <- if (is.infinite(ax) && ax < 0) 0 else pt(ax, df)
  ub <- if (is.infinite(bx) && bx > 0) 1 else pt(bx, df)

  eps <- 1e-8
  ua <- max(ua, eps)
  ub <- min(ub, 1 - eps)

  if (ub <= ua) return(0)

  g <- function(u) {
    x <- qt(u, df)

    sx <- sqrt((1 - rho^2) * (df + x^2) / (df + 1))

    ly <- (ay - rho * x) / sx
    uy <- (by - rho * x) / sx

    py <- tdiff(ly, uy, df + 1)
    py[!is.finite(py)] <- 0
    py
  }

  val <- integrate(
    g,
    lower = ua,
    upper = ub,
    rel.tol = rel.tol,
    subdivisions = 1000
  )$value

  max(val, 0)
}
