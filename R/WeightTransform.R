#' weight transformation methods
#'
#' @param x input weight vector
#' @param transform transformation methods, default:'neg_log_normal': -1*log10(x/max(x));
#' 'neg_log': -1*log10(x);'exp_neg_log':exp(1)^(-1*log10(x));'log':log10(x)
#'
#' @return transformed weight vector
#'
#' @examples
#' \donttest{
#'
#' }
#'
#' @export

WeightTransform=function(
    x,
    transform='neg_log_normal'
) {
  if (transform=='neg_log') {
    y <- -1*log10(x)
  }
  if (transform=='exp_neg_log') {
    y <- exp(1)^(-1*log10(x))
  }
  if (transform=='log') {
    y <- log10(x)
  }
  if (transform=='neg_log_normal') {
    y <- -1*log10( x/max(x) )
  }
  return(y)
}
