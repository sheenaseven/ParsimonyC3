#' Add transformed weights into cci
#'
#' @param cci The dataframe extracted from  DataTransformCPD function
#' @param features The selected features that can be included into the final weights, default: c('means_e','FC_ec')
#' @param w The weights controls features' contribution, default: c(1,1)
#'
#' @return The new cci with added weights
#'
#' @examples
#' \donttest{
#' cci <- AddWeights(cci,features=c('means_e','FC_ec'),w = c(1,1))
#' }
#'
#' @import dplyr
#' @export

AddWeights = function(
    cci,
    features=c('means_e','FC_ec'),
    w = c(1,1)
) {
  if ( !all(features %in% colnames(cci)) ) {
    print('error: columns "means_e","FC_ec" not in cci')
  } else {
    for (i in seq_along(features)) {
      column_name <- paste0("score_", features[i])
      cci[column_name] <- w[i]*WeightTransform(cci[features[i]])
    }
    if (length(features)==1) {
      cci$weight <- cci[, paste0("score_", features)]
    } else {
      cci$weight <- rowSums(cci[, paste0("score_", features)])
    }
    cci$negweight <- -1*cci$weight
    cci$weight_d <- (cci$negweight - min(cci$negweight)) / (max(cci$negweight)-min(cci$negweight))
    return(cci)
  }
}
