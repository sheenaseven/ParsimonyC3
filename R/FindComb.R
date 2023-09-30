#' find the ligand-receptor(LR) combinations
#'
#' @param cci The dataframe extracted from  DataTransformCPD function
#' @param n Maximum number of LR pairs in one combination, default:2
#'
#' @return the ligand-receptor(LR) combination list
#'
#' @examples
#' \donttest{
#' lr_combination <- FindComb(cci,n=2)
#' }
#'
#' @import dplyr
#' @export
FindComb = function(
    cci,
    n = 2
) {
  if ('inter_lr' %in% colnames(cci)) {
    lr <- unique(cci$inter_lr)
    if (n==1) {
      lr_combination <- combn(lr, 1)
      print(paste0('Total LR combinations: ', dim(lr_combination)[2]))
    } else if (n>1) {
      print('Current combination: n=1')
      lr_combination <- combn(lr, 1)
      num <- dim(lr_combination)[2]
      for (i in 2:n) {
        print(paste0('Current combination: n=',i))
        combni <- as.list(as.data.frame(combn(lr, i)))
        num <- num+dim(combn(lr, i))[2]
        lr_combination <- c(lr_combination,combni)
      }
      print(paste0('Total LR combinations: ', num))
    }
    else{
      print('error: n should be > 1')
    }
  }
  else{
    print('error: there is no inter_lr column in your input cci, please check it!')
  }
  return(lr_combination)
}
