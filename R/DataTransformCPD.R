#' merge the expr and ctrl CellPhoneDB data and get the differential cross-talking cell clusters
#'
#'
#' This function extracts the results from CellPhoneDB; merges the expr and ctrl results based
#' on cell-cell L-R interactions; calculates the foldchange (expr/ctrl) based on the LR expression;
#' and get the differential cross-talking cell clusters in the expr group compared to ctrl.
#'
#' @param exp_dir CellPhoneDB output directory of expr group, like cellphonedb/expr/out
#' @param ct_dir CellPhoneDB output directory of ctrl group, like cellphonedb/ctrl/out
#' @param celltypes the cell types of interest
#' @param SelfInter whether self-self communications are included. Default: FALSE, not included
#' @param FC_cutoff: foldchange (expr/out) cutoff for the differential communication, default:1.5
#' @param pvalue_cutoff: pvalue cutoff (extracted from CellPhoneDB) of expr group, default: 0.05
#' @param FC_constant small constant to avoid issues associated with division by zero, default: 0.01
#'
#' @return A dataframe that includes the differential ligand-receptor interactions in expr vs ctrl
#'
#' @examples
#' \donttest{
#' cci <- DataTransformCPD('cellphonedb/ccm/out','cellphonedb/ep/out',celltypes=c('Tcell','Bcell'))
#' }
#'
#' @import dplyr
#' @import stringr
#' @export

DataTransformCPD <- function(
    exp_dir,
    ct_dir,
    celltypes=c("Vn_ACKR1","Mph_C2_GPNMB",'PB_MZB1',"CD4_C5_FOXP3","Mo_C2_EREG","Mast_TPSB2","FB_C2_POSTN"),
    SelfInter=FALSE,
    FC_cutoff=1.5,
    pvalue_cutoff=0.05,
    FC_constant=0.01
) {
  pvals_e <- read.delim(paste0(exp_dir,'/pvalues.txt'), check.names = FALSE)
  means_e <- read.delim(paste0(exp_dir,'/means.txt'), check.names = FALSE)
  pvals_c <- read.delim(paste0(ct_dir,'/pvalues.txt'), check.names = FALSE)
  means_c <- read.delim(paste0(ct_dir,'/means.txt'), check.names = FALSE)
  fea <- celltypes
  means_e %>% dplyr::select("interacting_pair",starts_with(fea),ends_with(fea))  %>%
    reshape2::melt() -> meansdf
  pvals_e %>% dplyr::select("interacting_pair",starts_with(fea),ends_with(fea))  %>%
    reshape2::melt() -> pvalsdf
  colnames(meansdf) <- c("interacting_pair","CC","means")
  colnames(pvalsdf) <- c("interacting_pair","CC","pvals")
  pvalsdf$id <- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$id <- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf_e <- merge(pvalsdf,meansdf,by = "id")
  colnames(pldf_e) <- c("id",paste0(colnames(pldf_e)[2:length(colnames(pldf_e))],"_","e"))

  means_c %>% dplyr::select("interacting_pair",starts_with(fea),ends_with(fea))  %>%
    reshape2::melt() -> meansdf
  colnames(meansdf)<- c("interacting_pair","CC","means")
  pvals_c %>% dplyr::select("interacting_pair",starts_with(fea),ends_with(fea))  %>%
    reshape2::melt()-> pvalsdf
  colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
  pvalsdf$id<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$id<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf_c <- merge(pvalsdf,meansdf,by = "id")
  colnames(pldf_c) <- c("id",paste0(colnames(pldf_c)[2:length(colnames(pldf_c))],"_","c"))

  pldf <- merge(pldf_e,pldf_c,by = "id", all = TRUE)
  a <- FC_constant
  pldf$FC_ec <- (pldf$means_e+a)/(pldf$means_c+a)
  pldf_f <- pldf %>% filter((pvals_e<pvalue_cutoff) & (FC_ec>FC_cutoff | is.na(pvals_c)))

  pldf_f[c('from0', 'to0')] <- stringr::str_split_fixed(as.character(pldf_f$CC.x_e), '\\|', 2)
  pldf_f <- pldf_f[pldf_f$from0 %in% fea,]
  pldf_f <- pldf_f[pldf_f$to0 %in% fea,]
  pldf_f <- pldf_f[pldf_f$from0!=pldf_f$to0,]

  lr_rl=pvals_e[,c("interacting_pair","receptor_b","secreted","is_integrin")]
  colnames(lr_rl)=c("interacting_pair.x_e","receptor_b","secreted","is_integrin")
  pldf_f=merge(pldf_f,lr_rl,by = "interacting_pair.x_e",all.x = TRUE)
  pldf_f[c('L0', 'R0')] <- stringr::str_split_fixed(as.character(pldf_f$interacting_pair.x_e), '_', 2)
  pldf_f1=pldf_f[pldf_f$receptor_b=='True',]
  pldf_f1$L=pldf_f1$L0;pldf_f1$R=pldf_f1$R0;pldf_f1$from=pldf_f1$from0;pldf_f1$to=pldf_f1$to0
  pldf_f2=pldf_f[pldf_f$receptor_b=='False',]
  pldf_f2$L=pldf_f2$R0;pldf_f2$R=pldf_f2$L0;pldf_f2$from=pldf_f2$to0;pldf_f2$to=pldf_f2$from0
  pldf_m=rbind(pldf_f1,pldf_f2)
  if (SelfInter==F){
    pldf_m <- pldf_m[pldf_m$L!=pldf_m$R,]
  }
  pldf_m$inter_lr=paste0(pldf_m$L,'_',pldf_m$R)
  pldf_m$inter_cc=paste0(pldf_m$from,'|',pldf_m$to)
  return(pldf_m)
}
